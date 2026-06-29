/* dwd_core.c -- distributed qutrit weight-enumerator hot core (C + OpenMP).
 *
 * Computes a partial Hamming-weight histogram (int64 weight enumerator) of an
 * [n, K] code over F_3 by enumerating a contiguous range of message *blocks*.
 *
 * NOTE ON LANGUAGE: the frozen design (DISTRIBUTED_WEIGHTDIST_DESIGN.md sec.1)
 * names Rust+rayon as primary with "C+OpenMP the documented fallback ... if a
 * node lacks a Rust toolchain".  This build host has gcc 13.3 + OpenMP but NO
 * Rust toolchain, so the fallback is the actual core.  Same algorithm (loopless
 * reflected ternary Gray code + trit-packed full-add), same .g0 / .partial
 * byte formats, so it interoperates with the Python harness identically.
 *
 * ----------------------------------------------------------------------------
 * ALGORITHM
 * ----------------------------------------------------------------------------
 * Message space F_3^K is partitioned by fixing the TOP t message trits -> 3^t
 * independent blocks, each enumerating the remaining Kfree = K-t trits (3^Kfree
 * messages).  Free rows = G rows [0 .. Kfree).  Block rows = G rows [Kfree .. K).
 * A block id b in [0, 3^t) has base-3 digits (v_0..v_{t-1}) giving the fixed
 * values of the block rows; the block's BASE codeword is
 *      c_base = sum_i v_i * G[Kfree+i]   (mod 3).
 * Within a block we enumerate all 3^Kfree combinations of the free rows with a
 * loopless reflected ternary Gray code (Knuth TAOCP 7.2.1.1 Alg.H): consecutive
 * messages differ in exactly ONE free trit by +-1, so the codeword changes by
 * + or - exactly one generator row.  Weight is recomputed each step word-parallel
 * over the packed codeword (cost ~ n/32, independent of row weight -- the dense
 * RM_3(4,7) rows make this ~22x faster than a support scan, measured).
 *
 * Histogram counts are int64.  Every bin <= 3^K <= 3^29 << INT64_MAX, so no add
 * can overflow (proof: correctness.assert_int64_safe).  No floating point.
 *
 * ----------------------------------------------------------------------------
 * PACKING
 * ----------------------------------------------------------------------------
 * The codeword is stored as one byte per trit (value 0/1/2), length padded up to
 * a multiple of 32 (np).  Padding bytes are 0 and every generator row is 0 in the
 * padding, so packed adds keep them 0 and they never contribute to the weight.
 * The "trit-packing" is the 32-trits-per-256-bit-AVX2-lane full add: add g (mod 3)
 * to the whole codeword and derive the nonzero count by popcount of the
 * is-zero byte mask.  A scalar fallback (no AVX2) is compiled when __AVX2__ is
 * absent; `dwd_core selfcheck` asserts packed == scalar bit-exact (design sec.4.2).
 *
 * ----------------------------------------------------------------------------
 * CLI
 * ----------------------------------------------------------------------------
 *   dwd_core enum <g0> <t> <block_start> <block_count> <out.partial> [nthreads]
 *   dwd_core selfcheck                 # packed-vs-scalar bit-exact self test
 *
 * Build: see build.sh  (gcc -O3 -march=native -fopenmp).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif

#define KMAX 32

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }

/* ---- .g0 parser -------------------------------------------------------- */
typedef struct {
    int q, K, n, np;          /* np = n padded up to multiple of 32 */
    uint8_t *G;               /* K*np bytes, row-major, padding=0 */
    uint8_t *Gneg;            /* K*np bytes, (3-g)%3 == -g, padding=0 */
} g0_t;

static uint8_t* row(uint8_t*base,int np,int r){ return base + (size_t)r*np; }

static int read_g0(const char*path, g0_t*g){
    FILE*fp=fopen(path,"rb");
    if(!fp){ fprintf(stderr,"dwd_core: cannot open %s\n",path); return 1; }
    char hdr[256];
    if(!fgets(hdr,sizeof hdr,fp)){ fprintf(stderr,"empty .g0\n"); fclose(fp); return 1; }
    char magic[64],ver[64];
    if(sscanf(hdr,"%63s %63s %d %d %d",magic,ver,&g->q,&g->K,&g->n)!=5
       || strcmp(magic,"qmsd-g0") || strcmp(ver,"v1")){
        fprintf(stderr,"bad .g0 header: %s\n",hdr); fclose(fp); return 1;
    }
    if(g->q!=3){ fprintf(stderr,"dwd_core requires q=3 (got %d)\n",g->q); fclose(fp); return 1; }
    if(g->K<1||g->K>KMAX){ fprintf(stderr,"K=%d out of [1,%d]\n",g->K,KMAX); fclose(fp); return 1; }
    g->np=(g->n+31)&~31;
    g->G   =aligned_alloc(32,(size_t)g->K*g->np);
    g->Gneg=aligned_alloc(32,(size_t)g->K*g->np);
    memset(g->G,0,(size_t)g->K*g->np);
    memset(g->Gneg,0,(size_t)g->K*g->np);
    char *line=malloc(g->n+8);
    for(int r=0;r<g->K;r++){
        if(!fgets(line,g->n+8,fp)){ fprintf(stderr,"row %d missing\n",r); fclose(fp); return 1; }
        /* strip trailing newline/cr */
        int L=(int)strlen(line);
        while(L>0 && (line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0;
        if(L!=g->n){ fprintf(stderr,"row %d len %d != n=%d\n",r,L,g->n); fclose(fp); return 1; }
        uint8_t *gr=row(g->G,g->np,r), *gn=row(g->Gneg,g->np,r);
        for(int p=0;p<g->n;p++){
            int v=line[p]-'0';
            if(v<0||v>=3){ fprintf(stderr,"row %d col %d bad digit\n",r,p); fclose(fp); return 1; }
            gr[p]=(uint8_t)v; gn[p]=(uint8_t)((3-v)%3);
        }
    }
    free(line); fclose(fp);
    return 0;
}

/* ---- packed full-add weight: c += g (mod 3) over np bytes, return wt ----- */
static inline long packed_add_weight(uint8_t*restrict c, const uint8_t*restrict g, int np){
    long wt=0;
#ifdef __AVX2__
    const __m256i three=_mm256_set1_epi8(3), two=_mm256_set1_epi8(2), zero=_mm256_setzero_si256();
    for(int b=0;b<np;b+=32){
        __m256i cv=_mm256_load_si256((const __m256i*)(c+b));
        __m256i gv=_mm256_load_si256((const __m256i*)(g+b));
        __m256i sv=_mm256_add_epi8(cv,gv);                 /* 0..4 */
        __m256i ge=_mm256_cmpgt_epi8(sv,two);              /* signed ok (<128) */
        sv=_mm256_sub_epi8(sv,_mm256_and_si256(ge,three)); /* mod 3 */
        _mm256_store_si256((__m256i*)(c+b),sv);
        unsigned m=(unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(sv,zero));
        wt += 32 - __builtin_popcount(m);
    }
#else
    for(int p=0;p<np;p++){
        int nv=c[p]+g[p]; if(nv>=3) nv-=3; c[p]=(uint8_t)nv; wt += (nv!=0);
    }
#endif
    return wt;
}

/* ---- enumerate one block into a thread-private histogram ----------------- */
/* free rows are G rows [0..Kfree); block id b fixes rows [Kfree..K). */
static void enum_block(const g0_t*g, int Kfree, long long block,
                       uint8_t*restrict c, long*restrict lh){
    int np=g->np, t=g->K-Kfree;
    memset(c,0,np);
    /* base codeword from the fixed block rows */
    long long bb=block;
    for(int i=0;i<t;i++){
        int v=(int)(bb%3); bb/=3;
        if(v){ const uint8_t*gr=row(g->G,np,Kfree+i);
            if(v==1) for(int p=0;p<g->n;p++){ int nv=c[p]+gr[p]; if(nv>=3)nv-=3; c[p]=(uint8_t)nv; }
            else     for(int p=0;p<g->n;p++){ int nv=c[p]+2*gr[p]; while(nv>=3)nv-=3; c[p]=(uint8_t)nv; }
        }
    }
    long w0=0; for(int p=0;p<g->n;p++) w0 += (c[p]!=0);
    lh[w0]++;                                   /* the all-free-zero message */
    if(Kfree==0) return;
    /* Knuth loopless reflected ternary Gray code over Kfree free digits */
    int a[KMAX+1], f[KMAX+1], dir[KMAX+1];
    for(int j=0;j<=Kfree;j++){ a[j]=0; f[j]=j; dir[j]=1; }
    for(;;){
        int d=f[0]; f[0]=0; if(d==Kfree) break;
        int s=dir[d]; a[d]+=s;
        if(a[d]==0||a[d]==2){ dir[d]=-dir[d]; f[d]=f[d+1]; f[d+1]=d+1; }
        const uint8_t*gg = (s>0)? row(g->G,np,d) : row(g->Gneg,np,d);
        long wt=packed_add_weight(c,gg,np);
        lh[wt]++;
    }
}

static long long ipow3(int e){ long long r=1; for(int i=0;i<e;i++) r*=3; return r; }

/* ---- .partial writer (format DWDP0001, little-endian) -------------------- */
static int write_partial(const char*path,const g0_t*g,
                         long long block_start,long long block_count,
                         const long*hist){
    char tmp[4096]; snprintf(tmp,sizeof tmp,"%s.tmp",path);
    FILE*fp=fopen(tmp,"wb");
    if(!fp){ fprintf(stderr,"cannot write %s\n",tmp); return 1; }
    uint64_t checksum=0; for(int w=0;w<=g->n;w++) checksum+=(uint64_t)hist[w];
    uint32_t q=g->q,K=g->K,n=g->n,nblocks=(uint32_t)block_count;
    fwrite("DWDP0001",1,8,fp);
    fwrite(&q,4,1,fp); fwrite(&K,4,1,fp); fwrite(&n,4,1,fp); fwrite(&nblocks,4,1,fp);
    fwrite(&checksum,8,1,fp);
    for(int w=0;w<=g->n;w++){ int64_t v=hist[w]; fwrite(&v,8,1,fp); }
    for(long long b=0;b<block_count;b++){ uint32_t id=(uint32_t)(block_start+b); fwrite(&id,4,1,fp); }
    fflush(fp);
#ifdef __unix__
    int fd=fileno(fp); if(fd>=0){ extern int fsync(int); fsync(fd); }
#endif
    fclose(fp);
    if(rename(tmp,path)){ fprintf(stderr,"rename %s -> %s failed\n",tmp,path); return 1; }
    return 0;
}

static int cmd_enum(int argc,char**argv){
    if(argc<7){ fprintf(stderr,
        "usage: dwd_core enum <g0> <t> <block_start> <block_count> <out.partial> [nthreads]\n");
        return 2; }
    const char*g0path=argv[2];
    int t=atoi(argv[3]);
    long long block_start=atoll(argv[4]);
    long long block_count=atoll(argv[5]);
    const char*outpath=argv[6];
    int nthreads = (argc>7)? atoi(argv[7]) : 0;
#ifdef _OPENMP
    if(nthreads>0) omp_set_num_threads(nthreads);
#endif
    g0_t g;
    if(read_g0(g0path,&g)) return 1;
    if(t<0||t>g.K){ fprintf(stderr,"t=%d out of [0,K=%d]\n",t,g.K); return 1; }
    int Kfree=g.K-t;
    long long nblk=ipow3(t);
    if(block_start<0||block_count<0||block_start+block_count>nblk){
        fprintf(stderr,"block range [%lld,%lld) out of [0,3^%d=%lld)\n",
                block_start,block_start+block_count,t,nblk); return 1; }

    long *hist=calloc((size_t)g.n+1,sizeof(long));
    double t0=now();
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        long *lh=calloc((size_t)g.n+1,sizeof(long));
        uint8_t *c=aligned_alloc(32,g.np);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
#endif
        for(long long blk=block_start; blk<block_start+block_count; blk++){
            enum_block(&g,Kfree,blk,c,lh);
        }
#ifdef _OPENMP
        #pragma omp critical
#endif
        { for(int w=0;w<=g.n;w++) hist[w]+=lh[w]; }
        free(lh); free(c);
    }
    double dt=now()-t0;

    long long msgs=block_count*ipow3(Kfree);
    uint64_t chk=0; for(int w=0;w<=g.n;w++) chk+=(uint64_t)hist[w];
    if((long long)chk != msgs){
        fprintf(stderr,"FATAL checksum %llu != expected messages %lld\n",
                (unsigned long long)chk,msgs); return 1; }
    if(write_partial(outpath,&g,block_start,block_count,hist)) return 1;
    fprintf(stderr,"[dwd_core] K=%d n=%d t=%d Kfree=%d blocks=[%lld,%lld) msgs=%lld "
            "%.2fs %.3e msg/s -> %s\n",
            g.K,g.n,t,Kfree,block_start,block_start+block_count,msgs,dt,
            dt>0?msgs/dt:0.0,outpath);
    free(hist); free(g.G); free(g.Gneg);
    return 0;
}

/* ---- selfcheck: packed enumeration == scalar brute, on a small code ------ */
static uint64_t rs=88172645463325252ULL;
static inline uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

static int cmd_selfcheck(void){
    int fails=0;
    for(int trial=0; trial<6; trial++){
        int K=4+(int)(xr()%6);          /* 4..9 */
        int n=20+(int)(xr()%40);        /* 20..59 */
        int np=(n+31)&~31;
        g0_t g={3,K,n,np,NULL,NULL};
        g.G=aligned_alloc(32,(size_t)K*np); g.Gneg=aligned_alloc(32,(size_t)K*np);
        memset(g.G,0,(size_t)K*np); memset(g.Gneg,0,(size_t)K*np);
        for(int r=0;r<K;r++){ uint8_t*gr=row(g.G,np,r),*gn=row(g.Gneg,np,r);
            for(int p=0;p<n;p++){ int v=(int)(xr()%3); gr[p]=v; gn[p]=(3-v)%3; } }

        /* packed gray enumeration over all 3^K (t=0, single block) */
        long *hp=calloc((size_t)n+1,sizeof(long));
        uint8_t *c=aligned_alloc(32,np);
        enum_block(&g,K,0,c,hp);
        free(c);

        /* scalar brute: iterate every message, full m@G mod 3, count nonzeros */
        long *hs=calloc((size_t)n+1,sizeof(long));
        long long total=ipow3(K);
        uint8_t *cw=malloc(n);
        for(long long m=0;m<total;m++){
            memset(cw,0,n);
            long long mm=m;
            for(int r=0;r<K;r++){ int dgt=(int)(mm%3); mm/=3;
                if(dgt){ uint8_t*gr=row(g.G,np,r);
                    for(int p=0;p<n;p++){ int nv=cw[p]+dgt*gr[p]; nv%=3; cw[p]=(uint8_t)nv; } } }
            int w=0; for(int p=0;p<n;p++) w+=(cw[p]!=0);
            hs[w]++;
        }
        free(cw);
        int ok=1; for(int w=0;w<=n;w++) if(hp[w]!=hs[w]) ok=0;
        long long sp=0,ss=0; for(int w=0;w<=n;w++){sp+=hp[w];ss+=hs[w];}
        printf("[selfcheck] trial %d K=%d n=%d  packed==scalar:%s  sum=%lld(=3^%d=%lld)\n",
               trial,K,n, ok?"YES":"NO!!", sp, K, total);
        if(!ok || sp!=total || ss!=total) fails++;
        free(hp); free(hs); free(g.G); free(g.Gneg);
    }
    if(fails){ printf("SELFCHECK FAILED (%d)\n",fails); return 1; }
    printf("SELFCHECK PASSED: packed gray-code weight == scalar brute on all trials\n");
    return 0;
}

int main(int argc,char**argv){
    if(argc>=2 && !strcmp(argv[1],"enum"))      return cmd_enum(argc,argv);
    if(argc>=2 && !strcmp(argv[1],"selfcheck")) return cmd_selfcheck();
    fprintf(stderr,
      "dwd_core -- distributed qutrit weight-enumerator hot core\n"
      "  dwd_core enum <g0> <t> <block_start> <block_count> <out.partial> [nthreads]\n"
      "  dwd_core selfcheck\n");
    return 2;
}
