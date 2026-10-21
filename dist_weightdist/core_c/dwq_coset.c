/* dwq_coset.c -- q-ary weight enumerator of a code A and its cosets A + h_i (C + OpenMP).
 *
 * Companion to dwd_core.c (which is q=3 only and left untouched).  Used by
 * qmsd.performance for exact per-output distillation thresholds: one Gray-code walk
 * over the message space F_q^K of A = rowspan(G0) serves every logical row h_i, since at
 * each step the running codeword c is in cache and wt(c + h_i) is a non-destructive
 * add-and-count.
 *
 * ALGORITHM (same partitioning as dwd_core.c)
 *   The top t message digits are fixed -> q^t independent blocks; each block walks the
 *   remaining Kfree = K - t digits with a loopless reflected q-ary Gray code (Knuth
 *   TAOCP 7.2.1.1 Alg. H): consecutive messages differ in one digit by +-1, so c changes
 *   by +-one generator row.  Byte-per-symbol packing, padded to 32; AVX2 add-mod-q
 *   (q <= 127) with an is-zero popcount.  int64 histograms: every bin <= q^K < 2^63.
 *
 * INPUT (.cos text):   line 1: "qmsd-coset v1 q K L n"
 *                      K lines: G0 rows, n space-separated ints in [0,q)
 *                      L lines: logical rows h_1..h_L
 * OUTPUT (.cpart binary, little-endian):
 *   "DWQC0001", u32 q, K, L, n, u64 block_start, block_count, t,
 *   int64 hist[(1+L)*(n+1)]   -- row 0 = A, row i = A + h_i
 *
 * CLI:
 *   dwq_coset enum <in.cos> <t> <block_start> <block_count> <out.cpart> <nthreads>
 *   dwq_coset selfcheck
 * nthreads is REQUIRED (no default) so a run never silently takes every core.
 *
 * Build: gcc -O3 -funroll-loops -march=native -fopenmp -o dwq_coset dwq_coset.c
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

#define KMAX 40

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }

typedef struct {
    int q, K, L, n, np;
    uint8_t *G, *Gneg, *H;   /* K*np, K*np, L*np bytes; padding = 0 */
} cos_t;

static uint8_t* row(uint8_t*base,int np,int r){ return base + (size_t)r*np; }

static uint8_t* alloc_rows(int R,int np){
    size_t sz=(size_t)(R>0?R:1)*np;
    uint8_t*p=aligned_alloc(32,sz); memset(p,0,sz); return p;
}

static int read_cos(const char*path, cos_t*g){
    FILE*fp=fopen(path,"r");
    if(!fp){ fprintf(stderr,"dwq_coset: cannot open %s\n",path); return 1; }
    char magic[64],ver[64];
    if(fscanf(fp,"%63s %63s %d %d %d %d",magic,ver,&g->q,&g->K,&g->L,&g->n)!=6
       || strcmp(magic,"qmsd-coset") || strcmp(ver,"v1")){
        fprintf(stderr,"bad .cos header\n"); fclose(fp); return 1; }
    if(g->q<2||g->q>127){ fprintf(stderr,"q=%d out of [2,127]\n",g->q); fclose(fp); return 1; }
    if(g->K<1||g->K>KMAX){ fprintf(stderr,"K=%d out of [1,%d]\n",g->K,KMAX); fclose(fp); return 1; }
    if(g->L<0){ fprintf(stderr,"L<0\n"); fclose(fp); return 1; }
    g->np=(g->n+31)&~31;
    g->G=alloc_rows(g->K,g->np); g->Gneg=alloc_rows(g->K,g->np); g->H=alloc_rows(g->L,g->np);
    for(int r=0;r<g->K+g->L;r++){
        uint8_t*dst = r<g->K ? row(g->G,g->np,r) : row(g->H,g->np,r-g->K);
        for(int p=0;p<g->n;p++){
            int v; if(fscanf(fp,"%d",&v)!=1||v<0||v>=g->q){
                fprintf(stderr,"bad entry row %d col %d\n",r,p); fclose(fp); return 1; }
            dst[p]=(uint8_t)v;
            if(r<g->K) row(g->Gneg,g->np,r)[p]=(uint8_t)((g->q-v)%g->q);
        }
    }
    fclose(fp);
    return 0;
}

/* c += g (mod q) in place; return weight of the new c */
static inline long add_weight(uint8_t*restrict c,const uint8_t*restrict g,int np,int q){
    long wt=0;
#ifdef __AVX2__
    const __m256i qv=_mm256_set1_epi8((char)q), qm1=_mm256_set1_epi8((char)(q-1)),
                  zero=_mm256_setzero_si256();
    for(int b=0;b<np;b+=32){
        __m256i s=_mm256_add_epi8(_mm256_load_si256((const __m256i*)(c+b)),
                                  _mm256_load_si256((const __m256i*)(g+b)));
        s=_mm256_sub_epi8(s,_mm256_and_si256(_mm256_cmpgt_epi8(s,qm1),qv));
        _mm256_store_si256((__m256i*)(c+b),s);
        wt += 32 - __builtin_popcount((unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(s,zero)));
    }
#else
    for(int p=0;p<np;p++){ int v=c[p]+g[p]; if(v>=q) v-=q; c[p]=(uint8_t)v; wt+=(v!=0); }
#endif
    return wt;
}

/* weight of (c + h) mod q, c unchanged */
static inline long coset_weight(const uint8_t*restrict c,const uint8_t*restrict h,int np,int q){
    long wt=0;
#ifdef __AVX2__
    const __m256i qv=_mm256_set1_epi8((char)q), qm1=_mm256_set1_epi8((char)(q-1)),
                  zero=_mm256_setzero_si256();
    for(int b=0;b<np;b+=32){
        __m256i s=_mm256_add_epi8(_mm256_load_si256((const __m256i*)(c+b)),
                                  _mm256_load_si256((const __m256i*)(h+b)));
        s=_mm256_sub_epi8(s,_mm256_and_si256(_mm256_cmpgt_epi8(s,qm1),qv));
        wt += 32 - __builtin_popcount((unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(s,zero)));
    }
#else
    for(int p=0;p<np;p++){ int v=c[p]+h[p]; if(v>=q) v-=q; wt+=(v!=0); }
#endif
    return wt;
}

static inline void record(const cos_t*g,const uint8_t*c,long wt,long*lh){
    int n1=g->n+1;
    lh[wt]++;
    for(int i=0;i<g->L;i++) lh[(size_t)(i+1)*n1 + coset_weight(c,row(g->H,g->np,i),g->np,g->q)]++;
}

/* free digits = G rows [0..Kfree); block id fixes rows [Kfree..K) */
static void enum_block(const cos_t*g,int Kfree,long long block,uint8_t*restrict c,long*restrict lh){
    int np=g->np, q=g->q, t=g->K-Kfree;
    memset(c,0,np);
    long long bb=block;
    for(int i=0;i<t;i++){
        int v=(int)(bb%q); bb/=q;
        const uint8_t*gr=row(g->G,np,Kfree+i);
        for(int p=0;p<g->n;p++) c[p]=(uint8_t)((c[p]+v*gr[p])%q);
    }
    long w0=0; for(int p=0;p<g->n;p++) w0+=(c[p]!=0);
    record(g,c,w0,lh);
    if(Kfree==0) return;
    int a[KMAX+1], f[KMAX+1], dir[KMAX+1];
    for(int j=0;j<=Kfree;j++){ a[j]=0; f[j]=j; dir[j]=1; }
    for(;;){
        int d=f[0]; f[0]=0; if(d==Kfree) break;
        int s=dir[d]; a[d]+=s;
        if(a[d]==0||a[d]==q-1){ dir[d]=-dir[d]; f[d]=f[d+1]; f[d+1]=d+1; }
        long wt=add_weight(c, s>0 ? row(g->G,np,d) : row(g->Gneg,np,d), np, q);
        record(g,c,wt,lh);
    }
}

static long long ipow(int b,int e){ long long r=1; for(int i=0;i<e;i++) r*=b; return r; }

static int write_cpart(const char*path,const cos_t*g,long long bs,long long bc,int t,const long*hist){
    char tmp[4096]; snprintf(tmp,sizeof tmp,"%s.tmp",path);
    FILE*fp=fopen(tmp,"wb");
    if(!fp){ fprintf(stderr,"cannot write %s\n",tmp); return 1; }
    uint32_t h4[4]={(uint32_t)g->q,(uint32_t)g->K,(uint32_t)g->L,(uint32_t)g->n};
    uint64_t h8[3]={(uint64_t)bs,(uint64_t)bc,(uint64_t)t};
    fwrite("DWQC0001",1,8,fp); fwrite(h4,4,4,fp); fwrite(h8,8,3,fp);
    size_t N=(size_t)(1+g->L)*(g->n+1);
    for(size_t i=0;i<N;i++){ int64_t v=hist[i]; fwrite(&v,8,1,fp); }
    fflush(fp); fclose(fp);
    if(rename(tmp,path)){ fprintf(stderr,"rename failed\n"); return 1; }
    return 0;
}

static int run_enum(const cos_t*g,int t,long long bs,long long bc,long*hist){
    int Kfree=g->K-t;
    size_t N=(size_t)(1+g->L)*(g->n+1);
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        long*lh=calloc(N,sizeof(long));
        uint8_t*c=aligned_alloc(32,g->np);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
#endif
        for(long long b=bs;b<bs+bc;b++) enum_block(g,Kfree,b,c,lh);
#ifdef _OPENMP
        #pragma omp critical
#endif
        { for(size_t i=0;i<N;i++) hist[i]+=lh[i]; }
        free(lh); free(c);
    }
    long long msgs=bc*ipow(g->q,Kfree);
    for(int i=0;i<=g->L;i++){
        long long s=0; for(int w=0;w<=g->n;w++) s+=hist[(size_t)i*(g->n+1)+w];
        if(s!=msgs){ fprintf(stderr,"FATAL row %d sum %lld != %lld\n",i,s,msgs); return 1; }
    }
    return 0;
}

static int cmd_enum(int argc,char**argv){
    if(argc<8){ fprintf(stderr,
        "usage: dwq_coset enum <in.cos> <t> <block_start> <block_count> <out.cpart> <nthreads>\n");
        return 2; }
    int t=atoi(argv[3]); long long bs=atoll(argv[4]), bc=atoll(argv[5]);
    int nth=atoi(argv[7]);
    if(nth<1){ fprintf(stderr,"nthreads must be >= 1\n"); return 2; }
#ifdef _OPENMP
    omp_set_num_threads(nth);
#endif
    cos_t g; if(read_cos(argv[2],&g)) return 1;
    if(t<0||t>g.K){ fprintf(stderr,"t out of range\n"); return 1; }
    long long nblk=ipow(g.q,t);
    if(bs<0||bc<0||bs+bc>nblk){ fprintf(stderr,"block range out of [0,%lld)\n",nblk); return 1; }
    long*hist=calloc((size_t)(1+g.L)*(g.n+1),sizeof(long));
    double t0=now();
    if(run_enum(&g,t,bs,bc,hist)) return 1;
    double dt=now()-t0;
    if(write_cpart(argv[6],&g,bs,bc,t,hist)) return 1;
    long long msgs=bc*ipow(g.q,g.K-t);
    fprintf(stderr,"[dwq_coset] q=%d K=%d L=%d n=%d t=%d blocks=[%lld,%lld) msgs=%lld %.2fs %.3e msg/s (%d threads)\n",
            g.q,g.K,g.L,g.n,t,bs,bs+bc,msgs,dt,dt>0?msgs/dt:0.0,nth);
    return 0;
}

/* selfcheck: Gray enumeration == brute force (all messages, explicit cosets) */
static uint64_t rs=88172645463325252ULL;
static inline uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

static int cmd_selfcheck(void){
    int qs[]={2,3,5,7,11,13}, fails=0;
    for(int trial=0;trial<12;trial++){
        cos_t g; g.q=qs[trial%6]; g.K=2+(int)(xr()%3); g.L=1+(int)(xr()%3); g.n=10+(int)(xr()%40);
        if(ipow(g.q,g.K)>200000) g.K=2;
        g.np=(g.n+31)&~31;
        g.G=alloc_rows(g.K,g.np); g.Gneg=alloc_rows(g.K,g.np); g.H=alloc_rows(g.L,g.np);
        for(int r=0;r<g.K;r++) for(int p=0;p<g.n;p++){ int v=(int)(xr()%g.q);
            row(g.G,g.np,r)[p]=v; row(g.Gneg,g.np,r)[p]=(g.q-v)%g.q; }
        for(int r=0;r<g.L;r++) for(int p=0;p<g.n;p++) row(g.H,g.np,r)[p]=(uint8_t)(xr()%g.q);
        size_t N=(size_t)(1+g.L)*(g.n+1);
        long*hg=calloc(N,sizeof(long)), *hb=calloc(N,sizeof(long));
        int t=g.K>1?1:0;
        if(run_enum(&g,t,0,ipow(g.q,t),hg)) fails++;
        int*cw=malloc(sizeof(int)*g.n);
        for(long long m=0;m<ipow(g.q,g.K);m++){
            long long mm=m; memset(cw,0,sizeof(int)*g.n);
            for(int r=0;r<g.K;r++){ int dg=(int)(mm%g.q); mm/=g.q;
                for(int p=0;p<g.n;p++) cw[p]=(cw[p]+dg*row(g.G,g.np,r)[p])%g.q; }
            for(int i=0;i<=g.L;i++){ int w=0;
                for(int p=0;p<g.n;p++){ int v=i? (cw[p]+row(g.H,g.np,i-1)[p])%g.q : cw[p]; w+=(v!=0); }
                hb[(size_t)i*(g.n+1)+w]++; }
        }
        int ok=!memcmp(hg,hb,N*sizeof(long));
        printf("[selfcheck] q=%d K=%d L=%d n=%d gray==brute:%s\n",g.q,g.K,g.L,g.n,ok?"YES":"NO!!");
        if(!ok) fails++;
        free(cw); free(hg); free(hb); free(g.G); free(g.Gneg); free(g.H);
    }
    printf(fails? "SELFCHECK FAILED\n" : "SELFCHECK PASSED\n");
    return fails?1:0;
}

int main(int argc,char**argv){
    if(argc>=2 && !strcmp(argv[1],"enum"))      return cmd_enum(argc,argv);
    if(argc>=2 && !strcmp(argv[1],"selfcheck")) return cmd_selfcheck();
    fprintf(stderr,"dwq_coset enum <in.cos> <t> <block_start> <block_count> <out.cpart> <nthreads>\n"
                   "dwq_coset selfcheck\n");
    return 2;
}
