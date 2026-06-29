/* bench_hotcore.c -- throughput benchmark for the qutrit weight-enumerator hot core.
 *
 * Compares two single-thread kernels over a realistic [n,K] F_3 generator:
 *   (S) support-scan incremental: Gray step touches only support(changed row),
 *       running weight updated by delta.  Cost ~ row weight per step.
 *   (P) AVX2 byte-packed full-add: each step adds +-g_j over the WHOLE codeword
 *       (mod 3) and recomputes weight by popcount of nonzero-byte mask.
 *       Cost ~ n/32 per step, INDEPENDENT of row weight.
 *
 * Both use the loopless reflected ternary Gray code (Knuth TAOCP 7.2.1.1 Alg H):
 * consecutive messages differ in exactly ONE trit by +-1, never wrapping, so the
 * codeword changes by + or - one generator row.  We benchmark raw steps/sec; the
 * actual run histograms 3^K of them.
 *
 * Build: gcc -O3 -march=native -funroll-loops bench_hotcore.c -o bench_hotcore
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>

#ifndef NMAX
#define NMAX 2048
#endif
#ifndef KMAX
#define KMAX 30
#endif

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }

/* xorshift rng for synthetic generator */
static uint64_t rs=88172645463325252ULL;
static inline uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

int main(int argc,char**argv){
    int n = argc>1?atoi(argv[1]):1939;
    int K = argc>2?atoi(argv[2]):26;
    long  steps = argc>3?atol(argv[3]):200000000L; /* benchmark step count */
    int   density = argc>4?atoi(argv[4]):0; /* 0=dense random GF3, else target row weight */
    if(n>NMAX||K>KMAX){fprintf(stderr,"n/K too big\n");return 1;}

    int np = (n+31)&~31;                 /* pad n up to multiple of 32 for AVX2 */
    /* generator rows as bytes, value in {0,1,2}; G[j] length np (padding=0) */
    uint8_t (*G)[NMAX] = aligned_alloc(32, (size_t)K*NMAX);
    uint8_t (*Gneg)[NMAX] = aligned_alloc(32, (size_t)K*NMAX); /* 2*g mod3 = -g */
    memset(G,0,(size_t)K*NMAX); memset(Gneg,0,(size_t)K*NMAX);
    int rowwt[KMAX];
    for(int j=0;j<K;j++){
        int wt=0;
        if(density<=0){
            for(int p=0;p<n;p++){ int v=xr()%3; G[j][p]=v; if(v)wt++; }
        } else {
            /* sparse: place ~density nonzeros */
            for(int p=0;p<n;p++) G[j][p]=0;
            int placed=0; while(placed<density){ int p=xr()%n; if(!G[j][p]){G[j][p]=1+xr()%2; placed++;} }
            wt=density;
        }
        for(int p=0;p<np;p++) Gneg[j][p]=(uint8_t)((3-G[j][p])%3);
        rowwt[j]=wt;
    }
    /* order rows by weight ascending so low Gray digits (most frequent) get sparsest rows */
    int order[KMAX]; for(int j=0;j<K;j++)order[j]=j;
    for(int a=0;a<K;a++)for(int b=a+1;b<K;b++) if(rowwt[order[b]]<rowwt[order[a]]){int t=order[a];order[a]=order[b];order[b]=t;}
    double avgwt=0; for(int j=0;j<K;j++)avgwt+=rowwt[j]; avgwt/=K;

    /* support lists for kernel S, indexed by gray-digit position d (0=fastest) */
    int  *supp_pos[KMAX]; uint8_t *supp_val[KMAX]; int supp_len[KMAX];
    for(int d=0;d<K;d++){ int j=order[d];
        supp_pos[d]=malloc(sizeof(int)*n); supp_val[d]=malloc(n); int c=0;
        for(int p=0;p<n;p++) if(G[j][p]){ supp_pos[d][c]=p; supp_val[d][c]=G[j][p]; c++; }
        supp_len[d]=c;
    }

    long *hist = calloc(n+1,sizeof(long));

    /* ---- Knuth loopless ternary Gray code state ---- */
    int a[KMAX+1], f[KMAX+1], dir[KMAX+1];
    /* ================= Kernel S: support-scan ================= */
    {
        uint8_t *c = aligned_alloc(32, np); memset(c,0,np);
        for(int j=0;j<=K;j++){a[j]=0;f[j]=j;dir[j]=1;}
        long wt=0; long cnt=0; double t0=now();
        hist[0]++; cnt++;
        while(cnt<steps){
            int d=f[0]; f[0]=0; if(d==K) break;
            int s=dir[d]; a[d]+=s;
            int oldp=(a[d]-s); (void)oldp;
            if(a[d]==0||a[d]==2){ dir[d]=-dir[d]; f[d]=f[d+1]; f[d+1]=d+1; }
            /* apply +-row(d): change = s * g  (mod 3); update support */
            int *sp=supp_pos[d]; uint8_t *sv=supp_val[d]; int L=supp_len[d];
            if(s>0){ for(int i=0;i<L;i++){int p=sp[i]; int old=c[p]; int nv=old+sv[i]; if(nv>=3)nv-=3;
                       wt += (nv!=0)-(old!=0); c[p]=(uint8_t)nv; } }
            else   { for(int i=0;i<L;i++){int p=sp[i]; int old=c[p]; int nv=old+ (3-sv[i]); if(nv>=3)nv-=3;
                       wt += (nv!=0)-(old!=0); c[p]=(uint8_t)nv; } }
            hist[wt]++; cnt++;
        }
        double dt=now()-t0;
        fprintf(stderr,"[S support-scan]  n=%d K=%d avg_rowwt=%.0f  steps=%ld  %.2fs  %.3e steps/s/core\n",
                n,K,avgwt,cnt,dt,cnt/dt);
        free(c);
    }

    /* checksum so the compiler can't elide */
    long chk=0; for(int w=0;w<=n;w++) chk+=hist[w]; memset(hist,0,(n+1)*sizeof(long));

    /* ================= Kernel P: AVX2 byte-packed full add ================= */
    {
        uint8_t *c = aligned_alloc(32, np); memset(c,0,np);
        for(int j=0;j<=K;j++){a[j]=0;f[j]=j;dir[j]=1;}
        long cnt=0; double t0=now();
        __m256i three=_mm256_set1_epi8(3), two=_mm256_set1_epi8(2), zero=_mm256_setzero_si256();
        int nv=np/32;
        hist[0]++; cnt++;
        while(cnt<steps){
            int d=f[0]; f[0]=0; if(d==K) break;
            int s=dir[d]; a[d]+=s;
            if(a[d]==0||a[d]==2){ dir[d]=-dir[d]; f[d]=f[d+1]; f[d+1]=d+1; }
            int j=order[d];
            uint8_t *g = (s>0)? G[j] : Gneg[j];
            long wt=0;
            for(int b=0;b<nv;b++){
                __m256i cv=_mm256_load_si256((__m256i*)(c+32*b));
                __m256i gv=_mm256_load_si256((__m256i*)(g+32*b));
                __m256i sv=_mm256_add_epi8(cv,gv);
                /* sv in 0..4 ; subtract 3 where sv>2  (cmpgt is signed, ok since <128) */
                __m256i ge=_mm256_cmpgt_epi8(sv,two);
                sv=_mm256_sub_epi8(sv,_mm256_and_si256(ge,three));
                _mm256_store_si256((__m256i*)(c+32*b),sv);
                /* weight: bytes != 0 */
                __m256i isz=_mm256_cmpeq_epi8(sv,zero);
                unsigned m=(unsigned)_mm256_movemask_epi8(isz);
                wt += 32 - __builtin_popcount(m);
            }
            hist[wt]++; cnt++;
        }
        double dt=now()-t0;
        fprintf(stderr,"[P avx2-packed]   n=%d K=%d avg_rowwt=%.0f  steps=%ld  %.2fs  %.3e steps/s/core\n",
                n,K,avgwt,cnt,dt,cnt/dt);
        free(c);
    }
    long chk2=0; for(int w=0;w<=n;w++) chk2+=hist[w];
    fprintf(stderr,"chk %ld %ld\n",chk,chk2);
    return 0;
}
