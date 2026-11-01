/* Spanning thresholds for the continuous bidirectional (B) and undirected (U)
   lattice neighbour graphs on the box {0..L-1}^d with free boundary.
   Monotone coupling in beta in [0,2d]: vertex x draws a uniform permutation of
   the 2d directions and U_x ~ Unif[0,1); the direction in position i (i=0,..,2d-1)
   is kept iff beta >= i + U_x. For fixed beta = k+eps this keeps the first k
   directions and, with probability eps, the (k+1)st: the law (1-eps)ups_k + eps ups_{k+1}.
   Edge {x,x+e_j} has threshold max (B) or min (U) of the two keep thresholds.
   For each sample we output the smallest beta at which the faces x_1=0 and x_1=L-1
   are connected (Newman-Ziff style: edges added in increasing threshold order).
   usage: span d L samples seed B|U   -> prints one threshold per line */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){const uint64_t r=rotl(s[1]*5,7)*9;const uint64_t t=s[1]<<17;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45);return r;}
static inline double unif(void){return (nxt()>>11)*0x1.0p-53;}
static int *par;
static int find(int x){while(par[x]!=x){par[x]=par[par[x]];x=par[x];}return x;}
int main(int argc,char**argv){
  int d=atoi(argv[1]),L=atoi(argv[2]),S=atoi(argv[3]); uint64_t seed=strtoull(argv[4],0,10); int B=(argv[5][0]=='B');
  uint64_t z=seed*0x9E3779B97F4A7C15ULL+12345; for(int i=0;i<4;i++){z+=0x9E3779B97F4A7C15ULL;uint64_t w=z;w=(w^(w>>30))*0xBF58476D1CE4E5B9ULL;w=(w^(w>>27))*0x94D049BB133111EBULL;s[i]=w^(w>>31);}
  long N=1; for(int i=0;i<d;i++) N*=L; int m=2*d;
  long *stride=malloc(sizeof(long)*d); stride[0]=1; for(int i=1;i<d;i++) stride[i]=stride[i-1]*L;
  float *t=malloc(sizeof(float)*N*m);   /* keep threshold of direction j at x: index x*m+j, j=2i (+e_i), 2i+1 (-e_i) */
  long E=0; for(int i=0;i<d;i++) E+= N/L*(L-1);
  uint32_t *key=malloc(4*E),*key2=malloc(4*E); uint32_t *eid=malloc(4*E),*eid2=malloc(4*E);
  int *ea=malloc(4*E),*eb=malloc(4*E);
  par=malloc(sizeof(int)*(N+2)); int SRC=N,TGT=N+1;
  /* edge list (fixed) */
  long e=0; for(long x=0;x<N;x++){ for(int i=0;i<d;i++){ int c=(x/stride[i])%L; if(c<L-1){ea[e]=x;eb[e]=x+stride[i];e++;} } }
  int perm[64];
  for(int smp=0;smp<S;smp++){
    for(long x=0;x<N;x++){ for(int j=0;j<m;j++) perm[j]=j; for(int j=m-1;j>0;j--){int r=nxt()%(j+1);int tmp=perm[j];perm[j]=perm[r];perm[r]=tmp;}
      double u=unif(); for(int pos=0;pos<m;pos++) t[x*m+perm[pos]]=(float)(pos+u); }
    for(long q=0;q<E;q++){ long a=ea[q],b=eb[q]; long diff=b-a; int i=0; while(stride[i]!=diff) i++;
      float ta=t[a*m+2*i], tb=t[b*m+2*i+1]; float th=B?(ta>tb?ta:tb):(ta<tb?ta:tb);
      uint32_t kb; memcpy(&kb,&th,4); key[q]=kb; eid[q]=q; }
    /* LSD radix sort, 2 passes of 16 bits (keys are non-negative floats) */
    { static long cnt[65537]; for(int pass=0;pass<2;pass++){ int sh=16*pass; memset(cnt,0,sizeof(cnt));
        for(long q=0;q<E;q++) cnt[((key[q]>>sh)&0xFFFF)+1]++;
        for(int r=0;r<65536;r++) cnt[r+1]+=cnt[r];
        for(long q=0;q<E;q++){ long p=cnt[(key[q]>>sh)&0xFFFF]++; key2[p]=key[q]; eid2[p]=eid[q]; }
        uint32_t *tk=key;key=key2;key2=tk; tk=eid;eid=eid2;eid2=tk; } }
    for(long x=0;x<N+2;x++) par[x]=x;
    for(long x=0;x<N;x++){ int c=x%L; if(c==0) par[find(x)]=find(SRC); else if(c==L-1) par[find(x)]=find(TGT);} 
    float ans=-1;
    for(long q=0;q<E;q++){ long id=eid[q]; int ra=find(ea[id]),rb=find(eb[id]); if(ra!=rb){par[ra]=rb; if(find(SRC)==find(TGT)){ memcpy(&ans,&key[q],4); break;} } }
    printf("%.6f\n",ans); fflush(stdout);
  }
  return 0;
}
