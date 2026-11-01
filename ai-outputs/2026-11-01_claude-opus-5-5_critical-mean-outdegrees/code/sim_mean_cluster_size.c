/* Mean cluster size chi(beta) of the origin in the B, U or D graph on Z^d (no boundary).
   Keep sets are generated lazily (uniform order of the 2d directions plus U_x; the direction in
   position i is kept iff beta >= i+U_x) and stored in an open-addressing hash table.
   For D the out-cluster is explored. Clusters reaching CAP vertices are counted as CAP.
   usage: chi d beta samples seed B|U|D CAP  -> d beta model mean stderr fraction_capped */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){const uint64_t r=rotl(s[1]*5,7)*9;const uint64_t t=s[1]<<17;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45);return r;}
#define MAXD 12
typedef struct { int16_t c[MAXD]; uint32_t keep; uint32_t stamp; int used; } Ent;
static Ent *tab; static uint64_t HS; static int d; static double beta; static uint32_t stamp=1;
static uint64_t hashc(const int16_t*c){uint64_t h=1469598103934665603ULL; for(int i=0;i<d;i++){h^=(uint16_t)c[i]; h*=1099511628211ULL;} return h^(h>>29);}
static Ent* get(const int16_t*c){ uint64_t h=hashc(c)&(HS-1);
  while(1){ Ent*e=&tab[h];
    if(e->stamp!=stamp){ e->stamp=stamp; memcpy(e->c,c,sizeof(int16_t)*d); e->used=0;
      int m=2*d,perm[2*MAXD]; for(int j=0;j<m;j++)perm[j]=j;
      for(int j=m-1;j>0;j--){int r=nxt()%(j+1);int t=perm[j];perm[j]=perm[r];perm[r]=t;}
      double u=(nxt()>>11)*0x1.0p-53; uint32_t k=0; for(int p=0;p<m;p++) if(beta>=p+u) k|=1u<<perm[p];
      e->keep=k; return e; }
    if(memcmp(e->c,c,sizeof(int16_t)*d)==0) return e; h=(h+1)&(HS-1);} }
int main(int argc,char**argv){
  d=atoi(argv[1]); beta=atof(argv[2]); int S=atoi(argv[3]); uint64_t seed=strtoull(argv[4],0,10); char M=argv[5][0]; long CAP=atol(argv[6]);
  for(int i=0;i<4;i++){seed+=0x9E3779B97F4A7C15ULL;uint64_t w=seed;w=(w^(w>>30))*0xBF58476D1CE4E5B9ULL;w=(w^(w>>27))*0x94D049BB133111EBULL;s[i]=w^(w>>31);}
  HS=1; while(HS<(uint64_t)(CAP*(2*d+1)*2)) HS<<=1; tab=calloc(HS,sizeof(Ent));
  int16_t (*queue)[MAXD]=malloc(sizeof(int16_t)*MAXD*(CAP+1));
  double sum=0,sum2=0; long capped=0;
  for(int smp=0;smp<S;smp++){ stamp++; long qh=0,qt=0; int16_t o[MAXD]={0};
    Ent*e=get(o); e->used=1; memcpy(queue[qt++],o,sizeof(int16_t)*d); long size=1; int hit=0;
    while(qh<qt && !hit){ int16_t x[MAXD]; memcpy(x,queue[qh++],sizeof(int16_t)*d); uint32_t kx=get(x)->keep;
      for(int j=0;j<2*d;j++){ int i=j/2, sg=(j%2)?-1:1; int16_t y[MAXD]; memcpy(y,x,sizeof(int16_t)*d); y[i]+=sg;
        int jx=(kx>>j)&1, open;
        if(M=='D') open=jx; else { Ent*ey=get(y); int jy=(ey->keep>>(j^1))&1; open=(M=='B')?(jx&&jy):(jx||jy); }
        if(!open) continue;
        Ent*ey=get(y); if(ey->used) continue; ey->used=1; size++;
        if(size>=CAP){hit=1;break;}
        memcpy(queue[qt++],y,sizeof(int16_t)*d); } }
    if(hit) capped++; sum+=size; sum2+=(double)size*size;
    if(smp>=50 && capped*20>smp){ printf("%d %.5f %c %.4f %.4f %.5f\n",d,beta,M,sum/(smp+1),0.0,(double)capped/(smp+1)); return 0; } }
  double mean=sum/S, var=sum2/S-mean*mean;
  printf("%d %.5f %c %.4f %.4f %.5f\n",d,beta,M,mean,sqrt(var/S),(double)capped/S); return 0; }
