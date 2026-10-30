// P(origin's out-cluster reaches l1-distance R) for the 2dp-neighbour graph on Z^d, mean degree m in [1,2]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint64_t s0=88172645463325252ULL;
static inline uint64_t rng(void){s0^=s0<<13;s0^=s0>>7;s0^=s0<<17;return s0;}
static inline double U(void){return (rng()>>11)*(1.0/9007199254740992.0);}
#define HB 24
static uint64_t *keys; static uint32_t *stamp; static size_t HS=1ULL<<HB;
static int d,R; static uint32_t it;
static inline uint64_t enc(int*x){uint64_t k=0;for(int i=0;i<d;i++)k=k*256+(uint64_t)(x[i]+128);return k+1;}
static int insert(uint64_t k){ uint64_t h=(k*0x9E3779B97F4A7C15ULL)>>(64-HB);
  while(stamp[h]==it){ if(keys[h]==k) return 0; h=(h+1)&(HS-1);} stamp[h]=it; keys[h]=k; return 1;}
int main(int argc,char**argv){
  d=atoi(argv[1]); double m=atof(argv[2]); R=atoi(argv[3]); int N=atoi(argv[4]); s0+=atoi(argv[5]);
  keys=malloc(HS*8); stamp=calloc(HS,4);
  int cap=1<<22; int *q=malloc((size_t)cap*d*sizeof(int)); int hit=0;
  for(it=1;it<=(uint32_t)N;it++){
    int head=0,tail=0; memset(q,0,d*sizeof(int)); insert(enc(q)); tail=1; int reached=0;
    while(head<tail){
      int *x=q+(size_t)head*d; head++;
      int dist=0; for(int i=0;i<d;i++) dist+=abs(x[i]);
      if(dist>=R){reached=1;break;}
      int kk=1+(U()<m-1.0); int idx[64]; for(int j=0;j<2*d;j++) idx[j]=j;
      for(int j=0;j<kk;j++){ int r=j+(int)(U()*(2*d-j)); int t=idx[j];idx[j]=idx[r];idx[r]=t;
        int y[16]; memcpy(y,x,d*sizeof(int)); int c=idx[j]/2; y[c]+= (idx[j]%2)?-1:1;
        if(insert(enc(y))){ if(tail>=cap){reached=1;head=tail;break;} memcpy(q+(size_t)tail*d,y,d*sizeof(int)); tail++; } }
      if(tail>=cap) {reached=1;break;}
    }
    hit+=reached;
  }
  printf("d=%d m=%.3f R=%d P=%.4f\n",d,m,R,(double)hit/N); return 0;}
