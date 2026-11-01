"""Rigorous bounds on the critical mean out-degrees beta_c^B, beta_c^D, beta_c^U.
phi(beta)=k(k-1+2eps), beta=k+eps: m(m-1) times the probability that two given directions are kept.
B: percolation if phi > m(m-1) pc(d); no percolation if c(d) phi < m(m-1).
D: percolation if beta > m pc(d) [BJKLRT25, Cor. 3.7]; no percolation if c(d) beta < m; upper bounds of the
   beta_c^D <= 2 (CHJK26, LXZ26) and the column-projection bounds BAR are also used.
U: percolation if 1-(1-q)^2 > pc(d), q^2=phi/(m(m-1)), or beta > beta_c^D; no percolation if c(d) lambda(beta) < 1
   (transfer matrix) or beta <= 1.
All upper bounds on pc and c are rigorous published values (see the paper); bisection is done in exact
rationals on monotone functions and rounded outward."""
from fractions import Fraction as Fr
import math
def gps(d,it=70):
    e=[(d+i)//3 for i in range(3)]
    def F(p):
        pr=Fr(1)
        for x in e: pr*=1-(1-p)**x
        return pr-(2-sum((1-p)**x for x in e))
    lo,hi=Fr(1,1000),Fr(999,1000)
    for _ in range(it):
        mid=(lo+hi)/2
        lo,hi=(mid,hi) if F(mid)>0 else (lo,mid)
    return hi
FVH={11:Fr('1.01306'),12:Fr('1.00857'),13:Fr('1.006244'),14:Fr('1.00484'),15:Fr('1.0039'),20:Fr('1.001777')}
def pcbar(d):
    if d==2: return Fr(1,2)
    b=None
    for dd in range(3,d+1):
        x=gps(dd)
        if dd in FVH: x=min(x,FVH[dd]/(2*dd-1))
        b=x if b is None else min(b,x)
    return b
CBAR={2:Fr('2.679192495'),3:Fr('4.7387'),4:Fr('6.8040'),5:Fr('8.8602')}
def cbar(d): return CBAR.get(d,Fr(2*d-1))
# certified upper bounds on beta_c^D from the column projection (verify_uniform_bounds_and_plane.py); 1.457 for all d>=31
BAR={3:'1.834',4:'1.711',5:'1.621',6:'1.553',7:'1.512',8:'1.477',9:'1.447',10:'1.422',11:'1.405',12:'1.389',13:'1.375',14:'1.362',15:'1.35',16:'1.339',17:'1.328',18:'1.319',19:'1.31',20:'1.301',21:'1.293',22:'1.286',23:'1.279',24:'1.272',25:'1.266',26:'1.26',27:'1.255',28:'1.249',29:'1.244',30:'1.241'}
def phi(b):
    k=math.floor(b); e=b-k; return k*(k-1+2*e)
def phiinv(y):  # smallest beta with phi(beta)=y, y>0
    k=1
    while (k+1)*k<=y: k+=1
    return k+(y-k*(k-1))/(2*k)
def bis(f,lo,hi,it=60):  # f monotone increasing, find root
    lo,hi=Fr(lo),Fr(hi)
    for _ in range(it):
        mid=(lo+hi)/2
        if f(mid)>0: hi=mid
        else: lo=mid
    return lo,hi
def lam(b,m):
    # increasing majorant of the Perron root of M_{d,beta}: lambda <= r - s/2 + sqrt(s), r=beta/m, s=phi/(m(m-1))
    b=float(b); s=phi(b)/(m*(m-1)); r=b/m
    return r-s/2+math.sqrt(s)
rows=[]
for d in [2,3,4,5,6,7,8,9,10,11,12,15,20]:
    m=2*d; pc=pcbar(d); c=cbar(d)
    Bl=phiinv(Fr(m*(m-1))/c); Bu=phiinv(Fr(m*(m-1))*pc)
    Dl=Fr(m)/c; Du=min([Fr(m)*pc,Fr(2)]+([Fr(BAR[d])] if d in BAR else []))
    # U lower: c*lambda(beta)=1 (lambda increasing); float with outward margin 1e-9
    lo,hi=bis(lambda b: float(c)*lam(b,m)-1,1,m-1)
    Ul=max(1.0,float(lo)-1e-9)
    # U upper via comparison: smallest beta with 1-(1-q)^2 > pc
    qs=1-math.sqrt(1-float(pc)); Uc=phiinv(qs*qs*m*(m-1))+1e-9
    Uu=min(float(Du),Uc)
    rows.append((d,float(Bl),float(Bu),float(Dl),float(Du),Ul,Uu,Uc))
    print(f"d={d:2d}  B:[{float(Bl):.4f},{float(Bu):.4f}]  D:[{float(Dl):.4f},{float(Du):.4f}]  U:[{Ul:.4f},{Uu:.4f}] (U-comp {Uc:.4f})  sqrt(2d)+1/2={math.sqrt(2*d)+0.5:.4f}")
