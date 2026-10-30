from fractions import Fraction as Fr
import math
K=2**120
def up(x):  # round up to grid 1/K
    return Fr(math.ceil(x*K),K)
def start(d,m):
    n=2*d; e=m-1
    return [up((1-e)*Fr(n-l,n) + e*Fr((n-l)*(n-l-1),n*(n-1))) for l in range(n+1)]
def F(t):
    n=len(t)-1; out=[]
    for l in range(n-1):
        psi = Fr(1) if l==0 else t[l+1]/(1-t[l]+t[l+1])
        out.append(up((1-psi)**2*t[l+2]+2*psi*(1-psi)*t[l+1]+psi**2*t[l]))
    return out
def ok(d,m):
    t=start(d,m)
    for j in range(d,2,-1): t=F(t)
    return all(t[l] < Fr(1,2**l) for l in range(1,5)), [float(x) for x in t[1:5]]
table=[(3,Fr(1834,1000)),(4,Fr(1712,1000)),(5,Fr(1621,1000)),(6,Fr(1553,1000)),(7,Fr(1512,1000)),(8,Fr(1477,1000)),(9,Fr(1447,1000)),(10,Fr(1422,1000)),(15,Fr(135,100)),(20,Fr(1301,1000)),(30,Fr(1241,1000)),(50,Fr(119,100)),(100,Fr(113,100))]
for d,m in table:
    r,t=ok(d,m); print(d,float(m),r,["%.4g"%x for x in t],flush=True)
