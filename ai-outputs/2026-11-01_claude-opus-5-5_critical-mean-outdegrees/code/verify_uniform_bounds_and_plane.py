"""Rigorous (outward-rounded rational) verification of the computer-assisted statements.
Avoidance vectors are rounded UP after each application of F (valid by monotonicity of F)."""
from fractions import Fraction as Fr
import math
K=2**100
def up(x): return Fr(math.ceil(x*K),K)
def down(x): return Fr(math.floor(x*K),K)
def tau(d,b,rnd=up):
    n=2*d; e=b-1
    return [rnd((1-e)*Fr(n-l,n)+e*Fr((n-l)*(n-l-1),n*(n-1))) for l in range(n+1)]
def F(t):
    n=len(t)-1; out=[Fr(1)]
    for l in range(1,n-1):
        psi=t[l+1]/(1-t[l]+t[l+1])
        out.append(up((1-psi)**2*t[l+2]+2*psi*(1-psi)*t[l+1]+psi**2*t[l]))
    return out
def C(t): t=list(t); t[-1]=t[-2]=t[-3]; return t
def planar_ok(t): return all(t[l] < Fr(1,2**l) for l in range(1,5))
# (a) explicit bounds
bars={3:'1.834',4:'1.711',5:'1.621',6:'1.553',7:'1.512',8:'1.477',9:'1.447',10:'1.422',11:'1.405',12:'1.389',13:'1.375',14:'1.362',15:'1.35',16:'1.339',17:'1.328',18:'1.319',19:'1.31',20:'1.301',21:'1.293',22:'1.286',23:'1.279',24:'1.272',25:'1.266',26:'1.26',27:'1.255',28:'1.249',29:'1.244',30:'1.241'}
for d,bs in bars.items():
    b=Fr(bs); t=tau(d,b)
    for _ in range(d-2): t=F(t)
    print("(a)",d,bs,planar_ok(t),flush=True)
# (b) envelope at D=31, beta=1.457 (closure)
b=Fr('1.457'); D=31; t=C(tau(D,b))
for _ in range(D-2): t=C(F(t))
print("(b) omega(31) at 1.457 ok:",planar_ok(t),[float(x) for x in t[1:5]],flush=True)
# (c) window d=2, s=5 on [L2, 1.512]
L2=Fr('1.4929'); hi=Fr('1.512'); steps=40; okall=True
for i in range(steps):
    b1=L2+(hi-L2)*i/steps; b2=L2+(hi-L2)*(i+1)/steps
    t=tau(7,b1)
    for _ in range(5): t=F(t)
    r=tau(2,b2,rnd=down)
    okall &= all(t[l]<=r[l] for l in range(1,4))
print("(c) window d=2,s=5 on [1.4929,1.512]:",okall)
