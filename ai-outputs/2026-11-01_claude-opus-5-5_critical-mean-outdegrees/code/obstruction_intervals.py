"""Section 6: comparison intervals for Theorem 2.3 on a beta-grid, and difference vectors."""
from fractions import Fraction as Fr
def tau(d,b):
    n=2*d; e=b-1
    return [(1-e)*Fr(n-l,n)+e*Fr((n-l)*(n-l-1),n*(n-1)) for l in range(n+1)]
def F(t):
    n=len(t)-1; out=[Fr(1)]
    for l in range(1,n-1):
        psi=t[l+1]/(1-t[l]+t[l+1]); out.append((1-psi)**2*t[l+2]+2*psi*(1-psi)*t[l+1]+psi**2*t[l])
    return out
def holds(d,s,b):
    t=tau(d+s,b)
    for _ in range(s): t=F(t)
    r=tau(d,b); return all(t[l]<=r[l] for l in range(1,2*d))
def interval(d,s,step=Fr(1,200)):
    g=[Fr(1)+i*step for i in range(int(1/step)+1) if holds(d,s,Fr(1)+i*step)]
    return (float(min(g)),float(max(g))) if g else None
for d,s in [(2,1),(3,1),(4,1),(6,1),(10,1),(2,3),(2,4),(2,5),(2,10)]:
    print("d=%d s=%d interval:"%(d,s),interval(d,s),flush=True)
for d,b in [(2,Fr('1.84')),(3,Fr('1.38'))]:
    l=F(tau(d+1,b)); r=tau(d,b); print("differences d=%d beta=%s:"%(d,b),[round(float(r[i]-l[i]),4) for i in range(1,2*d)])
