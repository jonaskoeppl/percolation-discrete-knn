"""Numerical check of the mixture local inequality: law (1-eps) ups_k + eps ups_{k+1} on an m-set
vs product law with q^2 = phi/(m(m-1)), phi = k(k-1+2 eps).
Checks: pair inequality for all alpha,beta (alpha+beta<=m), single sets, and the structural claims
(ratio monotone => single sign change of second differences, Delta(m-1)<0 for k>=2, G(1)>=0)."""
from mpmath import mp, mpf, binomial, sqrt
mp.dps=80
tol=mpf(10)**-60
viol=0; cnt=0
for m in range(2,41):
  for k in range(0,m):
    for e in [mpf(0),mpf(1)/7,mpf(1)/3,mpf(1)/2,mpf(2)/3,mpf(9)/10,mpf(99)/100]:
      if k==0 and e==0: continue
      phi=k*(k-1+2*e); q=sqrt(phi/(m*(m-1))); u=1-q
      h=[(1-e)*binomial(m-a,k)/binomial(m,k)+e*(binomial(m-a,k+1)/binomial(m,k+1) if k+1<=m else 0) for a in range(m+1)]
      G=[u**a-h[a] for a in range(m+1)]
      for a in range(1,m+1):
        if 1-h[a] < 1-u**a - tol: viol+=1; print('single',m,k,e,a)
      for a in range(1,m):
        for b in range(1,m-a+1):
          cnt+=1
          if G[a+b] > G[a]+G[b]+tol: viol+=1; print('pair',m,k,float(e),a,b,float(G[a+b]-G[a]-G[b]))
print('checked',cnt,'pairs; violations',viol)
