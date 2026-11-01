"""Lemma 4.1: identity (4.1) and Bernstein certificate.
For n=2d+2, a=n-l, eps=beta-1:
  tau^(d)_l - F(tau^(d+1))_l = 2(n-a) S(a,n,eps) / ( n(n-1)(n-2)(n-3) ((n-1)^2 + eps(n+1-2a))^2 ),
and S>=0 for beta in [1.45,2], n>=8, 4<=a<=n-1."""
import sympy as sp
a,n,e,x,y,s=sp.symbols('a n e x y s')
N=n*(n-1); M=(n-2)*(n-3)
T=lambda b,nn,NN: (1-e)*b/nn + e*b*(b-1)/NN
t0,t1,t2=T(a,n,N),T(a-1,n,N),T(a-2,n,N)
psi=t1/(1-t0+t1)
Fl=(1-psi)**2*t2+2*psi*(1-psi)*t1+psi**2*t0
diff=T(a-2,n-2,M)-Fl
den=n*(n-1)*(n-2)*(n-3)*((n-1)**2+e*(n+1-2*a))**2
S=sp.factor(sp.cancel(diff*den/(2*(n-a))))
S=sp.expand(S); assert S.is_polynomial(a,n,e)  # exact cancellation, so the identity holds by construction
import random
for _ in range(20):
    v={a:random.randint(4,30),e:sp.Rational(random.randint(0,100),100)}; v[n]=v[a]+random.randint(1,20)
    assert diff.subs(v)==(2*(n-a)*S/den).subs(v)
print("S(a,n,eps) =",S)
print("degrees (a,n,eps):",sp.degree(S,a),sp.degree(S,n),sp.degree(S,e))
e0=sp.Rational(9,20); n0=8
Ss=sp.expand(S.subs(e,e0+(1-e0)*s)); deg=sp.degree(Ss,s)
regions=[sp.expand(Ss.subs({n:a+1+y}).subs(a,n0-1+x))]+[sp.expand(Ss.subs({a:A,n:n0+y})) for A in range(4,n0-1)]
neg=0
for R in regions:
    pe=sp.Poly(R,s); c=[pe.coeff_monomial(s**k) for k in range(deg+1)]
    for i in range(deg+1):
        b=sp.expand(sum(sp.binomial(i,k)/sp.binomial(deg,k)*c[k] for k in range(i+1)))
        gens=[g for g in (x,y) if b.has(g)] or [x]
        neg+=sum(1 for co in sp.Poly(b,*gens).coeffs() if co<0)
print("negative Bernstein coefficients:",neg, "(certificate valid iff 0)")
