import sympy as sp
a,n=sp.symbols('a n',positive=True)
N=n*(n-1); M=(n-2)*(n-3)
t0=a*(a-1)/N; t1=(a-1)*(a-2)/N; t2=(a-2)*(a-3)/N
psi=t1/(1-t0+t1)
L=(1-psi)**2*t2+2*psi*(1-psi)*t1+psi**2*t0
D=sp.together((a-2)*(a-3)/M - L)
num,den=sp.fraction(sp.factor(D))
print("den:",sp.factor(den)); print("num:",sp.factor(num))
num_e=sp.expand(num); print(sp.Poly(num_e,a).all_coeffs())
