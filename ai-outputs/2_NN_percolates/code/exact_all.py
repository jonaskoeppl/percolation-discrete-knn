from fractions import Fraction as Fr
def start(d,m=Fr(2)):
    n=2*d; e=m-1
    return [(1-e)*Fr(n-l,n) + e*Fr((n-l)*(n-l-1),n*(n-1)) for l in range(n+1)]
def C(t): t=list(t); t[-1]=t[-2]=t[-3]; return t
def F(t):
    n=len(t)-1; out=[]
    for l in range(n-1):
        psi = Fr(1) if l==0 else t[l+1]/(1-t[l]+t[l+1])
        out.append((1-psi)**2*t[l+2]+2*psi*(1-psi)*t[l+1]+psi**2*t[l])
    return out
def true(d,m=Fr(2)):
    t=start(d,m)
    for j in range(d,2,-1): t=F(t)
    return t
def W(d,m=Fr(2)):
    t=C(start(d,m))
    for j in range(d,2,-1): t=C(F(t))
    return t
p=Fr(13,25)  # 0.52
for d in (3,4,5):
    t=true(d); print("true d=%d"%d,[float(x) for x in t[1:]], all(t[l]<=(1-p)**l for l in range(1,5)), [len(str(x)) for x in t[1:]])
t=W(6); print("W(6)",[float(x) for x in t[1:]], all(t[l]<=(1-p)**l for l in range(1,5)))
print("(1-p)^l:",[float((1-p)**l) for l in range(1,5)])
print("d=3 exact:",[str(x) for x in true(3)])
