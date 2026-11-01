from itertools import combinations
import math
q=1/math.sqrt(2); D=range(4)
U=[set(D)-{r} for r in D]
def unif(f): return sum(f(N) for N in U)/4
def bern(f):
    s=0
    for m in range(16):
        N={i for i in D if m>>i&1}; s+=f(N)*q**len(N)*(1-q)**(4-len(N))
    return s
res=[]
for a in range(1,5):
    A=set(range(a)); f=lambda N:len(N&A)>0
    res.append((('S',a),unif(f),bern(f)))
for a in range(1,4):
  for b in range(1,5-a):
    A=set(range(a));B=set(range(a,a+b)); f=lambda N:len(N&A)>0 and len(N&B)>0
    res.append(((a,b),unif(f),bern(f)))
for r in res: print(r, r[1]-r[2])
