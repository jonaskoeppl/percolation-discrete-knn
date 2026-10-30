from trec import F
def start(d,m):
    n=2*d; e=m-1
    return [(1-e)*(n-l)/n + e*(n-l)*(n-l-1)/(n*(n-1)) for l in range(n+1)]
def C(t): t=list(t); t[-1]=t[-2]=t[-3]; return t
def good(t): return all(t[l] < 0.5**l for l in range(1,5))
def W(d,m):
    t=C(start(d,m))
    for j in range(d,2,-1): t=C(F(t))
    return t
def true(d,m):
    t=start(d,m)
    for j in range(d,2,-1): t=F(t)
    return t
def mstar(fun,d):
    lo,hi=1.0,2.0
    if not good(fun(d,2.0)): return None
    for _ in range(40):
        mid=(lo+hi)/2
        if good(fun(d,mid)): hi=mid
        else: lo=mid
    return hi
# check star_m
viol=[]
for m in [1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]:
    for d in range(2,60):
        a=F(start(d+1,m)); b=start(d,m)
        v=[l for l in range(2*d-1) if a[l]>b[l]+1e-15]
        if v: viol.append((m,d,v[:3])); break
print("star_m violations (first per m):",viol)
for d in [3,4,5,6,7,8,10,12,15,20,30,50,100]:
    print(d,"m*_true=%s"%(mstar(true,d)),"m*_W=%s"%(mstar(W,d)))
