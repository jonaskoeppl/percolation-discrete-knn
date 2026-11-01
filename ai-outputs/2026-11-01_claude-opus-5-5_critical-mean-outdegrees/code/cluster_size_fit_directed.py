"""Directed graph: mean out-cluster size chi(beta) on Z^d. Scan upwards from beta=1 with step h,
reducing the step by a factor 3 whenever chi exceeds 150, 400, 1000; stop at chi>2500 or capping.
Fit chi^{-1/gamma} linearly on the last points with chi>=40 (gamma as for Bernoulli percolation)."""
import json
from concurrent.futures import ThreadPoolExecutor
from cluster_size_fit_bidirectional import chi,fit,GAM
def scan(d):
    g=GAM.get(d,1.0); b=1.0; h={3:0.04,4:0.02}.get(d,0.01); pts=[]; thr=[150,400,1000]
    while True:
        c,_,cap=chi(d,b,'D',3000,int(b*1e5)+7*d)
        if cap>0.001:
            if h<1e-4: break
            b-=h; h/=3; b+=h; continue
        pts.append((b,c))
        if c>2500: break
        while thr and c>thr[0]: thr.pop(0); h/=3
        b+=h
    use=[p for p in pts if p[1]>=40][-6:]
    if len(use)<4: use=pts[-4:]
    est=fit(use,g); alts=[fit(use[1:],g),fit(use[:-1],g)]
    if d<6: alts+=[fit(use,g-0.05),fit(use,g+0.05)]
    r=('D',d,est,max(abs(a-est) for a in alts),pts)
    print('D',d,f'{est:.4f} +- {r[3]:.4f}',' '.join(f'({x:.4f},{y:.0f})' for x,y in pts),flush=True)
    return r
with ThreadPoolExecutor(2) as ex: res=list(ex.map(scan,range(3,11)))
json.dump(res,open('../data/chi_D.json','w'))
