"""Undirected graph: scan beta upwards from 1 with step h(d) until chi>2600 or clusters are capped;
then fit chi^{-1/gamma} linearly on the points with chi>=40 (gamma as in cluster_size_fit_bidirectional.py). Also reruns B for d=9."""
import subprocess,json
from concurrent.futures import ThreadPoolExecutor
from cluster_size_fit_bidirectional import chi,fit,GAM
def scan(case):
    M,d,b,h=case; g=GAM.get(d,1.0); pts=[]
    while True:
        c,_,cap=chi(d,b,M,3000,int(b*1e5)+d)
        if cap>0.001: break
        pts.append((b,c))
        if c>2600: break
        b+=h
    use=[p for p in pts if p[1]>=40][-6:]
    if len(use)<4: use=pts[-4:]
    est=fit(use,g); alts=[fit(use[1:],g),fit(use[:-1],g)]
    if d<6: alts+=[fit(use,g-0.05),fit(use,g+0.05)]
    r=(M,d,est,max(abs(a-est) for a in alts),pts)
    print(M,d,f'{est:.4f} +- {r[3]:.4f}',' '.join(f'({b:.4f},{c:.0f})' for b,c in pts),flush=True)
    return r
cases=[('U',3,1.0,0.006),('U',4,1.0,0.0025),('U',5,1.0,0.0015),('U',6,1.0,0.0012),('U',7,1.0,0.001),('U',8,1.0,0.0006),('U',9,0.997,0.0005),('U',10,0.997,0.0004),('B',9,4.70,0.012)]
with ThreadPoolExecutor(2) as ex: res=list(ex.map(scan,cases))
json.dump(res,open('../data/chi_scan.json','w'))
