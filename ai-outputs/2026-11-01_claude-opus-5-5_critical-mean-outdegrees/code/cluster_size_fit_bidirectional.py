"""Estimate beta_c from the mean cluster size: chi^{-1/gamma} is asymptotically linear in beta_c-beta.
gamma: d=3 1.80, d=4 1.44, d=5 1.18, d>=6 1 (mean field; logarithmic corrections in d=6 are ignored).
Adaptive: successive beta values aim at chi ~ 20,40,...,2560; then a linear fit of chi^{-1/gamma}
on the points with chi >= 40 is extrapolated to zero. Error: spread of the fits under gamma +- 0.05
(d<6) and dropping the first/last fitted point."""
import subprocess,sys,math,json
from concurrent.futures import ThreadPoolExecutor
GAM={3:1.80,4:1.44,5:1.18}
NUMPC={3:0.24881,4:0.16013,5:0.11817,6:0.09420,7:0.07868,8:0.06770,9:0.05950,10:0.05309}
def phiinv(y):
    k=1
    while (k+1)*k<=y: k+=1
    return k+(y-k*(k-1))/(2*k)
def chi(d,b,M,S,seed):
    out=subprocess.run(['./sim_mean_cluster_size',str(d),f'{b:.6f}',str(S),str(seed),M,'400000'],capture_output=True,text=True).stdout.split()
    return float(out[3]),float(out[4]),float(out[5])
def fit(pts,g):
    xs=[b for b,c in pts]; ys=[c**(-1/g) for b,c in pts]
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx; return mx-my/a
def run(case):
    M,d=case; g=GAM.get(d,1.0); m=2*d
    b0 = phiinv(m*(m-1)*0.75*NUMPC[d]) if M=='B' else 0.95
    S=3000; pts=[]; seed=d*100+ord(M)
    c,_,_=chi(d,b0,M,S,seed); pts.append((b0,c))
    b1=b0+(0.05 if M=='B' else 0.01); c1,_,_=chi(d,b1,M,S,seed+1); pts.append((b1,c1))
    target=max(c1,20)*2; it=0
    while target<=2600 and it<14:
        (xa,ca),(xb,cb)=pts[-2],pts[-1]
        ya,yb=ca**(-1/g),cb**(-1/g)
        if yb>=ya: bn=xb+(xb-xa)
        else:
            bc=xb+yb*(xb-xa)/(ya-yb); yt=target**(-1/g)
            bn=bc-(bc-xb)*yt/yb
        bn=min(bn,xb+3*(xb-xa))
        cn,_,cap=chi(d,bn,M,S,seed+2+it); it+=1
        if cap>0.001:
            pts_bad=bn; target=10**9; continue
        pts.append((bn,cn)); target=cn*2
    use=[p for p in pts if p[1]>=40]
    if len(use)<3: use=pts[-3:]
    est=fit(use,g)
    alts=[fit(use[1:],g),fit(use[:-1],g)]
    if d<6: alts+= [fit(use,g-0.05),fit(use,g+0.05)]
    err=max(abs(a-est) for a in alts)
    return (M,d,est,err,pts)
if __name__=="__main__":
  cases=[(M,d) for d in range(3,11) for M in "BU"]
  with ThreadPoolExecutor(2) as ex:
      res=list(ex.map(run,cases))
  for M,d,est,err,pts in res: print(M,d,f'{est:.4f} +- {err:.4f}',' '.join(f'({b:.4f},{c:.0f})' for b,c in pts))
  json.dump([(M,d,est,err,pts) for M,d,est,err,pts in res],open('../data/chi_estimates.json','w'))
