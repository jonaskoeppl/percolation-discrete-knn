"""Estimate critical points from spanning thresholds: for consecutive box sizes L<L', the empirical
distribution functions F_L and F_L' of the spanning threshold cross near beta_c. Bootstrap errors."""
import glob,re,random,statistics as st,sys,json
random.seed(1)
def load(f): return sorted(float(x) for x in open(f) if x.strip() and float(x)>=0)
def F(xs,b):
    import bisect; return bisect.bisect_right(xs,b)/len(xs)
def cross(a,b):
    lo=min(a[len(a)//10],b[len(b)//10]); hi=max(a[9*len(a)//10],b[9*len(b)//10])
    grid=[lo+(hi-lo)*i/400 for i in range(401)]
    diff=[F(a,g)-F(b,g) for g in grid]
    # last sign change from >0 (small box spans earlier) to <=0; use the point where diff crosses 0 nearest the medians
    best=None
    for i in range(400):
        if diff[i]>0 and diff[i+1]<=0 or diff[i]<0 and diff[i+1]>=0:
            x=grid[i]+(grid[i+1]-grid[i])*diff[i]/(diff[i]-diff[i+1]) if diff[i]!=diff[i+1] else grid[i]
            best=x if best is None or abs(x-st.median(b))<abs(best-st.median(b)) else best
    return best
res={}
for M in 'BU':
    for d in range(2,5):
        fs=sorted(glob.glob(f'../data/{M}_d{d}_L*.txt'),key=lambda f:int(re.search(r'L(\d+)',f).group(1)))
        data=[(int(re.search(r'L(\d+)',f).group(1)),load(f)) for f in fs]
        data=[(L,x) for L,x in data if len(x)>=100]
        out=[]
        for (L1,a),(L2,b) in zip(data,data[1:]):
            c=cross(a,b)
            if c is None: continue
            bs=[]
            for _ in range(100):
                aa=sorted(random.choices(a,k=len(a))); bb=sorted(random.choices(b,k=len(b)))
                cc=cross(aa,bb)
                if cc is not None: bs.append(cc)
            err=st.pstdev(bs) if len(bs)>5 else float('nan')
            out.append((L1,L2,len(a),len(b),c,err,st.median(b)))
        res[f'{M}{d}']=out
        for o in out: print(M,d,'L=%d/%d n=%d/%d crossing %.4f +- %.4f  median(L2) %.4f'%o)
json.dump(res,open('../data/estimates.json','w'),indent=1)
