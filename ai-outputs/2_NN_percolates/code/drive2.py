import subprocess, sys
def run(d,b,N,seed):
    out=subprocess.run(["./simr64",str(d),"%.4f"%b,str(N),str(seed)],capture_output=True,text=True).stdout.split()
    a,bb,c=int(out[3]),int(out[4]),int(out[5]); return bb/max(a,1), c/max(bb,1)
d=int(sys.argv[1]); lo=float(sys.argv[2]); hi=float(sys.argv[3]); N=int(sys.argv[4])
pts=[]; b=lo; s=1
while b<=hi+1e-9:
    r1,r2=run(d,b,N,100*d+s); pts.append((round(b,3),round(r2-r1,4))); s+=1; b+=0.01
cr=None
for (b0,f0),(b1,f1) in zip(pts,pts[1:]):
    if f0<0<=f1: cr=b0+(b1-b0)*(-f0)/(f1-f0); break
print(d, "crossing", None if cr is None else round(cr,3), pts, flush=True)
