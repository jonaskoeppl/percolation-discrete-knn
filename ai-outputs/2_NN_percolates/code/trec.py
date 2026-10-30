import sys
def F(t):
    n=len(t)-1   # t_0..t_n on n directions; output on n-2 directions
    out=[]
    for l in range(n-1):
        psi = 1.0 if l==0 else t[l+1]/(1-t[l]+t[l+1])
        out.append((1-psi)**2*t[l+2]+2*psi*(1-psi)*t[l+1]+psi**2*t[l])
    return out
def start(d,k=2):
    from math import comb
    return [comb(2*d-l,k)/comb(2*d,k) for l in range(2*d+1)]
def run(d):
    t=start(d)
    for j in range(d,2,-1): t=F(t)
    p=min(1-t[l]**(1/l) for l in range(1,5))
    return p,t[1:5]
if __name__=="__main__":
    for d in [3,4,5,6,7,8,10,15,20,30,50,100,200,500,1000,2000]:
        p,t=run(d); print(d,"p_max=%.5f"%p,["%.3e"%x for x in t])
