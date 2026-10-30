import random
def run(d,T=100000,cap=20000,seed=3):
    random.seed(seed); m=d-2
    fib=[(i,s) for i in range(m) for s in (1,-1)]
    D=['a','-a','b','-b']+fib
    ip=['a','-a','b','-b']
    hits=[0]*5; capped=0
    for _ in range(T):
        start=(0,)*m; seen={start}; stack=[start]; N=set(); 
        while stack:
            z=stack.pop(); P=random.sample(D,2)
            for x in P:
                if isinstance(x,str): N.add(x)
                else:
                    i,s=x; y=list(z); y[i]+=s; y=tuple(y)
                    if y not in seen: seen.add(y); stack.append(y)
            if len(seen)>cap: capped+=1; N=set(ip); break
        if 'a' in N: hits[1]+=1
        if N&{'a','b'}: hits[2]+=1
        if N&{'a','b','-a'}: hits[3]+=1
        if N: hits[4]+=1
    h=[x/T for x in hits]
    print(d,'capped',capped/T,'hit',[round(x,4) for x in h[1:]],'p allowed',[round(1-(1-h[l])**(1/l),4) for l in range(1,5)])
for d in (3,4,5,6): run(d)
