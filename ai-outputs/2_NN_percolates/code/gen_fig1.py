import random, sys
R=4  # box [-R,R]^2
dirs=[(1,0),(-1,0),(0,1),(0,-1)]
def sample(seed):
    rnd=random.Random(seed)
    N={(x,y):rnd.sample(dirs,2) for x in range(-R,R+1) for y in range(-R,R+1)}
    # forward cluster of origin inside box (and record whether it exits)
    C={(0,0)}; st=[(0,0)]; exits=False
    while st:
        v=st.pop()
        for d in N[v]:
            w=(v[0]+d[0],v[1]+d[1])
            if max(abs(w[0]),abs(w[1]))>R: exits=True; continue
            if w not in C: C.add(w); st.append(w)
    return N,C,exits
best=None
for seed in range(2000):
    N,C,ex=sample(seed)
    if ex and 18<=len(C)<=30:
        best=seed; break
N,C,ex=sample(best)
print("% seed",best,"cluster size",len(C),"exits",ex, file=sys.stderr)
out=[]
off=0.09
for (x,y),ds in N.items():
    for (dx,dy) in ds:
        tx,ty=x+dx,y+dy
        inside=max(abs(tx),abs(ty))<=R
        # offset perpendicular if reverse edge exists
        rev = inside and (-dx,-dy) in N[(tx,ty)]
        ox,oy=(-dy*off,dx*off) if rev else (0,0)
        sx,sy=x+ox+0.18*dx, y+oy+0.18*dy
        if inside: ex_,ey_=tx+ox-0.18*dx, ty+oy-0.18*dy
        else: ex_,ey_=x+ox+0.6*dx, y+oy+0.6*dy
        style="clus" if (x,y) in C else "edge"
        out.append(f"\\draw[{style}] ({sx:.2f},{sy:.2f}) -- ({ex_:.2f},{ey_:.2f});")
for x in range(-R,R+1):
    for y in range(-R,R+1):
        st="cv" if (x,y) in C else "v"
        out.append(f"\\node[{st}] at ({x},{y}) {{}};")
out.append("\\node[orig] at (0,0) {};")
open("fig1_body.tex","w").write("\n".join(out)+"\n")
