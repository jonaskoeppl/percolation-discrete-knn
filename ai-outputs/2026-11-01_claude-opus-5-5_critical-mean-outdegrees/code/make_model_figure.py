"""TikZ figure of the three graphs for one random configuration on a 6x6 patch of Z^2 with beta=2.5:
each vertex keeps 2 uniformly chosen directions and, with probability 1/2, a third one."""
import random
random.seed(11)
L=6; beta=2.5
dirs=[(1,0),(-1,0),(0,1),(0,-1)]
keep={}
for x in range(L):
    for y in range(L):
        order=random.sample(range(4),4); u=random.random()
        keep[(x,y)]={dirs[order[i]] for i in range(4) if beta>=i+u}
def inside(p): return 0<=p[0]<L and 0<=p[1]<L
def panel(kind,shift):
    out=[]
    S=0.66
    def P(p): return f'({shift+p[0]*S:.3f},{p[1]*S:.3f})'
    for p,K in keep.items():
        for d in dirs:
            q=(p[0]+d[0],p[1]+d[1])
            if not inside(q): continue
            if kind=='D':
                if d in K:
                    m=(p[0]+0.44*d[0],p[1]+0.44*d[1])
                    out.append(rf'\draw[serD,line width=0.9pt,-{{Stealth[length=3.4pt]}}] {P(p)}--{P(m)};')
            elif p<q:
                a=d in K; b=(-d[0],-d[1]) in keep[q]
                if kind=='U' and (a or b): out.append(rf'\draw[serU,line width=1.4pt] {P(p)}--{P(q)};')
                if kind=='B' and (a and b): out.append(rf'\draw[serB,line width=1.4pt] {P(p)}--{P(q)};')
    for p in keep: out.append(rf'\fill[black!70] {P(p)} circle (1.3pt);')
    return out
lines=[r'\begin{figure}[ht]',r'\centering',r'\begin{tikzpicture}']
titles={'D':'directed (D)','U':'undirected (U)','B':'bidirectional (B)'}
for i,k in enumerate('DUB'):
    sh=i*4.6
    lines+=panel(k,sh)
    lines.append(rf'\node[font=\small] at ({sh+2.5*0.66:.3f},{-0.55:.3f}) {{{titles[k]}}};')
lines+=[r'\end{tikzpicture}',
 r'\caption{One configuration of the keep sets on a $6\times6$ patch of $\Z^2$ with mean out-degree $\beta=2.5$: every vertex keeps two uniformly chosen directions and, with probability $1/2$, a third one. Left: the directed graph, with an arrow along half of the edge for every kept direction, so that an edge carries two arrows if and only if it is open in the bidirectional graph. Middle: the undirected graph, in which an edge is open if at least one endpoint keeps it. Right: the bidirectional graph, in which an edge is open if both endpoints keep it.}\label{fig:model}',
 r'\end{figure}']
open(__import__('os').path.join(__import__('os').path.dirname(__file__),'..','figure_model.tex'),'w').write('\n'.join(lines)+'\n')
print('ok')
