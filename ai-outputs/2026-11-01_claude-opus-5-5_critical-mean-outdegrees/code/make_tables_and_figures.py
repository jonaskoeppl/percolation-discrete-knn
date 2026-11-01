"""Generate the simulation table, the bounds table and the two figures (pgfplots) from
../data/estimates.json and rigorous_bounds.py; writes the table_*.tex and figure_*.tex files."""
import json,math,io,contextlib,importlib.util
from fractions import Fraction
from decimal import Decimal,ROUND_HALF_UP
def rhu(x,n):
    return str((Decimal(x.numerator)/Decimal(x.denominator)).quantize(Decimal(1).scaleb(-n),rounding=ROUND_HALF_UP))
spec=importlib.util.spec_from_file_location('b','rigorous_bounds.py'); b=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(b)
fin=json.load(open('../data/sim_final.json'))
SIMD={2:1.84,3:1.38,4:1.20}
for _M,_d,_e,_r,_p in json.load(open('../data/chi_D.json')):
    if _d>=5: SIMD[_d]=_e
sim={}
for M in 'BU':
    for d in range(2,11):
        if f'{M}{d}' in fin: sim[(M,d)]=tuple(fin[f'{M}{d}'][:2])+(0,0)
# --- simulation table
t=[r'\begin{table}[ht]',r'\centering',r'\renewcommand{\arraystretch}{1.15}',r'\begin{tabular}{c'+'c'*9+'}',r'\toprule',
   r'$d$ & '+' & '.join(str(d) for d in range(2,11))+r'\\',r'\midrule']
def cell(M,d):
    c=sim[(M,d)][0]
    if M=='U' and d>=5: return f'{c:.4f}'
    return f'{c:.3f}'
t.append(r'$\bc^{\mathrm B}(d)$ & '+' & '.join(cell('B',d) for d in range(2,11))+r'\\')
t.append(r'$\bc^{\mathrm D}(d)$ & '+' & '.join((f'{SIMD[d]:.2f}' if d<5 else f'{SIMD[d]:.3f}') for d in range(2,11))+r'\\')
t.append(r'$\bc^{\mathrm U}(d)$ & '+' & '.join(cell('U',d) for d in range(2,11))+r'\\')
t.append(r'\midrule')
t.append(r'$\varphi^{-1}(2d)$ & '+' & '.join(f'{float(b.phiinv(2*d)):.3f}' for d in range(2,11))+r'\\')
t.append(r'$1+1/(2d)$ & '+' & '.join(rhu(1+Fraction(1,2*d),3) for d in range(2,11))+r'\\')
t.append(r'$1+1/(8d^2)$ & '+' & '.join(rhu(1+Fraction(1,8*d*d),4) for d in range(2,11))+r'\\')
t+= [r'\bottomrule',r'\end{tabular}']
t.append(r'\caption{Simulated critical mean out-degrees. For B and U: $d=2$ from the crossing of the spanning-threshold distributions for $L\in\{512,1024\}$, $d\ge3$ from the extrapolation of the mean cluster size; the uncertainty is about $\pm0.005$ for B and $\pm0.001$ for U; for $d\in\{9,10\}$ only the points with $\beta\ge1$ are used for U. For D: \cite[Table~1]{CHJK26} for $d\le4$, and the out-cluster extrapolation for $d\ge5$, with an uncertainty of about $0.02$. The last three rows are the rigorous lower bound for $\bc^{\mathrm B}$ and the leading terms for $\bc^{\mathrm D}$ and $\bc^{\mathrm U}$ from \cref{thm:asymptotics} and \cref{conj:undirected}.}\label{tab:sim}')
t.append(r'\end{table}')
open('../table_simulated_values.tex','w').write('\n'.join(t)+'\n')
# --- errors for text
with open('../data/sim_errors.txt','w') as f:
    for k,v in sorted(sim.items()): f.write(f'{k} {v}\n')
# --- bounds table
rows=[]
for d in [2,3,4,5,6,7,8,9,10,11,12,15,20]:
    m=2*d; pc=b.pcbar(d); c=b.cbar(d)
    Bl=float(b.phiinv(b.Fr(m*(m-1))/c)); Bu=float(b.phiinv(b.Fr(m*(m-1))*pc))
    Dl=float(b.Fr(m)/c); Du=float(min([b.Fr(m)*pc,b.Fr(2)]+([b.Fr(b.BAR[d])] if d in b.BAR else [])))
    lo,hi=b.bis(lambda x: float(c)*b.lam(x,m)-1,1,m-1); Ul=max(1.0,float(lo))
    qs=1-math.sqrt(1-float(pc)); Uu=min(Du,b.phiinv(qs*qs*m*(m-1)))
    rows.append((d,Bl,Bu,Dl,Du,Ul,Uu))
def fl(x): return math.floor(x*1e4)/1e4
def ce(x): return math.ceil(x*1e4)/1e4
t=[r'\begin{table}[ht]',r'\centering',r'\begin{tabular}{crrrrrr}',r'\toprule',
   r'$d$ & \multicolumn{2}{c}{$\bc^{\mathrm B}(d)$} & \multicolumn{2}{c}{$\bc^{\mathrm D}(d)$} & \multicolumn{2}{c}{$\bc^{\mathrm U}(d)$}\\',
   r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}',r' & lower & upper & lower & upper & lower & upper\\',r'\midrule']
for d,Bl,Bu,Dl,Du,Ul,Uu in rows:
    bu=f'{ce(Bu):.4f}' if d!=2 else r'$<3$'
    du=f'{ce(Du):.4f}' if d!=2 else r'$<2$'
    t.append(f'{d} & {fl(Bl):.4f} & {bu} & {fl(Dl):.4f} & {du} & {fl(Ul):.4f} & {ce(Uu):.4f}\\\\')
t+=[r'\bottomrule',r'\end{tabular}',r'\caption{Rigorous bounds from \cref{thm:absence,thm:bounds,thm:planar,thm:directed}, rounded outwards. For $d=2$ the upper bounds are \cref{thm:planar} and \cite[Theorem~2.3]{CHJK26}. For $d\ge6$ the lower bounds use only $c(d)<2d-1$.}\label{tab:bounds}',r'\end{table}']
open('../table_rigorous_bounds.tex','w').write('\n'.join(t)+'\n')
# --- figure 1: simulated curves
def coords(pairs): return ' '.join(f'({x},{y:.4f})' for x,y in pairs)
Bp=[(d,sim[('B',d)][0]) for d in range(2,11)]
Up=[(d,sim[('U',d)][0]) for d in range(2,11)]
Dp=[(d,v) for d,v in SIMD.items() if d<=10]
f1=r"""\begin{figure}[ht]
\centering
\begin{tikzpicture}
\begin{axis}[width=0.82\textwidth,height=0.5\textwidth,xmin=1.6,xmax=10.9,ymin=0.85,ymax=5.4,
  xlabel={dimension $d$},ylabel={critical mean out-degree $\bc(d)$},xtick={2,...,10},
  axis lines=left,tick style={black!50},grid=major,grid style={black!8},
  every axis plot/.append style={line width=1pt},legend style={draw=none,fill=none,font=\small},legend pos=north west,legend cell align=left]
\addplot[black!45,dashed,domain=2:10.5,samples=60,line width=0.6pt] {sqrt(2*x)+0.5};
\addlegendentry{$\sqrt{2d}+1/2$}
\addplot[black!45,dotted,domain=2:10.5,samples=60,line width=0.8pt] {1+1/(2*x)};
\addlegendentry{$1+1/(2d)$}
\addplot[serB,mark=*,mark size=2.4pt] coordinates {"""+coords(Bp)+r"""};
\addlegendentry{bidirectional (B)}
\addplot[serD,mark=square*,mark size=2.2pt] coordinates {"""+coords(Dp)+r"""};
\addlegendentry{directed (D)}
\addplot[serU,mark=triangle*,mark size=2.8pt] coordinates {"""+coords(Up)+r"""};
\addlegendentry{undirected (U)}
\node[serB,anchor=west,font=\small] at (axis cs:10.15,5.02) {B};
\node[serD,anchor=west,font=\small] at (axis cs:10.15,1.13) {D};
\node[serU,anchor=west,font=\small] at (axis cs:10.15,0.95) {U};
\end{axis}
\end{tikzpicture}
\caption{Simulated critical mean out-degrees of the bidirectional (circles), directed (squares) and undirected (triangles) lattice neighbour graphs, see \cref{tab:sim}. The dashed and dotted curves are the leading asymptotics of $\bc^{\mathrm B}$ and $\bc^{\mathrm D}$ from \cref{thm:asymptotics}. The uncertainties are smaller than the markers, except for D with $d\ge5$ (about $0.02$).}\label{fig:curves}
\end{figure}
"""
open('../figure_critical_curves.tex','w').write(f1)
# --- figure 2: rigorous bands
ds=[r[0] for r in rows]
def band(i,j): return ' '.join(f'({r[0]},{(r[i]+r[j])/2:.4f}) +- (0,{(r[j]-r[i])/2:.4f})' for r in rows)
f2=r"""\begin{figure}[ht]
\centering
\begin{tikzpicture}
\begin{axis}[name=left,width=0.5\textwidth,height=0.46\textwidth,xmin=1,xmax=21,ymin=2,ymax=8,
  xlabel={$d$},title={\small bidirectional},axis lines=left,tick style={black!50},grid=major,grid style={black!8},
  xtick={2,5,10,15,20}]
\addplot[black!45,dashed,domain=2:20.5,samples=80,line width=0.6pt] {sqrt(2*x)+0.5};
\addplot[serB,draw=none,mark=none,line width=1.2pt,error bars/.cd,y dir=both,y explicit,error mark=-,error bar style={line width=1.2pt,serB}] coordinates {"""+band(1,2)+r"""};
\addplot[black,only marks,mark=*,mark size=1.4pt] coordinates {"""+coords(Bp)+r"""};
\end{axis}
\begin{axis}[at={(left.east)},anchor=west,xshift=1.2cm,width=0.5\textwidth,height=0.46\textwidth,xmin=1,xmax=21,ymin=0.95,ymax=2.1,
  xlabel={$d$},title={\small directed and undirected},axis lines=left,tick style={black!50},grid=major,grid style={black!8},
  xtick={2,5,10,15,20}]
\addplot[black!45,dotted,domain=2:20.5,samples=80,line width=0.8pt] {1+1/(2*x)};
\addplot[serD,draw=none,mark=none,line width=1.2pt,error bars/.cd,y dir=both,y explicit,error mark=-,error bar style={line width=1.2pt,serD}] coordinates {"""+' '.join(f'({r[0]-0.18},{(r[3]+r[4])/2:.4f}) +- (0,{(r[4]-r[3])/2:.4f})' for r in rows)+r"""};
\addplot[serU,draw=none,mark=none,line width=1.2pt,error bars/.cd,y dir=both,y explicit,error mark=-,error bar style={line width=1.2pt,serU}] coordinates {"""+' '.join(f'({r[0]+0.18},{(r[5]+r[6])/2:.4f}) +- (0,{(r[6]-r[5])/2:.4f})' for r in rows)+r"""};
\addplot[black,only marks,mark=square*,mark size=1.3pt] coordinates {"""+' '.join(f'({d-0.18},{v})' for d,v in SIMD.items())+r"""};
\addplot[black,only marks,mark=triangle*,mark size=1.7pt] coordinates {"""+' '.join(f'({d+0.18},{v:.4f})' for d,v in Up)+r"""};
\node[serD,font=\small,anchor=west] at (axis cs:11,1.85) {D (left bars)};
\node[serU,font=\small,anchor=west] at (axis cs:11,1.7) {U (right bars)};
\end{axis}
\end{tikzpicture}
\caption{Rigorous bounds of \cref{tab:bounds} (vertical bars) and simulated values (black markers) for $d\in\{2,\dots,12,15,20\}$. Left: bidirectional graph with $\sqrt{2d}+1/2$ (dashed). Right: directed (left bars, squares) and undirected (right bars, triangles) graph with $1+1/(2d)$ (dotted). For $d\ge11$ the bounds use the computer-assisted estimates on $\pc(d)$ of \cite{FvdH17}.}\label{fig:bounds}
\end{figure}
"""
open('../figure_rigorous_bounds.tex','w').write(f2)
print('ok', sim)

# (the integer-case table table_integer_case.tex is written by hand; see the consistency check classify_check in the paper's code folder)
