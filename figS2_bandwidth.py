import pickle, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, gamma
rows=pickle.load(open("data/master.pkl","rb"))
DON=['Asgardarchaeota','Planctomycetota','Myxococcota','Alphaproteobacteria']
col={'Alphaproteobacteria':'#c0552c','Myxococcota':'#2e8b2e','Planctomycetota':'#7a3a9a','Asgardarchaeota':'#1f78b4'}
dsl={d:np.array([float(x['sl']) for x in rows if x['donor']==d and x['db']=='A' and x['three']=='True' and x['sl'] and float(x['sl'])>0]) for d in DON}
def pipe_mode(v,adj):
    v=v[(v>=np.quantile(v,0.01))&(v<=np.quantile(v,0.99))]
    k=gaussian_kde(v); k.set_bandwidth(k.factor*adj)
    xs=np.linspace(v.min(),v.max(),1024); y=k(xs); pk=np.argmax(y); dy=np.diff(y); cut=v.max()
    for i in range(pk+1,len(dy)):
        if dy[i-1]<0 and dy[i]>=0: cut=xs[i]; break
    w=v[v<=cut]
    if len(w)<8: return np.nan
    a,_,sc=gamma.fit(w,floc=0); return (a-1)*sc if a>1 else 0.0
adjs=np.linspace(0.5,2.5,17)
fig,ax=plt.subplots(figsize=(8,4.6))
curves={d:[pipe_mode(dsl[d],a) for a in adjs] for d in DON}
for d in DON: ax.plot(adjs,curves[d],'-o',ms=3,color=col[d],label=d)
ax.set_xlabel('kernel bandwidth (× default)')
ax.set_ylabel('reported “wave” mode (sl)')
ax.legend(fontsize=8,loc='upper right'); ax.spines[['right','top']].set_visible(False)
ax.axvline(1.0,color='#bbb',ls=':',lw=1); ax.text(1.02,ax.get_ylim()[1]*0.96,'default',fontsize=7,color='#888')
plt.tight_layout()
for e in ('png','svg'): plt.savefig(f"figures/figS2_bandwidth.{e}",dpi=150)
# print the reordering explicitly
print("bandwidth  bacterial order (earliest->latest, excluding host)")
for a in [0.5,0.75,1.0,1.5,2.0,2.5]:
    ms={d:pipe_mode(dsl[d],a) for d in DON if d!='Asgardarchaeota'}
    print(f"  x{a}: "+" > ".join(sorted(ms,key=lambda d:-ms[d])))
