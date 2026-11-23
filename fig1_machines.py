import pickle, numpy as np, csv, sys
DB=sys.argv[1] if len(sys.argv)>1 else 'A'
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from collections import defaultdict
rows=pickle.load(open("data/master.pkl","rb"))
K="data"; name={}   # KEGG KO list expected at data/ko_list.txt (https://rest.kegg.jp/list/ko); see README
for line in open(f"{K}/ko_list.txt"):
    p=line.rstrip('\n').split('\t')
    if len(p)>=2: name[p[0].split(':')[-1]]=p[1].lower()
def by_name(*kw): return set(k for k,v in name.items() if any(w in v for w in kw))
DON=['Asgardarchaeota','Planctomycetota','Myxococcota','Alphaproteobacteria']
wcol={'Alphaproteobacteria':'#c0552c','Myxococcota':'#2e8b2e','Planctomycetota':'#7a3a9a','Asgardarchaeota':'#1f78b4'}
dsl={d:np.array([float(x['sl']) for x in rows if x['donor']==d and x['db']==DB and x['three']=='True' and x['sl'] and float(x['sl'])>0]) for d in DON}
BAY="data"
rmode_tab={}
for r in csv.DictReader(open("data/reported_modes.tsv"),delimiter='\t'):
    rmode_tab[(r['database'],r['clade'])]=float(r['reported_mode'])
rmode={d:rmode_tab[('TOLDB'+DB,d)] for d in DON}
dmed={d:float(np.median(dsl[d])) for d in DON}
comp=[('Mitoribosome, large subunit',by_name('large subunit ribosomal protein')),
      ('Mitoribosome, small subunit',by_name('small subunit ribosomal protein')),
      ('Complex I (NADH dehydrogenase)',by_name('nadh dehydrogenase','nadh:ubiquinone','nadh-ubiquinone')),
      ('Complex II (succinate dehydrogenase)',by_name('succinate dehydrogenase')),
      ('Complex III (cytochrome bc1)',by_name('ubiquinol-cytochrome c','cytochrome c1','cytochrome b;')),
      ('Complex IV (cytochrome c oxidase)',by_name('cytochrome c oxidase'))]
alpha=[(d['KO'],d['db'],float(d['sl'])) for d in rows if d['donor']=='Alphaproteobacteria' and d['db']==DB and d['three']=='True' and d['KO'] and d['sl'] and float(d['sl'])>0]
units=[(c,[(ko,db,sl) for ko,db,sl in alpha if ko in kos]) for c,kos in comp]; units=[(c,p) for c,p in units if len(p)>=4]

fig=plt.figure(figsize=(9.6,7.6)); gs=fig.add_gridspec(2,1,height_ratios=[1.05,2.05],hspace=0.72)
axT=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1],sharex=axT); xs=np.linspace(np.log(0.004),np.log(3.2),500)
for d in DON:
    v=dsl[d]; k=gaussian_kde(np.log(v)); y=k(xs)
    axT.fill_between(np.exp(xs),y,color=wcol[d],alpha=0.18); axT.plot(np.exp(xs),y,color=wcol[d],lw=1.2)
    axT.plot([dmed[d]],[np.interp(np.log(dmed[d]),xs,y)],marker='o',color=wcol[d],ms=7,mec='k',mew=0.5,zorder=6,label=d)
    axT.axvline(rmode[d],color=wcol[d],ls=':',lw=1.3,zorder=3)
axT.legend(fontsize=8,loc='upper left',framealpha=0.9,handletextpad=0.5)
axT.set_ylabel('density',fontsize=9); axT.set_yticks([]); axT.spines[['right','top']].set_visible(False)
for d in DON: axB.axvline(dmed[d],color=wcol[d],ls='--',lw=1.2,zorder=1)
rng=np.random.default_rng(0)
for yi,(c,pairs) in enumerate(units):
    v=np.array([sl for ko,db,sl in pairs])
    y=yi+rng.uniform(-0.17,0.17,len(v)); axB.scatter(v,y,s=18,color='#8a8f96',alpha=0.6,edgecolor='none',zorder=2)
    lo,hi=np.percentile(v,2.5),np.percentile(v,97.5); axB.plot([lo,hi],[yi,yi],color='#555',lw=1.2,zorder=3)
    axB.scatter([np.median(v)],[yi],marker='o',s=55,color='k',zorder=5,edgecolor='w',linewidth=0.6)
    if 'Complex II' in c:
        sb=sorted(sl for ko,db,sl in pairs if ko=='K00235' and db==DB)
        sc=sorted(sl for ko,db,sl in pairs if ko=='K00236' and db==DB)
        if sb and sc:
            xb,xc=sb[0],sc[0]
            axB.scatter([xb],[yi],marker='D',s=62,color='#c0392b',zorder=7,edgecolor='w',linewidth=0.6)
            axB.scatter([xc],[yi],marker='D',s=62,color='#1f6fb0',zorder=7,edgecolor='w',linewidth=0.6)
            axB.text(xb,yi-0.24,'SdhB',ha='center',va='bottom',fontsize=7.2,color='#c0392b')
            axB.text(xc,yi-0.24,'SdhC',ha='center',va='bottom',fontsize=7.2,color='#1f6fb0')
            axB.annotate('',xy=(xc,yi+0.28),xytext=(xb,yi+0.28),arrowprops=dict(arrowstyle='<->',color='#333',lw=1.0))
axB.set_yticks(range(len(units))); axB.set_yticklabels([c for c,_ in units],fontsize=8.8)
axB.set_ylim(-0.9,len(units)-0.4); axB.invert_yaxis(); axB.set_xscale('log'); axB.set_xlim(0.004,3.2)
axB.set_xlabel('normalized stem length  sl = rsl/ebl  (log scale)',fontsize=9)
axB.spines[['right','top']].set_visible(False); plt.tight_layout()
for e in ('png','svg'): plt.savefig(f"figures/Fig1_machines_{DB}.{e}",dpi=160,bbox_inches='tight')
cm=[float(np.median([sl for ko,db,sl in p])) for c,p in units]
print("complex medians:",{c.split(" (")[0]:round(float(np.median([sl for ko,db,sl in p])),3) for c,p in units})
print(f"[{DB}] complex median span {max(cm)/min(cm):.1f}x  | donor median span {max(dmed.values())/min(dmed.values()):.1f}x")
print("donor medians:",{d[:5]:round(dmed[d],3) for d in DON}," reported:",{d[:5]:round(rmode[d],3) for d in DON})
