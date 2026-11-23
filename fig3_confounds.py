#!/usr/bin/env python3
"""Figure 3 — euk-neighbour enrichment rises with (a) genome sampling and (b) branch length.
Run: python fig3_confounds.py   Inputs: data/{tax.pkl, phylum_rate.pkl, rate_enrich_merge.pkl}"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from _common import rate_enrich_genomes
keys,xs,ys,gc=rate_enrich_genomes(); idx={k:i for i,k in enumerate(keys)}
# combined R^2 (reported in caption)
le=np.log(ys); X=np.column_stack([np.ones(len(le)),np.log(gc),np.log(xs)]); bb,_,_,_=np.linalg.lstsq(X,le,rcond=None)
R2=1-np.sum((le-X@bb)**2)/np.sum((le-le.mean())**2)
print(f"combined genomes+branch R2 = {R2:.2f}")
cc={'Alphaproteobacteria':'#1f6fb0','Asgardarchaeota':'#1f6fb0','Myxococcota':'#c0392b','Planctomycetota':'#c0392b','Chlamydiota':'#7a3a9a'}
ll={'Alphaproteobacteria':'Alphaproteo.','Asgardarchaeota':'Asgardarch.','Myxococcota':'Myxococcota','Planctomycetota':'Planctomycetota','Chlamydiota':'Chlamydiota'}
fig,(a1,a2)=plt.subplots(1,2,figsize=(9.0,4.4))
def panel(ax,Xv,off,xlabel,title):
    ax.scatter(Xv,ys,s=16,c='#c2c7ce',edgecolor='none',zorder=2)
    lx=np.log(Xv); m=np.polyfit(lx,np.log(ys),1); xx=np.linspace(lx.min(),lx.max(),50)
    ax.plot(np.exp(xx),np.exp(m[0]*xx+m[1]),color='#555',lw=1.3,zorder=3)
    for k,c in cc.items():
        i=idx[k]; ax.scatter([Xv[i]],[ys[i]],s=64,c=c,edgecolor='w',lw=0.7,zorder=6)
        dx,dy,ha=off.get(k,(1,1.1,'center')); ax.annotate(ll[k],(Xv[i]*dx,ys[i]*dy),fontsize=7.6,color=c,ha=ha)
    ax.axhline(1,ls=':',color='#aaa',lw=0.9); ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(xlabel,fontsize=9); ax.set_ylabel('euk-neighbour enrichment',fontsize=9); ax.set_title(title,fontsize=9.2)
    ax.spines[['right','top']].set_visible(False)
panel(a1,gc,{'Asgardarchaeota':(0.9,1.12,'right'),'Myxococcota':(1.05,1.13,'left'),'Planctomycetota':(1.05,1.1,'left'),'Alphaproteobacteria':(1,1.11,'center'),'Chlamydiota':(1,1.13,'center')},
      'genomes sampled (log)','(a)')
panel(a2,xs,{'Asgardarchaeota':(1.06,0.89,'left'),'Myxococcota':(0.94,1.12,'right'),'Planctomycetota':(1,0.88,'center'),'Alphaproteobacteria':(1,1.11,'center'),'Chlamydiota':(0.95,1.12,'right')},
      'mean terminal branch length (log)','(b)')
plt.tight_layout()
for e in ('png','svg'): plt.savefig(f"figures/Fig3_confounds.{e}",dpi=170,bbox_inches='tight')
print("wrote figures/Fig3_confounds.png/svg")
