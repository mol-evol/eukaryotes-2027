#!/usr/bin/env python3
"""Figure 2 — observed euk-neighbour enrichment vs the level predicted from genome
sampling + branch length. Run from the code/ folder: python fig2_enrichment.py
Inputs: data/{tax.pkl, phylum_rate.pkl, rate_enrich_merge.pkl, lgt_A/B/C.pkl}"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from _common import rate_enrich_genomes, enrichment_ci
keys,xs,ys,gc=rate_enrich_genomes(); idx={k:i for i,k in enumerate(keys)}
le=np.log(ys); X=np.column_stack([np.ones(len(le)),np.log(gc),np.log(xs)])
b,_,_,_=np.linalg.lstsq(X,le,rcond=None); pred=np.exp(X@b)      # confound prediction
ci=enrichment_ci()
foc=[('Chlamydiota','unproposed'),('Alphaproteobacteria','established'),('Myxococcota','proposed'),
     ('Asgardarchaeota','established'),('Gammaproteobacteria','unproposed'),('Cyanobacteria','unproposed'),
     ('Verrucomicrobiota','unproposed'),('Planctomycetota','proposed'),('Bacteroidota','unproposed'),('Chloroflexota','unproposed')]
foc=sorted(foc,key=lambda r:ci[r[0]][0])
col={'established':'#1f6fb0','proposed':'#c0392b','unproposed':'#9aa0a8'}
lab={}
fig,ax=plt.subplots(figsize=(8.4,4.9))
for i,(k,cat) in enumerate(foc):
    o,lo,hi=ci[k]; p=pred[idx[k]]; c=col[cat]
    ax.plot([p,o],[i,i],color=c,lw=1.5,alpha=0.45,zorder=2)
    ax.plot([lo,hi],[i,i],color=c,lw=1.0,alpha=0.8,zorder=3)
    ax.scatter([p],[i],s=42,facecolor='white',edgecolor=c,lw=1.6,zorder=4)
    ax.scatter([o],[i],s=72,color=c,zorder=5)
    ax.text(0.30,i,lab.get(k,k),ha='right',va='center',fontsize=8.3,color=c if cat!='unproposed' else '#444')
ax.axvline(1,ls=':',color='#999',lw=1)
ax.set_xscale('log'); ax.set_xlim(0.30,7); ax.set_ylim(-1.1,len(foc)-0.3); ax.set_yticks([])
ax.set_xticks([0.5,1,2,3,5]); ax.set_xticklabels(['0.5','1','2','3','5'])
ax.set_xlabel('euk-neighbour enrichment (95% CI)')
ax.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor='#444',markersize=8,label='observed'),
                   Line2D([0],[0],marker='o',color='w',markerfacecolor='white',markeredgecolor='#444',markersize=8,label='predicted')],
          loc='lower right',fontsize=8,framealpha=0.95)
ax.spines[['right','top']].set_visible(False); plt.tight_layout()
for e in ('png','svg'): plt.savefig(f"figures/Fig2_enrichment.{e}",dpi=170,bbox_inches='tight')
print("wrote figures/Fig2_enrichment.png/svg")
