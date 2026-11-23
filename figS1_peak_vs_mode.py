import pickle, numpy as np, csv
from collections import defaultdict
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
rows=pickle.load(open("data/master.pkl","rb"))
rmode={}
for r in csv.DictReader(open("data/reported_modes.tsv"),delimiter='\t'):
    rmode[(r['database'],r['clade'])]=float(r['reported_mode'])
DON=['Alphaproteobacteria','Myxococcota','Planctomycetota','Asgardarchaeota']
fig,axs=plt.subplots(2,2,figsize=(10,6.4)); axs=axs.ravel()
xs=np.linspace(np.log(0.004),np.log(3.2),600)
for ax,d in zip(axs,DON):
    v=np.array([float(x['sl']) for x in rows if x['donor']==d and x['db']=='A' and x['three']=='True' and x['sl'] and float(x['sl'])>0])
    k=gaussian_kde(np.log(v)); y=k(xs); sl=np.exp(xs)
    ax.fill_between(sl,y,color='#9ab',alpha=0.35); ax.plot(sl,y,color='#4a6785',lw=1.3)
    peak=np.exp(xs[np.argmax(y)]); pm=rmode[('TOLDBA',d)]
    med=np.median(v); mean=v.mean()
    ax.axvline(peak,color='#1a7a1a',lw=2,label='peak')
    ax.axvline(pm,color='#c0392b',lw=2,ls='--',label='reported mode')
    ax.axvline(med,color='#666',lw=1,ls=':',label='median')
    ax.axvline(mean,color='#666',lw=1,ls='-.',label='mean')
    ax.set_xscale('log'); ax.set_xlim(0.004,3.2); ax.set_yticks([])
    ax.set_title(d,fontsize=10); ax.legend(fontsize=7.6,loc='upper right',framealpha=0.9)
    ax.set_xlabel('normalized stem length sl (log)',fontsize=8)
plt.tight_layout()
for e in ('png','svg'): plt.savefig(f"figures/figS1_peak_vs_mode.{e}",dpi=150)
print("saved")
