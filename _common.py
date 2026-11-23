"""Shared helpers: load per-phylum enrichment, sampling and branch-length data."""
import pickle, numpy as np
from collections import Counter
DATA="data/"
def load_tax(): return pickle.load(open(DATA+"tax.pkl","rb"))      # (sp2, prok_phyla)
def key_of(t):  return t[2] if t[1]=='Proteobacteria' else t[1]    # class for Proteobacteria, else phylum
def genome_counts():
    sp2,_=load_tax(); g=Counter()
    for sp,t in sp2.items():
        if t[0] in ('Bacteria','Archaea') and t[1] and t[1]!='NA': g[key_of(t)]+=1
    return g
def wilson(k,n,z=1.96):
    p=k/n; den=1+z*z/n; c=(p+z*z/2/n)/den; h=z*np.sqrt(p*(1-p)/n+z*z/4/n/n)/den; return c-h,c+h
def enrichment_ci():
    """observed enrichment + 95% CI per phylum, from lgt_{A,B,C}.pkl (euk/baseline neighbour counts)."""
    ec=Counter(); ea=Counter(); bc=Counter(); ba=Counter()
    for db in "ABC":
        d=pickle.load(open(DATA+f"lgt_{db}.pkl","rb")); ec+=d['ec']; ea+=d['ea']; bc+=d['bc']; ba+=d['ba']
    out={}
    for P in set(list(ec)+list(bc)):
        if ea[P]<150 or ba[P]<150: continue
        elo,ehi=wilson(ec[P],ea[P]); blo,bhi=wilson(bc[P],ba[P])
        out[P]=((ec[P]/ea[P])/(bc[P]/ba[P]), elo/bhi, ehi/blo)
    return out
def rate_enrich_genomes(min_leaves=2000):
    """aligned arrays: keys, mean terminal branch length (xs), enrichment (ys), genome count (gc)."""
    rate=pickle.load(open(DATA+"phylum_rate.pkl","rb"))            # key -> (geomean_bl, mean_bl, n_leaves)
    keys,xs,ys=pickle.load(open(DATA+"rate_enrich_merge.pkl","rb"))
    g=genome_counts()
    xs=np.array(xs); ys=np.array(ys); gc=np.array([g.get(k,0) for k in keys])
    m=(gc>0)
    keys=[k for k,o in zip(keys,m) if o]
    return keys, xs[m], ys[m], gc[m]
