#!/usr/bin/env python3
"""BUILDER: per-phylum mean terminal branch length (data/phylum_rate.pkl) and the
rate+enrichment merge (data/rate_enrich_merge.pkl), from the deposited gene trees.
INPUT: Zenodo 13897946 gene trees. Set DEPOSIT_TREES to the folder of TOLDB{A,B,C}.nwk
(tab-separated; newick in the last field). Needs data/tax.pkl and data/lgt_{A,B,C}.pkl.
The figure scripts run from the bundled .pkl outputs without re-running this."""
import math, pickle, numpy as np
from collections import defaultdict, Counter
import fig1a
from fig1a import parse, get_sp, sp2
DEPOSIT_TREES="trees_full"   # <-- set to the deposited gene-tree folder
def key_of(t): return t[2] if t[1]=='Proteobacteria' else t[1]
slog=defaultdict(float); cnt=defaultdict(int); lin=defaultdict(float)
for db in "ABC":
    for line in open(f"{DEPOSIT_TREES}/TOLDB{db}.nwk"):
        nw=line.rstrip('\n').split('\t')[-1]
        if '(' not in nw: continue
        try: par,dist,name,kids=parse(nw)
        except Exception: continue
        for i in range(len(par)):
            if kids[i] or name[i] is None: continue
            tax=sp2.get(get_sp(name[i]))
            if not tax or tax[0] not in ('Bacteria','Archaea') or not tax[1] or tax[1]=='NA': continue
            bl=dist[i]
            if bl and bl>0:
                k=key_of(tax); slog[k]+=math.log(bl); lin[k]+=bl; cnt[k]+=1
rate={k:(math.exp(slog[k]/cnt[k]), lin[k]/cnt[k], cnt[k]) for k in cnt}
pickle.dump(rate,open("data/phylum_rate.pkl","wb"))
# merge with availability-controlled enrichment
ec=Counter(); ea=Counter(); bc=Counter(); ba=Counter()
for db in "ABC":
    d=pickle.load(open(f"data/lgt_{db}.pkl","rb")); ec+=d['ec']; ea+=d['ea']; bc+=d['bc']; ba+=d['ba']
keys=[]; xs=[]; ys=[]
for P in set(list(ec)+list(bc)):
    if ea[P]<150 or ba[P]<150 or P not in rate or rate[P][2]<2000: continue
    keys.append(P); xs.append(rate[P][0]); ys.append((ec[P]/ea[P])/(bc[P]/ba[P]))
pickle.dump((keys,xs,ys),open("data/rate_enrich_merge.pkl","wb"))
print(f"wrote data/phylum_rate.pkl ({len(rate)} phyla) and data/rate_enrich_merge.pkl ({len(keys)} phyla)")
