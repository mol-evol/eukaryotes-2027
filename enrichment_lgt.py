#!/usr/bin/env python3
"""BUILDER: regenerate data/lgt_{A,B,C}.pkl (euk / baseline nearest-neighbour counts)
from the deposited gene trees (Zenodo 13897946). The figure scripts run from the
bundled data/ and do NOT need this. Set TF to the gene-tree folder."""
import pickle, sys, random, fig1a
sp2,prok=pickle.load(open("data/tax.pkl","rb"))
TF="trees_full"             # <-- deposited gene-tree folder (TOLDB{A,B,C}.nwk)
random.seed(1)
def anc(x,par):
    a=[]
    while x!=-1: a.append(x); x=par[x]
    return a
def mrca(nodes,par):
    c=set(anc(nodes[0],par))
    for n in nodes[1:]: c&=set(anc(n,par))
    return max(c,key=lambda z:len(anc(z,par))) if c else None
def nearest_phylum(src,par,dist,kids,ltax,leaves,own):
    INF=float('inf'); dd=[INF]*len(par); dd[src]=0.0; st=[src]
    while st:
        u=st.pop()
        for v in kids[u]:
            if dd[v]==INF: dd[v]=dd[u]+dist[v]; st.append(v)
        p=par[u]
        if p!=-1 and dd[p]==INF: dd[p]=dd[u]+dist[u]; st.append(p)
    best=None;bd=INF
    for i in leaves:
        t=ltax.get(i)
        if not t: continue
        dom,phy,cls=t
        if dom not in ('Bacteria','Archaea') or not phy or phy=='NA': continue
        g=cls if phy=='Proteobacteria' else phy
        if g==own: continue
        if dd[i]<bd: bd=dd[i]; best=g
    return best
from collections import Counter
db=sys.argv[1]; NF=int(sys.argv[2]) if len(sys.argv)>2 else 8
euk_chosen=Counter(); euk_avail=Counter(); base_chosen=Counter(); base_avail=Counter()
for line in open(f"{TF}/TOLDB{db}.nwk"):
    p=line.rstrip('\n').split('\t')
    if len(p)<2: continue
    try: par,dist,name,kids=fig1a.parse(p[-1])
    except: continue
    leaves=[i for i in range(len(par)) if not kids[i]]
    ltax={i:sp2.get(fig1a.get_sp(name[i])) for i in leaves if name[i]}
    ltax={i:t for i,t in ltax.items() if t}
    # prok phyla present with >=2 leaves
    byphy={}
    for i in leaves:
        t=ltax.get(i)
        if t and t[0] in ('Bacteria','Archaea') and t[1] and t[1]!='NA':
            g=t[2] if t[1]=='Proteobacteria' else t[1]
            byphy.setdefault(g,[]).append(i)
    present=[ph for ph,ls in byphy.items() if len(ls)>=2]
    if len(present)<3: continue
    # euk clade
    euk=[i for i in leaves if ltax.get(i) and ltax[i][0]=='Eukaryota']
    if len(euk)>=2:
        em=mrca(euk,par)
        if em is not None:
            nb=nearest_phylum(em,par,dist,kids,ltax,leaves,own=None)
            for ph in present: euk_avail[ph]+=1
            if nb: euk_chosen[nb]+=1
    # baseline: sample focal prok phyla
    foc=random.sample(present,min(NF,len(present)))
    for f in foc:
        fm=mrca(byphy[f],par)
        if fm is None: continue
        nb=nearest_phylum(fm,par,dist,kids,ltax,leaves,own=f)
        for ph in present:
            if ph!=f: base_avail[ph]+=1
        if nb: base_chosen[nb]+=1
pickle.dump(dict(ec=euk_chosen,ea=euk_avail,bc=base_chosen,ba=base_avail),open(f"data/lgt_{db}.pkl","wb"))
print(f"db {db}: euk assignments={sum(euk_chosen.values())} baseline assignments={sum(base_chosen.values())}")
