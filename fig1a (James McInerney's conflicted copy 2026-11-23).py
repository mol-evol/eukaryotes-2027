#!/usr/bin/env python3
"""Helper + BUILDER: newick parser, species->taxon map, and the Fig 1a nearest-phylum
counts from the deposited gene trees (Zenodo 13897946). Imported by the other builders.
The figure scripts run from the bundled data/ and do NOT need the trees. Set EX to the
deposited LECA_origins folder."""
import pickle, re, sys, time
sp2,prok_phyla=pickle.load(open("data/tax.pkl","rb"))
EX="LECA_origins"           # <-- deposited LECA_origins folder (gene trees)
NHX=re.compile(r'\[&&NHX[^\]]*\]')
def get_sp(h):
    if h.count('_')>1 and '.' not in h: return h.rsplit('_',1)[0]
    return h.split('_',1)[0]

def parse(nw):
    nw=NHX.sub('',nw).strip().rstrip(';')
    par=[-1]; dist=[0.0]; name=[None]; kids=[[]]
    cur=0; i=0; n=len(nw)
    def newnode(p):
        par.append(p); dist.append(0.0); name.append(None); kids.append([])
        idx=len(par)-1; kids[p].append(idx); return idx
    while i<n:
        c=nw[i]
        if c=='(':
            cur=newnode(cur); i+=1
        elif c==',':
            cur=newnode(par[cur]); i+=1
        elif c==')':
            cur=par[cur]; i+=1
            # optional internal label
            j=i
            while j<n and nw[j] not in '(),;': j+=1
            lab=nw[i:j]
            if ':' in lab:
                nm,_,dd=lab.partition(':')
                try: dist[cur]=float(dd)
                except: pass
            i=j
        else:
            j=i
            while j<n and nw[j] not in '(),;': j+=1
            tok=nw[i:j]
            nm,_,dd=tok.partition(':')
            name[cur]=nm
            try: dist[cur]=float(dd)
            except: pass
            i=j
    return par,dist,name,kids

def analyze_db(db, target_class, start=0, count=10**9):
    fn=f"{EX}/TOLDB{db}_LECA_trees.nwk"
    qual=0; nearest=[]; total=0; seen=0
    for line in open(fn):
        seen+=1
        if seen<=start: continue
        if total>=count: break
        total+=1
        parts=line.rstrip('\n').split('\t')
        if len(parts)<4: continue
        try: par,dist,name,kids=parse(parts[3])
        except Exception: continue
        leaves=[i for i in range(len(par)) if not kids[i]]
        # taxonomy per leaf
        ltax={}
        for i in leaves:
            if name[i] is None: continue
            t=sp2.get(get_sp(name[i]))
            if t: ltax[i]=t
        focal=[i for i in leaves if i in ltax and (ltax[i][2]==target_class or ltax[i][1]==target_class)]
        gen=set(get_sp(name[i]) for i in focal)
        if len(gen)<2: continue
        own_phylum=ltax[focal[0]][1]  # Proteobacteria for alpha
        # MRCA of focal (alpha) leaves; distance from that clade node to every leaf
        def ancestors(x):
            a=[]
            while x!=-1: a.append(x); x=par[x]
            return a
        common=set(ancestors(focal[0]))
        for f in focal[1:]: common&=set(ancestors(f))
        depth={}
        for cn in common:
            d=0; x=cn
            while x!=-1: d+=1; x=par[x]
            depth[cn]=d
        mrca=max(common,key=lambda c:depth[c])
        INF=float('inf'); dd=[INF]*len(par); dd[mrca]=0.0
        stack=[mrca]
        while stack:
            u=stack.pop()
            for v in kids[u]:
                if dd[v]==INF: dd[v]=dd[u]+dist[v]; stack.append(v)
            p=par[u]
            if p!=-1 and dd[p]==INF: dd[p]=dd[u]+dist[u]; stack.append(p)
        best=None;bestd=INF
        for i in leaves:
            t=ltax.get(i)
            if not t: continue
            dom,phy,cls=t
            if dom not in ('Bacteria','Archaea'): continue
            if phy==own_phylum or not phy or phy=='NA': continue
            if dd[i]<bestd: bestd=dd[i]; best=phy
        if best is not None:
            qual+=1; nearest.append(best)
    return total,qual,nearest

if __name__=="__main__":
    import sys
    from collections import Counter
    db=sys.argv[1]; tcls=sys.argv[2]; start=int(sys.argv[3]); count=int(sys.argv[4])
    t0=time.time()
    total,qual,nearest=analyze_db(db,tcls,start,count)
    pickle.dump((total,qual,nearest),open(f"res_{db}_{tcls}_{start}.pkl","wb"))
    print(f"db {db} {tcls} slice[{start}:+{count}] processed={total} qual={qual} time={time.time()-t0:.1f}s phyla={len(set(nearest))}")
