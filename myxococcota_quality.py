#!/usr/bin/env python3
"""BUILDER/ANALYSIS: quality of the authors' Myxococcota donor assignments, from their
own deposited per-family tables (sister sequences + branch support). Reproduces the
numbers in the reply to W. Martin and the Myxococcota_assigned_genes.xlsx spreadsheet.
INPUT: Zenodo 13897946 LECA_origins tables. Set DEPOSIT to the mLECA_OGs/LECA_origins
folder (TOLDB{A,B,C}.tsv with columns donor, Three_supergroups, leca_sister_ufboot,
MABT_prop, donor_seqs). Needs data/tax.pkl."""
import csv, pickle, numpy as np
from collections import Counter
sp2,_=pickle.load(open("data/tax.pkl","rb"))
DEPOSIT="LECA_origins"   # <-- set to the deposited LECA_origins folder
def get_sp(h): return h.rsplit('_',1)[0] if (h.count('_')>1 and '.' not in h) else h.split('_',1)[0]
EUK_SG={'Diaphoretickes','Amorphea','Discoba','CRuMs','Malawimonadida','Ancyromonadida','Metamonada','Opimoda','Provora','Hemimastigophora'}
def contam(seqs):
    e=v=0
    for s in seqs.split(';'):
        if not s: continue
        t=sp2.get(get_sp(s))
        if t and (t[0]=='Eukaryota' or t[1] in EUK_SG): e+=1
        if t and ('viric' in (t[1] or '').lower() or t[1]=='Nucleocytoviricota'): v+=1
    return e+v
def ff(x):
    try: return float(x)
    except: return None
rows=[]
for db in "ABC":
    for r in csv.DictReader(open(f"{DEPOSIT}/TOLDB{db}.tsv"),delimiter='\t'):
        if r['donor']=='Myxococcota' and r['Three_supergroups']=='True': rows.append(r)
N=len(rows)
uf=[ff(r['leca_sister_ufboot']) for r in rows]; pr=[ff(r['MABT_prop']) for r in rows]
nseq=[len([s for s in r['donor_seqs'].split(';') if s]) for r in rows]
cont=sum(1 for r in rows if contam(r['donor_seqs'])>0)
both=sum(1 for r in rows if (ff(r['leca_sister_ufboot']) or 0)>=75 and (ff(r['MABT_prop']) or 0)>=0.75)
print(f"Myxococcota (3-supergroup) family assignments: N={N}")
print(f"  sister-branch UFBoot < 75:        {100*sum(1 for u in uf if u and u<75)/N:.0f}%")
print(f"  donor-proportion MABT_prop < 0.75: {100*sum(1 for p in pr if p is not None and p<0.75)/N:.0f}%")
print(f"  sister defined by <= 2 sequences:  {100*sum(1 for n in nseq if n<=2)/N:.0f}%")
print(f"  eukaryotes/viruses in donor seqs:  {100*cont/N:.0f}%")
print(f"  meet BOTH UFBoot>=75 AND prop>=0.75: {both}/{N} = {100*both/N:.0f}%")
