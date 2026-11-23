#!/usr/bin/env python3
"""Excess-over-confounds ("positive controls") analysis behind the Fig. 2 / Results
numbers. Prints, for each lineage, observed enrichment and the ratio of observed to
the value predicted from genome sampling and/or branch length.
Run: python excess_residuals.py   Inputs: data/{tax.pkl, phylum_rate.pkl, rate_enrich_merge.pkl}"""
import numpy as np
from _common import rate_enrich_genomes
keys,xs,ys,gc=rate_enrich_genomes(); idx={k:i for i,k in enumerate(keys)}; le=np.log(ys)
def fit(cols):
    X=np.column_stack([np.ones(len(le))]+cols); b,_,_,_=np.linalg.lstsq(X,le,rcond=None); return np.exp(X@b)
mod={'genomes':fit([np.log(gc)]),'branch':fit([np.log(xs)]),'both':fit([np.log(gc),np.log(xs)])}
foc=['Alphaproteobacteria','Asgardarchaeota','Myxococcota','Planctomycetota','Chlamydiota']
print(f"{'phylum':20s} {'enrich':>7s} {'exc(genomes)':>12s} {'exc(branch)':>12s} {'exc(both)':>10s}")
for k in foc:
    i=idx[k]
    print(f"{k:20s} {ys[i]:7.2f} {ys[i]/mod['genomes'][i]:11.2f}x {ys[i]/mod['branch'][i]:11.2f}x {ys[i]/mod['both'][i]:9.2f}x")
