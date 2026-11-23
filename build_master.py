# INPUT: the deposited Zenodo 13897946 data. Set DEPOSIT to its mLECA_OGs path.\n# This is a BUILDER (regenerates data/*.pkl); the figure scripts run from the bundled data/ without it.\n
import csv, sys
EX="/sessions/charming-modest-babbage/work/ex/mLECA_OGs"
GAB="/tmp/gab"
csv.field_size_limit(1<<30)

# COG -> category letters
cog2cat={}
for r in csv.reader(open(f"{GAB}/cog-20.def.tab",encoding='latin-1'),delimiter='\t'):
    if len(r)>=2: cog2cat[r[0]]=r[1]
INFO=set("JAKLB")   # information storage & processing
# operational = metabolism + cellular processing (exclude R,S poorly characterized)
OPER=set("DYVTMNZWUOXCGEFHIPQ")

def load_anno(db, kind):
    d={}
    fn=f"{EX}/annotations/TOLDB{db}_LECA_{kind}.tsv"
    for i,r in enumerate(csv.reader(open(fn),delimiter='\t')):
        if i==0 or len(r)<2: continue
        d[r[0]]=r[1]   # key 'Orthogroup|LECA' -> KO or COG
    return d

rows=[]
for db in "ABC":
    ko=load_anno(db,"KOs"); cog=load_anno(db,"COGs")
    with open(f"{EX}/LECA_origins/TOLDB{db}.tsv") as f:
        rd=csv.reader(f,delimiter='\t'); H=next(rd)
        idx={h:i for i,h in enumerate(H)}
        for r in rd:
            key=r[idx['Orthogroup']]+"|"+r[idx['LECA']]
            def g(c):
                v=r[idx[c]]; return v
            rsl=g('rsl'); sl=g('sl'); ebl=g('ebl')
            rows.append(dict(db=db, key=key, og=g('Orthogroup'), leca=g('LECA'),
                family_type=g('family_type'), three=g('Three_supergroups'),
                five=g('Five_supergroups'), evotype=g('LECA_evotype'),
                rsl=rsl, ebl=ebl, sl=sl, donor=g('donor'), mediator=g('mediator'),
                KO=ko.get(key,''), COG=cog.get(key,'')))
# attach categories
for d in rows:
    d['cats']=cog2cat.get(d['COG'],'')
import pickle
pickle.dump(rows,open("master.pkl","wb"))
print("total families (rows):", len(rows))
from collections import Counter
print("per db:", Counter(d['db'] for d in rows))
print("with KO:", sum(1 for d in rows if d['KO']))
print("Alpha-assigned (all):", sum(1 for d in rows if d['donor']=='Alphaproteobacteria'))
