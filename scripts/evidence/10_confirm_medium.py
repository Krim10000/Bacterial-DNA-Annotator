"""Confirma candidatos de dominio medio con homología y evidencia EC/Rhea."""
import csv,json,subprocess,tempfile
from pathlib import Path
from collections import defaultdict,Counter
from concurrent.futures import ThreadPoolExecutor
D=Path(__file__).resolve().parents[2]
rows=[r for r in csv.DictReader(open(D/'data/metabolism/evidence/medium_review.tsv'),delimiter='\t') if r['category']=='strong_no_domain']
print("strong_no_domain:",len(rows),flush=True)
ref={}
for g in csv.DictReader(open(D/'out_fl/refined/refined_gaps.tsv'),delimiter='\t'): ref[(g['genome'],g['gap'])]=g
r2a=defaultdict(set)
for line in open(D/'data/metabolism/rhea2uniprot_sprot.tsv'):
    p=line.rstrip().split('\t')
    if len(p)>=4 and p[0]!='RHEA_ID': r2a['RHEA:'+p[0]].add(p[3])
pf=json.load(open(D/'data/metabolism/uniprot_pfam_cache.json'))
gap_dom={}; wanted=set()
for r in rows:
    ps=[x for x in ref[(r['genome'],r['gap'])].get('producers','').split(';') if x]
    accs=set().union(*(r2a.get(p,set()) for p in ps)) if ps else set()
    doms=set()
    for a in accs: doms |= set(pf.get(a,[]))
    gap_dom[(r['genome'],r['gap'])]=doms; wanted|=doms
print("dominios:",len(wanted),flush=True)
hmm=Path('/tmp/domain70.hmm'); block=[]; name=None; written=0
with open(D/'.cache/pfam/Pfam-A.hmm','r',errors='ignore') as fh, open(hmm,'w') as oh:
    for line in fh:
        if line.startswith('HMMER3'): block=[line]; name=None
        elif block:
            block.append(line)
            if line.startswith('NAME'): name=line.split(None,1)[1].strip()
            if line.startswith('//'):
                if name in wanted: oh.writelines(block); written+=1
                block=[]
print("HMMs",written,"bytes",hmm.stat().st_size,flush=True)
byg=defaultdict(list)
for r in rows: byg[r['genome']].append(r)
def run(g):
    faa=D/'data/annotations_ab_initio'/g/'proteins.faa'
    if not faa.exists(): return g,set()
    with tempfile.TemporaryDirectory() as td:
        o=Path(td)/'d.tbl'
        subprocess.run(['hmmsearch','--domtblout',str(o),'--noali','-E','1e-5',str(hmm),str(faa)],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
        hits=set()
        for line in open(o):
            if line.startswith('#'): continue
            f=line.split()
            if len(f)>6 and float(f[6])<=1e-5: hits.add(f[3])
        return g,hits
res={}
with ThreadPoolExecutor(max_workers=4) as ex:
    for i,(g,h) in enumerate(ex.map(run,byg),1):
        res[g]=h
        if i%20==0: print(f"  {i}/{len(byg)}",flush=True)
out=[]
for r in rows:
    k=(r['genome'],r['gap']); doms=gap_dom.get(k,set()); hits=doms & res.get(r['genome'],set())
    out.append({**r,'domain_confirmation':'domain_confirmed' if hits else ('no_domain_data' if not doms else 'structure_needed'),'domain_hits':';'.join(sorted(hits))})
with open(D/'data/metabolism/evidence/domain_confirm_medium.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(out)
print("estados:",dict(Counter(r['domain_confirmation'] for r in out)),flush=True)
