#!/usr/bin/env python3
"""Ranking de candidatos por dominio, homología y contexto."""
import csv, json, re, os
from collections import defaultdict, Counter
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=Path(os.environ.get("TARGET_DIR", D/"data"/"metabolism"/"targets")); S=O/"search"
targets={r["target_id"]:r for r in csv.DictReader(open(O/"targets.tsv"),delimiter="\t")}
acc_t=defaultdict(set); dom_t=defaultdict(set)
for r in csv.DictReader(open(O/"refs.tsv"),delimiter="\t"):
    acc_t[r["acc"]].add(r["target_id"])
    for d in r["domains"].split(";"):
        if d: dom_t[d].add(r["target_id"])
def genes(g):
    out=[]; p=D/"data"/"annotations_ab_initio"/g/"annotation.gff3"
    if p.exists():
        for l in open(p):
            if l.startswith("#"): continue
            f=l.split("\t")
            if len(f)<9 or f[2]!="CDS": continue
            mid=re.search(r"ID=([^;]+)",f[8])
            if mid: out.append({"pid":mid.group(1),"contig":f[0],"start":int(f[3]),"end":int(f[4]),"strand":f[6]})
    return out
def ann(g):
    m={}; p=D/"data"/"annotations_ab_initio"/g/"functional.tsv"
    if p.exists():
        for line in open(p):
            f=line.rstrip("\n").split("\t")
            if len(f)>=13: m[f[0]]=f[12]
    return m
rows=[]
genomes=sorted({r["genome"] for r in csv.DictReader(open(O/"target_genomes.tsv"),delimiter="\t")})
for gi,g in enumerate(genomes,1):
    anns=ann(g); gl=genes(g); cords={x["pid"]:(x["contig"],x["start"],x["end"],x["strand"]) for x in gl}; local=defaultdict(lambda:{"hom":None,"dom":None})
    for line in open(S/f"{g}.diamond.tsv",errors="ignore"):
        f=line.split("\t")
        if len(f)<6: continue
        pid,acc,pident,qcov,evalue,bits=f[:6]
        for t in acc_t.get(acc,[]):
            v=local[(pid,t)]; s=float(pident)*float(qcov)
            if v["hom"] is None or s>v["hom"][0]: v["hom"]=(s,acc,pident,qcov,evalue)
    for line in open(S/f"{g}.domtbl",errors="ignore"):
        if line.startswith("#"): continue
        f=line.split()
        if len(f)<22: continue
        pid,dom=f[0],f[3]
        try: ie=float(f[12]); cov=(int(f[16])-int(f[15])+1)/int(f[5])
        except: continue
        if ie>1e-5: continue
        for t in dom_t.get(dom,[]):
            v=local[(pid,t)]
            if v["dom"] is None or cov>v["dom"][0]: v["dom"]=(cov,dom,ie)
    for (pid,t),v in local.items():
        h=v["hom"]; d=v["dom"]
        if not h and not d: continue
        hs=(h[0]/10000) if h else 0; ds=d[0] if d else 0
        a=anns.get(pid,""); un=0 if a and not re.search(r"hypothetical|uncharacterized|DUF",a,re.I) else 1
        c=cords.get(pid,("","","","")); ctx=0.0; neigh=""
        if c[0]:
            ns=[x for x in gl if x["contig"]==c[0] and x["pid"]!=pid and not (x["end"]<c[1]-10000 or x["start"]>c[2]+10000)]
            ns.sort(key=lambda x:min(abs(x["start"]-c[2]),abs(c[1]-x["end"])))
            if ns:
                near=ns[0]; dist=min(abs(near["start"]-c[2]),abs(c[1]-near["end"]))
                ann_near=anns.get(near["pid"],"")
                neigh=f"{near['pid']}({dist}bp)"
                if ann_near: ctx+=0.5 if dist<=5000 else 0.25
                if near["strand"]==c[3]: ctx+=0.25
                ctx=min(ctx,1.0)
        score=0.4*ds+0.25*hs+0.15*un+0.2*ctx
        if score<0.15: continue
        rows.append({"target_id":t,"map_id":targets[t]["map_id"],"pathway":targets[t]["pathway"],"target":targets[t]["target"],
                     "genome":g,"protein_id":pid,"contig":c[0],"start":c[1],"end":c[2],"strand":c[3],
                     "annotation":a,"context_score":round(ctx,3),"nearest_gene":neigh,"score":round(score,4),"domain_score":round(ds,3),"homology_score":round(hs,3),
                     "best_acc":h[1] if h else "","pident":h[2] if h else "","qcov":h[3] if h else "","evalue":h[4] if h else "",
                     "known_exception":targets[t].get("known_exception",""),"example_taxon":targets[t].get("example_taxon",""),
                     "domains":d[1] if d else "","evidence":"domain+homology" if h and d else "domain" if d else "homology"})
    if gi%10==0: print(gi,"/",len(genomes),flush=True)
rows.sort(key=lambda x:(x["target_id"],-x["score"]))
top=[]; seen=Counter()
for r in rows:
    if seen[r["target_id"]]>=10: continue
    seen[r["target_id"]]+=1; top.append(r)
with open(O/"target_candidates.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(top[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(top)
want={(r["genome"],r["protein_id"]) for r in top}; seq={}
for g in {x[0] for x in want}:
    keep=False; hdr=None; buf=[]
    for line in open(D/"data"/"annotations_ab_initio"/g/"proteins.faa"):
        if line.startswith(">"):
            if hdr and (g,hdr) in want: seq[(g,hdr)]="".join(buf)
            hdr=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if hdr and (g,hdr) in want: seq[(g,hdr)]="".join(buf)
with open(O/"candidate_proteins.faa","w") as fh:
    for (g,pid),s in seq.items(): fh.write(f">{g}|{pid}\n{s}\n")
summary={"targets":len(targets),"genomes":len({r['genome'] for r in top}),"candidates":len(rows),"top":len(top),
         "by_target":dict(Counter(r["target_id"] for r in top)),"by_evidence":dict(Counter(r["evidence"] for r in top))}
json.dump(summary,open(O/"summary.json","w"),indent=2)
print(summary)
