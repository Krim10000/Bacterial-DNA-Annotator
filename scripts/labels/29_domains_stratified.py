#!/usr/bin/env python3
"""Piloto estratificado de dominios Pfam para no_data (20/genoma, --cut_ga)."""
from __future__ import annotations
import argparse, csv, json, subprocess
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--indir",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--per-genome",type=int,default=20)
    ap.add_argument("--threads",type=int,default=6)
    a=ap.parse_args(); work=a.indir/"domains_stratified"; work.mkdir(parents=True,exist_ok=True)
    # no_data actual
    nodata=set()
    with (a.data_root/"labels/gene_evidence.tsv").open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            if r.get("source")=="no_data": nodata.add((r["genome"],r["gene_id"]))
    mapping={r["protein_id"]:r for r in read_tsv(a.indir/"proteins_map.tsv")}
    byg=defaultdict(list)
    with (a.indir/"unannotated.faa").open(encoding="utf-8",errors="replace") as fh:
        name=None; seq=[]
        for line in fh:
            line=line.strip()
            if line.startswith(">"):
                if name: byg[name.split("|",1)[0]].append((name,"".join(seq)))
                name=line[1:].split()[0]; seq=[]
            elif line and name is not None: seq.append(line)
        if name: byg[name.split("|",1)[0]].append((name,"".join(seq)))
    chosen=[]
    for g,items in sorted(byg.items()):
        n=0
        for pid,seq in items:
            m=mapping.get(pid)
            if not m: continue
            if (m["genome"],m["gene_id"]) not in nodata: continue
            chosen.append((pid,seq)); n+=1
            if n>=a.per_genome: break
    q=work/"query.faa"
    with q.open("w") as out:
        for pid,seq in chosen: out.write(f">{pid}\n{seq}\n")
    dom=work/"domains.domtbl"
    subprocess.run(["hmmscan","--cpu",str(a.threads),"--cut_ga","--domtblout",str(dom),
                    "--noali",str(a.data_root/"db/bakta_db/pfam"),str(q)],check=True,
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    names={}; ns={}; cur=None
    for line in (a.data_root/"labels/external/go/go.obo").read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
    pfam=defaultdict(list)
    for line in (a.data_root/"labels/external/go/pfam2go").read_text(encoding="utf-8",errors="replace").splitlines():
        if line.startswith("!"): continue
        left,right=line.split(">",1); acc=left.split()[0].split(":",1)[1].split(".")[0]
        go=right.split(";")[-1].strip()
        if go.startswith("GO:"): pfam[acc].append(go)
    pmap=mapping; rows=[]; seen=set()
    with dom.open(errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f=line.split()
            if len(f)<22: continue
            qid=f[3]; pf=f[1].split(".")[0]; m=pmap.get(qid)
            if not m: continue
            for go in pfam.get(pf,[]):
                key=(m["genome"],m["gene_id"],go)
                if key in seen: continue
                seen.add(key)
                rows.append({"genome":m["genome"],"gene_id":m["gene_id"],"protein_id":qid,
                             "pfam":pf,"go_id":go,"go_name":names.get(go,""),
                             "aspect":ns.get(go,"")[:1].upper() if ns.get(go) else "",
                             "evidence":"IEA","source":"Pfam2GO","confidence":"medium"})
    out=work/"gene_domain_evidence.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"sampled":len(chosen),"genomes":len({p.split('|',1)[0] for p,_ in chosen}),
             "domain_rows":len(rows),"genes_with_domain_go":len({(r["genome"],r["gene_id"]) for r in rows}),
             "pfams":len({r["pfam"] for r in rows}),"terms":len({r["go_id"] for r in rows}),
             "per_genome":a.per_genome,"mode":"cut_ga"}
    (work/"domain_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
