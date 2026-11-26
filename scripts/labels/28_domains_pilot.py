#!/usr/bin/env python3
"""Piloto de anotación por dominio Pfam → GO para proteínas no_data."""
from __future__ import annotations
import argparse, csv, json, subprocess
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def fasta_head(src,dst,n):
    count=0
    with src.open(encoding="utf-8",errors="replace") as fh, dst.open("w") as out:
        for line in fh:
            if line.startswith(">"):
                if count>=n: break
                count+=1
            out.write(line)
    return count


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--indir",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--n",type=int,default=10000)
    ap.add_argument("--threads",type=int,default=4)
    a=ap.parse_args(); work=a.indir/"domains_pilot"; work.mkdir(parents=True,exist_ok=True)
    n=fasta_head(a.indir/"unannotated.faa",work/"query.faa",a.n)
    dom=work/"domains.domtbl"
    subprocess.run(["hmmscan","--cpu",str(a.threads),"--domtblout",str(dom),"--noali",
                    str(a.data_root/"db/bakta_db/pfam"),str(work/"query.faa")],check=True,
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    # GO names/aspect
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
    pmap={r["protein_id"]:r for r in read_tsv(a.indir/"proteins_map.tsv")}
    rows=[]; seen=set()
    for line in dom.read_text(errors="replace").splitlines():
        if not line or line.startswith("#"): continue
        f=line.split()
        if len(f)<22: continue
        try: ie=float(f[12])
        except ValueError: continue
        if ie>1e-5: continue
        pf=f[1].split(".")[0]; q=f[3]; m=pmap.get(q)
        if not m: continue
        for go in pfam.get(pf,[]):
            key=(m["genome"],m["gene_id"],go)
            if key in seen: continue
            seen.add(key)
            rows.append({"genome":m["genome"],"gene_id":m["gene_id"],"protein_id":q,
                         "pfam":pf,"go_id":go,"go_name":names.get(go,""),
                         "aspect":ns.get(go,"")[:1].upper() if ns.get(go) else "",
                         "evidence":"IEA","source":"Pfam2GO","ievalue":ie,
                         "confidence":"medium" if ie<1e-20 else "low"})
    out=work/"gene_domain_evidence.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"queried":n,"domain_rows":len(rows),"genes_with_domain_go":len({(r["genome"],r["gene_id"]) for r in rows}),
             "pfams":len({r["pfam"] for r in rows}),"terms":len({r["go_id"] for r in rows})}
    (work/"domain_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
