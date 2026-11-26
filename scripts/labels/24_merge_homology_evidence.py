#!/usr/bin/env python3
"""Fusiona homología ISS en gene_evidence.tsv, sin duplicar no_data."""
from __future__ import annotations
import csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    data=ROOT/"data"; src=data/"labels/gene_evidence.tsv"
    hom=data/"labels/go_indirect/gene_go_homology_batches.tsv"
    out=data/"labels/gene_evidence.tsv"; tmp=data/"labels/gene_evidence.merge.tmp"
    names={}; cur=None
    for line in (data/"labels/external/go/go.obo").read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
    best={}
    with hom.open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            key=(r["genome"],r["gene_id"],r["go_id"])
            score=(float(r.get("identity") or 0),float(r.get("coverage") or 0))
            if key not in best or score>best[key][0]: best[key]=(score,r)
    hom_genes={(k[0],k[1]) for k in best}
    with src.open(encoding="utf-8",errors="replace") as fh, tmp.open("w",newline="") as of:
        rd=csv.DictReader(fh,delimiter="\t"); fields=rd.fieldnames
        w=csv.DictWriter(of,fieldnames=fields,delimiter="\t"); w.writeheader()
        n=0
        for r in rd:
            if r.get("source")=="no_data" and (r["genome"],r["gene_id"]) in hom_genes:
                continue
            w.writerow(r); n+=1
        for (g,gid,go),(score,r) in best.items():
            w.writerow({"genome":g,"gene_id":gid,"protein_id":gid,
                        "source":"homology_transfer","evidence":"ISS","aspect":r.get("aspect",""),
                        "term_id":go,"term_name":names.get(go,""),"confidence":"low",
                        "with_from":r.get("with_from",""),"coverage":"homology",
                        "notes":f"identity={score[0]:.1f}; qcov={score[1]:.1f}; donor={r.get('donor_source','')}"}); n+=1
    tmp.replace(out)
    genes=set(); nd=0
    with out.open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            if r["source"]=="no_data": nd+=1
            else: genes.add((r["genome"],r["gene_id"]))
    summary={"rows":n,"no_data_rows":nd,"genes_with_evidence":len(genes),
             "homology_rows":len(best),"homology_genes":len(hom_genes)}
    (data/"labels/gene_evidence_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
