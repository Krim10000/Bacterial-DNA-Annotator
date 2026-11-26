#!/usr/bin/env python3
"""Añade una fila no_data por gen sin evidencia en gene_evidence.tsv."""
from __future__ import annotations
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    data=ROOT/"data"; p=data/"labels/gene_evidence.tsv"; tmp=p.with_suffix(".tmp")
    have=set()
    with p.open(encoding="utf-8",errors="replace") as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        fields=rd.fieldnames
        for r in rd:
            if r.get("source")!="no_data": have.add((r.get("genome",""),r.get("gene_id","")))
    n=0
    with p.open(encoding="utf-8",errors="replace") as fh, tmp.open("w",newline="") as out:
        rd=csv.DictReader(fh,delimiter="\t"); w=csv.DictWriter(out,fieldnames=fields,delimiter="\t"); w.writeheader()
        for r in rd: w.writerow(r)
        for d in sorted((data/"regulation").iterdir()):
            if not d.is_dir() or not (d/"genes.tsv").exists(): continue
            g=d.name
            with (d/"genes.tsv").open(encoding="utf-8",errors="replace") as gf:
                for gene in csv.DictReader(gf,delimiter="\t"):
                    gid=gene.get("gene_id") or gene.get("locus_tag") or ""
                    if gid and (g,gid) not in have:
                        w.writerow({"genome":g,"gene_id":gid,"protein_id":gid,
                                    "source":"no_data","evidence":"no_data","aspect":"",
                                    "term_id":"","term_name":"","confidence":"unknown",
                                    "with_from":"","coverage":"none","notes":"sin evidencia GO"})
                        n+=1
    tmp.replace(p)
    print("added",n,"have",len(have))
if __name__=="__main__": main()
