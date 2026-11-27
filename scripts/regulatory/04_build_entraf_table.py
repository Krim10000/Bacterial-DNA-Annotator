#!/usr/bin/env python3
"""Normaliza ENTRAF con hit_type/confidence y fuente explícita."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inp",type=Path,default=ROOT/"data/regulatory/external/entraf/entraf_entries.tsv")
    ap.add_argument("--out",type=Path,default=ROOT/"data/regulatory/entraf_entries.tsv")
    a=ap.parse_args()
    rows=[]
    for r in csv.DictReader(a.inp.open(encoding="utf-8",errors="replace"),delimiter="\t"):
        conf="high" if r.get("pubmed") and r.get("pdb") else ("medium" if r.get("pdb") else "low")
        rows.append({**r,"hit_type":"direct_species" if r.get("genome_match") else "no_match",
                     "confidence":conf,"inference":"none"})
    with a.out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["source_db"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"entries":len(rows),"matched":sum(1 for r in rows if r["hit_type"]=="direct_species"),
             "no_match":sum(1 for r in rows if r["hit_type"]=="no_match"),
             "by_genome":{},"by_role":{},"by_family":{}}
    for r in rows:
        if r.get("genome_match"): summary["by_genome"][r["genome_match"]]=summary["by_genome"].get(r["genome_match"],0)+1
        summary["by_role"][r.get("role","")]=summary["by_role"].get(r.get("role",""),0)+1
        summary["by_family"][r.get("family","")]=summary["by_family"].get(r.get("family",""),0)+1
    (a.out.parent/"entraf_summary_normalized.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1)[:2000])
if __name__=="__main__": main()
