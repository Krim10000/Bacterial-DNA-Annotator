#!/usr/bin/env python3
"""Anota rRNA/tRNA/ncRNA con Rfam2GO, independiente de proteínas."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from webapp.ncrna_go import build_for_genome

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/ncrna_go")
    a=ap.parse_args()
    genomes=[p.name[:-3] if p.name.endswith(".gz") else p.name for p in sorted((a.data_root/"annotations").glob("*.gff*"))]
    genomes=[g[:-4] if g.endswith(".gff") else g for g in genomes]
    for g in ("ecoli_k12","subtilis"):
        if g not in genomes: genomes.append(g)
    summary={}
    for g in genomes:
        r=build_for_genome(a.data_root,g,a.out/g); summary[g]=r["summary"]
    (a.out/"ncrna_index_summary.json").write_text(json.dumps({"genomes":len(summary),"per_genome":summary},indent=1),encoding="utf-8")
    print("genomes",len(summary),"rows",sum(v["n_rows"] for v in summary.values()),"features",sum(v["n_features"] for v in summary.values()))

if __name__=="__main__": raise SystemExit(main())
