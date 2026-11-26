#!/usr/bin/env python3
"""Compara etiquetas GO con y sin evidencia IEA."""
from __future__ import annotations
import argparse, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from webapp.labels_engine import read_obo, read_gaf, ancestors, GO_PATTERNS


def labels_for(gaf, names, parents):
    cache={}; out=defaultdict(set)
    for acc, gos in gaf.items():
        for go in gos:
            for g in ancestors(go, parents, cache):
                nm=names.get(g,"").lower()
                for lab,pat in GO_PATTERNS.items():
                    if nm and re.search(pat,nm): out[acc].add(lab)
    return out


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ext", type=Path, default=ROOT/"data/labels/external")
    ap.add_argument("--out", type=Path, default=ROOT/"reports/labels_pilot")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    names,parents=read_obo(a.ext/"go/go-basic.obo")
    files={"ecoli_k12":"ECOLI-uniprot.gaf.gz","subtilis":"BACSU-uniprot.gaf.gz"}
    out={}
    for genome,f in files.items():
        p=a.ext/"goa"/f
        g0=read_gaf(p,allow_iea=False); g1=read_gaf(p,allow_iea=True)
        l0=labels_for(g0,names,parents); l1=labels_for(g1,names,parents)
        added_genes={k for k,v in l1.items() if v - l0.get(k,set())}
        only=Counter(l for k,v in l1.items() for l in (v-l0.get(k,set())))
        out[genome]={"gaf_without_iea":len(g0),"gaf_with_iea":len(g1),
                     "genes_with_new_label":len(added_genes),
                     "new_labels":dict(only),
                     "label_genes_without":{k:sum(1 for v in l0.values() if k in v) for k in GO_PATTERNS},
                     "label_genes_with":{k:sum(1 for v in l1.values() if k in v) for k in GO_PATTERNS}}
    (a.out/"go_iea_diff.json").write_text(json.dumps(out,indent=1),encoding="utf-8")
    print(json.dumps(out,indent=1))
    return 0

if __name__=="__main__": raise SystemExit(main())
