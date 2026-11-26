#!/usr/bin/env python3
"""Deriva un sigmulón σ→TU de B. subtilis desde TSS gold + operones.

No es un gold independiente: usa la misma anotación BSGatlas/DBTBS. Se marca
como `derived_bsu_sigmulon` y se documenta.
"""
from __future__ import annotations
import argparse, csv, json
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read(p):
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"reports/labels_pilot/bsu_derived_sigmulon.tsv")
    ap.add_argument("--max-up",type=int,default=1500)
    a=ap.parse_args()
    gold=read(a.data_root/"sigma/gold/subtilis/tss_sigma.tsv")
    ops=read(a.data_root/"dbtbs/gt_operons.tsv")
    by_strand=defaultdict(list)
    for op in ops:
        try: by_strand[op.get("strand","")].append((int(op["start"]),int(op["end"]),op))
        except Exception: pass
    rows=[]; unmapped=0
    for t in gold:
        try: pos=int(t["pos"])
        except Exception: unmapped+=1; continue
        cand=[]
        for s,e,op in by_strand.get(t.get("strand",""),()):
            if s <= pos <= e or 0 <= s-pos <= a.max_up:
                cand.append((abs(s-pos),op))
        if not cand:
            unmapped+=1; continue
        op=sorted(cand,key=lambda x:x[0])[0][1]
        rows.append({"sigma":t.get("sigma",""),"operon_id":op.get("operon_id",""),
                     "operon_name":op.get("operon_name",""),"strand":t.get("strand",""),
                     "tss":pos,"operon_start":op.get("start",""),
                     "operon_end":op.get("end",""),"source":"derived_bsu_sigmulon",
                     "note":"TSS gold + BSGatlas operon; no independente"})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"n_tss":len(gold),"n_mapped":len(rows),"n_unmapped":unmapped,
             "by_sigma":dict(sorted((k,sum(1 for r in rows if r["sigma"]==k)) for k in {r["sigma"] for r in rows})),
             "n_operons":len({r["operon_id"] for r in rows}),
             "status":"derived_not_independent"}
    (a.out.with_suffix(".json")).write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
    return 0

if __name__=="__main__": raise SystemExit(main())
