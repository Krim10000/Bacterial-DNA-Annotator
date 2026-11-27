#!/usr/bin/env python3
"""Agrega las predicciones P2Rank de p2rank_lowcpu en una tabla."""
from __future__ import annotations
import csv, json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    root=ROOT/"data/sigma/structures/p2rank_lowcpu"
    q={r["accession"]:r for r in csv.DictReader(open(root/"queue.tsv"),delimiter="\t")}
    rows=[]; withp=set(); no=set()
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        preds=list(d.glob("*_predictions.csv"))
        acc=d.name
        if not preds:
            no.add(acc); continue
        withp.add(acc)
        for line in preds[0].read_text(errors="replace").splitlines():
            if not line.strip() or line.lower().startswith("name"): continue
            f=[x.strip() for x in line.split(",")]
            if len(f)<10: continue
            rows.append({"accession":acc,"origin":q.get(acc,{}).get("origin",""),
                         "pdb":q.get(acc,{}).get("pdb",""),"pocket_id":f[0],
                         "pocket_rank":f[1],"score":f[2],"probability":f[3],
                         "sas_points":f[4],"surf_atoms":f[5],"residues":f[9]})
    out=root/"p2rank_lowcpu_pockets.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["accession"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"structures_done":len(withp|no),"with_pocket":len(withp),"without_pocket":len(no),
             "pocket_rows":len(rows),"origins":dict(Counter(r["origin"] for r in rows)),"out":str(out)}
    (root/"p2rank_lowcpu_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
