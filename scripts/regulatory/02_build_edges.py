#!/usr/bin/env python3
"""Construye la tabla común TF–promotor (RegulonDB + PRODORIC)."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/regulatory/evidence_edges.tsv")
    a=ap.parse_args(); a.out.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    # RegulonDB
    rp=a.data_root/"regulondb/TF-RISet.tsv"
    for r in csv.DictReader(rp.open(encoding="utf-8",errors="replace"),delimiter="\t"):
        conf=(r.get("20)confidenceLevel") or "").strip().upper()
        rows.append({"genome":"ecoli_k12","source_db":"RegulonDB","source_version":"local",
                     "source_organism":"Escherichia coli K-12","tf_name":r.get("4)regulatorName",""),
                     "matrix_id":"","target":r.get("19)targetTuOrGene",""),
                     "target_type":"TU_or_gene","regulation":r.get("11)riFunction",""),
                     "binding_sequence":r.get("10)tfrsSeq",""),"evidence":r.get("21)tfrsEvidence",""),
                     "confidence":{"S":"high","C":"medium","W":"low"}.get(conf,conf.lower() or "unknown"),
                     "match_type":"curated_same_species","source_url":"regulondb",
                     "notes":f"riEvidence={r.get('22)riEvidence','')}; sigma={r.get('15)sigmaF','')}"})
    # PRODORIC
    pp=a.data_root/"regulatory/external/prodoric/prodoric_binding_sites.tsv"
    for r in csv.DictReader(pp.open(encoding="utf-8",errors="replace"),delimiter="\t"):
        for g in (r.get("genome_match","") or "").split(";"):
            if not g: continue
            rows.append({"genome":g,"source_db":"PRODORIC","source_version":r.get("source_version",""),
                         "source_organism":r.get("source_organism",""),"tf_name":r.get("tf_name",""),
                         "matrix_id":r.get("matrix_id",""),"target":r.get("target_operon",""),
                         "target_type":"operon","regulation":r.get("regulation",""),
                         "binding_sequence":r.get("binding_sequence",""),"evidence":r.get("experiments",""),
                         "confidence":r.get("confidence","low"),"match_type":r.get("match_type","same_species"),
                         "source_url":"https://www.prodoric.de","notes":f"source_acc={r.get('source_acc','')}; evidence_code={r.get('evidence','')}"})
    fields=["genome","source_db","source_version","source_organism","tf_name","matrix_id","target","target_type","regulation","binding_sequence","evidence","confidence","match_type","source_url","notes"]
    with a.out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"rows":len(rows),"by_source":{}, "by_genome":{}}
    for r in rows:
        summary["by_source"][r["source_db"]]=summary["by_source"].get(r["source_db"],0)+1
        summary["by_genome"][r["genome"]]=summary["by_genome"].get(r["genome"],0)+1
    (a.out.parent/"evidence_edges_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
