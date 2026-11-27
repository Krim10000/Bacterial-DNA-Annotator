#!/usr/bin/env python3
"""Descarga matrices/sitios PRODORIC y los mapea a nuestras especies."""
from __future__ import annotations
import argparse, csv, json, re, time, urllib.request
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"dna-token-bac/1.0"})
    with urllib.request.urlopen(req,timeout=60) as r: return r.read()

def norm(s):
    s=(s or "").lower(); s=re.sub(r"\b(candidatus|subsp\.?|str\.?|strain)\b"," ",s)
    s=re.sub(r"\([^)]*\)"," ",s); return " ".join(re.sub(r"[^a-z0-9 ]+"," ",s).split()[:2])

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/regulatory/external/prodoric")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    base="https://www.prodoric.de/api"
    mats=json.loads(get(base+"/matrix"))
    (a.out/"matrices.json").write_text(json.dumps(mats,indent=1),encoding="utf-8")
    details=[]
    for i,m in enumerate(mats,1):
        mx=m["mx"]; p=a.out/"matrices"/f"MX{mx}.json"
        p.parent.mkdir(exist_ok=True)
        if p.exists(): d=json.loads(p.read_text())
        else:
            d=json.loads(get(f"{base}/matrix/MX{mx}")); p.write_text(json.dumps(d,indent=1),encoding="utf-8")
        details.append(d)
        if i%50==0: print(i,"/",len(mats),flush=True)
        time.sleep(0.02)
    # species map
    manual={"pseudomonas aeruginosa":["paeruginosa"],"bacillus subtilis":["subtilis"],
            "escherichia coli":["ecoli_k12"]}
    by_species=defaultdict(list)
    for d in (a.data_root/"taxonomy/full").glob("*.json"):
        try: sp=json.loads(d.read_text()).get("species","")
        except Exception: continue
        if sp: by_species[norm(sp)].append(d.stem)
    rows=[]
    for d in details:
        org=d.get("organism",{}) or {}; sp=org.get("name","")
        matches=by_species.get(norm(sp), manual.get(norm(sp), []))
        match=";".join(matches) if matches else ""
        for b in d.get("bindingsites",[]) or []:
            rows.append({"source_db":"PRODORIC","source_version":"2021.9",
                         "source_organism":sp,"source_acc":org.get("accession",""),
                         "tf_name":d.get("name",""),"matrix_id":f"MX{d.get('mx','')}",
                         "target_operon":b.get("operon",""),"regulation":b.get("regulation",""),
                         "binding_sequence":b.get("sequence",""),"context":b.get("context",""),
                         "evidence":b.get("evidence",""),"experiments":",".join(b.get("experiments",[]) or []),
                         "genome_match":match,"match_type":"same_species" if match else "no_match",
                         "confidence":"medium" if b.get("experiments") else "low"})
    out=a.out/"prodoric_binding_sites.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["source_db"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"matrices":len(details),"binding_sites":len(rows),"species_matched":len({r['genome_match'] for r in rows if r['genome_match']}),
             "rows_matched":sum(1 for r in rows if r['genome_match']),"out":str(out)}
    (a.out/"prodoric_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8"); print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
