#!/usr/bin/env python3
"""Mapea GAF GOA a genomas por nombre de especie (taxdump NCBI)."""
from __future__ import annotations
import csv, gzip, json, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
GAF_DIR=ROOT/"data/labels/external/goa/by_code"
TAX=ROOT/"data/labels/external/goa/taxdump/names.dmp"


def norm(s):
    s=(s or "").lower()
    s=re.sub(r"\b(candidatus|subsp\.?|str\.?|strain)\b"," ",s)
    s=re.sub(r"\([^)]*\)"," ",s)
    return " ".join(re.sub(r"[^a-z0-9 ]+"," ",s).split()[:2])


def main():
    # taxid -> scientific name
    wanted=set()
    gaf_tax={}
    for p in sorted(GAF_DIR.glob("*.gaf.gz")):
        try:
            with gzip.open(p,"rt",errors="replace") as fh:
                for line in fh:
                    if line.startswith("!"): continue
                    f=line.rstrip("\n").split("\t")
                    m=next((x for x in f if x.startswith("taxon:")),None)
                    if m:
                        tid=m.split(":",1)[1].split("|")[0]; gaf_tax[p.name]=tid; wanted.add(tid); break
        except Exception: pass
    names={}
    with TAX.open(encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if "scientific name" not in line: continue
            f=[x.strip() for x in line.split("|")]
            if len(f)>=2 and f[0] in wanted: names[f[0]]=f[1]
    by_species={}
    for fname,tid in gaf_tax.items():
        nm=names.get(tid,"")
        if nm: by_species.setdefault(norm(nm),[]).append((fname,tid,nm))
    rows=[]
    for d in sorted((ROOT/"data/regulation").iterdir()):
        if not d.is_dir() or not (d/"genes.tsv").exists(): continue
        g=d.name; t=ROOT/"data/taxonomy/full"/f"{g}.json"; species=""
        if t.exists():
            try: species=json.loads(t.read_text()).get("species") or ""
            except Exception: pass
        cands=by_species.get(norm(species),[])
        # prefer uniprot over mod
        cands=sorted(cands,key=lambda x:(not x[0].endswith("-uniprot.gaf.gz"),x[0]))
        fname,tid,nm=(cands[0] if cands else ("","",""))
        rows.append({"genome":g,"species":species,"goa_file":fname,
                     "goa_taxid":tid,"goa_species":nm,
                     "status":"ok" if fname else "no_data"})
    out=ROOT/"data/labels/external/goa/goa_map.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"gaf_files":len(gaf_tax),"taxids_named":len(names),
             "genomes":len(rows),"mapped":sum(1 for r in rows if r["status"]=="ok"),
             "gaf_taxids":len(set(gaf_tax.values()))}
    (out.with_suffix(".json")).write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
    return 0

if __name__=="__main__": raise SystemExit(main())
