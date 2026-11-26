#!/usr/bin/env python3
"""Añade accesiones B. subtilis (functional hits) al índice GO SQLite."""
from __future__ import annotations
import gzip, json, sqlite3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    data=ROOT/"data"; db=data/"labels/go_mf_cc/goa.sqlite"
    hits=json.loads((ROOT/"results/subtilis/subtilis.functional_annotation.json").read_text()).get("hits",{})
    accs={str(v.get("uniprot","")).split("-")[0] for v in hits.values() if v.get("uniprot")}
    names={}; ns={}; cur=None
    for line in (data/"labels/external/go/go-basic.obo").read_text().splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
    con=sqlite3.connect(db); n=0
    with gzip.open(data/"labels/external/goa/BACSU-uniprot.gaf.gz","rt",errors="replace") as fh:
        batch=[]
        for line in fh:
            if line.startswith("!"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<15 or f[1].split("-")[0] not in accs or f[8] not in ("F","C"): continue
            ev=f[6]; cls="electronic" if ev=="IEA" else ("experimental" if ev in {"EXP","IDA","IPI","IMP","IGI","IEP"} else "curated")
            tax=next((x.split(":",1)[1] for x in f if x.startswith("taxon:")),"")
            batch.append((f[1].split("-")[0],f[4],names.get(f[4],""),f[8],ns.get(f[4],""),ev,cls,"BACSU-uniprot.gaf.gz",tax,f[14]))
        con.executemany("INSERT OR IGNORE INTO go VALUES (?,?,?,?,?,?,?,?,?,?)",batch); n+=len(batch)
    con.commit(); print("inserted",n,"accs",len(accs),"total",con.execute("select count(*) from go").fetchone()[0])
    return 0
if __name__=="__main__": raise SystemExit(main())
