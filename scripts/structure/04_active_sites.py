#!/usr/bin/env python3
"""Residuos catalíticos/de unión del homólogo UniProt (no del ORF exacto)."""
import csv, json, urllib.request
from collections import defaultdict
from pathlib import Path
D=Path(__file__).resolve().parents[2]; T=D/"data"/"metabolism"/"targets"
rows=list(csv.DictReader(open(T/"STRUCTURAL_CANDIDATES.tsv"),delimiter="\t"))
cache={}
out=[]
for r in rows:
    acc=r.get("best_acc","").strip()
    if not acc: continue
    if acc not in cache:
        try: cache[acc]=json.load(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json",timeout=20))
        except Exception: cache[acc]={}
    d=cache[acc]
    name=(d.get("proteinDescription",{}).get("recommendedName",{}) or {}).get("fullName",{}).get("value","")
    ecs=d.get("proteinDescription",{}).get("recommendedName",{}) or {}
    ec=";".join(x.get("value","") for x in ecs.get("ecNumbers",[]) or [])
    for f in d.get("features",[]) or []:
        if f.get("type") not in ("Active site","Binding site","Metal binding","Site"): continue
        loc=f.get("location",{}) or {}; s=(loc.get("start",{}) or {}).get("value",""); e=(loc.get("end",{}) or {}).get("value","")
        out.append({"target":r["target"],"orf":r["protein_id"],"acc":acc,"protein":name,"ec":ec,
                    "organism":d.get("organism",{}).get("scientificName",""),
                    "feature":f.get("type"),"description":f.get("description",""),"start":s,"end":e,
                    "note":"del homólogo; transferir por alineación"})
with open(T/"ACTIVE_SITES.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
by=defaultdict(list)
for r in out: by[(r["target"],r["acc"])].append(r)
with open(T/"ACTIVE_SITES.md","w") as fh:
    fh.write("# Residuos catalíticos de homólogos\n\n")
    fh.write("Son residuos del **homólogo UniProt**, no del ORF exacto. Transferir solo tras alineación estructural.\n\n")
    for (t,acc),rs in by.items():
        r0=rs[0]; fh.write(f"## {t} · {acc}\n\n{r0['protein']} ({r0['organism']}) · EC {r0['ec'] or '—'}\n\n")
        fh.write("| tipo | descripción | inicio | fin |\n|---|---|---:|---:|\n")
        for x in rs[:12]: fh.write(f"| {x['feature']} | {x['description']} | {x['start']} | {x['end']} |\n")
        fh.write("\n")
    fh.write("## Límites\n\n- UniProt no siempre anota sitio activo.\n- El residuo del homólogo no equivale al del ORF.\n- Sin alineación + conservación + experimento, no es prueba.\n")
print("homologs",len(by),"features",len(out))
