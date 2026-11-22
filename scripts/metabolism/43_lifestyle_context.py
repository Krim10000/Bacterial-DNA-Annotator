#!/usr/bin/env python3
"""Contexto tamaño/taxón para interpretar huecos: no es una llamada de estilo de vida."""
import csv, json
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parents[2]
meta = {r["genome"]: r for r in csv.DictReader(open(D/"data/taxonomy/full/genome_meta.tsv"), delimiter="\t")}
asm = {}
with open(D/"data/assembly_summary.txt", errors="ignore") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 25: asm[f[0]] = f
HOST_GENERA = {
    "chlamydia","chlamydophila","mycoplasma","ureaplasma","phytoplasma","spiroplasma",
    "rickettsia","orientia","ehrlichia","anaplasma","neorickettsia","wolbachia",
    "buchnera","blochmannia","wigglesworthia","sulcia","carsonella","tremblaya",
    "nasuia","hodgkinia","profftia","zinderia","kinetoplastibacterium","borrelia",
}
HOST_WORDS = ("endosymbiont","symbiont","endoparasite","intracellular","parasite")
rows=[]
for g,m in meta.items():
    acc=m.get("accession") or g
    size=int(float(asm.get(acc,[None]*26)[25])) if asm.get(acc) and asm[acc][25].isdigit() else None
    name=(m.get("organism_name") or "").lower()
    genus=(m.get("genus") or "").lower()
    kw=[w for w in HOST_WORDS if w in name]
    gen=genus in HOST_GENERA or any(x in name for x in ("candidatus","endosymbiont"))
    signals=[]
    if gen: signals.append("taxon_huésped")
    if kw: signals.append("nombre_huésped")
    if size is not None and size < 1_500_000: signals.append("genoma_pequeño")
    if size is not None and size >= 3_000_000: signals.append("genoma_grande")
    if gen or kw:
        exp="host_associated"; conf="high" if size and size<2_000_000 else "medium"
        msg="Taxón/nombre compatible con endosimbionte o endoparásito; se espera reducción de rutas. Un hueco no es automáticamente un error."
    elif size is not None and size < 1_500_000:
        exp="host_associated"; conf="medium"
        msg="Genoma pequeño (<1,5 Mb), compatible con vida dependiente de huésped; se espera reducción de rutas."
    elif size is not None and size >= 3_000_000:
        exp="free_living"; conf="medium"
        msg="Genoma grande sin señal de huésped; se esperan rutas core completas. Los huecos en core son prioritarios."
    else:
        exp="uncertain"; conf="low"
        msg="Sin señal concluyente de tamaño/taxón; interpretar huecos con cautela."
    rows.append({"genome":g,"organism":m.get("organism_name",""),"phylum":m.get("phylum",""),
                 "size_mb":round(size/1e6,3) if size else "",
                 "expectation":exp,"confidence":conf,"signals":";".join(signals),
                 "message":msg})
with open(D/"data/metabolism/essential/lifestyle_context.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
s={"rows":len(rows),"expectation":dict(Counter(r["expectation"] for r in rows)),
   "confidence":dict(Counter(r["confidence"] for r in rows)),
   "host_genera_rows":sum(1 for r in rows if "taxon_huésped" in r["signals"])}
json.dump(s,open(D/"data/metabolism/essential/lifestyle_summary.json","w"),indent=2)
print(json.dumps(s,indent=2))
