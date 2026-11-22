#!/usr/bin/env python3
"""Shortlist de ORFs no anotados para estructura (AFDB homólogo o ESMFold exacto)."""
import csv, json, urllib.request
from collections import Counter
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=D/"data"/"metabolism"/"targets"
rows=[r for r in csv.DictReader(open(O/"target_candidates_comparative.tsv"),delimiter="\t")
      if not (r.get("annotation") or "").strip()]
rows.sort(key=lambda r:(-float(r["score"]), r["target"]))
# uno por target, top 10
out=[]; seen=set()
for r in rows:
    if r["target"] in seen: continue
    seen.add(r["target"]); out.append(r)
    if len(out)>=10: break
for r in out:
    acc=r.get("best_acc","")
    r["afdb_entry"]=""; r["afdb_url"]=""; r["structural_route"]="ESMFold/ColabFold del ORF exacto"
    if acc:
        try:
            d=json.load(urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",timeout=10))
            if d:
                r["afdb_entry"]=d[0].get("entryId",""); r["afdb_url"]=d[0].get("pdbUrl","") or f"https://alphafold.ebi.ac.uk/entry/{acc}"
                r["structural_route"]="AFDB del homólogo + ESMFold del ORF exacto"
        except Exception:
            pass
with open(O/"STRUCTURAL_CANDIDATES.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
with open(O/"STRUCTURAL_CANDIDATES.md","w") as fh:
    fh.write("# Candidatos no anotados para estructura\n\n")
    fh.write("AlphaFold DB solo contiene proteínas UniProt conocidas; el ORF exacto de un genoma nuevo no está. Usar AFDB del homólogo + ESMFold/ColabFold del ORF.\n\n")
    fh.write("| target | ORF | dominio | homólogo | pident | score | AFDB | ruta |\n|---|---|---|---|---:|---:|---|---|\n")
    for r in out:
        fh.write(f"| {r['target']} | {r['protein_id']} | {r['domains'] or '—'} | {r['best_acc'] or '—'} | {r['pident']} | {r['score']} | {r['afdb_entry'] or 'no'} | {r['structural_route']} |\n")
    fh.write("\n## Límites\n\n- AFDB del homólogo no valida la función del ORF.\n- Identidad baja (<35%) exige cautela.\n- Sin estructura experimental ni ensayo.\n")
print("structural",len(out),"with_afdb",sum(1 for r in out if r["afdb_entry"]))
