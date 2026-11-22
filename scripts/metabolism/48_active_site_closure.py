#!/usr/bin/env python3
"""Readiness: qué gaps tienen sitio activo + candidato + estructura para cerrar."""
import csv, json
from collections import Counter, defaultdict
from pathlib import Path
D=Path(__file__).resolve().parents[2]; E=D/"data"/"metabolism"/"essential"
pri=list(csv.DictReader(open(E/"essential_priority_targets.tsv"),delimiter="\t"))
asite=defaultdict(int)
for r in csv.DictReader(open(E/"ACTIVE_SITES_ESSENTIAL.tsv"),delimiter="\t"):
    if r["feature"]!="none": asite[(r["map_id"],r["target"])]+=1
cands={}
for p in [D/"data/metabolism/targets/target_candidates_comparative.tsv", D/"data/metabolism/targets_p3/target_candidates.tsv"]:
    if not p.exists(): continue
    for r in csv.DictReader(open(p),delimiter="\t"):
        k=(r["map_id"],r["target"])
        if k not in cands or float(r["score"])>float(cands[k]["score"]): cands[k]=r
struct=set()
sp=D/"data/metabolism/targets/STRUCTURAL_CANDIDATES.tsv"
if sp.exists():
    for r in csv.DictReader(open(sp),delimiter="\t"): struct.add(r["target"])
rows=[]
for r in pri:
    k=(r["map_id"],r["target"]); c=cands.get(k); a=asite.get(k,0)
    if a and c and r["target"] in struct: st="structural_ready"
    elif a and c: st="needs_structure"
    elif c: st="needs_active_site"
    elif a: st="needs_candidate"
    else: st="no_reference_or_candidate"
    rows.append({"map_id":r["map_id"],"pathway":r["pathway"],"target":r["target"],"priority":r["priority"],
                 "active_site_features":a,"candidate":c["protein_id"] if c else "",
                 "candidate_score":c["score"] if c else "","candidate_evidence":c["evidence"] if c else "",
                 "best_acc":c["best_acc"] if c else "","readiness":st,"exception":r.get("known_exception","")})
with open(E/"ACTIVE_SITE_CLOSURE_READINESS.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
s=Counter(r["readiness"] for r in rows)
with open(E/"ACTIVE_SITE_CLOSURE_READINESS.md","w") as fh:
    fh.write("# Readiness de cierre por sitio activo\n\n")
    fh.write("No cierra gaps: indica qué falta para poder validar estructuralmente.\n\n")
    for k,v in s.most_common(): fh.write(f"- `{k}`: {v}\n")
    fh.write("\n| target | prioridad | features | candidato | evidencia | readiness |\n|---|---|---:|---|---|---|\n")
    for r in sorted(rows,key=lambda x:(-x["active_site_features"],x["priority"]))[:30]:
        fh.write(f"| {r['target']} | {r['priority']} | {r['active_site_features']} | {r['candidate'] or '—'} | {r['candidate_evidence'] or '—'} | {r['readiness']} |\n")
    fh.write("\n## Regla\n\n1. sitio activo del homólogo;\n2. ORF candidato;\n3. ESMFold;\n4. alineación;\n5. residuos conservados;\n6. bolsillo/docking;\n7. experimento.\n")
json.dump(dict(s),open(E/"ACTIVE_SITE_CLOSURE_READINESS.json","w"),indent=2)
print(dict(s))
