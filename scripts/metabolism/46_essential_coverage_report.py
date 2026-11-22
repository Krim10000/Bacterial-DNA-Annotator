#!/usr/bin/env python3
"""Informe de cobertura esencial y objetivos prioritarios, con excepciones."""
import csv, json, statistics
from collections import Counter, defaultdict
from pathlib import Path

D=Path(__file__).resolve().parents[2]; E=D/"data"/"metabolism"/"essential"
paths=list(csv.DictReader(open(E/"essential_pathways.tsv"),delimiter="\t"))
comp={r["map_id"]:r for r in csv.DictReader(open(E/"essential_completeness.tsv"),delimiter="\t")}
pc=defaultdict(list)
for r in csv.DictReader(open(D/"data/metabolism/pathway_completeness.tsv"),delimiter="\t"):
    pc[r["pathway"].split(" ")[0]].append(r)
q=list(csv.DictReader(open(E/"essential_gap_queue.tsv"),delimiter="\t"))
exc=defaultdict(list)
for r in csv.DictReader(open(E/"gap_exception_flags.tsv"),delimiter="\t"):
    exc[(r["genome"],r["gap"])].append(r)
ref={}
for r in csv.DictReader(open(D/"out_fl/refined/refined_gaps.tsv"),delimiter="\t"):
    ref.setdefault((r["genome"],r["gap"]),r)
ec=defaultdict(set)
for line in open(D/"data/reactions/tsv/rhea2ec.tsv"):
    f=line.rstrip().split("\t")
    if len(f)>=4 and f[0]!="RHEA_ID": ec["RHEA:"+f[0]].add(f[3])
# por ruta
by_path=defaultdict(list)
for r in q: by_path[r["map_id"]].append(r)
path_rows=[]
for p in paths:
    mid=p["map_id"]; rs=by_path.get(mid,[])
    cov=[float(x["coverage"]) for x in pc.get(mid,[])] or [0]
    uniq={(r["gap"]) for r in rs}
    uexc={r["gap"] for r in rs if exc.get((r["genome"],r["gap"]))}
    path_rows.append({**p,"genomes":len(pc.get(mid,[])),
        "mean_coverage":round(sum(cov)/len(cov),4),"median_coverage":round(statistics.median(cov),4),
        "p10_coverage":round(sorted(cov)[max(0,int(len(cov)*0.1)-1)],4),
        "gap_rows":len(rs),"unique_targets":len(uniq),"targets_with_exception":len(uexc),
        "targets_without_exception":len(uniq-uexc)})
with open(E/"essential_coverage_by_pathway.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(path_rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(path_rows)
# objetivos únicos
t=defaultdict(lambda: {"genomes":0,"supports":[],"exception":set(),"examples":set(),"transport":set(),"class":set(),"rhea":set(),"ec":set()})
for r in q:
    x=t[(r["map_id"],r["gap"])]; x["genomes"]+=1
    try: x["supports"].append(int(r["support"]))
    except: pass
    x["transport"].add(r["transport"]); x["class"].add(r["classification"])
    for e in exc.get((r["genome"],r["gap"]),[]):
        x["exception"].add(e["exception"]); x["examples"].add(e["example_taxon"])
    rr=ref.get((r["genome"],r["gap"]))
    if rr:
        for rid in rr.get("producers","").split(";"):
            if rid:
                x["rhea"].add(rid); x["ec"] |= ec.get(rid,set())
rows=[]
for (mid,gap),x in t.items():
    p=next((p for p in paths if p["map_id"]==mid),{})
    mx=max(x["supports"] or [0]); mean=round(sum(x["supports"])/len(x["supports"]),1) if x["supports"] else 0
    has_exc=bool(x["exception"]); core=p.get("tier")=="core"
    if core and not has_exc and mx>=100 and "no_transport" in x["transport"]: pri="P1"
    elif core and not has_exc and mx>=30: pri="P2"
    elif core and has_exc: pri="P3_exception"
    elif not core and not has_exc: pri="P4_conditional"
    else: pri="P5_conditional_exception"
    rows.append({"map_id":mid,"pathway":p.get("pathway",mid),"tier":p.get("tier",""),
                 "target":gap,"genomes":x["genomes"],"support_max":mx,"support_mean":mean,
                 "transport":";".join(sorted(x["transport"])),"classification":";".join(sorted(x["class"])),
                 "rhea":";".join(sorted(x["rhea"])),"ec":";".join(sorted(x["ec"])) or "",
                 "known_exception":";".join(sorted(x["exception"])),"example_taxon":";".join(sorted(x["examples"])),
                 "priority":pri})
order={"P1":1,"P2":2,"P3_exception":3,"P4_conditional":4,"P5_conditional_exception":5}
rows.sort(key=lambda r:(order[r["priority"]],-r["support_max"],r["target"]))
with open(E/"essential_priority_targets.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
# markdown
lines=["# Cobertura de rutas esenciales y objetivos prioritarios","",
       f"- Rutas: **{len(paths)}** (core {sum(1 for p in paths if p['tier']=='core')}, conditional {sum(1 for p in paths if p['tier']=='conditional')}).",
       f"- Filas genoma–gap: **{len(q)}**; objetivos únicos: **{len(rows)}**.",
       f"- Objetivos con excepción documentada: **{sum(1 for r in rows if r['known_exception'])}**.",
       "- Prioridad: P1 core, sin excepción, soporte ≥100, sin transporte; P2 core; P3 core con excepción; P4/P5 conditional.","",
       "## Cobertura por ruta","",
       "| ruta | tier | cobertura media | p10 | gaps | objetivos | con excepción |","|---|---|---:|---:|---:|---:|---:|"]
for r in sorted(path_rows,key=lambda x:(x["tier"]!="core",x["mean_coverage"])):
    lines.append(f"| {r['pathway']} | {r['tier']} | {r['mean_coverage']} | {r['p10_coverage']} | {r['gap_rows']} | {r['unique_targets']} | {r['targets_with_exception']} |")
lines += ["","## Objetivos prioritarios P1/P2 (buscar enzima)","",
          "| objetivo | ruta | genomas | soporte max | Rhea | EC | excepción |","|---|---|---:|---:|---|---|---|"]
for r in [x for x in rows if x["priority"] in ("P1","P2")][:40]:
    lines.append(f"| {r['target']} | {r['pathway']} | {r['genomes']} | {r['support_max']} | {r['rhea']} | {r['ec'] or '—'} | {r['known_exception'] or '—'} |")
lines += ["","## Objetivos con excepción documentada (no perseguir sin revisar)","",
          "| objetivo | ruta | genomas | excepción | ejemplo |","|---|---|---:|---|---|"]
for r in [x for x in rows if x["known_exception"]][:30]:
    lines.append(f"| {r['target']} | {r['pathway']} | {r['genomes']} | {r['known_exception']} | {r['example_taxon']} |")
lines += ["","## Límites","",
          "- Cobertura es Rhea/KEGG, no validación biológica.",
          "- Una excepción documentada es un ejemplo, no una exención automática.",
          "- `known_exception` se asigna por ruta, no por especie; revisar `lifestyle_context.tsv`.",
          "- No se marcó nada como `confirmed`."]
open(E/"COVERAGE_REPORT.md","w").write("\n".join(lines)+"\n")
print({"pathways":len(path_rows),"targets":len(rows),"P1":sum(1 for r in rows if r['priority']=='P1'),"P2":sum(1 for r in rows if r['priority']=='P2'),"exception_targets":sum(1 for r in rows if r['known_exception'])})
