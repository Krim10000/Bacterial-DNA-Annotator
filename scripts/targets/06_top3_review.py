#!/usr/bin/env python3
"""Revisión top 3 por target (in silico, no experimental)."""
import csv
from collections import OrderedDict
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=D/"data"/"metabolism"/"targets"
rows=list(csv.DictReader(open(O/"target_candidates_comparative.tsv"),delimiter="\t"))
by=OrderedDict()
for r in rows: by.setdefault(r["target_id"],[]).append(r)
lines=["# Revisión top 3 por target","","Criterios: dominio + homología + contexto + co-ocurrencia. `annotation=\"\"` = ORF no anotado.",""]
for t,rs in by.items():
    lines.append(f"## {rs[0]['target']} · `{rs[0]['map_id']}`")
    lines.append(f"Ruta: {rs[0]['pathway']}")
    lines.append("")
    lines.append("| # | genoma | ORF | score | evidencia | dominio | homólogo | contexto | co-ocurrencia |")
    lines.append("|---|---|---|---:|---|---|---|---|---|")
    for i,r in enumerate(rs[:3],1):
        lines.append(f"| {i} | {r['genome']} | {r['protein_id']} | {r['score']} | {r['evidence']} | {r['domains'] or '—'} | {r['best_acc'] or '—'} | {r['context_score']} | {r['comparative_specificity']} ({r['comparative_support']}) |")
    lines.append("")
lines += ["## Límites","","- Revisión in silico; no experimental.","- Co-ocurrencia = presencia del candidato en no-gap vs gap.","- No se marca `confirmed`."]
open(O/"TOP3_REVIEW.md","w").write("\n".join(lines)+"\n")
print("targets",len(by),"rows",len(rows))
