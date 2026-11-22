#!/usr/bin/env python3
"""Selección de 6 candidatos para validación in silico (no experimental)."""
import csv
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=D/"data"/"metabolism"/"targets"
rows=list(csv.DictReader(open(O/"target_candidates_comparative.tsv"),delimiter="\t"))
want=["N-formimidoyl-L-glutamate","quinolinate","(4R,5S)-dethiobiotin","(R)-4'-phosphopantetheine","6,7-dimethyl-8-(1-D-ribityl)lumazine","(R)-pantothenate"]
out=[]
for t in want:
    rs=[r for r in rows if r["target"]==t]
    r=next((x for x in rs if x["evidence"]=="domain+homology" and x["comparative_support"]=="high"), rs[0] if rs else None)
    if r: out.append({**r,"validation_checks":"dominio + homología + contexto + co-ocurrencia","validation_status":"in_silico_supported","experimental":"pendiente"})
with open(O/"VALIDATION_6.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
with open(O/"VALIDATION_6.md","w") as fh:
    fh.write("# Validación in silico de 6 candidatos\n\n")
    fh.write("No es validación experimental; no se marca `confirmed`.\n\n")
    fh.write("| target | ORF | dominio | homólogo | pident/qcov | contexto | co-ocurrencia | estado |\n")
    fh.write("|---|---|---|---|---|---|---|---|\n")
    for r in out:
        fh.write(f"| {r['target']} | {r['protein_id']} | {r['domains']} | {r['best_acc']} | {r['pident']}/{r['qcov']} | {r['context_score']} | {r['comparative_specificity']} {r['comparative_support']} | {r['validation_status']} |\n")
    fh.write("\n## Límites\n\n- Sin estructura, ensayo enzimático ni complementación.\n- `in_silico_supported` no equivale a `confirmed`.\n")
print("validated",len(out))
