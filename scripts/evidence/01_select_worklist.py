#!/usr/bin/env python3
"""01_select_worklist.py — selecciona el lote de brechas `unknown`.

Criterio (declarado, reproducible):
  1. Solo filas `classification == unknown` de la revisión vigente
     (`EVIDENCE_REVIEW`, por defecto `review/gap_review.tsv`).
  2. Prioridad de `REVIEW_QUEUE.tsv`: entran todas las `alta` y se completa el
     lote con las `media` de mayor `support` (desempate genome, gap).
  3. Se excluyen las accionables (`annotation_gap` + `domain_present`) y las
     claves del lote previo (`EVIDENCE_EXCLUDE`, p. ej. las 500 del lote 01).
  4. Tamaño del lote: `EVIDENCE_N` (500 por defecto; 1.000 en el lote 02).

Productores: reacciones de `out_fl/refined/refined_gaps.tsv` (columna
`producers`) y accesiones UniProt de `rhea2uniprot_sprot.tsv` (todas las filas
de cada RHEA productora, sin recortar). No se inventa ningún productor.

Salida: `data/metabolism/evidence/work/selected_gaps.tsv` + `selection.json`.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ACTIONABLE, BATCH, EXCLUDE, KEY, N_BATCH, R2A, REAL_GAPS,
                    REFINED, REVIEW, REVIEW_QUEUE, WORK, clean, read_tsv,
                    write_json, write_tsv)


def main() -> int:
    import os
    if (os.environ.get("EVIDENCE_REUSE_SELECTION") == "1"
            and (WORK / "selected_gaps.tsv").exists()):
        print("reutilizando selección congelada del lote "
              f"({sum(1 for _ in open(WORK / 'selected_gaps.tsv')) - 1} filas)")
        return 0
    review = read_tsv(REVIEW)
    queue = read_tsv(REVIEW_QUEUE)
    real = read_tsv(REAL_GAPS)
    refined = read_tsv(REFINED)
    r2a = read_tsv(R2A)

    # --- validaciones de entrada ------------------------------------------------
    # Las tablas base pueden evolucionar (el lead las actualiza); el lote se
    # construye sobre la intersección de claves que existe en las tres tablas.
    kr = [tuple(r[k] for k in KEY) for r in review]
    kq = [tuple(r[k] for k in KEY) for r in queue]
    kreal = [tuple(r[k] for k in KEY) for r in real]
    assert len(set(kr)) == len(kr), "gap_review con claves duplicadas"
    assert len(set(kq)) == len(kq), "REVIEW_QUEUE con claves duplicadas"
    assert len(set(kreal)) == len(kreal), "real_gaps con claves duplicadas"
    common = set(kr) & set(kq) & set(kreal)
    excluded_base = len(set(kr) - common)

    classes = Counter(r["classification"] for r in review)
    actionable = {tuple(r[k] for k in KEY)
                  for r in review if r["classification"] in ACTIONABLE}
    previous: set[tuple[str, str]] = set()
    if EXCLUDE is not None and Path(EXCLUDE).exists():
        previous = {(r["genome"], r["gap"]) for r in read_tsv(Path(EXCLUDE))}
    print(f"entrada: review={len(review)} queue={len(queue)} real={len(real)} "
          f"comunes={len(common)}; unknown={classes['unknown']} "
          f"accionables={len(actionable)} lote_previo={len(previous)}")

    # --- productores por gap (consistentes entre genomas) ----------------------
    producers: dict[str, set[str]] = defaultdict(set)
    chebi: dict[str, str] = {}
    role: dict[str, str] = {}
    for r in refined:
        g = r["gap"]
        ps = {p.strip() for p in r["producers"].split(";") if p.strip()}
        producers[g] |= ps
        chebi.setdefault(g, r["chebi"])
        role.setdefault(g, r.get("role", ""))
    # un mismo gap no puede tener conjuntos de productores distintos
    by_gap: dict[str, set[frozenset]] = defaultdict(set)
    for r in refined:
        by_gap[r["gap"]].add(frozenset(
            p.strip() for p in r["producers"].split(";") if p.strip()))
    multi = {g: v for g, v in by_gap.items() if len(v) > 1}
    assert not multi, f"productores inconsistentes por gap: {list(multi)[:3]}"

    rid_accs: dict[str, set[str]] = defaultdict(set)
    for r in r2a:
        if r["RHEA_ID"].isdigit() and r["ID"]:
            rid_accs["RHEA:" + r["RHEA_ID"]].add(r["ID"])

    # --- selección ---------------------------------------------------------------
    real_support = {tuple(r[k] for k in KEY): int(r["support_400"]) for r in real}
    priority = {tuple(r[k] for k in KEY): r["priority"] for r in queue}
    by_priority: dict[str, list[dict]] = defaultdict(list)
    for r in review:
        if r["classification"] != "unknown":
            continue
        key = tuple(r[k] for k in KEY)
        if key not in common or key in previous or key in actionable:
            continue
        by_priority[priority[key]].append({
            "genome": r["genome"], "gap": r["gap"], "family": r["family"],
            "support": real_support[key], "priority": priority[key],
            "transport": r["transport"], "role": role.get(r["gap"], ""),
            "chebi": chebi.get(r["gap"], ""),
        })

    alta = sorted(by_priority["alta"],
                  key=lambda r: (-r["support"], r["genome"], r["gap"]))
    media = sorted(by_priority["media"],
                   key=lambda r: (-r["support"], r["genome"], r["gap"]))
    baja = sorted(by_priority["baja"],
                  key=lambda r: (-r["support"], r["genome"], r["gap"]))
    rest = media + baja
    if len(alta) >= N_BATCH:
        selected = alta[:N_BATCH]
    else:
        selected = alta + rest[: N_BATCH - len(alta)]
    available = len(alta) + len(rest)
    assert len(selected) == min(N_BATCH, available), (
        f"selección: {len(selected)} (disponibles {available})")
    assert not ({tuple(r[k] for k in KEY) for r in selected} & previous), \
        "la selección incluye claves del lote previo"
    assert not ({(r["genome"], r["gap"]) for r in selected} & actionable), \
        "la selección incluye filas accionables"

    # --- productores/accesiones ---------------------------------------------------
    no_producer = 0
    no_acc = 0
    for r in selected:
        ps = sorted(producers.get(r["gap"], set()))
        accs: set[str] = set()
        for rid in ps:
            accs |= rid_accs.get(rid, set())
        r["producers"] = ";".join(ps)
        r["producer_accs"] = ";".join(sorted(accs))
        r["n_producers"] = len(ps)
        r["n_accs"] = len(accs)
        no_producer += not ps
        no_acc += not accs

    fields = ["genome", "gap", "chebi", "family", "support", "priority",
              "transport", "role", "producers", "producer_accs",
              "n_producers", "n_accs"]
    rows = [{k: clean(r[k]) for k in fields} for r in selected]
    out = WORK / "selected_gaps.tsv"
    write_tsv(out, rows, fields)

    sup = sorted(int(r["support"]) for r in rows)
    summary = {
        "batch": BATCH,
        "batch_size": len(rows),
        "requested_size": N_BATCH,
        "unknown_total": classes["unknown"],
        "excluded_previous_batch": len(previous),
        "excluded_actionable": len(actionable),
        "priorities": dict(Counter(r["priority"] for r in rows)),
        "unique_metabolites": len({r["gap"] for r in rows}),
        "genomes": len({r["genome"] for r in rows}),
        "families": len({r["family"] for r in rows}),
        "support_min": sup[0],
        "support_max": sup[-1],
        "support_median": sup[len(sup) // 2],
        "rows_without_producer": no_producer,
        "rows_without_accs": no_acc,
        "unique_producer_rhea": len({p for r in rows for p in r["producers"].split(";") if p}),
        "unique_producer_accs": len({a for r in rows for a in r["producer_accs"].split(";") if a}),
        "input_rows": {"review": len(review), "queue": len(queue),
                       "real": len(real), "common_keys": len(common),
                       "excluded_missing_base": excluded_base},
        "input": {
            "gap_review": str(REVIEW.relative_to(REVIEW.parents[2])),
            "review_queue": str(REVIEW_QUEUE.relative_to(REVIEW_QUEUE.parents[2])),
            "real_gaps": str(REAL_GAPS.relative_to(REAL_GAPS.parents[2])),
            "refined_gaps": str(REFINED.relative_to(REFINED.parents[3])),
            "rhea2uniprot_sprot": str(R2A.relative_to(R2A.parents[2])),
        },
        "rule": (f"todas las unknown alta + media de mayor support hasta "
                 f"{N_BATCH}; accionables y lote previo excluidos; "
                 "productores Rhea completos"),
    }
    write_json(WORK / "selection.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"\n-> {out.relative_to(WORK.parents[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
