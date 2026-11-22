#!/usr/bin/env python3
"""09_prepare_batch03.py — universo y exclusiones del lote 03.

El universo del lote 03 es `real_gaps.tsv` (3.667 `enzymatic_missing`), no la
revisión histórica de 4.021. Se excluyen, por clave:

  * 1.500 brechas ya cubiertas por los lotes 01+02 (`gap_evidence_all.tsv`);
  * las 77 accionables históricas (`annotation_gap` + `domain_present`);
  * las 343 promovidas a `homology_present` (`review_status=evidence_high`,
    `origin=homology_present`);
  * las 23 `reannotated` (`origin=reannotated`).

Las 343 promovidas ya están entre las 1.500 cubiertas y las 23 reannotated
entre las 77 accionables: la unión es 1.577 claves, no 1.943. Por eso el pool
real es 2.456, no ~1.927; queda declarado en `work_batch_03/universe.json`.

Genera:
  work_batch_03/universe.tsv       (3.667 filas, classification=unknown)
  work_batch_03/exclude_keys.tsv   (unión de exclusiones)
  work_batch_03/universe.json      (conteos y procedencia)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ACTIONABLE, META, REAL_GAPS, REVIEW_QUEUE, ROOT, WORK,
                    read_tsv, write_json, write_tsv)

EVID = META / "evidence"
ALL = EVID / "gap_evidence_all.tsv"
HISTORIC = META / "review" / "gap_review.tsv"
GI = META / "gap_interpretation.tsv"
KEY = ["genome", "gap"]


def _snapshots_from_universe() -> None:
    """Reconstruye los espejos desde el universo congelado (base ya cambiada)."""
    uni = read_tsv(WORK / "universe.tsv")
    real_now = {(r["genome"], r["gap"]): r for r in read_tsv(REAL_GAPS)}
    queue_now = {(r["genome"], r["gap"]): r for r in read_tsv(REVIEW_QUEUE)}
    sel = ({(r["genome"], r["gap"]): r for r in
            read_tsv(WORK / "selected_gaps.tsv")}
           if (WORK / "selected_gaps.tsv").exists() else {})
    write_tsv(WORK / "real_snapshot.tsv", [{
        "genome": r["genome"], "gap": r["gap"],
        "support_400": r["support"], "transport_status": r["transport"],
    } for r in uni], ["genome", "gap", "support_400", "transport_status"])
    qrows = []
    for r in uni:
        k = (r["genome"], r["gap"])
        prio = (queue_now.get(k, {}).get("priority")
                or sel.get(k, {}).get("priority") or "media")
        qrows.append({"genome": k[0], "gap": k[1], "priority": prio,
                      "family": r["family"]})
    write_tsv(WORK / "queue_snapshot.tsv", qrows,
              ["genome", "gap", "priority", "family"])


def keys(rows: list[dict]) -> set[tuple[str, str]]:
    return {(r["genome"], r["gap"]) for r in rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="rehacer el universo desde la base actual")
    a = ap.parse_args()
    if not a.force and (WORK / "universe.tsv").exists():
        if not (WORK / "real_snapshot.tsv").exists() or not (
                WORK / "queue_snapshot.tsv").exists():
            _snapshots_from_universe()
        st = json.loads((WORK / "universe.json").read_text())
        print(f"reutilizando snapshot del lote 03: {st['enzymatic_missing']} "
              f"filas, pool {st['pool']}")
        return 0

    gi = read_tsv(GI)
    real = read_tsv(REAL_GAPS)
    queue = read_tsv(REVIEW_QUEUE)
    allrows = read_tsv(ALL)
    historic = read_tsv(HISTORIC)

    universe = keys(real)
    print(f"base actual: real_gaps {len(universe)} filas")
    assert all(r["origin"] == "enzymatic_missing"
               for r in gi if (r["genome"], r["gap"]) in universe), \
        "real_gaps contiene filas que no son enzymatic_missing"

    covered = keys(allrows)
    actionable = {(r["genome"], r["gap"]) for r in historic
                  if r["classification"] in ACTIONABLE}
    promoted = keys([r for r in gi if r["review_status"] == "evidence_high"])
    reannotated = keys([r for r in gi if r["origin"] == "reannotated"])
    excluded = covered | actionable | promoted | reannotated
    pool = universe - excluded

    prio = {(r["genome"], r["gap"]): r["priority"] for r in queue}
    sup = {(r["genome"], r["gap"]): int(r["support_400"]) for r in real}
    fam = {(r["genome"], r["gap"]): r["family"] for r in queue}
    trans = {(r["genome"], r["gap"]): r["transport_status"] for r in real}

    missing = [k for k in universe if k not in prio or k not in sup or k not in fam]
    assert not missing, f"universo sin prioridad/support/familia: {missing[:3]}"

    # universo para 01_select_worklist.py (schema de gap_review, todo unknown)
    urows = [{
        "genome": g, "gap": p, "family": fam[(g, p)],
        "support": str(sup[(g, p)]), "transport": trans[(g, p)],
        "homolog_status": "not_checked", "domain_status": "not_checked",
        "classification": "unknown", "action": "",
        "evidence": "universo batch_03: enzymatic_missing vigente",
    } for g, p in sorted(universe)]
    write_tsv(WORK / "universe.tsv", urows,
              ["genome", "gap", "family", "support", "transport",
               "homolog_status", "domain_status", "classification", "action",
               "evidence"])

    erows = [{"genome": g, "gap": p} for g, p in sorted(excluded)]
    write_tsv(WORK / "exclude_keys.tsv", erows, ["genome", "gap"])

    # espejo de las tablas base para que un rerun no dependa de cambios externos
    rmap = {(r["genome"], r["gap"]): r for r in real}
    qmap = {(r["genome"], r["gap"]): r for r in queue}
    write_tsv(WORK / "real_snapshot.tsv", [{
        "genome": g, "gap": p,
        "support_400": rmap[(g, p)]["support_400"],
        "transport_status": rmap[(g, p)]["transport_status"],
    } for g, p in sorted(universe)],
        ["genome", "gap", "support_400", "transport_status"])
    write_tsv(WORK / "queue_snapshot.tsv", [{
        "genome": g, "gap": p,
        "priority": qmap[(g, p)]["priority"],
        "family": qmap[(g, p)]["family"],
    } for g, p in sorted(universe)],
        ["genome", "gap", "priority", "family"])

    stats = {
        "batch": "03",
        "snapshot": True,
        "enzymatic_missing": len(universe),
        "covered_total": len(covered),
        "covered_in_universe": len(covered & universe),
        "promoted_total": len(promoted),
        "promoted_in_universe": len(promoted & universe),
        "promoted_in_covered": len(promoted & covered),
        "actionable_total": len(actionable),
        "actionable_in_universe": len(actionable & universe),
        "reannotated_total": len(reannotated),
        "reannotated_in_actionable": len(reannotated & actionable),
        "excluded_union": len(excluded),
        "excluded_in_universe": len(excluded & universe),
        "pool": len(pool),
        "pool_priorities": {
            p: sum(1 for k in pool if prio[k] == p)
            for p in ("alta", "media", "baja")
        },
        "pool_genomes": len({g for g, _ in pool}),
        "note": ("El pool real es 2.456, no ~1.927: las 343 promovidas ya "
                 "estaban en las 1.500 cubiertas y las 23 reannotated entre "
                 "las 77 accionables; la unión de exclusiones es 1.577."),
    }
    write_json(WORK / "universe.json", stats)
    print(f"universo {len(universe)} · excluidas {len(excluded)} "
          f"({len(excluded & universe)} dentro) · pool {len(pool)}")
    print(f"prioridades pool: {stats['pool_priorities']} · "
          f"genomas {stats['pool_genomes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
