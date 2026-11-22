#!/usr/bin/env python3
"""11_prepare_batch04.py — universo y exclusiones del lote 04 (cierre).

Universo: `real_gaps.tsv` (3.433 `enzymatic_missing`). Exclusiones por clave:

  * 2.500 brechas ya cubiertas por los lotes 01-03 (`gap_evidence_all.tsv`);
  * 77 accionables históricas (`annotation_gap` + `domain_present`);
  * 23 `reannotated` (`origin=reannotated`);
  * las promovidas a `homology_present` (`origin=homology_present`), que
    incluyen las 343 + 213 de contexto y otras ya cerradas.

El pool resultante es 1.456, no 1.510: las 54 accionables que siguen dentro
del universo explican la diferencia. Se deja congelado en `universe.json` y en
los espejos `real_snapshot.tsv` / `queue_snapshot.tsv` para que el rerun no
dependa de cambios posteriores de la base.

Genera:
  work_batch_04/universe.tsv
  work_batch_04/exclude_keys.tsv
  work_batch_04/real_snapshot.tsv
  work_batch_04/queue_snapshot.tsv
  work_batch_04/universe.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ACTIONABLE, META, REAL_GAPS, REVIEW_QUEUE, WORK, read_tsv, write_json, write_tsv

EVID = META / "evidence"
ALL = EVID / "gap_evidence_all.tsv"
HISTORIC = META / "review" / "gap_review.tsv"
GI = META / "gap_interpretation.tsv"


def keys(rows: list[dict]) -> set[tuple[str, str]]:
    return {(r["genome"], r["gap"]) for r in rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="rehacer el universo desde la base actual")
    a = ap.parse_args()
    if (not a.force and (WORK / "universe.json").exists()
            and (WORK / "universe.tsv").exists()
            and (WORK / "real_snapshot.tsv").exists()
            and (WORK / "queue_snapshot.tsv").exists()):
        st = json.loads((WORK / "universe.json").read_text())
        print(f"reutilizando snapshot del lote 04: "
              f"{st['enzymatic_missing']} filas, pool {st['pool']}")
        return 0
    gi = read_tsv(GI)
    # base real, no los espejos que este mismo script va a crear
    real = read_tsv(META / "real_gaps.tsv")
    queue = read_tsv(META / "REVIEW_QUEUE.tsv")
    allrows = read_tsv(ALL)
    historic = read_tsv(HISTORIC)

    universe = keys(real)
    covered = keys(allrows)
    actionable = {(r["genome"], r["gap"]) for r in historic
                  if r["classification"] in ACTIONABLE}
    reannotated = keys([r for r in gi if r["origin"] == "reannotated"])
    promoted = keys([r for r in gi if r["origin"] == "homology_present"])
    excluded = covered | actionable | reannotated | promoted
    pool = universe - excluded

    prio = {(r["genome"], r["gap"]): r["priority"] for r in queue}
    sup = {(r["genome"], r["gap"]): int(r["support_400"]) for r in real}
    fam = {(r["genome"], r["gap"]): r["family"] for r in queue}
    trans = {(r["genome"], r["gap"]): r["transport_status"] for r in real}
    missing = [k for k in universe if k not in prio or k not in sup or k not in fam]
    assert not missing, f"universo sin prioridad/support/familia: {missing[:3]}"

    write_tsv(WORK / "universe.tsv", [{
        "genome": g, "gap": p, "family": fam[(g, p)],
        "support": str(sup[(g, p)]), "transport": trans[(g, p)],
        "homolog_status": "not_checked", "domain_status": "not_checked",
        "classification": "unknown", "action": "",
        "evidence": "universo batch_04: enzymatic_missing vigente",
    } for g, p in sorted(universe)],
        ["genome", "gap", "family", "support", "transport", "homolog_status",
         "domain_status", "classification", "action", "evidence"])

    write_tsv(WORK / "exclude_keys.tsv",
              [{"genome": g, "gap": p} for g, p in sorted(excluded)],
              ["genome", "gap"])
    write_tsv(WORK / "real_snapshot.tsv", [{
        "genome": g, "gap": p, "support_400": str(sup[(g, p)]),
        "transport_status": trans[(g, p)],
    } for g, p in sorted(universe)],
        ["genome", "gap", "support_400", "transport_status"])
    write_tsv(WORK / "queue_snapshot.tsv", [{
        "genome": g, "gap": p, "priority": prio[(g, p)],
        "family": fam[(g, p)],
    } for g, p in sorted(universe)],
        ["genome", "gap", "priority", "family"])

    stats = {
        "batch": "04",
        "snapshot": True,
        "enzymatic_missing": len(universe),
        "covered_total": len(covered),
        "covered_in_universe": len(covered & universe),
        "promoted_total": len(promoted),
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
        "note": ("Pool real 1.456, no 1.510: 54 accionables siguen dentro del "
                 "universo. Las promovidas y reannotated ya están excluidas por "
                 "las otras listas."),
    }
    write_json(WORK / "universe.json", stats)
    print(f"universo {len(universe)} · excluidas {len(excluded)} "
          f"({len(excluded & universe)} dentro) · pool {len(pool)}")
    print(f"prioridades pool: {stats['pool_priorities']} · "
          f"genomas {stats['pool_genomes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
