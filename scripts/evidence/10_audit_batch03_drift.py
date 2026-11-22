#!/usr/bin/env python3
"""10_audit_batch03_drift.py — deriva de la base después de congelar el lote 03.

El lote 03 se construyó sobre 3.667 `enzymatic_missing` (snapshot). Al cerrar,
la base ya tenía 3.433: 234 brechas pasaron a `homology_present`/`spontaneous`
durante la corrida. Este script mide el efecto y lo deja escrito en
`summary_batch_03.json` (`post_close_audit`) y en `README_batch_03.md`, sin
tocar la tabla de evidencia ni las tablas base.

No es una validación de evidencia: es trazabilidad del universo.

Uso:
    EVIDENCE_BATCH=03 python3 10_audit_batch03_drift.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BATCH, EVID, META, REAL_GAPS, README_OUT, SUMMARY_OUT,
                    WORK,
                    read_tsv)

ALL = META / "evidence" / "gap_evidence_all.tsv"
TABLE = EVID / f"gap_evidence_batch_{BATCH}.tsv"


def keys(rows: list[dict]) -> set[tuple[str, str]]:
    return {(r["genome"], r["gap"]) for r in rows}


def main() -> int:
    # base actual, no el snapshot: el objetivo es medir la deriva
    real = read_tsv(META / "real_gaps.tsv")
    gi = read_tsv(META / "gap_interpretation.tsv")
    allrows = read_tsv(ALL)
    table = read_tsv(TABLE)
    universe = json.loads((WORK / "universe.json").read_text())

    now = keys(real)
    b03 = keys(table)
    promoted = {(r["genome"], r["gap"]) for r in gi
                if r["origin"] == "homology_present"}
    spontaneous = {(r["genome"], r["gap"]) for r in gi
                   if r["origin"] == "spontaneous"}
    covered = keys(allrows) | b03
    audit = {
        "snapshot_enzymatic_missing": universe["enzymatic_missing"],
        "current_enzymatic_missing": len(now),
        "delta": universe["enzymatic_missing"] - len(now),
        "batch_rows": len(b03),
        "batch_still_unknown": len(b03 & now),
        "batch_promoted": len(b03 & promoted),
        "batch_spontaneous": len(b03 & spontaneous),
        "batch_other_removed": len(b03 - now - promoted - spontaneous),
        "current_unknown_covered_by_batches": len(covered & now),
        "current_unknown_remaining": len(now - covered),
        "note": (f"snapshot {universe['enzymatic_missing']} vs base actual "
                 f"{len(now)} (Δ {len(now) - universe['enzymatic_missing']}); "
                 f"el lote {BATCH} queda congelado sobre su snapshot y se "
                 "documenta qué filas ya fueron promovidas"),
    }
    if BATCH == "04":
        prior = 0
        for name in ("summary.json", "summary_batch_02.json",
                     "summary_batch_03.json"):
            path = EVID / name
            if path.exists():
                prior += json.loads(path.read_text())["evidence_counts"].get(
                    "not_checked_homolog", 0)
        batch_nbr = json.loads(SUMMARY_OUT.read_text()).get(
            "no_bacterial_reference", 0)
        audit["no_bacterial_reference_prior"] = prior
        audit["no_bacterial_reference_batch"] = batch_nbr
        audit["no_bacterial_reference_total"] = prior + batch_nbr
    (WORK / "drift_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n")

    summary = json.loads(SUMMARY_OUT.read_text())
    summary["post_close_audit"] = audit
    SUMMARY_OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    readme = README_OUT.read_text()
    marker = "## Deriva post-cierre"
    if marker in readme:
        readme = readme.split(marker)[0].rstrip() + "\n"
    if marker not in readme:
        block = "\n".join([
            "",
            marker,
            "",
            f"- Snapshot del lote: {audit['snapshot_enzymatic_missing']} "
            f"`enzymatic_missing`; base al cerrar: "
            f"{audit['current_enzymatic_missing']} (Δ {audit['delta']}).",
            f"- De las {audit['batch_rows']} filas del lote, "
            f"{audit['batch_still_unknown']} siguen `unknown` y "
            f"{audit['batch_promoted']} ya pasaron a `homology_present` "
            "durante la corrida.",
            f"- Cobertura sobre la base actual: "
            f"{audit['current_unknown_covered_by_batches']} unknown cubiertos; "
            f"quedan {audit['current_unknown_remaining']}.",
            *([f"- `no_bacterial_reference`: "
               f"{audit['no_bacterial_reference_batch']} nuevas en el lote 04; "
               f"{audit['no_bacterial_reference_prior']} acumuladas en 01-03; "
               f"total {audit['no_bacterial_reference_total']}."]
              if "no_bacterial_reference_total" in audit else []),
            "",
        ])
        README_OUT.write_text(readme.rstrip() + "\n" + block)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
