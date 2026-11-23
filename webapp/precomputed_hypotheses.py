"""Carga hipótesis estructurales precalculadas desde TSV."""
from __future__ import annotations

import csv
from pathlib import Path

from structural_hypotheses import build

DEFAULT = Path("data/metabolism/targets/UNCHAR_AFDB_PRESENT_filters.tsv")


def _b(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "si", "sí")


def load(path: str | Path | None = None, genome: str = "") -> list[dict]:
    p = Path(path) if path else DEFAULT
    if not p.exists():
        return []
    out = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        if genome and r.get("genome") != genome:
            continue
        row = dict(r)
        for k in ("f1_foldseek", "f2_enzymm", "f3_pocket", "f4_alphafill", "f5_active_site"):
            row[k] = _b(row.get(k))
        meta = {
            "genome": r.get("genome", ""),
            "protein_id": r.get("uniprot", "") or r.get("protein_id", ""),
            "uniprot": r.get("uniprot", ""),
            "pdb": r.get("pdb", ""),
        }
        h = build(row, meta)
        # protein_id único para el selector de la GUI
        h["protein_id"] = r.get("uniprot") or r.get("protein_id", "")
        h["source_dataset"] = "UNCHAR_AFDB_PRESENT"
        h["review_decision"] = "pending"
        out.append(h)
    return out
