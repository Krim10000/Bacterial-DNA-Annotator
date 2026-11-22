"""Hipótesis estructurales revisables a partir de los 5 filtros.

No sustituyen la anotación base ni cierran brechas.
"""
from __future__ import annotations

import json
from typing import Any

SCHEMA = [
    "genome", "protein_id", "uniprot", "pdb", "hypothesis", "ec", "cofactor",
    "pocket", "template_mcsa", "foldseek_top", "f1_foldseek", "f2_enzymm",
    "f3_pocket", "f4_alphafill", "f5_active_site", "n_sources", "confidence",
    "review_decision", "limitation", "evidence_json",
]


def _uniq(items) -> list[str]:
    out = []
    for x in items:
        x = str(x).strip()
        if x and x not in out:
            out.append(x)
    return out


def build(row: dict[str, Any], meta: dict[str, Any] | None = None) -> dict[str, Any]:
    meta = meta or {}
    fs = row.get("f2_enzymm")
    f3 = row.get("f3_pocket")
    f4 = row.get("f4_alphafill")
    f5 = row.get("f5_active_site")
    ecs = _uniq((row.get("enzymm_ecs") or "").split(";"))
    templates = _uniq(str(row.get("enzymm_templates") or "").split(";"))
    ligands = _uniq(str(row.get("alphafill_ligands") or "").split(";"))
    pockets = row.get("p2rank_pockets", "")
    if ecs:
        hypothesis = "posible EC " + ", ".join(ecs[:3])
    elif fs:
        hypothesis = "posible enzima (plantilla M-CSA " + ", ".join(templates[:3]) + ")"
    elif f4 and ligands:
        hypothesis = "posible unión a " + ", ".join(x.split("(")[0] for x in ligands[:3])
    else:
        hypothesis = "sin hipótesis estructural"
    if fs and f3 and f4:
        confidence = "supported_hypothesis"
    elif fs and f5:
        confidence = "high_priority"
    elif fs:
        confidence = "orientative"
    elif row.get("f1_foldseek") and (f3 or f4):
        confidence = "structural_only"
    else:
        confidence = "uncharacterized"
    sources = []
    if row.get("f1_foldseek"): sources.append("Foldseek")
    if fs: sources.append("EnzyMM/M-CSA")
    if f3: sources.append("P2Rank")
    if f4: sources.append("AlphaFill")
    if f5: sources.append("UniProt/estructura")
    evidence = {k: row.get(k, "") for k in (
        "foldseek_top", "foldseek_bits", "enzymm_templates", "enzymm_ecs",
        "p2rank_pockets", "pocket_active_hits", "alphafill_ligands",
        "active_positions", "error", "source")}
    return {
        "genome": meta.get("genome", row.get("genome", "")),
        "protein_id": meta.get("protein_id", ""),
        "uniprot": meta.get("uniprot", row.get("uniprot", "")),
        "pdb": meta.get("pdb", row.get("pdb", "")),
        "hypothesis": hypothesis,
        "ec": ";".join(ecs),
        "cofactor": ";".join(ligands[:5]),
        "pocket": str(pockets),
        "template_mcsa": ";".join(templates[:5]),
        "foldseek_top": row.get("foldseek_top", ""),
        "f1_foldseek": row.get("f1_foldseek", ""),
        "f2_enzymm": row.get("f2_enzymm", ""),
        "f3_pocket": row.get("f3_pocket", ""),
        "f4_alphafill": row.get("f4_alphafill", ""),
        "f5_active_site": row.get("f5_active_site", ""),
        "n_sources": row.get("n_sources", len(sources)),
        "confidence": confidence,
        "review_decision": "pending",
        "limitation": "Hipótesis computacional; no confirma función ni cierra brechas.",
        "evidence_json": json.dumps({"row": row, "sources": sources, "meta": meta},
                                    ensure_ascii=False)[:200000],
    }


def from_manifest(rows: list[dict], filters: list[dict]) -> list[dict]:
    """Une manifiesto y resultados de filtros por pdb."""
    by_pdb = {str(x.get("pdb", "")): x for x in filters}
    out = []
    for m in rows:
        f = by_pdb.get(str(m.get("pdb", "")))
        if not f:
            continue
        out.append(build(f, m))
    return out
