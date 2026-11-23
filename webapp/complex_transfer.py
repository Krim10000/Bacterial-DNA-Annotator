"""Carga las transferencias K-12 → genoma como sospecha fuerte, no cierre.

`transfer_<genoma>.tsv` viene de `scripts/complexes/transfer_orthology.py`.
Se usa cuando el genoma ya está en los 400. Un genoma nuevo no tiene este
archivo: se anota con el modelo K-12 (`complex_model.py`) y las firmas.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def load(data_root: Path, genome: str, features: list[dict] | None = None) -> dict[str, Any]:
    base = Path(genome).stem
    candidates = [genome, base, base.replace("_genomic", "")]
    path = next((Path(data_root) / "complexes" / f"transfer_{g}.tsv"
                 for g in candidates
                 if (Path(data_root) / "complexes" / f"transfer_{g}.tsv").exists()),
                Path(data_root) / "complexes" / f"transfer_{genome}.tsv")
    if not path.exists():
        return {"ok": False, "reason": f"falta {path}", "complexes": [],
                "by_protein": {}}
    by_id = {str(f.get("id", "")): f for f in (features or [])
             if f.get("type") == "CDS"}
    # Para un genoma de los 400, el ID de transferencia es el de genes.tsv,
    # no el cds_* del pipeline. Se agrega ese mapa de coordenadas.
    for cand in (base, base.replace("_genomic", "")):
        gp = Path(data_root) / "regulation" / cand / "genes.tsv"
        if not gp.exists():
            continue
        with gp.open(encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gid = r.get("gene_id") or r.get("locus_tag")
                if gid:
                    by_id[str(gid)] = {
                        "id": str(gid), "type": "CDS",
                        "contig": r.get("contig", ""),
                        "start": int(float(r.get("start") or 0)),
                        "end": int(float(r.get("end") or 0)),
                        "strand": r.get("strand", "+"),
                    }
        break
    feature_by_coord = {(str(f.get("contig","")), int(f.get("start",0)),
                         int(f.get("end",0))): f
                        for f in (features or []) if f.get("type") == "CDS"}
    complexes, by_protein = [], {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ids = [x for x in (r.get("genes_genoma") or "").split(",")
                   if x in by_id]
            if len(ids) < 2:
                continue
            coords = [(int(by_id[x]["start"]), int(by_id[x]["end"])) for x in ids]
            frac = float(r.get("frac_mapped") or 0)
            expected = int(r.get("expected_n") or 0)
            observed = len(ids)
            status = ("complex_suspected" if frac >= 0.8
                      else "complex_candidate")
            cid = "K12:" + str(r.get("complex_id", ""))
            c = {
                "complex_id": cid,
                "name": "K-12 transferido: " + str(r.get("complex_name", "")),
                "contig": by_id[ids[0]].get("contig", ""),
                "start": min(a for a, _ in coords),
                "end": max(b for _, b in coords),
                "subunits": ";".join(ids), "expected": expected,
                "missing": (f"{max(0, expected - observed)} subunidades"
                            if expected else ""),
                "status": status, "flag": "yellow",
                "confidence": "high" if frac >= 0.8 else "medium",
                "evidence": (f"RBH K-12; frac={frac:.2f}; "
                             f"identidad={r.get('identity_mean','')}%; "
                             f"phylo={r.get('phylo_jacc','')}; "
                             f"masa={r.get('mass_kda','')} kDa"),
                "source": "transferencia por ortología K-12",
                "ecs": r.get("ecs", ""), "score": r.get("score", ""),
                "limitation": ("Sospecha fuerte por ortología, no interacción "
                               "demostrada. No cierra brechas."),
            }
            complexes.append(c)
            for pid in ids:
                payload = {
                    "complex_id": cid, "complex_name": c["name"],
                    "complex_status": status, "complex_role": "candidate",
                    "complex_confidence": c["confidence"],
                }
                by_protein[pid] = payload
                g = by_id.get(pid, {})
                f = feature_by_coord.get((str(g.get("contig","")),
                                          int(g.get("start",0)),
                                          int(g.get("end",0))))
                if f:
                    by_protein[str(f.get("id", ""))] = payload
    return {"ok": bool(complexes), "reason": "", "complexes": complexes,
            "by_protein": by_protein, "n_candidate": len(complexes),
            "n_complete": 0, "n_incomplete": 0, "n_ambiguous": 0}
