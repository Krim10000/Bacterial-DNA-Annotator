"""Estado estructural del camino B.

No pliega ni filtra: busca modelos AFDB locales y deja un manifiesto. Un modelo
disponible no es evidencia de función ni de complejo; solo habilita el siguiente
paso (filtros/estructura).
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def _model_dir(data_root: Path) -> Path:
    return Path(data_root) / "metabolism" / "targets" / "afdb25k"


def available(data_root: Path) -> tuple[bool, str]:
    d = _model_dir(data_root)
    if d.exists() and any(d.glob("*.pdb")):
        return True, ""
    return False, f"sin modelos AFDB en {d}"


def run(rows: list[dict], outdir: Path, genome: str,
        data_root: Path) -> dict[str, Any]:
    ok, why = available(data_root)
    models = _model_dir(data_root)
    manifest = []
    for r in rows:
        if r.get("path") != "B":
            continue
        acc = str(r.get("uniprot") or "").split("|")[-1]
        pdb = models / f"{acc}.pdb"
        if acc and pdb.exists():
            r["structural_status"] = "afdb_model_available"
            r["structural_evidence"] = "local_afdb_model"
            r["structural_note"] = ("Modelo AFDB disponible; no evaluado con "
                                    "filtros, no confirma función ni complejo.")
            manifest.append({"genome": genome, "protein_id": r.get("protein_id", ""),
                             "uniprot": acc, "pdb": str(pdb)})
        else:
            r["structural_status"] = "no_afdb_model"
            r["structural_evidence"] = "none"
            r["structural_note"] = ("Sin modelo AFDB local; se necesita descarga "
                                    "o plegado para evaluar estructura.")
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{genome}.structure_b_manifest.tsv"
    cols = ["genome", "protein_id", "uniprot", "pdb"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(manifest)
    return {"ok": ok, "reason": why, "n_available": len(manifest),
            "manifest": str(path)}
