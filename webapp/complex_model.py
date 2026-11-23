"""Puente CPU al buscador de complejos ya existente en `scripts/complexes`.

No reimplementa el modelo: escribe `genes.tsv`/`operons.tsv` desde la
anotación y llama a `build_candidates.py --apply-only` con el modelo K-12.
El resultado es una **sospecha** (`complex_candidate`), no una confirmación.
"""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "complexes" / "build_candidates.py"
sys = __import__("sys")


def available(data_root: Path) -> tuple[bool, str]:
    model = Path(data_root) / "complexes" / "model_k12.json"
    if not SCRIPT.exists():
        return False, f"falta {SCRIPT}"
    if not model.exists():
        return False, f"falta {model} (regenerable con fetch_gold_k12.py --fit)"
    try:
        import sklearn  # noqa: F401
    except Exception:
        return False, "falta scikit-learn (extra [full])"
    return True, ""


def _write_genes(features: list[dict], hits: dict[str, dict], path: Path) -> None:
    cols = ["locus_tag", "gene_name", "ec", "contig", "start", "end",
            "strand", "uniprot"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for f in features:
            if f.get("type") != "CDS":
                continue
            h = hits.get(f.get("id", ""), {})
            w.writerow({
                "locus_tag": f.get("id", ""),
                "gene_name": h.get("gene", ""),
                "ec": h.get("ec", ""),
                "contig": f.get("contig", ""),
                "start": f.get("start", 0),
                "end": f.get("end", 0),
                "strand": f.get("strand", "+"),
                "uniprot": h.get("uniprot", ""),
            })


def _write_operons(features: list[dict], path: Path, max_gap: int = 300) -> None:
    try:  # paquete webapp
        from .complexes import neighborhoods
    except ImportError:
        from complexes import neighborhoods
    cols = ["operon_id", "genes"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for i, group in enumerate(neighborhoods(features, max_gap=max_gap), 1):
            ids = [str(f.get("id", "")) for f in group if f.get("type") == "CDS"]
            if len(ids) >= 2:
                w.writerow({"operon_id": f"cx_operon_{i}", "genes": ",".join(ids)})


def run(features: list[dict], hits: dict[str, dict], outdir: Path, genome: str,
        data_root: Path) -> dict[str, Any]:
    """Ejecuta el modelo K-12 sobre un genoma nuevo y devuelve sospechas."""
    ok, why = available(data_root)
    if not ok:
        return {"ok": False, "reason": why, "complexes": [], "by_protein": {}}
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cx_model_") as td:
        tmp = Path(td)
        genes = tmp / "genes.tsv"
        operons = tmp / "operons.tsv"
        _write_genes(features, hits, genes)
        _write_operons(features, operons)
        model = Path(data_root) / "complexes" / "model_k12.json"
        cmd = [
            shutil.which("python") or "python", str(SCRIPT),
            "--apply-only", "--genome", genome,
            "--genes", str(genes), "--operons", str(operons),
            "--model", str(model), "--out-dir", str(outdir),
            "--synteny", "/dev/null", "--string", "/dev/null",
            "--string-species", "0",
        ]
        import sys
        cmd[0] = sys.executable
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if res.returncode != 0:
            return {"ok": False, "reason": (res.stderr or res.stdout)[-500:],
                    "complexes": [], "by_protein": {}}
        cluster_path = outdir / f"clusters_{genome}.tsv"
        rows = list(csv.DictReader(cluster_path.open(encoding="utf-8"),
                                   delimiter="\t")) if cluster_path.exists() else []
    genes_by_id = {f.get("id"): f for f in features if f.get("type") == "CDS"}
    complexes, by_protein = [], {}
    for r in rows:
        ids = [x for x in (r.get("genes") or "").split(",") if x in genes_by_id]
        if len(ids) < 2:
            continue
        coords = [(int(genes_by_id[x]["start"]), int(genes_by_id[x]["end"]))
                  for x in ids]
        cid = r.get("cluster_id", "")
        name = "complejo candidato: " + (r.get("gene_names") or "genes")
        c = {
            "complex_id": cid, "name": name,
            "contig": genes_by_id[ids[0]].get("contig", ""),
            "start": min(a for a, _ in coords), "end": max(b for _, b in coords),
            "subunits": ";".join(ids), "expected": "", "missing": "",
            "status": "complex_candidate",
            "evidence": f"modelo K-12; {r.get('evidence','')}",
            "confidence": r.get("confidence", ""),
            "source": "modelo K-12 (build_candidates.py)",
            "ecs": "", "flag": "yellow",
            "limitation": ("Sospecha por modelo genómico K-12; no demuestra "
                           "interacción ni estequiometría. No cierra brechas."),
            "score": r.get("mean_score", ""),
        }
        complexes.append(c)
        for pid in ids:
            by_protein[pid] = {"complex_id": cid, "complex_name": name,
                               "complex_status": "complex_candidate",
                               "complex_role": "candidate",
                               "complex_confidence": r.get("confidence", "")}
    return {"ok": True, "reason": "", "complexes": complexes,
            "by_protein": by_protein, "n_candidate": len(complexes),
            "raw_clusters": rows}
