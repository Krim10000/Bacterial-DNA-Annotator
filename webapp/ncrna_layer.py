"""ncRNA con Infernal/Rfam (CPU, opcional).

Usa un subconjunto Rfam local (`data/rfam_models/bac_ncRNA.cm`) con cmsearch.
No es Rfam completo y no sustituye la anotación experimental. Si faltan
cmsearch o el panel, devuelve `not_available`; no inventa ncRNA.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

PANEL = Path("rfam_models") / "bac_ncRNA.cm"


def available(data_root: Path) -> tuple[bool, str]:
    if shutil.which("cmsearch") is None:
        return False, "cmsearch no está en PATH"
    panel = Path(data_root) / PANEL
    if not panel.exists():
        return False, f"falta {panel}"
    return True, ""


def _parse_tblout(path: Path, genome: str) -> list[dict[str, Any]]:
    rows = []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 17:
                continue
            try:
                start, end = int(f[7]), int(f[8])
                score, evalue = float(f[14]), float(f[15])
            except ValueError:
                continue
            rows.append({
                "contig": f[0], "start": min(start, end),
                "end": max(start, end), "strand": f[9],
                "type": "ncRNA", "family": f[2], "accession": f[3],
                "score": score, "evalue": evalue,
                "evidence": "cmsearch_rfam", "genome": genome,
            })
    return rows


def run(sequences: dict[str, str], features: list[dict], outdir: Path,
        genome: str, data_root: Path, cpu: int = 4) -> dict[str, Any]:
    ok, why = available(data_root)
    if not ok:
        return {"ok": False, "reason": why, "rows": [], "n": 0}
    if not sequences:
        return {"ok": True, "reason": "", "rows": [], "n": 0}
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ncrna_") as td:
        tmp = Path(td)
        fa = tmp / "genome.fa"
        with fa.open("w", encoding="utf-8") as fh:
            for name, seq in sequences.items():
                fh.write(f">{name}\n{seq}\n")
        tbl = tmp / "hits.tbl"
        cmd = ["cmsearch", "--rfam", "--cut_ga", "--noali",
               "--cpu", str(max(1, cpu)), "--tblout", str(tbl),
               str(Path(data_root) / PANEL), str(fa)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            return {"ok": False, "reason": (proc.stderr or proc.stdout)[-400:],
                    "rows": [], "n": 0}
        rows = _parse_tblout(tbl, genome)
    # Agrega features ncRNA al resultado y evita IDs colisionando.
    existing = {str(f.get("id", "")) for f in features}
    for i, row in enumerate(rows, 1):
        fid = f"ncrna_{i}"
        while fid in existing:
            i += 1
            fid = f"ncrna_{i}"
        existing.add(fid)
        row["id"] = fid
        row["product"] = row["family"]
        row["source"] = "cmsearch_rfam"
        features.append({**row, "gene_id": fid})
    path = outdir / f"{genome}.ncrna.tsv"
    cols = ["genome", "contig", "start", "end", "strand", "type",
            "family", "accession", "score", "evalue", "evidence"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    return {"ok": True, "reason": "", "rows": rows, "n": len(rows),
            "tsv": str(path), "features": features}
