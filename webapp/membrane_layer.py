"""Membrana/transporte para un FASTA nuevo.

Envoltorio CPU sobre `scripts/membrane/03_annotate_membrane.py`: escribe GFF,
proteínas y functional.tsv en un directorio temporal y ejecuta hmmscan/Pfam.
No sustituye TMHMM/SignalP; son heurísticas declaradas.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "membrane" / "03_annotate_membrane.py"


def available(data_root: Path) -> tuple[bool, str]:
    if not SCRIPT.exists():
        return False, f"falta {SCRIPT}"
    if shutil.which("hmmscan") is None:
        return False, "hmmscan no está en PATH"
    panel = Path(data_root) / "membrane" / "pfam_panel"
    if not (panel / "membrane_panel.hmm").exists():
        return False, f"falta {panel}/membrane_panel.hmm"
    return True, ""


def _write_functional(res: dict[str, Any], path: Path) -> None:
    hits = (res.get("functional") or {}).get("hits") or {}
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        for pid, h in hits.items():
            w.writerow([pid, h.get("uniprot", ""), h.get("pident", 0),
                        h.get("aln_len", 0), 0, 0, 0, h.get("aln_len", 0),
                        0, h.get("aln_len", 0), h.get("evalue", 1.0),
                        h.get("bitscore", 0.0), h.get("product", "")])


def run(res: dict[str, Any], outdir: Path, data_root: Path,
        cpu: int = 4) -> dict[str, Any]:
    ok, why = available(data_root)
    if not ok:
        return {"ok": False, "reason": why, "rows": [], "tsv": "", "gff": ""}
    genome = str(res.get("genome", "genome"))
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="membrane_") as td:
        tmp = Path(td)
        annot = tmp / genome
        annot.mkdir(parents=True)
        (annot / "annotation.gff3").write_text(res.get("gff", ""), encoding="utf-8")
        (annot / "proteins.faa").write_text(res.get("proteins_fasta", ""),
                                            encoding="utf-8")
        _write_functional(res, annot / "functional.tsv")
        panel = Path(data_root) / "membrane" / "pfam_panel"
        cmd = [sys.executable, str(SCRIPT), "--genome", genome,
               "--annot-root", str(tmp), "--panel-dir", str(panel),
               "--outdir", str(tmp / "out"), "--cpu", str(max(1, cpu))]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if proc.returncode != 0:
            return {"ok": False, "reason": (proc.stderr or proc.stdout)[-500:],
                    "rows": [], "tsv": "", "gff": ""}
        tsv = tmp / "out" / f"{genome}.tsv"
        gff = tmp / "out" / "gff" / f"{genome}.gff3"
        if not tsv.exists():
            return {"ok": False, "reason": "sin TSV de membrana",
                    "rows": [], "tsv": "", "gff": ""}
        final_tsv = outdir / f"{genome}.membrane.tsv"
        final_gff = outdir / f"{genome}.membrane.gff3"
        shutil.copyfile(tsv, final_tsv)
        if gff.exists():
            shutil.copyfile(gff, final_gff)
        rows = list(csv.DictReader(tsv.open(encoding="utf-8"), delimiter="\t"))
    return {"ok": True, "reason": "", "rows": rows, "n_rows": len(rows),
            "tsv": str(final_tsv), "gff": str(final_gff)}
