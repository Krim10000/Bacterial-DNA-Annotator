"""Anotación estructural automática con EnzyMM/PyJess (plantillas M-CSA).

EnzyMM busca arreglos geométricos de residuos catalíticos en un PDB.
Un hit es una hipótesis de plantilla, no una función confirmada.
"""
from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "metabolism" / "targets"
CACHE = DATA / "structures" / "annotation_cache"


def _exe() -> str:
    p = Path(sys.executable).parent / "enzymm"
    if p.exists():
        return str(p)
    return shutil.which("enzymm") or "enzymm"


def _cache_key(pdb: Path) -> str:
    return "enzymm_" + hashlib.sha256(pdb.read_bytes()).hexdigest()[:20]


def scan_pdbs(pdbs: list[Path | str], outdir: Path | str, jobs: int = 4) -> dict:
    """Escanea uno o varios PDB con EnzyMM. Devuelve {hits, no_template, error, tsv, residues_tsv}."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = [Path(p) for p in pdbs if Path(p).exists()]
    if not paths:
        return {"hits": [], "no_template": [], "error": "sin PDB", "tsv": "", "residues_tsv": ""}
    # EnzyMM list mode loads templates once
    list_file = out / "query_pdbs.list"
    list_file.write_text("\n".join(str(p) for p in paths) + "\n")
    tsv = out / "enzymm.tsv"
    cmd = [_exe(), "-l", str(list_file), "-o", str(tsv), "-j", str(max(1, jobs)),
           "--per-residue-results", "-v"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if not tsv.exists():
        return {"hits": [], "no_template": [str(p) for p in paths],
                "error": f"EnzyMM {r.returncode}: {r.stderr[-300:]}", "tsv": "", "residues_tsv": ""}
    hits = []
    with open(tsv) as fh:
        rows = [line for line in fh if not line.startswith("#")]
    if rows:
        for row in csv.DictReader(rows, delimiter="\t"):
            hits.append(row)
    queried = {h.get("query_id", "") for h in hits}
    no_template = [str(p) for p in paths if str(p) not in queried and p.name not in queried and p.stem not in queried]
    return {"hits": hits, "no_template": no_template, "error": "",
            "tsv": str(tsv), "residues_tsv": str(tsv).replace(".tsv", ".residues.tsv")}


def summarize(hits: list[dict]) -> list[dict]:
    """Resumen legible por candidato/plantilla."""
    out = []
    for h in hits:
        out.append({
            "query_id": h.get("query_id", ""),
            "template_pdb_id": h.get("template_pdb_id", ""),
            "template_mcsa_id": h.get("template_mcsa_id", ""),
            "template_uniprot_id": h.get("template_uniprot_id", ""),
            "template_ec": h.get("template_ec", ""),
            "rmsd": h.get("rmsd", ""),
            "log_evalue": h.get("log_evalue", ""),
            "matched_residues": h.get("matched_residues", ""),
            "template_effective_size": h.get("template_effective_size", ""),
            "predicted_correct": h.get("predicted_correct", ""),
            "source": "EnzyMM/PyJess · M-CSA",
            "limitation": "Match de plantilla geométrica; no confirma reacción, sustrato ni cierre.",
        })
    return out
