#!/usr/bin/env python3
"""00 — Comprueba que herramientas y datos de referencia estan disponibles.

Escribe `data/taxonomy/tooling.json` con lo que hay y lo que falta. No descarga
nada: solo detecta. La descarga puntual de fastANI la hace 00_setup_fastani.sh y
la de la taxonomia GTDB la hace 03_assign_taxonomy.py.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

TOOLS = {
    "fastANI": ["fastANI", "fastani"],
    "skani": ["skani"],
    "gtdb-tk": ["gtdb-tk", "gtdbtk"],
    "ncbi-datasets": ["datasets"],
    "mash": ["mash"],
    "sourmash": ["sourmash"],
    "blastn": ["blastn"],
    "makeblastdb": ["makeblastdb"],
}


def which_any(names: list[str]) -> str | None:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def version_of(path: str) -> str | None:
    for flag in ("--version", "version"):
        try:
            out = subprocess.run([path, flag], capture_output=True, text=True, timeout=30)
        except Exception:
            continue
        text = (out.stdout + out.stderr).strip().splitlines()
        if out.returncode == 0 and text:
            return text[0][:200]
    return None


def main() -> int:
    L.ensure_dirs()
    found: dict[str, dict] = {}
    for label, names in TOOLS.items():
        path = which_any(names)
        found[label] = {
            "available": bool(path),
            "path": path,
            "version": version_of(path) if path else None,
        }

    # fastANI descargado a mano en el cache del modulo
    local_fastani = L.WORK / "bin" / "fastANI"
    if local_fastani.exists():
        found["fastANI_cached"] = {
            "available": True,
            "path": str(local_fastani),
            "version": version_of(str(local_fastani)),
        }

    data = {
        "assembly_summary": {
            "path": str(L.ASSEMBLY_SUMMARY),
            "exists": L.ASSEMBLY_SUMMARY.exists(),
            "bytes": L.ASSEMBLY_SUMMARY.stat().st_size if L.ASSEMBLY_SUMMARY.exists() else 0,
        },
        "taxdump": {
            "path": str(L.TAXDUMP),
            "exists": (L.TAXDUMP / "nodes.dmp").exists() and (L.TAXDUMP / "names.dmp").exists(),
        },
        "pilot_tsv": {"path": str(L.PILOT_TSV), "exists": L.PILOT_TSV.exists()},
        "local_fasta_dirs": {
            str(d): len(list(d.glob("*.fna.gz"))) if d.exists() else 0 for d in L.RAW_DIRS
        },
        "gtdb_taxonomy_cached": {
            "path": str(L.GTDB_DIR / "gtdb_taxonomy.tsv"),
            "exists": (L.GTDB_DIR / "gtdb_taxonomy.tsv").exists(),
        },
        "work_dir": str(L.WORK),
    }

    payload = {
        "tools": found,
        "data": data,
        "missing": [k for k in ("skani", "gtdb-tk", "ncbi-datasets") if not found[k]["available"]],
    }
    L.write_json(L.OUT / "tooling.json", payload)

    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
