#!/usr/bin/env python3
"""Foldseek-Multimer: busca un modelo de complejo en ensamblajes PDB.

No crea la base de ensamblajes ni inventa templates. Si faltan el binario o la
DB, devuelve `not_configured`. Uso:

    python scripts/complexes/foldseek_multimer.py \
        --model complejo.pdb --db data/db/pdb_assemblies --out hits.tsv
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def available(db: Path) -> tuple[bool, str]:
    exe = shutil.which("foldseek")
    if exe is None:
        return False, "foldseek no está en PATH"
    if not db.exists():
        return False, f"falta la base de ensamblajes: {db}"
    return True, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--db", type=Path, default=Path("data/db/pdb_assemblies"))
    ap.add_argument("--out", type=Path, default=Path("foldseek_complex.tsv"))
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    ok, why = available(a.db)
    if not ok:
        print(f"not_configured: {why}")
        return 2
    if not a.model.exists():
        print(f"not_configured: falta el modelo {a.model}")
        return 2
    a.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [shutil.which("foldseek") or "foldseek", "easy-complexsearch",
           str(a.model), str(a.db), str(a.out), "--threads", str(a.threads)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("foldseek falló:", (res.stderr or res.stdout)[-400:])
        return 1
    print("escrito:", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
