#!/usr/bin/env python3
"""00_make_exclusions.py — claves a excluir de un lote nuevo.

Escribe un TSV `genome, gap` con las claves de las tablas indicadas (p. ej. el
`gap_evidence.tsv` del lote anterior). `01_select_worklist.py` lo consume vía
`EVIDENCE_EXCLUDE` para no repetir brechas ya cubiertas.

Uso:
    python 00_make_exclusions.py --from /ruta/gap_evidence.tsv \
        --out /ruta/work_batch_02/exclude_keys.tsv [--from otra.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="sources", action="append", required=True,
                    help="TSV con columnas genome y gap (repetible)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    keys: set[tuple[str, str]] = set()
    for src in a.sources:
        path = Path(src)
        if not path.exists():
            print(f"aviso: no existe {path}", file=sys.stderr)
            continue
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("genome") and r.get("gap"):
                    keys.add((r["genome"], r["gap"]))
    rows = [{"genome": g, "gap": p} for g, p in sorted(keys)]
    write_tsv(Path(a.out), rows, ["genome", "gap"])
    print(f"exclusiones: {len(rows)} claves -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
