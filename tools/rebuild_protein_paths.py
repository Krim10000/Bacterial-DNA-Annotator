#!/usr/bin/env python3
"""Reconstruye protein_paths.tsv (A/B/C) desde functional.tsv ya calculados.

No vuelve a correr DIAMOND. Sirve para los 400, donde functional.tsv existe
pero protein_paths.tsv no estaba persistido.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "webapp"))
from functional import _row_for, parse_hit  # noqa: E402


def proteins(path: Path) -> list[dict]:
    out, name, seq = [], None, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                out.append({"id": name, "aa": "".join(seq)})
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name:
        out.append({"id": name, "aa": "".join(seq)})
    return out


def hits(path: Path) -> dict[str, dict]:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        c = line.split("\t")
        if len(c) < 7:
            continue
        try:
            h = parse_hit(c[0], c[1], c[-1] if len(c) >= 13 else "",
                          float(c[2]), int(c[3]),
                          float(c[10]) if len(c) > 10 else 1.0,
                          float(c[11]) if len(c) > 11 else 0.0)
        except ValueError:
            continue
        if c[0] not in out or h.bitscore > out[c[0]].bitscore:
            out[c[0]] = h
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    ann = a.data_root / "annotations_ab_initio"
    dirs = sorted(p for p in ann.iterdir() if p.is_dir())
    if a.limit:
        dirs = dirs[:a.limit]
    total = Counter()
    for i, d in enumerate(dirs, 1):
        out = d / "protein_paths.tsv"
        if out.exists() and not a.force:
            continue
        fa = d / "proteins.faa"
        ft = d / "functional.tsv"
        if not fa.exists() or not ft.exists():
            continue
        prots = proteins(fa)
        by_id = hits(ft)
        rows = [_row_for(p, by_id.get(p["id"])) for p in prots]
        cols = list(rows[0]) if rows else ["protein_id"]
        with out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
            w.writeheader(); w.writerows(rows)
        total.update(r["path"] for r in rows)
        if i % 50 == 0:
            print(f"{i}/{len(dirs)} A={total['A']} B={total['B']} C={total['C']}",
                  flush=True)
    print("total", dict(total))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
