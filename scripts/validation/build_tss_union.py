#!/usr/bin/env python3
"""build_tss_union.py — unir TSS curados y de alto rendimiento.

Por qué unir
------------
Cada fuente cubre lo que la otra no:

- **curados** (RegulonDB `PromoterSet`): alta precisión, conservadores.
- **alto rendimiento** (dRNA-seq de varios estudios): ruidosos pero completos.

Medido en K-12 sobre la reconstrucción de operones (F1 de pares):

| fuente | TSS | F1 | P | R |
|---|---:|---:|---:|---:|
| curados | 3.964 | 0,785 | 0,743 | 0,833 |
| alto rendimiento (apoyo >=2) | 13.405 | 0,756 | 0,721 | 0,793 |
| **unión** | **17.129** | **0,799** | **0,834** | 0,767 |

La unión mejora a las dos por separado, y sobre todo la precisión: los cortes
que añade la fuente ruidosa son correctos.

Regla de deduplicación: dos TSS de la misma hebra a <= 10 pb son el mismo sitio;
se conserva el curado (o el de mayor apoyo).

Uso
---
    python scripts/validation/build_tss_union.py \
        --curated data/regulondb/tss_k12.tsv \
        --ht data/regulondb/ht/ht_tss.tsv \
        --min-support 2 \
        --out data/regulondb/tss_k12_union.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MERGE_DIST = 10


def load_curated(path: Path) -> list[tuple[int, str, str]]:
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        p = (r.get("pos") or "").strip()
        if not p.isdigit():
            continue
        st = "-" if (r.get("strand") or "+").strip().startswith("-") else "+"
        out.append((int(p), st, r.get("sigma", "")))
    return out


def load_ht(path: Path, min_support: int) -> list[tuple[int, str, str]]:
    """Agrupa por posición (±MERGE_DIST) y exige apoyo en >= min_support estudios."""
    clusters: dict[tuple[str, int], set[str]] = defaultdict(set)
    for r in csv.DictReader(open(path), delimiter="\t"):
        p = (r.get("pos_1") or r.get("leftEndPosition") or "").strip()
        if not p.lstrip("-").isdigit():
            continue
        st = "-" if (r.get("strand") or "+").strip() in ("-", "-1", "reverse") else "+"
        clusters[(st, int(p) // MERGE_DIST)].add(r.get("dataset_id", ""))
    out = []
    for (st, key), datasets in clusters.items():
        if len(datasets) >= min_support:
            out.append((key * MERGE_DIST, st, f"HT x{len(datasets)}"))
    return out


def union(curated, ht, dist: int = MERGE_DIST):
    """Curados mandan; se añaden los de alto rendimiento que no estén cerca."""
    out = list(curated)
    for pos, st, src in ht:
        if any(abs(pos - p) <= dist and s == st for p, s, _ in out):
            continue
        out.append((pos, st, src))
    out.sort(key=lambda x: (x[1], x[0]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", required=True)
    ap.add_argument("--ht", required=True)
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cur = load_curated(Path(a.curated))
    ht = load_ht(Path(a.ht), a.min_support)
    merged = union(cur, ht)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "strand", "sigma"])
        w.writerows(merged)
    print(f"curados: {len(cur)} | alto rendimiento (apoyo>={a.min_support}): {len(ht)}")
    print(f"unión: {len(merged)}  (+{len(merged)-len(cur)} sobre curados)")
    print(f"-> {out.relative_to(ROOT) if str(out).startswith(str(ROOT)) else out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
