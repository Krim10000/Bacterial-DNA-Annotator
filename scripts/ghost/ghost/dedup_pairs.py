#!/usr/bin/env python3
"""Análisis de pares para el experimento de deduplicación.

Usa el all-vs-all de DIAMOND para definir la verdad: dos proteínas "deberían"
agruparse si son el mismo gen (misma clave GN/acceso) y su identidad real es
mayor o igual que el umbral. Con eso calcula precisión y recall reales por
umbral, en lugar del recall engañoso sobre ortólogos lejanos.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from .dedup_eval import collect_annotations, sample_rows


def parse_allvall(path: Path) -> dict[tuple[str, str], float]:
    """Identidad por par no ordenado, con cobertura >=80% en ambas secuencias."""
    pairs: dict[tuple[str, str], float] = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            q, s, pid, aln, ql, sl = p[0], p[1], float(p[2]), int(p[3]), int(p[4]), int(p[5])
            if q == s:
                continue
            if aln / max(1, ql) < 0.8 or aln / max(1, sl) < 0.8:
                continue
            key = (q, s) if q < s else (s, q)
            if pid > pairs.get(key, 0.0):
                pairs[key] = pid
    return pairs


def cluster_of(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out[p[1]] = p[0]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--sample", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    rows = collect_annotations(a.data_root)
    sample = sample_rows(rows, a.sample, a.seed)
    key_of = {r["pid"]: r["key"] for r in sample}

    allvall = parse_allvall(a.out_dir / "allvall.tsv")
    print(f"pares con identidad y cobertura >=80%: {len(allvall):,}")

    results = []
    for tsv in sorted(a.out_dir.glob("clusters_*.tsv"),
                      key=lambda p: float(p.stem.split("_")[1])):
        t = float(tsv.stem.split("_")[1])
        cmap = cluster_of(tsv)
        grouped = set()
        by_cluster: dict[str, list[str]] = defaultdict(list)
        for pid, c in cmap.items():
            by_cluster[c].append(pid)
        for members in by_cluster.values():
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    x, y = members[i], members[j]
                    grouped.add((x, y) if x < y else (y, x))

        grouped_same = sum(1 for (x, y) in grouped if key_of.get(x) == key_of.get(y))
        precision = grouped_same / len(grouped) if grouped else 1.0

        truth = {(x, y) for (x, y), pid in allvall.items()
                 if pid >= t and key_of.get(x) == key_of.get(y)}
        hit = len(grouped & truth)
        recall = hit / len(truth) if truth else 1.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

        results.append({
            "threshold": t,
            "grouped_pairs": len(grouped),
            "same_gene_grouped": grouped_same,
            "precision": round(precision, 4),
            "truth_pairs": len(truth),
            "recall": round(recall, 4),
            "f1": round(f1, 4),
        })
        print(f"  {t:>3.0f}%: agrupados={len(grouped):>6} precisión={precision:.3f} "
              f"recall={recall:.3f} F1={f1:.3f}")

    ok = [r for r in results if r["precision"] >= 0.95]
    chosen = max(ok, key=lambda r: r["recall"]) if ok else max(results, key=lambda r: r["f1"])
    print(f"\nelegido: {chosen['threshold']:.0f}% (precisión>={0.95}, máximo recall)")

    with open(a.out_dir / "dedup_pairs.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    (a.out_dir / "dedup_pairs_choice.json").write_text(
        json.dumps({"chosen": chosen, "criterion": "precisión>=0.95, máximo recall",
                    "results": results}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
