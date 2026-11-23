#!/usr/bin/env python3
"""Cobertura de los nombres 'subunit' que quedan fuera del catálogo K-12.

No decide si son subunidades reales. Mide si el modelo de complejos las cubre
como candidatas y cuáles no tienen ninguna cobertura.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAT = re.compile(r"\bsubunit\b|\bcomponent\b", re.I)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", type=Path,
                    default=ROOT / "data" / "regulation" / "ecoli_k12" / "genes.tsv")
    ap.add_argument("--gold", type=Path,
                    default=ROOT / "data" / "complexes" / "gold_k12.tsv")
    ap.add_argument("--clusters", type=Path,
                    default=ROOT / "data" / "complexes" / "clusters_ecoli_k12_evidencia.tsv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "reports" / "subunit_coverage.json")
    a = ap.parse_args()
    genes = {}
    with a.genes.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("locus_tag"):
                genes[r["locus_tag"]] = r
    in_gold = set()
    with a.gold.open(encoding="utf-8", errors="replace") as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            for b in (c.get("locus_tags") or "").split(","):
                if b:
                    in_gold.add(b)
    outside = [b for b, r in genes.items()
               if PAT.search(r.get("product", "") or "") and b not in in_gold]
    cluster_of = {}
    if a.clusters.exists():
        with a.clusters.open(encoding="utf-8", errors="replace") as fh:
            for c in csv.DictReader(fh, delimiter="\t"):
                for b in (c.get("genes") or "").split(","):
                    if b:
                        cluster_of[b] = c
    covered = [b for b in outside if b in cluster_of]
    uncovered = [b for b in outside if b not in cluster_of]
    report = {
        "outside_catalog": len(outside),
        "covered_by_model": len(covered),
        "covered_pct": round(100 * len(covered) / max(1, len(outside)), 1),
        "uncovered": len(uncovered),
        "uncovered_pct": round(100 * len(uncovered) / max(1, len(outside)), 1),
        "model_evidence": dict(Counter(cluster_of[b].get("evidencia", "")
                                       for b in covered)),
        "covered_examples": [
            {"locus": b, "gene": genes[b].get("gene_name", ""),
             "product": genes[b].get("product", ""),
             "evidencia": cluster_of[b].get("evidencia", "")}
            for b in covered[:20]],
        "uncovered_examples": [
            {"locus": b, "gene": genes[b].get("gene_name", ""),
             "product": genes[b].get("product", "")}
            for b in uncovered],
        "caveat": ("Cubierto por el modelo = candidato, no miembro confirmado. "
                   "Fuera del catálogo no significa falso."),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report, ensure_ascii=False, indent=2),
                     encoding="utf-8")
    print(json.dumps({k: report[k] for k in
                      ("outside_catalog", "covered_by_model", "covered_pct",
                       "uncovered", "uncovered_pct", "model_evidence")},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
