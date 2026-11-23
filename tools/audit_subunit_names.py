#!/usr/bin/env python3
"""Mide si "subunit/component" en el nombre predice pertenencia a complejo.

Usa el gold curado de K-12. No es verdad experimental: "fuera del catálogo"
no significa "no es subunidad"; puede ser un complejo no curado.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAT = re.compile(r"\bsubunit\b|\bcomponent\b", re.I)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", type=Path,
                    default=ROOT / "data" / "regulation" / "ecoli_k12" / "genes.tsv")
    ap.add_argument("--gold", type=Path,
                    default=ROOT / "data" / "complexes" / "gold_k12.tsv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "reports" / "subunit_name_audit.json")
    a = ap.parse_args()
    genes = {}
    with a.genes.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("locus_tag"):
                genes[r["locus_tag"]] = r
    in_complex = set()
    with a.gold.open(encoding="utf-8", errors="replace") as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            for b in (c.get("locus_tags") or "").split(","):
                if b:
                    in_complex.add(b)
    named = [b for b, r in genes.items()
             if PAT.search(r.get("product", "") or "")]
    inside = [b for b in named if b in in_complex]
    outside = [b for b in named if b not in in_complex]
    complex_not_named = [b for b in in_complex
                         if b in genes and not PAT.search(genes[b].get("product", "") or "")]
    report = {
        "genes": len(genes), "complex_subunits": len(in_complex),
        "named": len(named), "named_in_complex": len(inside),
        "named_outside_catalog": len(outside),
        "named_outside_pct": round(100 * len(outside) / max(1, len(named)), 1),
        "name_precision_vs_catalog": round(100 * len(inside) / max(1, len(named)), 1),
        "name_recall_vs_catalog": round(100 * len(inside) / max(1, len(in_complex)), 1),
        "complex_subunits_not_named_pct": round(100 * len(complex_not_named)
                                               / max(1, len(in_complex)), 1),
        "examples_outside": [
            {"locus": b, "product": genes[b].get("product", "")}
            for b in outside[:15]],
        "caveat": ("'Fuera del catálogo' no significa 'no es subunidad': "
                   "puede ser un complejo no curado o no mapeado."),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report, ensure_ascii=False, indent=2),
                     encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
