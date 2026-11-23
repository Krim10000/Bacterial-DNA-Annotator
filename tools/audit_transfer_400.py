#!/usr/bin/env python3
"""Audita transferencias K-12 contra las brechas de los 400.

No cierra brechas. Solo mide si una EC del complejo transferido produce el
metabolito de la brecha según Rhea.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_ec_products(data_root: Path) -> dict[str, set[str]]:
    with open(data_root / "metabolism" / "rhea_reactions.pkl", "rb") as fh:
        reactions = pickle.load(fh)
    out: dict[str, set[str]] = defaultdict(set)
    with open(data_root / "reactions" / "tsv" / "rhea2ec.tsv",
              encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or not c[3] or not c[0].isdigit():
                continue
            out[c[3]] |= set(reactions.get(int(c[0]), {}).get("products", set()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "reports" / "transfer_gap_audit.tsv")
    a = ap.parse_args()
    ec_products = load_ec_products(a.data_root)
    gaps_by_genome: dict[str, set[str]] = defaultdict(set)
    with open(a.data_root / "metabolism" / "gap_interpretation.tsv",
              encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("chebi"):
                gaps_by_genome[r.get("genome", "")].add(r["chebi"])
    rows = []
    details = []
    for path in sorted((a.data_root / "complexes").glob("transfer_*.tsv")):
        genome = path.stem.removeprefix("transfer_")
        gaps = gaps_by_genome.get(genome, set())
        with path.open(encoding="utf-8", errors="replace") as fh:
            tr = list(csv.DictReader(fh, delimiter="\t"))
        matched: set[str] = set()
        strong: set[str] = set()
        for r in tr:
            products = set()
            for ec in str(r.get("ecs", "")).split(";"):
                if ec:
                    products |= ec_products.get(ec.strip(), set())
            hit = gaps & products
            if hit:
                matched |= hit
                if float(r.get("frac_mapped") or 0) >= 0.8:
                    strong |= hit
                for chebi in hit:
                    details.append({"genome": genome, "chebi": chebi,
                                    "complex_id": r.get("complex_id", ""),
                                    "complex_name": r.get("complex_name", ""),
                                    "frac_mapped": r.get("frac_mapped", ""),
                                    "status": ("complex_suspected"
                                               if float(r.get("frac_mapped") or 0) >= 0.8
                                               else "complex_candidate")})
        rows.append({"genome": genome, "gaps": len(gaps),
                     "transferred": len(tr),
                     "strong": sum(1 for r in tr
                                   if float(r.get("frac_mapped") or 0) >= 0.8),
                     "gaps_with_transfer": len(matched),
                     "gaps_with_strong_transfer": len(strong),
                     "gaps_unassociated": len(gaps - matched)})
    a.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0]) if rows else ["genome"]
    with a.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    detail_path = a.out.with_name("transfer_gap_details.tsv")
    with detail_path.open("w", newline="", encoding="utf-8") as fh:
        cols = list(details[0]) if details else ["genome"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(details)
    summary = {k: sum(r[k] for r in rows) for k in
               ("gaps", "transferred", "strong", "gaps_with_transfer",
                "gaps_with_strong_transfer", "gaps_unassociated")}
    (a.out.with_suffix(".json")).write_text(json.dumps({"summary": summary,
                                                        "genomes": len(rows)},
                                                       indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
