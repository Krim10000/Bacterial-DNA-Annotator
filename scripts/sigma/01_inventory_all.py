#!/usr/bin/env python3
"""Inventario de sigma para todas las especies del catálogo.

No asigna sigma ni usa gold: solo cuenta candidatos por especie a partir de
`data/regulation/<genoma>/genes.tsv`. Es la primera capa de la extensión a 400.
"""
from __future__ import annotations
import argparse, csv, json, sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from webapp.sigma_inventory import load_inventory, summary  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "sigma" / "inventory")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    reg = a.data_root / "regulation"
    genomes = sorted(d.name for d in reg.iterdir()
                     if d.is_dir() and (d / "genes.tsv").exists())
    if a.limit:
        genomes = genomes[:a.limit]
    a.out.mkdir(parents=True, exist_ok=True)
    rows, summaries = [], {}
    for i, genome in enumerate(genomes, 1):
        factors = load_inventory(a.data_root, genome)
        s = summary(factors)
        summaries[genome] = s
        for f in factors:
            rows.append({"genome": genome, **f.as_dict()})
        if i % 25 == 0:
            print(f"{i}/{len(genomes)}", flush=True)
    cols = ["genome", "species", "native_name", "protein_id", "locus_tag", "product",
            "family", "evidence", "confidence"]
    with (a.out / "all_sigma_inventory.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    (a.out / "all_sigma_inventory.json").write_text(
        json.dumps({"n_genomes": len(genomes), "summaries": summaries},
                   indent=1, ensure_ascii=False), encoding="utf-8")
    counts = Counter(s["n_factors"] for s in summaries.values())
    fam = Counter()
    for s in summaries.values():
        fam.update(s["families"])
    print(json.dumps({"genomes": len(genomes), "rows": len(rows),
                      "n_factors_hist": dict(sorted(counts.items())),
                      "families": dict(fam)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
