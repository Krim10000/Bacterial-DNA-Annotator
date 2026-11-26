#!/usr/bin/env python3
"""Inventario de sigma por especie; no impone seis factores globales."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from webapp.sigma_inventory import load_inventory, summary, to_tsv

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", type=Path, default=Path("data"))
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    factors = load_inventory(a.data_root, a.genome)
    if not factors:
        print(f"{a.genome}: no hay genes.tsv o no se detectaron sigma", file=sys.stderr)
        return 1
    s = summary(factors)
    print(f"{a.genome}: {s['n_factors']} sigma detectados {s['families']}")
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(to_tsv(factors), encoding="utf-8")
        a.out.with_suffix(".json").write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"-> {a.out}")
    else:
        print(to_tsv(factors))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
