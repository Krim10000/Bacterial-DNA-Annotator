#!/usr/bin/env python3
"""Ejecuta los 5 filtros estructurales sobre uno o varios PDB."""
import argparse, csv, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
from structural_filters import run_filters  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", nargs="+", type=Path, required=True)
    ap.add_argument("--genome", default="")
    ap.add_argument("--uniprot", default="")
    ap.add_argument("--out", type=Path, default=ROOT / "data/metabolism/targets/STRUCTURAL_FILTERS.tsv")
    a = ap.parse_args()
    rows = [run_filters(p, a.genome, a.uniprot) for p in a.pdb]
    cols = ["pdb", "genome", "uniprot", "f1_foldseek", "foldseek_top", "foldseek_bits",
            "f2_enzymm", "enzymm_hits", "enzymm_templates", "f3_pocket", "p2rank_pockets",
            "pocket_active_hits", "f4_alphafill", "alphafill_ligands", "f5_active_site",
            "active_positions", "n_pass", "priority", "source", "limitation", "error"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(json.dumps({"n": len(rows), "high": sum(r["priority"] == "high" for r in rows),
                      "review": sum(r["priority"] == "review" for r in rows), "out": str(a.out)}))


if __name__ == "__main__":
    main()
