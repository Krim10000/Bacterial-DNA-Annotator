#!/usr/bin/env python3
"""Escanea PDBs con EnzyMM/PyJess y escribe un TSV de plantillas M-CSA."""
import argparse, csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
from structural_annotation import scan_pdbs, summarize  # noqa: E402

TARGETS = ROOT / "data" / "metabolism" / "targets"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=TARGETS / "ENZYMM_SCAN.tsv")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    res = scan_pdbs(a.pdb, a.out.parent, jobs=a.jobs)
    rows = summarize(res["hits"])
    cols = ["query_id", "template_pdb_id", "template_mcsa_id", "template_uniprot_id",
            "template_ec", "rmsd", "log_evalue", "matched_residues",
            "template_effective_size", "predicted_correct", "source", "limitation"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print({"hits": len(rows), "no_template": len(res["no_template"]), "error": res["error"], "out": str(a.out)})


if __name__ == "__main__":
    main()
