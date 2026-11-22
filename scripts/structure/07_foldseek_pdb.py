#!/usr/bin/env python3
"""Foldseek API contra PDB100 para los PDB ESMFold del proyecto."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
from structure_evidence import foldseek_search  # noqa: E402

TARGETS = ROOT / "data" / "metabolism" / "targets"


def main():
    results = list(csv.DictReader(open(TARGETS / "NEUROSNAP_4_RESULTS.tsv"), delimiter="\t"))
    rows = []
    for r in results:
        pdb = Path(r["PDB"])
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        if not pdb.exists():
            continue
        fs = foldseek_search(pdb, "pdb100", 10)
        for h in (fs.get("hits") or [])[:10]:
            rows.append({"genoma": r["genoma"], "pdb_query": str(pdb), "source": fs["source"],
                         "pdb_hit": h.get("pdb_id", ""), "target": h.get("target", ""),
                         "seq_id": h.get("seq_id", ""), "bits": h.get("bits", ""),
                         "evalue": h.get("evalue", ""), "aln_len": h.get("aln_len", ""),
                         "error": fs.get("error", "")})
        print(r["genoma"], fs["source"], len(fs.get("hits") or []), fs.get("error", ""))
    cols = ["genoma", "pdb_query", "source", "pdb_hit", "target", "seq_id", "bits", "evalue", "aln_len", "error"]
    with open(TARGETS / "FOLDSEEK_PDB_HITS.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print("wrote", TARGETS / "FOLDSEEK_PDB_HITS.tsv", len(rows))


if __name__ == "__main__":
    main()
