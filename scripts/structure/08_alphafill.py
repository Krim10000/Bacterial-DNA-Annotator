#!/usr/bin/env python3
"""AlphaFill: ligandos/cofactores transplantados para los UniProt esperados."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
from structure_evidence import alphafill  # noqa: E402

TARGETS = ROOT / "data" / "metabolism" / "targets"


def main():
    results = list(csv.DictReader(open(TARGETS / "NEUROSNAP_4_RESULTS.tsv"), delimiter="\t"))
    rows = []
    for r in results:
        acc = (r.get("esperado") or "").strip()
        if not acc:
            continue
        af = alphafill(acc)
        for h in (af.get("hits") or [])[:20]:
            rows.append({"genoma": r["genoma"], "uniprot": acc, "source": af["source"],
                         "pdb_id": h.get("pdb_id", ""), "ligand": h.get("ligand", ""),
                         "rmsd": h.get("rmsd", ""), "clash": h.get("clash", ""),
                         "error": af.get("error", "")})
        print(r["genoma"], acc, len(af.get("hits") or []), af.get("error", ""))
    cols = ["genoma", "uniprot", "source", "pdb_id", "ligand", "rmsd", "clash", "error"]
    with open(TARGETS / "ALPHAFILL_HITS.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print("wrote", TARGETS / "ALPHAFILL_HITS.tsv", len(rows))


if __name__ == "__main__":
    main()
