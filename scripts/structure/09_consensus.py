#!/usr/bin/env python3
"""Consenso multi-fuente: Foldseek + AlphaFill + sitios activos locales."""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
from structure_evidence import summarize  # noqa: E402

TARGETS = ROOT / "data" / "metabolism" / "targets"


def main():
    results = list(csv.DictReader(open(TARGETS / "NEUROSNAP_4_RESULTS.tsv"), delimiter="\t"))
    rows = []
    for r in results:
        pdb = Path(r["PDB"])
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        uniprot = (r.get("esperado") or "").strip()
        if not pdb.exists():
            continue
        s = summarize(r["genoma"], str(pdb), uniprot)
        s["genoma"] = s["genome"]
        s.pop("genome", None)
        s["plddt"] = r.get("pLDDT", "")
        rows.append(s)
        print(s["genoma"], s["verdict"], s["n_sources"], s["active_site_conserved"], s["error"])
    cols = ["genoma", "pdb", "uniprot", "plddt", "foldseek_source", "foldseek_top",
            "foldseek_bits", "foldseek_seq_id", "alphafill_ligands", "active_site_rows",
            "active_site_conserved", "n_sources", "verdict", "closure", "error", "limitation"]
    with open(TARGETS / "STRUCTURAL_EVIDENCE.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    md = ["# Evidencia estructural multi-fuente", "",
          "Foldseek PDB100 + AlphaFill + sitios activos locales. **No es cierre.**", "",
          "| genoma | pLDDT | Foldseek top | bits | AlphaFill (clash≤20) | active cons | fuentes | veredicto | cierre |",
          "|---|---:|---|---:|---|---:|---:|---|---|"]
    for s in rows:
        md.append(f"| {s['genoma']} | {s['plddt']} | {s['foldseek_top']} | {s['foldseek_bits']} | "
                  f"{s['alphafill_ligands'] or '—'} | {s['active_site_conserved']}/{s['active_site_rows']} | "
                  f"{s['n_sources']} | {s['verdict']} | {s['closure']} |")
    md += ["", "## Límite", "",
           "- `multi_source_support` significa que ≥2 fuentes coinciden; no prueba actividad.",
           "- `closure=no` en todos; `confirmed` requiere experimento.",
           "- AlphaFill transplantado puede tener clash; solo se listan ≤20.",
           "- Foldseek puede dar un top distinto del homólogo esperado; se conserva el top real."]
    (TARGETS / "STRUCTURAL_EVIDENCE.md").write_text("\n".join(md) + "\n")
    print("wrote", TARGETS / "STRUCTURAL_EVIDENCE.tsv")


if __name__ == "__main__":
    main()
