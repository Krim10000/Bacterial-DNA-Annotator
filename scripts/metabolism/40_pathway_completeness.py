#!/usr/bin/env python3
"""Completitud de rutas KEGG: reacciones esperadas vs presentes por genoma."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


def main() -> int:
    pathway_rx: dict[str, set[str]] = defaultdict(set)
    for r in csv.DictReader(open(DATA / "kegg" / "pathway_reactions.tsv"), delimiter="\t"):
        pathway_rx[f"{r['pathway_id']} {r['pathway_name']}"].add(r["rhea"])
    gap_count: dict[tuple[str, str], int] = defaultdict(int)
    essential_count: dict[tuple[str, str], int] = defaultdict(int)
    real = DATA / "metabolism" / "real_gaps.tsv"
    if real.exists():
        for r in csv.DictReader(open(real), delimiter="\t"):
            essential = (int(r.get("support_400", 0)) >= 100
                         and r.get("role") == "anabolic"
                         and r.get("transport_status") == "no_transport")
            for p in r["pathways"].split("; "):
                if p:
                    gap_count[(r["genome"], p)] += 1
                    if essential:
                        essential_count[(r["genome"], p)] += 1
    rows = []
    for net_file in sorted((DATA / "metabolism").glob("*_network.json")):
        genome = net_file.name.replace("_network.json", "")
        d = json.load(open(net_file))
        net = {l.get("target") for l in d.get("links", []) if l.get("target")}
        for pathway, expected in pathway_rx.items():
            present = len(expected & net)
            missing = len(expected - net)
            total = len(expected)
            rows.append({
                "genome": genome,
                "pathway": pathway,
                "expected": total,
                "present": present,
                "missing": missing,
                "coverage": round(present / total, 3) if total else 0,
                "gaps": gap_count.get((genome, pathway), 0),
                "essential_gaps": essential_count.get((genome, pathway), 0),
            })
    out = DATA / "metabolism" / "pathway_completeness.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # resumen por ruta
    by_path: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_path[r["pathway"]].append(r)
    summary = []
    for pathway, rs in by_path.items():
        n = len(rs)
        summary.append({
            "pathway": pathway,
            "genomes": n,
            "mean_coverage": round(sum(x["coverage"] for x in rs) / n, 3),
            "mean_missing": round(sum(x["missing"] for x in rs) / n, 1),
            "gaps": sum(x["gaps"] for x in rs),
            "essential_gaps": sum(x["essential_gaps"] for x in rs),
        })
    for r in summary:
        r["essentiality"] = ("alta" if r["essential_gaps"] >= 10
                             else "media" if r["essential_gaps"] > 0 else "baja")
    summary.sort(key=lambda x: (-x["essential_gaps"], -x["gaps"], -x["mean_missing"]))
    with open(DATA / "metabolism" / "pathway_completeness_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    lines = ["# Completitud de rutas KEGG", "",
             "| ruta | genomas | cobertura media | missing medio | gaps | gaps esenciales |",
             "|---|---|---|---|---|---|"]
    for r in summary[:25]:
        lines.append(f"| {r['pathway']} | {r['genomes']} | {r['mean_coverage']} | "
                     f"{r['mean_missing']} | {r['gaps']} | {r['essential_gaps']} |")
    lines += ["", "## Límites", "",
              "- Cobertura = reacciones Rhea de la ruta presentes en la red.",
              "- `gaps` = brechas del genoma mapeadas a esa ruta.",
              "- `essential_gaps` = soporte ≥100, anabólico y sin transporte; proxy, no FBA.",
              "- Rutas Rhea-only no aparecen en KEGG."]
    (DATA / "metabolism" / "PATHWAY_COMPLETENESS.md").write_text("\n".join(lines) + "\n")
    print(f"filas={len(rows)} rutas={len(summary)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
