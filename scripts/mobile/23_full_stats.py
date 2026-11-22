#!/usr/bin/env python3
"""23_full_stats.py — tablas y cifras para `README_full.md`.

Lee `summary_full.json`, `full_genomes.tsv` y los intermedios de
`_work/full/` y escribe un informe de texto (Markdown) por stdout. No modifica
nada: es solo para construir el README y revisar cobertura.

Uso: python scripts/mobile/23_full_stats.py > /tmp/full_stats.md
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "mobile_typing"
WORK = BASE / "_work" / "full"
TYPES = ["transposon", "virulence", "cas", "defense", "crispr", "amr", "prophage"]


def main() -> int:
    summary = json.loads((BASE / "summary_full.json").read_text())
    meta = {r["genome"]: r for r in csv.DictReader((BASE / "full_genomes.tsv").open(), delimiter="\t")}
    counts = summary["counts"]
    cov = summary.get("coverage", {})
    print("## Cobertura\n")
    print(f"- genomas esperados: {cov.get('genomes_expected')}")
    print(f"- procesados (TSV): {cov.get('genomes_processed')}")
    print(f"- fallidos: {cov.get('genomes_failed')}")
    print(f"- con FASTA: {cov.get('genomes_with_fasta')} · sin FASTA: {cov.get('genomes_without_fasta')}")
    print(f"- con Pfam: {cov.get('pfam_scanned')} · con CRISPR: {cov.get('crispr_scanned')}")
    print(f"- filas: {counts['rows_total']} · arrays: {counts['crispr_arrays']} · espaciadores: {counts['crispr_spacers']}")
    print("\n## Tipos\n")
    print("| tipo | filas | high | medium | low |")
    print("|---|---:|---:|---:|---:|")
    conf_by_type = defaultdict(Counter)
    for genome in summary.get("pilot_genomes", []):
        if genome not in summary["counts"]["by_genome"]:
            continue
        path = BASE / "full" / f"{genome}.tsv"
        if not path.exists():
            continue
        for row in csv.DictReader(path.open(), delimiter="\t"):
            conf_by_type[row["element_type"]][row["confidence"]] += 1
    for t in TYPES:
        c = conf_by_type[t]
        print(f"| {t} | {sum(c.values())} | {c['high']} | {c['medium']} | {c['low']} |")
    print("\n## Por filo\n")
    phy = defaultdict(Counter)
    processed = set(summary["counts"]["by_genome"])
    for genome in processed:
        p = meta.get(genome, {}).get("phylum", "unknown")
        phy[p]["genomas"] += 1
        phy[p]["filas"] += summary["counts"]["by_genome"][genome]["rows"]
        for t, n in summary["counts"]["by_genome"][genome]["by_type"].items():
            phy[p][t] += n
    print("| filo | genomas | transposon | cas | defense | amr | virulence | crispr | profago | filas |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for p in sorted(phy, key=lambda x: -phy[x]["filas"]):
        c = phy[p]
        print(f"| {p} | {c['genomas']} | {c['transposon']} | {c['cas']} | {c['defense']} | "
              f"{c['amr']} | {c['virulence']} | {c['crispr']} | {c['prophage']} | {c['filas']} |")
    print("\n## Subtipos más frecuentes\n")
    print("| tipo:subtipo | filas |")
    print("|---|---:|")
    for key, n in sorted(counts["by_subtype"].items(), key=lambda kv: -kv[1])[:30]:
        print(f"| {key} | {n} |")
    print("\n## IS (subtipos)\n")
    is_counts = {k.split(":", 1)[1]: v for k, v in counts["by_subtype"].items()
                 if k.startswith("transposon:IS")}
    print(", ".join(f"{k} {v}" for k, v in sorted(is_counts.items(), key=lambda kv: -kv[1])[:15]))
    print("\n## Demisiones / reclasificaciones\n")
    ev = summary["evidence"]
    print(f"- previos: {ev['prior_items_total']} · mapeados ab initio: {ev['prior_mapped_to_ab_initio']} "
          f"· ncbi_only: {ev['prior_ncbi_only']} · fuera de contig: {ev['prior_out_of_contig']}")
    print(f"- descartados/reclasificados: {ev['prior_demoted_or_reclassified']}")
    for section, n in sorted(ev.get("prior_demoted_by_section", {}).items(), key=lambda kv: -kv[1]):
        total = ev.get("prior_by_section", {}).get(section, 0)
        print(f"  - {section}: {n}/{total} ({100*n/max(total,1):.1f} %)")
    print("\n## Fallos\n")
    for genome, reason in sorted(cov.get("failures", {}).items()):
        print(f"- `{genome}`: {reason}")
    print("\n## Sin Pfam (primeros 20)\n")
    print(", ".join(cov.get("without_pfam", [])[:20]) or "ninguno")
    print("\n## Sin CRISPR (sin FASTA)\n")
    print(", ".join(cov.get("without_crispr", [])[:20]) or "ninguno")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
