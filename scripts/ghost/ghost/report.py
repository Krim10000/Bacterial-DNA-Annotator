#!/usr/bin/env python3
"""Genera un reporte Markdown a partir de las salidas del pipeline."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .gap_link import contig_to_genome, sp_for_protein


def _fmt(v, default="-"):
    if v is None or v == "":
        return default
    return str(v)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--genomes", nargs="+", required=True)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()

    c2g = contig_to_genome(a.data_root)
    L: list[str] = ["# Reporte de genes fantasma\n"]
    L.append(
        "Pipeline: regiones intergénicas → ORFs en 6 marcos con ambigüedad → "
        "scoring ligero (codones, RBS, composición, ESM-2) → conservación "
        "(DIAMOND vs 400 proteomas) → ranking.\n"
    )

    for g in a.genomes:
        s = json.load(open(a.out_dir / f"{g}_ghost_summary.json"))
        L.append(f"\n## {g}\n")
        L.append(
            f"- Genoma: {s['genome_bp']:,} bp · {s['n_regions']:,} regiones "
            f"intergénicas ({s['region_frac']*100:.1f}% del genoma) · "
            f"{s['n_orfs']:,} ORFs"
        )
        vpath = a.out_dir / f"{g}_validation.json"
        if vpath.exists():
            v = json.load(open(vpath))
            L.append(
                f"- Validación (100 CDS ocultados): recall **{v['recall']}**, "
                f"identidad {v['mean_identity']}, FPR@0,75 "
                f"{v['negative_control']['frac_rank_ge_075']}, FPR@0,55 "
                f"{v['negative_control']['frac_rank_ge_055']}"
            )
        L.append("\n| locus | heb. | aa | LLR | RBS | ESM | ann | rank | homólogo | filo propio |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        rows = list(csv.DictReader(open(a.out_dir / f"{g}_ghost_candidates.tsv"), delimiter="\t"))
        rows.sort(key=lambda r: float(r["rank"]), reverse=True)
        for r in rows[: a.top]:
            L.append(
                f"| {r['start']}-{r['end']} | {r['strand']} | {r['aa_len']} | "
                f"{float(r['codon_llr']):+.2f} | {float(r['rbs']):.2f} | "
                f"{_fmt(r.get('esm_llr'))} | {float(r['ann_frac']):.2f} | "
                f"{float(r['rank']):.3f} | {_fmt(r.get('hit_genome'))} | "
                f"{_fmt(r.get('same_phylum'))} |"
            )
        # función inferida del top con hit
        L.append("\nFunción inferida de los 5 mejores con homólogo:\n")
        n = 0
        for r in rows:
            if not r.get("hit") or r.get("self_hit") == "True":
                continue
            donor = c2g.get(r["hit"].split("__")[0], "")
            acc, desc = sp_for_protein(a.data_root, donor, r["hit"]) if donor else ("", "")
            L.append(
                f"- `{r['start']}-{r['end']}{r['strand']}` ({r['aa_len']} aa, "
                f"rank {float(r['rank']):.3f}): {desc or '(sin hit Swiss-Prot)'}"
            )
            n += 1
            if n >= 5:
                break

    L.append("\n## Limitaciones y siguiente paso\n")
    L.append("- Los ORFs cortos generan hits composicionales: exigir sintenía y homología de longitud completa.")
    L.append("- El cierre de brechas metabólicas requiere función (Swiss-Prot/Rhea) o estructura (Foldseek).")
    L.append("- Solo el top pasa a plegamiento (ESMFold/AlphaFold).")
    a.report.write_text("\n".join(L) + "\n")
    print(f"reporte -> {a.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
