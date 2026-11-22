#!/usr/bin/env python3
"""Refina las brechas: solo pasos enzimáticos, no transporte y con demanda.

Una brecha "de verdad" debe cumplir:
  1. Tener al menos una reacción productora **enzimática** (con enzima Swiss-Prot)
     y **no de transporte**.
  2. El metabolito debe ser **consumido** por alguna reacción presente en el
     genoma (es un intermedio, no un callejón sin salida).
  3. No ser un ión inorgánico, ROS ni metabolito currency.
  4. Estar clasificada como `candidate_loss` por conservación de filo.

Uso:
    python3 -m ghost.gap_refine --data-root <data> --out out_fl/refined
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
from collections import Counter, defaultdict

from .labels import SOURCES, label_row
from pathlib import Path

# Iones, ROS y currency que no son objetivos de gen.
EXCLUDE = re.compile(
    r"molybdate|tungstate|superoxide|proton|water|oxygen|hydrogen peroxide|"
    r"hydroxide|metal|cobalt|nickel|zinc|\biron\b|copper|manganese|calcium|"
    r"magnesium|sodium|potassium|chloride|ammonium|ammonia|nitrate|nitrite|"
    r"sulfate|sulfite|phosphate|selenate|selenite|cadmium|mercury",
    re.I,
)

CURRENCY = {
    "CHEBI:15377",  # water
    "CHEBI:15378",  # proton
    "CHEBI:15379",  # dioxygen
    "CHEBI:16526",  # CO2
    "CHEBI:30616",  # ATP
    "CHEBI:456216",  # ADP
    "CHEBI:456215",  # AMP
    "CHEBI:57945",  # NADH
    "CHEBI:57783",  # NADPH
    "CHEBI:57287",  # CoA
    "CHEBI:58349",  # NADP
    "CHEBI:57540",  # NAD
}


def load_reactions(cache: Path) -> dict[str, dict]:
    raw = pickle.load(open(cache, "rb"))
    return {f"RHEA:{k}": v for k, v in raw.items()}


def load_rhea_to_acc(path: Path) -> dict[str, set[str]]:
    m: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3]:
                m[f"RHEA:{p[0]}"].add(p[3])
    return m


def network_reactions(path: Path) -> set[str]:
    try:
        d = json.load(open(path))
    except Exception:
        return set()
    return {
        n["id"] for n in d.get("nodes", [])
        if n.get("type") == "reaction" and str(n.get("id", "")).startswith("RHEA:")
    }


def refine(genome: str, data_root: Path, reactions: dict, r2a: dict) -> list[dict]:
    trace = json.load(open(data_root / "metabolism" / f"{genome}_network_trace.json"))
    net_path = data_root / "metabolism" / f"{genome}_network.json"
    net_rx = network_reactions(net_path)
    cons_path = Path("out_fl") / f"{genome}_gap_conservation.tsv"
    conserved = {}
    if cons_path.exists():
        conserved = {
            r["chebi"]: True for r in csv.DictReader(open(cons_path), delimiter="\t")
            if r["verdict"] == "candidate_loss"
        }

    rows = []
    for g in trace.get("gaps_missing_enzyme", []):
        if g.get("triage") != "genuine_missing_enzyme":
            continue
        chebi = g["chebi"]
        if chebi in CURRENCY or EXCLUDE.search(g.get("name", "")):
            continue
        producing = g.get("producing_rhea", [])
        enz = [r for r in producing
               if r in reactions and not reactions[r].get("transport")
               and r2a.get(r)]
        if not enz:
            continue
        consumers = [r for r in net_rx
                     if chebi in reactions.get(r, {}).get("substrates", set())]
        if not consumers:
            continue
        if conserved and chebi not in conserved:
            continue
        row = {
            "genome": genome, "gap": g.get("name", chebi), "chebi": chebi,
            "n_enz_producers": len(enz), "n_consumers": len(consumers),
            "producers": ";".join(sorted(enz)[:5]),
            "consumers": ";".join(sorted(consumers)[:5]),
            "producers_total": g.get("producers_total", 0),
        }
        label_row(
            row, label="gap_biosynthetic",
            source=" + ".join([SOURCES["rhea"], SOURCES["diamond_swissprot"],
                               SOURCES["phylum"]]),
            evidence=(f"enz_producers={len(enz)}; consumers={len(consumers)}; "
                      f"producers={row['producers']}; conserved=candidate_loss"),
            confidence="high" if len(enz) >= 2 else "medium",
            interpretation=("Brecha biosintética: existe al menos una enzima que produce "
                            "el metabolito, no es transporte y la red del genoma lo consume."),
        )
        rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cache", type=Path, default=Path("cache/rhea_reactions.pkl"))
    ap.add_argument("--genome-list", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    reactions = load_reactions(a.cache)
    r2a = load_rhea_to_acc(a.data_root / "metabolism" / "rhea2uniprot_sprot.tsv")
    if a.genome_list and a.genome_list.exists():
        genomes = json.load(open(a.genome_list))
    else:
        genomes = sorted(p.name.replace("_gap_conservation.tsv", "")
                         for p in Path("out_fl").glob("*_gap_conservation.tsv"))
    if a.limit:
        genomes = genomes[: a.limit]
    print(f"genomas: {len(genomes)} · reacciones: {len(reactions):,}")

    all_rows: list[dict] = []
    per_genome: dict[str, int] = {}
    for i, g in enumerate(genomes, 1):
        rows = refine(g, a.data_root, reactions, r2a)
        if rows:
            per_genome[g] = len(rows)
            all_rows.extend(rows)
        if i % 50 == 0:
            print(f"  {i}/{len(genomes)} · brechas refinadas: {len(all_rows):,}")

    with open(a.out / "refined_gaps.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()) if all_rows else ["genome"], delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    met = Counter(r["gap"] for r in all_rows)
    summary = {
        "genomes": len(per_genome),
        "refined_gaps": len(all_rows),
        "distinct_metabolites": len(met),
        "top": met.most_common(25),
        "per_genome_max": max(per_genome.values()) if per_genome else 0,
    }
    (a.out / "refined_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
