#!/usr/bin/env python3
"""Filtro de demanda biológica: ¿la brecha está en una ruta hacia biomasa?

Calcula, sobre la red completa de Rhea, el conjunto de metabolitos que son
ancestros (precursores) de un conjunto de precursores de biomasa. Una brecha
que no es ancestro de nada esencial no es un objetivo de gen faltante.

Reutiliza el parser de Rhea del proyecto (misma versión y reglas).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import pickle
from collections import defaultdict
from pathlib import Path

# Precursores de biomasa bacteriana (nombres tal como aparecen en ChEBI).
BIOMASS_NAMES = [
    # aminoácidos
    "L-Alanine", "L-Arginine", "L-Asparagine", "L-Aspartate", "L-Cysteine",
    "L-Glutamate", "L-Glutamine", "Glycine", "L-Histidine", "L-Isoleucine",
    "L-Leucine", "L-Lysine", "L-Methionine", "L-Phenylalanine", "L-Proline",
    "L-Serine", "L-Threonine", "L-Tryptophan", "L-Tyrosine", "L-Valine",
    # nucleótidos
    "ATP", "GTP", "CTP", "UTP", "dATP", "dGTP", "dCTP", "dTTP",
    # cofactores
    "NAD", "NADP", "FAD", "FMN", "CoA", "tetrahydrofolate",
    "pyridoxal 5'-phosphate", "thiamine diphosphate", "biotin",
    "lipoic acid", "heme", "siroheme", "menaquinone", "ubiquinone",
    # lípidos / pared
    "phosphatidate", "CDP-diacylglycerol",
]


def load_name_map(data_root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in open(data_root / "reactions" / "tsv" / "chebiId_name.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            out[p[1].strip().lower()] = p[0].strip()
    return out


def load_reactions(data_root: Path, cache: Path) -> dict:
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    path = data_root.parent / "scripts" / "metabolism" / "04_network_trace.py"
    spec = importlib.util.spec_from_file_location("dna_trace", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    reactions, _eq, _master = mod.parse_rhea_rdf(
        data_root / "metabolism" / "tsv" / "rhea.rdf.gz"
    )
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_bytes(pickle.dumps(reactions))
    return reactions


def ancestors(reactions: dict, targets: set[str], max_depth: int = 8, cap: int = 20000) -> set[str]:
    """Cierra hacia atrás desde los targets (tratando Rhea como reversible)."""
    producers: dict[str, list[set[str]]] = defaultdict(list)
    for r in reactions.values():
        if r.get("transport"):
            continue
        subs, prods = r.get("substrates", set()), r.get("products", set())
        if not subs or not prods:
            continue
        for p in prods:
            producers[p].append(subs)
        for s in subs:
            producers[s].append(prods)
    seen = set(targets)
    frontier = set(targets)
    for _ in range(max_depth):
        nxt: set[str] = set()
        for m in frontier:
            for partners in producers.get(m, ()):
                for x in partners:
                    if x not in seen:
                        seen.add(x)
                        nxt.add(x)
        frontier = nxt
        if not frontier or len(seen) > cap:
            break
    return seen


def classify(data_root: Path, genome: str, cache: Path, out_dir: Path,
             max_depth: int = 8) -> None:
    name2chebi = load_name_map(data_root)
    targets = set()
    missing = []
    for n in BIOMASS_NAMES:
        c = name2chebi.get(n.lower())
        if c:
            targets.add(c)
        else:
            missing.append(n)
    reactions = load_reactions(data_root, cache)
    anc = ancestors(reactions, targets, max_depth)

    trace = json.load(open(data_root / "metabolism" / f"{genome}_network_trace.json"))
    gaps = [g for g in trace.get("gaps_missing_enzyme", [])
            if g.get("triage") == "genuine_missing_enzyme"]

    # rareza (si existe el TSV de gap_demand)
    demand_path = out_dir / f"{genome}_gap_demand.tsv"
    demand = {}
    if demand_path.exists():
        demand = {r["chebi"]: r for r in csv.DictReader(open(demand_path), delimiter="\t")}

    rows = []
    for g in gaps:
        d = demand.get(g["chebi"], {})
        on_path = g["chebi"] in anc
        specific = d.get("verdict") == "specific"
        rows.append({
            "gap": g.get("name", g["chebi"]),
            "chebi": g["chebi"],
            "on_biomass_path": on_path,
            "phylum_frac": d.get("phylum_frac", ""),
            "global_frac": d.get("global_frac", ""),
            "specific": specific,
            "priority": "alta" if (on_path and specific) else (
                "media" if on_path else "baja"),
        })
    rows.sort(key=lambda r: ({"alta": 0, "media": 1, "baja": 2}[r["priority"]],
                             -len(str(r["gap"]))))

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{genome}_gap_priority.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gap"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    counts = {p: sum(1 for r in rows if r["priority"] == p) for p in ("alta", "media", "baja")}
    print(f"{genome} | targets biomasa mapeados {len(targets)}/{len(BIOMASS_NAMES)} | ancestros {len(anc)}")
    print(f"  brechas={len(rows)} | prioridad alta={counts['alta']} media={counts['media']} baja={counts['baja']}")
    if missing:
        print(f"  sin mapear: {', '.join(missing[:6])}{'...' if len(missing) > 6 else ''}")
    for r in rows:
        if r["priority"] == "alta":
            print(f"    ALTA  {r['gap'][:45]:45s} filo={r['phylum_frac']} global={r['global_frac']}")
    print(f"  -> {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--cache", type=Path, default=None)
    ap.add_argument("--max-depth", type=int, default=8)
    a = ap.parse_args()
    cache = a.cache or (a.out_dir / "rhea_reactions.pkl")
    classify(a.data_root, a.genome, cache, a.out_dir, a.max_depth)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
