#!/usr/bin/env python3
"""Filtro de demanda: ¿qué brechas son objetivos reales de gen faltante?

Una brecha que aparece en 379/400 genomas no es una pérdida de gen: o el
metabolito no se necesita, o falta cobertura en Rhea. Una brecha que aparece
en 1–2 genomas mientras sus parientes sí producen el metabolito es el perfil
de un gen perdido/no anotado.

Este módulo mide la rareza de cada brecha (global y dentro del filo) y
clasifica las brechas de un genoma en `specific` (objetivo) o `common`.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


def _acc(name: str) -> str:
    import re

    m = re.match(r"(GC[AF]_\d+\.\d+)", name or "")
    return m.group(1) if m else (name or "")


def load_phylum(data_root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    path = data_root / "species_phylum.tsv"
    if not path.exists():
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ip = header.index("phylum")
        except ValueError:
            return out
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > ip:
                out[_acc(p[0])] = p[ip]
    return out


def build_rarity(data_root: Path) -> dict:
    """Cuenta en cuántos genomas (y en cuántos de su filo) falta cada metabolito."""
    phylum = load_phylum(data_root)
    traces = sorted(glob.glob(str(data_root / "metabolism" / "*_network_trace.json")))
    total = 0
    by_chebi: Counter = Counter()
    name: dict[str, str] = {}
    by_phylum_genomes: Counter = Counter()
    chebi_phylum: dict[str, Counter] = defaultdict(Counter)

    for path in traces:
        genome = os.path.basename(path).replace("_network_trace.json", "")
        ph = phylum.get(_acc(genome), "?")
        total += 1
        by_phylum_genomes[ph] += 1
        seen = set()
        try:
            d = json.load(open(path))
        except Exception:
            continue
        for g in d.get("gaps_missing_enzyme", []):
            if g.get("triage") != "genuine_missing_enzyme":
                continue
            chebi = g["chebi"]
            name[chebi] = g.get("name", chebi)
            seen.add(chebi)
            chebi_phylum[chebi][ph] += 1
        for chebi in seen:
            by_chebi[chebi] += 1

    return {
        "total_genomes": total,
        "by_chebi": by_chebi,
        "name": name,
        "by_phylum_genomes": by_phylum_genomes,
        "chebi_phylum": chebi_phylum,
        "phylum": phylum,
    }


def classify(data_root: Path, genome: str, out_dir: Path, max_frac: float = 0.25) -> dict:
    r = build_rarity(data_root)
    ph = r["phylum"].get(_acc(genome), "?")
    n_ph = r["by_phylum_genomes"].get(ph, 0)
    trace = json.load(open(data_root / "metabolism" / f"{genome}_network_trace.json"))
    gaps = [g for g in trace.get("gaps_missing_enzyme", [])
            if g.get("triage") == "genuine_missing_enzyme"]

    rows = []
    for g in gaps:
        chebi = g["chebi"]
        global_n = r["by_chebi"].get(chebi, 0)
        ph_n = r["chebi_phylum"][chebi].get(ph, 0)
        frac_ph = ph_n / n_ph if n_ph else 1.0
        rows.append({
            "gap": g.get("name", chebi),
            "chebi": chebi,
            "global_n": global_n,
            "global_frac": round(global_n / r["total_genomes"], 3) if r["total_genomes"] else 0,
            "phylum": ph,
            "phylum_n": ph_n,
            "phylum_frac": round(frac_ph, 3),
            "producers_total": g.get("producers_total", 0),
            "verdict": "specific" if frac_ph <= max_frac else "common",
        })
    rows.sort(key=lambda x: (x["phylum_frac"], x["global_frac"], -x["producers_total"]))

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{genome}_gap_demand.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gap"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    specific = [x for x in rows if x["verdict"] == "specific"]
    print(f"{genome} | filo={ph} ({n_ph} genomas) | brechas={len(rows)}")
    print(f"  objetivos específicos (raro en el filo): {len(specific)}")
    print(f"  comunes (probable artefacto / no requerido): {len(rows) - len(specific)}")
    for x in specific[:10]:
        print(f"    {x['gap'][:45]:45s} filo={x['phylum_n']}/{n_ph} global={x['global_n']}/{r['total_genomes']} producers={x['producers_total']}")
    print(f"  -> {out}")
    return {"rows": rows, "specific": specific, "out": out}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--max-frac", type=float, default=0.25)
    a = ap.parse_args()
    classify(a.data_root, a.genome, a.out_dir, a.max_frac)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
