#!/usr/bin/env python3
"""Conservación de la reacción: ¿el filo puede producir el metabolito?

Para cada brecha de un genoma, mira si algún genoma del mismo filo tiene en su
red alguna de las reacciones Rhea que producen ese metabolito.

  - Si ≥1 pariente la tiene y este genoma no -> pérdida de gen (objetivo real).
  - Si ningún pariente la tiene -> el metabolito no forma parte del metabolismo
    del filo (artefacto de Rhea / química no bacteriana).

Esto separa los alcaloides vegetales de las enzimas realmente perdidas.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
from collections import Counter
from pathlib import Path


def _acc(name: str) -> str:
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


def network_reactions(path: Path) -> set[str]:
    try:
        d = json.load(open(path))
    except Exception:
        return set()
    return {
        n["id"] for n in d.get("nodes", [])
        if n.get("type") == "reaction" and str(n.get("id", "")).startswith("RHEA:")
    }


def genome_of(path: Path) -> str:
    return os.path.basename(path).replace("_network.json", "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--max-relatives", type=int, default=40)
    a = ap.parse_args()

    phylum = load_phylum(a.data_root)
    my_ph = phylum.get(_acc(a.genome), "?")
    nets = sorted(glob.glob(str(a.data_root / "metabolism" / "*_network.json")))
    relatives = [p for p in nets
                 if genome_of(p) != a.genome
                 and phylum.get(_acc(genome_of(p)), "?") == my_ph][: a.max_relatives]
    rel_rx = {genome_of(p): network_reactions(Path(p)) for p in relatives}

    trace = json.load(open(a.data_root / "metabolism" / f"{a.genome}_network_trace.json"))
    gaps = [g for g in trace.get("gaps_missing_enzyme", [])
            if g.get("triage") == "genuine_missing_enzyme"]

    rows = []
    for g in gaps:
        producing = set(g.get("producing_rhea", []))
        with_rx = [name for name, rx in rel_rx.items() if producing & rx]
        rows.append({
            "gap": g.get("name", g["chebi"]),
            "chebi": g["chebi"],
            "producers_total": g.get("producers_total", 0),
            "n_relatives": len(rel_rx),
            "n_relatives_with_reaction": len(with_rx),
            "verdict": "candidate_loss" if with_rx else "phylum_absent",
            "examples": ";".join(with_rx[:3]),
        })
    rows.sort(key=lambda r: (r["verdict"] != "candidate_loss",
                             -r["n_relatives_with_reaction"],
                             -r["producers_total"]))

    a.out_dir.mkdir(parents=True, exist_ok=True)
    out = a.out_dir / f"{a.genome}_gap_conservation.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gap"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    loss = [r for r in rows if r["verdict"] == "candidate_loss"]
    absent = [r for r in rows if r["verdict"] == "phylum_absent"]
    print(f"{a.genome} | filo={my_ph} | parientes={len(rel_rx)} | brechas={len(rows)}")
    print(f"  pérdida de gen (el filo sí puede): {len(loss)}")
    print(f"  ausente en el filo (no objetivo):  {len(absent)}")
    for r in loss[:12]:
        print(f"    {r['gap'][:44]:44s} parientes_con_reaccion={r['n_relatives_with_reaction']}/{r['n_relatives']}")
    if absent:
        print("  descartadas:", ", ".join(r["gap"][:22] for r in absent[:6]))
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
