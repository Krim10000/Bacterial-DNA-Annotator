#!/usr/bin/env python3
"""Cierre de brechas refinadas con candidatos del mismo genoma (cis).

Para cada brecha refinada:
  - enzimas objetivo = accesos Swiss-Prot de las reacciones productoras;
  - candidatos T1 del genoma con hit Swiss-Prot;
  - si el hit está entre las enzimas objetivo -> cierre candidato.

Uso:
    python3 -m ghost.gap_closure --data-root <data> --refined out_fl/refined/refined_gaps.tsv \
        --candidates-dir out_fl --out out_fl/closures.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def load_rhea_to_acc(path: Path) -> dict[str, set[str]]:
    m: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3]:
                m[f"RHEA:{p[0]}"].add(p[3])
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--refined", required=True, type=Path)
    ap.add_argument("--candidates-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    r2a = load_rhea_to_acc(a.data_root / "metabolism" / "rhea2uniprot_sprot.tsv")
    gaps_by_genome: dict[str, list[dict]] = defaultdict(list)
    for r in csv.DictReader(open(a.refined), delimiter="\t"):
        gaps_by_genome[r["genome"]].append(r)

    rows: list[dict] = []
    for genome, gaps in gaps_by_genome.items():
        cand_path = a.candidates_dir / f"{genome}_ghost_candidates.json"
        if not cand_path.exists():
            continue
        cands = json.load(open(cand_path)).get("candidates", [])
        # índice acc -> mejor candidato
        by_acc: dict[str, dict] = {}
        for c in cands:
            acc = c.get("hit")
            if not acc or c.get("self_hit"):
                continue
            if acc not in by_acc or c.get("rank", 0) > by_acc[acc].get("rank", 0):
                by_acc[acc] = c
        for g in gaps:
            producers = [x for x in g["producers"].split(";") if x]
            targets: set[str] = set()
            for r in producers:
                targets |= r2a.get(r, set())
            hits = targets & set(by_acc)
            for acc in sorted(hits):
                c = by_acc[acc]
                rows.append({
                    "genome": genome, "gap": g["gap"], "chebi": g["chebi"],
                    "n_target_accs": len(targets),
                    "candidate_locus": f"{c['contig']}:{c['start']}-{c['end']}{c['strand']}",
                    "rank": c.get("rank"), "aa_len": c.get("aa_len"),
                    "hit_acc": acc, "hit_pident": c.get("hit_pident"),
                    "hit_qcov": c.get("hit_qcov"),
                    "hit_desc": (c.get("hit_title") or "")[:140],
                })

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["genome"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"cierres cis: {len(rows)}")
    for r in rows[:15]:
        print(f"  {r['genome'][:28]} {r['gap'][:38]:38s} {r['candidate_locus']} rank={r['rank']} {r['hit_acc']} {r['hit_desc'][:40]}")
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
