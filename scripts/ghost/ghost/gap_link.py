#!/usr/bin/env python3
"""Vincula candidatos con brechas metabólicas (posible cierre de gap).

Para cada candidato con hit en otro genoma:
  hit (proteína del genoma donante) -> acceso Swiss-Prot -> reacciones Rhea
  -> ¿alguna reacción produce un metabolito que el genoma receptor no puede
     producir (gaps_missing_enzyme)?

Salida: TSV con candidato, genoma donante, proteína, Rhea y metabolitos
faltantes que cerraría.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from pathlib import Path


def contig_to_genome(data_root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for gff in glob.glob(str(data_root / "annotations_ab_initio" / "*" / "annotation.gff3")):
        genome = os.path.basename(os.path.dirname(gff))
        with open(gff) as fh:
            for line in fh:
                if line.startswith("##sequence-region"):
                    out[line.split()[1]] = genome
                    break
    return out


def load_rhea_map(path: Path) -> dict[str, set[str]]:
    m: dict[str, set[str]] = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            rhea, _dir, _master, acc = p[:4]
            m.setdefault(acc, set()).add(f"RHEA:{rhea}")
    return m


def sp_for_protein(data_root: Path, genome: str, protein_id: str) -> tuple[str, str]:
    """Devuelve (acceso Swiss-Prot, descripción) del hit de esa proteína."""
    path = data_root / "annotations_ab_initio" / genome / "functional.tsv"
    if not path.exists():
        return "", ""
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 13 and p[0] == protein_id:
                sseqid = p[1]
                acc = sseqid.split("|")[1] if "|" in sseqid else ""
                return acc, " ".join(p[12:])[:160]
    return "", ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    c2g = contig_to_genome(a.data_root)
    rhea_map = load_rhea_map(a.data_root / "metabolism" / "rhea2uniprot_sprot.tsv")
    trace = json.load(
        open(a.data_root / "metabolism" / f"{a.genome}_network_trace.json")
    )
    gap_by_rhea: dict[str, list[str]] = {}
    for g in trace.get("gaps_missing_enzyme", []):
        for r in g.get("producing_rhea", []):
            gap_by_rhea.setdefault(r, []).append(g["name"])

    rows = list(csv.DictReader(open(a.candidates), delimiter="\t"))
    out = []
    for r in rows:
        if not r.get("hit") or r.get("self_hit") == "True":
            continue
        subj = r["hit"]
        contig = subj.split("__")[0]
        genome = c2g.get(contig)
        if not genome:
            continue
        acc, desc = sp_for_protein(a.data_root, genome, subj)
        if not acc:
            continue
        rheas = rhea_map.get(acc, set())
        closes = sorted({g for rx in rheas for g in gap_by_rhea.get(rx, [])})
        if not closes:
            continue
        out.append({
            "start": r["start"], "end": r["end"], "strand": r["strand"],
            "aa_len": r["aa_len"], "rank": r["rank"], "ann_frac": r["ann_frac"],
            "donor_genome": genome, "donor_protein": subj,
            "swissprot": acc, "desc": desc,
            "rhea": ";".join(sorted(rheas)),
            "closes_gaps": ";".join(closes),
        })

    out.sort(key=lambda x: float(x["rank"]), reverse=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()) if out else ["start"], delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"candidatos con posible cierre de gap: {len(out)} -> {a.out}")
    for x in out[:8]:
        print(f"  {x['start']}-{x['end']}{x['strand']} rank={x['rank']} {x['swissprot']} {x['desc'][:50]} -> {x['closes_gaps'][:60]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
