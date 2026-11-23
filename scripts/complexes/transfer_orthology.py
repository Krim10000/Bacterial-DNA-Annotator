#!/usr/bin/env python3
"""Transfiere complejos K-12 a otros genomas por ortología recíproca.

Insumos:
  data/complexes/gold_k12.tsv          catálogo curado (273 complejos)
  data/terminators/k12_vs_all.tsv      K-12 -> genomas (b-number, target, ...)
  .cache/synteny/hits.tsv              genomas -> K-12 (target, b-number, ...)

No hay cobertura de alineamiento en los TSV actuales. Se aplica RBH por
bitscore + identidad ≥30 % + evalue ≤1e-5, y se declara esa limitación.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CX = DATA / "complexes"


def read_hits(path: Path, target_first: bool = True) -> dict[tuple[str, str], tuple[float, float, float]]:
    """(target_key, b_number) -> (pident, evalue, bitscore).

    ``k12_vs_all.tsv`` empieza por b-number; ``hits.tsv`` por target.
    """
    out: dict[tuple[str, str], tuple[float, float, float]] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        c = line.split("\t")
        if len(c) < 5:
            continue
        try:
            target, b = (c[0], c[1]) if target_first else (c[1], c[0])
            out[(target, b)] = (float(c[2]), float(c[3]), float(c[4]))
        except ValueError:
            continue
    return out


def rbh(k12: dict, target: dict) -> dict[tuple[str, str], str]:
    """(genome, b-number) -> target_key con best hit recíproco por genoma.

    El RBH se calcula dentro de cada genoma: un mismo gen K-12 puede tener su
    mejor ortólogo en muchos genomas distintos.
    """
    best_by_genome: dict[str, dict[str, tuple[str, float]]] = {}
    for (t, b), (pid, ev, bits) in k12.items():
        if pid < 30 or ev > 1e-5:
            continue
        genome = t.split("|", 1)[0]
        cur = best_by_genome.setdefault(genome, {})
        if b not in cur or bits > cur[b][1]:
            cur[b] = (t, bits)
    # reverse best b para cada target
    rev_best: dict[str, tuple[str, float]] = {}
    for (t, b), (_pid, ev, bits) in target.items():
        if ev > 1e-5:
            continue
        if t not in rev_best or bits > rev_best[t][1]:
            rev_best[t] = (b, bits)
    out: dict[tuple[str, str], str] = {}
    for genome, bs in best_by_genome.items():
        for b, (t, _bits) in bs.items():
            if rev_best.get(t, ("", 0))[0] == b:
                out[(genome, b)] = t
    return out



def expected_n(props: str) -> int | None:
    m = re.search(r"Hetero\s*(\d+)-mer", props or "", re.I)
    return int(m.group(1)) if m else None


def complex_map() -> dict[str, dict]:
    out = {}
    with open(CX / "gold_k12.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tags = [x for x in (r.get("locus_tags") or "").split(",") if x]
            if len(tags) < 2:
                continue
            out[r["complex_id"]] = {
                "complex_id": r["complex_id"], "complex_name": r["complex_name"],
                "locus_tags": tags, "gene_names": r.get("gene_names", ""),
                "expected_n": expected_n(r.get("properties", "")),
                "evidence_class": r.get("evidence_class", ""),
                "properties": r.get("properties", ""),
            }
    return out


def read_lengths(genome: str) -> dict[str, int]:
    path = DATA / "annotations_ab_initio" / genome / "proteins.faa"
    out = {}
    if not path.exists():
        return out
    name = None
    n = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                out[name] = n
            name = line[1:].split()[0]
            n = 0
        else:
            n += len(line.strip())
    if name:
        out[name] = n
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k12", type=Path, default=DATA / "terminators" / "k12_vs_all.tsv")
    ap.add_argument("--target", type=Path, default=ROOT / ".cache" / "synteny" / "hits.tsv")
    ap.add_argument("--out", type=Path, default=CX)
    ap.add_argument("--summary", type=Path,
                    default=ROOT / "reports" / "complex_transfer_summary.tsv")
    args = ap.parse_args()

    k12 = read_hits(args.k12, target_first=False)
    target = read_hits(args.target, target_first=True)
    gold = complex_map()
    # EC de los genes K-12 para poder asociar transferencias a reacciones.
    ec_by_b: dict[str, set[str]] = defaultdict(set)
    k12_genes = DATA / "regulation" / "ecoli_k12" / "genes.tsv"
    if k12_genes.exists():
        with k12_genes.open(encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                b = (r.get("locus_tag") or "").strip()
                for ec in re.split(r"[;|,]", r.get("ec") or ""):
                    if ec.strip():
                        ec_by_b[b].add(ec.strip())
    mapped = rbh(k12, target)
    genomes = sorted({g for g, _b in mapped})
    print(f"RBH: {len(mapped)} pares; {len(genomes)} genomas; {len(gold)} complejos")

    # presencia por b-number para perfil filogenético
    presence: dict[str, set[str]] = defaultdict(set)
    for (t, b), (_p, ev, _bits) in k12.items():
        if ev <= 1e-5:
            presence[b].add(t.split("|", 1)[0])

    def jaccard(a: str, b: str) -> float | None:
        pa, pb = presence.get(a, set()), presence.get(b, set())
        if not pa and not pb:
            return None
        return len(pa & pb) / len(pa | pb) if (pa | pb) else None

    args.out.mkdir(parents=True, exist_ok=True)
    summary = []
    for g in genomes:
        b2t = {b: t for (gg, b), t in mapped.items() if gg == g}
        if not b2t:
            continue
        lens = read_lengths(g)
        rows = []
        for cid, c in gold.items():
            found = [b for b in c["locus_tags"] if b in b2t]
            if len(found) < 2:
                continue
            ids = [b2t[b] for b in found]
            ids = [x.split("|", 1)[1] for x in ids]
            identities = [k12.get((b2t[b], b), (0, 1, 0))[0] for b in found]
            n_exp = c["expected_n"]
            mass = sum(lens.get(x, 0) for x in ids) * 110 / 1000
            js = [j for a, b in __import__("itertools").combinations(c["locus_tags"], 2)
                  if (j := jaccard(a, b)) is not None]
            rows.append({
                "complex_id": cid, "complex_name": c["complex_name"],
                "n_subunits": len(c["locus_tags"]), "n_mapped": len(found),
                "frac_mapped": round(len(found) / len(c["locus_tags"]), 3),
                "identity_mean": round(sum(identities) / len(identities), 1),
                "score": round(len(found) / len(c["locus_tags"])
                               * sum(identities) / len(identities) / 100, 4),
                "genes_k12": ",".join(found), "genes_genoma": ",".join(ids),
                "observed_n": len(ids), "expected_n": n_exp or "",
                "stoich_status": ("complete" if n_exp and len(ids) >= n_exp
                                  else "incomplete" if n_exp else "unknown"),
                "mass_kda": round(mass, 1),
                "phylo_jacc": round(sum(js) / len(js), 3) if js else "",
                "ecs": ";".join(sorted({e for b in found for e in ec_by_b.get(b, set())})),
                "evidence_class": "documentado_inferido",
                "eco": "ECO:0005544",
                "limitation": "RBH sin cobertura de alineamiento en los TSV actuales.",
            })
        path = args.out / f"transfer_{g}.tsv"
        cols = list(rows[0]) if rows else ["complex_id"]
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
            w.writeheader(); w.writerows(rows)
        summary.append({"genome": g, "complexes": len(rows),
                        "high": sum(1 for r in rows if r["frac_mapped"] >= 0.8),
                        "complete": sum(1 for r in rows if r["stoich_status"] == "complete"),
                        "with_phylo": sum(1 for r in rows if r["phylo_jacc"] != "")})
        print(f"  {g}: {len(rows)} complejos", flush=True)

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="", encoding="utf-8") as fh:
        cols = list(summary[0]) if summary else ["genome"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(summary)
    (args.summary.with_suffix(".json")).write_text(
        json.dumps({"genomes": len(summary),
                    "complexes": sum(r["complexes"] for r in summary),
                    "high": sum(r["high"] for r in summary),
                    "complete": sum(r["complete"] for r in summary)},
                   indent=2), encoding="utf-8")
    print("listo:", args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
