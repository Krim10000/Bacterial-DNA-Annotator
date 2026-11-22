#!/usr/bin/env python3
"""Cierre de brechas dirigido por enzima.

Para cada metabolito que el genoma no puede producir (brecha "genuine"):
  1. Toma las reacciones Rhea que podrían producirlo.
  2. Expande a los accesos Swiss-Prot que catalizan esas reacciones.
  3. Comprueba si el genoma ya tiene alguno anotado (brecha de mapeo) o no
     (brecha de gen).
  4. Busca en los ORFs fantasma un hit Swiss-Prot que sea uno de esos accesos.

Salida: TSV por brecha + resumen.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
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


def load_annotated_accs(functional_tsv: Path) -> set[str]:
    accs: set[str] = set()
    if not functional_tsv.exists():
        return accs
    with open(functional_tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and "|" in p[1]:
                accs.add(p[1].split("|")[1])
    return accs


def run_diamond(query: Path, db: Path, out: Path, threads: int = 4) -> None:
    subprocess.run(
        [
            "diamond", "blastp", "--query", str(query), "--db", str(db),
            "--out", str(out), "--outfmt", "6", "qseqid", "sseqid", "pident",
            "length", "qlen", "evalue", "bitscore",
            "--evalue", "1e-5", "--max-target-seqs", "5",
            "--threads", str(threads), "--quiet",
        ],
        check=True,
    )


def parse_hits(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            acc = p[1].split("|")[1] if "|" in p[1] else p[1]
            hits[p[0]].append({
                "acc": acc, "pident": float(p[2]),
                "qcov": round(int(p[3]) / max(1, int(p[4])), 3),
                "evalue": float(p[5]),
            })
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--demand-tsv", type=Path, default=None,
                    help="TSV de gap_demand para restringir a brechas 'specific'")
    ap.add_argument("--conservation-tsv", type=Path, default=None,
                    help="TSV de gap_conservation para restringir a 'candidate_loss'")
    a = ap.parse_args()

    trace = json.load(open(a.data_root / "metabolism" / f"{a.genome}_network_trace.json"))
    gaps = [g for g in trace.get("gaps_missing_enzyme", [])
            if g.get("triage") == "genuine_missing_enzyme"]
    r2a = load_rhea_to_acc(a.data_root / "metabolism" / "rhea2uniprot_sprot.tsv")
    ann = load_annotated_accs(
        a.data_root / "annotations_ab_initio" / a.genome / "functional.tsv"
    )

    faa = a.out_dir / f"{a.genome}_ghost_orfs.faa"
    hits_path = a.out_dir / f"{a.genome}_ghost_vs_swissprot.tsv"
    if not hits_path.exists() and faa.exists():
        run_diamond(faa, a.db, hits_path, a.threads)
    hits = parse_hits(hits_path) if hits_path.exists() else {}
    by_acc: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for q, hs in hits.items():
        for h in hs:
            by_acc[h["acc"]].append((q, h))

    # coordenadas de cada query ORF
    coord: dict[str, dict] = {}
    cand_tsv = a.out_dir / f"{a.genome}_ghost_candidates.tsv"
    if cand_tsv.exists():
        for r in csv.DictReader(open(cand_tsv), delimiter="\t"):
            qid = f"{a.genome}|{r['contig']}|{r['start']}|{r['end']}|{r['strand']}"
            coord[qid] = r

    only_specific: set[str] | None = None
    if a.conservation_tsv and a.conservation_tsv.exists():
        only_specific = {
            r["chebi"] for r in csv.DictReader(open(a.conservation_tsv), delimiter="\t")
            if r.get("verdict") == "candidate_loss"
        }
    elif a.demand_tsv and a.demand_tsv.exists():
        only_specific = {
            r["chebi"] for r in csv.DictReader(open(a.demand_tsv), delimiter="\t")
            if r.get("verdict") == "specific"
        }

    rows = []
    for g in gaps:
        if only_specific is not None and g["chebi"] not in only_specific:
            continue
        targets: set[str] = set()
        for r in g.get("producing_rhea", []):
            targets |= r2a.get(r, set())
        if not targets:
            continue
        annotated = targets & ann
        closers = []
        for acc in targets:
            for q, h in by_acc.get(acc, []):
                c = coord.get(q)
                closers.append({
                    "acc": acc,
                    "locus": f"{c['start']}-{c['end']}{c['strand']}" if c else q,
                    "pident": h["pident"],
                    "qcov": h["qcov"],
                    "evalue": h["evalue"],
                })
        rows.append({
            "gap": g["name"],
            "chebi": g["chebi"],
            "producers_total": g.get("producers_total"),
            "n_target_enzymes": len(targets),
            "n_already_annotated": len(annotated),
            "n_candidate_closers": len(closers),
            "closers": ";".join(
                f"{x['locus']}={x['acc']}({x['pident']:.0f}%,{x['qcov']:.2f})"
                for x in closers[:5]
            ),
        })

    a.out_dir.mkdir(parents=True, exist_ok=True)
    out = a.out_dir / f"{a.genome}_gap_targets.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["gap"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    mapping_gaps = sum(1 for r in rows if r["n_already_annotated"] > 0)
    gene_gaps = sum(1 for r in rows if r["n_already_annotated"] == 0)
    with_closers = sum(1 for r in rows if r["n_candidate_closers"] > 0)
    print(f"brechas genuine analizadas: {len(rows)}")
    print(f"  ya tienen enzima anotada (brecha de mapeo): {mapping_gaps}")
    print(f"  sin enzima anotada (brecha de gen):        {gene_gaps}")
    print(f"  con candidato fantasma que las cierra:     {with_closers}")
    for r in sorted(rows, key=lambda x: x["n_candidate_closers"], reverse=True)[:6]:
        print(f"  {r['gap'][:42]:42s} targets={r['n_target_enzymes']:>3} anot={r['n_already_annotated']:>2} closers={r['n_candidate_closers']} {r['closers'][:50]}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
