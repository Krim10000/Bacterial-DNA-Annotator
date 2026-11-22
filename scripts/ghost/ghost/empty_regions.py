#!/usr/bin/env python3
"""Regiones sospechosamente vacías: desiertos de anotación.

Fusiona TODAS las señales disponibles del genoma (CDS, tRNA, rRNA, ncRNA,
oriC, transposones, defensa, ARN) y busca tramos largos sin ninguna.

Para cada tramo reporta longitud, GC, fracción de N, ORFs fantasma dentro,
su mejor rank y si alguno pertenece a un cluster conservado (≥2 genomas).

Uso:
    python3 -m ghost.empty_regions --data-root <data> --out-fl out_fl \
        --out out_fl/empty_regions --min-len 3000
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
from pathlib import Path

from .labels import SOURCES, label_row
from .orfs import merge_intervals
from .pipeline import _find_fasta
from .seqio import read_fasta, read_gff_features


def _largest_contig(contigs: dict[str, str]) -> str:
    return max(contigs, key=lambda c: len(contigs[c]))


def load_signals(genome: str, data_root: Path, contigs: dict[str, str],
                 feats: list[tuple]) -> dict[str, list[tuple[int, int]]]:
    """Todas las señales por contig."""
    sig: dict[str, list[tuple[int, int]]] = {}
    for c, s, e, _, _ in feats:
        sig.setdefault(c, []).append((s, e))
    main = _largest_contig(contigs)

    def add(s: int, e: int) -> None:
        sig.setdefault(main, []).append((s, e))

    # oriC
    arch = data_root / "architecture_annotations.json"
    if arch.exists():
        d = json.load(open(arch)).get(genome, {})
        ori = d.get("oriC")
        if ori:
            add(int(ori["start"]), int(ori["end"]))
    orca = data_root / "orca_predictions" / f"{genome}.json"
    if orca.exists():
        try:
            pos = json.load(open(orca)).get("orca_position")
            if pos:
                add(max(1, int(pos) - 500), int(pos) + 500)
        except Exception:
            pass

    # móviles y defensa
    for fname, key in (("mobile_elements.json", None), ("defense_systems.json", None)):
        p = data_root / fname
        if not p.exists():
            continue
        d = json.load(open(p)).get(genome, {})
        items = d if isinstance(d, list) else [x for v in d.values() for x in v]
        for it in items:
            if isinstance(it, dict) and "start" in it and "end" in it:
                add(int(it["start"]), int(it["end"]))

    # ARN anotado
    p = data_root / "rna_annotations" / "rna_annotations.json"
    if p.exists():
        d = json.load(open(p)).get(genome, {})
        for items in d.values():
            for it in items:
                add(int(it["start"]), int(it["end"]))
    return sig


def gc_n(seq: str) -> tuple[float, float]:
    s = seq.upper()
    if not s:
        return 0.0, 0.0
    return (s.count("G") + s.count("C")) / len(s), s.count("N") / len(s)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-fl", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-len", type=int, default=3000)
    ap.add_argument("--genome-list", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    if a.genome_list and a.genome_list.exists():
        genomes = json.load(open(a.genome_list))
    else:
        genomes = [Path(p).name.replace("_ghost_candidates.json", "")
                   for p in a.out_fl.glob("*_ghost_candidates.json")]
    if a.limit:
        genomes = genomes[: a.limit]

    # clusters conservados
    cl_ng = {r["cluster"]: int(r["n_genomes"])
             for r in csv.DictReader(open(a.out_fl / "dedup" / "dedup_clusters.tsv"), delimiter="\t")}
    locus_cluster = {}
    for r in csv.DictReader(open(a.out_fl / "dedup" / "dedup_members.tsv"), delimiter="\t"):
        locus_cluster[(r["genome"], r["locus"])] = r["cluster"]

    rows: list[dict] = []
    for i, genome in enumerate(genomes, 1):
        gff = a.data_root / "annotations_ab_initio" / genome / "annotation.gff3"
        fasta = _find_fasta(a.data_root, genome)
        if not gff.exists() or fasta is None:
            continue
        contigs = dict(read_fasta(fasta))
        _, feats = read_gff_features(gff)
        sig = load_signals(genome, a.data_root, contigs, feats)

        cand_path = a.out_fl / f"{genome}_ghost_candidates.json"
        cands = json.load(open(cand_path)).get("candidates", []) if cand_path.exists() else []

        for contig, intervals in sig.items():
            seq = contigs.get(contig, "")
            merged = merge_intervals(intervals)
            gaps = []
            cur = 1
            for s, e in merged:
                if s - cur >= a.min_len:
                    gaps.append((cur, s - 1))
                cur = max(cur, e + 1)
            if len(seq) - cur + 1 >= a.min_len:
                gaps.append((cur, len(seq)))
            for s, e in gaps:
                in_c = [c for c in cands if c["contig"] == contig
                        and c["start"] >= s and c["end"] <= e]
                g, n = gc_n(seq[s - 1:e])
                best = max(in_c, key=lambda c: c["rank"], default=None)
                conserved = 0
                for c in in_c:
                    loc = f"{c['contig']}:{c['start']}-{c['end']}{c['strand']}"
                    cl = locus_cluster.get((genome, loc))
                    if cl and cl_ng.get(cl, 1) >= 2:
                        conserved += 1
                row = {
                    "genome": genome, "contig": contig, "start": s, "end": e,
                    "length": e - s + 1, "gc": round(g, 3), "n_frac": round(n, 3),
                    "ghost_orfs": len(in_c),
                    "max_rank": round(best["rank"], 3) if best else 0.0,
                    "conserved_orfs": conserved,
                }
                if conserved > 0:
                    label, conf = "desert_candidate", "high"
                    interp = ("Hueco de anotación con ORFs presentes en ≥2 genomas: "
                              "candidato a gen/operón no anotado.")
                elif len(in_c) > 0:
                    label, conf = "desert_region", "medium"
                    interp = ("Hueco de anotación con ORFs de novo, pero sin conservación "
                              "entre genomas; candidato débil.")
                else:
                    label, conf = "desert_region", "low"
                    interp = ("Hueco sin ORFs; probablemente no codificante o hueco de "
                              "ensamblaje, no un gen escondido.")
                label_row(
                    row, label=label,
                    source=" + ".join([SOURCES["pyrodigal"], SOURCES["infernal"],
                                       SOURCES["trnascan"], SOURCES["orca"],
                                       SOURCES["mobile"], SOURCES["defense"],
                                       SOURCES["rna"], SOURCES["orf"],
                                       SOURCES["conservation"]]),
                    evidence=(f"len={row['length']}; gc={row['gc']}; n={row['n_frac']}; "
                              f"orfs={row['ghost_orfs']}; conserved={conserved}; "
                              f"max_rank={row['max_rank']}"),
                    confidence=conf, interpretation=interp,
                )
                rows.append(row)
        if i % 50 == 0:
            print(f"  {i}/{len(genomes)} · regiones: {len(rows):,}")

    rows.sort(key=lambda r: (r["conserved_orfs"] > 0, r["ghost_orfs"] > 0, r["length"]), reverse=True)
    with open(a.out / "empty_regions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["genome"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summary = {
        "genomes": len({r["genome"] for r in rows}),
        "regions": len(rows),
        "with_ghost_orfs": sum(1 for r in rows if r["ghost_orfs"] > 0),
        "with_conserved_orfs": sum(1 for r in rows if r["conserved_orfs"] > 0),
        "longest": rows[:10],
    }
    (a.out / "empty_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
