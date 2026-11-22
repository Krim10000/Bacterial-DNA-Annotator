#!/usr/bin/env python3
"""Tier 3: modelos alternativos para genes anotados de forma imprecisa.

Selecciona CDS con evidencia débil (baja confianza, parcial, corto, inicio
raro o varios inicios en fase) y re-enumera la ventana ±300 nt en los 6 marcos.
Conserva los ORFs que solapan el gen pero no son idénticos a la proteína
anotada: posibles extensiones, truncamientos, cambios de marco o antisentido.

Uso:
    python3 -m ghost.tier3 --data-root <data> --out out/tier3
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

from .orfs import find_orfs_oriented
from .seqio import read_fasta, revcomp, translate

STARTS = frozenset({"ATG", "GTG", "TTG"})


def parse_cds(gff: Path) -> list[dict]:
    out: list[dict] = []
    for line in open(gff):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9 or p[2] != "CDS":
            continue
        a = p[8]
        m = re.search(r"ID=([^;]+)", a)
        conf = re.search(r"conf=([\d.]+)", a)
        part = re.search(r"partial=(\d+)", a)
        st = re.search(r"start_type=(\w+)", a)
        out.append({
            "contig": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6],
            "id": m.group(1) if m else f"cds_{len(out)}",
            "conf": float(conf.group(1)) if conf else 100.0,
            "partial": part.group(1) if part else "00",
            "start_type": st.group(1) if st else "?",
        })
    return out


def alt_start_within(seq: str, start: int, end: int, strand: str, window: int = 30) -> bool:
    """¿Hay otro inicio en fase dentro de los primeros `window` nt?"""
    sub = seq[start - 1:end]
    if strand == "-":
        sub = revcomp(sub)
    head = sub[:window]
    return any(head[i:i + 3] in STARTS for i in range(3, len(head) - 2, 3))


def is_imprecise(g: dict, seq: str) -> bool:
    length = g["end"] - g["start"] + 1
    if g["partial"] != "00" or g["conf"] < 80 or length < 150:
        return True
    if g["start_type"] not in STARTS:
        return True
    return alt_start_within(seq, g["start"], g["end"], g["strand"])


def orfs_in_window(seq: str, ws: int, we: int, min_aa: int = 30) -> list[dict]:
    sub = seq[ws - 1:we]
    out: list[dict] = []
    for strand in ("+", "-"):
        oriented = sub if strand == "+" else revcomp(sub)
        for hit in find_orfs_oriented(oriented, min_aa=min_aa):
            i0, i1 = hit["start_idx"], hit["end_idx"]
            if strand == "+":
                gs, ge = ws + i0, ws + i1 - 1
            else:
                gs, ge = we - i1 + 1, we - i0
            out.append({
                "start": gs, "end": ge, "strand": strand, "aa": hit["aa"],
                "start_codon": hit["start_codon"], "partial": hit["partial"],
            })
    return out


def gene_protein(seq: str, gene: dict) -> str:
    sub = seq[gene["start"] - 1:gene["end"]]
    if gene["strand"] == "-":
        sub = revcomp(sub)
    return translate(sub, to_stop=True)


def classify(gene: dict, orf: dict) -> str:
    if orf["strand"] != gene["strand"]:
        return "antisense"
    if gene["strand"] == "+":
        delta = orf["start"] - gene["start"]
    else:
        delta = gene["end"] - orf["end"]
    if delta < 0:
        return "extended_5p"
    if delta > 0:
        return "truncated_5p" if delta % 3 == 0 else "frameshift"
    return "same_start_alt_frame"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--min-aa", type=int, default=30)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    candidates: list[dict] = []
    per_genome: dict[str, int] = {}
    genomes = sorted(glob.glob(str(a.data_root / "annotations_ab_initio" / "*")))
    genomes = genomes[a.offset:]
    if a.limit:
        genomes = genomes[: a.limit]
    print(f"genomas: {len(genomes)}")

    for gdir in genomes:
        gdir = Path(gdir)
        genome = gdir.name
        gff = gdir / "annotation.gff3"
        if not gff.exists():
            continue
        raw = sorted((a.data_root / "raw").glob(f"{genome}*.fna.gz")) + \
              sorted((a.data_root / "raw_excluded").glob(f"{genome}*.fna.gz"))
        if not raw:
            continue
        seqs = dict(read_fasta(raw[0]))
        n = 0
        for gene in parse_cds(gff):
            seq = seqs.get(gene["contig"])
            if not seq or not is_imprecise(gene, seq):
                continue
            annotated = gene_protein(seq, gene)
            ws = max(1, gene["start"] - 100)
            we = min(len(seq), gene["end"] + 100)
            for orf in orfs_in_window(seq, ws, we, a.min_aa):
                if orf["partial"]:
                    continue
                ov = max(0, min(orf["end"], gene["end"]) - max(orf["start"], gene["start"]) + 1)
                if ov / (gene["end"] - gene["start"] + 1) < 0.5:
                    continue
                if len(orf["aa"]) > (gene["end"] - gene["start"]) // 3 + 50:
                    continue
                relation = classify(gene, orf)
                if relation not in ("extended_5p", "truncated_5p", "frameshift"):
                    continue
                if orf["aa"] == annotated:
                    continue
                candidates.append({
                    "genome": genome, "contig": gene["contig"],
                    "start": orf["start"], "end": orf["end"], "strand": orf["strand"],
                    "aa": orf["aa"], "aa_len": len(orf["aa"]), "rank": 0.5,
                    "tier": 3, "source_gene": gene["id"],
                    "relation": relation,
                    "gene_conf": gene["conf"],
                })
                n += 1
        if n:
            per_genome[genome] = n

    print(f"candidatos Tier 3: {len(candidates):,} de {len(per_genome)} genomas")
    import csv
    with open(a.out / "tier3_candidates.tsv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh, delimiter="\t", extrasaction="ignore",
            fieldnames=["genome", "contig", "start", "end", "strand", "aa_len",
                        "relation", "source_gene", "gene_conf"],
        )
        w.writeheader()
        for c in candidates:
            w.writerow(c)
    (a.out / "tier3_candidates.json").write_text(
        json.dumps({"stats": {"n": len(candidates)}, "candidates": candidates})
    )
    with open(a.out / "tier3_candidates.faa", "w") as fh:
        for c in candidates:
            fh.write(f">{c['genome']}|{c['contig']}|{c['start']}|{c['end']}|{c['strand']}\n")
            for j in range(0, len(c["aa"]), 60):
                fh.write(c["aa"][j:j + 60] + "\n")

    rel: dict[str, int] = {}
    for c in candidates:
        rel[c["relation"]] = rel.get(c["relation"], 0) + 1
    summary = {
        "genomes": len(per_genome),
        "candidates": len(candidates),
        "relations": rel,
        "per_genome_max": max(per_genome.values()) if per_genome else 0,
    }
    (a.out / "tier3_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
