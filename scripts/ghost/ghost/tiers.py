#!/usr/bin/env python3
"""Tier 2: proteínas anotadas de forma poco específica como candidatas.

No hace falta releer ADN: hay ~53.600 proteínas ya predichas cuyo mejor hit
Swiss-Prot es "hypothetical protein", "uncharacterized protein" o un dominio
DUF. Cualquiera de ellas podría ser la enzima que falta mal nombrada.

Genera un catálogo con coordenadas y secuencia, en el mismo formato que
`ghost_candidates.json` para poder pasarlo por `ghost.dedup`.

Uso:
    python3 -m ghost.tiers --data-root <data> --out out/tier2
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
from pathlib import Path

PATTERNS = re.compile(
    r"hypothetical|uncharacterized|\bDUF\d+|unknown function|putative uncharacterized",
    re.IGNORECASE,
)


def parse_functional(path: Path) -> dict[str, str]:
    """proteína -> descripción del hit, solo si es poco específica."""
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            desc = " ".join(p[12:]).replace("\t", " ")
            if PATTERNS.search(desc):
                out[p[0]] = desc[:200]
    return out


def parse_gff_coords(path: Path) -> dict[str, tuple[str, int, int, str]]:
    """ID de CDS -> (contig, inicio, fin, hebra)."""
    out: dict[str, tuple[str, int, int, str]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", p[8])
            if m:
                out[m.group(1)] = (p[0], int(p[3]), int(p[4]), p[6])
    return out


def read_selected_proteins(faa: Path, wanted: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    out[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name in wanted:
        out[name] = "".join(buf)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    candidates: list[dict] = []
    per_genome: dict[str, int] = {}
    genomes = sorted(glob.glob(str(a.data_root / "annotations_ab_initio" / "*")))
    if a.limit:
        genomes = genomes[: a.limit]
    print(f"genomas: {len(genomes)}")

    for gdir in genomes:
        gdir = Path(gdir)
        genome = gdir.name
        functional = gdir / "functional.tsv"
        faa = gdir / "proteins.faa"
        gff = gdir / "annotation.gff3"
        if not (functional.exists() and faa.exists() and gff.exists()):
            continue
        hyp = parse_functional(functional)
        if not hyp:
            continue
        coords = parse_gff_coords(gff)
        seqs = read_selected_proteins(faa, set(hyp))
        n = 0
        for pid, desc in hyp.items():
            aa = seqs.get(pid, "")
            if not aa:
                continue
            c = coords.get(pid)
            if not c:
                continue
            contig, start, end, strand = c
            candidates.append({
                "genome": genome, "contig": contig, "start": start, "end": end,
                "strand": strand, "aa": aa, "aa_len": len(aa),
                "rank": 0.5, "tier": 2, "desc": desc,
            })
            n += 1
        if n:
            per_genome[genome] = n

    print(f"candidatos Tier 2: {len(candidates):,} de {len(per_genome)} genomas")

    # TSV
    import csv
    with open(a.out / "tier2_candidates.tsv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh, delimiter="\t",
            fieldnames=["genome", "contig", "start", "end", "strand", "aa_len", "desc"],
            extrasaction="ignore",
        )
        w.writeheader()
        for c in candidates:
            w.writerow(c)

    # JSON compatible con ghost.dedup
    (a.out / "tier2_candidates.json").write_text(
        json.dumps({"stats": {"n": len(candidates)}, "candidates": candidates})
    )

    # FASTA
    with open(a.out / "tier2_candidates.faa", "w") as fh:
        for i, c in enumerate(candidates):
            fh.write(f">{c['genome']}|{c['contig']}|{c['start']}|{c['end']}|{c['strand']}\n")
            for j in range(0, len(c["aa"]), 60):
                fh.write(c["aa"][j:j + 60] + "\n")

    summary = {
        "genomes": len(per_genome),
        "candidates": len(candidates),
        "mean_len": round(sum(c["aa_len"] for c in candidates) / max(1, len(candidates)), 1),
        "per_genome_top": sorted(per_genome.items(), key=lambda x: -x[1])[:10],
        "per_genome_min": min(per_genome.values()) if per_genome else 0,
        "per_genome_max": max(per_genome.values()) if per_genome else 0,
    }
    (a.out / "tier2_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"  -> {a.out/'tier2_candidates.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
