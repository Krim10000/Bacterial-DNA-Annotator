#!/usr/bin/env python3
"""CLI del pipeline de genes fantasma.

Ejemplo:
    python -m ghost.run --genome GCF_000006985.1_ASM698v1 \
        --data-root data \
        --out ./out --diamond swissprot --db .../data/db/swissprot.dmnd
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .pipeline import (
    add_esm_scores,
    attach_hits,
    contig_to_genome,
    load_genome_phylum,
    build_candidates,
    parse_diamond,
    run_diamond,
    write_faa,
    write_outputs,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Descubre ORFs no anotados y los prioriza.")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-region-len", type=int, default=150)
    ap.add_argument("--min-aa", type=int, default=60)
    ap.add_argument("--max-variants", type=int, default=8)
    ap.add_argument("--alt-starts", action="store_true")
    ap.add_argument("--max-regions", type=int, default=None)
    ap.add_argument("--pad", type=int, default=30)
    ap.add_argument("--exclude-partial", action="store_true")
    ap.add_argument("--diamond", choices=["none", "swissprot", "proteomes"], default="none")
    ap.add_argument("--db", type=Path, default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--esm", action="store_true", help="puntuar con ESM-2 8M (requiere .venv)")
    ap.add_argument("--esm-top", type=int, default=None)
    args = ap.parse_args(argv)

    t0 = time.time()
    res = build_candidates(
        args.data_root,
        args.genome,
        min_region_len=args.min_region_len,
        min_aa=args.min_aa,
        max_variants=args.max_variants,
        alt_starts=args.alt_starts,
        max_regions=args.max_regions,
        pad=args.pad,
        exclude_partial=args.exclude_partial,
    )
    print(
        f"[{args.genome}] contigs={res.n_contigs} bp={res.genome_bp:,} "
        f"regiones={res.n_regions} ({res.region_bp:,} bp) ORFs={res.n_orfs}"
    )

    if args.esm:
        n = add_esm_scores(res.candidates, top_n=args.esm_top)
        print(f"  ESM-2: {n} péptidos puntuados")

    if args.diamond != "none":
        if args.db is None:
            print("ERROR: --diamond requiere --db", file=sys.stderr)
            return 2
        args.out.mkdir(parents=True, exist_ok=True)
        faa = args.out / f"{args.genome}_ghost_orfs.faa"
        n = write_faa(res.candidates, faa)
        print(f"  DIAMOND: {n} péptidos -> {args.db.name}")
        raw = args.out / f"{args.genome}_diamond_{args.diamond}.tsv"
        run_diamond(faa, args.db, raw, threads=args.threads)
        hits = parse_diamond(raw)
        c2g = contig_to_genome(args.data_root) if args.diamond == "proteomes" else {}
        phy = load_genome_phylum(args.data_root) if args.diamond == "proteomes" else {}
        attach_hits(res.candidates, hits, res.contigs, c2g, phy, args.genome)
        print(f"  hits: {len(hits)}")

    paths = write_outputs(res, args.out)
    s = res.stats()
    print(
        f"  candidatos={s['n_candidates']} "
        f"(alto={sum(1 for c in res.candidates if c['rank'] >= 0.75)}, "
        f"medio={sum(1 for c in res.candidates if 0.55 <= c['rank'] < 0.75)}) "
        f"en {time.time() - t0:.1f}s"
    )
    for k, p in paths.items():
        print(f"  {k}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
