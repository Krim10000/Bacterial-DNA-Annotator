#!/usr/bin/env python3
"""build_synteny_k12.py — ¿cuánto se conserva cada par de genes de K-12?

Para cada par de genes contiguos de *E. coli* K-12 se pregunta: en cuántos de
los 400 genomas del proyecto sus ortólogos siguen **juntos y en el mismo orden**.

Eso es sintenia, y es una señal que no mira el hueco ni el promotor: mide la
conservación a lo largo de la evolución. Sirve justo donde el hueco no decide
(la banda ambigua de 43–100 pb).

Método
------
1. Proteínas de K-12 desde el GBFF (id = locus_tag) → base DIAMOND.
2. Proteínas de los 400 genomas → consultas, con id ``<genoma>|<locus_tag>``.
3. Mejor hit de cada consulta contra K-12 (``--max-target-seqs 1``).
4. Para cada par contiguo de K-12, en cada genoma: ¿hay hits a ambos genes
   adyacentes (mismo contig, misma hebra, hueco <= 50 pb)?

Salida: ``data/operon_finder/k12_synteny.tsv`` (left, right, n_genomes, frac).

Uso
---
    python scripts/validation/build_synteny_k12.py [--threads 8] [--limit 40]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import find_operons as F  # noqa: E402

GBFF = ROOT / "data" / "test_ecoli.gbff"
REG = ROOT / "data" / "regulation"
WORK = ROOT / ".cache" / "synteny"
OUT = ROOT / "data" / "operon_finder" / "k12_synteny.tsv"


def k12_genes() -> list[dict]:
    _, genes = F.read_gbff(GBFF)
    return F.ordered([g for g in genes if g.get("locus_tag")])


def write_k12_proteins(genes: list[dict], path: Path) -> None:
    """Traducciones del GBFF con Biopython (locus_tag como id)."""
    from Bio import SeqIO
    want = {g["locus_tag"] for g in genes}
    n = 0
    with open(path, "w") as fh:
        for rec in SeqIO.parse(str(GBFF), "genbank"):
            for f in rec.features:
                if f.type != "CDS":
                    continue
                tag = (f.qualifiers.get("locus_tag") or [""])[0]
                tr = (f.qualifiers.get("translation") or [""])[0]
                if tag in want and tr:
                    fh.write(f">{tag}\n{tr}\n")
                    n += 1
    print(f"K-12 proteínas escritas: {n}")


def collect_queries(limit: int | None) -> tuple[Path, dict]:
    """Escribe las proteínas de los genomas y devuelve sus coordenadas."""
    genomes = sorted(d for d in REG.iterdir()
                     if d.is_dir() and (d / "genes.tsv").exists())
    if limit:
        genomes = genomes[:limit]
    coords: dict[str, list[dict]] = {}
    path = WORK / "queries.faa"
    n = 0
    with open(path, "w") as out:
        for gdir in genomes:
            gname = gdir.name
            table = [r for r in csv.DictReader(open(gdir / "genes.tsv"), delimiter="\t")
                     if (r.get("gene_id") or r.get("locus_tag"))
                     and r.get("start", "").isdigit()]
            if not table:
                continue
            coords[gname] = [{"tag": r.get("gene_id") or r["locus_tag"],
                              "contig": r.get("contig", ""),
                              "start": int(r["start"]), "end": int(r["end"]),
                              "strand": r.get("strand", "1")} for r in table]
            faa = ROOT / "data" / "annotations_ab_initio" / gname / "proteins.faa"
            if not faa.exists():
                continue
            wrote = False
            for line in open(faa):
                if line.startswith(">"):
                    tag = line[1:].strip().split()[0]
                    out.write(f">{gname}|{tag}\n")
                    n += 1
                    wrote = True
                elif wrote:
                    out.write(line)
    print(f"consultas: {n:,} proteínas de {len(coords)} genomas -> {path}")
    return path, coords


def run_diamond(queries: Path, k12db: Path, threads: int) -> Path:
    hits = WORK / "hits.tsv"
    if hits.exists() and hits.stat().st_size > 0:
        print(f"hits ya calculados: {hits}")
        return hits
    cmd = ["diamond", "blastp", "-q", str(queries), "-d", str(k12db),
           "-o", str(hits), "-f", "6", "qseqid", "sseqid", "pident", "evalue",
           "bitscore", "--max-target-seqs", "1", "--fast", "--threads", str(threads),
           "--quiet"]
    print("diamond:", " ".join(cmd[:6]), "…")
    subprocess.run(cmd, check=True)
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None,
                    help="usar solo los primeros N genomas (prueba rápida)")
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)

    genes = k12_genes()
    print(f"K-12: {len(genes)} CDS")
    k12fa = WORK / "k12.faa"
    write_k12_proteins(genes, k12fa)
    db = WORK / "k12db"
    if not Path(str(db) + ".dmnd").exists():
        subprocess.run(["diamond", "makedb", "--in", str(k12fa), "-d", str(db)],
                       check=True)
        print("base DIAMOND creada")

    queries, coords = collect_queries(a.limit)
    hits = run_diamond(queries, db, a.threads)

    # hits por genoma y por gen de K-12
    by_genome: dict[str, dict] = defaultdict(lambda: defaultdict(list))
    n_hits = 0
    for line in open(hits):
        q, s = line.split("\t")[:2]
        if "|" not in q:
            continue
        gname, tag = q.split("|", 1)
        by_genome[gname][s].append(tag)
        n_hits += 1
    print(f"hits: {n_hits:,} en {len(by_genome)} genomas")

    # índice de coordenadas por genoma
    idx = {}
    for gname, table in coords.items():
        idx[gname] = {r["tag"]: (r["contig"], r["start"], r["end"], r["strand"])
                      for r in table}

    pares = [(g["locus_tag"], h["locus_tag"]) for g, h in zip(genes, genes[1:])
             if g.get("seqid") == h.get("seqid") and g.get("strand") == h.get("strand")]

    total = len(by_genome)
    rows = []
    for left, right in pares:
        conservados = 0
        for gname, hitmap in by_genome.items():
            coords_g = idx.get(gname, {})
            if not coords_g:
                continue
            for ql in hitmap.get(left, []):
                cl = coords_g.get(ql)
                if not cl:
                    continue
                for qr in hitmap.get(right, []):
                    cr = coords_g.get(qr)
                    if not cr or cl[0] != cr[0] or cl[3] != cr[3]:
                        continue
                    gap = cr[1] - cl[2] - 1 if cl[3] in ("1", "+") else cl[1] - cr[2] - 1
                    if -20 <= gap <= 50:
                        conservados += 1
                        break
                else:
                    continue
                break
        rows.append({"left": left, "right": right, "n_genomes": conservados,
                     "frac": round(conservados / max(1, total), 4)})

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["left", "right", "n_genomes", "frac"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"pares de K-12: {len(rows)} -> {OUT.relative_to(ROOT)}")
    altos = sum(1 for r in rows if r["frac"] >= 0.8)
    print(f"pares conservados en >=80% de los {total} genomas: {altos}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
