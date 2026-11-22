#!/usr/bin/env python3
"""dense_synteny.py — sintenia con ortología recíproca (best hits bidireccionales).

Por qué "densa"
---------------
La sintenia anterior (build_synteny_k12.py) usaba solo el mejor hit de cada
proteína del genoma contra K-12. Eso deja fuera todo lo que no es el mejor hit y
la cobertura se queda en ~22% de los genes.

Aquí se exige **reciprocidad**: A es el mejor hit de B *y* B es el mejor hit de A.
Es el criterio estándar de ortología (RBH) y elimina parecidos espurios, pero
también es más estricto: la cobertura no sube, sube la **calidad** de la señal.

Entradas
--------
- `--forward`: hits de los 400 genomas contra K-12  (qseqid <genoma>|<tag>, sseqid b####)
- `--reverse`: hits de K-12 contra los 400 genomas  (qseqid b####, sseqid <genoma>|<tag>)

Salida
------
    data/operon_finder/k12_synteny_dense.tsv   left, right, n_genomes, frac

Uso
---
    python scripts/validation/dense_synteny.py \
        --forward .cache/synteny/hits.tsv \
        --reverse data/terminators/k12_vs_all.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))


def best_per_query(path: Path) -> dict[str, tuple[str, float]]:
    """Mejor hit global (mayor bitscore) por consulta."""
    best: dict[str, tuple[str, float]] = {}
    for line in open(path):
        f = line.split("\t")
        if len(f) < 5:
            continue
        q, s, bits = f[0], f[1], float(f[4])
        if q not in best or bits > best[q][1]:
            best[q] = (s, bits)
    return best


def best_per_query_genome(path: Path) -> dict[str, dict[str, tuple[str, float]]]:
    """Mejor hit por (consulta, genoma): lo que hace falta para sintenia."""
    best: dict[str, dict[str, tuple[str, float]]] = defaultdict(dict)
    for line in open(path):
        f = line.split("\t")
        if len(f) < 5:
            continue
        q, s, bits = f[0], f[1], float(f[4])
        genome, tag = s.split("|", 1)
        cur = best[q].get(genome)
        if cur is None or bits > cur[1]:
            best[q][genome] = (tag, bits)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--forward", required=True)
    ap.add_argument("--reverse", required=True)
    ap.add_argument("--gbff", default=str(ROOT / "data" / "test_ecoli.gbff"))
    ap.add_argument("--out", default=str(ROOT / "data" / "operon_finder"
                                        / "k12_synteny_dense.tsv"))
    a = ap.parse_args()

    fwd = best_per_query(Path(a.forward))          # genoma|tag -> b####
    rev = best_per_query_genome(Path(a.reverse))   # b#### -> {genoma: (tag, bits)}
    print(f"hits forward: {len(fwd):,} | genes K-12 con hits: {len(rev):,}")

    # recíproco: para cada gen de K-12 y cada genoma, su mejor hit es ortólogo
    # si ese hit devuelve al mismo gen de K-12 en la dirección contraria.
    rbh: dict[str, list[str]] = defaultdict(list)   # b#### -> [genoma|tag]
    for k12, per_genome in rev.items():
        for genome, (tag, _) in per_genome.items():
            hit = fwd.get(f"{genome}|{tag}")
            if hit and hit[0] == k12:
                rbh[k12].append(f"{genome}|{tag}")
    print(f"pares recíprocos: {sum(len(v) for v in rbh.values()):,}")

    # coordenadas de los genes en cada genoma
    import find_operons as F
    _, genes = F.read_gbff(Path(a.gbff))
    genes = F.ordered([g for g in genes if g.get("locus_tag")])

    coords: dict[str, dict[str, tuple]] = {}
    for gdir in sorted((ROOT / "data" / "regulation").iterdir()):
        table = gdir / "genes.tsv"
        if not table.exists():
            continue
        d = {}
        for r in csv.DictReader(open(table), delimiter="\t"):
            tag = r.get("gene_id") or r.get("locus_tag")
            if tag and r.get("start", "").isdigit():
                d[tag] = (r.get("contig", ""), int(r["start"]), int(r["end"]),
                          r.get("strand", "1"))
        coords[gdir.name] = d
    print(f"genomas con coordenadas: {len(coords)}")

    pares = [(x["locus_tag"], y["locus_tag"]) for x, y in zip(genes, genes[1:])
             if x.get("seqid") == y.get("seqid") and x.get("strand") == y.get("strand")]
    total = len(coords)
    # índice: b#### -> {genoma: tag}
    idx: dict[str, dict[str, str]] = {}
    for k12, lst in rbh.items():
        d = {}
        for item in lst:
            genome, tag = item.split("|", 1)
            d.setdefault(genome, tag)
        idx[k12] = d

    rows = []
    for left, right in pares:
        conservados = 0
        for genome, tag_l in idx.get(left, {}).items():
            tag_r = idx.get(right, {}).get(genome)
            if not tag_r:
                continue
            cl = coords.get(genome, {}).get(tag_l)
            cr = coords.get(genome, {}).get(tag_r)
            if not cl or not cr or cl[0] != cr[0] or cl[3] != cr[3]:
                continue
            gap = (cr[1] - cl[2] - 1) if cl[3] in ("1", "+") else (cl[1] - cr[2] - 1)
            if -20 <= gap <= 50:
                conservados += 1
        rows.append({"left": left, "right": right, "n_genomes": conservados,
                     "frac": round(conservados / max(1, total), 4)})

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    altos = sum(1 for r in rows if r["n_genomes"] >= 5)
    print(f"{len(rows)} pares -> {out.relative_to(ROOT)}")
    print(f"pares conservados en >=5 genomas: {altos}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
