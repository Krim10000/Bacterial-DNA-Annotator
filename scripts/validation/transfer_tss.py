#!/usr/bin/env python3
"""transfer_tss.py — transferir TSS por ortología a un genoma sin anotación.

Para qué
--------
Los 400 genomas del proyecto no tienen TSS medidos. Si un gen tiene ortólogo en
K-12 (o en otra especie documentada), su TSS se puede **transferir** colocándolo
a la misma distancia del inicio del gen.

Se valida en *B. subtilis*, donde sí hay verdad: si la transferencia desde K-12
se acerca al resultado con TSS propios, sirve como puente para los genomas sin
documentación; si no, no sirve y hay que decirlo.

Método
------
1. Traducir los CDS del genoma destino (el GBFF de BSGatlas no trae traducciones).
2. DIAMOND contra la base de K-12; mejor hit por proteína.
3. Para cada gen de K-12 con TSS aguas arriba, calcular el desfase
   `offset = inicio_gen − TSS`.
4. Transferir ese desfase al ortólogo del genoma destino.
5. Segmentar por los TSS transferidos y medir contra el ground truth.

Uso
---
    python scripts/validation/transfer_tss.py
"""
from __future__ import annotations

import csv
import subprocess
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import bisect

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import find_operons as F  # noqa: E402

GBFF = ROOT / "data" / "dbtbs" / "subtilis.gbff"
FAA = ROOT / "data" / "dbtbs" / "subtilis.faa"
K12DB = ROOT / ".cache" / "synteny" / "k12db"
OUT = ROOT / "data" / "dbtbs" / "bsu_transferred_tss.tsv"


def translate_genome() -> dict[str, dict]:
    """Traduce los CDS y devuelve sus metadatos."""
    from Bio import SeqIO
    rec = next(SeqIO.parse(str(GBFF), "genbank"))
    seq = rec.seq
    genes = {}
    with open(FAA, "w") as fh:
        for f in rec.features:
            if f.type != "CDS":
                continue
            tag = (f.qualifiers.get("locus_tag") or [""])[0]
            if not tag:
                continue
            sub = f.extract(seq)
            if f.location.strand == -1:
                sub = sub.reverse_complement()
            prot = str(sub.translate(to_stop=False))
            if len(prot) < 20:
                continue
            fh.write(f">{tag}\n{prot}\n")
            genes[tag] = {"start": int(f.location.start), "end": int(f.location.end),
                          "strand": "+" if f.location.strand == 1 else "-"}
    print(f"proteínas traducidas: {len(genes)}")
    return genes


def k12_orthologs() -> dict[str, tuple[str, float]]:
    """Mejor hit de cada proteína del destino contra K-12."""
    out = ROOT / "data" / "dbtbs" / "bsu_vs_k12.tsv"
    if not out.exists():
        subprocess.run(["diamond", "blastp", "-q", str(FAA), "-d", str(K12DB),
                        "-o", str(out), "-f", "6", "qseqid", "sseqid",
                        "pident", "evalue", "bitscore",
                        "--max-target-seqs", "1", "--fast", "--threads", "8",
                        "--quiet"], check=True)
    best = {}
    for line in open(out):
        f = line.split("\t")
        if len(f) >= 5:
            q, s, bits = f[0], f[1], float(f[4])
            if q not in best or bits > best[q][1]:
                best[q] = (s, bits)
    print(f"genes del destino con ortólogo en K-12: {len(best)}")
    return best


def k12_tss_offsets() -> dict[str, int]:
    """Para cada gen de K-12 con TSS aguas arriba: offset = inicio − TSS."""
    _, genes = F.read_gbff(ROOT / "data" / "test_ecoli.gbff")
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    tss = {"+": [], "-": []}
    for r in csv.DictReader(open(ROOT / "data" / "regulondb" / "tss_k12_union.tsv"),
                            delimiter="\t"):
        if (r.get("pos") or "").isdigit():
            tss["-" if (r.get("strand") or "+").startswith("-") else "+"].append(int(r["pos"]))
    for v in tss.values():
        v.sort()
    off = {}
    for g in genes:
        st = "+" if str(g["strand"]) in ("1", "+") else "-"
        arr = tss[st]
        start = int(g["start"]) if st == "+" else int(g["end"])
        # TSS mas cercano aguas arriba, hasta 600 pb
        lo = bisect.bisect_left(arr, start - 600)
        hi = bisect.bisect_right(arr, start)
        if hi > lo:
            p = arr[hi - 1]
            off[g["locus_tag"]] = start - p
    print(f"genes de K-12 con TSS aguas arriba: {len(off)}")
    return off


def segment(genes, tss_of) -> list[list[str]]:
    ops, cur = [], []
    order = sorted(genes.items(), key=lambda kv: kv[1]["start"])
    cur = [order[0][0]]
    for (ta, ga), (tb, gb) in zip(order, order[1:]):
        cut = ga["strand"] != gb["strand"]
        if not cut:
            # ¿hay TSS transferido en el hueco?
            arr = tss_of[ga["strand"]]
            i = bisect.bisect_left(arr, ga["end"] + 1)
            cut = i < len(arr) and arr[i] <= gb["start"] - 1
        if cut:
            ops.append(cur); cur = [tb]
        else:
            cur.append(tb)
    ops.append(cur)
    return ops


def main() -> int:
    genes = translate_genome()
    ortho = k12_orthologs()
    off = k12_tss_offsets()

    # transferir
    tss_of = {"+": [], "-": []}
    n = 0
    for tag, g in genes.items():
        hit = ortho.get(tag)
        if not hit or hit[0] not in off:
            continue
        d = off[hit[0]]
        if g["strand"] == "+":
            p = g["start"] - d
        else:
            p = g["end"] + d
        if p > 0:
            tss_of[g["strand"]].append(p)
            n += 1
    for v in tss_of.values():
        v.sort()
    print(f"TSS transferidos: {n}")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "strand", "source"])
        for st in "+-":
            for p in tss_of[st]:
                w.writerow([p, st, "transferido_K12"])

    # evaluar contra operones de BSGatlas
    gt = {}
    for r in csv.DictReader(open(ROOT / "data/dbtbs/gt_operons.tsv"), delimiter="\t"):
        for x in (r.get("genes") or "").split(","):
            if x.strip():
                gt[x.strip()] = r["operon_id"]
    ops = segment(genes, tss_of)
    pred, real = set(), set()
    for op in ops:
        tags = [t for t in op if t in gt]
        for a, b in combinations(sorted(tags), 2):
            pred.add(frozenset((a, b)))
    by = defaultdict(set)
    for g, o in gt.items():
        by[o].add(g)
    for tags in by.values():
        for a, b in combinations(sorted(tags), 2):
            real.add(frozenset((a, b)))
    tp = len(pred & real)
    P = tp / len(pred) if pred else 0
    R = tp / len(real) if real else 0
    F1 = 2 * P * R / (P + R) if P + R else 0
    print(f"\ncontra operones BSGatlas: F1={F1:.3f} P={P:.3f} R={R:.3f} "
          f"({len(ops)} operones)")
    print("referencia: TSS propios 0,799 | motor de distancia 0,752")
    return 0


if __name__ == "__main__":
    sys.exit(main())
