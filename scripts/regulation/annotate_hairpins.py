#!/usr/bin/env python3
"""annotate_hairpins.py — anota horquillas de terminación en cualquier genoma.

Una horquilla es una **pista**, no un decisor de frontera: 70 de los 375
terminadores curados de RegulonDB están *dentro* del operón, y 8 de los 37
dependientes de Rho también tienen horquilla. Pero la estructura sí indica que
algo puede ocurrir en la maquinaria (pausa de la polimerasa), así que se anota
con su clase y su evidencia.

Método: en cada región intergénica se deslizan ventanas de tamaño fijo y se
toma la mejor. Se miden tallo, GC, bucle y poli-U (proxy) y ΔG (ViennaRNA).

Salida (`data/operon_finder/<genoma>_hairpins.tsv`):
    genome contig start end strand class stem gc loop polyU dG gap
    left_gene right_gene

Uso
---
    python scripts/regulation/annotate_hairpins.py \
        --gff data/annotations_ab_initio/<G>/annotation.gff3 \
        --fasta data/raw/<G>_genomic.fna.gz --genome <G>

    python scripts/regulation/annotate_hairpins.py --k12
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402
import terminator_classes as T  # noqa: E402

OUT = ROOT / "data" / "operon_finder"
WINDOW = 44
STEP = 16


def best_hairpin(seq: str, use_dg: bool) -> dict | None:
    """Estructura sobre la región completa (una pasada), no ventanas deslizantes.

    Deslizar ventanas de 44 nt con paso 16 era 16x más caro y no cambiaba la
    clasificación: `parse_structure` ya devuelve la posición del mejor tallo.
    """
    if len(seq) < 20:
        return None
    region = seq[:400]                      # las intergénicas enormes se truncan
    st = T.parse_structure(region)
    if st["stem"] < 4:
        return None
    best = dict(st)
    best["cls"] = T.classify(best)
    if use_dg:
        try:
            import RNA
            lo = max(0, st["start"] - 10)
            win = region[lo:st["start"] + 2 * st["stem"] + st["loop"] + 20]
            best["dG"] = round(-RNA.fold(win)[1], 2) if len(win) >= 20 else ""
        except Exception:
            best["dG"] = ""
    else:
        best["dG"] = ""
    best["window"] = region[st["start"]:st["start"] + 2 * 12 + 10]
    return best


def _runs(seq: str) -> list[int]:
    out, run = [], 0
    for b in seq.upper() + "A":
        run = run + 1 if b == "T" else 0
        if run:
            out.append(run)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gff", type=Path)
    ap.add_argument("--fasta", type=Path)
    ap.add_argument("--gbff", type=Path)
    ap.add_argument("--genome", default="ecoli_k12")
    ap.add_argument("--k12", action="store_true", help="usa data/test_ecoli.gbff")
    ap.add_argument("--no-dg", action="store_true", help="omite ViennaRNA")
    ap.add_argument("--out", type=Path, default=None,
                    help="carpeta de salida (default data/operon_finder); permite "
                         "separar genomas excluidos en data/operon_finder/excluded")
    a = ap.parse_args()

    if a.k12:
        seq, genes = F.read_gbff(ROOT / "data" / "test_ecoli.gbff")
        genes = F.ordered([g for g in genes if g.get("locus_tag")])
        seqs = {"": seq}
        a.genome = "ecoli_k12"
    else:
        if not (a.gff and a.fasta):
            print("hacen falta --gff y --fasta (o --k12)", file=sys.stderr)
            return 1
        seqs = F.read_fasta(a.fasta)
        genes = F.ordered(F.read_gff_genes(a.gff))

    rows = []
    by_contig: dict[str, list[dict]] = {}
    for g in genes:
        by_contig.setdefault(g.get("seqid", ""), []).append(g)

    for contig, glist in by_contig.items():
        s = seqs.get(contig, "")
        if not s:
            continue
        glist = F.ordered(glist)
        for x, y in zip(glist, glist[1:]):
            gap = y["start"] - x["end"] - 1
            if gap < 20:
                continue
            inter = s[x["end"]:y["start"]]
            h = best_hairpin(inter, not a.no_dg)
            if h is None:
                continue
            rows.append({
                "genome": a.genome, "contig": contig,
                "start": x["end"] + h["start"],
                "end": x["end"] + h["start"] + len(h["window"]),
                "strand": "intergénica", "class": h["cls"],
                "stem": h["stem"], "gc": f"{h['gc']:.2f}", "loop": h["loop"],
                "polyU": h["polyU"], "dG": h["dG"], "gap": gap,
                "left_gene": x.get("locus_tag", ""), "right_gene": y.get("locus_tag", ""),
            })

    out_dir = a.out if a.out else OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{a.genome}_hairpins.tsv"
    # escritura atómica: si el proceso muere a mitad, no queda un .tsv truncado
    # que el driver confunda con "ya hecho"
    tmp = path.with_suffix(".tsv.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)

    c = Counter(r["class"] for r in rows)
    dg = [r["dG"] for r in rows if r["dG"] != "" and r["class"] == "A_fuerte"]
    try:
        destino = path.relative_to(ROOT)
    except ValueError:
        destino = path
    print(f"{a.genome}: {len(rows)} regiones intergénicas con horquilla -> {destino}")
    for k, v in c.most_common():
        print(f"  {k:10s} {v:5d} ({100 * v / len(rows):4.0f}%)")
    if dg:
        print(f"  ΔG mediana clase A: {statistics.median(dg):.1f} kcal/mol")
    return 0


if __name__ == "__main__":
    sys.exit(main())
