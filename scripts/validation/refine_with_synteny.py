#!/usr/bin/env python3
"""refine_with_synteny.py — segunda pasada: redecidir fronteras con sintenia.

La primera pasada decide cada frontera con señales locales (hueco, promotor,
terminador). Sus errores se concentran en la banda de 43–100 pb, donde el hueco
no distingue. La sintenia —¿siguen juntos los ortólogos de estos dos genes en
otros genomas?— es una señal independiente que sí puede decidir ahí.

Regla de segunda pasada
-----------------------
- Si la primera pasada **cortó** y la sintenia del par es alta (>= merge_at),
  se **deshace el corte** (los genes van juntos en la evolución).
- Si **unió** y la sintenia es muy baja (<= split_at), se **corta**.
- En el resto, se respeta la primera pasada.

Se barren los umbrales y se reporta el F1 de pares contra RegulonDB.

Uso
---
    python scripts/validation/refine_with_synteny.py \
        --boundary data/operon_finder/k12_p16b/ecoli_k12_boundary_rules.tsv \
        --synteny data/operon_finder/k12_synteny.tsv \
        --gbff data/test_ecoli.gbff
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import find_operons as F  # noqa: E402


def load_synteny(path: Path) -> dict[tuple[str, str], float]:
    out = {}
    if not Path(path).exists():
        return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            out[(r["left"], r["right"])] = float(r["frac"])
        except (KeyError, ValueError):
            continue
    return out


def build(genes: list[dict], cut_of: dict) -> list[list[dict]]:
    """Encadena genes según las decisiones de corte."""
    ops: list[list[dict]] = []
    cur: list[dict] = []
    for x, y in zip(genes, genes[1:]):
        same = x["strand"] == y["strand"] and x.get("seqid") == y.get("seqid")
        if not same:
            if cur:
                ops.append(cur)
            cur = [y]
            continue
        if not cur:
            cur = [x]
        if cut_of.get((x["locus_tag"], y["locus_tag"]), False):
            ops.append(cur)
            cur = [y]
        else:
            cur.append(y)
    if cur:
        ops.append(cur)
    return ops


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", required=True)
    ap.add_argument("--synteny", required=True)
    ap.add_argument("--gbff", default=str(ROOT / "data" / "test_ecoli.gbff"))
    ap.add_argument("--out", default=None,
                    help="escribe los operones refinados (mejor configuracion) en TSV")
    ap.add_argument("--ec", default=str(ROOT / "data" / "regulation" / "ecoli_k12" / "genes.tsv"),
                    help="tabla con locus_tag y ec para vetar cortes por coherencia")
    a = ap.parse_args()

    syn = load_synteny(Path(a.synteny))
    if not syn:
        print(f"sin datos de sintenia en {a.synteny}")
        return 1
    cut0 = {}
    for r in csv.DictReader(open(a.boundary), delimiter="\t"):
        cut0[(r["left"], r["right"])] = str(r.get("cut", "0")) == "1"

    _, genes = F.read_gbff(Path(a.gbff))
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    op_of = F.load_operon_of()

    ec = {}
    if Path(a.ec).exists():
        for row in csv.DictReader(open(a.ec), delimiter="\t"):
            if row.get("ec") and row.get("locus_tag"):
                ec[row["locus_tag"]] = row["ec"]

    def ec_same(x: str, y: str) -> bool:
        a1, a2 = ec.get(x, ""), ec.get(y, "")
        return bool(a1 and a2 and ".".join(a1.split(".")[:2]) == ".".join(a2.split(".")[:2]))

    base_cut = dict(cut0)
    p, r, f = F.pair_f1(build(genes, base_cut), op_of)
    print(f"{'configuración':34s} {'P':>6s} {'R':>6s} {'F1':>6s}  cambios")
    print(f"{'primera pasada':34s} {p:6.3f} {r:6.3f} {f:6.3f}      0")

    # EC como segunda pasada independiente
    cut = dict(cut0); n = 0
    for key, c in cut0.items():
        if c and ec_same(*key):
            cut[key] = False; n += 1
    p, r, f = F.pair_f1(build(genes, cut), op_of)
    print(f"{'2a pasada: EC (nivel 2)':34s} {p:6.3f} {r:6.3f} {f:6.3f}  {n:5d}")

    # EC + sintenia combinadas
    for merge_at in (0.05, 0.1, 0.2):
        cut = dict(cut0); n = 0
        for key, c in cut0.items():
            if not c:
                continue
            if ec_same(*key) or syn.get(key, 0) >= merge_at:
                cut[key] = False; n += 1
        ops = build(genes, cut)
        p, r, f = F.pair_f1(ops, op_of)
        print(f"{f'2a: EC + sintenia>={merge_at}':34s} {p:6.3f} {r:6.3f} {f:6.3f}  {n:5d}")
        if a.out and merge_at == 0.05:
            out = Path(a.out)
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["operon_id", "genome", "contig", "strand", "start", "end",
                            "n_genes", "genes", "gene_names", "confidence", "status"])
                for i, op in enumerate(ops, 1):
                    w.writerow([f"ref:{i}", "refined", op[0].get("seqid", ""),
                                op[0]["strand"], min(g["start"] for g in op),
                                max(g["end"] for g in op), len(op),
                                ",".join(g["locus_tag"] for g in op),
                                ",".join(g.get("gene", "") or "" for g in op),
                                "refined", "candidate"])
            print(f"-> {out}")

    # solo la dirección "deshacer corte": cortar por sintenia baja fragmenta
    # el genoma, porque la mayoría de los pares no están conservados.
    for merge_at in (0.05, 0.1, 0.2, 0.3):
        for split_at in (0.0,):
            cut = dict(cut0)
            cambios = 0
            for key, c in cut0.items():
                s = syn.get(key)
                if s is None:
                    continue
                if c and s >= merge_at:
                    cut[key] = False
                    cambios += 1
                elif split_at and (not c) and s <= split_at:
                    cut[key] = True
                    cambios += 1
            p, r, f = F.pair_f1(build(genes, cut), op_of)
            print(f"{f'merge>={merge_at} split<={split_at}':34s} "
                  f"{p:6.3f} {r:6.3f} {f:6.3f}  {cambios:5d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
