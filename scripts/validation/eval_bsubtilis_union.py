#!/usr/bin/env python3
"""eval_bsubtilis_union.py — ¿la unión de TSS también mejora en B. subtilis?

En K-12, unir TSS curados y de alto rendimiento subió el F1 de 0,785 a 0,798.
Aquí se repite en *B. subtilis* con sus dos fuentes:

- **BSGatlas** (`gt_tss.tsv`): integra DBTBS, BsubCyc, SubtiWiki y datos propios.
- **DBTBS** (`dbtbs_tss.tsv`): raspado de sus páginas de promotor (TSS estimado
  desde la posición absoluta de la región −42:+4).

Se evalúa contra **dos ground truths** para no depender de uno solo:
- operones de BSGatlas (`gt_operons.tsv`, 2.266)
- operones de DBTBS (`gt_operons_dbtbs.tsv`, 1.093), mapeados a locus tags BSU

Uso
---
    python scripts/validation/eval_bsubtilis_union.py
"""
from __future__ import annotations

import csv
import re
import sys
from itertools import combinations
from collections import defaultdict
from pathlib import Path
import bisect

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import find_operons as F  # noqa: E402

GBFF = ROOT / "data" / "dbtbs" / "subtilis.gbff"
GFF = ROOT / "data" / "dbtbs" / "BSGatlas_v1.1.gff"


def genes_bsu():
    _, genes = F.read_gbff(GBFF)
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    return genes


def load_tss(path: Path) -> dict[str, list[int]]:
    out = {"+": [], "-": []}
    if not path.exists():
        return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        p = (r.get("pos") or r.get("pos_1") or "").strip()
        if not p.isdigit():
            continue
        st = (r.get("strand") or "+").strip()
        out["-" if st.startswith("-") or st == "reverse" else "+"].append(int(p))
    for v in out.values():
        v.sort()
    return out


def union(a: dict, b: dict, dist: int = 10) -> dict[str, list[int]]:
    out = {}
    for st in "+-":
        merged = sorted(set(a.get(st, [])) | set(b.get(st, [])))
        kept = []
        for p in merged:
            if not kept or p - kept[-1] > dist:
                kept.append(p)
        out[st] = kept
    return out


def gt_bsgatlas() -> dict[str, str]:
    out = {}
    for r in csv.DictReader(open(ROOT / "data/dbtbs/gt_operons.tsv"), delimiter="\t"):
        for g in (r.get("genes") or "").split(","):
            if g.strip():
                out[g.strip()] = r["operon_id"]
    return out


def gt_dbtbs() -> dict[str, str]:
    """DBTBS da nombres de gen; se mapean a locus tags por Name del GFF."""
    byname: dict[str, str] = {}
    for line in open(GFF, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"locus_tag=([^;]+)", f[8])
        n = re.search(r"Name=([^;]+)", f[8])
        if m and n:
            tag = m.group(1).split(",")[0].strip()
            byname.setdefault(n.group(1).strip(), tag)
    out = {}
    for r in csv.DictReader(open(ROOT / "data/dbtbs/gt_operons_dbtbs.tsv"), delimiter="\t"):
        for nm in (r.get("genes") or "").split(","):
            tag = byname.get(nm.strip())
            if tag:
                out[tag] = r["operon"]
    return out


def segment(genes, tss) -> list[list[dict]]:
    ops, cur = [], [genes[0]]
    strand_of = {g["locus_tag"]: ("+" if str(g["strand"]) in ("1", "+") else "-")
                 for g in genes}
    for a, b in zip(genes, genes[1:]):
        if a["strand"] != b["strand"] or a.get("seqid") != b.get("seqid"):
            ops.append(cur); cur = [b]; continue
        arr = tss[strand_of[a["locus_tag"]]]
        i = bisect.bisect_left(arr, int(a["end"]) + 1)
        cut = i < len(arr) and arr[i] <= int(b["start"]) - 1
        if cut:
            ops.append(cur); cur = [b]
        else:
            cur.append(b)
    ops.append(cur)
    return ops


def score(ops, gt) -> tuple[float, float, float, int]:
    pred = set()
    for op in ops:
        tags = [g["locus_tag"] for g in op if g["locus_tag"] in gt]
        for a, b in combinations(sorted(tags), 2):
            pred.add(frozenset((a, b)))
    by = defaultdict(set)
    for g, o in gt.items():
        by[o].add(g)
    real = set()
    for tags in by.values():
        for a, b in combinations(sorted(tags), 2):
            real.add(frozenset((a, b)))
    tp = len(pred & real)
    P = tp / len(pred) if pred else 0
    R = tp / len(real) if real else 0
    F = 2 * P * R / (P + R) if P + R else 0
    return F, P, R, len(ops)


def main() -> int:
    genes = genes_bsu()
    bsg = load_tss(ROOT / "data/dbtbs/gt_tss.tsv")
    dbt = load_tss(ROOT / "data/dbtbs/dbtbs_tss.tsv")
    un = union(bsg, dbt)
    print(f"genes: {len(genes)}")
    print(f"TSS BSGatlas: {len(bsg['+'])+len(bsg['-'])} | "
          f"DBTBS: {len(dbt['+'])+len(dbt['-'])} | "
          f"unión: {len(un['+'])+len(un['-'])}")
    for nombre, gt in (("BSGatlas", gt_bsgatlas()), ("DBTBS", gt_dbtbs())):
        print(f"\nGT {nombre} ({len(set(gt.values()))} operones, {len(gt)} genes)")
        for etiqueta, tss in (("BSGatlas", bsg), ("DBTBS", dbt), ("unión", un)):
            f, p, r, n = score(segment(genes, tss), gt)
            print(f"  TSS {etiqueta:8s} F1={f:.3f} P={p:.3f} R={r:.3f} ops={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
