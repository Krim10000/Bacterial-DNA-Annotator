#!/usr/bin/env python3
"""operon_components.py — ¿qué piezas identificadas tiene cada operón?

Antes de integrar más criterios hace falta saber qué está realmente
identificado en cada operón, no solo qué genes lo forman. Este script cruza
`OperonSet` × `TuPromOperonTFSet` × `PromoterSet` × `TerminatorSet` ×
`EvidenceSet` y deja una fila por operón con sus componentes.

Columnas de `gt_operon_components.tsv`:

    operon_id, operon_name, n_genes, n_tus, confidence, evidence_codes,
    has_promoter, has_tss, has_minus10, has_minus35, has_sigma, has_tf,
    has_terminator, n_pmids, n_components, complete

`n_components` cuenta de 0 a 6: promotor, TSS, -10, -35, sigma, terminador.
`complete` = tiene las 6.

Uso
---
    python scripts/validation/operon_components.py
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "regulondb"
OUT = RAW / "gt" / "gt_operon_components.tsv"
DETAIL = RAW / "gt" / "gt_operon_tu_detail.tsv"


def load(name: str) -> list[dict]:
    path = RAW / f"{name}.tsv"
    rows = list(csv.reader(open(path), delimiter="\t"))
    header = [c.split(")")[-1].strip() for c in rows[0]]
    return [dict(zip(header, r)) for r in rows[1:] if len(r) >= 3]


def main() -> int:
    ops = load("OperonSet")
    tus = load("TuPromOperonTFSet")
    prom = {p["id"]: p for p in load("PromoterSet")}
    evidence = {e["code"]: e for e in load("EvidenceSet")}

    by_operon: dict[str, list[dict]] = defaultdict(list)
    for t in tus:
        by_operon[t["operonId"]].append(t)

    detail = []
    rows = []
    for op in ops:
        oid = op["operonId"]
        group = by_operon.get(oid, [])
        flags = Counter()
        pmids: set[str] = set()
        for t in group:
            p = prom.get(t.get("promoterId", ""), {})
            flags["has_promoter"] += bool(t.get("promoterId", "").strip())
            flags["has_tss"] += bool(p.get("posTSS", "").strip())
            flags["has_minus10"] += bool(p.get("boxMinus10seq", "").strip())
            flags["has_minus35"] += bool(p.get("boxMinus35seq", "").strip())
            flags["has_sigma"] += bool(t.get("sigmaFactor", "").strip())
            flags["has_tf"] += bool(t.get("tfsNames", "").strip())
            flags["has_terminator"] += bool(
                (t.get("terminatorIds", "") or "").strip()
                or (t.get("terminatorPositions", "") or "").strip())
            detail.append({
                "operon_id": oid, "operon_name": op.get("operonName", ""),
                "tu_id": t.get("tuId", ""), "tu_name": t.get("tuName", ""),
                "genes": t.get("tuGenesBnumber", "").replace("|", ","),
                "promoter_id": t.get("promoterId", ""),
                "promoter_name": t.get("promoterName", ""),
                "tss": p.get("posTSS", ""),
                "minus10_pos": p.get("boxMinus10pos", ""),
                "minus10_seq": p.get("boxMinus10seq", ""),
                "minus35_pos": p.get("boxMinus35pos", ""),
                "minus35_seq": p.get("boxMinus35seq", ""),
                "sigma": t.get("sigmaFactor", "") or p.get("sigmaFactor", ""),
                "tfs": t.get("tfsNames", ""),
                "terminators": (t.get("terminatorIds", "") or "").replace("|", ","),
                "terminator_positions": t.get("terminatorPositions", ""),
                "pmids": p.get("pmids", ""),
                "evidence": t.get("tuEvidence", ""),
            })
            for src in (p.get("pmids", ""), t.get("tuEvidence", "")):
                pmids.update(x for x in src.replace(";", ",").split(",") if x.strip().isdigit())
        codes = [c for c in op.get("operonEvidence", "").replace(",", ":").split(":")
                 if c.strip()]
        comps = sum(1 for k in ("has_promoter", "has_tss", "has_minus10",
                                "has_minus35", "has_sigma", "has_terminator")
                    if flags[k])
        rows.append({
            "operon_id": oid, "operon_name": op.get("operonName", ""),
            "n_genes": op.get("numberOfGenes", ""), "n_tus": len(group),
            "confidence": op.get("confidenceLevel", ""),
            "evidence_codes": ";".join(sorted(set(codes))),
            # booleano por operón: ¿alguna de sus TUs tiene la pieza?
            **{k: int(bool(flags[k])) for k in ("has_promoter", "has_tss",
                                                "has_minus10", "has_minus35",
                                                "has_sigma", "has_tf",
                                                "has_terminator")},
            "n_pmids": len(pmids), "n_components": comps,
            "complete": "yes" if comps == 6 else "no",
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if detail:
        with open(DETAIL, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(detail[0]), delimiter="\t")
            w.writeheader()
            w.writerows(detail)

    n = len(rows)
    multi = [r for r in rows if str(r["n_genes"]).isdigit() and int(r["n_genes"]) >= 2]
    print(f"operones: {n} (>=2 genes: {len(multi)}) -> {OUT.relative_to(ROOT)}")
    print(f"detalle por TU: {len(detail)} -> {DETAIL.relative_to(ROOT)}")
    print("\ncobertura de componentes por operón:")
    for k in ("has_promoter", "has_tss", "has_minus10", "has_minus35",
              "has_sigma", "has_tf", "has_terminator"):
        c = sum(r[k] for r in rows)
        print(f"  {k:16s} {c:5d}/{n} ({100 * c / n:4.0f}%)")
    print("\ncuántas piezas tiene cada operón (0-6):")
    for k, v in sorted(Counter(r["n_components"] for r in rows).items()):
        print(f"  {k} piezas: {v:5d} ({100 * v / n:4.0f}%)")
    print(f"\ncompletos (6/6): {sum(1 for r in rows if r['complete'] == 'yes')}")
    print(f"con >=1 PMID: {sum(1 for r in rows if r['n_pmids'])} "
          f"({100 * sum(1 for r in rows if r['n_pmids']) / n:.0f}%)")
    print(f"con códigos de evidencia: {sum(1 for r in rows if r['evidence_codes'])}")
    if evidence:
        print(f"códigos en el catálogo EvidenceSet: {len(evidence)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
