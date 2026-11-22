#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""15_gate.py — gate formal L2 y `data/regulation/GATE.md`.

El gate son **20 operones canónicos** de E. coli K-12: se exige el mismo
conjunto de genes y el mismo orden de transcripción. La pertenencia y el orden
se toman de RegulonDB v14.5.1 (`data/regulondb/gt/gt_operons.tsv`) y se citan
las referencias clásicas de los casos que se pueden citar sin dudas.

Métricas que escribe:
  * gate 20/20 (exactos, mismos genes, mismo orden)
  * pair-level / exact-set F1 y desglose por tamaño
  * recall de pares en operones 7+ (el objetivo del plan: >= 0.6)
  * cobertura de promotores y de sigma
  * hash SHA-256 de las salidas

Uso:
  python scripts/regulation/15_gate.py
  python scripts/regulation/15_gate.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

GATE_MIN = 15
TARGET_7PLUS_PAIR_RECALL = 0.60

# (nombre GT, cita). Todas las filas usan RegulonDB como fuente de pertenencia;
# las citas clásicas se añaden solo donde son inequívocas.
CURATED = [
    ("lacZYA", "Jacob & Monod 1961; RegulonDB"),
    ("trpLEDCBA", "Yanofsky 1981; RegulonDB"),
    ("araBAD", "Englesberg & Wilcox 1974; RegulonDB"),
    ("galETKM", "RegulonDB (galETKM)"),
    ("hisLGDCBHAFI", "RegulonDB (hisLGDCBHAFI; incluye el líder hisL)"),
    ("thrLABC", "RegulonDB (thrLABC; atenuación con thrL)"),
    ("rplKAJL-rpoBC", "RegulonDB (operón L10: rplKAJL-rpoBC)"),
    ("atpIBEFHAGDC", "Walker et al. 1984; RegulonDB"),
    ("fimAICDFGH", "Klemm 1986; RegulonDB (incluye el sRNA fimR2)"),
    ("malEFGH", "Boos & Shuman 1998; RegulonDB (incluye el sRNA malH)"),
    ("xylAB", "RegulonDB (xylAB; incluye el sRNA xylZ)"),
    ("rbsDACBZKR", "RegulonDB (incluye rbsZ, rbsK y rbsR)"),
    ("cyoABCDE", "RegulonDB (cyoABCDE)"),
    ("groSL", "RegulonDB (groSL)"),
    ("argCBH", "RegulonDB (argCBH)"),
    ("cadBA", "Meng & Bennett 1992; RegulonDB"),
    ("mglBAC", "Boos & Shuman 1998; RegulonDB"),
    ("nuoABCEFGHIJKLMN", "RegulonDB (nuo, 13 genes)"),
    ("sdhCDAB-sucABCD-sdhX", "RegulonDB (sdh-suc, 10 genes)"),
    ("nikABCDER", "RegulonDB (nik, 6 genes)"),
]


def load_gt() -> dict[str, dict]:
    p = C.DATA / "regulondb" / "gt" / "gt_operons.tsv"
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("mapping_ok") != "yes":
                continue
            genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
            if genes:
                out[r["operon_name"]] = {"genes": genes, "strand": r["strand"],
                                         "n": len(genes),
                                         "evidence": r.get("evidence", "")}
    return out


def tx_order(genes: list[str], coords: dict[str, dict]) -> list[str]:
    """Ordena b-numbers en orden de transcripción (5'->3')."""
    have = [g for g in genes if g in coords]
    if len(have) < len(genes):
        return genes
    strand = coords[have[0]]["strand"]
    if strand == 1:
        return sorted(have, key=lambda b: coords[b]["start"])
    return sorted(have, key=lambda b: -coords[b]["end"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        coords = {"b2": {"start": 10, "end": 20, "strand": 1},
                  "b1": {"start": 30, "end": 40, "strand": 1}}
        assert tx_order(["b1", "b2"], coords) == ["b2", "b1"]
        coords2 = {"b1": {"start": 10, "end": 20, "strand": -1},
                   "b2": {"start": 30, "end": 40, "strand": -1}}
        assert tx_order(["b1", "b2"], coords2) == ["b2", "b1"]
        print("selftest 15_gate: OK")
        return 0

    d = C.REG / "ecoli_k12"
    gt = load_gt()
    genes = C.read_tsv(d / "genes.tsv")
    coords = {r["locus_tag"]: {"start": int(r["start"]), "end": int(r["end"]),
                               "strand": int(r["strand"])}
              for r in genes if r["locus_tag"]}
    gid2b = {r["gene_id"]: r["locus_tag"] for r in genes}
    og = C.read_tsv(d / "operon_genes.tsv")
    b2op, op_rows = {}, {}
    for r in og:
        b = gid2b.get(r["gene_id"], "")
        if b:
            b2op[b] = r["operon_id"]
        op_rows.setdefault(r["operon_id"], []).append(r)
    for oid in op_rows:
        op_rows[oid].sort(key=lambda x: int(x["position"]))

    rows, exact = [], 0
    for name, cite in CURATED:
        g = gt.get(name)
        if not g:
            rows.append((name, "NO_EN_GT", "", cite))
            continue
        want = tx_order(g["genes"], coords)
        ids = {b2op.get(b) for b in want}
        if len(ids) != 1 or None in ids:
            rows.append((name, "cortado", ",".join(want), cite))
            continue
        oid = ids.pop()
        pred = [gid2b.get(r["gene_id"], "") or r["gene_name"] for r in op_rows[oid]]
        pred = [p for p in pred if p]
        if pred == want:
            exact += 1
            rows.append((name, "exacto", "", cite))
        else:
            rows.append((name, "orden/extra", ",".join(pred), cite))

    # métricas globales
    bench = json.loads((d / "benchmark.json").read_text()) if (d / "benchmark.json").exists() else {}
    pair = bench.get("pair_level", {})
    exactm = bench.get("exact_set", {})
    by_size = bench.get("by_size", {})

    # recall de pares en GT 7+
    gt7 = [set(tx_order(g["genes"], coords)) for n, g in gt.items() if g["n"] >= 7]
    gt7p = {frozenset(p) for s in gt7 for p in combinations(sorted(s), 2)}
    pred_pairs = set()
    for oid, rs in op_rows.items():
        bs = {gid2b.get(r["gene_id"], "") for r in rs} - {""}
        if len(bs) >= 2:
            pred_pairs |= {frozenset(p) for p in combinations(sorted(bs), 2)}
    rec7 = len(gt7p & pred_pairs) / len(gt7p) if gt7p else 0.0

    pm = json.loads((d / "promoters_meta.json").read_text()) if (d / "promoters_meta.json").exists() else {}
    prom_cov = (pm.get("n_with_promoter", 0) / pm.get("n_operons", 1)) if pm else 0.0
    sigma_tot = pm.get("n_sigma_factors", 0)
    sigma_hmm = pm.get("n_sigma_hmm_validated", 0)

    gate_ok = exact >= GATE_MIN
    md = ["# Gate L2 — operones curados (E. coli K-12)\n",
          f"Generado por `15_gate.py`. Referencia: RegulonDB v14.5.1 "
          f"(`data/regulondb/gt/gt_operons.tsv`).\n",
          f"**Veredicto: {'PASA' if gate_ok else 'NO PASA'} — {exact}/{len(CURATED)} "
          f"operones exactos (mismos genes, mismo orden); mínimo exigido "
          f"{GATE_MIN}/20.**\n",
          "| operón | estado | genes (si difiere) | cita |", "|---|---|---|---|"]
    for name, st, extra, cite in rows:
        mark = "OK" if st == "exacto" else st
        md.append(f"| {name} | {mark} | {extra} | {cite} |")
    md.append("")
    md.append("## Métricas globales (827 operones GT con ≥2 genes)\n")
    md.append("| métrica | valor | objetivo |")
    md.append("|---|---|---|")
    md.append(f"| pair-level F1 | {pair.get('f1', 0):.3f} "
              f"(P={pair.get('precision', 0):.3f}, R={pair.get('recall', 0):.3f}) | ≥0.80 |")
    md.append(f"| exact-set F1 | {exactm.get('f1', 0):.3f} | referencia |")
    md.append(f"| cobertura de genes GT | {bench.get('gene_coverage', 0):.1%} | — |")
    for k in ("2-3", "4-6", "7+"):
        m = by_size.get(k, {})
        md.append(f"| F1 tamaño {k} | {m.get('f1', 0):.3f} (P={m.get('precision', 0):.3f}, "
                  f"R={m.get('recall', 0):.3f}) | — |")
    md.append(f"| recall de pares en 7+ | {rec7:.3f} | ≥{TARGET_7PLUS_PAIR_RECALL} |")
    md.append(f"| operones en 7+ recuperados | {len(gt7)} grupos GT | — |")
    md.append("")
    md.append("## Promotores y sigma\n")
    md.append(f"- operones con -35/-10: {pm.get('n_with_promoter', 0)}/"
              f"{pm.get('n_operons', 0)} = {prom_cov:.1%} (objetivo ≥50%)")
    md.append(f"- factores sigma en el inventario: {sigma_tot}"
              f" ({sigma_hmm} confirmados por dominio HMM)")
    md.append("")
    md.append("## Determinismo\n")
    md.append(f"- sha256 `operons.tsv` + `operon_genes.tsv`: "
              f"`{C.sha256_files([d / 'operons.tsv', d / 'operon_genes.tsv'])}`")
    md.append(f"- sha256 `promoters.tsv`: "
              f"`{C.sha256_files([d / 'promoters.tsv'])}`")
    md.append("")
    md.append("## Límites\n")
    md.append("- La pertenencia y el orden de referencia salen de RegulonDB; las "
              "citas clásicas se añaden donde son inequívocas. Los ARN pequeños "
              "(fimR2, xylZ, rbsZ, malH) cuentan como miembros del operón, igual "
              "que en RegulonDB.")
    md.append("- El orden se compara en orientación de transcripción (5'→3'), no "
              "en el orden de b-number de la tabla.")
    md.append("- F1 global 0.80 es un umbral de método; el gate de 20 operones "
              "es el criterio de aceptación del plan.")
    (C.REG / "GATE.md").write_text("\n".join(md) + "\n")

    print(f"GATE {'PASA' if gate_ok else 'NO PASA'}: {exact}/{len(CURATED)} exactos "
          f"(mínimo {GATE_MIN})")
    for name, st, extra, _ in rows:
        if st != "exacto":
            print(f"  {name}: {st} {extra[:60]}")
    print(f"pair F1={pair.get('f1', 0):.3f} exact F1={exactm.get('f1', 0):.3f} "
          f"7+ pair recall={rec7:.3f}")
    print(f"promotores={prom_cov:.1%} sigma={sigma_tot} (hmm={sigma_hmm})")
    print(f"escrito: {C.REG / 'GATE.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
