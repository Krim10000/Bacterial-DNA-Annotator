#!/usr/bin/env python3
"""build_candidates.py — de un genoma anotado a hipótesis de complejos.

Dos pasos
---------
1. `--fit`: aprende a puntuar pares con el gold de K-12. Los negativos son todos
   los pares entre genes que sí están en algún complejo curado (mundo cerrado),
   así que la prevalencia es realista para ese universo.
2. `--apply`: puntúa los pares candidatos del genoma, los filtra y agrupa en
   candidatos a complejo con top-k mutuo.

Salidas
-------
- `pairs_<genome>.tsv`: cada par con su score y la evidencia declarada.
- `clusters_<genome>.tsv`: candidatos a complejo, con confianza y evidencia.

Uso
---
    # K-12: ajustar y aplicar
    python scripts/complexes/build_candidates.py --fit --apply

    # solo aplicar un genoma nuevo con un modelo ya ajustado
    python scripts/complexes/build_candidates.py --apply \
        --genome GCF_000006985.1 --genes data/regulation/<GCF>/genes.tsv \
        --operons data/operon_finder/<GCF>_operons.tsv \
        --string data/complexes/raw/string_<taxid>_full.txt.gz
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import complexes_lib as L  # noqa: E402


def load_context(a: argparse.Namespace) -> L.Context:
    genes = L.read_genes_tsv(a.genes)
    if not genes:
        raise SystemExit(f"sin genes en {a.genes}")
    op_of, _ = L.read_operons(a.operons)
    synteny = L.read_synteny(a.synteny)
    string = (L.read_string_full(a.string, a.string_species)
              if a.string and Path(a.string).is_file() else {})
    print(f"  genes={len(genes)} operones={len(set(op_of.values()))} "
          f"STRING={len(string)} sintenia={len(synteny)}")
    return L.build_context(genes, op_of, synteny, string)


def load_gold(path: Path) -> list[list[str]]:
    out: list[list[str]] = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tags = [g for g in (r.get("locus_tags") or "").split(",") if g]
            if len(tags) >= 2:
                out.append(sorted(set(tags)))
    return out


def fit(a: argparse.Namespace, ctx: L.Context) -> dict:
    try:
        import sklearn  # noqa: F401
    except ImportError:  # pragma: no cover
        raise SystemExit("--fit requiere scikit-learn: pip install '.[full]'")
    gold = load_gold(a.gold)
    if not gold:
        raise SystemExit(f"gold vacío: {a.gold}")
    universe = sorted({g for c in gold for g in c if g in ctx.genes})
    pos = {L.pair_key(g, h) for c in gold for g, h in itertools.combinations(sorted(c), 2)}
    pairs = list(itertools.combinations(universe, 2))
    print(f"  ajuste: {len(universe)} genes, {len(pairs)} pares, "
          f"{len(pos)} positivos")
    names = [f for f in L.FEATURES_ALL
             if a.include_string_db or f != "string_db"]
    names_gen = L.FEATURES_GENOMIC
    X: list[list[float]] = []
    Xg: list[list[float]] = []
    for g, h in pairs:
        fr = L.features_row(g, h, ctx)
        X.append([fr[n] for n in names])
        Xg.append([fr[n] for n in names_gen])
    y = [1 if p in pos else 0 for p in pairs]
    model = L.fit_model(X, y, names, C=a.C)
    model.update({
        "trained_genome": a.genome,
        "gold": str(a.gold),
        "gold_complexes": len(gold),
        "universe_genes": len(universe),
        "include_string_db": bool(a.include_string_db),
    })
    L.save_model(model, a.model)
    print(f"  modelo: {a.model}  (db={'sí' if a.include_string_db else 'no'})")

    # Modelo de respaldo para genomas sin STRING: solo señales genómicas.
    mgen = L.fit_model(Xg, y, names_gen, C=a.C)
    mgen.update({
        "trained_genome": a.genome,
        "gold": str(a.gold),
        "gold_complexes": len(gold),
        "universe_genes": len(universe),
        "genomic_only": True,
    })
    gen_path = genomic_path(a.model)
    L.save_model(mgen, gen_path)
    print(f"  respaldo sin STRING: {gen_path}")
    return model


def genomic_path(model_path: Path) -> Path:
    return model_path.with_name(model_path.stem + ".genomic" + model_path.suffix)


def apply(a: argparse.Namespace, ctx: L.Context, model: dict) -> tuple[list[dict], list[dict]]:
    gen_path = genomic_path(a.model)
    if (not ctx.string and gen_path.is_file()
            and any(f.startswith("string_") for f in model["features"])):
        print(f"  sin STRING: uso el modelo genómico {gen_path.name}")
        model = L.load_model(gen_path)
    cand = L.candidate_pairs(ctx, min_string=a.min_string)
    print(f"  candidatos: {len(cand)} pares")
    names = model["features"]
    rows: list[dict] = []
    edges: dict[tuple[str, str], float] = {}
    for g, h in cand:
        try:
            f = L.features_row(g, h, ctx)
            p = L.predict(model, [f[n] for n in names])
        except (KeyError, ValueError):
            continue
        edges[L.pair_key(g, h)] = p
        if p < a.min_pair:
            continue
        geo = L.pair_geometry(g, h, ctx)
        rows.append({
            "gene_a": geo["upstream"], "gene_b": geo["downstream"],
            "name_a": ctx.genes[geo["upstream"]].gene_name,
            "name_b": ctx.genes[geo["downstream"]].gene_name,
            "score": round(p, 4),
            "gap_bp": geo["gap"], "dist_bp": geo["dist"],
            "same_operon": int(f["same_operon"]), "adjacent": int(f["adjacent"]),
            "synteny": round(f["synteny"], 3), "same_ec": int(f["same_ec"]),
            "string_exp": f["string_exp"], "string_nbr": f["string_nbr"],
            "string_fus": f["string_fus"], "string_cooc": f["string_cooc"],
            "string_coex": f["string_coex"],
            "string_db": f["string_db"],
            "evidence": ";".join(L.evidence_tags(g, h, ctx)),
        })
    rows.sort(key=lambda r: (-r["score"], r["gene_a"], r["gene_b"]))

    clusters: list[dict] = []
    for i, cl in enumerate(L.cluster_edges(edges, top_k=a.top_k, min_p=a.min_p), 1):
        s = L.cluster_summary(cl, edges, ctx)
        clusters.append({
            "cluster_id": f"{a.genome}:cx{i}",
            "genome": a.genome,
            "n_genes": s["n_genes"],
            "genes": ",".join(cl),
            "gene_names": ",".join(ctx.genes[g].gene_name or g for g in cl),
            "mean_score": s["mean_p"], "min_score": s["min_p"],
            "confidence": s["confidence"], "evidence": s["evidence"],
        })
    clusters.sort(key=lambda r: ({"alta": 0, "media": 1, "baja": 2}[r["confidence"]],
                                 -r["mean_score"]))
    return rows, clusters


def write_tsv(rows: list[dict], cols: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  escrito: {path} ({len(rows)} filas)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", default="ecoli_k12")
    ap.add_argument("--genes", type=Path, default=L.CX / "k12_genes.tsv")
    ap.add_argument("--operons", type=Path,
                    default=L.DATA / "operon_finder" / "ecoli_k12_operons_rules.tsv")
    ap.add_argument("--string", type=Path,
                    default=L.CX / "raw" / "string_511145_full.txt.gz")
    ap.add_argument("--string-species", default="511145")
    ap.add_argument("--synteny", type=Path,
                    default=L.DATA / "operon_finder" / "k12_synteny.tsv")
    ap.add_argument("--gold", type=Path, default=L.CX / "gold_k12.tsv")
    ap.add_argument("--model", type=Path, default=L.CX / "model_k12.json")
    ap.add_argument("--out-dir", type=Path, default=L.CX)
    ap.add_argument("--fit", action="store_true", help="(re)ajustar el modelo")
    ap.add_argument("--apply", action="store_true", help="puntuar y agrupar")
    ap.add_argument("--fit-only", action="store_true")
    ap.add_argument("--apply-only", action="store_true")
    ap.add_argument("--include-string-db", action="store_true",
                    help="usa el canal database de STRING (circular con el gold)")
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--min-string", type=float, default=0.2,
                    help="señal STRING mínima para considerar un par")
    ap.add_argument("--min-pair", type=float, default=0.1,
                    help="probabilidad mínima para escribir un par")
    ap.add_argument("--min-p", type=float, default=0.1,
                    help="probabilidad mínima para una arista al agrupar")
    ap.add_argument("--top-k", type=int, default=2,
                    help="vecinos mutuos para el clustering")
    a = ap.parse_args(argv)

    do_fit = a.fit or a.fit_only or not (a.apply or a.apply_only)
    do_apply = a.apply or a.apply_only or not (a.fit or a.fit_only)

    print("contexto")
    ctx = load_context(a)

    model = None
    if do_fit:
        print("ajuste")
        model = fit(a, ctx)
    if do_apply:
        print("aplicación")
        if model is None:
            if not a.model.exists():
                raise SystemExit(f"falta {a.model}; corre --fit primero")
            model = L.load_model(a.model)
        rows, clusters = apply(a, ctx, model)
        write_tsv(rows, list(rows[0]) if rows else [], a.out_dir / f"pairs_{a.genome}.tsv")
        write_tsv(clusters, list(clusters[0]) if clusters else [],
                  a.out_dir / f"clusters_{a.genome}.tsv")
        n_high = sum(1 for c in clusters if c["confidence"] == "alta")
        print(f"  clusters: {len(clusters)} (alta: {n_high})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
