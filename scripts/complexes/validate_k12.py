#!/usr/bin/env python3
"""validate_k12.py — ¿cuánto se recupera de los complejos documentados de K-12?

Protocolo
---------
Gold: ComplexPortal 83333 (complejos curados) mapeado a b-numbers con el GenBank.
Universo: genes que aparecen en algún complejo con >=2 genes distintos. Los
negativos son todos los pares del universo que no comparten complejo (mundo
cerrado). La prevalencia resultante es la real de ese universo.

Modelo: regresión logística con validación cruzada 5-fold **agrupada por
complejo**; se reportan las predicciones fuera de fold (OOF), nunca las de
entrenamiento.

Tres conjuntos de features, para separar el mérito:
  - `genomica`: solo operón, adyacencia, distancia, sintenia y EC;
  - `sin_db`: + canales de STRING *excepto* `database` (evita el círculo con
    ComplexPortal);
  - `con_db`: + `database`, que es el techo si se permite la vía curada.

Uso
---
    python scripts/complexes/validate_k12.py
    python scripts/complexes/validate_k12.py --save-oof
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np

try:  # extra `full` de pyproject; sin él no hay validación posible
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import GroupKFold
except ImportError:  # pragma: no cover
    raise SystemExit("validate_k12.py requiere scikit-learn: pip install '.[full]'")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import complexes_lib as L  # noqa: E402

# Casos canónicos para la tabla legible (nombre exacto en ComplexPortal).
EXAMPLES = [
    ("ATP synthase complex", "ATP sintasa, operón atp"),
    ("30S small ribosomal subunit", "subunidad ribosomal 30S"),
    ("50S large ribosomal subunit", "subunidad ribosomal 50S"),
    ("DNA-directed RNA polymerase holoenzyme complex, Sigma70 variant",
     "ARN polimerasa σ70"),
    ("GroEL-GroES complex", "chaperonina GroEL/GroES"),
    ("Divisome complex", "divisoma"),
    ("Chemotaxis phosphorelay complex CheA-CheY", "quimiotaxis CheA-CheY"),
    ("AcrAB-TolC multidrug efflux transport complex", "bomba de eflujo AcrAB-TolC"),
    ("BAM complex", "complejo BAM"),
    ("2-oxoglutarate dehydrogenase complex", "2-oxoglutarato deshidrogenasa"),
    ("DNA polymerase III core complex", "ADN polimerasa III, núcleo"),
]


def load_data(a: argparse.Namespace):
    genes = L.read_genes_tsv(a.genes)
    op_of, _ = L.read_operons(a.operons)
    synteny = L.read_synteny(a.synteny)
    string = L.read_string_full(a.string, a.string_species)
    ctx = L.build_context(genes, op_of, synteny, string)
    gold: list[list[str]] = []
    with open(a.gold, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tags = [g for g in (r.get("locus_tags") or "").split(",") if g]
            if len(tags) >= 2:
                gold.append(sorted(set(tags)))
    return ctx, gold, string


def make_pairs(gold: list[list[str]]):
    universe = sorted({g for c in gold for g in c})
    pair_of_cx: dict[tuple[str, str], list[int]] = defaultdict(list)
    for i, c in enumerate(gold):
        for a, b in itertools.combinations(sorted(c), 2):
            pair_of_cx[L.pair_key(a, b)].append(i)
    pairs = list(itertools.combinations(universe, 2))
    y = np.array([1 if p in pair_of_cx else 0 for p in pairs], dtype=int)
    groups = np.array([
        f"cx{pair_of_cx[p][0]}" if p in pair_of_cx
        else f"neg{zlib.crc32((p[0] + p[1]).encode()) % 101}"
        for p in pairs
    ])
    return universe, pairs, y, groups, pair_of_cx


def cv_scores(X: list[list[float]], y: np.ndarray, groups: np.ndarray,
              names: list[str], folds: int = 5) -> np.ndarray:
    oof = np.zeros(len(y))
    gkf = GroupKFold(n_splits=folds)
    for tr, te in gkf.split(np.zeros(len(y)), y, groups):
        model = L.fit_model([X[i] for i in tr], [int(y[i]) for i in tr], names)
        for i in te:
            oof[i] = L.predict(model, X[i])
    return oof


def ranking_metrics(oof: np.ndarray, pairs: list[tuple[str, str]],
                    gold_pairs: set[tuple[str, str]],
                    universe: list[str], ks=(1, 3, 5, 10)) -> dict:
    fam: dict[str, list[tuple[float, str]]] = defaultdict(list)
    for i, (a, b) in enumerate(pairs):
        fam[a].append((oof[i], b))
        fam[b].append((oof[i], a))
    out: dict[str, float] = {}
    ranks: list[int] = []
    for g in universe:
        ordered = sorted(fam[g], reverse=True)
        for rank, (_, h) in enumerate(ordered, 1):
            if L.pair_key(g, h) in gold_pairs:
                ranks.append(rank)
                break
        else:
            ranks.append(0)
    for k in ks:
        hit = sum(1 for g in universe
                  if any(L.pair_key(g, h) in gold_pairs
                         for _, h in sorted(fam[g], reverse=True)[:k]))
        out[f"recall@{k}"] = round(hit / max(1, len(universe)), 4)
    out["mrr"] = round(sum(1.0 / r for r in ranks if r) / max(1, len(ranks)), 4)
    return out


def top1_by_evidence(oof: np.ndarray, pairs: list[tuple[str, str]],
                     gold_rows: list[dict]) -> dict:
    """Top-1 partner dentro del mismo complejo, separando evidencia A/B."""
    best: dict[str, tuple[float, str]] = {}
    for i, (a, b) in enumerate(pairs):
        if oof[i] > best.get(a, (-1.0, ""))[0]:
            best[a] = (oof[i], b)
        if oof[i] > best.get(b, (-1.0, ""))[0]:
            best[b] = (oof[i], a)
    in_cx: dict[str, set[str]] = defaultdict(set)
    gene_cls: dict[str, set[str]] = defaultdict(set)
    for r in gold_rows:
        gs = set((r.get("locus_tags") or "").split(",")) - {""}
        for g in gs:
            in_cx[g] |= gs
            gene_cls[g].add(r.get("evidence_class") or "otra")
    out = {}
    for cls in ("experimental", "inferida"):
        tot = hit = 0
        for g, classes in gene_cls.items():
            if cls not in classes:
                continue
            tot += 1
            hit += best.get(g, (-1.0, ""))[1] in in_cx[g]
        out[cls] = {"n": tot, "top1": round(hit / tot, 4) if tot else 0.0}
    return out


def pair_metrics(oof: np.ndarray, y: np.ndarray) -> dict:
    out = {
        "pr_auc": round(float(average_precision_score(y, oof)), 4),
        "roc_auc": round(float(roc_auc_score(y, oof)), 4),
        "prevalence": round(float(y.mean()), 5),
        "n_pairs": int(len(y)),
        "n_pos": int(y.sum()),
    }
    for t in (0.5, 0.9):
        sel = oof >= t
        tp = int((sel & (y == 1)).sum())
        fp = int((sel & (y == 0)).sum())
        fn = int((~sel & (y == 1)).sum())
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        out[f"precision@{t}"] = round(prec, 4)
        out[f"recall@{t}"] = round(rec, 4)
    return out


def cluster_metrics(clusters: list[list[str]], gold: list[list[str]],
                    edges: dict | None = None) -> dict:
    cl = [set(c) for c in clusters if len(c) > 1]
    gold_sets = [set(c) for c in gold]
    gold_pairs = {L.pair_key(a, b) for c in gold for a, b in itertools.combinations(sorted(c), 2)}
    pure = 0
    for c in cl:
        ps = [L.pair_key(a, b) for a, b in itertools.combinations(sorted(c), 2)]
        if all(p in gold_pairs for p in ps):
            pure += 1
    jac: list[float] = []
    by_size: dict[int, list[float]] = defaultdict(list)
    for gs in gold_sets:
        j = max((len(gs & c) / len(gs | c) for c in cl), default=0.0)
        jac.append(j)
        by_size[len(gs)].append(j)
    by_conf: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    if edges is not None:
        for c in cl:
            ps = [L.pair_key(a, b) for a, b in itertools.combinations(sorted(c), 2)]
            vals = [edges.get(p, 0.0) for p in ps]
            t = L.confidence(len(c), min(vals), sum(vals) / len(vals))
            by_conf[t][0] += 1
            by_conf[t][1] += int(all(p in gold_pairs for p in ps))
    return {
        "by_confidence": {k: {"n": v[0], "purity": round(v[1] / v[0], 4)}
                          for k, v in sorted(by_conf.items())},
        "n_clusters": len(cl),
        "n_pure": pure,
        "purity": round(pure / max(1, len(cl)), 4),
        "n_gold": len(gold_sets),
        "coverage_j50": round(sum(1 for j in jac if j >= 0.5) / max(1, len(jac)), 4),
        "exact": round(sum(1 for j in jac if j >= 0.999) / max(1, len(jac)), 4),
        "mean_jaccard": round(float(np.mean(jac)) if jac else 0.0, 4),
        "by_size": {int(k): {"n": len(v), "mean_jaccard": round(float(np.mean(v)), 4)}
                    for k, v in sorted(by_size.items())},
    }


def find_complex(gold_rows: list[dict], needle: str) -> dict | None:
    for r in gold_rows:
        if r["complex_name"].lower() == needle.lower():
            return r
    for r in gold_rows:
        if needle.lower() in r["complex_name"].lower():
            return r
    return None


def evaluate_examples(clusters: list[list[str]], gold_rows: list[dict],
                      genes: dict[str, L.Gene]) -> list[dict]:
    cl = [set(c) for c in clusters if len(c) > 1]
    out = []
    for needle, label in EXAMPLES:
        r = find_complex(gold_rows, needle)
        if r is None:
            out.append({"complex": needle, "label": label, "status": "ausente"})
            continue
        gs = set(r["locus_tags"].split(","))
        best, best_c = 0.0, set()
        for c in cl:
            j = len(gs & c) / len(gs | c)
            if j > best:
                best, best_c = j, c
        out.append({
            "complex": r["complex_name"],
            "label": label,
            "n_genes": len(gs),
            "genes": ",".join(genes[g].gene_name or g for g in sorted(gs)),
            "best_jaccard": round(best, 3),
            "recovered": best >= 0.5,
            "cluster_n": len(best_c),
            "cluster": ",".join(genes[g].gene_name or g for g in sorted(best_c)),
        })
    return out


def markdown(res: dict) -> str:
    lines = ["# Validación de complejos en *E. coli* K-12", ""]
    lines += [
        "Gold: **ComplexPortal** (taxid 83333), mapeado a b-numbers con el "
        "`db_xref` UniProt del GenBank de K-12.",
        "",
        f"- complejos usables (≥2 genes distintos): **{res['n_complexes']}**",
        f"- genes en el universo: **{res['n_genes']}**",
        f"- pares del universo: **{res['n_pairs']}**, positivos: "
        f"**{res['n_pos']}** (prevalencia {res['prevalence']:.2%})",
        "",
        "Validación cruzada 5-fold agrupada por complejo; métricas sobre "
        "predicciones fuera de fold. `sin_db` excluye el canal `database` de "
        "STRING, que contiene complejos curados: es la medición honesta.",
        "",
        "## Recuperación de pares y de partner",
        "",
        "| features | PR-AUC | ROC-AUC | R@1 | R@3 | R@5 | MRR |",
        "|---|---|---|---|---|---|---|",
    ]
    for key, label in (("genomica", "genómica"), ("sin_db", "genómica+STRING sin db"),
                       ("con_db", "genómica+STRING con db")):
        m, r = res[key], res[key + "_rank"]
        lines.append(f"| {label} | {m['pr_auc']} | {m['roc_auc']} | "
                     f"{r['recall@1']} | {r['recall@3']} | {r['recall@5']} | {r['mrr']} |")
    lines += [
        "",
        "`R@k`: para cada proteína del universo, fracción en que alguno de sus "
        "**k mejores partners** pertenece a su complejo curado.",
        "",
        "## ¿Importa si el gold es experimental?",
        "",
        "Top-1 partner dentro del mismo complejo, separando complejos con "
        "evidencia experimental (`ECO:0000353`) de los curados por inferencia:",
        "",
        "| evidencia del gold | genes | top-1 correcto |",
        "|---|---|---|",
    ]
    for cls in ("experimental", "inferida"):
        d = res.get("top1_by_evidence", {}).get(cls)
        if d:
            lines.append(f"| {cls} | {d['n']} | {d['top1']:.1%} |")
    lines += [
        "",
        "El modelo acierta más donde hay experimento que donde hay inferencia "
        "curatorial. No es un detalle: mitad del gold es inferencia.",
        "",
        "## Pares individuales (modelo sin db)",
        "",
        f"- precisión con score ≥0.5: **{res['sin_db']['precision@0.5']:.1%}**; "
        f"recall **{res['sin_db']['recall@0.5']:.1%}**",
        f"- precisión con score ≥0.9: **{res['sin_db']['precision@0.9']:.1%}**; "
        f"recall **{res['sin_db']['recall@0.9']:.1%}**",
        "",
        "La precisión por par es baja a propósito: STRING conecta proteínas de "
        "la misma ruta que no forman complejo. La especificidad no viene del par, "
        "viene del agrupamiento.",
        "",
        "## Clusters (top-2 mutuo sobre scores OOF sin db)",
        "",
        f"- clusters predichos: **{res['clusters']['n_clusters']}**",
        f"- pureza (todos sus pares son del mismo complejo): "
        f"**{res['clusters']['purity']:.1%}** ({res['clusters']['n_pure']})",
        f"- cobertura Jaccard ≥0.5: **{res['clusters']['coverage_j50']:.1%}**",
        f"- recuperación exacta: **{res['clusters']['exact']:.1%}**",
        f"- Jaccard medio: **{res['clusters']['mean_jaccard']:.3f}**",
        "",
        "| confianza | n | pureza |",
        "|---|---|---|",
    ]
    for k, d in res["clusters"].get("by_confidence", {}).items():
        lines.append(f"| {k} | {d['n']} | {d['purity']:.1%} |")
    lines += [
        "",
        "| tamaño del complejo | n | Jaccard medio |",
        "|---|---|---|",
    ]
    for size, d in res["clusters"]["by_size"].items():
        lines.append(f"| {size} | {d['n']} | {d['mean_jaccard']} |")
    lines += ["", "## Casos documentados", "",
              "| complejo | genes | J | recuperado | cluster predicho |",
              "|---|---|---|---|---|"]
    for e in res["examples"]:
        if e.get("status") == "ausente":
            continue
        lines.append(f"| {e['label']} | {e['n_genes']} | {e['best_jaccard']} | "
                     f"{'sí' if e['recovered'] else 'no'} | {e['cluster'][:70]} |")
    lines += [
        "",
        "## Límites",
        "",
        "- El universo es de mundo cerrado: solo genes que ya están en algún "
        "complejo curado. Un gen fuera de ComplexPortal no se evalúa.",
        "- La CV agrupa por complejo, pero un gen puede estar en dos complejos y "
        "aparecer en *train* y *test*. Mide recuperación, no generalización a "
        "genes nuevos.",
        "- `score` no es una probabilidad calibrada: el modelo entrena con clases "
        "balanceadas y la prevalencia real es ~0.6%.",
        "- Complejos de más de ~8 subunidades se fragmentan (el top-2 mutuo no "
        "puede mantener un anillo grande unido).",
        "- Las features genómicas dependen de operones predichos: si el operón "
        "está mal, la evidencia también.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=Path, default=L.CX / "k12_genes.tsv")
    ap.add_argument("--operons", type=Path,
                    default=L.DATA / "operon_finder" / "ecoli_k12_operons_rules.tsv")
    ap.add_argument("--string", type=Path,
                    default=L.CX / "raw" / "string_511145_full.txt.gz")
    ap.add_argument("--string-species", default="511145")
    ap.add_argument("--synteny", type=Path,
                    default=L.DATA / "operon_finder" / "k12_synteny.tsv")
    ap.add_argument("--gold", type=Path, default=L.CX / "gold_k12.tsv")
    ap.add_argument("--out-dir", type=Path, default=L.CX)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--save-oof", action="store_true")
    a = ap.parse_args(argv)

    for p in (a.genes, a.gold, a.operons):
        if not Path(p).exists():
            raise SystemExit(f"falta {p}; corre fetch_gold_k12.py primero")

    print("cargando")
    ctx, gold, _ = load_data(a)
    universe, pairs, y, groups, pair_of_cx = make_pairs(gold)
    gold_pairs = set(pair_of_cx)
    print(f"  {len(gold)} complejos, {len(universe)} genes, {len(pairs)} pares, "
          f"{int(y.sum())} positivos")

    sets = {
        "genomica": L.FEATURES_GENOMIC,
        "sin_db": [f for f in L.FEATURES_ALL if f != "string_db"],
        "con_db": L.FEATURES_ALL,
    }
    res: dict = {
        "n_complexes": len(gold), "n_genes": len(universe),
        "n_pairs": len(pairs), "n_pos": int(y.sum()),
        "prevalence": round(float(y.mean()), 5),
    }
    oof = {}
    for key, names in sets.items():
        print(f"  CV {key} ({len(names)} features)")
        X = [[L.features_row(g, h, ctx)[n] for n in names] for g, h in pairs]
        oof[key] = cv_scores(X, y, groups, names, a.folds)
        res[key] = pair_metrics(oof[key], y)
        res[key + "_rank"] = ranking_metrics(oof[key], pairs, gold_pairs, universe)
        print(f"    PR-AUC={res[key]['pr_auc']} R@1={res[key + '_rank']['recall@1']}")

    print("  clusters (OOF sin db)")
    edges = {L.pair_key(g, h): float(oof["sin_db"][i])
             for i, (g, h) in enumerate(pairs) if oof["sin_db"][i] >= 0.1}
    clusters = L.cluster_edges(edges, top_k=2, min_p=0.1)
    res["clusters"] = cluster_metrics(clusters, gold, edges)

    gold_rows = []
    with open(a.gold, newline="") as fh:
        gold_rows = list(csv.DictReader(fh, delimiter="\t"))
    res["examples"] = evaluate_examples(clusters, gold_rows, ctx.genes)
    res["top1_by_evidence"] = top1_by_evidence(oof["sin_db"], pairs, gold_rows)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "validation.json").write_text(
        json.dumps(res, indent=2, sort_keys=True) + "\n")
    (a.out_dir / "VALIDATION.md").write_text(markdown(res))
    if a.save_oof:
        out = a.out_dir / "oof_k12.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene_a", "gene_b", "score_sin_db", "score_con_db", "y"])
            for i, (g, h) in enumerate(pairs):
                w.writerow([g, h, round(oof["sin_db"][i], 4),
                            round(oof["con_db"][i], 4), int(y[i])])
        print(f"  escrito: {out}")

    print(f"  escrito: {a.out_dir / 'VALIDATION.md'}")
    print(f"  R@1 sin db = {res['sin_db_rank']['recall@1']}, "
          f"cobertura clusters = {res['clusters']['coverage_j50']}, "
          f"pureza = {res['clusters']['purity']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
