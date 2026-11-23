#!/usr/bin/env python3
"""annotate_evidence.py — ¿este candidato ya está documentado o es sospecha?

Toma los clusters del buscador y los compara con un catálogo curado
(ComplexPortal mapeado a locus_tag). Para cada candidato declara:

- `documentado`   : reproduce el complejo curado (Jaccard >= 0.8);
- `subcomplejo`   : todos sus genes pertenecen a un mismo complejo curado,
                    pero el candidato es un trozo (Jaccard < 0.8);
- `parcial`       : al menos la mitad de sus genes están en un complejo curado
                    (típico de operones que codifican dos enzimas distintas);
- `nuevo`         : no se parece a nada del catálogo (sospecha).

Y separa la evidencia del catálogo: `experimental` (ECO:0000353) vs `inferida`
(ortología, parálogos, conocimiento previo). No todo lo curado es experimental.

Esto es **anotación, no validación**: si el catálogo se usó para entrenar el
modelo, el solapamiento está garantizado. Para validar, `validate_k12.py`.

Uso
---
    python scripts/complexes/annotate_evidence.py
    python scripts/complexes/annotate_evidence.py \
        --clusters data/complexes/clusters_ecoli_k12.tsv \
        --gold data/complexes/gold_k12.tsv \
        --out data/complexes/clusters_ecoli_k12_evidencia.tsv
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import complexes_lib as L  # noqa: E402


def read_gold(path: Path) -> list[dict]:
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tags = {g for g in (r.get("locus_tags") or "").split(",") if g}
            if len(tags) < 2:
                continue
            out.append({
                "id": r["complex_id"],
                "name": r["complex_name"],
                "genes": tags,
                "evidence": r.get("evidence_class") or "otra",
                "eco": r.get("evidence_code") or "",
            })
    return out


def classify(cluster: set[str], gold: list[dict]) -> dict:
    """Mejor solapamiento; `within` = fracción del cluster dentro del curado."""
    best = {"jaccard": 0.0, "within": 0.0, "id": "", "name": "",
            "evidence": "", "eco": "", "n_shared": 0}
    for g in gold:
        inter = len(cluster & g["genes"])
        if not inter:
            continue
        union = len(cluster | g["genes"])
        jac = inter / union
        within = inter / len(cluster)
        # se prioriza Jaccard, luego cobertura del cluster
        if (jac, within) > (best["jaccard"], best["within"]):
            best = {"jaccard": jac, "within": within, "id": g["id"],
                    "name": g["name"], "evidence": g["evidence"],
                    "eco": g["eco"], "n_shared": inter}
    j, w = best["jaccard"], best["within"]
    if j >= 0.8:
        status = "documentado"
    else:
        # ¿es un trozo de un complejo curado, aunque el mejor Jaccard sea otro?
        full = next((g for g in gold if cluster <= g["genes"]), None)
        if full is not None:
            best = {"jaccard": j, "within": 1.0, "id": full["id"],
                    "name": full["name"], "evidence": full["evidence"],
                    "eco": full["eco"], "n_shared": len(cluster)}
            status = "subcomplejo"
        elif w >= 0.5:
            status = "parcial"
        else:
            status = "nuevo"
    best["status"] = status
    return best


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", type=Path,
                    default=L.CX / "clusters_ecoli_k12.tsv")
    ap.add_argument("--gold", type=Path, default=L.CX / "gold_k12.tsv")
    ap.add_argument("--out", type=Path,
                    default=L.CX / "clusters_ecoli_k12_evidencia.tsv")
    ap.add_argument("--summary", type=Path,
                    default=L.CX / "evidence_summary.json")
    a = ap.parse_args(argv)

    if not a.clusters.is_file():
        raise SystemExit(f"falta {a.clusters}; corre build_candidates.py --apply")
    if not a.gold.is_file():
        raise SystemExit(f"falta {a.gold}; corre fetch_gold_k12.py")

    gold = read_gold(a.gold)
    print(f"catálogo: {len(gold)} complejos "
          f"({Counter(g['evidence'] for g in gold)})")

    rows: list[dict] = []
    with open(a.clusters, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = {g for g in (r.get("genes") or "").split(",") if g}
            if len(genes) < 2:
                continue
            m = classify(genes, gold)
            r.update({
                "evidencia": m["status"],
                "curado_id": m["id"],
                "curado_nombre": m["name"],
                "jaccard": round(m["jaccard"], 3),
                "cobertura_cluster": round(m["within"], 3),
                "curado_evidencia": m["evidence"],
            })
            rows.append(r)

    cols = list(rows[0]) if rows else []
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # Resumen: evidencia × confianza del modelo
    table: dict[str, Counter] = {}
    for r in rows:
        table.setdefault(r["evidencia"], Counter())[r.get("confidence", "?")] += 1
    summary = {
        "n_clusters": len(rows),
        "por_evidencia": {k: dict(v) for k, v in sorted(table.items())},
        "total_por_evidencia": {k: sum(v.values()) for k, v in sorted(table.items())},
    }
    a.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"escrito: {a.out} ({len(rows)} clusters)")
    print(f"{'evidencia':<14}{'alta':>6}{'media':>7}{'baja':>6}{'total':>7}")
    for k, v in sorted(table.items()):
        print(f"{k:<14}{v.get('alta', 0):>6}{v.get('media', 0):>7}"
              f"{v.get('baja', 0):>6}{sum(v.values()):>7}")
    nuevos = [r for r in rows if r["evidencia"] == "nuevo"]
    print(f"nuevos: {len(nuevos)} (de ellos, alta: "
          f"{sum(1 for r in nuevos if r.get('confidence') == 'alta')})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
