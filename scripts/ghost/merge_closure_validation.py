#!/usr/bin/env python3
"""merge_closure_validation.py — clasifica acuerdos y desacuerdos entre cierres.

La pregunta no es "cuántos cierres nuevos hay" sino "qué evidencia sostiene cada
uno". Este script cruza el lote de dominio nuevo con `closure_annotations.tsv`
(el catálogo previo) y etiqueta cada brecha con el vocabulario del proyecto
(`docs/PIPELINE_LOGIC.md`):

    reproduced_same_protein   mismo dominio y misma proteína  -> consensus
    reproduced_other_protein  mismo dominio, otra proteína    -> consensus
    disagreement_domain       dominio distinto                -> anomalous
    new                       sin cierre previo               -> candidate
    prior_only                solo el método previo la cerró  -> (se cuenta aparte)

Por qué `consensus` y no `validated`: los cierres previos salen del **mismo**
método de dominio (`source=dominio` en 373 de 387). Reproducirlo con un panel
Pfam reconstruido (2026) demuestra determinismo y equivalencia de panel, no
independencia de evidencia. `validated` exige dos métodos independientes; eso
sería cruzar con un método no-dominio (p. ej. homología DIAMOND directa).

Los desacuerdos son el subproducto más valioso: no se descartan, se listan para
revisión con el dominio de cada lado y el e-value que lo sostiene.

Uso
---
    python scripts/ghost/merge_closure_validation.py --top 1000
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data" / "metabolism"
NEW = META / "domain_closure_candidates_top{top}.tsv"
OLD = META / "closure_annotations.tsv"
OUT = META / "closure_validation_top{top}.tsv"
CANDIDATE = {"candidate", "candidate_weak"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=1000)
    a = ap.parse_args()

    new_path = Path(str(NEW).format(top=a.top))
    rows = [r for r in csv.DictReader(open(new_path), delimiter="\t")
            if r["status"] in CANDIDATE]
    old = list(csv.DictReader(open(OLD), delimiter="\t"))

    old_map: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in old:
        old_map[(r["genome"], r["gap"])].append(r)
    new_map: dict[tuple[str, str], dict] = {}
    for r in rows:
        new_map[(r["genome"], r["gap"])] = r  # una fila por brecha (mejor e-value)

    out = []
    for key, n in sorted(new_map.items()):
        o = old_map.get(key)
        if not o:
            verdict, detail = "new", ""
        else:
            same_dom = [x for x in o if x["domain"] == n["domain"]]
            same_all = [x for x in same_dom if x["query"] == n["query"]]
            if same_all:
                verdict, detail = "reproduced_same_protein", same_all[0]["evalue"]
            elif same_dom:
                verdict = "reproduced_other_protein"
                detail = f"{same_dom[0]['query']}|{same_dom[0]['evalue']}"
            else:
                verdict = "disagreement_domain"
                detail = f"{o[0]['domain']}|{o[0]['query']}|{o[0]['evalue']}"
        out.append({
            "genome": n["genome"], "gap": n["gap"], "role": n["role"],
            "support_400": n["support_400"],
            "new_domain": n["domain"], "new_query": n["query"],
            "new_evalue": n["evalue"], "direction": n["direction"],
            "prior": detail, "verdict": verdict,
            "status_map": {"new": "candidate",
                           "disagreement_domain": "anomalous"}.get(verdict, "consensus"),
            "evidence": (f"dominio={n['domain']}; query={n['query']}; "
                         f"evalue={n['evalue']}; direccion={n['direction']}; "
                         f"veredicto={verdict}"),
        })

    path = Path(str(OUT).format(top=a.top))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    c = Counter(r["verdict"] for r in out)
    prior_only = set(old_map) - set(new_map)
    print(f"brechas con candidato nuevo: {len(out)} -> {path.relative_to(ROOT)}")
    for k, v in c.most_common():
        print(f"  {k:26s} {v}")
    print(f"  {'prior_only':26s} {len(prior_only)}")
    print(f"\nmapeo al vocabulario:")
    for k, v in Counter(r["status_map"] for r in out).most_common():
        print(f"  {k:10s} {v}")
    print("\ndesacuerdos de dominio (revisión):")
    for r in [x for x in out if x["verdict"] == "disagreement_domain"][:8]:
        print(f"  {r['gap'][:34]:34s} nuevo={r['new_domain'][:22]:22s} "
              f"previo={r['prior'].split('|')[0][:22]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
