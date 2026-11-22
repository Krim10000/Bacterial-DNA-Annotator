#!/usr/bin/env python3
"""F2 — análisis del co-folding: ¿el signo de Δ separa inductor de correpresor?

Entrada: data/allostery/boltz_results.json (salida de 04_boltz_modal.py).

Definiciones (Boltz-2):
    affinity_pred_value = log10(IC50) con IC50 en µM  -> más bajo = une mejor
    Δ = af(estado con ADN) − af(estado libre)
        Δ < 0  ⇒ el ligando prefiere el estado unido  ⇒ correpresor
        Δ > 0  ⇒ prefiere el estado libre             ⇒ inductor

Se reporta:
  - signo de Δ por TF y ligando (verdadero y señuelo)
  - acierto del signo contra el modo real
  - especificidad: Δ_verdadero − Δ_señuelo
  - test de signos (binomial) y sensibilidad con p(bind)

Uso: python scripts/allostery/05_boltz_analyze.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"


def binom_two_sided(k: int, n: int, p: float = 0.5) -> float:
    """p-valor exacto de un test de signos (dos colas)."""
    def pmf(i):
        return math.comb(n, i) * p ** i * (1 - p) ** (n - i)
    obs = pmf(k)
    return min(1.0, sum(pmf(i) for i in range(n + 1) if pmf(i) <= obs + 1e-12))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=str(OUT / "boltz_results.json"))
    ap.add_argument("--states", default="cognate,scrambled")
    ap.add_argument("--out", default="boltz_analysis")
    a = ap.parse_args()
    st_a, st_b = a.states.split(",")
    data = json.load(open(a.results))
    by = {}
    for r in data:
        m = r.get("meta") or {}
        aff = (r.get("affinity") or {}).get("affinity_pred_value")
        p = (r.get("affinity") or {}).get("affinity_probability_binary")
        by.setdefault((m.get("tf"), m.get("kind")), {})[m.get("state")] = {
            "aff": aff, "p": p, "rc": r.get("returncode")}

    rows = []
    for (tf, kind), st in sorted(by.items()):
        if st_a not in st or st_b not in st:
            continue
        va, vb = st[st_a]["aff"], st[st_b]["aff"]
        if va is None or vb is None:
            continue
        rows.append({"tf": tf, "kind": kind,
                     "mode": next(r["meta"]["mode"] for r in data
                                  if (r.get("meta") or {}).get("tf") == tf),
                     "aff_a": va, "aff_b": vb, "delta": va - vb,
                     "p_a": st[st_a]["p"], "p_b": st[st_b]["p"]})

    tfs = sorted({r["tf"] for r in rows})
    summary, detail = [], []
    for tf in tfs:
        t = next(r for r in rows if r["tf"] == tf and r["kind"] == "true")
        d = next((r for r in rows if r["tf"] == tf and r["kind"] == "decoy"), None)
        expected = "negativo" if t["mode"] == "correpresor" else "positivo"
        got = "negativo" if t["delta"] < 0 else "positivo"
        spec = (t["delta"] - d["delta"]) if d else None
        detail.append({"tf": tf, "mode": t["mode"],
                       "delta_true": round(t["delta"], 3),
                       "delta_decoy": round(d["delta"], 3) if d else None,
                       "especificidad": round(spec, 3) if spec is not None else None,
                       "esperado": expected, "observado": got,
                       "acierta": expected == got,
                       "p_true": [round(t["p_a"], 3), round(t["p_b"], 3)]})
    # Discriminación absoluta verdadero vs señuelo (la pregunta "¿este bolsillo
    # distingue el ligando correcto?"), independiente del estado.
    abs_rows = []
    for tf in tfs:
        t = next(r for r in rows if r["tf"] == tf and r["kind"] == "true")
        d = next((r for r in rows if r["tf"] == tf and r["kind"] == "decoy"), None)
        if not d:
            continue
        for est, lab in (("a", st_a), ("b", st_b)):
            ta, da = t[f"aff_{est}"], d[f"aff_{est}"]
            abs_rows.append({"tf": tf, "estado": est, "mejor": ta < da,
                             "gap": ta - da})
    hits = sum(1 for d in detail if d["acierta"])
    n = len(detail)
    pval = binom_two_sided(hits, n) if n else None
    summary = {"n_tf": n, "aciertos_signo": hits,
               "discrimina_abs": {
                   "n": len(abs_rows),
                   "aciertos": sum(1 for a in abs_rows if a["mejor"]),
                   "gap_medio": round(sum(a["gap"] for a in abs_rows) / max(1, len(abs_rows)), 3)},
               "p_signos": pval,
               "delta_medio_true": round(sum(abs(d["delta_true"]) for d in detail) / max(1, n), 3),
               "especificidad_media": round(
                   sum(d["especificidad"] for d in detail if d["especificidad"] is not None)
                   / max(1, sum(1 for d in detail if d["especificidad"] is not None)), 3)}

    (OUT / f"{a.out}.json").write_text(
        json.dumps({"summary": summary, "por_tf": detail}, indent=1, ensure_ascii=False))
    L = ["# F2 — co-folding Boltz-2: ¿Δ separa inductor de correpresor?", "",
         f"n={n} TFs · estados {st_a} vs {st_b} · aciertos de signo={hits}/{n} · p(signos)={pval}", "",
         "| TF | modo | Δ verdadero | Δ señuelo | especificidad | esperado | observado |",
         "|---|---|---:|---:|---:|---|---|"]
    for d in detail:
        L.append(f"| {d['tf']} | {d['mode']} | {d['delta_true']:+.3f} | "
                 f"{d['delta_decoy']:+.3f} | {d['especificidad']:+.3f} | "
                 f"{d['esperado']} | {d['observado']} |")
    L += ["", "## Discriminación absoluta (verdadero vs señuelo)", "",
          f"- aciertos: {summary['discrimina_abs']['aciertos']}/{summary['discrimina_abs']['n']}"
          f" · gap medio: {summary['discrimina_abs']['gap_medio']:+.3f} log10(IC50)",
          "", "Δ = afinidad(con ADN) − afinidad(libre); más bajo = une mejor. "
          "Δ<0 ⇒ correpresor; Δ>0 ⇒ inductor.", ""]
    (OUT / f"{a.out}.md").write_text("\n".join(L) + "\n")
    print(json.dumps({"summary": summary, "por_tf": detail}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
