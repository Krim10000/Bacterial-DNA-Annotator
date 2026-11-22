#!/usr/bin/env python3
"""rule_operons.py — buscador de operones DETERMINISTA, sin modelo estadístico.

Reemplaza al clasificador logístico por una **tabla de decisión** con reglas
explícitas, cada una con su motivo auditable. La idea central viene de la
evidencia: una horquilla **dentro** de un operón es un atenuador, no una
frontera, y su firma determinista es un **ORF corto de péptido líder** delante.

Reglas (por prioridad; la primera que aplica decide)
---------------------------------------------------
  R1  atenuador     el gen previo es un péptido líder (<=30 aa)  -> UNIR
                    (la horquilla corta el transcrito pero el operón sigue;
                     trpLEDCBA, hisLGDCBHAFI, leuLABCD son un solo operón)
  R2  terminador    horquilla clase A con poli-U>=4 a <=80 pb del fin del gen
                    previo, y NO seguida de ORF líder               -> CORTAR
  R3  promotor      PWM sigma70 en posición plausible (20-120 pb del inicio
                    del gen siguiente) con score alto               -> CORTAR
                    OJO: con el umbral por defecto (8.0) R3 dispara de más y
                    cuesta ~0,16 de F1 también en K-12; un promotor interno
                    fragmenta operones. Ver docs/PWM_RECALIBRATION.md.
  R4  distancia     hueco intergénico >= umbral (calibrado en K-12) -> CORTAR
  R5  defecto                                                        -> UNIR

Cada frontera sale con `rule` y `why`, así que la decisión es reproducible y
explicable sin pesos ni umbrales opacos.

Uso
---
    python scripts/regulation/rule_operons.py
    python scripts/regulation/rule_operons.py --gap 60 --out data/operon_finder
"""
from __future__ import annotations

import argparse
import bisect
import csv
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402

HAIRPINS = ROOT / "data" / "operon_finder" / "ecoli_k12_hairpins.tsv"
LEADER_MAX_AA = 35
STRING_COUPLED = 0.70
TERM_MAX_DIST = 80
PROM_WINDOW = (20, 120)
PROM_MIN = 8.0


def load_hairpins(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            out.append({"start": int(r["start"]), "end": int(r["end"]),
                        "class": r["class"], "polyU": int(r["polyU"] or 0),
                        "dG": float(r["dG"]) if r["dG"] else None})
        except (ValueError, KeyError):
            continue
    return sorted(out, key=lambda x: x["start"])


def hairpin_between(hps: list[dict], lo: int, hi: int) -> dict | None:
    """Horquilla cuyo inicio cae en [lo, hi]."""
    best = None
    for h in hps:
        if lo <= h["start"] <= hi:
            if best is None or h["class"] < best["class"]:
                best = h
    return best


def calibrate_gap(genes: list[dict], bins: int = 60) -> float:
    """Umbral de hueco por Otsu sobre log(hueco+1), SIN usar etiquetas.

    La distribución no tiene un valle limpio (no es bimodal), así que se usa
    Otsu: el corte que minimiza la varianza intra-clase. Es determinista y se
    calcula con el propio genoma, no con la verdad de referencia.
    """
    vals = [math.log1p(y["start"] - x["end"] - 1)
            for x, y in zip(genes, genes[1:])
            if y["strand"] == x["strand"] and y["start"] - x["end"] - 1 > 0]
    if len(vals) < 20:
        return 60.0
    lo, hi = min(vals), max(vals)
    w = (hi - lo) / bins
    hist = [0] * bins
    for v in vals:
        hist[min(int((v - lo) / w), bins - 1)] += 1
    tot = len(vals)
    best_var, best_t = -1.0, bins // 2
    for t in range(1, bins):
        w0 = sum(hist[:t])
        w1 = tot - w0
        if not w0 or not w1:
            continue
        m0 = sum((i + 0.5) * hist[i] for i in range(t)) / w0
        m1 = sum((i + 0.5) * hist[i] for i in range(t, bins)) / w1
        var = w0 * w1 * (m0 - m1) ** 2
        if var > best_var:
            best_var, best_t = var, t
    return round(math.expm1(lo + best_t * w), 1)


def string_score(a: dict, b: dict) -> float:
    return F.STRING_SCORES.get((a["locus_tag"], b["locus_tag"]), 0.0)


def decide(a: dict, b: dict, seq: str, hps: list[dict], boxes: list[dict],
           gap_max: int) -> tuple[bool, str, str]:
    """(¿cortar?, regla, motivo) para la frontera a|b."""
    gap = b["start"] - a["end"] - 1
    aa_prev = (a["end"] - a["start"]) // 3
    hp = hairpin_between(hps, a["end"], max(a["end"], b["start"]))

    # R1: péptido líder -> atenuador, el operón continúa
    if aa_prev <= LEADER_MAX_AA and gap < 400:
        why = f"gen previo {aa_prev} aa (péptido líder)"
        if hp:
            why += f"; horquilla {hp['class']} detrás = atenuador"
        return False, "R1_atenuador", why

    # R1b: horquilla + genes funcionalmente acoplados -> atenuación, no frontera
    if hp and hp["class"] == "A_fuerte" and string_score(a, b) >= STRING_COUPLED:
        return (False, "R1b_atenuador_funcional",
                f"horquilla {hp['class']} con STRING {string_score(a, b):.2f}: "
                f"terminación interna, el operón sigue")

    # R2: terminador fuerte pegado al gen previo
    if hp and hp["class"] == "A_fuerte" and hp["polyU"] >= 4 \
            and string_score(a, b) < STRING_COUPLED:
        d = hp["start"] - a["end"]
        if d <= TERM_MAX_DIST:
            return True, "R2_terminador", f"horquilla A a {d} pb, poli-U {hp['polyU']}"

    # R3: promotor en posición plausible
    up = seq[max(0, b["start"] - 150):max(0, b["start"] - 1)]
    pscore, pdist = F.promoter_scan(up, boxes)
    if pscore >= PROM_MIN and PROM_WINDOW[0] <= pdist <= PROM_WINDOW[1]:
        return True, "R3_promotor", f"sigma70 score {pscore:.1f} a {pdist} pb"

    # R4: distancia
    if gap >= gap_max:
        return True, "R4_distancia", f"hueco {gap} pb >= {gap_max}"

    return False, "R5_defecto", f"hueco {gap} pb, sin señal fuerte"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gbff", type=Path, default=ROOT / "data" / "test_ecoli.gbff")
    ap.add_argument("--genome", default="ecoli_k12")
    ap.add_argument("--gap", type=float, default=0,
                    help="umbral de distancia (pb); 0 = calibrar por Otsu")
    # R3 como frontera no ayuda: con el umbral por defecto (8) cuesta 0,19 de
    # F1 en K-12, y con 14 (5 cortes) sigue por debajo de apagarla del todo.
    # Los pocos cortes que sobreviven rompen operones largos. Default 16 =
    # R3 apagada; los promotores quedan como anotación, no como frontera.
    # Medido en docs/PWM_RECALIBRATION.md.
    ap.add_argument("--prom", type=float, default=16.0,
                    help="score mínimo de promotor (16 = R3 apagada)")
    ap.add_argument("--hairpins", type=Path, default=HAIRPINS)
    ap.add_argument("--tss", type=Path, default=None,
                    help="TSS (columnas pos,strand). Si se da, la frontera la decide "
                         "el promotor: R6_tss corta si hay un TSS en el hueco y las "
                         "horquillas dejan de usarse como frontera (medido: no "
                         "aportan cuando hay TSS)")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "operon_finder")
    ap.add_argument("--pwm", type=Path, default=None,
                    help="PWM de promotor alternativo (p. ej. B. subtilis); "
                         "por defecto usa el de E. coli")
    ap.add_argument("--grid", action="store_true", help="barre umbrales y reporta F1")
    a = ap.parse_args()
    a.out = a.out.resolve()          # permite --out relativo
    a.gbff = a.gbff.resolve()
    a.hairpins = a.hairpins.resolve()

    if a.pwm is not None:
        F.PROM_PWM = a.pwm
    boxes = F.load_promoter_pwm()
    global PROM_MIN
    # Con TSS documentados, la regla del promotor (R3) queda sustituida por
    # R6_tss: predecir el promotor cuando ya lo tenemos medido solo añade ruido.
    PROM_MIN = float("inf") if a.tss else a.prom
    F.load_string(a.out / "string_scores.tsv")
    seq, genes = F.read_gbff(a.gbff)
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    hps = {} if a.tss else load_hairpins(a.hairpins)
    tss: dict[str, list[int]] = {}
    if a.tss:
        for row in csv.DictReader(open(a.tss.resolve()), delimiter="\t"):
            try:
                pos = int(row.get("pos") or row.get("start") or 0)
            except ValueError:
                continue
            st = (row.get("strand") or "+").strip()
            tss.setdefault(st if st in ("+", "-") else "+", []).append(pos)
        for v in tss.values():
            v.sort()
        print(f"TSS cargados: {sum(len(v) for v in tss.values())} "
              f"(+ {len(tss.get('+', []))}, - {len(tss.get('-', []))})")
    if not a.gap:
        a.gap = calibrate_gap(genes)
        print(f"umbral calibrado por Otsu (sin etiquetas): {a.gap} pb")
    print(f"{a.genome}: {len(genes)} CDS, {len(hps)} horquillas, "
          f"umbral de hueco {a.gap} pb")

    ops: list[list[dict]] = []
    cur: list[dict] = []
    rules = Counter()
    detail = []
    for x, y in zip(genes, genes[1:]):
        same = x["strand"] == y["strand"] and x.get("seqid") == y.get("seqid")
        if not same:
            if cur:
                ops.append(cur)
            cur = [y]
            continue
        if not cur:
            cur = [x]
        if tss:
            stx = "+" if str(x.get("strand")) in ("1", "+") else "-"
            pos = tss.get(stx, [])
            i = bisect.bisect_left(pos, int(x["end"]) + 1)
            if i < len(pos) and pos[i] <= int(y["start"]) - 1:
                cut, rule, why = True, "R6_tss", f"TSS en el hueco ({pos[i]} pb)"
            else:
                cut, rule, why = decide(x, y, seq, hps, boxes, a.gap)
        else:
            cut, rule, why = decide(x, y, seq, hps, boxes, a.gap)
        rules[rule] += 1
        detail.append({"left": x["locus_tag"], "right": y["locus_tag"],
                       "rule": rule, "why": why, "cut": int(cut)})
        if cut:
            ops.append(cur)
            cur = [y]
        else:
            cur.append(y)
    if cur:
        ops.append(cur)

    # el ground truth de RegulonDB es de E. coli: en otra especie no aplica
    op_of = F.load_operon_of() if a.genome == "ecoli_k12" else {}
    if op_of and detail:
        p_, r_, f_ = F.pair_f1(ops, op_of)
        print(f"F1 de pares vs RegulonDB: {f_:.3f} (P {p_:.3f} R {r_:.3f})")

    out = a.out / f"{a.genome}_operons_rules.tsv"
    a.out.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["operon_id", "genome", "contig", "strand", "start", "end",
                    "n_genes", "genes", "gene_names", "confidence", "status"])
        for i, op in enumerate(ops, 1):
            w.writerow([f"{a.genome}:r{i}", a.genome, op[0].get("seqid", ""),
                        op[0]["strand"], min(g["start"] for g in op),
                        max(g["end"] for g in op), len(op),
                        ",".join(g["locus_tag"] for g in op),
                        ",".join(g.get("gene", "") or "" for g in op),
                        "rule", "candidate"])
    with open(a.out / f"{a.genome}_boundary_rules.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(detail[0]), delimiter="\t")
        w.writeheader()
        w.writerows(detail)

    try:
        destino = out.relative_to(ROOT)
    except ValueError:
        destino = out
    print(f"operones: {len(ops)} -> {destino}")
    print("decisiones por regla:")
    for r, n in rules.most_common():
        print(f"  {r:16s} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
