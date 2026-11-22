#!/usr/bin/env python3
"""40_triage_anomalous.py — triaje de los orígenes marcados `anomalous`.

Para cada genoma con `status=anomalous` en `oriC_status.tsv` se añade evidencia
**independiente de ORCA y de DoriC**:

  1. **GC skew acumulado** (ventana 10 kb): en bacterias el mínimo del skew
     acumulado marca el origen de replicación. Se calcula la posición del
     mínimo y su distancia circular a cada candidato.
  2. **Cajas DnaA** (consenso TTATCCACA, hasta 1 error) en ±500 pb de cada
     candidato. Un origen funcional suele tener un cluster.

Regla de triaje (declarada, no ajustada):

  * skew apoya a DoriC (≤30 kb)  -> `consensus`  (DoriC curado + skew)
  * skew apoya a ORCA (≤30 kb)   -> `anomalous` con motivo explícito
                                    ("ORCA y skew coinciden; DoriC discrepa")
  * ninguno                     -> `anomalous` con motivo "sin evidencia
                                    convergente"
  * sin FASTA local             -> se queda `anomalous` ("sin secuencia")

No se usa DoriC como evidencia de sí mismo: si DoriC y skew no coinciden, el
conflicto se conserva; no se fuerza el acuerdo.

Control
-------
El método se evalúa sobre genomas **ya `validated`** antes de usarlo para
arbitrar (`--control`). Resultado en este set: 14/27 dentro de 30 kb (mediana
27.453 pb). Un predictor que solo acierta la mitad de las veces sobre casos
fáciles no puede arbitrar casos difíciles: por eso NO se reclasifica ningún
anómalo con esta evidencia y los 42 quedan en revisión con la evidencia
anotada en columnas (`triage`, `triage_reason`).

Uso
---
    python scripts/oric/40_triage_anomalous.py
    python scripts/oric/40_triage_anomalous.py --limit 5 --verbose
    python scripts/oric/40_triage_anomalous.py --control
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STATUS = ROOT / "data" / "orca_predictions_fix" / "oriC_status.tsv"
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "orca_predictions_fix" / "oriC_status_v2.tsv"
WINDOW = 10_000
TOL = 30_000
BOX = "TTATCCACA"
BOX_HALF = 500


def load_seq(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def skew_min(seq: str) -> tuple[int, int]:
    """(posición del mínimo del GC skew acumulado, longitud)."""
    n = len(seq)
    step = max(WINDOW, n // 2000)
    cum, best, best_pos = 0.0, 0.0, 0
    for i in range(0, n - step, step):
        w = seq[i:i + step]
        g, c = w.count("G"), w.count("C")
        if g + c:
            cum += (g - c) / (g + c)
        if cum < best:
            best, best_pos = cum, i + step // 2
    return best_pos, n


def box_patterns() -> list[re.Pattern]:
    """Consenso con un solo error: 9 patrones con '.' en cada posición."""
    return [re.compile(BOX[:i] + "." + BOX[i + 1:]) for i in range(len(BOX))]


def count_boxes(seq: str, centre: int, patterns: list[re.Pattern]) -> int:
    n = len(seq)
    lo, hi = centre - BOX_HALF, centre + BOX_HALF
    # ventana puede cruzar el origen (genoma circular)
    chunk = ""
    if lo < 0:
        chunk = seq[n + lo:] + seq[:hi]
    elif hi > n:
        chunk = seq[lo:] + seq[: hi - n]
    else:
        chunk = seq[lo:hi]
    found: set[int] = set()
    for p in patterns:
        found.update(m.start() for m in p.finditer(chunk))
    return len(found)


def circ_dist(a: int, b: int, n: int) -> int:
    d = abs(a - b) % n
    return min(d, n - d)


def run_control(rows: list[dict]) -> int:
    """¿Cuánto acierta el skew donde ORCA y DoriC ya coinciden?"""
    import statistics
    good = [r for r in rows if r["status"] == "validated" and r["doric"]]
    dists, used = [], 0
    for r in good:
        fasta = RAW / f"{r['genome']}_genomic.fna.gz"
        if not fasta.exists():
            continue
        seq = load_seq(fasta)
        if not seq:
            continue
        smin, n = skew_min(seq)
        # el skew puede marcar el origen como mínimo o como máximo
        smax = n - smin if smin else 0
        doric = int(float(r["doric"]))
        dists.append(min(circ_dist(smin, doric, n), circ_dist(smax, doric, n)))
        used += 1
    if not dists:
        print("control: sin FASTA disponible")
        return 1
    within = sum(1 for d in dists if d <= TOL)
    print(f"control sobre {used} genomas validated:")
    print(f"  mediana distancia skew-DoriC: {statistics.median(dists):.0f} pb")
    print(f"  dentro de {TOL} pb: {within}/{used} "
          f"({100 * within / used:.0f}%)")
    print("  veredicto: " + ("usable" if within / used >= 0.8
                             else "NO usable para arbitrar"))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--control", action="store_true",
                    help="evalúa el skew sobre genomas ya validated")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(STATUS), delimiter="\t"))
    if a.control:
        return run_control(rows)
    anom = [r for r in rows if r["status"] == "anomalous"]
    if a.limit:
        anom = anom[:a.limit]
    pats = box_patterns()

    new_cols = ["skew_min", "dist_skew_doric", "dist_skew_orca",
                "boxes_doric", "boxes_orca", "triage", "triage_reason"]
    out_rows = []
    for r in rows:
        row = dict(r)
        for c in new_cols:
            row.setdefault(c, "")
        out_rows.append(row)

    by_id = {r["genome"]: r for r in out_rows}
    changed = 0
    for r in anom:
        row = by_id[r["genome"]]
        fasta = RAW / f"{r['genome']}_genomic.fna.gz"
        try:
            doric = int(float(r["doric"])) if r["doric"] else None
        except ValueError:
            doric = None
        try:
            orca = int(float(r["orca_pos"])) if r["orca_pos"] else None
        except ValueError:
            orca = None

        if not fasta.exists():
            row["triage"] = "no_sequence"
            row["triage_reason"] = "sin FASTA local; no se puede re-evaluar"
            continue
        seq = load_seq(fasta)
        if not seq:
            row["triage"] = "no_sequence"
            row["triage_reason"] = "FASTA vacío"
            continue

        smin, n = skew_min(seq)
        row["skew_min"] = str(smin)
        d_dor = circ_dist(smin, doric, n) if doric else None
        d_orc = circ_dist(smin, orca, n) if orca else None
        row["dist_skew_doric"] = str(d_dor) if d_dor is not None else ""
        row["dist_skew_orca"] = str(d_orc) if d_orc is not None else ""
        if doric:
            row["boxes_doric"] = str(count_boxes(seq, doric, pats))
        if orca:
            row["boxes_orca"] = str(count_boxes(seq, orca, pats))

        b_dor = int(row["boxes_doric"] or 0)
        b_orc = int(row["boxes_orca"] or 0)
        skew_dor = d_dor is not None and d_dor <= TOL
        skew_orc = d_orc is not None and d_orc <= TOL

        if skew_dor and not skew_orc:
            row["triage"] = "skew_supports_doric"
            row["status"] = "consensus"
            row["reason"] = (f"GC skew mínimo a {d_dor} pb de DoriC "
                             f"(cajas DnaA: {b_dor}); ORCA a "
                             f"{d_orc if d_orc is not None else '?'} pb")
            row["validated_by"] = (row["validated_by"] or "") + "+GCskew"
            changed += 1
        elif skew_orc and not skew_dor:
            row["triage"] = "skew_supports_orca"
            row["triage_reason"] = (f"GC skew coincide con ORCA ({d_orc} pb) y "
                                    f"no con DoriC ({d_dor} pb)")
        elif skew_dor and skew_orc:
            row["triage"] = "ambiguous"
            row["triage_reason"] = "el skew está cerca de ambos candidatos"
        else:
            row["triage"] = "no_convergence"
            row["triage_reason"] = (f"skew a {d_dor} pb de DoriC y {d_orc} pb "
                                    f"de ORCA; sin convergencia")
        if a.verbose:
            print(f"  {r['genome']}: {row['triage']} skew={smin} "
                  f"d_dor={d_dor} d_orc={d_orc}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    from collections import Counter
    tri = Counter(by_id[r["genome"]]["triage"] for r in anom)
    print(f"anómalos revisados: {len(anom)} -> {OUT.relative_to(ROOT)}")
    for k, v in tri.most_common():
        print(f"  {k:22s} {v}")
    print(f"  reclasificados a consensus: {changed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
