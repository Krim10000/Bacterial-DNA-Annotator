#!/usr/bin/env python3
"""build_tf_pwms.py — PWMs de factores de transcripción desde sitios reales.

Fuente: `data/regulondb/RISet.tsv`, que trae por interacción el **sitio de
unión** (`tfrsSeq`), sus coordenadas y la función (activador/represor). Son
sitios curados de RegulonDB, no motivos de manual.

Salida (`data/regulondb/gt/tf_pwms.tsv`): una fila por TF con
  tf, n_sites, width, consensus, pwm (log-odds serializado), info (bits)

Por qué una PWM y no un motivo con huecos: la capa L2 usa hoy un patrón
`TTGACA`/`TATAAT` con umbral de mismatches, y su precisión posicional medida
contra `gt_promoters.tsv` a ±10 nt es P=0.31 / R=0.19. Una PWM puntúa cada
posición y permite elegir el máximo local, que es lo que hace falta para saber
si un operador solapa el -10 (de eso depende el signo de una arista de L3).

Filtros (declarados):
  - >= 8 sitios por TF (menos no da matriz estable)
  - sitios solo ACGT, longitud entre 6 y 40
  - se recorta/centra a la longitud modal del TF

Además escribe `promoter_pwms.tsv`: PWM de las cajas -10 y -35 por factor
sigma, construidas desde `gt_promoters.tsv` (posiciones y secuencias curadas).
Eso es lo que permite sustituir el escaneo por umbral de mismatches por una
puntuación posicional, que es la métrica donde L2 pierde (P=0.31 a ±10 nt).

Uso
---
    python scripts/validation/build_tf_pwms.py
    python scripts/validation/build_tf_pwms.py --min-sites 15 --test
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RISET = ROOT / "data" / "regulondb" / "RISet.tsv"
OUT = ROOT / "data" / "regulondb" / "gt" / "tf_pwms.tsv"
BASES = "ACGT"
PSEUDO = 0.5


def read_sites(path: Path) -> dict[str, list[str]]:
    sites: dict[str, list[str]] = defaultdict(list)
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = [c.split(")")[-1].strip() for c in rows[0]]
    idx_tf = header.index("regulatorName") if "regulatorName" in header else None
    idx_seq = header.index("tfrsSeq") if "tfrsSeq" in header else None
    if idx_tf is None or idx_seq is None:
        raise SystemExit(f"faltan columnas en {path}: {header}")
    for row in rows[1:]:
        if len(row) <= max(idx_tf, idx_seq):
            continue
        tf = row[idx_tf].strip()
        seq = row[idx_seq].strip().upper()
        if not tf or not seq:
            continue
        if any(b not in BASES for b in seq) or not (6 <= len(seq) <= 40):
            continue
        sites[tf].append(seq)
    return sites


def modal_width(seqs: list[str]) -> int:
    return Counter(len(s) for s in seqs).most_common(1)[0][0]


def centred(seq: str, width: int) -> str:
    """Centra el sitio a `width`, recortando o rellenando con N."""
    if len(seq) > width:
        extra = len(seq) - width
        lo = extra // 2
        return seq[lo:lo + width]
    padr = width - len(seq)
    lo = padr // 2
    return "N" * lo + seq + "N" * (padr - lo)


def build_pwm(seqs: list[str], width: int) -> list[dict[str, float]]:
    counts = [Counter() for _ in range(width)]
    n = 0
    for s in seqs:
        c = centred(s, width)
        if "N" in c:
            continue
        n += 1
        for i, b in enumerate(c):
            counts[i][b] += 1
    if not n:
        return []
    pwm = []
    for i in range(width):
        row = {}
        for b in BASES:
            p = (counts[i][b] + PSEUDO) / (n + 4 * PSEUDO)
            row[b] = round(math.log2(p / 0.25), 4)
        pwm.append(row)
    return pwm


def consensus(pwm: list[dict[str, float]]) -> str:
    return "".join(max(row, key=row.get) for row in pwm)


def information(pwm: list[dict[str, float]]) -> float:
    total = 0.0
    for row in pwm:
        for b in BASES:
            p = 0.25 * 2 ** row[b]
            total += p * row[b]
    return round(total, 3)


def score(pwm: list[dict[str, float]], seq: str) -> float:
    if len(seq) != len(pwm):
        return float("-inf")
    return sum(pwm[i].get(b, -10.0) for i, b in enumerate(seq))


def self_test(pwm: list[dict[str, float]], seqs: list[str], n: int = 500) -> float:
    """AUC de sitios reales vs aleatorios de igual composición."""
    width = len(pwm)
    rnd = random.Random(42)
    pool = "".join(seqs)
    real = [score(pwm, centred(s, width)) for s in seqs if "N" not in centred(s, width)]
    fake = [score(pwm, "".join(rnd.choice(pool) for _ in range(width)))
            for _ in range(n)]
    if not real or not fake:
        return float("nan")
    wins = sum(1 for r in real for f in fake if r > f)
    ties = sum(1 for r in real for f in fake if r == f)
    return (wins + 0.5 * ties) / (len(real) * len(fake))


def build_promoter_pwms() -> None:
    """PWM de -10 y -35 por factor sigma desde las cajas curadas."""
    src = ROOT / "data" / "regulondb" / "gt" / "gt_promoters.tsv"
    out = ROOT / "data" / "regulondb" / "gt" / "promoter_pwms.tsv"
    if not src.exists():
        print("  (sin gt_promoters.tsv; se omite promoter_pwms)")
        return
    boxes: dict[tuple[str, str], list[str]] = defaultdict(list)
    with open(src, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sigma = (row.get("sigma_factor") or "unknown").strip()
            for name, key in (("-10", "minus10_seq"), ("-35", "minus35_seq")):
                seq = (row.get(key) or "").strip().upper()
                if seq and all(b in BASES for b in seq):
                    boxes[(sigma, name)].append(seq)
    rows = []
    for (sigma, name), seqs in sorted(boxes.items()):
        if len(seqs) < 8:
            continue
        width = modal_width(seqs)
        pwm = build_pwm(seqs, width)
        if not pwm:
            continue
        rows.append({"sigma": sigma, "box": name, "n_sites": len(seqs),
                     "width": width, "consensus": consensus(pwm),
                     "info_bits": information(pwm),
                     "pwm": json.dumps([[p[b] for b in BASES] for p in pwm])})
    if not rows:
        return
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  PWMs de promotor: {len(rows)} -> {out.relative_to(ROOT)}")
    for r in rows:
        print(f"    {r['sigma']:10s} {r['box']:4s} n={r['n_sites']:4d} "
              f"{r['consensus']} ({r['info_bits']} bits)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--min-sites", type=int, default=8)
    ap.add_argument("--test", action="store_true", help="AUC sitios vs aleatorio")
    ap.add_argument("--promoters", action="store_true",
                    help="construye también las PWM de cajas -10/-35 por sigma")
    a = ap.parse_args()

    if not RISET.exists():
        print(f"falta {RISET}", file=sys.stderr)
        return 1
    sites = read_sites(RISET)

    rows = []
    for tf, seqs in sorted(sites.items()):
        if len(seqs) < a.min_sites:
            continue
        width = modal_width(seqs)
        pwm = build_pwm(seqs, width)
        if not pwm:
            continue
        row = {"tf": tf, "n_sites": len(seqs), "width": width,
               "consensus": consensus(pwm), "info_bits": information(pwm),
               "pwm": json.dumps([[p[b] for b in BASES] for p in pwm])}
        if a.test:
            row["auc_selftest"] = round(self_test(pwm, seqs), 3)
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"TFs con PWM: {len(rows)} (>= {a.min_sites} sitios) -> {OUT.relative_to(ROOT)}")
    print(f"  sitios usados: {sum(r['n_sites'] for r in rows)} "
          f"de {sum(len(v) for v in sites.values())}")
    if a.test:
        aucs = [r["auc_selftest"] for r in rows if r.get("auc_selftest") == r.get("auc_selftest")]
        if aucs:
            print(f"  AUC self-test: mediana {statistics.median(aucs):.3f} "
                  f"(min {min(aucs):.3f})")
    if a.promoters:
        build_promoter_pwms()

    cap = [r for r in rows if r["info_bits"] >= 6]
    print(f"  con >=6 bits de información: {len(cap)}")
    for r in rows[:5]:
        print(f"    {r['tf']:12s} n={r['n_sites']:3d} w={r['width']:2d} "
              f"{r['consensus']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
