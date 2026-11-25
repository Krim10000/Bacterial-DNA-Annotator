#!/usr/bin/env python3
"""PWMs de sigma específicas de B. subtilis desde TSS curados.

No transfiere matrices de E. coli. Para cada etiqueta de DBTBS/BSGatlas
(A/B/D/I/L) construye una matriz de una caja con las ventanas aguas arriba
del TSS. Las etiquetas con pocos TSS se marcan con su n y no se presentan
como equivalentes funcionales perfectos de los sigma de E. coli.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TSS = ROOT / "data" / "dbtbs" / "gt_tss.tsv"
FASTA = ROOT / "data" / "dbtbs" / "genome.fna"
OUT = ROOT / "data" / "dbtbs" / "promoter_pwms_subtilis_all.tsv"
BASES = "ACGT"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
LABELS = {"A": "sigma70", "B": "sigma38", "D": "sigma28",
          "I": "sigma32", "L": "sigma54"}


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">")).upper()


def read_tss(sigma: str) -> list[tuple[int, str]]:
    out = []
    for r in csv.DictReader(open(TSS), delimiter="\t"):
        if (r.get("sigma") or "").strip() != sigma:
            continue
        try:
            pos = int(r["pos"])
        except (KeyError, ValueError):
            continue
        strand = (r.get("strand") or "+").strip()
        if pos >= 1 and strand in ("+", "-"):
            out.append((pos, strand))
    return out


def upstream(genome: str, pos: int, strand: str, width: int) -> str | None:
    if strand == "+":
        a, b = pos - 1 - width, pos - 1
        return genome[a:b] if a >= 0 else None
    a, b = pos, pos + width
    return genome[a:b].translate(COMP)[::-1] if b <= len(genome) else None


def build(windows: list[str], lo: int, hi: int, bg: dict[str, float]) -> dict:
    width = hi - lo
    counts = [[0] * 4 for _ in range(width)]
    for w in windows:
        for i, ch in enumerate(w[lo:hi]):
            if ch in BASES:
                counts[i][BASES.index(ch)] += 1
    n = len(windows)
    pwm, probs = [], []
    for i in range(width):
        col = [(counts[i][j] + 0.5) / (n + 2.0) for j in range(4)]
        probs.append(col)
        pwm.append([round(math.log2(col[j] / bg[BASES[j]]), 4) for j in range(4)])
    cons = "".join(BASES[max(range(4), key=lambda j: counts[i][j])] for i in range(width))
    info = round(sum(2.0 + sum(p * math.log2(p) for p in col if p > 0)
                     for col in probs), 3)
    return {"width": width, "n_sites": n, "consensus": cons,
            "info_bits": info, "pwm": pwm}


def best_box(windows: list[str], bg: dict[str, float], up: int) -> tuple[dict, int]:
    best = None
    # Región del promotor: no buscar en cualquier parte del gen.
    for width in (6, 7, 8, 9, 10, 12):
        for lo in range(max(0, up - 30), max(1, up - 8)):
            if lo + width > up:
                continue
            box = build(windows, lo, lo + width, bg)
            if best is None or box["info_bits"] > best[0]["info_bits"]:
                best = (box, lo)
    return best  # type: ignore[return-value]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--up", type=int, default=45)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    genome = read_fasta(FASTA)
    bg = {b: genome.count(b) / len(genome) for b in BASES}
    rows = []
    for label, canonical in LABELS.items():
        tss = read_tss(label)
        windows = [w for w in (upstream(genome, p, s, a.up) for p, s in tss) if w]
        if len(windows) < 5:
            print(f"{label}/{canonical}: omitido, solo {len(windows)} ventanas")
            continue
        box, lo = best_box(windows, bg, a.up)
        row = {"sigma": canonical, "box": "-10", "n_sites": box["n_sites"],
               "width": box["width"], "consensus": box["consensus"],
               "info_bits": box["info_bits"], "pwm": json.dumps(box["pwm"])}
        rows.append(row)
        print(f"{label}/{canonical}: n={len(windows)} "
              f"box={box['consensus']} info={box['info_bits']}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
