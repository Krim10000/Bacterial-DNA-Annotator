#!/usr/bin/env python3
"""summarize_hairpins.py — resumen agregado y clase normalizada por GC.

El lote anotó horquillas en los 311 genomas con FASTA. Pero la clase A se
definió con `GC del tallo >= 0.60`, un umbral calibrado en *E. coli* (50,8% GC).
La correlación entre GC del genoma y proporción de clase A es r=0.87: en un
genoma de 27% GC casi ningún tallo llega a 0,60, así que la clase A
prácticamente desaparece. No es que no haya terminadores; es que el umbral no
porta.

Este script recalcula la clase con un umbral **relativo al fondo del genoma**
(`GC_tallo >= GC_genoma + 0.10`) y reporta ambas versiones. No re-escanea: las
medidas crudas (stem, gc, polyU) ya están guardadas, así que solo se re-etiqueta.

Salida: `data/operon_finder/summary.json` + tabla por consola.

Uso
---
    python scripts/regulation/summarize_hairpins.py
"""
from __future__ import annotations

import csv
import gzip
import json
import statistics
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "operon_finder"
RAW = ROOT / "data" / "raw"
MARGEN = 0.10          # el tallo debe superar el GC del genoma en 10 puntos


def gc_of(genome: str) -> float | None:
    path = RAW / f"{genome}_genomic.fna.gz"
    if not path.exists():
        return None
    try:
        s = "".join(l.strip() for l in gzip.open(path, "rt")
                    if not l.startswith(">")).upper()
        return (s.count("G") + s.count("C")) / len(s) if s else None
    except Exception:
        return None


def classify(row: dict, thr: float) -> str:
    stem = int(row.get("stem") or 0)
    gc = float(row.get("gc") or 0)
    polyu = int(row.get("polyU") or 0)
    if stem >= 6 and gc >= thr and polyu >= 4:
        return "A_fuerte"
    if stem >= 5 and polyu >= 3:
        return "B_media"
    return "C_debil"


def main() -> int:
    files = sorted(OUT.glob("*_hairpins.tsv"))
    genomes = [f.name.replace("_hairpins.tsv", "") for f in files]
    print(f"genomas con horquillas: {len(genomes)}")

    with Pool(8) as pool:
        gcs = dict(zip(genomes, pool.map(gc_of, genomes)))

    rows = []
    for f in files:
        g = f.name.replace("_hairpins.tsv", "")
        data = list(csv.DictReader(open(f), delimiter="\t"))
        if not data:
            continue
        gc = gcs.get(g)
        thr = (gc + MARGEN) if gc else 0.60
        norm = Counter(classify(r, thr) for r in data)
        orig = Counter(r["class"] for r in data)
        rows.append({"genome": g, "n": len(data),
                     "gc": round(gc, 4) if gc else None,
                     "thr": round(thr, 3),
                     "A_orig": orig["A_fuerte"],
                     "A_norm": norm["A_fuerte"],
                     "B_orig": orig["B_media"], "B_norm": norm["B_media"],
                     "C_orig": orig["C_debil"], "C_norm": norm["C_debil"]})

    total = sum(r["n"] for r in rows)
    a_orig = sum(r["A_orig"] for r in rows)
    a_norm = sum(r["A_norm"] for r in rows)
    print(f"horquillas totales: {total:,}")
    print(f"  clase A con umbral fijo 0.60 : {a_orig:,} ({100*a_orig/total:.1f}%)")
    print(f"  clase A normalizada por GC    : {a_norm:,} ({100*a_norm/total:.1f}%)")

    # correlación GC vs clase A (la que motivó la corrección)
    xs = [r["gc"] for r in rows if r["gc"]]
    ys = [100 * r["A_orig"] / r["n"] for r in rows if r["gc"]]
    if len(xs) > 2:
        mx, my = statistics.mean(xs), statistics.mean(ys)
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
        print(f"  correlación GC ↔ clase A (fija): r = {cov/den:.2f}")

    print(f"\n{'genoma':30s} {'GC':>5s} {'A fija%':>8s} {'A norm%':>8s}")
    con_gc = [r for r in rows if r["gc"]]
    for r in sorted(con_gc, key=lambda x: -x["gc"])[:6]:
        print(f"{r['genome'][:30]:30s} {100*(r['gc'] or 0):5.1f} "
              f"{100*r['A_orig']/r['n']:7.1f}% {100*r['A_norm']/r['n']:7.1f}%")
    print("  ...")
    for r in sorted(con_gc, key=lambda x: x["gc"])[:6]:
        print(f"{r['genome'][:30]:30s} {100*(r['gc'] or 0):5.1f} "
              f"{100*r['A_orig']/r['n']:7.1f}% {100*r['A_norm']/r['n']:7.1f}%")

    summary = {"genomes": len(rows), "hairpins": total,
               "margin": MARGEN,
               "A_fixed": a_orig, "A_normalized": a_norm,
               "pct_A_fixed": round(100 * a_orig / total, 2),
               "pct_A_normalized": round(100 * a_norm / total, 2),
               "per_genome": rows}
    if "--rewrite" in __import__("sys").argv:
        n = 0
        for f in files:
            g = f.name.replace("_hairpins.tsv", "")
            gc = gcs.get(g)
            thr = (gc + MARGEN) if gc else 0.60
            data = list(csv.DictReader(open(f), delimiter="\t"))
            if not data or "class_norm" in data[0]:
                continue
            for r in data:
                r["class_fixed"] = r.get("class", "")
                r["class"] = classify(r, thr)
            tmp = f.with_suffix(".tsv.tmp")
            with open(tmp, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t")
                w.writeheader(); w.writerows(data)
            tmp.replace(f); n += 1
        print(f"re-etiquetados {n} genomas (columna class = normalizada; "
              f"class_fixed conserva la original)")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print(f"\n-> {(OUT / 'summary.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
