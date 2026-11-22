#!/usr/bin/env python3
"""fetch_dbtbs_promoters.py — TSS de DBTBS desde sus páginas de operón.

Por qué
-------
Para validar que **unir dos fuentes de TSS** mejora el resultado hacía falta una
segunda fuente para *B. subtilis*. BSGatlas ya está integrado; DBTBS publica sus
promotores en páginas por operón (`COG/prom/<operón>.html`) con:

- factor de transcripción (SigA, SigF…)
- `Location` relativa, p. ej. `-42:+4`
- **`Absolute position`** en el genoma, p. ej. `3849749..3849794`
- secuencia y referencia experimental

De ahí se obtiene el TSS: si la región cubre `a:b` respecto al inicio del
transcrito (+1), entonces `TSS = absoluto_inicio - a`. Con `-42:+4`:
`3849749 + 42 = 3849791`.

Es una **estimación** (±2 pb) y así se etiqueta en la salida: no sustituye a un
TSS medido, pero sí es una fuente independiente de BSGatlas.

Salida: `data/dbtbs/dbtbs_tss.tsv`  (pos, strand, sigma, evidence, operon)

Uso
---
    python scripts/validation/fetch_dbtbs_promoters.py [--limit 20]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "dbtbs" / "dbtbs_tss.tsv"
BASE = "https://dbtbs.hgc.jp/"
UA = {"User-Agent": "Mozilla/5.0 (academic use; operon benchmark)"}


def fetch(url: str) -> str:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=45) as fh:
        return fh.read().decode("latin-1", errors="replace")


def operon_pages() -> list[str]:
    """Enlaces a páginas de operón desde la tabla de promotores."""
    html = fetch(BASE + "promtable.html")
    links = re.findall(r'href="(COG/prom/[^"]+\.html)"', html)
    seen, out = set(), []
    for l in links:
        if l not in seen:
            seen.add(l)
            out.append(l)
    return out


def parse_page(html: str, operon: str):
    """Extrae (pos, strand, sigma, evidence) de cada promotor de la página."""
    # hebra: columna Direction de la tabla de genes ("+" o "-")
    m = re.search(r">Direction<.*?>\s*([+-])\s*<", html, re.S)
    strand = m.group(1) if m else "+"
    # bloques de fila con posición absoluta y localización relativa
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", html, re.S):
        txt = re.sub(r"<[^>]*>", "|", tr)
        txt = re.sub(r"\|+", "|", txt)
        cells = [c.strip() for c in txt.split("|") if c.strip()]
        abspos = None
        loc = None
        for c in cells:
            if re.fullmatch(r"\d{4,8}\.\.\d{4,8}", c):
                abspos = c
            if re.fullmatch(r"-?\d+:[+-]?\d+", c):
                loc = c
        if not (abspos and loc):
            continue
        a, b = abspos.split("..")
        rel_a = int(loc.split(":")[0])
        tss = int(a) - rel_a
        # factor: primera celda tipo SigA/SigF o nombre de proteína
        sigma = next((c for c in cells if re.fullmatch(r"Sig[A-Za-z0-9]+", c)), "")
        ev = next((c for c in cells if "(" in c and ")" in c and len(c) > 12), "")
        rows.append((tss, strand, sigma, ev[:80]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--delay", type=float, default=0.7)
    a = ap.parse_args()
    pages = operon_pages()
    if a.limit:
        pages = pages[:a.limit]
    print(f"páginas de operón: {len(pages)}")
    rows = []
    for i, page in enumerate(pages, 1):
        operon = page.split("/")[-1].replace(".html", "")
        try:
            rows += [(p, s, sig, ev, operon) for p, s, sig, ev in parse_page(fetch(BASE + page), operon)]
        except Exception as exc:  # noqa: BLE001
            print(f"  {operon}: {exc}")
        if i % 100 == 0:
            print(f"  {i}/{len(pages)} páginas, {len(rows)} promotores")
        time.sleep(a.delay)
    # dedupe por posición
    seen, uniq = set(), []
    for r in rows:
        if r[0] not in seen:
            seen.add(r[0])
            uniq.append(r)
    uniq.sort(key=lambda r: r[0])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "strand", "sigma", "evidence", "operon"])
        w.writerows(uniq)
    print(f"{len(uniq)} TSS de DBTBS -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
