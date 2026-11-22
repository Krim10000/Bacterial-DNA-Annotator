#!/usr/bin/env python3
"""fetch_string.py — scores de asociación funcional de STRING por par adyacente.

La literatura (Taboada et al. 2010; Operon-mapper) llega a 94,6% de accuracy en
E. coli usando solo **distancia intergénica + score funcional de STRING**. Es el
feature que a nuestro modelo le falta: `same_ec` es un proxy muy pobre.

Se descargan los scores de todos los pares de genes adyacentes en la misma
hebra de K-12 y se cachean para no repetir red.

Salida (`data/operon_finder/string_scores.tsv`):
    pair_a, pair_b, score, nscore, fscore, pscore, ascore, escore, dscore, tscore

`score` es el combined score de STRING (0-1).

Uso
---
    python scripts/regulation/fetch_string.py
    python scripts/regulation/fetch_string.py --chunk 150 --force
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402

OUT = ROOT / "data" / "operon_finder" / "string_scores.tsv"
URL = "https://string-db.org/api/tsv/network"
SPECIES = "511145"
COLS = ["score", "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]


def sid(locus: str) -> str:
    return f"{SPECIES}.{locus}"


def fetch_chunk(ids: list[str], tries: int = 3) -> dict[tuple[str, str], dict]:
    params = urllib.parse.urlencode({"identifiers": "\r".join(ids),
                                     "species": SPECIES, "required_score": "0"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(f"{URL}?{params}", timeout=90) as fh:
                text = fh.read().decode()
            break
        except Exception:
            if attempt == tries - 1:
                return {}
            time.sleep(2 * (attempt + 1))
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    out = {}
    for r in rows:
        a = r["stringId_A"].split(".", 1)[-1]
        b = r["stringId_B"].split(".", 1)[-1]
        out[(a, b)] = {c: r.get(c, "") for c in COLS}
        out[(b, a)] = out[(a, b)]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gbff", type=Path, default=ROOT / "data" / "test_ecoli.gbff")
    ap.add_argument("--chunk", type=int, default=150)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    if OUT.exists() and not a.force:
        n = sum(1 for _ in open(OUT)) - 1
        print(f"cache existente: {n} pares en {OUT.relative_to(ROOT)} (usa --force)")
        return 0

    seq, genes = F.read_gbff(a.gbff)
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    pairs = [(x["locus_tag"], y["locus_tag"]) for x, y in F.pairs(genes)]
    print(f"{len(genes)} CDS, {len(pairs)} pares adyacentes misma hebra")

    scores: dict[tuple[str, str], dict] = {}
    for i in range(0, len(pairs), a.chunk):
        block = pairs[i:i + a.chunk]
        ids = sorted({sid(g) for pr in block for g in pr})
        got = fetch_chunk(ids)
        for pr in block:
            if pr in got:
                scores[pr] = got[pr]
        done = min(i + a.chunk, len(pairs))
        print(f"  [{done}/{len(pairs)}] pares con score: {len(scores)}", flush=True)
        time.sleep(0.3)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair_a", "pair_b"] + COLS)
        for (x, y), sc in sorted(scores.items()):
            w.writerow([x, y] + [sc[c] for c in COLS])

    cover = len(scores) / len(pairs) if pairs else 0
    print(f"\npares con score STRING: {len(scores)}/{len(pairs)} ({cover:.1%})")
    print(f"-> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
