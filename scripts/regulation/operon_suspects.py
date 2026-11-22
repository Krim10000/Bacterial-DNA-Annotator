#!/usr/bin/env python3
"""operon_suspects.py — sospecha de operón por proximidad + red metabólica.

Idea
----
Si N enzimas están **muy cerca en el ADN** y sus reacciones **pertenecen a la
misma red de reacciones** (comparten metabolitos, sin contar los ubicuos como
agua o ATP), eso es una pista de que forman un operón metabólico. No es prueba:
muchos operones no tienen ruta asociada y muchas rutas contiguas no son
operones. Por eso se marca como **sospecha**, nunca como hecho.

Se apoya en datos que ya existen: los `gene_map.tsv` (gen → reacción Rhea) y el
grafo de reacciones. Coste: segundos por genoma.

Salida (`data/metabolism/operon_suspects.tsv`):
    genome contig strand start end n_genes genes reactions shared_metabolites
    gap_mean confidence evidence

Uso
---
    python scripts/regulation/operon_suspects.py
    python scripts/regulation/operon_suspects.py --gap 60 --limit 5
"""
from __future__ import annotations

import argparse
import csv
import pickle
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENEMAPS = ROOT / "data" / "metabolism" / "genemaps"
REACTIONS = ROOT / "data" / "metabolism" / "rhea_reactions.pkl"
OUT = ROOT / "data" / "metabolism" / "operon_suspects.tsv"

# metabolitos ubicuos: compartirlos no indica relación funcional
UBIQUITOUS = {
    "CHEBI:15377", "CHEBI:15378", "CHEBI:15379", "CHEBI:16526", "CHEBI:30616",
    "CHEBI:456216", "CHEBI:456215", "CHEBI:57945", "CHEBI:57783", "CHEBI:57287",
    "CHEBI:58349", "CHEBI:57540",
}


def load_reactions() -> dict:
    with open(REACTIONS, "rb") as fh:
        return pickle.load(fh)


def genes_of(path: Path) -> list[dict]:
    """Una fila por gen (un gen con 2 reacciones se une)."""
    by_gene: dict[str, dict] = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        gid = r.get("gene_id") or r.get("locus_tag") or ""
        if not gid:
            continue
        g = by_gene.setdefault(gid, {
            "gene": gid, "product": r.get("product", ""),
            "contig": r.get("contig", ""), "start": int(r.get("start") or 0),
            "end": int(r.get("end") or 0), "strand": r.get("strand", "+"),
            "rhea": set(),
        })
        rid = r.get("rhea_id", "")
        if rid:
            g["rhea"].add(rid)
    return sorted(by_gene.values(), key=lambda g: (g["contig"], g["start"]))


def components(genes: list[dict], rx: dict) -> list[tuple[set[str], set[str]]]:
    """Componentes conexas del grafo de reacciones de un grupo de genes.

    Dos reacciones se conectan si comparten un metabolito no ubicuo.
    """
    all_rx = set().union(*(g["rhea"] for g in genes)) if genes else set()
    if len(all_rx) < 2:
        return []
    # metabolito -> reacciones del grupo que lo tocan
    by_metabolite: dict[str, set[str]] = defaultdict(set)
    for rid in all_rx:
        try:
            key = int(rid.split(":")[1])
        except (IndexError, ValueError):
            continue
        r = rx.get(key) or rx.get(str(key)) or {}
        for m in set(r.get("products", ())) | set(r.get("substrates", ())):
            if m not in UBIQUITOUS:
                by_metabolite[m].add(rid)
    # unión-find sencilla
    parent = {rid: rid for rid in all_rx}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for mets, rids in by_metabolite.items():
        rids = sorted(rids)
        for a, b in zip(rids, rids[1:]):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
    groups: dict[str, set[str]] = defaultdict(set)
    for rid in all_rx:
        groups[find(rid)].add(rid)
    out = []
    for root, rids in groups.items():
        if len(rids) < 2:
            continue
        shared = {m for m, rs in by_metabolite.items() if len(rs & rids) >= 2}
        out.append((rids, shared))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gap", type=int, default=60, help="distancia máxima entre genes (pb)")
    ap.add_argument("--min-genes", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    rx = load_reactions()
    maps = sorted(GENEMAPS.glob("*_gene_map.tsv"))
    if a.limit:
        maps = maps[:a.limit]
    print(f"genomas: {len(maps)} | umbral de hueco {a.gap} pb")

    rows = []
    per_genome = {}
    for path in maps:
        genome = path.name.replace("_gene_map.tsv", "")
        genes = genes_of(path)
        found = 0
        # correr por contig+hebra y agrupar por proximidad
        for i, g in enumerate(genes):
            if i and (genes[i - 1]["contig"] != g["contig"]
                      or genes[i - 1]["strand"] != g["strand"]):
                continue
            run = [g]
            j = i + 1
            while j < len(genes):
                prev, nxt = genes[j - 1], genes[j]
                if (nxt["contig"] != g["contig"] or nxt["strand"] != g["strand"]):
                    break
                if nxt["start"] - prev["end"] - 1 > a.gap:
                    break
                run.append(nxt)
                j += 1
            if len(run) < a.min_genes:
                continue
            if not all(x["rhea"] for x in run):
                continue
            comps = components(run, rx)
            if not comps:
                continue
            rids, shared = max(comps, key=lambda c: len(c[0]))
            involved = [x for x in run if x["rhea"] & rids]
            if len(involved) < 2:
                continue
            gaps = [run[k + 1]["start"] - run[k]["end"] - 1 for k in range(len(run) - 1)]
            found += 1
            rows.append({
                "genome": genome, "contig": g["contig"], "strand": g["strand"],
                "start": min(x["start"] for x in run),
                "end": max(x["end"] for x in run),
                "n_genes": len(run),
                "genes": ",".join(x["gene"] for x in run),
                "products": "|".join(x["product"][:40] for x in run),
                "reactions": ",".join(sorted(rids)),
                "shared_metabolites": ",".join(sorted(shared)[:6]),
                "gap_mean": round(statistics.mean(gaps), 1) if gaps else 0,
                "confidence": "high" if (len(involved) >= 3 and shared) else "medium",
                "evidence": (f"{len(involved)} genes contiguos (hueco medio "
                             f"{round(statistics.mean(gaps),1) if gaps else 0} pb) "
                             f"con {len(shared)} metabolitos compartidos"),
            })
        per_genome[genome] = found

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"sospechas: {len(rows)} en {len(per_genome)} genomas")
    if per_genome:
        vals = [v for v in per_genome.values() if v]
        print(f"  genomas con al menos una: {len(vals)} | mediana por genoma: "
              f"{statistics.median(vals) if vals else 0:.0f}")
    print(f"-> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
