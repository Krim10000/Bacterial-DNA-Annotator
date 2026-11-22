#!/usr/bin/env python3
"""build_domain_worklist.py — alcance del lote de cierre por dominio (top N).

Prepara, para las brechas de mayor soporte, todo lo que NO requiere el escaneo
HMM: productores Rhea, accesiones UniProt, dominios Pfam y **dirección**
(LR/RL/UN/BI + master) desde `rhea2uniprot_sprot.tsv`.

Por qué existe la separación
----------------------------
El escaneo de dominios necesita Pfam-A.hmm (658 MB) y en esta máquina no había
ninguna copia: el panel local (`data/db/pfam/panel_focus.hmm`) son 18 modelos
de sulfatasas/timidilato, ninguno de los 12 dominios que pide este lote. Este
script deja el lote listo y marca cada fila con su estado real, en vez de
inventar cierres.

Estados
-------
  pending_hmm   tiene dominios pero faltan los HMM (bloqueo externo)
  no_pfam       los productores no tienen dominio Pfam anotado
  no_accession  ninguna reacción productora tiene accesión UniProt
  no_producer   la brecha no tiene reacción productora en Rhea

Uso
---
    python scripts/ghost/build_domain_worklist.py --top 1000 --resolve
    python scripts/ghost/build_domain_worklist.py --top 50        # sin red
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
import re
import sys
from math import ceil
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLIST = ROOT / "data" / "metabolism" / "gap_worklist.tsv"
REACTIONS = ROOT / "data" / "metabolism" / "rhea_reactions.pkl"
R2A = ROOT / "data" / "metabolism" / "rhea2uniprot_sprot.tsv"
CACHE = ROOT / "data" / "metabolism" / "uniprot_pfam_cache.json"
OUT = ROOT / "data" / "metabolism" / "domain_worklist_top{top}.tsv"
PRIORITY = {"LR": 3, "RL": 3, "UN": 2, "BI": 1, "": 0}


def load_producers() -> dict[str, set[int]]:
    rx = pickle.load(open(REACTIONS, "rb"))
    out: dict[str, set[int]] = defaultdict(set)
    for rid, r in rx.items():
        if r.get("transport"):
            continue
        for p in r.get("products", ()):
            out[p].add(int(rid))
    return out


def load_r2a() -> tuple[dict[int, set[str]], dict[int, tuple[str, str]]]:
    accs: dict[int, set[str]] = defaultdict(set)
    meta: dict[int, tuple[str, str]] = {}
    for line in open(R2A):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0].isdigit() and p[3]:
            accs[int(p[0])].add(p[3])
            meta[int(p[0])] = (p[1], p[2])  # direction, master_id
    return accs, meta


def resolve_missing(keys: list[str], cache: dict, limit: int | None = None) -> int:
    """Completa el cache UniProt -> Pfam por REST. Devuelve cuántas resolvió."""
    todo = keys[:limit] if limit else keys
    ok = 0
    for i, acc in enumerate(todo, 1):
        try:
            with urllib.request.urlopen(
                    f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=25) as fh:
                d = json.load(fh)
            cache[acc] = sorted({
                x["properties"][0]["value"]
                for x in d.get("uniProtKBCrossReferences", [])
                if x.get("database") == "Pfam" and x.get("properties")})
            ok += 1
        except Exception:
            cache.setdefault(acc, [])
        if i % 25 == 0:
            print(f"  resolve {i}/{len(todo)} ok={ok}", file=sys.stderr)
            CACHE.write_text(json.dumps(cache))
        time.sleep(0.035)
    CACHE.write_text(json.dumps(cache))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=1000)
    ap.add_argument("--resolve", action="store_true",
                    help="consulta UniProt para las accesiones sin dominios")
    ap.add_argument("--max-producers", type=int, default=5)
    ap.add_argument("--max-accs", type=int, default=10)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(WORKLIST), delimiter="\t"))

    def support(r: dict) -> int:
        v = r.get("support_400") or ""
        return int(v) if v.isdigit() else 0

    rows.sort(key=lambda r: (-support(r), r["genome"], r["gap"]))
    top = rows[:a.top]
    producers, accs, meta = load_producers(), *load_r2a()
    cache = json.load(open(CACHE))

    def pick_accs(chebi: str) -> list[str]:
        out: set[str] = set()
        for rid in sorted(producers.get(chebi, ()))[:a.max_producers]:
            out |= accs.get(rid, set())
        return sorted(out)[:a.max_accs]

    # 1) alcance sin red
    plan = []
    for r in top:
        acc_list = pick_accs(r["chebi"])
        rids = sorted(producers.get(r["chebi"], ()))
        direction, master = "", ""
        best = -1
        for rid in rids:
            d, m = meta.get(rid, ("", ""))
            if PRIORITY.get(d, 0) > best:
                best, direction, master = PRIORITY.get(d, 0), d, m
        plan.append({"row": r, "accs": acc_list, "rids": rids,
                     "direction": direction, "master": master})

    if a.resolve:
        missing = sorted({x for p in plan for x in p["accs"] if x not in cache})
        print(f"accesiones sin dominios en cache: {len(missing)}", file=sys.stderr)
        if missing:
            resolve_missing(missing, cache)

    # 2) dominios y estado
    out = []
    for p in plan:
        r = p["row"]
        sets = [set(cache.get(x, [])) for x in p["accs"]]
        union = sorted(set().union(*sets)) if sets else []
        # El dominio que define la enzima suele estar en TODOS los productores;
        # los accesorios (y los falsos positivos) aparecen solo en algunos.
        core: list[str] = []
        majority: list[str] = []
        if len(sets) >= 2:
            cnt = Counter(d for s in sets for d in s)
            core = sorted(d for d, c in cnt.items() if c == len(sets))
            thr = max(2, ceil(0.6 * len(sets)))
            majority = sorted(d for d, c in cnt.items() if c >= thr)
        else:
            core = majority = sorted(sets[0]) if sets else []
        if core:
            policy, allowed = "core", core
        elif majority:
            policy, allowed = "majority", majority
        else:
            policy, allowed = "any", union
        doms = allowed
        if not p["accs"]:
            status, reason = "no_accession", "sin accesión UniProt en productores"
        elif not doms:
            status, reason = "no_pfam", "productores sin dominio Pfam anotado"
        else:
            status, reason = "pending_hmm", "faltan los HMM de Pfam (658 MB)"
        if not p["rids"]:
            status, reason = "no_producer", "sin reacción productora en Rhea"
        out.append({
            "genome": r["genome"], "gap": r["gap"], "chebi": r["chebi"],
            "role": r.get("role", ""), "origin": r.get("origin", ""),
            "support_400": support(r), "n_producers": len(p["rids"]),
            "producers": ";".join(map(str, p["rids"][:5])),
            "producer_accs": ";".join(p["accs"]),
            "pfam_domains": ";".join(doms),
            "domain_policy": policy,
            "direction": p["direction"], "rhea_master": p["master"],
            "status": status, "reason": reason,
        })

    path = Path(str(OUT).format(top=a.top))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    c = Counter(o["status"] for o in out)
    dom = {d for o in out for d in o["pfam_domains"].split(";") if d}
    print(f"lote top {a.top}: {len(out)} filas -> {path.relative_to(ROOT)}")
    for k, v in c.most_common():
        print(f"  {k:14s} {v}")
    print(f"  dominios Pfam distintos: {len(dom)}")
    print(f"  dirección LR/RL: {sum(1 for o in out if o['direction'] in ('LR','RL'))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
