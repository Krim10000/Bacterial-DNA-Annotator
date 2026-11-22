#!/usr/bin/env python3
"""03_resolve_domains.py — dominios Pfam esperados por brecha.

Para las accesiones productoras **no eucariotas** del lote:
  * resuelve sus dominios Pfam (accesiones PFxxxxx) por lotes contra UniProt
    REST (`fields=accession,xref_pfam`), con caché local reanudable en
    `work/pfam_acc_cache.json`;
  * por brecha aplica la política del proyecto: `core` (en todos los
    productores con Pfam) > `majority` (>=60%) > `any` (unión). Elegir el
    dominio de mejor e-value entre todos producía falsos positivos;
  * escribe `work/gap_domains.tsv`.

Si UniProt no responde para una accesión, queda sin dominios: la brecha se
marcará `not_checked` en el eje de dominio, nunca `no_domain`.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import WORK, read_tsv, write_json, write_tsv

CACHE = WORK.parent / "work" / "pfam_acc_cache.json"
BATCH = 80
TRIES = 4


def fetch_batch(accs: list[str]) -> dict[str, list[str]]:
    query = " OR ".join(f"accession:{a}" for a in accs)
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({"query": query,
                                     "fields": "accession,xref_pfam",
                                     "format": "tsv", "size": 500}))
    last = None
    for attempt in range(TRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                text = fh.read().decode()
            out: dict[str, list[str]] = {}
            for line in text.splitlines()[1:]:
                p = line.split("\t")
                if len(p) < 2:
                    continue
                pfams = [x.strip() for x in p[1].replace(",", ";").split(";")
                         if x.strip().startswith("PF")]
                out[p[0]] = pfams
            return out
        except Exception as exc:  # red intermitente
            last = exc
            time.sleep(1.5 * (attempt + 1))
    print(f"  fallo lote ({len(accs)} accs): {last}", file=sys.stderr)
    return {}


def main() -> int:
    rows = read_tsv(WORK / "selected_gaps.tsv")
    acc_info = {r["acc"]: r for r in read_tsv(WORK / "acc_info.tsv")}
    prok = {a for a, r in acc_info.items() if r["euk"] == "0"}
    all_accs = sorted({a for r in rows for a in r["producer_accs"].split(";") if a})
    todo = [a for a in all_accs if a in prok]
    print(f"accesiones del lote: {len(all_accs)}; no eucariotas: {len(todo)}")

    cache: dict[str, list[str]] = {}
    if CACHE.exists():
        cache = json.loads(CACHE.read_text())
    missing = [a for a in todo if a not in cache]
    print(f"sin consultar: {len(missing)}")
    for i in range(0, len(missing), BATCH):
        chunk = missing[i:i + BATCH]
        got = fetch_batch(chunk)
        for a in chunk:
            cache[a] = got.get(a, [])
        CACHE.write_text(json.dumps(cache))
        print(f"  {min(i + BATCH, len(missing))}/{len(missing)} "
              f"(resueltas {sum(1 for a in chunk if cache[a])}/{len(chunk)})",
              file=sys.stderr)
        time.sleep(0.2)
    unresolved = [a for a in todo if a not in cache]
    no_pfam = [a for a in todo if a in cache and not cache[a]]

    out_rows = []
    for r in rows:
        accs = [a for a in r["producer_accs"].split(";") if a and a in prok]
        sets = [set(cache.get(a, [])) for a in accs if cache.get(a)]
        n_with = len(sets)
        union = set().union(*sets) if sets else set()
        core: list[str] = []
        majority: list[str] = []
        if n_with >= 2:
            cnt = Counter(d for s in sets for d in s)
            core = sorted(d for d, c in cnt.items() if c == n_with)
            thr = max(2, (60 * n_with + 99) // 100)
            majority = sorted(d for d, c in cnt.items() if c >= thr)
        elif n_with == 1:
            core = majority = sorted(sets[0])
        if core:
            policy, domains = "core", core
        elif majority:
            policy, domains = "majority", majority
        else:
            policy, domains = "any", sorted(union)
        if not accs:
            status, reason = "not_checked", "sin accesiones no eucariotas"
        elif not n_with:
            status, reason = "not_checked", "ningún productor con Pfam resuelto"
        elif not domains:
            status, reason = "not_checked", "productores sin dominio Pfam anotado"
        else:
            status, reason = "ok", f"pfam en {n_with}/{len(accs)} accs"
        out_rows.append({
            "genome": r["genome"], "gap": r["gap"],
            "n_producers": r["n_producers"], "n_accs": r["n_accs"],
            "n_accs_prok": len(accs), "n_accs_with_pfam": n_with,
            "domains_core": ";".join(core),
            "domains_majority": ";".join(majority),
            "domains": ";".join(domains),
            "domain_policy": policy, "status": status, "reason": reason,
        })

    write_tsv(WORK / "gap_domains.tsv", out_rows,
              ["genome", "gap", "n_producers", "n_accs", "n_accs_prok",
               "n_accs_with_pfam", "domains_core", "domains_majority",
               "domains", "domain_policy", "status", "reason"])
    stats = {
        "accs_total": len(all_accs),
        "accs_prok": len(todo),
        "accs_with_pfam": len(todo) - len(no_pfam) - len(unresolved),
        "accs_without_pfam": len(no_pfam),
        "accs_unresolved": len(unresolved),
        "unique_domains": len({d for r in out_rows for d in r["domains"].split(";") if d}),
        "gap_domain_status": dict(Counter(r["status"] for r in out_rows)),
        "policies": dict(Counter(r["domain_policy"] for r in out_rows)),
    }
    write_json(WORK / "domains.json", stats)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
