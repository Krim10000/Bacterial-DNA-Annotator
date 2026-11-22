#!/usr/bin/env python3
"""Cierre de brechas con evidencia dominio → EC → Rhea, etiquetado.

Para cada cluster con un dominio Pfam mapeado a EC:
  - EC → reacciones Rhea;
  - para cada genoma miembro, si una brecha refinada tiene esa reacción entre
    sus productoras, se emite un `gap_closure_candidate`.

La confianza es `high` si la reacción es direccional (LR/RL) y el cluster es
conservado; `medium` si la reacción es reversible (UN/BI) o el cluster es único.
"""
from __future__ import annotations

import argparse
import csv
import pickle
from collections import defaultdict
from pathlib import Path

from .labels import SOURCES, label_row

# Mapa curado dominio Pfam -> EC.
DOMAIN_EC = {
    "RibD_C": "1.1.1.193", "GTP_cyclohydro2": "3.5.4.25", "GTP_cyclohydroI": "3.5.4.16",
    "Glyoxalase": "4.4.1.5", "Glyoxalase_2": "4.4.1.5", "Glyoxalase_3": "4.4.1.5",
    "Glyoxalase_4": "4.4.1.5", "DHFR_1": "1.5.1.3", "URO-D": "4.1.1.37",
    "Coprogen_oxidas": "1.3.3.3", "FolB": "4.1.2.25", "FolK": "2.7.6.3",
    "FolP": "2.5.1.15", "FolE": "3.5.4.16", "ThiE": "2.5.1.3", "TPK": "2.7.6.2",
    "Transketolase": "2.2.1.1", "MobB": "2.7.7.75",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--tier2-eval", required=True, type=Path)
    ap.add_argument("--tier2-members", required=True, type=Path)
    ap.add_argument("--tier1-panel", required=True, type=Path)
    ap.add_argument("--tier1-members", required=True, type=Path)
    ap.add_argument("--refined", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cache", type=Path, default=Path("cache/rhea_reactions.pkl"))
    a = ap.parse_args()

    ec2r: dict[str, set[str]] = defaultdict(set)
    for line in open(a.data_root / "reactions" / "tsv" / "rhea2ec.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0] != "RHEA_ID":
            ec2r[p[3]].add(f"RHEA:{p[0]}")
    rx_type = {}
    if a.cache.exists():
        raw = pickle.load(open(a.cache, "rb"))
        rx_type = {f"RHEA:{k}": v.get("type", "?") for k, v in raw.items()}

    domains: dict[str, set[str]] = defaultdict(set)
    for r in csv.DictReader(open(a.tier2_eval), delimiter="\t"):
        if r.get("pfam"):
            domains[r["rep"].split("|")[0]].add(r["pfam"])
    for line in open(a.tier1_panel):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) >= 23:
            domains[p[0]].add(p[3])

    members: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for f in (a.tier1_members, a.tier2_members):
        for r in csv.DictReader(open(f), delimiter="\t"):
            members[r["cluster"]].append((r["genome"], r["locus"]))
    ng = {c: len({g for g, _ in v}) for c, v in members.items()}

    gaps: dict[str, list[tuple[str, set[str]]]] = defaultdict(list)
    for r in csv.DictReader(open(a.refined), delimiter="\t"):
        gaps[r["genome"]].append((r["gap"], set(r["producers"].split(";"))))

    rows: list[dict] = []
    for cid, ds in domains.items():
        for dom in ds:
            ec = DOMAIN_EC.get(dom)
            if not ec:
                continue
            for rhea in ec2r.get(ec, set()):
                for genome, locus in members.get(cid, []):
                    for gap, prod in gaps.get(genome, []):
                        if rhea not in prod:
                            continue
                        directional = rx_type.get(rhea) in ("LR", "RL")
                        conf = "high" if (directional and ng.get(cid, 1) >= 2) else "medium"
                        row = {
                            "cluster": cid, "domain": dom, "ec": ec, "rhea": rhea,
                            "genome": genome, "locus": locus, "gap": gap,
                            "n_genomes": ng.get(cid, 1),
                        }
                        label_row(
                            row, label="gap_closure_candidate",
                            source="Pfam-A → EC → Rhea + " + SOURCES["conservation"],
                            evidence=(f"pfam={dom}; EC={ec}; {rhea}; genoma={genome}; "
                                      f"ng={ng.get(cid, 1)}; direccion={rx_type.get(rhea, '?')}"),
                            confidence=conf,
                            interpretation=("Proteína conservada con el dominio de la enzima que "
                                            "produce el metabolito faltante; verificar dirección "
                                            "y especificidad antes de confirmar."),
                        )
                        rows.append(row)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["cluster"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"cierres etiquetados: {len(rows)}")
    for r in rows[:10]:
        print(f"  {r['cluster']} {r['domain']} EC={r['ec']} {r['rhea']} {r['genome']} "
              f"gap={r['gap']} ng={r['n_genomes']} conf={r['label_confidence']}")
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
