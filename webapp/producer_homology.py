"""¿La enzima productora está presente pero fuera de la red mapeada?

Busca homólogos de las accs de los productores Rhea de cada brecha.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

EUK = ("Homo sapiens", "Mus musculus", "Rattus", "Arabidopsis", "Saccharomyces")


def load_r2a(data_root: Path) -> dict[str, set[str]]:
    r2a: dict[str, set[str]] = defaultdict(set)
    path = Path(data_root) / "metabolism" / "rhea2uniprot_sprot.tsv"
    if not path.exists():
        return r2a
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0] != "RHEA_ID":
            r2a["RHEA:" + p[0]].add(p[3])
    return r2a


def load_rhea_ec(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not Path(path).exists():
        return out
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0] != "RHEA_ID":
            out["RHEA:" + p[0]].add(p[3])
    return out


def load_db_accs(path: Path) -> set[str]:
    out = set()
    if not Path(path).exists():
        return out
    for line in open(path):
        if line.startswith('>') and '|' in line:
            out.add(line.split('|')[1])
    return out


def annotate(gaps: list[dict], hits: dict[str, dict], r2a: dict[str, set[str]],
             db_accs: set[str] | None = None,
             rhea_ec: dict[str, set[str]] | None = None) -> list[dict]:
    for g in gaps:
        prod = g.get("producers", [])
        if isinstance(prod, str):
            prod = [x for x in prod.split(";") if x]
        accs = set().union(*(r2a.get(p, set()) for p in prod)) if prod else set()
        checked = accs & db_accs if db_accs else accs
        if not checked:
            g["homolog_status"] = "not_checked"
            continue
        cand = [hits[a] for a in checked if a in hits]
        if not cand:
            g["homolog_status"] = "no_homolog"
            continue
        b = sorted(cand, key=lambda h: -h["bits"])[0]
        euk = any(x in b.get("desc", "") for x in EUK)
        g["homolog_acc"] = b["acc"]
        g["homolog_pident"] = b["pident"]
        g["homolog_qcov"] = round(b["qcov"], 3)
        g["homolog_description"] = b.get("desc", "")
        cand_rhea = [p for p in prod if b["acc"] in r2a.get(p, set())]
        g["homolog_rhea"] = cand_rhea[0] if cand_rhea else ""
        if rhea_ec and g["homolog_rhea"]:
            g["homolog_ec"] = ";".join(sorted(rhea_ec.get(g["homolog_rhea"], set())))
        if (not euk) and b["pident"] >= 35 and b["qcov"] >= 0.7:
            g["homolog_status"] = "homolog_strong"
            g["origin"] = "homology_present"
            g["status"] = "consensus"
            g["status_reason"] = f"homólogo {b['acc']} pident={b['pident']} qcov={g['homolog_qcov']}; reacción no mapeada"
            g["review_status"] = "homology_mapping"
            g["review_reason"] = "enzima presente; corregir mapeo Rhea/EC"
            g["mapping_suggestion"] = "suggested"
        elif b["pident"] >= 30 and b["qcov"] >= 0.6:
            g["homolog_status"] = "homolog_weak"
            g["review_status"] = "homology_weak"
            g["review_reason"] = "homólogo débil; revisar antes de mapear"
            g["mapping_suggestion"] = "review"
        else:
            g["homolog_status"] = "homolog_low"
    return gaps
