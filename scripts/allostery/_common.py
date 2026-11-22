#!/usr/bin/env python3
"""Utilidades compartidas del módulo de alosteria (F0/F1).

- Set decisivo: represores de RegulonDB con efector en UNA sola conformación.
- Secuencias de proteína desde el GenBank de E. coli K-12.
- Familia real por dominio Pfam (hmmsearch) — no por nombre de gen.
"""
from __future__ import annotations

import collections
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RD = ROOT / "data" / "regulondb"
GBFF = ROOT / "data" / "test_ecoli.gbff"
WORK = ROOT / "data" / "allostery" / "work"

# Dominio Pfam -> familia de TF. Se completa con lo que devuelve hmmsearch.
FAMILY_MAP = {
    "LacI": "LacI/GalR",
    "LacI_substrate": "LacI/GalR",
    "HTH_1": "LysR",
    "LysR_substrate": "LysR",
    "TetR_N": "TetR",
    "GntR": "GntR",
    "UTRA": "GntR",
    "HTH_DeoR": "DeoR",
    "DeoR": "DeoR",
    "HTH_MARR": "MarR",
    "MarR": "MarR",
    "HTH_AraC": "AraC",
    "AraC_binding": "AraC",
    "HTH_CRP": "Crp/FNR",
    "HTH_3": "HTH_3",
    "HTH_XRE": "XRE",
    "HTH_5": "ArsR/SmtB",
    "ArsR": "ArsR/SmtB",
    "MerR": "MerR",
    "MerR_DNA": "MerR",
    "HTH_9": "MerR",
    "HTH_17": "MerR",
    "DtxR": "DtxR",
    "Fe_dep_repress": "DtxR",
    "Fur": "Fur",
    "HTH_24": "Fur",
    "HTH_29": "Fur",
    "HTH_26": "PadR",
    "PadR": "PadR",
    "Rrf2": "Rrf2",
    "AsnC_trans_reg": "AsnC",
    "IclR": "IclR",
    "RpiR": "RpiR",
    "SorC": "SorC",
    "HTH_8": "PadR",
    "NrdR": "NrdR",
    "MetJ": "MetJ",
    "HxlR": "HxlR",
    "NikR": "NikR",
    "CopY": "CopY",
    "CsoR": "CsoR",
    "HTH_31": "MarR",
    "HTH_34": "MarR",
    "HTH_42": "MarR",
}

# Dominios que no son de unión a ADN: se excluyen al elegir la familia.
NOT_DBD = re.compile(r"(substrate|Peripla|PAS|GAF|SBP_bac|HAMP|Response_reg|"
                     r"Trans_reg|AAA|ATPase|kinase|peptidase|TM|DUF)", re.I)


def load_proteins():
    txt = GBFF.read_text(errors="ignore")
    prot = {}
    for chunk in txt.split("     CDS     ")[1:]:
        g = re.search(r'/gene="([^"]+)"', chunk)
        t = re.search(r'/translation="([^"]+)"', chunk)
        if g and t:
            prot[g.group(1).lower()] = "".join(t.group(1).split())
    return prot


def load_roles():
    roles = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(RD / "NetworkRegulatorGene.tsv"), delimiter="\t"):
        roles[r["2)regulatorName"]][r["6)function"]] += 1
    out = {}
    for t, c in roles.items():
        pos, neg = c.get("+", 0), c.get("-", 0)
        out[t] = ("dual" if pos and neg else
                  "activador" if pos else "represor" if neg else "sin_aristas")
    return out


def load_decisive():
    """Represores (no duales) con efector en una sola conformación."""
    prot = load_proteins()
    roles = load_roles()
    out = []
    for r in csv.DictReader(open(RD / "TFSet.tsv"), delimiter="\t"):
        tf = r["2)name"]
        ea = r["10)tfActiveConformationsEffector"].strip()
        ei = r["11)tfInactiveConformationsEffector"].strip()
        if not (ea or ei) or roles.get(tf) != "represor" or (ea and ei):
            continue
        out.append({
            "tf": tf,
            "mode": "correpresor" if ea else "inductor",
            "effector": (ea or ei).split("|")[0],
            "conf": r["17)confidenceLevel"],
            "pmid": r["18)tfConformationPMID"].strip(),
            "seq": prot.get(tf.lower(), ""),
        })
    return out


def load_families(domtbl: Path | None = None):
    """Familia por el mejor dominio Pfam de unión a ADN (hmmsearch --domtblout)."""
    domtbl = domtbl or (WORK / "decisive.domtbl")
    best = collections.defaultdict(list)
    for line in domtbl.read_text(errors="ignore").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        # target(seq) acc query(HMM) acc ... evalue score ... hmmfrom hmmto ...
        seq_name, hmm_name = f[0], f[3]
        evalue = float(f[12])
        score = float(f[13])
        if evalue > 1e-5:
            continue
        best[seq_name].append((score, hmm_name, evalue))
    out = {}
    for tf, hits in best.items():
        dbd = [(s, n) for s, n, _ in hits if n in FAMILY_MAP]
        if not dbd:
            dbd = [(s, n) for s, n, _ in hits if not NOT_DBD.search(n)]
        if not dbd:
            out[tf] = "sin_dominio"
            continue
        score, name = max(dbd)
        out[tf] = FAMILY_MAP.get(name, name)
    return out
