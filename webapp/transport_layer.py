"""Capa de transporte para genomas subidos.

Busca homólogos de transportadores en la base Rhea/UniProt y los asocia a la
reacción de transporte exacta del metabolito. Solo marca importación con
homólogo no eucariota, pident >=35% y qcov >=0,7.
"""
from __future__ import annotations

import pickle
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

EUK = re.compile(r"Homo sapiens|Mus musculus|Rattus|Arabidopsis|Saccharomyces|"
                 r"Drosophila|Caenorhabditis|Eptatretus|Danio|Xenopus|Bos taurus|"
                 r"Sus scrofa|Gallus|eukaryot|human|mouse|rat\b", re.I)


def load_transport_by_chebi(cache: Path) -> dict[str, set[str]]:
    rx = pickle.load(open(cache, "rb"))
    out: dict[str, set[str]] = defaultdict(set)
    for rid, d in rx.items():
        if d.get("transport"):
            for c in d.get("substrates", set()) | d.get("products", set()):
                out[c].add(f"RHEA:{rid}")
    return out


def search(proteins: list[dict], db: Path, threads: int = 2) -> dict[str, dict]:
    if not db.exists() or not proteins:
        return {}
    with tempfile.TemporaryDirectory() as td:
        faa = Path(td) / "proteins.faa"
        with open(faa, "w") as fh:
            for p in proteins:
                fh.write(f">{p['id']}\n")
                aa = p["aa"]
                for i in range(0, len(aa), 60):
                    fh.write(aa[i:i + 60] + "\n")
        out = Path(td) / "hits.tsv"
        subprocess.run(
            ["diamond", "blastp", "--query", str(faa), "--db", str(db),
             "--out", str(out), "--outfmt", "6", "qseqid", "sseqid", "pident",
             "length", "qlen", "evalue", "bitscore", "stitle",
             "--max-target-seqs", "5", "--evalue", "1e-5",
             "--threads", str(threads), "--quiet"],
            check=True,
        )
        best: dict[str, dict] = {}
        for line in open(out):
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            try:
                acc = f[1].split("|")[1]
            except Exception:
                continue
            qcov = float(f[3]) / float(f[4]) if float(f[4]) else 0.0
            h = {"query": f[0], "acc": acc, "pident": float(f[2]),
                 "qcov": qcov, "evalue": float(f[5]), "bits": float(f[6]),
                 "desc": (f[7] if len(f) > 7 else "")[:100]}
            if acc not in best or h["bits"] > best[acc]["bits"]:
                best[acc] = h
        return best


def annotate_gaps(gaps: list[dict], hits: dict[str, dict],
                  transport_by_chebi: dict[str, set[str]],
                  r2a: dict[str, set[str]]) -> list[dict]:
    for g in gaps:
        chebi = g.get("chebi", "")
        trs = transport_by_chebi.get(chebi, set())
        g["transport_rhea"] = ";".join(sorted(trs)[:2])
        if not trs:
            g["transport_status"] = "no_transport"
            g["transport_reason"] = "sin reacción de transporte para el metabolito"
            continue
        accs = set().union(*(r2a.get(t, set()) for t in trs)) if trs else set()
        found = []
        for acc in accs:
            h = hits.get(acc)
            if not h:
                continue
            found.append(h)
        good = [h for h in found if not EUK.search(h.get("desc", ""))
                and h["pident"] >= 35 and h["qcov"] >= 0.7]
        if good:
            b = sorted(good, key=lambda h: h["evalue"])[0]
            g["transport_status"] = "t1_bacterial"
            g["transporter_acc"] = b["acc"]
            g["transport_pident"] = b["pident"]
            g["transport_qcov"] = round(b["qcov"], 3)
            g["transport_evalue"] = b["evalue"]
            g["transport_reason"] = f"transportador {b['acc']} pident={b['pident']} qcov={g['transport_qcov']}"
        elif found:
            b = sorted(found, key=lambda h: h["evalue"])[0]
            g["transport_status"] = "t1_review"
            g["transporter_acc"] = b["acc"]
            g["transport_pident"] = b["pident"]
            g["transport_qcov"] = round(b["qcov"], 3)
            g["transport_evalue"] = b["evalue"]
            g["transport_reason"] = "homólogo eucariota o identidad/cobertura baja; revisar"
        else:
            g["transport_status"] = "transport_reaction_only"
            g["transport_reason"] = "reacción de transporte en Rhea; sin homólogo en el genoma"
    return gaps
