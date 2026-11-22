"""Red metabólica dinámica para un genoma subido.

- Si la especie es conocida (GCF en el proyecto), reutiliza su red precalculada.
- Si no, anota las proteínas con DIAMOND vs Swiss-Prot, las mapea a Rhea y
  construye la red.
- Detecta brechas y las prioriza por cuántos de los 400 genomas de referencia
  tienen una enzima que las produciría (reutilización inteligente).
"""
from __future__ import annotations

import csv
import json
import pickle
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SOURCES = {
    "CHEBI:15377", "CHEBI:15378", "CHEBI:15379", "CHEBI:16526", "CHEBI:30616",
    "CHEBI:456216", "CHEBI:456215", "CHEBI:57945", "CHEBI:57783", "CHEBI:57287",
    "CHEBI:58349", "CHEBI:57540",
}


def load_reactions(cache: Path) -> dict:
    return pickle.load(open(cache, "rb"))


def load_r2a(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3]:
                out[f"RHEA:{p[0]}"].add(p[3])
    return out


def load_names(path: Path) -> dict[str, str]:
    out = {}
    if not path.exists():
        return out
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            out[p[0].strip()] = p[1].strip()
    return out


def load_reference_counts(path: Path) -> dict[str, int]:
    try:
        return json.load(open(path))
    except Exception:
        return {}


def _diamond(proteins: list[dict], db: Path, threads: int) -> set[str]:
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
             "length", "qlen", "evalue", "bitscore",
             "--evalue", "1e-5", "--max-target-seqs", "1",
             "--threads", str(threads), "--quiet"],
            check=True,
        )
        accs = set()
        with open(out) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and "|" in p[1]:
                    accs.add(p[1].split("|")[1])
        return accs


def build_network(proteins: list[dict], data_root: Path, ref_gcf: str | None = None,
                  threads: int = 4, panel_hmm: Path | None = None,
                  transport_db: Path | None = None,
                  producer_db: Path | None = None) -> dict:
    r2a = load_r2a(data_root / "metabolism" / "rhea2uniprot_sprot.tsv")
    acc2r: dict[str, set[str]] = defaultdict(set)
    for r, accs in r2a.items():
        for acc in accs:
            acc2r[acc].add(r)
    reactions = load_reactions(data_root / "metabolism" / "rhea_reactions.pkl")
    ref_counts = load_reference_counts(data_root / "metabolism" / "reference_reactions.json")
    names = load_names(data_root / "reactions" / "tsv" / "chebiId_name.tsv")

    ref_net = data_root / "metabolism" / f"{ref_gcf}_network.json" if ref_gcf else None
    if ref_net and ref_net.exists():
        d = json.load(open(ref_net))
        net = {n["id"] for n in d.get("nodes", []) if n.get("type") == "reaction"}
        source = f"red precalculada de {ref_gcf}"
    else:
        db = data_root / "db" / "swissprot.dmnd"
        if not db.exists():
            raise FileNotFoundError(
                "Falta Swiss-Prot. Descárgalo con: python -m webapp.setup_data --data-root data"
            )
        accs = _diamond(proteins, db, threads)
        net = {r for acc in accs for r in acc2r.get(acc, ())}
        source = "DIAMOND vs Swiss-Prot → Rhea"

    gaps = detect_gaps(net, reactions, r2a, ref_counts, names)
    try:
        from gap_logic import interpret_gaps, annotate_pathways
        gaps = interpret_gaps(gaps)
        gaps = annotate_pathways(gaps, data_root)
    except Exception:
        pass
    if transport_db and Path(transport_db).exists():
        try:
            from transport_layer import load_transport_by_chebi, search, annotate_gaps
            rx_cache = data_root / "metabolism" / "rhea_reactions.pkl"
            tbc = load_transport_by_chebi(rx_cache)
            hits = search(proteins, Path(transport_db), threads=4)
            gaps = annotate_gaps(gaps, hits, tbc, r2a)
        except Exception:
            pass
    if producer_db and Path(producer_db).exists():
        try:
            from producer_homology import annotate as annotate_producers, load_db_accs, load_rhea_ec
            from transport_layer import search as search_homologs
            hits = search_homologs(proteins, Path(producer_db), threads=4)
            db_accs = load_db_accs(data_root / "db" / "producer_family.fasta")
            gaps = annotate_producers(gaps, hits, r2a, db_accs, load_rhea_ec(data_root / 'reactions' / 'tsv' / 'rhea2ec.tsv'))
        except Exception:
            pass
    try:
        from gap_logic import classify_review
        gaps = classify_review(gaps, data_root)
    except Exception:
        pass
    hmm_hits = []
    if panel_hmm and Path(panel_hmm).exists():
        try:
            from hmm_panel import confirm
            hmm_hits = confirm(gaps, proteins, Path(panel_hmm))
        except Exception:
            hmm_hits = []
    return {
        "source": source,
        "reactions": len(net),
        "proteins": len(proteins),
        "gaps": gaps,
        "n_gaps": len(gaps),
        "net": sorted(net),
        "gap_producers": sorted({r for g in gaps for r in g["producers"].split(";") if r}),
        "hmm_hits": hmm_hits,
    }


def _acc_to_rhea(acc: str, r2a: dict[str, set[str]]) -> set[str]:
    return {r for r, accs in r2a.items() if acc in accs}


def detect_gaps(net: set[str], reactions: dict, r2a: dict[str, set[str]],
                ref_counts: dict[str, int], names: dict[str, str],
                max_gaps: int = 500) -> list[dict]:
    producers: dict[str, set[str]] = defaultdict(set)
    consumers: dict[str, list[str]] = defaultdict(list)
    for rid, r in reactions.items():
        if r.get("transport"):
            continue
        key = f"RHEA:{rid}"
        for p in r.get("products", set()):
            producers[p].add(key)
        for s in r.get("substrates", set()):
            consumers[s].append(key)

    metabolites: set[str] = set()
    for rid in net:
        try:
            r = reactions[int(rid.split(":")[1])]
        except (KeyError, ValueError, IndexError):
            continue
        metabolites |= r.get("products", set()) | r.get("substrates", set())

    gaps = []
    for m in metabolites:
        if m in SOURCES:
            continue
        prods = producers.get(m, set())
        if prods & net:
            continue
        enz = [r for r in prods if r2a.get(r)]
        if not enz:
            continue
        cons = [r for r in consumers.get(m, []) if r in net]
        if not cons:
            continue
        support = max((ref_counts.get(r, 0) for r in enz), default=0)
        gaps.append({
            "chebi": m, "gap": names.get(m, m),
            "n_producers": len(enz), "support": support,
            "producers": ";".join(sorted(enz)[:4]),
            "consumers": ";".join(sorted(cons)[:3]),
            "label": "gap_supported" if support > 0 else "gap_novel",
        })
    gaps.sort(key=lambda g: (g["support"], g["n_producers"]), reverse=True)
    return gaps[:max_gaps]
