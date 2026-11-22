#!/usr/bin/env python3
"""Utilidades compartidas del módulo `transport_t1`.

Reglas del módulo:
  - Solo lectura de las tablas base de `data/metabolism/` y de
    `data/annotations_ab_initio/`.
  - Escritura exclusiva en `data/metabolism/transport_t1/`.
  - Sin commits.

Definición estricta de reacción de transporte (clave del módulo):
  una reacción Rhea se considera transporte *del metabolito* solo si
  `transport=True` **y** el ChEBI aparece a la vez en `substrates` y en
  `products` (misma molécula a ambos lados). El flag `transport` por sí solo
  también marca reacciones que comparten un ion (H+, Na+) y que no son
  transporte del metabolito (p. ej. acil-CoA sintetasas o decarboxilasas).
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
MET = ROOT / "data" / "metabolism"
OUT = MET / "transport_t1"
REF = OUT / "ref"
DIAMOND_DIR = OUT / "diamond"
ANN = ROOT / "data" / "annotations_ab_initio"
DB_DIR = ROOT / "data" / "db"
TAXDUMP = ROOT / "data" / "taxdump"

TRANSPORT_AUDIT = MET / "transport_audit.tsv"
GAP_TABLE = MET / "gap_interpretation.tsv"
RHEA2ACC = MET / "rhea2uniprot_sprot.tsv"
RHEA_PKL = MET / "rhea_reactions.pkl"
TRANSPORT_FASTA = DB_DIR / "transport.fasta"
TRANSPORT_DMND = DB_DIR / "transport.dmnd"

# Umbrales declarados (idénticos al intento previo de `webapp/transport_layer.py`,
# para que la comparación sea justa; ver README).
MIN_PIDENT = 35.0
MIN_QCOV = 0.70
MAX_EVALUE = 1e-5
# Búsqueda DIAMOND: se conservan hits débiles para la cola de revisión.
SEARCH_EVALUE = 1e-3
MAX_TARGET_SEQS = 0  # 0 = sin límite (todos los hits por consulta)

STATUSES = [
    "t1_bacterial",
    "t1_review",
    "transport_reaction_only",  # T2
    "no_transport",  # T3
]
STATUS_RANK = {s: i for i, s in enumerate(STATUSES)}
STATUS_LABEL = {
    "t1_bacterial": "T1",
    "t1_review": "T1_review (no importable)",
    "transport_reaction_only": "T2",
    "no_transport": "T3",
}

TSV_COLUMNS = [
    "genome",
    "gap",
    "chebi",
    "rhea",
    "acc",
    "pident",
    "qcov",
    "evalue",
    "status",
    "evidence",
]

HEADER_RE = re.compile(r"^>sp\|([^|]+)\|(\S+)\s+(.+)$")
OS_RE = re.compile(r"OS=(.*?) OX=(\d+)")
SUPERKINGDOM_IDS = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota"}


def clean(value: object, default: str = "") -> str:
    """Normaliza un campo para TSV/JSON: sin tabs, saltos ni extremos."""
    if value is None:
        return default
    text = str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return " ".join(text.split()).strip()


def read_tsv(path: Path) -> Iterator[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield row


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)
        fh.write("\n")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def count_rows(path: Path) -> int:
    with open(path, newline="", encoding="utf-8") as fh:
        return max(0, sum(1 for _ in fh) - 1)


def load_taxonomy_parents(nodes_path: Path = TAXDUMP / "nodes.dmp") -> Dict[int, int]:
    """Devuelve taxid -> parent taxid desde `nodes.dmp` de NCBI."""
    parent: Dict[int, int] = {}
    with open(nodes_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            head = line.split("\t|\t", 2)
            if len(head) < 2:
                continue
            try:
                parent[int(head[0])] = int(head[1])
            except ValueError:
                continue
    return parent


def superkingdom(taxid: int, parent: Dict[int, int]) -> str:
    """Clasifica un taxid en Bacteria/Archaea/Eukaryota/other.

    `other` es deliberadamente conservador: taxids ausentes o con linaje
    incompleto no se consideran bacterias y nunca producen T1.
    """
    seen = set()
    t = taxid
    while t and t not in seen and t != 1:
        if t in SUPERKINGDOM_IDS:
            return SUPERKINGDOM_IDS[t]
        seen.add(t)
        t = parent.get(t, 0)
    return "other"


def load_rhea(rhea_pkl: Path = RHEA_PKL) -> Tuple[dict, Dict[str, List[int]], Dict[str, List[int]]]:
    """Carga Rhea y devuelve (rx, estricto_por_chebi, un_lado_por_chebi).

    `estricto_por_chebi[chebi]`: reacciones de transporte con el ChEBI a ambos
    lados (transporte real del metabolito).
    `un_lado_por_chebi[chebi]`: reacciones con `transport=True` donde el ChEBI
    aparece solo en un lado (p. ej. acil-CoA sintetasas); no son T2.
    """
    import pickle

    with open(rhea_pkl, "rb") as fh:
        rx = pickle.load(fh)
    strict: Dict[str, set] = defaultdict(set)
    one_side: Dict[str, set] = defaultdict(set)
    for rid, data in rx.items():
        if not data.get("transport"):
            continue
        rid = int(rid)
        subs = data.get("substrates", set()) or set()
        prods = data.get("products", set()) or set()
        for chebi in subs & prods:
            strict[chebi].add(rid)
        for chebi in subs ^ prods:
            one_side[chebi].add(rid)
    return (
        rx,
        {k: sorted(v) for k, v in strict.items()},
        {k: sorted(v) for k, v in one_side.items()},
    )


def load_rhea2acc(
    path: Path = RHEA2ACC, transport_rheas: set | None = None
) -> Dict[int, set]:
    """Carga Rhea -> accs Swiss-Prot (opcionalmente solo Rhea de transporte)."""
    out: Dict[int, set] = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] == "RHEA_ID":
                continue
            try:
                rid = int(f[0])
            except ValueError:
                continue
            if transport_rheas is not None and rid not in transport_rheas:
                continue
            out[rid].add(f[3])
    return dict(out)


def parse_transport_fasta(path: Path = TRANSPORT_FASTA) -> List[dict]:
    """Parsea `transport.fasta`; aborta si una cabecera no es `>sp|acc|entry`."""
    records = []
    seen = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line.rstrip("\n"))
            if not m:
                raise ValueError(f"cabecera no parseable en {path}: {line[:120]!r}")
            acc, entry, desc = m.group(1), m.group(2), m.group(3)
            if acc in seen:
                raise ValueError(f"acc duplicada en {path}: {acc}")
            seen.add(acc)
            om = OS_RE.search(desc)
            if not om:
                raise ValueError(f"cabecera sin OS/OX: {acc}")
            records.append(
                {
                    "acc": acc,
                    "entry": entry,
                    "organism": clean(om.group(1)),
                    "taxid": int(om.group(2)),
                    "description": clean(desc),
                }
            )
    return records


def fasta_ids(path: Path) -> List[str]:
    ids = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def fmt_pident(value: float) -> str:
    return f"{value:.1f}"


def fmt_qcov(value: float) -> str:
    return f"{value:.3f}"


def fmt_evalue(value: float) -> str:
    return f"{value:.2e}"


def load_ref() -> Tuple[Dict[str, dict], Dict[str, List[int]], Dict[int, set], Dict[str, set]]:
    """Carga las tablas de referencia generadas por `01_build_reference.py`."""
    db: Dict[str, dict] = {}
    for row in read_tsv(REF / "transport_db.tsv"):
        db[row["acc"]] = row
    chebi2rhea: Dict[str, List[int]] = defaultdict(list)
    for row in read_tsv(REF / "chebi_transport_rhea.tsv"):
        chebi2rhea[row["chebi"]].append(int(row["rhea"]))
    rhea2acc: Dict[int, set] = defaultdict(set)
    acc2rhea: Dict[str, set] = defaultdict(set)
    for row in read_tsv(REF / "rhea2acc_transport.tsv"):
        rid, acc = int(row["rhea"]), row["acc"]
        rhea2acc[rid].add(acc)
        acc2rhea[acc].add(rid)
    return (
        db,
        {k: sorted(v) for k, v in chebi2rhea.items()},
        dict(rhea2acc),
        dict(acc2rhea),
    )
