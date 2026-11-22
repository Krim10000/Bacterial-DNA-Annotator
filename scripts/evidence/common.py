#!/usr/bin/env python3
"""Constantes y utilidades comunes del módulo de evidencia por brecha.

Solo lectura sobre las tablas base; todas las escrituras van a
`data/metabolism/evidence/` (salida + `work/` intermedio) o a este directorio
de scripts.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data" / "metabolism"
REVIEW_QUEUE = Path(os.environ.get("EVIDENCE_QUEUE",
                                   str(META / "REVIEW_QUEUE.tsv")))
REAL_GAPS = Path(os.environ.get("EVIDENCE_REAL",
                                str(META / "real_gaps.tsv")))
REFINED = ROOT / "out_fl" / "refined" / "refined_gaps.tsv"
R2A = META / "rhea2uniprot_sprot.tsv"
SPROT = ROOT / "data" / "db" / "uniprot_sprot.fasta"
ANNOT = ROOT / "data" / "annotations"          # <genome>.faa.gz (NCBI)
ABINITIO = ROOT / "data" / "annotations_ab_initio"  # <genome>/proteins.faa
TAXDUMP = ROOT / "data" / "taxdump"

EVID = META / "evidence"

# --- configuración por lote -------------------------------------------------
# El lote 01 usa las rutas históricas (work/, gap_evidence.tsv, summary.json).
# Los lotes siguientes se aíslan en work_batch_<NN>/ y salidas con sufijo.
BATCH = os.environ.get("EVIDENCE_BATCH", "01")
WORK = EVID / ("work" if BATCH == "01" else f"work_batch_{BATCH}")
WORK.mkdir(parents=True, exist_ok=True)
REVIEW = Path(os.environ.get(
    "EVIDENCE_REVIEW", str(META / "review" / "gap_review.tsv")))
EXCLUDE = (Path(os.environ["EVIDENCE_EXCLUDE"])
           if os.environ.get("EVIDENCE_EXCLUDE") else None)
N_BATCH = int(os.environ.get("EVIDENCE_N", "500"))
_SUFFIX = "" if BATCH == "01" else f"_batch_{BATCH}"
TABLE_OUT = EVID / f"gap_evidence{_SUFFIX}.tsv"
SUMMARY_OUT = EVID / f"summary{_SUFFIX}.json"
README_OUT = EVID / f"README{_SUFFIX}.md"

KEY = ["genome", "gap"]
# Umbrales declarados (los mismos del resto del proyecto, sin ajustar al lote):
#   fuerte  : pident >= 35 %, qcov >= 0.70, evalue <= 1e-5, productor no eucariota
#   débil   : pident >= 30 %, qcov >= 0.60, evalue <= 1e-5, productor no eucariota
MIN_PIDENT_STRONG, MIN_QCOV_STRONG = 35.0, 0.70
MIN_PIDENT_WEAK, MIN_QCOV_WEAK = 30.0, 0.60
MAX_EVALUE = 1e-5
# Dominios genéricos: no sostienen por sí solos "dominio presente".
GENERIC = {
    "FAD_binding_4", "FAD_binding_2", "NAD_binding_1", "NAD_binding_2",
    "Rossmann-like", "Aldedh", "Aminotran_1_2", "Aminotran_3", "Aminotran_4",
    "Aminotran_5", "Metallophos", "Radical_SAM", "GAF_2", "AAA_26",
}
# Acciones del gap-reviewer que ya cubre `reannotation-builder` (77 filas).
ACTIONABLE = {"annotation_gap", "domain_present"}


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"sin filas para {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=False) + "\n")


def clean(value: object, default: str = "") -> str:
    """Normaliza a una celda TSV sin tabs ni saltos."""
    if value is None:
        return default
    text = str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return " ".join(text.split()).strip() or default


def fasta_ids(path: Path):
    """Itera (id, header) de un FASTA plano o .gz."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                yield header.split()[0], header


def sprot_acc(header_id: str) -> str:
    """`sp|P12345|NAME` -> `P12345`; deja igual un id sin pipes."""
    parts = header_id.split("|")
    return parts[1] if len(parts) >= 3 and parts[0] == "sp" else header_id
