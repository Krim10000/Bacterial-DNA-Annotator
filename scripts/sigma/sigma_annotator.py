#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Núcleo de anotación de sigma para genomas bacterianos.

Lee las tablas de regulación ya calculadas (`genes.tsv`, `operons.tsv`,
`operon_genes.tsv`, `promoters.tsv`, `regulators_inventory.tsv`), escanea
promotores con las PWMs de sigma de RegulonDB y cruza la evidencia con:

  * arquitectura del promotor (cajas y espaciador),
  * inventario de factores sigma del genoma,
  * contexto del operón (función de los genes),
  * expresión, si el usuario aporta una tabla,
  * ubicación (promotor aguas arriba del primer gen).

No usa Streamlit, no descarga datos y no escribe fuera de `out_dir`.

Salida por genoma:
  sigma_factors.tsv       una fila por factor sigma del inventario
  sigma_promoters.tsv     una fila por operón/promotor con su sigma asignado
  sigma_summary.json      conteos, cobertura y metadatos
  sigma_report.md         informe legible con evidencia
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

# --------------------------------------------------------------------------
# Constantes
# --------------------------------------------------------------------------
SIGMA_ORDER = ("sigma70", "sigma54", "sigma32", "sigma38", "sigma28", "sigma24")
SIGMA_CLASS = {
    "sigma70": "housekeeping",
    "sigma54": "alternative",
    "sigma32": "stress",
    "sigma38": "stationary",
    "sigma28": "flagellar",
    "sigma24": "ecf",
    "sigma19": "ecf",
    "sigmaE": "ecf",
    "sigmaB": "stress",
    "sigmaH": "sporulation",
    "sigmaK": "sporulation",
    "sigmaF": "sporulation",
    "sigma_unresolved": "unknown",
}
# Rango de espaciador entre la caja aguas arriba y la caja aguas abajo.
# sigma54 usa -24/-12, no -35/-10; se documenta en la evidencia.
SIGMA_SPACER = {
    "sigma70": (15, 19),
    "sigma54": (10, 11),
    "sigma32": (12, 16),
    "sigma38": (15, 19),
    "sigma28": (15, 19),
    "sigma24": (14, 18),
}
SIGMA_ARCH = {
    "sigma70": "-35/-10",
    "sigma54": "-24/-12",
    "sigma32": "-35/-10",
    "sigma38": "-35/-10",
    "sigma28": "-35/-10",
    "sigma24": "-35/-10",
}
SIGMA_KEYWORDS = {
    "sigma28": ("flagel", "flif", "flg", "flh", "fli", "che", "mot", "hook",
                "basal body", "chemotaxis"),
    "sigma32": ("heat", "chaper", "dnak", "grol", "htpg", "ibp", "lon",
                "clp", "proteas", "htpx"),
    "sigma38": ("stationary", "osmo", "trehalose", "glycogen", "catalase",
                "kate", "bola", "biofilm", "starvation"),
    "sigma54": ("nitrogen", "nif", "formate", "hydrogenase", "transport",
                "pili", "gluconate", "dicarboxylate", "c4"),
    "sigma24": ("envelope", "ecf", "iron", "zinc", "copper", "efflux",
                "membrane", "cell wall", "periplasm"),
    "sigma70": ("ribosom", "trna", "rrna", "translation", "housekeeping",
                "polymerase", "replication"),
}
OUTPUT_FACTOR_COLS = [
    "genome", "gene_id", "locus_tag", "gene_name", "product", "sigma_type",
    "sigma_class", "status", "confidence", "confidence_label", "confidence_score",
    "evidence", "detected_by", "pipeline_step", "promoter_hits",
    "source_status", "source_confidence", "source_evidence", "factor_note",
    "countable", "quality_flag", "ortholog_group", "copy_number", "is_expansion",
    "provenance",
]
OUTPUT_PROMOTER_COLS = [
    "genome", "operon_id", "operon_label", "contig", "strand", "first_gene",
    "first_gene_name", "sigma", "sigma_class", "sigma_arch", "spacer",
    "score", "zscore", "pvalue", "score_margin", "dist_to_start",
    "upstream_box", "downstream_box", "sigma_status", "confidence",
    "confidence_label", "confidence_score", "status", "evidence", "detected_by",
    "pipeline_step", "alternatives", "pathway_support", "pathway_hits",
    "expression_support", "expression_detail", "inventory_support",
    "inventory_detail", "location_support", "provenance",
]
OUTPUT_SUMMARY_KEYS = [
    "genome", "n_operons", "n_scanned", "n_assigned", "n_ambiguous",
    "n_unresolved", "n_factors", "n_factors_validated", "sigma_counts",
    "confidence_counts", "type_counts", "window", "null_windows", "seed",
    "version", "sha256",
]
VERSION = "1.0.0"
DEFAULT_WINDOW = 180
DEFAULT_NULL_WINDOWS = 400
DEFAULT_SEED = 17
DEFAULT_MAX_P = 0.05
DEFAULT_MIN_Z = 2.0
DEFAULT_MIN_MARGIN = 1.0

_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMP = str.maketrans("ACGT", "TGCA")


# --------------------------------------------------------------------------
# Utilidades de E/S
# --------------------------------------------------------------------------
def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            vals = []
            for col in cols:
                value = row.get(col, "")
                if isinstance(value, (list, dict)):
                    value = json.dumps(value, ensure_ascii=False, sort_keys=True)
                text = str(value).replace("\t", " ").replace("\n", " ")
                vals.append(text)
            fh.write("\t".join(vals) + "\n")


def _sha256(paths: Iterable[Path]) -> str:
    h = hashlib.sha256()
    for path in sorted(paths):
        if path.exists():
            h.update(path.name.encode())
            h.update(path.read_bytes())
    return h.hexdigest()


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").strip().lower())


# --------------------------------------------------------------------------
# Secuencia
# --------------------------------------------------------------------------
def _load_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def load_sequence(data_root: Path, genome: str) -> dict[str, str]:
    """Carga la secuencia usando `genes_meta.json` (raw o gbff)."""
    base = data_root / "regulation" / genome
    meta_path = base / "genes_meta.json"
    if not meta_path.exists():
        meta_path = base / "meta.json"
    if not meta_path.exists():
        return {}
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    src = meta.get("seq_source") or {}
    kind, path = src.get("kind"), src.get("path")
    candidates: list[Path] = []
    if kind and path:
        candidates.append(Path(path))
    # Fallback: algunas corridas no guardaron seq_source, pero el GBFF/FASTA
    # de NCBI sigue en data/annotations/.
    candidates += [
        data_root / "annotations" / f"{genome}_genomic.gbff.gz",
        data_root / "annotations" / f"{genome}_genomic.gbff",
        data_root / "annotations" / f"{genome}_genomic.fna.gz",
        data_root / "annotations" / f"{genome}_genomic.fna",
        data_root / "raw" / f"{genome}_genomic.fna.gz",
    ]
    for p in candidates:
        if not p.exists():
            continue
        if p.suffix == ".gz" and ".gbff" in p.name:
            try:
                from Bio import SeqIO
            except ImportError:
                continue
            out = {}
            with gzip.open(p, "rt") as fh:
                for rec in SeqIO.parse(fh, "genbank"):
                    out[rec.id] = str(rec.seq).upper()
            if out:
                return out
        if ".gbff" in p.name:
            try:
                from Bio import SeqIO
            except ImportError:
                continue
            out = {}
            for rec in SeqIO.parse(str(p), "genbank"):
                out[rec.id] = str(rec.seq).upper()
            if out:
                return out
        if ".fna" in p.name:
            out = _load_fasta(p)
            if out:
                return out
    return {}


def upstream_seq(seq: str, start: int, end: int, strand: int, n: int) -> str:
    """n nt aguas arriba del 5' del gen, en orientación de transcripción."""
    if not seq:
        return ""
    if strand == 1:
        a = max(0, start - 1 - n)
        return seq[a:start - 1]
    a = min(len(seq), end)
    return revcomp(seq[a:a + n])


# --------------------------------------------------------------------------
# PWMs
# --------------------------------------------------------------------------
def load_pwm_library(data_root: Path, pwm_path: Path | str | None = None) -> dict[str, dict[str, Any]]:
    """Carga una biblioteca de PWMs; acepta una o dos cajas por sigma."""
    path = Path(pwm_path) if pwm_path else (data_root / "regulondb" / "gt" / "promoter_pwms.tsv")
    rows = _read_tsv(path)
    out: dict[str, dict[str, Any]] = {}
    for row in rows:
        sigma = (row.get("sigma") or "").strip()
        box = (row.get("box") or "").strip()
        if sigma not in SIGMA_ORDER or box not in ("-35", "-10"):
            continue
        try:
            pwm = json.loads(row.get("pwm") or "[]")
            width = int(row.get("width") or len(pwm))
        except (TypeError, ValueError, json.JSONDecodeError):
            continue
        if not pwm or width <= 0:
            continue
        rec = out.setdefault(sigma, {
            "sigma": sigma,
            "class": SIGMA_CLASS.get(sigma, "unknown"),
            "arch": SIGMA_ARCH.get(sigma, "unknown"),
            "source": str(path),
            "n_sites": 0,
        })
        rec["n_sites"] = max(int(rec.get("n_sites") or 0),
                             int(row.get("n_sites") or 0))
        rec["up" if box == "-35" else "down"] = {
            "width": width,
            "pwm": np.asarray(pwm, dtype=float),
            "consensus": row.get("consensus", ""),
            "info_bits": float(row.get("info_bits") or 0.0),
            "n_sites": int(row.get("n_sites") or 0),
        }
    # sigma54 se etiqueta -35/-10 en la fuente, pero su geometría es -24/-12.
    for sigma, rec in out.items():
        if sigma == "sigma54":
            rec["arch"] = "-24/-12"
    return {k: v for k, v in out.items()
            if "up" in v or "down" in v}


def _encode(seq: str) -> np.ndarray:
    arr = np.full(len(seq), -1, dtype=np.int16)
    for base, idx in _BASE_IDX.items():
        arr[np.frombuffer(seq.encode("ascii", errors="ignore"), dtype=np.uint8) == ord(base)] = idx
    return arr


def _pwm_scores(seq: str, pwm: np.ndarray) -> np.ndarray:
    """Score de la PWM en cada posición de `seq` (ventana deslizante)."""
    width = pwm.shape[0]
    if len(seq) < width:
        return np.empty(0, dtype=float)
    codes = _encode(seq)
    n = len(codes) - width + 1
    scores = np.zeros(n, dtype=float)
    for k in range(width):
        block = codes[k:k + n]
        valid = block >= 0
        safe = np.where(valid, block, 0)
        vals = pwm[k][safe]
        scores += np.where(valid, vals, -6.0)
    return scores


def _best_geometry(up_scores: np.ndarray, down_scores: np.ndarray,
                   up_w: int, down_w: int, spacers: tuple[int, int],
                   window_len: int) -> dict[str, Any] | None:
    best = None
    for spacer in range(spacers[0], spacers[1] + 1):
        offset = up_w + spacer
        n = min(len(up_scores) - offset, len(down_scores) - offset)
        if n <= 0:
            continue
        totals = up_scores[:n] + down_scores[offset:offset + n]
        j = int(np.argmax(totals))
        score = float(totals[j])
        if best is None or score > best["score"]:
            up_pos = j
            down_pos = j + offset
            best = {
                "score": score,
                "spacer": spacer,
                "up_pos": up_pos,
                "down_pos": down_pos,
                "dist_to_start": window_len - (down_pos + down_w),
                "up_box": "",  # se rellena en el llamador con la secuencia
                "down_box": "",
            }
    return best


def score_sigma(seq: str, rec: dict[str, Any], window_len: int) -> dict[str, Any] | None:
    """Mejor arquitectura de un sigma; soporta una caja única."""
    up, down = rec.get("up"), rec.get("down")
    if up is None and down is not None:
        scores = _pwm_scores(seq, down["pwm"])
        if len(scores) == 0:
            return None
        j = int(np.argmax(scores))
        return {
            "score": float(scores[j]), "spacer": None,
            "up_pos": None, "down_pos": j,
            "dist_to_start": window_len - (j + down["width"]),
            "up_box": "", "down_box": seq[j:j + down["width"]],
            "sigma": rec["sigma"], "class": rec["class"],
            "arch": "single-box", "upstream_box": "",
            "downstream_box": seq[j:j + down["width"]],
            "upstream_pwm": "", "downstream_pwm": down["consensus"],
        }
    if down is None and up is not None:
        scores = _pwm_scores(seq, up["pwm"])
        if len(scores) == 0:
            return None
        j = int(np.argmax(scores))
        return {
            "score": float(scores[j]), "spacer": None,
            "up_pos": j, "down_pos": None,
            "dist_to_start": window_len - (j + up["width"]),
            "up_box": seq[j:j + up["width"]], "down_box": "",
            "sigma": rec["sigma"], "class": rec["class"],
            "arch": "single-box", "upstream_box": seq[j:j + up["width"]],
            "downstream_box": "", "upstream_pwm": up["consensus"],
            "downstream_pwm": "",
        }
    if up is None or down is None:
        return None
    up_scores = _pwm_scores(seq, up["pwm"])
    down_scores = _pwm_scores(seq, down["pwm"])
    if len(up_scores) == 0 or len(down_scores) == 0:
        return None
    best = _best_geometry(up_scores, down_scores, up["width"], down["width"],
                          SIGMA_SPACER.get(rec["sigma"], (14, 19)), window_len)
    if best is None:
        return None
    seq_up = seq[best["up_pos"]:best["up_pos"] + up["width"]]
    seq_down = seq[best["down_pos"]:best["down_pos"] + down["width"]]
    best.update({
        "sigma": rec["sigma"],
        "class": rec["class"],
        "arch": rec["arch"],
        "upstream_box": seq_up,
        "downstream_box": seq_down,
        "upstream_pwm": up["consensus"],
        "downstream_pwm": down["consensus"],
    })
    return best


# --------------------------------------------------------------------------
# Fondo nulo
# --------------------------------------------------------------------------
def _random_windows(seqs: dict[str, str], length: int, n: int,
                    seed: int) -> list[str]:
    rng = np.random.default_rng(seed)
    usable = [(contig, seq) for contig, seq in seqs.items() if len(seq) >= length]
    if not usable:
        return []
    weights = np.array([len(seq) for _, seq in usable], dtype=float)
    weights /= weights.sum()
    out = []
    for _ in range(n):
        idx = int(rng.choice(len(usable), p=weights))
        seq = usable[idx][1]
        start = int(rng.integers(0, len(seq) - length + 1))
        out.append(seq[start:start + length])
    return out


def build_null(seqs: dict[str, str], pwms: dict[str, dict[str, Any]],
               window: int, n: int, seed: int) -> dict[str, dict[str, Any]]:
    """Distribución nula del mejor score de cada sigma en ventanas aleatorias."""
    windows = _random_windows(seqs, window, n, seed)
    raw: dict[str, list[float]] = defaultdict(list)
    for w in windows:
        for sigma, rec in pwms.items():
            hit = score_sigma(w, rec, len(w))
            if hit is not None:
                raw[sigma].append(hit["score"])
    out = {}
    for sigma, vals in raw.items():
        arr = np.asarray(vals, dtype=float)
        out[sigma] = {
            "n": int(arr.size),
            "mean": float(arr.mean()) if arr.size else 0.0,
            "sd": float(arr.std(ddof=1)) if arr.size > 1 else 1.0,
            "sorted": np.sort(arr),
        }
    return out


def _pvalue(score: float, null: dict[str, Any] | None) -> float:
    if not null or null["n"] == 0:
        return 1.0
    return float((1 + int(np.sum(null["sorted"] >= score))) / (1 + null["n"]))


def _zscore(score: float, null: dict[str, Any] | None) -> float:
    if not null or null["n"] < 2:
        return 0.0
    sd = null["sd"] or 1.0
    return float((score - null["mean"]) / sd)


# --------------------------------------------------------------------------
# Evidencia contextual
# --------------------------------------------------------------------------
def _gene_text(row: dict[str, str]) -> str:
    return _norm(" ".join([
        row.get("gene_name", ""), row.get("product", ""), row.get("fam", ""),
    ]))


def pathway_support(sigma: str, genes: list[dict[str, str]]) -> tuple[bool, list[str]]:
    """Coherencia entre el tipo de sigma y la función del operón."""
    words = SIGMA_KEYWORDS.get(sigma, ())
    hits: list[str] = []
    for gene in genes:
        txt = _gene_text(gene)
        for word in words:
            if word and word in txt:
                label = (gene.get("gene_name") or gene.get("locus_tag")
                         or gene.get("gene_id") or "gen")
                hits.append(f"{label}:{word}")
                break
    return bool(hits), hits[:8]


def load_expression(data_root: Path, genome: str) -> dict[str, float]:
    """Carga una tabla de expresión opcional.

    Formatos aceptados: TSV con columnas `gene_id`/`locus_tag`/`gene` y
    `tpm`/`expression`/`value`. Si no existe, devuelve {}.
    """
    candidates = [
        data_root / "expression" / f"{genome}.tsv",
        data_root / "expression" / f"{genome}_expression.tsv",
        data_root / "regulation" / genome / "expression.tsv",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        return {}
    rows = _read_tsv(path)
    if not rows:
        return {}
    keys = ("gene_id", "locus_tag", "gene", "id")
    vals = ("tpm", "expression", "value", "counts")
    key = next((k for k in keys if k in rows[0]), None)
    val = next((k for k in vals if k in rows[0]), None)
    if not key or not val:
        return {}
    out: dict[str, float] = {}
    for row in rows:
        ident = (row.get(key) or "").strip()
        if not ident:
            continue
        try:
            out[ident] = float(row.get(val) or 0.0)
        except ValueError:
            continue
    return out


def expression_support(genes: list[dict[str, str]],
                       expression: dict[str, float]) -> tuple[str, str]:
    if not expression:
        return "not_available", "sin tabla de expresión"
    measured = []
    for gene in genes:
        ident = (gene.get("gene_id") or gene.get("locus_tag") or "").strip()
        if ident in expression:
            measured.append(expression[ident])
    if not measured:
        return "not_available", "genes del operón sin medida"
    frac = sum(1 for v in measured if v > 0) / len(measured)
    if frac >= 0.5:
        return "supported", f"{sum(1 for v in measured if v > 0)}/{len(measured)} con TPM>0"
    return "weak", f"{sum(1 for v in measured if v > 0)}/{len(measured)} con TPM>0"


def inventory_support(inventory: dict[str, dict[str, str]], sigma: str,
                      first_gene: dict[str, str]) -> tuple[bool, str]:
    """¿El genoma tiene un factor de ese tipo, y está cerca del operón?"""
    candidates = []
    for row in inventory.values():
        if (row.get("family") or "").strip() != sigma:
            continue
        candidates.append(row)
    if not candidates:
        return False, "sin factor sigma de ese tipo en el inventario"
    # Preferir el factor del mismo contig/hebra y más cercano.
    best = candidates[0]
    best_d = None
    for row in candidates:
        if row.get("contig") and first_gene.get("contig") and row["contig"] != first_gene["contig"]:
            continue
        try:
            d = abs(int(row.get("start", 0)) - int(first_gene.get("start", 0)))
        except (TypeError, ValueError):
            d = None
        if d is not None and (best_d is None or d < best_d):
            best, best_d = row, d
    label = best.get("gene_name") or best.get("gene_id") or "sigma"
    if best_d is None:
        return True, f"{label} ({best.get('status', '')})"
    return True, f"{label} a {best_d} pb ({best.get('status', '')})"


# --------------------------------------------------------------------------
# Confianza
# --------------------------------------------------------------------------
def _confidence(pvalue: float, z: float, margin: float,
                inv: bool, pathway: bool, expr: str,
                max_p: float, min_z: float, min_margin: float) -> tuple[str, str, float]:
    """Devuelve (confidence, confidence_label, score 0-1)."""
    # Score continuo: mezcla de evidencia estadística y contextual.
    stat = max(0.0, min(1.0, (z / 5.0) * 0.6 + (1.0 - min(pvalue / max_p, 1.0)) * 0.4))
    context = 0.0
    if inv:
        context += 0.15
    if pathway:
        context += 0.25
    if expr == "supported":
        context += 0.20
    elif expr == "weak":
        context += 0.08
    score = min(1.0, stat * 0.6 + context)
    # La confianza media/fuerte exige contexto funcional o expresión: el PWM
    # solo no basta para distinguir sigma70 de sigma38/sigma24.
    context_ok = pathway or expr == "supported"
    strong = (pvalue <= max_p / 5 and z >= max(3.0, min_z + 1)
              and margin >= max(2.0, min_margin) and inv and context_ok)
    # "media" = PWM + arquitectura + inventario; el contexto funcional es lo
    # que la promueve a "fuerte". No es una confirmación experimental.
    medium = (pvalue <= max_p and z >= min_z
              and margin >= max(1.0, min_margin) and inv)
    if strong:
        return "high", "sospecha_fuerte", round(score, 3)
    if medium:
        return "medium", "sospecha_media", round(score, 3)
    if pvalue <= max_p * 2 and z >= 1.0:
        return "low", "sospecha_baja", round(score, 3)
    return "low", "sin_sospecha", round(score, 3)


def _status_for(confidence: str, margin: float, pvalue: float,
                max_p: float) -> str:
    if confidence in ("high", "medium"):
        return "consensus"
    if pvalue <= max_p and margin < DEFAULT_MIN_MARGIN:
        return "anomalous"
    if pvalue <= max_p * 2:
        return "novel"
    return "unresolved"


# --------------------------------------------------------------------------
# Anotación
# --------------------------------------------------------------------------
def _first_gene_map(data_root: Path, genome: str) -> dict[str, dict[str, str]]:
    rows = _read_tsv(data_root / "regulation" / genome / "operon_genes.tsv")
    out = {}
    for row in rows:
        if (row.get("position") or "").strip() == "1":
            out[row["operon_id"]] = row
    return out


def _genes_by_id(data_root: Path, genome: str) -> dict[str, dict[str, str]]:
    rows = _read_tsv(data_root / "regulation" / genome / "genes.tsv")
    return {(r.get("gene_id") or "").strip(): r for r in rows if r.get("gene_id")}


def _op_genes(data_root: Path, genome: str) -> dict[str, list[str]]:
    rows = _read_tsv(data_root / "regulation" / genome / "operon_genes.tsv")
    out: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        out[row["operon_id"]].append((row.get("gene_id") or "").strip())
    return out


def _inventory(data_root: Path, genome: str,
               genes_by_id: dict[str, dict[str, str]] | None = None) -> dict[str, dict[str, str]]:
    rows = _read_tsv(data_root / "regulation" / genome / "regulators_inventory.tsv")
    out = {}
    for row in rows:
        key = (row.get("gene_id") or row.get("gene_name") or "").strip()
        if not key:
            continue
        gene = (genes_by_id or {}).get(key, {})
        for field in ("locus_tag", "contig", "start", "end", "strand", "gene_name"):
            if not row.get(field) and gene.get(field):
                row[field] = gene[field]
        out[key] = row
    return out


def factor_quality(row: dict[str, str]) -> tuple[str, str, bool]:
    """Devuelve (quality_flag, factor_note, countable) para un factor.

    Un hit de un solo dominio sin producto no es un factor sigma contable.
    """
    family = (row.get("sigma_type") or row.get("family") or "sigma_unresolved").strip()
    status = (row.get("source_status") or row.get("status") or "").strip()
    product = (row.get("product") or "").strip()
    evidence = (row.get("source_evidence") or row.get("evidence") or "")
    m = re.search(r"domains:([^;]*)", evidence)
    doms = []
    if m:
        doms = [d.strip() for d in m.group(1).split(",") if d.strip() and d.strip() != "-"]
    has_product = bool(product) and "hypothetical" not in product.lower()
    if family == "sigma_unresolved":
        flag = "partial_domain" if not has_product and len(doms) <= 1 else "unresolved"
        return flag, "sin familia asignable", False
    if status == "validated":
        return "validated", "", True
    if status == "consensus":
        return "consensus", "", bool(has_product)
    if status == "novel":
        if has_product and len(doms) >= 2:
            return "novel_two_domains", "familia inferida sin nombre canónico", True
        return "partial_domain", "hit parcial sin producto", False
    if not has_product:
        return "no_product", "sin nombre/producto", False
    return "fragment", "evidencia incompleta", False


def annotate_genome(data_root: Path, genome: str, out_dir: Path | None = None,
                    *, window: int = DEFAULT_WINDOW,
                    null_windows: int = DEFAULT_NULL_WINDOWS,
                    seed: int = DEFAULT_SEED,
                    max_p: float = DEFAULT_MAX_P,
                    min_z: float = DEFAULT_MIN_Z,
                    min_margin: float = DEFAULT_MIN_MARGIN,
                    use_expression: bool = True,
                    pwm_path: Path | str | None = None,
                    overwrite: bool = True) -> dict[str, Any]:
    """Anota sigma para un genoma y devuelve rutas + métricas."""
    base = data_root / "regulation" / genome
    if not base.exists():
        raise FileNotFoundError(f"no existe data/regulation/{genome}")
    out = Path(out_dir) if out_dir else data_root / "sigma" / genome
    out.mkdir(parents=True, exist_ok=True)

    operons = _read_tsv(base / "operons.tsv")
    first = _first_gene_map(data_root, genome)
    genes_by_id = _genes_by_id(data_root, genome)
    op_genes = _op_genes(data_root, genome)
    inventory = _inventory(data_root, genome, genes_by_id)
    seqs = load_sequence(data_root, genome)
    pwms = load_pwm_library(data_root, pwm_path)
    expression = load_expression(data_root, genome) if use_expression else {}
    if not operons:
        raise ValueError(f"sin operones para {genome}")
    if not pwms:
        raise ValueError("sin PWMs de sigma en data/regulondb/gt/promoter_pwms.tsv")

    sequence_available = bool(seqs)
    null = build_null(seqs, pwms, window, null_windows, seed) if seqs else {}

    promoter_rows: list[dict[str, Any]] = []
    sigma_counts: Counter[str] = Counter()
    conf_counts: Counter[str] = Counter()
    assigned = ambiguous = unresolved = 0

    for op in operons:
        oid = op["operon_id"]
        fg = first.get(oid)
        if not fg:
            unresolved += 1
            promoter_rows.append({
                "genome": genome, "operon_id": oid, "operon_label": op.get("label", ""),
                "contig": op.get("contig", ""), "strand": op.get("strand", ""),
                "sigma_status": "unresolved", "confidence": "low",
                "confidence_label": "sin_sospecha", "status": "unresolved",
                "evidence": "sin primer gen del operón",
                "detected_by": "operon_genes", "pipeline_step": "sigma.scan.v1",
                "provenance": json.dumps({"operon": oid, "genome": genome}),
            })
            continue
        contig = fg.get("contig") or op.get("contig") or ""
        seq = seqs.get(contig, "")
        if not seq:
            unresolved += 1
            promoter_rows.append({
                "genome": genome, "operon_id": oid, "operon_label": op.get("label", ""),
                "contig": contig, "strand": fg.get("strand", ""),
                "first_gene": fg.get("gene_id", ""),
                "first_gene_name": fg.get("gene_name", ""),
                "sigma_status": "unresolved", "confidence": "low",
                "confidence_label": "sin_sospecha", "status": "unresolved",
                "evidence": "sin secuencia del contig",
                "detected_by": "sequence", "pipeline_step": "sigma.scan.v1",
                "provenance": json.dumps({"operon": oid, "contig": contig}),
            })
            continue
        try:
            start, end, strand = int(fg["start"]), int(fg["end"]), int(fg["strand"])
        except (KeyError, ValueError):
            unresolved += 1
            continue
        w = upstream_seq(seq, start, end, strand, window)
        if len(w) < 30:
            unresolved += 1
            continue

        hits = {}
        for sigma, rec in pwms.items():
            hit = score_sigma(w, rec, len(w))
            if hit is None:
                continue
            hit["z"] = round(_zscore(hit["score"], null.get(sigma)), 3)
            hit["pvalue"] = round(_pvalue(hit["score"], null.get(sigma)), 5)
            hits[sigma] = hit
        if not hits:
            unresolved += 1
            continue

        ranked = sorted(hits.values(),
                        key=lambda h: (h["z"], h["pvalue"], h["score"]),
                        reverse=True)
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        margin = round(best["z"] - second["z"], 3) if second else 99.0

        op_gene_ids = op_genes.get(oid, [])
        op_gene_rows = [genes_by_id[g] for g in op_gene_ids if g in genes_by_id]
        p_support, p_hits = pathway_support(best["sigma"], op_gene_rows)
        e_status, e_detail = expression_support(op_gene_rows, expression)
        inv, inv_detail = inventory_support(inventory, best["sigma"], fg)
        confidence, label, cscore = _confidence(
            best["pvalue"], best["z"], margin, inv, p_support, e_status,
            max_p, min_z, min_margin)
        status = _status_for(confidence, margin, best["pvalue"], max_p)
        if confidence in ("high", "medium"):
            sigma_status = "assigned"
            assigned += 1
        elif best["pvalue"] <= max_p * 2 or margin < DEFAULT_MIN_MARGIN:
            sigma_status = "ambiguous"
            ambiguous += 1
        else:
            sigma_status = "unresolved"
            unresolved += 1
        sigma_counts[best["sigma"]] += 1
        conf_counts[label] += 1

        detected = ["promoter_pwm:" + best["sigma"],
                    "architecture:" + best["arch"],
                    "spacer:" + str(best["spacer"])]
        if inv:
            detected.append("inventory:" + best["sigma"])
        if p_support:
            detected.append("pathway:" + best["sigma"])
        if e_status == "supported":
            detected.append("expression")
        detected.append("location:upstream")

        evidence = (
            f"PWM {best['sigma']} z={best['z']} p={best['pvalue']} "
            f"spacer={best['spacer']} ({best['arch']}); "
            f"margen={margin}; inventario={inv_detail}; "
            f"contexto={'sí' if p_support else 'no'}; "
            f"expresión={e_status}"
        )
        alternatives = [
            {"sigma": h["sigma"], "z": h["z"], "pvalue": h["pvalue"],
             "score": round(h["score"], 3), "spacer": h["spacer"]}
            for h in ranked[:3]
        ]
        promoter_rows.append({
            "genome": genome, "operon_id": oid, "operon_label": op.get("label", ""),
            "contig": contig, "strand": strand, "first_gene": fg.get("gene_id", ""),
            "first_gene_name": fg.get("gene_name", ""),
            "sigma": best["sigma"], "sigma_class": best["class"],
            "sigma_arch": best["arch"], "spacer": best["spacer"],
            "score": round(best["score"], 3), "zscore": best["z"],
            "pvalue": best["pvalue"], "score_margin": margin,
            "dist_to_start": best["dist_to_start"],
            "upstream_box": best.get("upstream_box", ""),
            "downstream_box": best.get("downstream_box", ""),
            "sigma_status": sigma_status, "confidence": confidence,
            "confidence_label": label, "confidence_score": cscore,
            "status": status, "evidence": evidence,
            "detected_by": ";".join(detected),
            "pipeline_step": "sigma.scan.v1",
            "alternatives": alternatives,
            "pathway_support": "yes" if p_support else "no",
            "pathway_hits": p_hits,
            "expression_support": e_status,
            "expression_detail": e_detail,
            "inventory_support": "yes" if inv else "no",
            "inventory_detail": inv_detail,
            "location_support": "upstream",
            "provenance": json.dumps({
                "genome": genome, "operon_id": oid, "contig": contig,
                "window": window, "seed": seed, "version": VERSION,
            }, ensure_ascii=False, sort_keys=True),
        })

    # Factores sigma del inventario, enriquecidos con sus promotores.
    factor_hits = Counter()
    for row in promoter_rows:
        if row.get("sigma"):
            factor_hits[row["sigma"]] += 1
    factor_rows = []
    for key, row in inventory.items():
        sigma_type = (row.get("family") or "sigma_unresolved").strip()
        src_status = (row.get("status") or "").strip()
        src_conf = (row.get("confidence") or "").strip()
        conf = ("high" if src_status == "validated" else
                "medium" if src_status == "consensus" else "low")
        label = ("sospecha_fuerte" if conf == "high" else
                 "sospecha_media" if conf == "medium" else "sospecha_baja")
        product = (row.get("product") or "").lower()
        note = ""
        if "leader peptide" in product:
            note = "posible péptido líder, no factor sigma"
        elif "anti-sigma" in product or "anti sigma" in product:
            note = "posible anti-sigma, no factor sigma"
        qflag, qnote, countable = factor_quality({
            "sigma_type": sigma_type, "source_status": src_status,
            "product": row.get("product", ""), "source_evidence": row.get("evidence", ""),
        })
        if qnote and not note:
            note = qnote
        factor_rows.append({
            "genome": genome, "gene_id": row.get("gene_id", ""),
            "locus_tag": row.get("locus_tag", ""),
            "gene_name": row.get("gene_name", "") or row.get("locus_tag", ""),
            "product": row.get("product", ""), "sigma_type": sigma_type,
            "sigma_class": SIGMA_CLASS.get(sigma_type, "unknown"),
            "status": src_status or "unresolved", "confidence": conf,
            "confidence_label": label,
            "confidence_score": 0.9 if conf == "high" else 0.6 if conf == "medium" else 0.3,
            "evidence": row.get("evidence", ""),
            "detected_by": "inventory:" + (row.get("attempted", "") or "nombre"),
            "pipeline_step": "sigma.inventory.v1",
            "promoter_hits": factor_hits.get(sigma_type, 0),
            "source_status": src_status, "source_confidence": src_conf,
            "source_evidence": row.get("evidence", ""),
            "factor_note": note,
            "countable": "yes" if countable else "no",
            "quality_flag": qflag,
            "ortholog_group": "",
            "copy_number": "",
            "is_expansion": "",
            "provenance": json.dumps({"genome": genome, "gene_id": row.get("gene_id", "")},
                                     ensure_ascii=False, sort_keys=True),
        })

    # Preservar la agrupación global si ya existe (20_group_sigma.py).
    old_factors = {r.get("gene_id", ""): r for r in _read_tsv(out / "sigma_factors.tsv")}
    for row in factor_rows:
        old = old_factors.get(row.get("gene_id", ""))
        if old:
            for col in ("ortholog_group", "copy_number", "is_expansion"):
                if old.get(col):
                    row[col] = old[col]

    # Orden estable
    promoter_rows.sort(key=lambda r: (str(r.get("contig", "")), str(r.get("operon_id", ""))))
    factor_rows.sort(key=lambda r: (str(r.get("sigma_type", "")), str(r.get("gene_id", ""))))

    _write_tsv(out / "sigma_promoters.tsv", promoter_rows, OUTPUT_PROMOTER_COLS)
    _write_tsv(out / "sigma_factors.tsv", factor_rows, OUTPUT_FACTOR_COLS)
    summary = {
        "genome": genome,
        "n_operons": len(operons),
        "n_scanned": len([r for r in promoter_rows if r.get("sigma")]),
        "n_assigned": assigned,
        "n_ambiguous": ambiguous,
        "n_unresolved": unresolved,
        "n_factors": len(factor_rows),
        "n_factors_validated": sum(1 for r in factor_rows if r["status"] == "validated"),
        "n_factors_countable": sum(1 for r in factor_rows if r["countable"] == "yes"),
        "n_factors_unresolved": sum(1 for r in factor_rows if r["sigma_type"] == "sigma_unresolved"),
        "n_ortholog_groups": len({r["ortholog_group"] for r in factor_rows if r.get("ortholog_group")}),
        "n_expansions": len({r["ortholog_group"] for r in factor_rows
                             if r.get("is_expansion") == "yes"}),
        "sigma_counts": dict(sigma_counts.most_common()),
        "confidence_counts": dict(conf_counts.most_common()),
        "type_counts": dict(Counter(r.get("sigma_type") for r in factor_rows).most_common()),
        "window": window,
        "null_windows": null_windows,
        "sequence_available": sequence_available,
        "seed": seed,
        "version": VERSION,
        "sha256": _sha256([out / "sigma_promoters.tsv", out / "sigma_factors.tsv"]),
    }
    (out / "sigma_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # Informe legible
    lines = [
        f"# Sigma — {genome}", "",
        f"Operones: {summary['n_operons']} · escaneados: {summary['n_scanned']} · "
        f"asignados: {assigned} · ambiguos: {ambiguous} · sin resolver: {unresolved}",
        "", "## Conteo por tipo", "",
        "| sigma | promotores |", "|---|---:|",
    ]
    for sigma, n in summary["sigma_counts"].items():
        lines.append(f"| {sigma} | {n} |")
    lines += ["", "## Confianza", "", "| etiqueta | n |", "|---|---:|"]
    for label, n in summary["confidence_counts"].items():
        lines.append(f"| {label} | {n} |")
    lines += ["", "## Top promotores", "",
              "| operón | gen | sigma | confianza | evidencia |", "|---|---|---|---|---|"]
    top = sorted([r for r in promoter_rows if r.get("sigma")],
                 key=lambda r: (-float(r.get("confidence_score") or 0),
                                -float(r.get("zscore") or 0)))[:20]
    for r in top:
        lines.append(
            f"| {r.get('operon_label') or r.get('operon_id')} | "
            f"{r.get('first_gene_name') or r.get('first_gene')} | "
            f"{r.get('sigma')} | {r.get('confidence_label')} | "
            f"{r.get('evidence')} |")
    lines += ["", "## Límites", "",
              "- `consensus` no es confirmación experimental; es acuerdo entre PWM, "
              "arquitectura y contexto.",
              "- `sin_sospecha`/`unresolved` no significa ausencia de promotor: "
              "significa que ninguna PWM superó el umbral.",
              "- La expresión es opcional; si no hay tabla, se declara "
              "`not_available`.",
              "- sigma54 se escanea con geometría -24/-12, no -35/-10.",
              ""]
    (out / "sigma_report.md").write_text("\n".join(lines), encoding="utf-8")
    return {
        "genome": genome, "out_dir": str(out), "summary": summary,
        "promoters": out / "sigma_promoters.tsv",
        "factors": out / "sigma_factors.tsv",
        "report": out / "sigma_report.md",
    }


# --------------------------------------------------------------------------
# Autoprueba / calibración
# --------------------------------------------------------------------------
def selftest() -> str:
    """Prueba rápida: PWM canónica, geometría y confianza."""
    pwm = np.array([
        [1.0, -1.0, -1.0, -1.0], [-1.0, 1.0, -1.0, -1.0],
        [-1.0, -1.0, 1.0, -1.0], [-1.0, -1.0, -1.0, 1.0],
        [1.0, -1.0, -1.0, -1.0], [-1.0, 1.0, -1.0, -1.0],
    ])
    rec = {"sigma": "sigma70", "class": "housekeeping", "arch": "-35/-10",
           "up": {"width": 6, "pwm": pwm, "consensus": "ACGTAC", "info_bits": 4.0},
           "down": {"width": 6, "pwm": pwm, "consensus": "ACGTAC", "info_bits": 4.0}}
    seq = "AAAAAA" + "ACGTAC" + "A" * 17 + "ACGTAC" + "A" * 20
    hit = score_sigma(seq, rec, len(seq))
    assert hit is not None
    assert hit["score"] > 0
    assert hit["spacer"] == 17
    # spacer inválido para sigma70 no debe dar el mejor hit con 17
    seq2 = "A" * 30 + "ACGTAC" + "A" * 30 + "ACGTAC"
    hit2 = score_sigma(seq2, rec, len(seq2))
    assert hit2 is None or hit2["spacer"] != 17
    conf, label, score = _confidence(0.001, 4.0, 2.5, True, True, "supported",
                                     0.05, 2.0, 1.0)
    assert conf == "high" and label == "sospecha_fuerte" and 0 <= score <= 1
    conf2, label2, _ = _confidence(0.04, 2.2, 1.6, True, False, "not_available",
                                   0.05, 2.0, 1.0)
    assert conf2 == "medium" and label2 == "sospecha_media"
    conf3, label3, _ = _confidence(0.2, 0.5, 0.5, False, False, "not_available",
                                   0.05, 2.0, 1.0)
    assert conf3 == "low" and label3 == "sin_sospecha"
    return "scanner sigma70, geometría y confianza OK"
