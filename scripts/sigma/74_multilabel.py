#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Candidatos sigma multietiqueta por TSS y evaluación top-k contra el gold.

El pipeline sigma clásico escribe una fila por operón. Eso oculta que un mismo
TSS puede tener más de un factor sigma curado y que el ranking top-k contiene
varias etiquetas plausibles. Este script construye una tabla de candidatos:

    una fila por (genoma, TSS/operón, sigma candidato)

y la evalúa como multietiqueta contra ``data/sigma/gold/<genoma>/tss_sigma.tsv``.

Fuentes (en orden de preferencia):

1. ``data/sigma/layers/<genoma>/layer_sigma_features.tsv``
   Una fila por TSS x sigma. Es la única fuente que aporta un TSS real.
2. ``data/sigma/<genoma>/sigma_promoters.tsv`` (respaldo)
   Una fila por operón con ``alternatives`` (top-3 JSON). El ancla del
   ``operon_id`` NO es un TSS: se escribe como ``coord_type=operon_anchor`` y
   no se evalúa contra un gold de TSS. No se inventa ninguna coordenada.

Reglas de construcción:

- ranking por ``z`` descendente (desempate: ``pvalue``, orden canónico, nombre);
- top-k configurable (por defecto 3);
- umbrales opcionales ``--min-z``, ``--max-pvalue`` y ``--min-pwm-score``;
- ``confidence_label`` heurística basada solo en señal (z, p, margen);
- ``refined_confidence`` opcional para sigma54/sigma38 usando margen y
  distancia, sin cambiar nunca el sigma asignado.

No usa las columnas ``true``, ``is_top1``, ``n_true``, ``gold_confidence`` ni
``source`` de las features: son metadatos del gold, no señales del pipeline.
El gold solo se lee después de construir los candidatos, para evaluar.

Uso:
    python scripts/sigma/74_multilabel.py --genome ecoli_k12 --genome subtilis
    python scripts/sigma/74_multilabel.py --all
    python scripts/sigma/74_multilabel.py --genome ecoli_k12 --top-k 4 \
        --min-z 0 --max-pvalue 0.10
    python scripts/sigma/74_multilabel.py --genome ecoli_k12 --source promoters

Salidas por genoma:
    data/sigma/multilabel/<genoma>/sigma_candidates.tsv
    data/sigma/multilabel/<genoma>/metrics.json
    data/sigma/multilabel/<genoma>/report.md
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]

VERSION = "1.0.0"
DEFAULT_TOP_K = 3
DEFAULT_FEATURES_ROOT = ROOT / "data" / "sigma" / "layers"
DEFAULT_PROMOTERS_ROOT = ROOT / "data" / "sigma"
DEFAULT_GOLD_ROOT = ROOT / "data" / "sigma" / "gold"
DEFAULT_OUT_ROOT = ROOT / "data" / "sigma" / "multilabel"

SIGMA_ORDER = ("sigma70", "sigma54", "sigma32", "sigma38", "sigma28", "sigma24")
SIGMA_INDEX = {s: i for i, s in enumerate(SIGMA_ORDER)}

# Equivalencias de letras de B. subtilis usadas por el gold DBTBS.
BSU_MAP = {
    "A": "sigma70", "SIGA": "sigma70",
    "B": "sigma38", "SIGB": "sigma38",
    "D": "sigma28", "SIGD": "sigma28",
    "I": "sigma32", "SIGI": "sigma32",
    "L": "sigma54", "SIGL": "sigma54",
    "E": "", "SIGE": "",
}

# Espaciadores canónicos (copia documentada de SIGMA_SPACER del núcleo; no se
# importa el anotador para que este script no dependa de sus efectos laterales).
SIGMA_SPACER_RANGE = {
    "sigma70": (15, 19),
    "sigma54": (10, 11),
    "sigma32": (12, 16),
    "sigma38": (15, 19),
    "sigma28": (15, 19),
    "sigma24": (14, 18),
}

#: Umbrales de ``confidence_label``. Son heurísticos y no están calibrados.
CONFIDENCE_THRESHOLDS = {
    "alta": "pvalue <= 0.05 y margin >= 1.0",
    "media": "pvalue <= 0.05, o (margin >= 1.0 y z >= 1.0)",
    "baja": "z >= 0.0, o pvalue <= 0.20",
    "muy_baja": "resto",
}

#: Reglas de refinamiento opcionales. No cambian el sigma, solo la confianza.
REFINEMENT_RULES = {
    "sigma54": {
        "name": "sigma54_margen_distancia",
        "description": "pvalue <= 0.05 AND margin >= 1.0 AND 0 <= dist_to_start <= 60",
        "pvalue_max": 0.05,
        "margin_min": 1.0,
        "dist_max": 60.0,
    },
    "sigma38": {
        "name": "sigma38_margen_distancia",
        "description": "margin >= 1.0 AND 0 <= dist_to_start <= 60",
        "pvalue_max": None,
        "margin_min": 1.0,
        "dist_max": 60.0,
    },
}

CANDIDATE_COLUMNS = [
    "genome",
    "contig",
    "tss",
    "strand",
    "locus_id",
    "coord_type",
    "operon_anchor",
    "rank",
    "sigma",
    "z",
    "pvalue",
    "margin",
    "score",
    "spacer",
    "spacer_dev",
    "dist_to_start",
    "inventory",
    "confidence_label",
    "refined_confidence",
    "refinement_rule",
    "evidence",
    "source",
]

COL_ALIASES = {
    "tss": ("tss", "pos", "position", "start", "coordinate"),
    "sigma": ("sigma", "sigma_factor", "factor", "sigmas", "sigma_type"),
    "contig": ("contig", "seqid", "chromosome", "replicon"),
    "strand": ("strand", "hebra"),
    "operon": ("operon", "operon_id", "target_operon", "tu"),
    "gene": ("gene", "first_gene", "gene_name", "target_gene"),
    "confidence": ("confidence", "confianza"),
}

LIMITATIONS = [
    "Las features de capa están ancladas al conjunto de TSS del gold: la tabla "
    "cubre esos TSS, no un escaneo genómico independiente.",
    "El gold es positivo: una predicción fuera del gold no se cuenta como falso "
    "positivo; se reporta aparte como 'fuera del gold'.",
    "La evaluación usa el mismo gold que inspiró los umbrales de confianza y "
    "las reglas de refinamiento; no es una validación held-out.",
    "La precisión top-k incluye k etiquetas por TSS; no es comparable con la "
    "precisión top-1 de una fila por operón.",
    "En B. subtilis se usan PWMs de E. coli; las métricas mezclan transferencia "
    "de motivo con asignación.",
    "El respaldo por sigma_promoters.tsv es a nivel de operón: su coordenada es "
    "un ancla de operón, no un TSS, y no se evalúa contra el gold de TSS.",
    "Las reglas de refinamiento son filtros de alta precisión, no una nueva "
    "llamada de sigma ni una validación experimental.",
]


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def _f(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    if math.isnan(out) or math.isinf(out):
        return float(default)
    return out


def _maybe_f(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def _int_or_none(value: Any) -> int | None:
    out = _maybe_f(value)
    if out is None:
        return None
    return int(round(out))


def _clean(value: Any) -> str:
    return str(value or "").replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: Iterable[dict[str, Any]],
               columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: _format_cell(row.get(c), c) for c in columns})


def _format_cell(value: Any, column: str = "") -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if column == "pvalue":
            return f"{value:.5f}"
        if column in {"z", "margin", "score"}:
            return f"{value:.3f}"
        if column in {"spacer_dev"}:
            return f"{value:g}"
        if value.is_integer():
            return str(int(value))
        return f"{value:g}"
    return _clean(value)


def _fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if math.isnan(value):
            return "—"
        return f"{value:.{digits}f}"
    return str(value)


def _margin_rows(rows: list[dict[str, Any]]) -> None:
    """Completa ``margin`` con z menos el mejor rival del mismo locus."""
    zs = sorted((_f(r.get("z")) for r in rows), reverse=True)
    for row in rows:
        if row.get("margin") is not None:
            continue
        z = _f(row.get("z"))
        others = list(zs)
        try:
            others.remove(z)
        except ValueError:
            pass
        row["margin"] = z - (others[0] if others else z)


def _norm_strand(value: Any) -> str:
    low = str(value or "").strip().lower()
    if low in {"+", "1", "forward", "fwd", "plus"}:
        return "+"
    if low in {"-", "-1", "reverse", "rev", "minus"}:
        return "-"
    return ""


def _canonical_sigma(raw: Any) -> str:
    """Normaliza sigma70/σ70/letras DBTBS. Devuelve '' si no se reconoce."""
    text = str(raw or "").strip()
    if not text or text.lower() in {"?", "na", "n/a", "none", "unlabeled", "-"}:
        return ""
    up = text.upper()
    if up in BSU_MAP:
        return BSU_MAP[up]
    low = text.lower().replace("σ", "sigma").replace(" ", "")
    if low in SIGMA_ORDER:
        return low
    match = re.fullmatch(r"sigma0*(70|54|32|38|28|24)", low)
    if match:
        return "sigma" + match.group(1)
    if low in {"sigmaa", "sigma_a"}:
        return "sigma70"
    return ""


def _pick(row: dict[str, str], key: str) -> str:
    for alias in COL_ALIASES.get(key, (key,)):
        value = row.get(alias)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _split_sigmas(value: str) -> list[str]:
    raw = str(value or "").strip()
    if not raw or raw.lower() in {"?", "na", "n/a", "none", "unlabeled", "-"}:
        return []
    return [p.strip() for p in re.split(r"[;,|]+", raw) if p.strip()]


def _rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def _rank_key(row: dict[str, Any]) -> tuple[Any, ...]:
    sigma = str(row.get("sigma") or "")
    return (-_f(row.get("z")), _f(row.get("pvalue"), 1.0),
            SIGMA_INDEX.get(sigma, 99), sigma)


def _confidence_label(z: float, pvalue: float, margin: float) -> str:
    if pvalue <= 0.05 and margin >= 1.0:
        return "alta"
    if pvalue <= 0.05 or (margin >= 1.0 and z >= 1.0):
        return "media"
    if z >= 0.0 or pvalue <= 0.20:
        return "baja"
    return "muy_baja"


def _rule_passes(row: dict[str, Any], sigma: str) -> bool:
    rule = REFINEMENT_RULES.get(sigma)
    if not rule:
        return False
    margin = _f(row.get("margin"), -999.0)
    pvalue = _f(row.get("pvalue"), 1.0)
    dist = _maybe_f(row.get("dist_to_start"))
    if margin < float(rule["margin_min"]):
        return False
    if dist is None or dist < 0.0 or dist > float(rule["dist_max"]):
        return False
    if rule.get("pvalue_max") is not None and pvalue > float(rule["pvalue_max"]):
        return False
    return True


def _refined_confidence(row: dict[str, Any], base: str, enabled: bool) -> str:
    if not enabled:
        return base
    sigma = str(row.get("sigma") or "")
    if sigma not in REFINEMENT_RULES:
        return base
    if _rule_passes(row, sigma):
        return "alta"
    # Una confianza alta sin el margen/distancia exigido se rebaja a media; no
    # se toca el sigma ni se sube ninguna etiqueta que no pase la regla.
    return "media" if base == "alta" else base


# ---------------------------------------------------------------------------
# Construcción de candidatos
# ---------------------------------------------------------------------------
@dataclass
class BuildResult:
    genome: str
    source: str
    source_path: Path
    candidates: list[dict[str, Any]]
    all_rows: list[dict[str, Any]]
    stats: dict[str, Any]
    top_k: int
    refine: bool
    filters: dict[str, Any] = field(default_factory=dict)


def _layer_key(row: dict[str, str]) -> tuple[str, int, str] | None:
    contig = (row.get("contig") or "").strip()
    tss = _int_or_none(row.get("tss"))
    strand = _norm_strand(row.get("strand"))
    if not contig or tss is None or tss <= 0 or not strand:
        return None
    return (contig, tss, strand)


def _assign_ranks(groups: dict[Any, list[dict[str, Any]]],
                  genome: str, top_k: int, refine: bool,
                  source: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    all_rows: list[dict[str, Any]] = []
    for _key, rows in groups.items():
        rows.sort(key=_rank_key)
        for rank, row in enumerate(rows, start=1):
            row["genome"] = genome
            row["rank"] = rank
            base = _confidence_label(
                _f(row.get("z")), _f(row.get("pvalue"), 1.0),
                _f(row.get("margin"), 0.0))
            row["confidence_label"] = base
            row["refined_confidence"] = _refined_confidence(row, base, refine)
            row["refinement_rule"] = (
                REFINEMENT_RULES[row["sigma"]]["description"]
                if refine and row.get("sigma") in REFINEMENT_RULES else "")
            row["source"] = source
            all_rows.append(row)

    def sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
        tss = row.get("tss")
        anchor = row.get("operon_anchor")
        coord = tss if tss not in (None, "") else anchor
        return (str(row.get("contig") or ""),
                _f(coord, 10**12),
                str(row.get("strand") or ""),
                int(row.get("rank") or 0),
                str(row.get("locus_id") or ""))

    all_rows.sort(key=sort_key)
    candidates = [r for r in all_rows if int(r.get("rank") or 0) <= top_k]
    return all_rows, candidates


def _make_stats(source: str, source_path: Path, genome: str,
                n_input_rows: int, n_input_groups: int,
                n_invalid: int, n_filtered: int,
                all_rows: list[dict[str, Any]],
                candidates: list[dict[str, Any]],
                top_k: int, filters: dict[str, Any],
                coordinate_note: str) -> dict[str, Any]:
    by_tss: dict[Any, int] = Counter()
    for row in all_rows:
        by_tss[row["_key"]] += 1
    rank_counts = Counter(int(r["rank"]) for r in candidates)
    return {
        "source": source,
        "source_path": _rel(source_path),
        "coordinate_note": coordinate_note,
        "n_input_rows": n_input_rows,
        "n_input_loci": n_input_groups,
        "n_invalid_rows": n_invalid,
        "n_filtered_rows": n_filtered,
        "n_candidates_all": len(all_rows),
        "n_candidates_topk": len(candidates),
        "n_tss_with_candidates": len(by_tss),
        "n_tss_without_candidates": max(0, n_input_groups - len(by_tss)),
        "n_tss_fewer_than_k": sum(1 for n in by_tss.values() if n < top_k),
        "mean_candidates_per_tss": (len(all_rows) / len(by_tss)) if by_tss else None,
        "rank_counts_topk": {str(k): rank_counts.get(k, 0) for k in range(1, top_k + 1)},
        "sigma_counts_topk": dict(Counter(r["sigma"] for r in candidates)),
        "top_k": top_k,
        "filters": filters,
        "genome": genome,
    }


def build_candidates_from_layers(path: Path, genome: str,
                                 top_k: int = DEFAULT_TOP_K,
                                 refine: bool = True,
                                 min_z: float | None = None,
                                 max_pvalue: float | None = None,
                                 min_pwm_score: float | None = None) -> BuildResult:
    """Construye candidatos TSS-level desde layer_sigma_features.tsv."""
    input_rows = _read_tsv(path)
    groups: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    valid_keys: set[tuple[str, int, str]] = set()
    n_invalid = 0
    n_filtered = 0

    for raw in input_rows:
        key = _layer_key(raw)
        sigma = _canonical_sigma(raw.get("sigma"))
        if key is None or not sigma:
            n_invalid += 1
            continue
        valid_keys.add(key)
        z = _f(raw.get("z"))
        pvalue = _f(raw.get("pvalue"), 1.0)
        score = _f(raw.get("score"))
        if min_z is not None and z < min_z:
            n_filtered += 1
            continue
        if max_pvalue is not None and pvalue > max_pvalue:
            n_filtered += 1
            continue
        if min_pwm_score is not None and score < min_pwm_score:
            n_filtered += 1
            continue
        margin = _maybe_f(raw.get("margin"))
        spacer = _int_or_none(raw.get("spacer"))
        spacer_dev = _maybe_f(raw.get("spacer_dev"))
        dist = _maybe_f(raw.get("dist_to_start"))
        inventory = (raw.get("inventory") or "").strip().lower()
        if inventory not in {"yes", "no"}:
            inventory = ""
        groups[key].append({
            "_key": key,
            "contig": key[0],
            "tss": key[1],
            "strand": key[2],
            "locus_id": f"{key[0]}:{key[1]}:{key[2]}",
            "coord_type": "tss",
            "operon_anchor": "",
            "sigma": sigma,
            "z": z,
            "pvalue": pvalue,
            "margin": margin,
            "score": score,
            "spacer": spacer,
            "spacer_dev": spacer_dev,
            "dist_to_start": dist,
            "inventory": inventory,
            "evidence": "",
        })

    # Margen de respaldo si la tabla no lo trae: z menos el mejor rival del TSS.
    for rows in groups.values():
        _margin_rows(rows)
        for row in rows:
            row["evidence"] = (
                f"PWM {row['sigma']} z={_f(row['z']):.3f} p={_f(row['pvalue'], 1.0):.5f}; "
                f"margen={_f(row['margin']):.3f}; "
                f"spacer={row.get('spacer') if row.get('spacer') is not None else 'NA'} "
                f"(dev={_fmt(row.get('spacer_dev'), 0)}); "
                f"dist_inicio={_fmt(row.get('dist_to_start'), 0)}; "
                f"inventario={row.get('inventory') or 'NA'}"
            )

    all_rows, candidates = _assign_ranks(groups, genome, top_k, refine, "layers")
    filters = {"min_z": min_z, "max_pvalue": max_pvalue,
               "min_pwm_score": min_pwm_score}
    stats = _make_stats(
        "layers", path, genome, len(input_rows), len(valid_keys), n_invalid,
        n_filtered, all_rows, candidates, top_k, filters,
        "TSS 1-based (primera base transcrita) tomado de layer_sigma_features.tsv; "
        "no se infiere ni se desplaza ninguna coordenada.")
    return BuildResult(genome, "layers", path, candidates, all_rows, stats,
                       top_k, refine, filters)


def _parse_operon_id(value: str) -> tuple[str, str, int | None, str]:
    parts = str(value or "").split(":")
    if len(parts) >= 4:
        return parts[0], parts[1], _int_or_none(parts[2]), _norm_strand(parts[3])
    return "", "", None, ""

def build_candidates_from_promoters(path: Path, genome: str,
                                    top_k: int = DEFAULT_TOP_K,
                                    refine: bool = True,
                                    min_z: float | None = None,
                                    max_pvalue: float | None = None,
                                    min_pwm_score: float | None = None) -> BuildResult:
    """Respaldo a nivel de operón desde sigma_promoters.tsv.

    La coordenada del ``operon_id`` es un ancla de operón, no un TSS. Se
    conserva como ``operon_anchor`` con ``coord_type=operon_anchor`` y ``tss``
    vacío. No se evalúa contra un gold de TSS.
    """
    input_rows = _read_tsv(path)
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    loci_seen: set[tuple[str, str]] = set()
    n_invalid = 0
    n_filtered = 0

    for raw in input_rows:
        operon_id = _clean(raw.get("operon_id") or "")
        _genome, contig, anchor, strand = _parse_operon_id(operon_id)
        if not contig:
            contig = _clean(raw.get("contig"))
        if not strand:
            strand = _norm_strand(raw.get("strand"))
        if not contig or not strand:
            n_invalid += 1
            continue
        locus = operon_id or f"{contig}:{anchor if anchor is not None else ''}:{strand}"
        key = ("__operon__", locus)
        loci_seen.add(key)

        alternatives: list[dict[str, Any]] = []
        try:
            parsed = json.loads(raw.get("alternatives") or "[]")
            if isinstance(parsed, list):
                alternatives = [x for x in parsed if isinstance(x, dict)]
        except (TypeError, ValueError, json.JSONDecodeError):
            alternatives = []
        if not alternatives:
            alternatives = [{
                "sigma": raw.get("sigma"),
                "z": raw.get("zscore"),
                "pvalue": raw.get("pvalue"),
                "score": raw.get("score"),
                "spacer": raw.get("spacer"),
            }]

        by_sigma: dict[str, dict[str, Any]] = {}
        for alt in alternatives:
            sigma = _canonical_sigma(alt.get("sigma"))
            if not sigma:
                n_invalid += 1
                continue
            z = _f(alt.get("z", alt.get("zscore")))
            pvalue = _f(alt.get("pvalue"), 1.0)
            score = _f(alt.get("score"))
            if min_z is not None and z < min_z:
                n_filtered += 1
                continue
            if max_pvalue is not None and pvalue > max_pvalue:
                n_filtered += 1
                continue
            if min_pwm_score is not None and score < min_pwm_score:
                n_filtered += 1
                continue
            prev = by_sigma.get(sigma)
            if prev is None or z > _f(prev.get("z")):
                spacer = _int_or_none(alt.get("spacer"))
                spacer_dev: float | None = None
                if spacer is not None:
                    lo, hi = SIGMA_SPACER_RANGE.get(sigma, (0, 0))
                    spacer_dev = (0.0 if lo <= spacer <= hi
                                  else float(min(abs(spacer - lo),
                                                 abs(spacer - hi))))
                by_sigma[sigma] = {
                    "_key": key,
                    "contig": contig,
                    "tss": "",
                    "strand": strand,
                    "locus_id": locus,
                    "coord_type": "operon_anchor",
                    "operon_anchor": anchor if anchor is not None else "",
                    "sigma": sigma,
                    "z": z,
                    "pvalue": pvalue,
                    "margin": None,
                    "score": score,
                    "spacer": spacer,
                    "spacer_dev": spacer_dev,
                    "dist_to_start": _maybe_f(raw.get("dist_to_start")),
                    "inventory": "",
                    "evidence": "",
                }
        rows = list(by_sigma.values())
        if not rows:
            continue
        _margin_rows(rows)
        for row in rows:
            context = _clean(raw.get("evidence"))
            row["evidence"] = (
                f"{context}; " if context else ""
            ) + (
                f"candidato de operón {row['sigma']} z={_f(row['z']):.3f} "
                f"p={_f(row['pvalue'], 1.0):.5f}; margen_top3={_f(row['margin']):.3f}; "
                f"spacer={row.get('spacer') if row.get('spacer') is not None else 'NA'}; "
                f"dist_inicio={_fmt(row.get('dist_to_start'), 0)}; "
                f"coordenada=ancla_de_operón (NO TSS)"
            )
        groups[key].extend(rows)

    all_rows, candidates = _assign_ranks(groups, genome, top_k, refine, "promoters")
    filters = {"min_z": min_z, "max_pvalue": max_pvalue,
               "min_pwm_score": min_pwm_score}
    stats = _make_stats(
        "promoters", path, genome, len(input_rows), len(loci_seen), n_invalid,
        n_filtered, all_rows, candidates, top_k, filters,
        "Coordenada de operón (ancla del operon_id), NO TSS. No se evalúa "
        "contra el gold de TSS ni se convierte a coordenada de inicio de "
        "transcripción.")
    return BuildResult(genome, "promoters", path, candidates, all_rows, stats,
                       top_k, refine, filters)


def resolve_source(genome: str, source: str, features_root: Path,
                   promoters_root: Path) -> tuple[str, Path]:
    layer_path = features_root / genome / "layer_sigma_features.tsv"
    promoter_path = promoters_root / genome / "sigma_promoters.tsv"
    if source == "layers":
        if not layer_path.exists():
            raise FileNotFoundError(f"sin features de capa: {layer_path}")
        return "layers", layer_path
    if source == "promoters":
        if not promoter_path.exists():
            raise FileNotFoundError(f"sin sigma_promoters.tsv: {promoter_path}")
        return "promoters", promoter_path
    if layer_path.exists():
        return "layers", layer_path
    if promoter_path.exists():
        return "promoters", promoter_path
    raise FileNotFoundError(
        f"{genome}: no hay layer_sigma_features.tsv ni sigma_promoters.tsv")


# ---------------------------------------------------------------------------
# Gold y evaluación
# ---------------------------------------------------------------------------
def load_gold(path: Path) -> tuple[dict[tuple[str, int, str], set[str]],
                                   dict[str, Any]]:
    """Carga el gold TSS y agrupa sus sigma por (contig, tss, strand)."""
    if not path.exists():
        raise FileNotFoundError(f"gold no encontrado: {path}")
    rows = _read_tsv(path)
    gold: dict[tuple[str, int, str], set[str]] = defaultdict(set)
    per_sigma: Counter[str] = Counter()
    n_invalid = 0
    n_unmapped = 0
    unmapped_values: Counter[str] = Counter()

    for raw in rows:
        contig = _pick(raw, "contig")
        tss = _int_or_none(_pick(raw, "tss"))
        strand = _norm_strand(_pick(raw, "strand"))
        if not contig or tss is None or tss <= 0 or not strand:
            n_invalid += 1
            continue
        labels = _split_sigmas(_pick(raw, "sigma"))
        for label in labels:
            sigma = _canonical_sigma(label)
            if not sigma:
                n_unmapped += 1
                unmapped_values[label] += 1
                continue
            gold[(contig, tss, strand)].add(sigma)

    per_sigma = Counter(s for labels in gold.values() for s in labels)
    n_tss = len(gold)
    n_multi = sum(1 for labels in gold.values() if len(labels) > 1)
    stats = {
        "path": _rel(path),
        "n_rows": len(rows),
        "n_tss": n_tss,
        "n_single": n_tss - n_multi,
        "n_multi": n_multi,
        "n_labels": sum(len(v) for v in gold.values()),
        "per_sigma": dict(sorted(per_sigma.items())),
        "n_invalid_rows": n_invalid,
        "n_unmapped_labels": n_unmapped,
        "unmapped_values": dict(unmapped_values),
        "coordinate_system": "1-based TSS; key=(contig, pos, strand)",
    }
    return dict(gold), stats


def _make_records(all_rows: list[dict[str, Any]],
                  gold: dict[tuple[str, int, str], set[str]]
                  ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_key: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in all_rows:
        by_key[row["_key"]].append(row)
    for rows in by_key.values():
        rows.sort(key=lambda r: int(r.get("rank") or 0))

    records: list[dict[str, Any]] = []
    for key, true in gold.items():
        rows = by_key.get(key, [])
        preds: list[str] = []
        seen: set[str] = set()
        for row in rows:
            sigma = row["sigma"]
            if sigma not in seen:
                seen.add(sigma)
                preds.append(sigma)
        records.append({
            "_key": key,
            "true": set(true),
            "preds_all": preds,
            "n_true": len(true),
            "n_pred_all": len(preds),
        })

    gold_keys = set(gold)
    outside_keys = [k for k in by_key if k not in gold_keys]
    outside_rows = sum(len(by_key[k]) for k in outside_keys)
    info = {
        "n_gold_tss": len(gold),
        "n_tss_with_candidates": sum(1 for r in records if r["preds_all"]),
        "n_gold_tss_without_candidates": sum(1 for r in records if not r["preds_all"]),
        "n_candidate_tss_outside_gold": len(outside_keys),
        "n_candidate_rows_outside_gold": outside_rows,
    }
    return records, info


def _per_sigma_counts(records: list[dict[str, Any]], k: int,
                      labels: Sequence[str]) -> dict[str, dict[str, int]]:
    out = {label: {"support": 0, "predicted": 0, "tp": 0, "fp": 0, "fn": 0}
           for label in labels}
    for rec in records:
        true = rec["true"]
        pred = set(rec["preds_all"][:k])
        for label in true:
            out.setdefault(label, {"support": 0, "predicted": 0,
                                   "tp": 0, "fp": 0, "fn": 0})
            out[label]["support"] += 1
        for label in pred:
            out.setdefault(label, {"support": 0, "predicted": 0,
                                   "tp": 0, "fp": 0, "fn": 0})
            out[label]["predicted"] += 1
        for label in pred & true:
            out[label]["tp"] += 1
        for label in pred - true:
            out[label]["fp"] += 1
        for label in true - pred:
            out[label]["fn"] += 1
    return out


def _prf(tp: int, fp: int, fn: int) -> dict[str, float | None]:
    precision = tp / (tp + fp) if (tp + fp) else None
    recall = tp / (tp + fn) if (tp + fn) else None
    if precision is None or recall is None:
        f1 = 0.0 if recall is not None else None
    else:
        f1 = (2 * precision * recall / (precision + recall)
              if (precision + recall) else 0.0)
    return {"precision": precision, "recall": recall, "f1": f1}


def _micro(per: dict[str, dict[str, int]]) -> dict[str, float | None]:
    tp = sum(v["tp"] for v in per.values())
    fp = sum(v["fp"] for v in per.values())
    fn = sum(v["fn"] for v in per.values())
    return _prf(tp, fp, fn)


def _macro(per: dict[str, dict[str, int]]) -> dict[str, Any]:
    ps: list[float] = []
    rs: list[float] = []
    f1s: list[float] = []
    for label, counts in per.items():
        if counts["support"] == 0:
            continue
        prf = _prf(counts["tp"], counts["fp"], counts["fn"])
        precision = prf["precision"] or 0.0
        recall = prf["recall"] or 0.0
        f1 = prf["f1"] if prf["f1"] is not None else 0.0
        ps.append(precision)
        rs.append(recall)
        f1s.append(f1)
    return {
        "n_labels": len(ps),
        "precision": (sum(ps) / len(ps)) if ps else None,
        "recall": (sum(rs) / len(rs)) if rs else None,
        "f1": (sum(f1s) / len(f1s)) if f1s else None,
    }


def _summarize_block(records: list[dict[str, Any]], k: int,
                     labels: Sequence[str]) -> dict[str, Any]:
    n = len(records)
    any_hits = 0
    all_hits = 0
    jaccards: list[float] = []
    n_true_labels = 0
    n_pred_labels = 0
    n_without = 0
    for rec in records:
        true = rec["true"]
        pred = set(rec["preds_all"][:k])
        n_true_labels += len(true)
        n_pred_labels += len(pred)
        if not pred:
            n_without += 1
        inter = true & pred
        if inter:
            any_hits += 1
        if true and true <= pred:
            all_hits += 1
        union = true | pred
        jaccards.append(len(inter) / len(union) if union else 1.0)
    per = _per_sigma_counts(records, k, labels)
    return {
        "k": k,
        "n_tss": n,
        "n_true_labels": n_true_labels,
        "n_pred_labels": n_pred_labels,
        "n_tss_without_candidates": n_without,
        "topk_any": any_hits,
        "topk_all": all_hits,
        "topk_any_rate": (any_hits / n) if n else None,
        "topk_all_rate": (all_hits / n) if n else None,
        "jaccard_mean": (sum(jaccards) / len(jaccards)) if jaccards else None,
        "micro": _micro(per),
        "macro": _macro(per),
    }


def _summary(blocks: dict[int, dict[str, Any]], top_k: int) -> dict[str, Any]:
    b1 = blocks.get(1)
    bk = blocks.get(top_k)
    b3 = blocks.get(3)

    def get(block: dict[str, Any] | None, field_name: str) -> Any:
        return block.get(field_name) if block else None

    return {
        "top_k": top_k,
        "n_tss": get(bk, "n_tss"),
        "n_true_labels": get(bk, "n_true_labels"),
        "n_pred_labels": get(bk, "n_pred_labels"),
        "n_tss_without_candidates": get(bk, "n_tss_without_candidates"),
        "top1_any": get(b1, "topk_any"),
        "top1_all": get(b1, "topk_all"),
        "top1_any_rate": get(b1, "topk_any_rate"),
        "top1_all_rate": get(b1, "topk_all_rate"),
        "top3_any": get(b3, "topk_any"),
        "top3_all": get(b3, "topk_all"),
        "top3_any_rate": get(b3, "topk_any_rate"),
        "top3_all_rate": get(b3, "topk_all_rate"),
        "topk_any": get(bk, "topk_any"),
        "topk_all": get(bk, "topk_all"),
        "topk_any_rate": get(bk, "topk_any_rate"),
        "topk_all_rate": get(bk, "topk_all_rate"),
        "jaccard_top1": get(b1, "jaccard_mean"),
        "jaccard_top3": get(b3, "jaccard_mean"),
        "jaccard_topk": get(bk, "jaccard_mean"),
        "micro": get(bk, "micro"),
        "macro": get(bk, "macro"),
        "top1_micro": get(b1, "micro"),
        "top1_macro": get(b1, "macro"),
        "by_k": [blocks[k] for k in sorted(blocks)],
    }


def _coverage(records: list[dict[str, Any]], k: int,
              labels: Sequence[str]) -> dict[str, Any]:
    multi = [r for r in records if r["n_true"] > 1]
    n = len(multi)
    any_hits = 0
    all_hits = 0
    labels_total = 0
    labels_covered = 0
    by_sigma = {label: {"support_multi": 0, "covered_multi": 0,
                        "outside_multi": 0} for label in labels}
    for rec in multi:
        true = rec["true"]
        pred = set(rec["preds_all"][:k])
        inter = true & pred
        if inter:
            any_hits += 1
        if true <= pred:
            all_hits += 1
        labels_total += len(true)
        labels_covered += len(inter)
        for label in true:
            slot = by_sigma.setdefault(
                label, {"support_multi": 0, "covered_multi": 0,
                        "outside_multi": 0})
            slot["support_multi"] += 1
            if label in pred:
                slot["covered_multi"] += 1
            else:
                slot["outside_multi"] += 1
    return {
        "k": k,
        "n_multi_tss": n,
        "multi_any": any_hits,
        "multi_all": all_hits,
        "multi_any_rate": (any_hits / n) if n else None,
        "multi_all_rate": (all_hits / n) if n else None,
        "multi_not_fully_covered": n - all_hits,
        "multi_missed_entirely": n - any_hits,
        "gold_labels_multi": labels_total,
        "gold_labels_covered": labels_covered,
        "gold_labels_outside": labels_total - labels_covered,
        "gold_labels_covered_rate": (labels_covered / labels_total) if labels_total else None,
        "by_sigma": by_sigma,
    }


def _per_sigma_table(records: list[dict[str, Any]], k: int,
                     labels: Sequence[str]) -> list[dict[str, Any]]:
    per_k = _per_sigma_counts(records, k, labels)
    per_1 = _per_sigma_counts(records, 1, labels)
    single = [r for r in records if r["n_true"] == 1]
    multi = [r for r in records if r["n_true"] > 1]
    single_counts = _per_sigma_counts(single, k, labels)
    multi_counts = _per_sigma_counts(multi, k, labels)
    out: list[dict[str, Any]] = []
    for label in labels:
        counts = per_k.get(label, {"support": 0, "predicted": 0,
                                   "tp": 0, "fp": 0, "fn": 0})
        prf = _prf(counts["tp"], counts["fp"], counts["fn"])
        c1 = per_1.get(label, {"support": 0, "predicted": 0,
                               "tp": 0, "fp": 0, "fn": 0})
        prf1 = _prf(c1["tp"], c1["fp"], c1["fn"])
        out.append({
            "sigma": label,
            "support": counts["support"],
            "predicted": counts["predicted"],
            "tp": counts["tp"],
            "fp": counts["fp"],
            "fn": counts["fn"],
            "precision": prf["precision"],
            "recall": prf["recall"],
            "f1": prf["f1"],
            "support_single": single_counts.get(label, {}).get("support", 0),
            "support_multi": multi_counts.get(label, {}).get("support", 0),
            "tp_multi": multi_counts.get(label, {}).get("tp", 0),
            "covered_multi": multi_counts.get(label, {}).get("tp", 0),
            "outside_multi": (multi_counts.get(label, {}).get("support", 0)
                              - multi_counts.get(label, {}).get("tp", 0)),
            "top1_precision": prf1["precision"],
            "top1_recall": prf1["recall"],
            "top1_f1": prf1["f1"],
        })
    return out


def _multi_examples(records: list[dict[str, Any]], k: int,
                    limit: int = 20) -> list[dict[str, Any]]:
    examples: list[dict[str, Any]] = []
    for rec in records:
        if rec["n_true"] <= 1:
            continue
        true = rec["true"]
        pred = rec["preds_all"][:k]
        inter = true & set(pred)
        union = true | set(pred)
        examples.append({
            "contig": rec["_key"][0],
            "tss": rec["_key"][1],
            "strand": rec["_key"][2],
            "true_sigmas": sorted(true, key=lambda s: SIGMA_INDEX.get(s, 99)),
            "pred_topk": pred,
            "missing": sorted(true - set(pred),
                              key=lambda s: SIGMA_INDEX.get(s, 99)),
            "jaccard": (len(inter) / len(union)) if union else 1.0,
        })
    examples.sort(key=lambda e: (len(e["missing"]) == 0, e["contig"], e["tss"]),
                  reverse=False)
    # Primero los que dejan más etiquetas fuera; luego por coordenada.
    examples.sort(key=lambda e: (-len(e["missing"]), e["contig"], e["tss"]))
    return examples[:limit]


def evaluate_candidates(result: BuildResult,
                        gold: dict[tuple[str, int, str], set[str]]) -> dict[str, Any]:
    if result.source != "layers":
        return {
            "status": "not_applicable",
            "reason": "candidatos a nivel de operón (operon_anchor), no de TSS; "
                      "no se fuerza un emparejamiento con el gold de TSS",
            "n_candidates": len(result.candidates),
            "n_loci": len({r["locus_id"] for r in result.candidates}),
        }
    if not gold:
        return {"status": "no_gold", "reason": "gold vacío o ausente"}

    records, info = _make_records(result.all_rows, gold)
    labels = [s for s in SIGMA_ORDER
              if s in {r["sigma"] for r in result.all_rows} or s in
              {x for true in gold.values() for x in true}]
    k_eval = result.top_k
    blocks = {k: _summarize_block(records, k, labels)
              for k in range(1, k_eval + 1)}
    single = [r for r in records if r["n_true"] == 1]
    multi = [r for r in records if r["n_true"] > 1]
    single_blocks = {k: _summarize_block(single, k, labels)
                     for k in range(1, k_eval + 1)} if single else {}
    multi_blocks = {k: _summarize_block(multi, k, labels)
                    for k in range(1, k_eval + 1)} if multi else {}

    overall = _summary(blocks, k_eval)
    overall["by_k"] = [blocks[k] for k in sorted(blocks)]
    single_summary = _summary(single_blocks, k_eval) if single_blocks else {
        "top_k": k_eval, "n_tss": 0}
    multi_summary = _summary(multi_blocks, k_eval) if multi_blocks else {
        "top_k": k_eval, "n_tss": 0}

    coverage = _coverage(records, k_eval, labels)
    max_sensitivity_k = max(
        k_eval,
        min(6, max((len(r["preds_all"]) for r in records), default=1)))
    sensitivity = [_coverage(records, k, labels)
                   for k in range(1, max_sensitivity_k + 1)]
    per_sigma = _per_sigma_table(records, k_eval, labels)
    examples = _multi_examples(records, k_eval, limit=20)

    return {
        "status": "ok",
        "k_eval": k_eval,
        "n_evaluated_tss": len(records),
        "n_single_tss": len(single),
        "n_multi_tss": len(multi),
        "labels": labels,
        "records_info": info,
        "overall": overall,
        "single": single_summary,
        "multi": multi_summary,
        "per_sigma": per_sigma,
        "multi_sigma_coverage": coverage,
        "topk_sensitivity": sensitivity,
        "multi_examples": examples,
    }


def evaluate_refinement(result: BuildResult,
                        gold: dict[tuple[str, int, str], set[str]]) -> dict[str, Any]:
    if not result.refine:
        return {"enabled": False, "rules": REFINEMENT_RULES}
    if result.source != "layers":
        return {"enabled": True, "status": "not_applicable",
                "reason": "las reglas se evalúan contra gold de TSS; los "
                          "candidatos de operón no comparten esa coordenada",
                "rules": REFINEMENT_RULES}
    if not gold:
        return {"enabled": True, "status": "no_gold",
                "rules": REFINEMENT_RULES}

    by_sigma: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in result.all_rows:
        by_sigma[row["sigma"]].append(row)

    out: dict[str, Any] = {"enabled": True, "scope": "all_signal_rows",
                           "rules": {}, "per_sigma": []}
    for sigma, rule in REFINEMENT_RULES.items():
        rows = by_sigma.get(sigma, [])
        pass_rows = [r for r in rows if _rule_passes(r, sigma)]
        tp = fp = outside = 0
        for row in pass_rows:
            true = gold.get(row["_key"])
            if true is None:
                outside += 1
            elif sigma in true:
                tp += 1
            else:
                fp += 1
        support = sum(1 for labels in gold.values() if sigma in labels)
        # Un gold con sigma y sin fila de features también es un FN de la regla.
        rows_by_key = {r["_key"]: r for r in rows}
        fn = 0
        for key, labels in gold.items():
            if sigma not in labels:
                continue
            row = rows_by_key.get(key)
            if row is None or not _rule_passes(row, sigma):
                fn += 1
        prf = _prf(tp, fp, fn)

        # Línea base: confianza alta antes de refinar.
        base_rows = [r for r in rows
                     if r.get("confidence_label") == "alta" and r["_key"] in gold]
        btp = sum(1 for r in base_rows if sigma in gold[r["_key"]])
        bfp = len(base_rows) - btp
        bfn = support - btp
        baseline = _prf(btp, bfp, max(0, bfn))

        n_refined_alta = sum(1 for r in rows
                             if r.get("refined_confidence") == "alta")
        upgraded = sum(1 for r in rows
                       if r.get("refined_confidence") == "alta"
                       and r.get("confidence_label") != "alta")
        downgraded = sum(1 for r in rows
                         if r.get("confidence_label") == "alta"
                         and r.get("refined_confidence") != "alta")
        n_in_table = sum(1 for r in result.candidates
                         if r["sigma"] == sigma)
        out["rules"][sigma] = rule
        out["per_sigma"].append({
            "sigma": sigma,
            "rule": rule["name"],
            "description": rule["description"],
            "n_rows": len(rows),
            "n_rows_in_candidate_table": n_in_table,
            "n_pass": len(pass_rows),
            "n_outside_gold_pass": outside,
            "support": support,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": prf["precision"],
            "recall": prf["recall"],
            "f1": prf["f1"],
            "baseline_alta_tp": btp,
            "baseline_alta_fp": bfp,
            "baseline_alta_precision": baseline["precision"],
            "baseline_alta_recall": baseline["recall"],
            "baseline_alta_f1": baseline["f1"],
            "n_refined_alta": n_refined_alta,
            "n_upgraded": upgraded,
            "n_downgraded": downgraded,
        })
    return out


# ---------------------------------------------------------------------------
# Reporte
# ---------------------------------------------------------------------------
def _md_table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> list[str]:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |",
             "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return lines


def _rate_text(value: Any, total: Any = None) -> str:
    if value is None:
        return "—"
    if total:
        return f"{value}/{total} ({_fmt(100.0 * value / total, 1)} %)"
    return str(value)


def build_report(genome: str, result: BuildResult,
                 gold_stats: dict[str, Any] | None,
                 evaluation: dict[str, Any],
                 refinement: dict[str, Any]) -> str:
    stats = result.stats
    cfg = {
        "top_k": result.top_k,
        "refine": result.refine,
        "filters": result.filters,
    }
    lines: list[str] = [
        f"# Candidatos sigma multietiqueta — `{genome}`",
        "",
        f"> Generado por `scripts/sigma/74_multilabel.py` v{VERSION}. "
        f"Fuente: `{stats['source_path']}`.",
        "",
        "## 1. Entrada y coordenada",
        "",
        f"- Fuente: **{result.source}**.",
        f"- Filas de entrada: {stats['n_input_rows']}.",
        f"- TSS/operones de entrada: {stats['n_input_loci']}.",
        f"- Candidatos generados (ranking completo): {stats['n_candidates_all']}.",
        f"- Candidatos en la tabla top-{result.top_k}: {stats['n_candidates_topk']}.",
        f"- TSS/operones con al menos un candidato: {stats['n_tss_with_candidates']}.",
        f"- TSS/operones sin candidatos tras umbrales: {stats['n_tss_without_candidates']}.",
        f"- Coordenada: {stats['coordinate_note']}",
        "",
        "La tabla no inventa TSS: usa la columna `tss` de las features de capa. "
        "Si se usa el respaldo de `sigma_promoters.tsv`, `tss` queda vacío y la "
        "coordenada se escribe en `operon_anchor` con "
        "`coord_type=operon_anchor`.",
        "",
        "## 2. Configuración",
        "",
        f"- `top_k`: {result.top_k}.",
        f"- Umbrales: `min_z={result.filters.get('min_z')}`, "
        f"`max_pvalue={result.filters.get('max_pvalue')}`, "
        f"`min_pwm_score={result.filters.get('min_pwm_score')}`.",
        f"- Refinamiento σ54/σ38: {'activo' if result.refine else 'desactivado'}.",
        "- `confidence_label`: " + "; ".join(
            f"`{label}` = {condition}"
            for label, condition in CONFIDENCE_THRESHOLDS.items()) + ".",
        "",
    ]

    if gold_stats:
        lines.extend([
            "## 3. Gold y cobertura",
            "",
            f"- Gold: `{gold_stats['path']}`.",
            f"- TSS gold: {gold_stats['n_tss']} "
            f"(single={gold_stats['n_single']}, multi={gold_stats['n_multi']}).",
            f"- Etiquetas gold: {gold_stats['n_labels']}.",
            f"- TSS evaluados: {evaluation.get('n_evaluated_tss', '—')}.",
            f"- TSS gold sin candidatos: "
            f"{evaluation.get('records_info', {}).get('n_gold_tss_without_candidates', '—')}.",
            f"- TSS candidatos fuera del gold: "
            f"{evaluation.get('records_info', {}).get('n_candidate_tss_outside_gold', '—')} "
            f"(no se cuentan como falsos positivos).",
            "",
        ])

        if evaluation.get("status") == "ok":
            ov = evaluation["overall"]
            k = evaluation["k_eval"]
            jaccard_rows: list[tuple[str, Any]] = [
                ("Jaccard top1", _fmt(ov["jaccard_top1"]))]
            if ov.get("jaccard_top3") is not None and k != 3:
                jaccard_rows.append(("Jaccard top3", _fmt(ov["jaccard_top3"])))
            if k != 1:
                jaccard_rows.append((f"Jaccard top{k}", _fmt(ov["jaccard_topk"])))
            lines.extend([
                "## 4. Métricas globales (top-k)",
                "",
                f"Evaluación con **top-{k}** como conjunto predicho. "
                "`top1_all` solo puede ser cierto en TSS con una etiqueta; "
                "`top3_all` es la métrica clave para multi-sigma.",
                "",
                *_md_table(
                    ["Métrica", "Valor"],
                    [
                        ("TSS evaluados", ov["n_tss"]),
                        ("Single / multi", f"{evaluation['n_single_tss']} / "
                                           f"{evaluation['n_multi_tss']}"),
                        ("top1_any", _rate_text(ov["top1_any"], ov["n_tss"])),
                        ("top1_all", _rate_text(ov["top1_all"], ov["n_tss"])),
                        ("top3_any", _rate_text(ov["top3_any"], ov["n_tss"])),
                        ("top3_all", _rate_text(ov["top3_all"], ov["n_tss"])),
                        *jaccard_rows,
                        ("Micro P / R / F1",
                         f"{_fmt(ov['micro']['precision'])} / "
                         f"{_fmt(ov['micro']['recall'])} / "
                         f"{_fmt(ov['micro']['f1'])}"),
                        ("Macro P / R / F1",
                         f"{_fmt(ov['macro']['precision'])} / "
                         f"{_fmt(ov['macro']['recall'])} / "
                         f"{_fmt(ov['macro']['f1'])}"),
                    ]),
                "",
                "### 4.1 Single vs multi",
                "",
                *_md_table(
                    ["Subconjunto", "n", "top1_any", "top1_all", "top3_any",
                     "top3_all", "Jaccard top3", "Micro F1"],
                    [
                        ("single", evaluation["single"].get("n_tss"),
                         _rate_text(evaluation["single"].get("top1_any"),
                                    evaluation["single"].get("n_tss")),
                         _rate_text(evaluation["single"].get("top1_all"),
                                    evaluation["single"].get("n_tss")),
                         _rate_text(evaluation["single"].get("top3_any"),
                                    evaluation["single"].get("n_tss")),
                         _rate_text(evaluation["single"].get("top3_all"),
                                    evaluation["single"].get("n_tss")),
                         _fmt(evaluation["single"].get("jaccard_top3")),
                         _fmt((evaluation["single"].get("micro") or {}).get("f1"))),
                        ("multi", evaluation["multi"].get("n_tss"),
                         _rate_text(evaluation["multi"].get("top1_any"),
                                    evaluation["multi"].get("n_tss")),
                         _rate_text(evaluation["multi"].get("top1_all"),
                                    evaluation["multi"].get("n_tss")),
                         _rate_text(evaluation["multi"].get("top3_any"),
                                    evaluation["multi"].get("n_tss")),
                         _rate_text(evaluation["multi"].get("top3_all"),
                                    evaluation["multi"].get("n_tss")),
                         _fmt(evaluation["multi"].get("jaccard_top3")),
                         _fmt((evaluation["multi"].get("micro") or {}).get("f1"))),
                    ]),
                "",
                "Lectura: en multi-sigma, `top1_all` es 0 por construcción; "
                "`top1_any` mide si el primer candidato es al menos uno de los "
                "sigma verdaderos. `top3_all` mide si el top-3 cubre todos.",
                "",
                "## 5. Por sigma (conjunto top-k)",
                "",
                *_md_table(
                    ["Sigma", "Support", "Pred", "TP", "FP", "FN", "P", "R",
                     "F1", "Support single", "Support multi", "Cubiertos multi"],
                    [
                        (r["sigma"], r["support"], r["predicted"], r["tp"],
                         r["fp"], r["fn"], _fmt(r["precision"]),
                         _fmt(r["recall"]), _fmt(r["f1"]),
                         r["support_single"], r["support_multi"],
                         f"{r['covered_multi']}/{r['support_multi']}")
                        for r in evaluation["per_sigma"]
                    ]),
                "",
                "La precisión se calcula sobre el conjunto top-k: por eso es "
                "menor que la precisión top-1. Los TSS fuera del gold no "
                "entran en FP.",
                "",
            ])

            lines.extend([
                "## 6. Cobertura multi-sigma por k",
                "",
                "Sensibilidad sobre el ranking completo de features, con los "
                "mismos umbrales; la fila "
                f"marcada con `*` es la tabla top-{k} entregada.",
                "",
                *_md_table(
                    ["k", "TSS multi", "Alguno cubierto", "Todos cubiertos",
                     "No cubiertos del todo", "Fuera del top-k (TSS)",
                     "Etiquetas cubiertas", "Etiquetas fuera"],
                    [
                        (f"{c['k']}{'*' if c['k'] == k else ''}",
                         c["n_multi_tss"],
                         _rate_text(c["multi_any"], c["n_multi_tss"]),
                         _rate_text(c["multi_all"], c["n_multi_tss"]),
                         c["multi_not_fully_covered"],
                         c["multi_missed_entirely"],
                         _rate_text(c["gold_labels_covered"],
                                    c["gold_labels_multi"]),
                         c["gold_labels_outside"])
                        for c in evaluation["topk_sensitivity"]
                    ]),
                "",
                "## 7. Refinamiento σ54/σ38",
                "",
            ])
            if refinement.get("enabled") and refinement.get("per_sigma"):
                lines.extend([
                    "El refinamiento **no cambia el sigma**: solo escribe "
                    "`refined_confidence`. Se evalúa sobre todas las filas de "
                    "señal, no solo sobre las que entran en top-k, para medir "
                    "la regla completa.",
                    "",
                    *_md_table(
                        ["Sigma", "Regla", "n pass", "TP", "FP", "FN", "P", "R",
                         "F1", "Alta base P/R", "Subidas", "Bajadas"],
                        [
                            (r["sigma"], r["description"], r["n_pass"],
                             r["tp"], r["fp"], r["fn"],
                             _fmt(r["precision"]), _fmt(r["recall"]),
                             _fmt(r["f1"]),
                             f"{_fmt(r['baseline_alta_precision'])} / "
                             f"{_fmt(r['baseline_alta_recall'])}",
                             r["n_upgraded"], r["n_downgraded"])
                            for r in refinement["per_sigma"]
                        ]),
                    "",
                ])
            elif refinement.get("enabled"):
                lines.extend([
                    "Refinamiento activo, pero no evaluable: "
                    f"{refinement.get('reason') or 'sin gold de TSS'}.",
                    "",
                ])
            else:
                lines.extend(["Refinamiento desactivado.", ""])

            lines.extend(["## 8. Ejemplos multi-sigma", ""])
            examples = evaluation.get("multi_examples", [])
            if examples:
                lines.extend(_md_table(
                    ["Contig", "TSS", "Hebra", "Verdad", "Top-k", "Faltan",
                     "Jaccard"],
                    [
                        (e["contig"], e["tss"], e["strand"],
                         ", ".join(e["true_sigmas"]),
                         ", ".join(e["pred_topk"]),
                         ", ".join(e["missing"]) or "—",
                         _fmt(e["jaccard"]))
                        for e in examples[:10]
                    ]))
                lines.append("")
            else:
                lines.extend(["No hay TSS multi-sigma en el gold de este "
                              "genoma.", ""])

        else:
            reason = evaluation.get("reason") or "evaluación no disponible"
            lines.extend([
                "## 4. Evaluación",
                "",
                f"No se calculan métricas de TSS: **{reason}**.",
                "",
            ])

    else:
        lines.extend([
            "## 3. Gold",
            "",
            "No se encontró gold para este genoma; solo se generaron "
            "candidatos, sin métricas supervisadas.",
            "",
        ])

    lines.extend([
        "## 9. Límites",
        "",
        *[f"- {item}" for item in LIMITATIONS],
        "",
        "## 10. Reproducir",
        "",
        "```bash",
        f"python scripts/sigma/74_multilabel.py --genome {genome} "
        f"--top-k {result.top_k}"
        + ("" if result.refine else " --no-refine"),
        "```",
        "",
        "Salidas: `data/sigma/multilabel/"
        f"{genome}/sigma_candidates.tsv`, `metrics.json`, `report.md`.",
        "",
    ])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Orquestación
# ---------------------------------------------------------------------------
def _discover_genomes(source: str, features_root: Path,
                      promoters_root: Path) -> list[str]:
    found: set[str] = set()
    if source in {"auto", "layers"}:
        if features_root.exists():
            found.update(p.parent.name for p in
                         features_root.glob("*/layer_sigma_features.tsv"))
    if source == "promoters":
        for p in promoters_root.glob("*/sigma_promoters.tsv"):
            found.add(p.parent.name)
    return sorted(found)


def _write_outputs(result: BuildResult, gold_stats: dict[str, Any] | None,
                   evaluation: dict[str, Any],
                   refinement: dict[str, Any],
                   out_root: Path) -> dict[str, Path]:
    out_dir = out_root / result.genome
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates_path = out_dir / "sigma_candidates.tsv"
    metrics_path = out_dir / "metrics.json"
    report_path = out_dir / "report.md"

    _write_tsv(candidates_path, result.candidates, CANDIDATE_COLUMNS)
    metrics = {
        "version": VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "genome": result.genome,
        "source": result.source,
        "config": {
            "top_k": result.top_k,
            "filters": result.filters,
            "refine": result.refine,
            "confidence_thresholds": CONFIDENCE_THRESHOLDS,
            "refinement_rules": REFINEMENT_RULES,
            "coordinate_system": "1-based TSS (layers) / operon anchor (fallback)",
        },
        "inputs": {
            "candidates": _rel(result.source_path),
            "gold": gold_stats.get("path") if gold_stats else None,
        },
        "candidate_stats": result.stats,
        "gold_stats": gold_stats,
        "evaluation": evaluation,
        "refinement": refinement,
        "limitations": LIMITATIONS,
        "outputs": {
            "candidates": _rel(candidates_path),
            "metrics": _rel(metrics_path),
            "report": _rel(report_path),
        },
    }
    metrics_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2),
                            encoding="utf-8")
    report = build_report(result.genome, result, gold_stats, evaluation,
                          refinement)
    report_path.write_text(report, encoding="utf-8")
    return {"candidates": candidates_path, "metrics": metrics_path,
            "report": report_path}


def run_genome(genome: str, args: argparse.Namespace,
               out_root: Path) -> dict[str, Any]:
    source, source_path = resolve_source(
        genome, args.source, args.features_root, args.promoters_root)
    if source == "layers":
        result = build_candidates_from_layers(
            source_path, genome, top_k=args.top_k, refine=args.refine,
            min_z=args.min_z, max_pvalue=args.max_pvalue,
            min_pwm_score=args.min_pwm_score)
    else:
        result = build_candidates_from_promoters(
            source_path, genome, top_k=args.top_k, refine=args.refine,
            min_z=args.min_z, max_pvalue=args.max_pvalue,
            min_pwm_score=args.min_pwm_score)

    gold: dict[tuple[str, int, str], set[str]] = {}
    gold_stats: dict[str, Any] | None = None
    gold_path = args.gold_root / genome / "tss_sigma.tsv"
    if gold_path.exists():
        gold, gold_stats = load_gold(gold_path)

    evaluation = evaluate_candidates(result, gold)
    refinement = evaluate_refinement(result, gold)
    paths = _write_outputs(result, gold_stats, evaluation, refinement, out_root)
    return {"genome": genome, "source": source, "result": result,
            "evaluation": evaluation, "refinement": refinement,
            "paths": paths}


def _print_summary(payload: dict[str, Any]) -> None:
    result: BuildResult = payload["result"]
    evaluation = payload["evaluation"]
    stats = result.stats
    tag = payload["genome"]
    if evaluation.get("status") == "ok":
        ov = evaluation["overall"]
        print(
            f"{tag}: {stats['n_candidates_topk']} candidatos en "
            f"{stats['n_tss_with_candidates']} TSS; "
            f"top1_any={ov['top1_any']}/{ov['n_tss']}; "
            f"top3_any={ov['top3_any']}/{ov['n_tss']}; "
            f"top3_all={ov['top3_all']}/{ov['n_tss']}; "
            f"multi top3_all={evaluation['multi'].get('top3_all', 0)}/"
            f"{evaluation['n_multi_tss']}; "
            f"-> {payload['paths']['candidates']}")
    else:
        print(
            f"{tag}: {stats['n_candidates_topk']} candidatos; "
            f"evaluación={evaluation.get('status')}; "
            f"-> {payload['paths']['candidates']}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Candidatos sigma multietiqueta por TSS y evaluación top-k")
    ap.add_argument("--genome", action="append", default=[],
                    help="Genoma a procesar (repetible).")
    ap.add_argument("--all", action="store_true",
                    help="Procesa todos los genomas con features de capa.")
    ap.add_argument("--source", choices=("auto", "layers", "promoters"),
                    default="auto",
                    help="Fuente de candidatos; auto prefiere layers.")
    ap.add_argument("--top-k", type=int, default=DEFAULT_TOP_K,
                    help=f"Etiquetas por TSS/operón (default {DEFAULT_TOP_K}).")
    ap.add_argument("--min-z", type=float, default=None,
                    help="Umbral mínimo de z para incluir un candidato.")
    ap.add_argument("--max-pvalue", type=float, default=None,
                    help="Umbral máximo de p-valor para incluir un candidato.")
    ap.add_argument("--min-pwm-score", type=float, default=None,
                    help="Umbral mínimo del score PWM crudo.")
    ap.add_argument("--features-root", type=Path,
                    default=DEFAULT_FEATURES_ROOT)
    ap.add_argument("--promoters-root", type=Path,
                    default=DEFAULT_PROMOTERS_ROOT)
    ap.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD_ROOT)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    ap.add_argument("--refine", dest="refine", action="store_true",
                    default=True, help="Aplica reglas σ54/σ38 (default).")
    ap.add_argument("--no-refine", dest="refine", action="store_false",
                    help="Desactiva refined_confidence.")
    return ap.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.top_k < 1:
        print("--top-k debe ser >= 1", file=sys.stderr)
        return 2
    if args.min_z is not None and not math.isfinite(args.min_z):
        print("--min-z inválido", file=sys.stderr)
        return 2
    if args.max_pvalue is not None and not (0.0 <= args.max_pvalue <= 1.0):
        print("--max-pvalue debe estar entre 0 y 1", file=sys.stderr)
        return 2

    if args.genome:
        genomes = list(dict.fromkeys(args.genome))
    elif args.all:
        genomes = _discover_genomes(args.source, args.features_root,
                                    args.promoters_root)
        if not genomes:
            print("No hay genomas con features para --all", file=sys.stderr)
            return 2
    else:
        print("Indicá --genome (repetible) o --all", file=sys.stderr)
        return 2

    out_root: Path = args.out_root
    payloads: list[dict[str, Any]] = []
    for genome in genomes:
        try:
            payload = run_genome(genome, args, out_root)
        except FileNotFoundError as exc:
            print(f"{genome}: {exc}", file=sys.stderr)
            return 2
        except (OSError, ValueError) as exc:
            print(f"{genome}: error: {exc}", file=sys.stderr)
            return 1
        payloads.append(payload)
        _print_summary(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
