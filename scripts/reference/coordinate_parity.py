#!/usr/bin/env python3
"""Paridad por coordenadas: GFF del sistema vs GFF de genomas de referencia.

Compara los features ``CDS``, ``tRNA`` y ``rRNA`` que el sistema anota en
``results/<id>/<id>.gff3`` contra los del GFF de referencia declarado en
``data/reference_genomes.json`` (campo ``annotation_path``).

Reglas de comparación:

- un feature multiparte (varias filas con el mismo ``ID`` en el mismo contig,
  tipo y hebra) se colapsa a un intervalo ``[min(inicio), max(fin)]``;
- el emparejamiento es 1:1 y exige mismo contig, mismo tipo y misma hebra;
- dos features son candidatos si sus intervalos se solapan recíprocamente
  (fracción de solape mínima configurable) o si ambos bordes difieren como
  mucho la tolerancia configurable (5 pb por defecto);
- la asignación 1:1 es de cardinalidad máxima y prefiere, en este orden,
  bordes exactos, menor diferencia de bordes y mayor solape.

Métricas por genoma y clase: TP, FP, FN, precisión, recall, F1, coincidencia
exacta de bordes y diferencia media/mediana de inicio y fin (sistema -
referencia).

No evalúa función, producto ni identidad de secuencia. No usa red, GPU ni
binarios externos y no modifica ``results/``, los GFF ni el manifiesto.

Salidas: ``reports/reference_coordinate_parity.tsv`` y
``reports/reference_coordinate_parity.md``.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import re
import statistics
import sys
from collections import Counter, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
ELIGIBLE_STATUSES = ("annotated", "available")
CLASSES = ("CDS", "tRNA", "rRNA")
REPORT_CLASSES = CLASSES + ("ALL",)
DEFAULT_TOLERANCE = 5
DEFAULT_MIN_RECIPROCAL_OVERLAP = 0.5
_ID_RE = re.compile(r"(?:^|;)\s*ID=([^;]+)")


# ---------------------------------------------------------------------------
# Modelo de datos
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Feature:
    """Un feature colapsado a un único intervalo genómico."""

    seqid: str
    kind: str
    start: int
    end: int
    strand: str
    feature_id: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class GffData:
    """Resultado de parsear un GFF3."""

    features: list[Feature]
    raw_counts: dict[str, int]
    invalid_lines: int
    rows: int

    def by_class(self) -> dict[str, list[Feature]]:
        grouped: dict[str, list[Feature]] = {kind: [] for kind in CLASSES}
        for feature in self.features:
            grouped[feature.kind].append(feature)
        return grouped


# ---------------------------------------------------------------------------
# Entradas
# ---------------------------------------------------------------------------
def open_text(path: Path | str):
    """Abre texto plano o gzip según la extensión."""
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def load_manifest(path: Path | str) -> list[dict[str, Any]]:
    """Lee ``reference_genomes.json`` y devuelve sus registros."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no existe el manifiesto: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"manifiesto ilegible ({path}): {exc}") from exc
    if not isinstance(data, list):
        raise ValueError(f"el manifiesto debe ser una lista: {path}")
    records = [rec for rec in data if isinstance(rec, dict)]
    if not records:
        raise ValueError(f"el manifiesto no contiene registros: {path}")
    return records


def _resolve(root: Path, raw: Any) -> Path | None:
    if raw is None or str(raw).strip() == "":
        return None
    path = Path(str(raw))
    return path if path.is_absolute() else root / path


def _feature_id(attributes: str, fallback: str) -> str:
    match = _ID_RE.search(attributes)
    if match:
        value = match.group(1).strip()
        if value:
            return value
    return fallback


def parse_gff(path: Path | str) -> GffData:
    """Parsea un GFF3 y colapsa features multiparte por ``ID``.

    Agrupa por ``(seqid, tipo, hebra, ID)`` para no fusionar identificadores
    repetidos en contigs o hebras distintas. Las filas sin ``ID`` reciben un
    identificador único basado en su número de línea. Las líneas malformadas se
    cuentan en ``invalid_lines`` en lugar de abortar el archivo.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no existe el GFF: {path}")

    grouped: dict[tuple[str, str, str, str], Feature] = {}
    raw_counts: Counter[str] = Counter()
    invalid_lines = 0
    rows = 0
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                invalid_lines += 1
                continue
            seqid = fields[0].strip()
            kind = fields[2].strip()
            if not seqid or not kind:
                invalid_lines += 1
                continue
            try:
                raw_start = int(fields[3])
                raw_end = int(fields[4])
            except (TypeError, ValueError):
                invalid_lines += 1
                continue
            if raw_start < 1 or raw_end < 1:
                invalid_lines += 1
                continue
            # Algunos GFF del sistema emiten las features de hebra negativa con
            # inicio > fin. El intervalo genómico es el mismo; se normaliza a
            # [menor, mayor] para comparar coordenadas.
            start, end = (
                (raw_start, raw_end)
                if raw_start <= raw_end
                else (raw_end, raw_start)
            )
            rows += 1
            raw_counts[kind] += 1
            if kind not in CLASSES:
                continue
            strand = fields[6].strip() or "."
            fid = _feature_id(fields[8], f"__row_{lineno}")
            key = (seqid, kind, strand, fid)
            previous = grouped.get(key)
            if previous is None:
                grouped[key] = Feature(seqid, kind, start, end, strand, fid)
            else:
                grouped[key] = Feature(
                    seqid,
                    kind,
                    min(previous.start, start),
                    max(previous.end, end),
                    strand,
                    fid,
                )

    features = sorted(
        grouped.values(),
        key=lambda f: (f.seqid, f.kind, f.start, f.end, f.strand, f.feature_id),
    )
    return GffData(
        features=features,
        raw_counts=dict(raw_counts),
        invalid_lines=invalid_lines,
        rows=rows,
    )


def _system_gff_candidates(record: dict[str, Any], root: Path) -> list[Path]:
    """Rutas convencionales del GFF del sistema, en orden de preferencia."""
    genome = str(record.get("id") or "").strip()
    base_raw = str(record.get("system_result") or "").strip()
    base = _resolve(root, base_raw) if base_raw else root / "results" / genome
    if base is None:
        base = root / "results" / genome
    if base.suffix.lower() in (".gff", ".gff3", ".gz") and base.is_file():
        return [base]
    names = (
        f"{genome}.gff3",
        f"{genome}.gff3.gz",
        f"{genome}.gff",
        "annotation.gff3",
        "annotation.gff3.gz",
    )
    return [base / name for name in names]


def system_gff_path(record: dict[str, Any], root: Path | str) -> tuple[Path | None, str]:
    """Devuelve el GFF del sistema existente y la ruta declarada/esperada."""
    root = Path(root)
    candidates = _system_gff_candidates(record, root)
    expected = candidates[0]
    for candidate in candidates:
        if candidate.is_file():
            display = (
                str(candidate.relative_to(root))
                if candidate.is_relative_to(root)
                else str(candidate)
            )
            return candidate, display
    display = (
        str(expected.relative_to(root))
        if expected.is_relative_to(root)
        else str(expected)
    )
    return None, display


# ---------------------------------------------------------------------------
# Emparejamiento 1:1
# ---------------------------------------------------------------------------
def _overlap(left: Feature, right: Feature) -> int:
    return max(0, min(left.end, right.end) - max(left.start, right.start) + 1)


def _pair_quality(
    system: Feature,
    reference: Feature,
    tolerance: int,
    min_reciprocal_overlap: float,
) -> tuple[int, int, int, int, int] | None:
    """Devuelve la clave de calidad del par, o ``None`` si no es válido.

    Un par es válido si cumple al menos una condición:

    1. solape recíproco: la fracción de solape respecto de cada feature es
       mayor o igual a ``min_reciprocal_overlap``;
    2. bordes próximos: ``|Δinicio| <= tolerancia`` y
       ``|Δfin| <= tolerancia``.

    La tolerancia rescata features pequeñas o desplazadas donde el solape
    recíproco por sí solo no alcanzaría el umbral. No relaja la exigencia de
    mismo contig, tipo y hebra.
    """
    overlap = _overlap(system, reference)
    delta_start = system.start - reference.start
    delta_end = system.end - reference.end
    within_edges = (
        abs(delta_start) <= tolerance and abs(delta_end) <= tolerance
    )
    reciprocal = (
        overlap > 0
        and overlap / system.length >= min_reciprocal_overlap
        and overlap / reference.length >= min_reciprocal_overlap
    )
    if not within_edges and not reciprocal:
        return None
    exact = 1 if delta_start == 0 and delta_end == 0 else 0
    max_edge = max(abs(delta_start), abs(delta_end))
    return (
        0 if exact else 1,
        max_edge,
        abs(delta_start) + abs(delta_end),
        -overlap,
        0,
    )


def _quality_sort_key(entry: tuple[tuple[int, int, int, int, int], int]) -> tuple:
    quality, ref_idx = entry
    return (quality, ref_idx)


def match_features(
    system_features: Sequence[Feature],
    reference_features: Sequence[Feature],
    tolerance: int = DEFAULT_TOLERANCE,
    min_reciprocal_overlap: float = DEFAULT_MIN_RECIPROCAL_OVERLAP,
) -> list[tuple[Feature, Feature]]:
    """Empareja 1:1 features del sistema y de referencia.

    Exige mismo contig, tipo y hebra. La asignación final es de cardinalidad
    máxima: se parte de un emparejamiento greedy que respeta el orden de
    calidad y se completa con caminos de aumento. Así se conservan primero los
    pares con bordes exactos sin sacrificar el número total de pares posibles.
    """
    if tolerance < 0:
        raise ValueError("la tolerancia debe ser >= 0")
    if not 0.0 <= min_reciprocal_overlap <= 1.0:
        raise ValueError("min_reciprocal_overlap debe estar entre 0 y 1")

    # Índice de referencia por la clave que debe coincidir exactamente. Se
    # guarda el índice global de cada feature para que el emparejamiento no
    # confunda posiciones locales de un grupo con posiciones del vector global.
    index: dict[tuple[str, str, str], list[tuple[int, Feature]]] = {}
    for ref_idx, reference in enumerate(reference_features):
        index.setdefault(
            (reference.seqid, reference.kind, reference.strand), []
        ).append((ref_idx, reference))
    for values in index.values():
        values.sort(key=lambda item: (item[1].start, item[1].end, item[1].feature_id))
    starts = {
        key: [feature.start for _ref_idx, feature in values]
        for key, values in index.items()
    }
    # ``prefix_max_end[i]`` es el mayor ``end`` entre las features 0..i. Como
    # es monótono, permite saltar por bisección todos los features que
    # terminan antes de la ventana de consulta y evita un barrido O(n²).
    prefix_max_end = {}
    for key, values in index.items():
        running = 0
        maxima: list[int] = []
        for _ref_idx, feature in values:
            running = max(running, feature.end)
            maxima.append(running)
        prefix_max_end[key] = maxima

    adjacency: list[list[tuple[tuple[int, int, int, int, int], int]]] = []
    for system in system_features:
        key = (system.seqid, system.kind, system.strand)
        candidates = index.get(key, [])
        entries: list[tuple[tuple[int, int, int, int, int], int]] = []
        if candidates:
            starts_for_key = starts[key]
            limit = system.end + tolerance
            hi = bisect.bisect_right(starts_for_key, limit)
            left = bisect.bisect_left(
                prefix_max_end[key], system.start - tolerance
            )
            for position in range(left, hi):
                ref_idx, reference = candidates[position]
                if reference.end < system.start - tolerance:
                    continue
                quality = _pair_quality(
                    system, reference, tolerance, min_reciprocal_overlap
                )
                if quality is not None:
                    entries.append((quality, ref_idx))
        entries.sort(key=_quality_sort_key)
        adjacency.append(entries)

    matched_system: dict[int, int] = {}
    matched_reference: dict[int, int] = {}

    # 1) Greedy global por calidad: fija primero los mejores pares y conserva
    #    los emparejamientos exactos.
    all_pairs: list[tuple[tuple[int, int, int, int, int], int, int]] = []
    for sys_idx, entries in enumerate(adjacency):
        for quality, ref_idx in entries:
            all_pairs.append((quality, sys_idx, ref_idx))
    all_pairs.sort(key=lambda item: (item[0], item[1], item[2]))
    for _quality, sys_idx, ref_idx in all_pairs:
        if sys_idx in matched_system or ref_idx in matched_reference:
            continue
        matched_system[sys_idx] = ref_idx
        matched_reference[ref_idx] = sys_idx

    # 2) Caminos de aumento hasta alcanzar cardinalidad máxima.
    def augment(start: int) -> bool:
        parent_system: dict[int, int | None] = {start: None}
        ref_source: dict[int, int] = {}
        seen_system = {start}
        queue: deque[int] = deque([start])
        found: int | None = None

        while queue and found is None:
            sys_idx = queue.popleft()
            for _quality, ref_idx in adjacency[sys_idx]:
                if ref_idx in ref_source:
                    continue
                ref_source[ref_idx] = sys_idx
                if ref_idx not in matched_reference:
                    found = ref_idx
                    break
                next_system = matched_reference[ref_idx]
                if next_system not in seen_system:
                    seen_system.add(next_system)
                    parent_system[next_system] = ref_idx
                    queue.append(next_system)

        if found is None:
            return False

        current_ref: int | None = found
        while current_ref is not None:
            current_system = ref_source[current_ref]
            previous_ref = parent_system[current_system]
            matched_reference[current_ref] = current_system
            matched_system[current_system] = current_ref
            current_ref = previous_ref
        return True

    while True:
        changed = False
        for sys_idx in range(len(system_features)):
            if sys_idx in matched_system:
                continue
            if augment(sys_idx):
                changed = True
        if not changed:
            break

    matches = [
        (system_features[sys_idx], reference_features[ref_idx])
        for sys_idx, ref_idx in sorted(matched_system.items())
    ]
    return matches


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
def _ratio(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return round(numerator / denominator, 4)


def _percent_ratio(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return round(100.0 * numerator / denominator, 2)


def _stats(values: Sequence[int]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    return round(statistics.fmean(values), 2), round(statistics.median(values), 2)


def _empty_metrics() -> dict[str, Any]:
    return {
        "state": "not_available",
        "system_features": None,
        "reference_features": None,
        "tp": None,
        "fp": None,
        "fn": None,
        "precision": None,
        "recall": None,
        "f1": None,
        "exact_edges": None,
        "exact_pct": None,
        "mean_start_diff": None,
        "median_start_diff": None,
        "mean_end_diff": None,
        "median_end_diff": None,
        "max_abs_edge_diff": None,
        "start_diffs": [],
        "end_diffs": [],
    }


def compare_feature_sets(
    system_features: Sequence[Feature],
    reference_features: Sequence[Feature],
    *,
    tolerance: int = DEFAULT_TOLERANCE,
    min_reciprocal_overlap: float = DEFAULT_MIN_RECIPROCAL_OVERLAP,
    label: str = "ALL",
) -> dict[str, Any]:
    """Calcula TP/FP/FN y estadísticas de bordes para un conjunto de features."""
    matches = match_features(
        system_features,
        reference_features,
        tolerance=tolerance,
        min_reciprocal_overlap=min_reciprocal_overlap,
    )
    tp = len(matches)
    fp = len(system_features) - tp
    fn = len(reference_features) - tp
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    if precision is None or recall is None:
        f1 = None
    elif precision + recall <= 0:
        f1 = 0.0
    else:
        f1 = round(2 * precision * recall / (precision + recall), 4)

    start_diffs = [system.start - reference.start for system, reference in matches]
    end_diffs = [system.end - reference.end for system, reference in matches]
    exact_edges = sum(
        1 for system, reference in matches
        if system.start == reference.start and system.end == reference.end
    )
    mean_start, median_start = _stats(start_diffs)
    mean_end, median_end = _stats(end_diffs)
    max_abs = max((abs(value) for value in start_diffs + end_diffs), default=None)

    if not system_features and not reference_features:
        state = "not_available"
        reason = "la clase objetivo está ausente en ambas anotaciones"
    elif fp == 0 and fn == 0 and tp > 0 and (max_abs is None or max_abs <= tolerance):
        state = "ok"
        reason = "paridad 1:1 dentro de la tolerancia de bordes"
    else:
        state = "partial"
        parts = [f"TP={tp}", f"FP={fp}", f"FN={fn}"]
        if max_abs is not None:
            parts.append(f"deriva máxima={max_abs} pb")
        if tp == 0:
            parts.append("sin pares 1:1")
        reason = "; ".join(parts)

    return {
        "label": label,
        "state": state,
        "system_features": len(system_features),
        "reference_features": len(reference_features),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "exact_edges": exact_edges,
        "exact_pct": _percent_ratio(exact_edges, tp),
        "mean_start_diff": mean_start,
        "median_start_diff": median_start,
        "mean_end_diff": mean_end,
        "median_end_diff": median_end,
        "max_abs_edge_diff": max_abs,
        "start_diffs": start_diffs,
        "end_diffs": end_diffs,
        "reason": reason,
    }


def _genome_state(classes: dict[str, dict[str, Any]]) -> str:
    present = [
        metrics
        for metrics in classes.values()
        if metrics["system_features"] is not None
        and (
            metrics["system_features"] > 0
            or metrics["reference_features"] > 0
        )
    ]
    if not present:
        return "not_available"
    if all(metrics["state"] == "ok" for metrics in present):
        return "ok"
    return "partial"


# ---------------------------------------------------------------------------
# Comparación por genoma
# ---------------------------------------------------------------------------
def _not_available_row(
    record: dict[str, Any],
    *,
    system_gff: str,
    reference_gff: str,
    reason: str,
    tolerance: int,
    min_reciprocal_overlap: float,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "genome": str(record.get("id") or ""),
        "species": str(record.get("species") or ""),
        "strain": str(record.get("strain") or ""),
        "manifest_status": str(record.get("status") or ""),
        "state": "not_available",
        "reason": reason,
        "system_gff": system_gff,
        "reference_gff": reference_gff,
        "tolerance": tolerance,
        "min_reciprocal_overlap": min_reciprocal_overlap,
        "classes": {kind: _empty_metrics() for kind in CLASSES},
        "all": _empty_metrics(),
    }
    for metrics in row["classes"].values():
        metrics["reason"] = reason
    row["all"]["reason"] = reason
    return row


def compare_record(
    record: dict[str, Any],
    root: Path | str,
    *,
    tolerance: int = DEFAULT_TOLERANCE,
    min_reciprocal_overlap: float = DEFAULT_MIN_RECIPROCAL_OVERLAP,
) -> dict[str, Any]:
    """Compara un genoma del manifiesto a nivel de coordenadas."""
    root = Path(root)
    genome = str(record.get("id") or "").strip()
    annotation_raw = str(record.get("annotation_path") or "").strip()
    annotation_path = _resolve(root, annotation_raw)
    system_path, system_display = system_gff_path(record, root)
    reference_display = annotation_raw

    if system_path is None:
        return _not_available_row(
            record,
            system_gff=system_display,
            reference_gff=reference_display,
            reason=f"falta el GFF del sistema: {system_display}",
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )
    if annotation_path is None:
        return _not_available_row(
            record,
            system_gff=system_display,
            reference_gff=reference_display,
            reason="el manifiesto no declara annotation_path",
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )
    if not annotation_path.is_file():
        return _not_available_row(
            record,
            system_gff=system_display,
            reference_gff=reference_display,
            reason=f"falta el GFF de referencia: {annotation_raw}",
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )

    try:
        system_gff = parse_gff(system_path)
        reference_gff = parse_gff(annotation_path)
    except (OSError, EOFError, ValueError) as exc:
        return _not_available_row(
            record,
            system_gff=system_display,
            reference_gff=reference_display,
            reason=f"GFF ilegible: {exc}",
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )

    system_by_class = system_gff.by_class()
    reference_by_class = reference_gff.by_class()
    if not reference_gff.features:
        return _not_available_row(
            record,
            system_gff=system_display,
            reference_gff=reference_display,
            reason="el GFF de referencia no contiene CDS/tRNA/rRNA",
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )

    classes: dict[str, dict[str, Any]] = {}
    for kind in CLASSES:
        classes[kind] = compare_feature_sets(
            system_by_class[kind],
            reference_by_class[kind],
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
            label=kind,
        )
    all_metrics = compare_feature_sets(
        system_gff.features,
        reference_gff.features,
        tolerance=tolerance,
        min_reciprocal_overlap=min_reciprocal_overlap,
        label="ALL",
    )

    state = _genome_state(classes)
    if state == "ok":
        max_abs = all_metrics["max_abs_edge_diff"]
        reason = (
            "paridad 1:1 en las clases con datos; "
            f"deriva máxima {max_abs} pb ≤ {tolerance} pb"
        )
    elif state == "partial":
        reason = all_metrics["reason"]
    else:
        reason = "ninguna clase CDS/tRNA/rRNA es comparable"
    all_metrics["reason"] = reason

    return {
        "genome": genome,
        "species": str(record.get("species") or ""),
        "strain": str(record.get("strain") or ""),
        "manifest_status": str(record.get("status") or ""),
        "state": state,
        "reason": reason,
        "system_gff": system_display,
        "reference_gff": reference_display,
        "tolerance": tolerance,
        "min_reciprocal_overlap": min_reciprocal_overlap,
        "classes": classes,
        "all": all_metrics,
    }


def _metrics_for(row: dict[str, Any], kind: str) -> dict[str, Any]:
    return row["all"] if kind == "ALL" else row["classes"][kind]


def _aggregate(rows: Sequence[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Agrega métricas por clase sobre los genomas comparables."""
    aggregate: dict[str, dict[str, Any]] = {}
    for kind in REPORT_CLASSES:
        metrics_rows = [
            _metrics_for(row, kind)
            for row in rows
            if _metrics_for(row, kind)["system_features"] is not None
        ]
        if not metrics_rows:
            aggregate[kind] = {
                "system_features": 0,
                "reference_features": 0,
                "tp": 0,
                "fp": 0,
                "fn": 0,
                "precision": None,
                "recall": None,
                "f1": None,
                "exact_edges": 0,
                "exact_pct": None,
                "mean_start_diff": None,
                "median_start_diff": None,
                "mean_end_diff": None,
                "median_end_diff": None,
                "max_abs_edge_diff": None,
            }
            continue
        tp = sum(item["tp"] for item in metrics_rows)
        fp = sum(item["fp"] for item in metrics_rows)
        fn = sum(item["fn"] for item in metrics_rows)
        precision = _ratio(tp, tp + fp)
        recall = _ratio(tp, tp + fn)
        if precision is None or recall is None:
            f1 = None
        elif precision + recall <= 0:
            f1 = 0.0
        else:
            f1 = round(2 * precision * recall / (precision + recall), 4)
        start_diffs = [
            value for item in metrics_rows for value in item["start_diffs"]
        ]
        end_diffs = [
            value for item in metrics_rows for value in item["end_diffs"]
        ]
        mean_start, median_start = _stats(start_diffs)
        mean_end, median_end = _stats(end_diffs)
        exact_edges = sum(item["exact_edges"] for item in metrics_rows)
        aggregate[kind] = {
            "system_features": sum(item["system_features"] for item in metrics_rows),
            "reference_features": sum(
                item["reference_features"] for item in metrics_rows
            ),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "exact_edges": exact_edges,
            "exact_pct": _percent_ratio(exact_edges, tp),
            "mean_start_diff": mean_start,
            "median_start_diff": median_start,
            "mean_end_diff": mean_end,
            "median_end_diff": median_end,
            "max_abs_edge_diff": max(
                (abs(value) for value in start_diffs + end_diffs), default=None
            ),
        }
    return aggregate


def build_report(
    manifest_path: Path | str,
    root: Path | str,
    *,
    generated_at: str | None = None,
    tolerance: int = DEFAULT_TOLERANCE,
    min_reciprocal_overlap: float = DEFAULT_MIN_RECIPROCAL_OVERLAP,
    statuses: Sequence[str] = ELIGIBLE_STATUSES,
) -> dict[str, Any]:
    """Construye el informe completo desde el manifiesto y los GFF."""
    if tolerance < 0:
        raise ValueError("la tolerancia debe ser >= 0")
    if not 0.0 <= min_reciprocal_overlap <= 1.0:
        raise ValueError("min_reciprocal_overlap debe estar entre 0 y 1")
    manifest_path = Path(manifest_path)
    root = Path(root)
    records = load_manifest(manifest_path)
    eligible = tuple(status.strip().lower() for status in statuses)

    selected: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for record in records:
        status = str(record.get("status") or "").strip().lower()
        if status not in eligible:
            skipped.append(
                {
                    "genome": str(record.get("id") or ""),
                    "status": status,
                    "reason": f"status '{status}' fuera de {', '.join(eligible)}",
                }
            )
        else:
            selected.append(record)

    rows = [
        compare_record(
            record,
            root,
            tolerance=tolerance,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )
        for record in selected
    ]
    state_counts = Counter(row["state"] for row in rows)
    comparable = sum(1 for row in rows if row["state"] != "not_available")
    return {
        "generated_at": generated_at
        or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifest": str(manifest_path),
        "root": str(root),
        "tolerance": tolerance,
        "min_reciprocal_overlap": min_reciprocal_overlap,
        "rows": rows,
        "skipped": skipped,
        "summary": {
            "candidates": len(selected),
            "comparable": comparable,
            "ok": state_counts.get("ok", 0),
            "partial": state_counts.get("partial", 0),
            "not_available": state_counts.get("not_available", 0),
            "skipped": len(skipped),
        },
        "aggregate": _aggregate(rows),
    }


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------
TSV_COLUMNS = (
    "genome",
    "species",
    "strain",
    "manifest_status",
    "class",
    "state",
    "system_features",
    "reference_features",
    "tp",
    "fp",
    "fn",
    "precision",
    "recall",
    "f1",
    "exact_edges",
    "exact_pct",
    "mean_start_diff",
    "median_start_diff",
    "mean_end_diff",
    "median_end_diff",
    "max_abs_edge_diff",
    "tolerance",
    "min_reciprocal_overlap",
    "system_gff",
    "reference_gff",
    "reason",
)


def _cell(value: Any) -> Any:
    return "" if value is None else value


def _rate_cell(value: float | None) -> str:
    return "" if value is None else f"{value:.4f}"


def _diff_cell(value: float | None) -> str:
    return "" if value is None else f"{value:.2f}"


def _pct_cell(value: float | None) -> str:
    return "" if value is None else f"{value:.2f}"


def write_tsv(rows: Iterable[dict[str, Any]], path: Path | str) -> None:
    """Escribe la tabla TSV máquina-legible (una fila por genoma y clase)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for row in rows:
            for kind in REPORT_CLASSES:
                metrics = _metrics_for(row, kind)
                writer.writerow(
                    [
                        row["genome"],
                        row["species"],
                        row["strain"],
                        row["manifest_status"],
                        kind,
                        metrics["state"],
                        _cell(metrics["system_features"]),
                        _cell(metrics["reference_features"]),
                        _cell(metrics["tp"]),
                        _cell(metrics["fp"]),
                        _cell(metrics["fn"]),
                        _rate_cell(metrics["precision"]),
                        _rate_cell(metrics["recall"]),
                        _rate_cell(metrics["f1"]),
                        _cell(metrics["exact_edges"]),
                        _pct_cell(metrics["exact_pct"]),
                        _diff_cell(metrics["mean_start_diff"]),
                        _diff_cell(metrics["median_start_diff"]),
                        _diff_cell(metrics["mean_end_diff"]),
                        _diff_cell(metrics["median_end_diff"]),
                        _cell(metrics["max_abs_edge_diff"]),
                        row["tolerance"],
                        row["min_reciprocal_overlap"],
                        row["system_gff"],
                        row["reference_gff"],
                        metrics["reason"],
                    ]
                )


def _md(value: Any) -> str:
    if value is None or value == "":
        return "—"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _md_rate(value: float | None) -> str:
    return "—" if value is None else f"{value:.4f}"


def _md_diff(value: float | None) -> str:
    return "—" if value is None else f"{value:+.2f}"


def _md_pct(value: float | None) -> str:
    return "—" if value is None else f"{value:.2f}%"


def render_markdown(result: dict[str, Any]) -> str:
    """Renderiza el informe legible en Markdown."""
    summary = result["summary"]
    tolerance = result["tolerance"]
    min_reciprocal = result["min_reciprocal_overlap"]
    rows = result["rows"]
    lines: list[str] = [
        "# Paridad por coordenadas: sistema vs GFF de referencia",
        "",
        f"**Generado:** `{result['generated_at']}`  ",
        f"**Manifiesto:** `{result['manifest']}`  ",
        f"**Raíz de datos:** `{result['root']}`  ",
        f"**Tolerancia de bordes:** {tolerance} pb  ",
        f"**Solape recíproco mínimo:** {min_reciprocal:.2f}",
        "",
        "## Alcance",
        "",
        "Compara los features **CDS, tRNA y rRNA** del GFF del sistema "
        "(`results/<id>/<id>.gff3`) contra el GFF de referencia "
        "(`annotation_path`) para cada registro con status `annotated` o "
        "`available` del manifiesto.",
        "",
        "El emparejamiento es **1:1** y exige mismo contig, mismo tipo y misma "
        "hebra. Dos features son candidatos si se solapan recíprocamente o si "
        f"sus dos bordes difieren como mucho **{tolerance} pb**. Si falta "
        "cualquiera de los dos GFF, el genoma queda como `not_available`; no "
        "se imputan coordenadas ni se completa con otros archivos.",
        "",
        "Esta versión **no evalúa función, producto, identidad de secuencia ni "
        "calidad biológica** de las predicciones. Tampoco modifica `results/`, "
        "los GFF ni el manifiesto.",
        "",
        "## Cobertura",
        "",
        f"- Registros con status elegible: **{summary['candidates']}**",
        f"- Comparables (ambos GFF legibles): **{summary['comparable']}** "
        f"(ok **{summary['ok']}**, partial **{summary['partial']}**)",
        f"- No disponibles: **{summary['not_available']}**",
        f"- Omitidos por status: **{summary['skipped']}**",
        "",
    ]

    unavailable = [row for row in rows if row["state"] == "not_available"]
    if unavailable:
        lines.extend(
            [
                "### Genomas no disponibles",
                "",
                "| genoma | status | motivo |",
                "|---|---|---|",
            ]
        )
        for row in unavailable:
            lines.append(
                f"| `{_md(row['genome'])}` | {_md(row['manifest_status'])} | "
                f"{_md(row['reason'])} |"
            )
        lines.append("")

    lines.extend(
        [
            "## Resultados por genoma (todas las clases)",
            "",
            "TP/FP/FN se calculan sobre CDS + tRNA + rRNA. Precisión, recall "
            "y F1 son fracciones (0–1); `% exactos` es un porcentaje.",
            "",
            "| genoma | estado | sistema | referencia | TP | FP | FN | "
            "precisión | recall | F1 |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        metrics = row["all"]
        cells = [
            f"`{_md(row['genome'])}`",
            _md(row["state"]),
            _md(metrics["system_features"]),
            _md(metrics["reference_features"]),
            _md(metrics["tp"]),
            _md(metrics["fp"]),
            _md(metrics["fn"]),
            _md_rate(metrics["precision"]),
            _md_rate(metrics["recall"]),
            _md_rate(metrics["f1"]),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    lines.extend(
        [
            "### Bordes de los pares emparejados",
            "",
            "Las diferencias son **sistema − referencia** en pb (GFF 1-based, "
            "inclusivo). `exactos` cuenta pares con inicio y fin idénticos.",
            "",
            "| genoma | pares | exactos | % exactos | Δinicio medio | "
            "Δinicio mediano | Δfin medio | Δfin mediano | deriva máxima |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        metrics = row["all"]
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(metrics['tp'])} | "
            f"{_md(metrics['exact_edges'])} | {_md_pct(metrics['exact_pct'])} | "
            f"{_md_diff(metrics['mean_start_diff'])} | "
            f"{_md_diff(metrics['median_start_diff'])} | "
            f"{_md_diff(metrics['mean_end_diff'])} | "
            f"{_md_diff(metrics['median_end_diff'])} | "
            f"{_md(metrics['max_abs_edge_diff'])} |"
        )
    lines.append("")

    lines.extend(
        [
            "## Resultados por clase",
            "",
            "| genoma | clase | estado | sistema | referencia | TP | FP | FN | "
            "precisión | recall | F1 | exactos |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        for kind in CLASSES:
            metrics = _metrics_for(row, kind)
            lines.append(
                f"| `{_md(row['genome'])}` | {kind} | {_md(metrics['state'])} | "
                f"{_md(metrics['system_features'])} | "
                f"{_md(metrics['reference_features'])} | {_md(metrics['tp'])} | "
                f"{_md(metrics['fp'])} | {_md(metrics['fn'])} | "
                f"{_md_rate(metrics['precision'])} | "
                f"{_md_rate(metrics['recall'])} | {_md_rate(metrics['f1'])} | "
                f"{_md(metrics['exact_edges'])} |"
            )
    lines.append("")

    aggregate = result["aggregate"]
    lines.extend(
        [
            "## Agregado (micro)",
            "",
            "Suma TP/FP/FN de todos los genomas comparables para cada clase.",
            "",
            "| clase | TP | FP | FN | precisión | recall | F1 | exactos | "
            "% exactos | Δinicio medio | Δfin medio |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for kind in REPORT_CLASSES:
        agg = aggregate[kind]
        lines.append(
            f"| {kind} | {agg['tp']} | {agg['fp']} | {agg['fn']} | "
            f"{_md_rate(agg['precision'])} | {_md_rate(agg['recall'])} | "
            f"{_md_rate(agg['f1'])} | {agg['exact_edges']} | "
            f"{_md_pct(agg['exact_pct'])} | "
            f"{_md_diff(agg['mean_start_diff'])} | "
            f"{_md_diff(agg['mean_end_diff'])} |"
        )
    lines.append("")

    lines.extend(
        [
            "## Estado y motivos",
            "",
            "| genoma | estado | motivo |",
            "|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['state'])} | "
            f"{_md(row['reason'])} |"
        )
    lines.append("")

    lines.extend(
        [
            "## Entradas",
            "",
            "| genoma | especie | cepa | status | GFF del sistema | "
            "GFF de referencia |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['species'])} | "
            f"{_md(row['strain'])} | {_md(row['manifest_status'])} | "
            f"`{_md(row['system_gff'])}` | `{_md(row['reference_gff'])}` |"
        )
    lines.append("")

    lines.extend(
        [
            "## Criterio de estado",
            "",
            "- **ok**: todas las clases con datos tienen `FP=0`, `FN=0`, al menos "
            "un par 1:1 y la deriva máxima de bordes no supera la tolerancia.",
            "- **partial**: la comparación fue posible, pero hay FP/FN o algún "
            "par con deriva de bordes mayor que la tolerancia.",
            "- **not_available**: falta o es ilegible el GFF del sistema o el de "
            "referencia, el manifiesto no declara `annotation_path`, o la "
            "referencia no contiene CDS/tRNA/rRNA. No se imputa nada.",
            "",
            "## Límites",
            "",
            "- **Depende del GFF del sistema**: se compara lo que el pipeline "
            "escribió en `results/<id>/<id>.gff3`, no el resumen de conteos.",
            "- **Solo tres clases**: `ncRNA`, `tmRNA`, `pseudogene`, `gene`, "
            "`exon` y otros tipos quedan fuera aunque estén en el GFF.",
            "- **Emparejamiento geométrico**: un par 1:1 no demuestra que sea el "
            "mismo gen ortólogo ni que la secuencia traducida coincida.",
            "- **Contig y hebra exactos**: diferencias de nomenclatura de contig "
            "o hebras no normalizadas producen FP/FN masivos, no un empate "
            "silencioso.",
            "- **Multiparte**: varias filas con el mismo `ID` en el mismo contig, "
            "tipo y hebra se colapsan a `[min(inicio), max(fin)]`; no se "
            "evalúan intrones ni fases.",
            "- **El gold no es verdad absoluta**: hereda criterios y posibles "
            "errores de la anotación RefSeq/NCBI o de la fuente usada.",
            "- **Tolerancia**: `ok` es paridad geométrica dentro de la "
            "tolerancia, no identidad exacta de bordes; el porcentaje de pares "
            "exactos se reporta aparte.",
            "- **Sin causalidad**: un FP o FN no indica por sí solo un error del "
            "sistema ni de la referencia; no se usa para entrenamiento.",
            "",
            "## Reproducir",
            "",
            "```bash",
            ".venv/bin/python scripts/reference/coordinate_parity.py",
            "```",
            "",
            "Opciones: `--root`, `--manifest`, `--tsv`, `--md`, "
            "`--tolerance`, `--min-reciprocal-overlap`.",
            "",
        ]
    )
    return "\n".join(lines)


def write_markdown(result: dict[str, Any], path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_markdown(result), encoding="utf-8")


def write_reports(
    result: dict[str, Any], tsv_path: Path | str, md_path: Path | str
) -> tuple[Path, Path]:
    """Escribe ambas salidas y devuelve sus rutas."""
    tsv_path, md_path = Path(tsv_path), Path(md_path)
    write_tsv(result["rows"], tsv_path)
    write_markdown(result, md_path)
    return tsv_path, md_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compara coordenadas CDS/tRNA/rRNA del GFF del sistema contra los "
            "GFF de reference_genomes.json."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="raíz del proyecto (default: la del script)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="manifiesto (default: <root>/data/reference_genomes.json)",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=None,
        help=(
            "salida TSV "
            "(default: <root>/reports/reference_coordinate_parity.tsv)"
        ),
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=None,
        help=(
            "salida Markdown "
            "(default: <root>/reports/reference_coordinate_parity.md)"
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=int,
        default=DEFAULT_TOLERANCE,
        help="tolerancia de bordes en pb (default: 5)",
    )
    parser.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=DEFAULT_MIN_RECIPROCAL_OVERLAP,
        help="fracción mínima de solape recíproco (default: 0.5)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.tolerance < 0:
        print("ERROR: --tolerance debe ser >= 0", file=sys.stderr)
        return 2
    if not 0.0 <= args.min_reciprocal_overlap <= 1.0:
        print(
            "ERROR: --min-reciprocal-overlap debe estar entre 0 y 1",
            file=sys.stderr,
        )
        return 2

    root = (args.root or ROOT).resolve()
    manifest = (args.manifest or root / "data" / "reference_genomes.json").resolve()
    tsv_path = (
        args.tsv or root / "reports" / "reference_coordinate_parity.tsv"
    ).resolve()
    md_path = (
        args.md or root / "reports" / "reference_coordinate_parity.md"
    ).resolve()

    try:
        result = build_report(
            manifest,
            root,
            tolerance=args.tolerance,
            min_reciprocal_overlap=args.min_reciprocal_overlap,
        )
    except (OSError, EOFError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    write_reports(result, tsv_path, md_path)
    summary = result["summary"]
    print(
        f"{summary['candidates']} registros elegibles, "
        f"{summary['comparable']} comparables: "
        f"{summary['ok']} ok, {summary['partial']} partial, "
        f"{summary['not_available']} not_available; "
        f"{summary['skipped']} omitidos por status"
    )
    print(f"TSV: {tsv_path}")
    print(f"MD:  {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
