"""Auditoría de CDS solapados y de su marco de lectura relativo.

Este módulo es **solo de lectura**: detecta pares de CDS que comparten contig
y coordenadas, describe el solapamiento y advierte para revisión manual. Nunca
elimina, recorta ni elige un gen "correcto"; devuelve filas de auditoría.

Convenciones
------------
- Coordenadas 1-based inclusivas (``start <= end``); ``end`` es el último
  nucleótido del CDS.
- Hebra ``+``/``-``. En la hebra ``-`` el inicio 5' del CDS es ``end`` y la
  transcripción avanza hacia coordenadas menores.
- Solapamiento: ``max(start_a, start_b) <= min(end_a, end_b)``. Dos CDS
  contiguos (``end_a + 1 == start_b``) **no** se consideran solapados.
- El marco relativo (``frame_shift``) se calcula únicamente si se aportan las
  secuencias. Sin genoma el valor es ``"unknown"``: no se inventa un marco a
  partir de coordenadas que no se han podido comprobar. En hebras opuestas el
  marco relativo tampoco es comparable y también queda ``"unknown"`` con
  ``frame_status = "opposite_strand"``.

``frame_shift`` es el desfase en pb entre los inicios 5' del primer y del
segundo CDS, medido en el sentido de transcripción de su hebra, módulo 3:
``0`` = mismo marco; ``1``/``2`` = marcos distintos.
"""
from __future__ import annotations

import csv
import io
from collections.abc import Mapping, Sequence
from typing import Any

FIELDS: tuple[str, ...] = (
    "contig",
    "gene_a",
    "gene_b",
    "start_a",
    "end_a",
    "strand_a",
    "start_b",
    "end_b",
    "strand_b",
    "strand_relation",
    "overlap_start",
    "overlap_end",
    "overlap_bp",
    "frame_reference",
    "frame_target",
    "frame_shift",
    "frame_status",
    "frame_reason",
    "severity",
    "warning",
)

FRAME_STATUSES: tuple[str, ...] = (
    "same_frame",
    "different_frame",
    "opposite_strand",
    "unknown",
)


def _as_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _normalise_strand(value: Any) -> str:
    text = str(value or "").strip()
    if text in ("+", "+1", "1", "plus", "forward"):
        return "+"
    if text in ("-", "-1", "minus", "reverse"):
        return "-"
    return ""


def _normalise(feature: Any, index: int) -> dict[str, Any] | None:
    """Devuelve una copia normalizada del CDS o ``None`` si no es auditable."""
    if not isinstance(feature, Mapping):
        return None
    if str(feature.get("type", "")).upper() != "CDS":
        return None
    contig = str(feature.get("contig", "") or "")
    if not contig:
        return None
    start, end = _as_int(feature.get("start")), _as_int(feature.get("end"))
    if start is None or end is None or start < 1 or end < 1:
        return None
    if end < start:
        start, end = end, start
    gene_id = str(
        feature.get("id")
        or feature.get("protein_id")
        or feature.get("gene_id")
        or feature.get("locus_tag")
        or f"cds_{index + 1}"
    )
    return {
        "contig": contig,
        "id": gene_id,
        "start": start,
        "end": end,
        "strand": _normalise_strand(feature.get("strand")),
        "_index": index,
    }


def _five_prime(cds: Mapping[str, Any]) -> int:
    """Coordenada 1-based del primer nucleótido traducido (inicio 5')."""
    return int(cds["start"]) if cds["strand"] == "+" else int(cds["end"])


def _ordered_pair(
    a: Mapping[str, Any], b: Mapping[str, Any], strand: str
) -> tuple[Mapping[str, Any], Mapping[str, Any], int]:
    """Ordena dos CDS de la misma hebra en sentido de transcripción.

    Devuelve ``(primero, segundo, desfase_5')``. En hebra ``-`` el CDS
    transcrito primero es el de mayor coordenada 5' (``end``).
    """
    five_a, five_b = _five_prime(a), _five_prime(b)
    if strand == "+":
        first, second = (a, b) if five_a <= five_b else (b, a)
        return first, second, abs(five_b - five_a)
    first, second = (a, b) if five_a >= five_b else (b, a)
    return first, second, abs(five_a - five_b)


def _frame_result(
    a: Mapping[str, Any],
    b: Mapping[str, Any],
    strand: str,
    sequences: Mapping[str, Any] | None,
) -> tuple[str, str, int | str, str, str]:
    """Calcula marco relativo y motivo.

    Devuelve ``(reference_id, target_id, frame_shift, frame_status,
    frame_reason)``.
    """
    first, second, offset = _ordered_pair(a, b, strand)
    reference, target = str(first["id"]), str(second["id"])
    if sequences is None:
        return (reference, target, "unknown", "unknown",
                "sin secuencias: no se calcula el marco relativo")
    seq = sequences.get(str(a["contig"]))
    if not isinstance(seq, str) or not seq:
        return (reference, target, "unknown", "unknown",
                f"contig {a['contig']} ausente en sequences")
    if max(int(a["end"]), int(b["end"])) > len(seq):
        return (reference, target, "unknown", "unknown",
                "coordenadas fuera de la secuencia aportada")
    shift = offset % 3
    status = "same_frame" if shift == 0 else "different_frame"
    reason = (
        f"desfase {offset} pb entre inicios 5' en hebra {strand}: "
        f"{'mismo marco' if shift == 0 else 'marcos distintos'} (mod 3 = {shift})"
    )
    return reference, target, shift, status, reason


def _warning(
    relation: str, frame_status: str, overlap_bp: int
) -> tuple[str, str]:
    """Devuelve ``(severity, texto)`` sin decidir qué gen es correcto."""
    if relation == "opposite":
        return ("info",
                f"CDS solapados {overlap_bp} pb en hebras opuestas (antisentido); "
                "puede ser regulación solapada. Solo auditoría, no implica error.")
    if relation == "unknown":
        return ("warning",
                f"CDS solapados {overlap_bp} pb con hebra desconocida; "
                "revisar coordenadas y hebra. Solo auditoría.")
    if frame_status == "same_frame":
        return ("warning",
                f"CDS solapados {overlap_bp} pb en el mismo marco; revisar si son "
                "dos genes, un inicio alternativo o coordenadas duplicadas. "
                "No se descarta ningún gen.")
    if frame_status == "different_frame":
        return ("warning",
                f"CDS solapados {overlap_bp} pb con marcos distintos; revisar "
                "posible error de coordenadas, fusión o CDS alternativo. "
                "No se descarta ningún gen.")
    return ("warning",
            f"CDS solapados {overlap_bp} pb; marco relativo no disponible. "
            "Revisar manualmente; no se descarta ningún gen.")


def analyze(
    features: Sequence[Mapping[str, Any]],
    sequences: Mapping[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Detecta todos los pares de CDS solapados y audita su marco relativo.

    Parámetros
    ----------
    features:
        Features del anotador. Solo se consideran ``type == "CDS"`` con
        ``contig`` no vacío y ``start``/``end`` enteros 1-based inclusivos.
        Se ignoran ``tRNA``, ``rRNA``, coordenadas inválidas y features que no
        sean mappings.
    sequences:
        Mapa opcional ``{contig: secuencia}`` en orientación ``+``. Si falta,
        o si el contig no está, el marco relativo queda ``"unknown"``.

    Returns
    -------
    list[dict]
        Una fila por par solapado, con las claves de :data:`FIELDS`. La lista
        se ordena por contig, inicio del solapamiento, fin e identificadores.
        Las features de entrada no se modifican.
    """
    if sequences is not None and not isinstance(sequences, Mapping):
        sequences = None

    by_contig: dict[str, list[dict[str, Any]]] = {}
    for index, feature in enumerate(features or []):
        cds = _normalise(feature, index)
        if cds is not None:
            by_contig.setdefault(cds["contig"], []).append(cds)

    rows: list[dict[str, Any]] = []
    for contig in sorted(by_contig):
        items = sorted(
            by_contig[contig],
            key=lambda c: (c["start"], c["end"], c["id"], c["_index"]),
        )
        for i, a in enumerate(items):
            for b in items[i + 1:]:
                # Ordenados por ``start``: si b empieza después del final de a,
                # ningún CDS posterior puede solapar con a.
                if int(b["start"]) > int(a["end"]):
                    break
                overlap_start = max(int(a["start"]), int(b["start"]))
                overlap_end = min(int(a["end"]), int(b["end"]))
                if overlap_start > overlap_end:
                    continue
                overlap_bp = overlap_end - overlap_start + 1

                strand_a = str(a["strand"])
                strand_b = str(b["strand"])
                if strand_a and strand_b:
                    relation = "same" if strand_a == strand_b else "opposite"
                else:
                    relation = "unknown"

                reference = target = ""
                if relation == "same":
                    reference, target, frame_shift, frame_status, reason = (
                        _frame_result(a, b, strand_a, sequences)
                    )
                elif relation == "opposite":
                    frame_shift = "unknown"
                    frame_status = "opposite_strand"
                    reason = ("hebras opuestas: no hay marco de lectura "
                              "relativo comparable")
                else:
                    frame_shift = "unknown"
                    frame_status = "unknown"
                    reason = "hebra desconocida: no se puede calcular el marco"

                severity, warning = _warning(relation, frame_status, overlap_bp)
                rows.append({
                    "contig": contig,
                    "gene_a": str(a["id"]),
                    "gene_b": str(b["id"]),
                    "start_a": int(a["start"]),
                    "end_a": int(a["end"]),
                    "strand_a": strand_a,
                    "start_b": int(b["start"]),
                    "end_b": int(b["end"]),
                    "strand_b": strand_b,
                    "strand_relation": relation,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_bp": overlap_bp,
                    "frame_reference": reference,
                    "frame_target": target,
                    "frame_shift": frame_shift,
                    "frame_status": frame_status,
                    "frame_reason": reason,
                    "severity": severity,
                    "warning": warning,
                })

    rows.sort(key=lambda r: (str(r["contig"]), int(r["overlap_start"]),
                             int(r["overlap_end"]), str(r["gene_a"]),
                             str(r["gene_b"])))
    return rows


def to_tsv(rows: Sequence[Mapping[str, Any]]) -> str:
    """Serializa filas de auditoría a TSV (cabecera siempre presente).

    Los valores ``None`` se escriben como cadena vacía y las claves extra se
    ignoran. Devuelve una cadena, no escribe archivos.
    """
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=FIELDS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows or []:
        if not isinstance(row, Mapping):
            continue
        writer.writerow({key: "" if row.get(key) is None else row.get(key)
                         for key in FIELDS})
    return buffer.getvalue()


def summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Resume las filas: pares, genes implicados, pb y estados de marco.

    No modifica ``rows``. ``audit_only`` es siempre ``True`` para dejar claro
    que ningún gen se ha eliminado ni elegido.
    """
    clean = [row for row in rows or [] if isinstance(row, Mapping)]
    genes: set[str] = set()
    overlaps: list[int] = []
    counts = {status: 0 for status in FRAME_STATUSES}
    same_strand = opposite = unknown_strand = 0
    warnings = 0
    for row in clean:
        for key in ("gene_a", "gene_b"):
            gene = str(row.get(key, "") or "")
            if gene:
                genes.add(gene)
        overlaps.append(int(row.get("overlap_bp", 0) or 0))
        relation = str(row.get("strand_relation", "unknown"))
        if relation == "same":
            same_strand += 1
        elif relation == "opposite":
            opposite += 1
        else:
            unknown_strand += 1
        status = str(row.get("frame_status", "unknown"))
        if status not in counts:
            status = "unknown"
        counts[status] += 1
        if str(row.get("severity", "")) == "warning":
            warnings += 1
    total = sum(overlaps)
    return {
        "n_pairs": len(clean),
        "n_same_strand": same_strand,
        "n_opposite_strand": opposite,
        "n_unknown_strand": unknown_strand,
        "n_genes": len(genes),
        "genes": sorted(genes),
        "overlap_bp_total": total,
        "overlap_bp_max": max(overlaps) if overlaps else 0,
        "frame_same": counts["same_frame"],
        "frame_different": counts["different_frame"],
        "frame_opposite_strand": counts["opposite_strand"],
        "frame_unknown": counts["unknown"],
        "n_warnings": warnings,
        "audit_only": True,
        "message": (
            f"{len(clean)} par(es) de CDS solapados; auditoría de solo lectura, "
            "ningún gen eliminado"
        ),
    }
