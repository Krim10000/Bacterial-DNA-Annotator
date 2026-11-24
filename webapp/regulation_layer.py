"""Operones y promotores putativos por reglas (CPU, sin bases externas).

Esta capa **no lee RegulonDB/DBTBS ni usa red**. Construye dos artefactos
declaradamente hipotéticos a partir de las features ya anotadas:

- **Operón putativo**: agrupa CDS del mismo contig y hebra cuando el hueco
  entre genes consecutivos es ``<= max_gap`` (por defecto 300 pb). Es una
  regla de distancia, no una prueba de cotranscripción.
- **Promotor putativo**: ventana de 50 pb inmediatamente upstream del primer
  gen del operón, en la orientación de la hebra (``+``/``-``). Si hay
  secuencia disponible, busca el mejor par de cajas de consenso sigma70
  ``-35 TTGACA`` / ``-10 TATAAT`` con espaciador 15-19 pb y puntúa las
  coincidencias (0-1). No es un TSS experimental.

Cada fila escrita lleva ``evidence="rule_based"``, ``confidence`` (solo
``low``/``medium``) y una ``limitation`` explícita. El estado
``hypothesis_not_experimental`` recuerda que ninguna de las dos tablas es
evidencia experimental.

La agrupación es determinista: se ordenan contigs/hebras y genes por
coordenadas, y los identificadores se numeran en ese orden.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Mapping

DEFAULT_MAX_GAP = 300
UPSTREAM_WINDOW = 50
MINUS10 = "TATAAT"
MINUS35 = "TTGACA"
SPACER_MIN = 15
SPACER_MAX = 19
MIN_BOX_MATCHES = 4
MIN_TOTAL_MATCHES = 8
TIGHT_GAP = 150

EVIDENCE = "rule_based"
STATUS = "hypothesis_not_experimental"
OPERON_LABEL = "operon_rule_based"
PROMOTER_LABEL = "promoter_rule_based"

OPERON_LIMITATION = (
    "Hipótesis no experimental: agrupa CDS del mismo contig y hebra por "
    "hueco <= max_gap. No demuestra cotranscripción, no usa TSS, terminador "
    "ni expresión, y puede unir o partir operones reales."
)
PROMOTER_LIMITATION = (
    "Hipótesis no experimental: ventana fija de 50 pb upstream del primer "
    "gen y coincidencia de consenso -35 TTGACA / -10 TATAAT (espaciador "
    "15-19). No es un TSS ni un promotor confirmado; no modela sigma "
    "alternativo, UP-element, -10 extendida, reguladores ni promotores "
    "internos. Un motivo puede ser un falso positivo."
)

OPERON_FIELDS = (
    "genome", "operon_id", "contig", "strand", "start", "end", "length",
    "n_genes", "genes", "gene_coords", "internal_gaps", "max_internal_gap",
    "promoter_id", "label", "evidence", "confidence", "confidence_reason",
    "status", "limitation",
)
PROMOTER_FIELDS = (
    "genome", "promoter_id", "operon_id", "contig", "strand", "gene_id",
    "start", "end", "length", "window", "sequence", "box35", "mm35",
    "box10", "mm10", "spacer", "score", "motif", "box_found",
    "sequence_status", "label", "evidence", "confidence",
    "confidence_reason", "status", "limitation",
)

_COMPLEMENT = str.maketrans(
    "ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverso complementario (A<->T, C<->G, N->N)."""
    return str(seq).translate(_COMPLEMENT)[::-1]


def _as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _sequence_for(sequences: Mapping[str, str] | None,
                  contig: str) -> str | None:
    if not sequences:
        return None
    seq = sequences.get(contig)
    return str(seq).upper() if seq is not None else None


def upstream_region(seq: str | None, start: int, end: int, strand: str,
                    window: int = UPSTREAM_WINDOW) -> dict[str, Any]:
    """Ventana upstream del gen en orientación de transcripción.

    ``start``/``end`` son 1-based inclusivos (convención de pyrodigal/GFF3).
    Para ``+`` toma ``[start-window, start-1]``; para ``-`` toma
    ``[end+1, end+window]`` y lo reverse-complementa. Si la secuencia del
    contig no está disponible, devuelve coordenadas nominales y
    ``sequence=""`` con ``sequence_status="no_sequence_available"``.
    """
    if strand == "+":
        a = max(1, start - window)
        b = start - 1
        if seq is not None:
            b = min(b, len(seq))
        if b < a:
            return {"sequence": "", "start": "", "end": "", "length": 0,
                    "sequence_status": "no_upstream"}
        region = "" if seq is None else seq[a - 1:b]
        status = "no_sequence_available" if seq is None else "available"
    else:
        a = end + 1
        b = end + window
        if seq is not None:
            b = min(b, len(seq))
        if a > b:
            return {"sequence": "", "start": "", "end": "", "length": 0,
                    "sequence_status": "no_upstream"}
        region = "" if seq is None else revcomp(seq[a - 1:b])
        status = "no_sequence_available" if seq is None else "available"
    return {"sequence": region, "start": a, "end": b, "length": len(region),
            "sequence_status": status}


def _matches(seq: str, motif: str) -> int:
    return sum(1 for a, b in zip(seq, motif) if a == b)


def promoter_motif(upstream: str) -> dict[str, Any]:
    """Mejor par -35/-10 de consenso en la ventana upstream.

    Devuelve ``box35``/``box10`` (secuencias encontradas), ``mm35``/``mm10``
    (mismatches), ``spacer``, ``score`` (coincidencias/12), ``motif`` y
    ``box_found``. Sin secuencia devuelve ``box_found="no"`` y
    ``sequence_status`` se resuelve en la fila.
    """
    up = (upstream or "").upper()
    n = len(up)
    empty: dict[str, Any] = {
        "box35": "", "mm35": "", "box10": "", "mm10": "", "spacer": "",
        "score": 0.0, "motif": "", "box_found": "no",
    }
    if n < SPACER_MIN + 12:
        return empty

    # La caja -35 se busca a -42/-26 y la -10 termina a -14/-4 respecto del
    # inicio del gen (posición n de la ventana). Así el motivo queda en una
    # posición promotora plausible y no en cualquier parte de la ventana.
    lo35 = max(0, n - 42)
    hi35 = max(lo35, n - 26)
    lo10_end = max(6, n - 14)
    hi10_end = min(n, n - 4)
    best: tuple[tuple[int, int, int], dict[str, Any]] | None = None
    for i in range(lo35, min(hi35, n - 6) + 1):
        box35 = up[i:i + 6]
        if len(box35) < 6:
            continue
        m35 = _matches(box35, MINUS35)
        if m35 < MIN_BOX_MATCHES:
            continue
        for spacer in range(SPACER_MIN, SPACER_MAX + 1):
            j = i + 6 + spacer
            e = j + 6
            if e < lo10_end or e > hi10_end or e > n:
                continue
            box10 = up[j:e]
            m10 = _matches(box10, MINUS10)
            if m10 < MIN_BOX_MATCHES:
                continue
            total = m35 + m10
            if total < MIN_TOTAL_MATCHES:
                continue
            # Determinista: más coincidencias, luego espaciador más cercano a
            # 17 pb, luego la caja -35 más upstream.
            key = (total, -abs(spacer - 17), -i)
            cand = {
                "box35": box35, "mm35": 6 - m35,
                "box10": box10, "mm10": 6 - m10,
                "spacer": spacer, "score": round(total / 12.0, 2),
                "motif": f"-35:{box35};-10:{box10};esp={spacer}",
                "box_found": "yes",
            }
            if best is None or key > best[0]:
                best = (key, cand)
    return best[1] if best else empty


def _confidence(n_genes: int, max_internal_gap: int | None) -> tuple[str, str]:
    if n_genes >= 2 and max_internal_gap is not None \
            and max_internal_gap <= TIGHT_GAP:
        return "medium", f"n_genes>=2 y hueco_interno_max<={TIGHT_GAP}"
    if n_genes == 1:
        return "low", "gen_unico_sin_vecino_en_la_hebra"
    return "low", f"hueco_interno_max>{TIGHT_GAP}"


def _write_tsv(path: Path, fields: tuple[str, ...],
               rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(features: list[dict[str, Any]], outdir: Path | str, genome: str,
        max_gap: int = DEFAULT_MAX_GAP,
        sequences: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Construye operones/promotores putativos y escribe dos TSV.

    Parámetros
    ----------
    features:
        Lista de features del anotador (solo se usan ``type == "CDS"`` con
        ``contig``, ``start``/``end`` 1-based inclusivos, ``strand`` y ``id``).
    outdir:
        Directorio de salida; se crea si no existe.
    genome:
        Nombre del genoma, usado en los nombres de archivo e identificadores.
    max_gap:
        Hueco máximo (pb) para unir CDS consecutivos en el mismo contig/hebra.
    sequences:
        Mapa opcional ``{contig: secuencia}`` para calcular el motivo -10/-35.
        Sin él, los promotores se escriben con coordenadas nominales y sin
        score (``sequence_status="no_sequence_available"``).

    Escribe ``<genome>.operons_auto.tsv`` y ``<genome>.promoters_auto.tsv``.
    Devuelve ``{"ok": True, "n_operons": ..., "n_promoters": ...,
    "rows": [operones], "promoter_rows": [promotores], ...}``.
    """
    genome = str(genome or "genome")
    outdir = Path(outdir)
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return {"ok": False, "n_operons": 0, "n_promoters": 0, "rows": [],
                "promoter_rows": [], "reason": f"no se pudo crear {outdir}: {exc}"}

    try:
        gap_limit = max(0, int(max_gap))
    except (TypeError, ValueError):
        gap_limit = DEFAULT_MAX_GAP

    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for idx, feature in enumerate(features or []):
        if not isinstance(feature, Mapping):
            continue
        if str(feature.get("type", "")).upper() != "CDS":
            continue
        strand = str(feature.get("strand", ""))
        if strand not in ("+", "-"):
            continue
        contig = str(feature.get("contig", ""))
        start = _as_int(feature.get("start"))
        end = _as_int(feature.get("end"))
        if start is None or end is None or start < 1 or end < 1:
            continue
        if end < start:
            start, end = end, start
        gene_id = str(feature.get("id") or feature.get("protein_id")
                      or f"cds_{idx + 1}")
        groups.setdefault((contig, strand), []).append({
            "id": gene_id, "start": start, "end": end,
        })

    operon_rows: list[dict[str, Any]] = []
    promoter_rows: list[dict[str, Any]] = []
    index = 0
    for contig, strand in sorted(groups):
        genes = sorted(groups[(contig, strand)],
                       key=lambda g: (g["start"], g["end"], g["id"]))
        clusters: list[list[dict[str, Any]]] = []
        current: list[dict[str, Any]] = []
        for gene in genes:
            if not current:
                current = [gene]
                continue
            gap = gene["start"] - current[-1]["end"] - 1
            if gap <= gap_limit:
                current.append(gene)
            else:
                clusters.append(current)
                current = [gene]
        if current:
            clusters.append(current)

        seq = _sequence_for(sequences, contig)
        for cluster in clusters:
            index += 1
            operon_id = f"{genome}:op:{index:05d}"
            promoter_id = f"{genome}:prom:{index:05d}"
            op_start = min(g["start"] for g in cluster)
            op_end = max(g["end"] for g in cluster)
            gaps = [b["start"] - a["end"] - 1
                    for a, b in zip(cluster, cluster[1:])]
            max_internal_gap = max(gaps) if gaps else None
            confidence, reason = _confidence(len(cluster), max_internal_gap)
            # Los genes se escriben en orden de transcripción (5'->3'): en la
            # hebra - eso es el orden inverso de las coordenadas crecientes.
            tx_cluster = (cluster if strand == "+"
                          else list(reversed(cluster)))
            tx_gaps = gaps if strand == "+" else list(reversed(gaps))
            first = tx_cluster[0]

            promoter = upstream_region(seq, first["start"], first["end"],
                                       strand)
            motif = promoter_motif(promoter["sequence"])
            if promoter["sequence_status"] != "available":
                motif = {"box35": "", "mm35": "", "box10": "", "mm10": "",
                         "spacer": "", "score": 0.0, "motif": "",
                         "box_found": "no"}
            if motif["box_found"] == "yes" and float(motif["score"]) >= 0.75:
                p_conf, p_reason = "medium", "motivo_-35/-10_con_score>=0.75"
            elif motif["box_found"] == "yes":
                p_conf, p_reason = "low", "motivo_debil_o_incompleto"
            elif promoter["sequence_status"] != "available":
                p_conf, p_reason = "low", "sin_secuencia_para_evaluar_motivo"
            else:
                p_conf, p_reason = "low", "sin_motivo_-35/-10_plausible"

            operon_rows.append({
                "genome": genome,
                "operon_id": operon_id,
                "contig": contig,
                "strand": strand,
                "start": op_start,
                "end": op_end,
                "length": op_end - op_start + 1,
                "n_genes": len(cluster),
                "genes": ";".join(g["id"] for g in tx_cluster),
                "gene_coords": ";".join(
                    f"{g['id']}:{g['start']}-{g['end']}" for g in tx_cluster),
                "internal_gaps": ";".join(str(g) for g in tx_gaps),
                "max_internal_gap": "" if max_internal_gap is None
                else max_internal_gap,
                "promoter_id": promoter_id,
                "label": OPERON_LABEL,
                "evidence": EVIDENCE,
                "confidence": confidence,
                "confidence_reason": reason,
                "status": STATUS,
                "limitation": OPERON_LIMITATION,
            })
            promoter_rows.append({
                "genome": genome,
                "promoter_id": promoter_id,
                "operon_id": operon_id,
                "contig": contig,
                "strand": strand,
                "gene_id": first["id"],
                "start": promoter["start"],
                "end": promoter["end"],
                "length": promoter["length"],
                "window": f"-{UPSTREAM_WINDOW}/-1",
                "sequence": promoter["sequence"],
                "box35": motif["box35"],
                "mm35": motif["mm35"],
                "box10": motif["box10"],
                "mm10": motif["mm10"],
                "spacer": motif["spacer"],
                "score": motif["score"],
                "motif": motif["motif"],
                "box_found": motif["box_found"],
                "sequence_status": promoter["sequence_status"],
                "label": PROMOTER_LABEL,
                "evidence": EVIDENCE,
                "confidence": p_conf,
                "confidence_reason": p_reason,
                "status": STATUS,
                "limitation": PROMOTER_LIMITATION,
            })

    operons_tsv = outdir / f"{genome}.operons_auto.tsv"
    promoters_tsv = outdir / f"{genome}.promoters_auto.tsv"
    try:
        _write_tsv(operons_tsv, OPERON_FIELDS, operon_rows)
        _write_tsv(promoters_tsv, PROMOTER_FIELDS, promoter_rows)
    except OSError as exc:
        return {"ok": False, "n_operons": len(operon_rows),
                "n_promoters": len(promoter_rows), "rows": operon_rows,
                "promoter_rows": promoter_rows,
                "reason": f"no se pudieron escribir los TSV: {exc}"}

    return {
        "ok": True,
        "n_operons": len(operon_rows),
        "n_promoters": len(promoter_rows),
        "rows": operon_rows,
        "promoter_rows": promoter_rows,
        "operons_tsv": str(operons_tsv),
        "promoters_tsv": str(promoters_tsv),
        "max_gap": gap_limit,
        "reason": "",
    }
