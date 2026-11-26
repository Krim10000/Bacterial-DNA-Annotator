"""Regiones oriC: intervalo DoriC y punto ORCA, sin estimar intervalos.

ORCA devuelve un **punto**; DoriC publica un **intervalo** para genomas con
accession. Este módulo separa ambas evidencias:

- si hay una región DoriC válida, ``start``/``end`` son sus coordenadas;
- si solo hay punto, ``start == end == pos`` y ``span_source`` es
  ``"ORCA (punto)"``; nunca se inventa un ancho.

DoriC usa convenciones de coordenadas heterogéneas (0/1-based según la fuente,
``end == 0`` para regiones que llegan al origen y regiones que cruzan el
origen). Por eso el intervalo se normaliza a **segmentos** sin reescribir los
endpoints publicados: ``length_reported`` y ``length_effective`` pueden diferir
en ±1 pb y esa diferencia se reporta, no se corrige a ciegas.

Funciones puras (sin Streamlit, pandas ni acceso a disco) para poder probarlas
con secuencias y regiones sintéticas.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

POINT_SPAN_SOURCE = "ORCA (punto)"
DORIC_SPAN_SOURCE = "DoriC"
INDEX_RELATIVE = Path("oric_regions") / "oriC_regions.json"

#: Categorías de recurrencia de DoriC. No son confianza experimental: el
#: manifiesto de la release lo dice explícitamente. Solo desempatan.
_RECURRENCE_RANK = {
    "Recurrent": 4,
    "Multi-accession": 3,
    "Cross-taxon": 2,
    "Singleton": 1,
    "Not analyzed": 0,
}


@dataclass(frozen=True)
class RegionSpan:
    """Intervalo oriC normalizado a segmentos 1-based inclusivos."""

    start: int
    end: int
    length: int
    segments: tuple[tuple[int, int], ...]
    genome_length: int | None = None
    wrap: bool = False
    end_effective: int | None = None
    length_reported: int | None = None
    length_delta: int | None = None

    def contains(self, pos: int) -> bool:
        """True si ``pos`` (1-based) cae en algún segmento."""
        return any(a <= pos <= b for a, b in self.segments)

    def distance_to(self, pos: int) -> int:
        """Distancia circular en pb al borde más cercano; 0 si está dentro."""
        if self.contains(pos):
            return 0
        if not self.segments:
            return 0
        length = self.genome_length
        best = None
        for a, b in self.segments:
            for edge in (a, b):
                if length and length > 0:
                    d = min((pos - edge) % length, (edge - pos) % length)
                else:
                    d = abs(pos - edge)
                best = d if best is None else min(best, d)
        return int(best or 0)

    def extract(self, sequence: str, *, strand: str = "+") -> str:
        """Extrae la secuencia de los segmentos; ``""`` si algo no cuadra."""
        if not sequence:
            return ""
        try:
            pieces = [sequence[a - 1:b] for a, b in self.segments]
        except (IndexError, TypeError):
            return ""
        if len(pieces) != len(self.segments) or any(
                len(p) != b - a + 1 for p, (a, b) in zip(pieces, self.segments)):
            return ""
        out = "".join(pieces).upper()
        if str(strand) == "-":
            out = out.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        return out

    def as_dict(self) -> dict[str, Any]:
        return {
            "start": self.start,
            "end": self.end,
            "end_effective": self.end_effective,
            "length": self.length,
            "length_reported": self.length_reported,
            "length_delta": self.length_delta,
            "wrap": self.wrap,
            "segments": [list(s) for s in self.segments],
        }


def normalize_region(start: int | str | None, end: int | str | None, *,
                     length: int | str | None = None,
                     genome_length: int | None = None) -> RegionSpan | None:
    """Normaliza una región publicada sin inventar endpoints.

    ``start`` es 1-based; ``end`` puede ser 0 (centinela de "hasta el origen")
    y puede ser menor que ``start`` (región que cruza el origen). Si
    ``length`` viene de DoriC y ``end == 0``, se usa la longitud declarada para
    reconstruir el final, porque ``end = 0`` no es una coordenada utilizable.
    Devuelve ``None`` si la región no es válida o si falta la longitud del
    genoma para resolver un cruce de origen.
    """
    try:
        s = int(start)  # type: ignore[arg-type]
        e = int(end)  # type: ignore[arg-type]
        rep = int(length) if length not in (None, "") else None  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if s < 1 or e < 0 or (rep is not None and rep < 1):
        return None
    if genome_length is not None:
        genome_length = int(genome_length)
        if genome_length < 1 or s > genome_length or e > genome_length:
            return None
        if rep is not None and rep > genome_length:
            return None

    if genome_length is None:
        # Sin longitud del genoma no se puede resolver un cruce de origen.
        if e == 0 or e < s:
            return None
        segments = ((s, e),)
        wrap = False
        end_effective = e
        effective = e - s + 1
    elif e == 0:
        # Centinela de DoriC: el intervalo llega al origen. La longitud
        # declarada (y ``len(oric_seq)``) es la fuente más fina.
        if rep is not None:
            end_effective = s + rep - 1
            if end_effective <= genome_length:
                segments = ((s, end_effective),)
                wrap = False
            else:
                remainder = end_effective - genome_length
                segments = ((s, genome_length), (1, remainder))
                wrap = True
                end_effective = remainder
        else:
            segments = ((s, genome_length),)
            wrap = True
            end_effective = genome_length
        effective = sum(b - a + 1 for a, b in segments)
    elif e < s:
        # Cruce de origen: [start, L] ∪ [1, end].
        segments = ((s, genome_length), (1, e))
        wrap = True
        end_effective = e
        effective = (genome_length - s + 1) + e
    else:
        segments = ((s, e),)
        wrap = False
        end_effective = e
        effective = e - s + 1

    delta = None if rep is None else effective - rep
    return RegionSpan(
        start=s, end=e, length=effective, segments=segments,
        genome_length=genome_length, wrap=wrap,
        end_effective=end_effective, length_reported=rep,
        length_delta=delta,
    )


def point_in_region(pos: int | str, region: RegionSpan) -> bool:
    """``True`` si el punto cae en la región (incluye cruces de origen)."""
    try:
        return region.contains(int(pos))
    except (TypeError, ValueError):
        return False


def point_to_region_bp(pos: int | str, region: RegionSpan) -> int:
    """Distancia circular del punto a la región; 0 si está dentro."""
    try:
        return region.distance_to(int(pos))
    except (TypeError, ValueError):
        return 0


def region_sequence(seqs: Mapping[str, str], region: RegionSpan, *,
                    contig: str = "", strand: str = "+") -> str:
    """Extrae la secuencia de una ``RegionSpan``.

    Si ``contig`` no existe y hay un único contig, se usa ese. No descarga ni
    rellena huecos: si el intervalo no cabe, devuelve ``""``.
    """
    if not seqs:
        return ""
    seq = seqs.get(contig) if contig else None
    if seq is None and len(seqs) == 1:
        seq = next(iter(seqs.values()))
    if not seq:
        return ""
    if region.genome_length is None or region.genome_length != len(seq):
        # El intervalo se normalizó con otra longitud: no mezclar coordinates.
        region = normalize_region(
            region.start, region.end, length=region.length_reported,
            genome_length=len(seq))
        if region is None:
            return ""
    return region.extract(seq, strand=strand)


def extract_region_sequence(seqs: Mapping[str, str], contig: str,
                            start: int | str | None, end: int | str | None,
                            *, length: int | str | None = None,
                            strand: str = "+") -> str:
    """Atajo puro para ``load_oric``/tests: normaliza y extrae en un paso."""
    if not seqs:
        return ""
    seq = seqs.get(contig) if contig else None
    if seq is None and len(seqs) == 1:
        seq = next(iter(seqs.values()))
    if not seq:
        return ""
    region = normalize_region(start, end, length=length,
                              genome_length=len(seq))
    if region is None:
        return ""
    return region.extract(seq, strand=strand)


def _candidate(record: Mapping[str, Any], genome_length: int) -> dict | None:
    """Normaliza un registro DoriC (dict con start/end/length) y lo anota."""
    if not isinstance(record, Mapping):
        return None
    span = normalize_region(
        record.get("start", record.get("oric_start")),
        record.get("end", record.get("oric_end")),
        length=record.get("length", record.get("oric_length")),
        genome_length=genome_length,
    )
    if span is None:
        return None
    support = str(record.get("comparative_support") or "")
    return {
        "start": span.start,
        "end": span.end,
        "length": span.length,
        "length_reported": span.length_reported,
        "length_delta": span.length_delta,
        "wrap": span.wrap,
        "segments": [list(s) for s in span.segments],
        "support": support,
        "recurrence_rank": _RECURRENCE_RANK.get(support, 0),
        "ori_id": str(record.get("ori_id") or ""),
        "organism": str(record.get("organism") or ""),
        "at_content": record.get("at_content", record.get("oric_at_content")),
        "span": span,
    }


def select_doric_region(records: Iterable[Mapping[str, Any]],
                        genome_length: int,
                        point: int | None = None) -> dict[str, Any] | None:
    """Elige una región DoriC entre varias publicadas para el mismo genoma.

    Regla explícita y auditable:

    1. si el punto (ancla ORCA/dnaA/dnaN) cae dentro de una candidata, se elige
       esa (concordancia);
    2. si no, la candidata DoriC más cercana al punto;
    3. sin punto, la de mayor recurrencia y luego la más larga.

    No se recorta ni se extiende ninguna candidata: todas son intervalos
    publicados. Las demás quedan en ``alternatives`` con su distancia.
    """
    candidates = [c for c in (_candidate(r, genome_length) for r in records)
                  if c is not None]
    if not candidates:
        return None
    # Deduplicar intervalos idénticos conservando la primera aparición.
    seen: set[tuple[int, int, int | None]] = set()
    unique = []
    for c in candidates:
        key = (c["start"], c["end"], c["length_reported"])
        if key not in seen:
            seen.add(key)
            unique.append(c)
    candidates = unique
    for c in candidates:
        c["distance_bp"] = c["span"].distance_to(int(point)) if point is not None else None
        c["point_inside"] = c["span"].contains(int(point)) if point is not None else None

    if point is not None:
        inside = [c for c in candidates if c["point_inside"]]
        if inside:
            selected = max(inside, key=lambda c: (c["recurrence_rank"], c["length"]))
            how = "concordante_con_el_punto"
        else:
            selected = min(
                candidates,
                key=lambda c: (c["distance_bp"], -c["recurrence_rank"], -c["length"]))
            how = "mas_cercana_al_punto"
    else:
        selected = max(candidates, key=lambda c: (c["recurrence_rank"], c["length"]))
        how = "mayor_recurrencia"

    alternatives = [c for c in candidates if c is not selected]
    alternatives.sort(key=lambda c: (
        c["distance_bp"] if c["distance_bp"] is not None else 0,
        -c["recurrence_rank"], -c["length"]))
    return {
        "selected": _public(selected),
        "alternatives": [_public(c) for c in alternatives],
        "selection": how,
        "n_candidates": len(candidates),
    }


def _public(dict_or_span: Mapping[str, Any]) -> dict[str, Any]:
    """Quita objetos ``RegionSpan`` para serializar a JSON."""
    out = {}
    for k, v in dict_or_span.items():
        if k == "span":
            continue
        out[k] = v
    return out


def resolve_index_entry(entry: Mapping[str, Any]) -> dict[str, Any] | None:
    """Convierte una entrada del índice en un dict compatible con ``load_oric``.

    Mantiene las claves históricas (``pos``, ``start``, ``end``, ``length``,
    ``confidence``, ``source``, ``span_source``, ``replicon``, ``doric``) y
    añade ``region`` (segmentos) y ``point`` (procedencia del ancla).
    """
    if not isinstance(entry, Mapping):
        return None
    try:
        pos = int(entry["pos"])
    except (KeyError, TypeError, ValueError):
        return None
    region = entry.get("region")
    span = None
    if isinstance(region, Mapping):
        span = normalize_region(
            region.get("start"), region.get("end"),
            length=region.get("length_reported"),
            genome_length=entry.get("genome_length"))
    start = span.start if span else pos
    end = span.end if span else pos
    out = {
        "genome": entry.get("genome", ""),
        "pos": pos,
        "start": start,
        "end": end,
        "length": entry.get("genome_length"),
        "confidence": entry.get("confidence"),
        "source": entry.get("source", ""),
        "span_source": entry.get("span_source")
        or (DORIC_SPAN_SOURCE if span else POINT_SPAN_SOURCE),
        "replicon": entry.get("replicon", ""),
        "doric": entry.get("doric"),
        "region": span.as_dict() if span else None,
        "point": entry.get("point") or {},
        "alternatives": entry.get("alternatives") or [],
        "discrepancy": entry.get("discrepancy") or {},
    }
    return out


@lru_cache(maxsize=4)
def _load_index_cached(path: str) -> dict[str, Any]:
    p = Path(path)
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def load_region_index(data_root: Path | str) -> dict[str, Any]:
    """Carga el índice compacto de regiones; ``{}`` si no existe.

    ``data_root`` es la carpeta de datos (la misma que recibe ``load_oric``,
    normalmente ``data/``). También se acepta la raíz del proyecto.
    """
    root = Path(data_root)
    # Fallback empaquetado para la release, que no versiona data/.
    bundled = Path(__file__).with_name("oriC_regions.json")
    for candidate in (root / INDEX_RELATIVE, root / "data" / INDEX_RELATIVE,
                      bundled):
        if candidate.exists():
            return _load_index_cached(str(candidate))
    return {}


def load_region(data_root: Path | str, genome: str) -> dict[str, Any] | None:
    """Entrada de región desde el índice; si no existe, cae a las predicciones.

    El respaldo mantiene el comportamiento histórico y añade el esquema de
    ``orca_predictions`` (``orca_position``), que ``load_oric`` no leía.
    """
    index = load_region_index(data_root)
    genomes = index.get("genomes") if isinstance(index, dict) else None
    entry = genomes.get(genome) if isinstance(genomes, dict) else None
    if isinstance(entry, dict):
        resolved = resolve_index_entry(entry)
        if resolved:
            return resolved
    return fallback_from_prediction(data_root, genome)


def fallback_from_prediction(data_root: Path | str, genome: str) -> dict[str, Any] | None:
    """Comportamiento histórico de ``load_oric`` como respaldo sin índice.

    Lee ``orca_predictions_fix`` y luego ``orca_predictions``. Es deliberadamente
    conservador: si el JSON no trae ``doric_truth.orics``, devuelve punto.
    """
    root = Path(data_root)
    if not (root / "orca_predictions_fix").exists() and \
            (root / "data" / "orca_predictions_fix").exists():
        root = root / "data"
    for folder in ("orca_predictions_fix", "orca_predictions"):
        path = root / folder / f"{genome}.json"
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        chosen = data.get("chosen") or {}
        pos = chosen.get("pos")
        if pos is None:
            pos = data.get("orca_position")
        if pos is None:
            continue
        doric = data.get("doric_truth") or data.get("doric") or {}
        orics = doric.get("orics") if isinstance(doric, dict) else None
        first = orics[0] if isinstance(orics, list) and orics else None
        span = None
        if isinstance(first, dict):
            span = normalize_region(
                first.get("start"), first.get("end"),
                length=first.get("length"),
                genome_length=data.get("L"))
        start, end = (span.start, span.end) if span else (int(pos), int(pos))
        confidence = chosen.get("conf", data.get("orca_confidence"))
        source = data.get("source") or (
            "ORCA" if data.get("orca_position") is not None else "")
        return {
            "genome": genome, "pos": int(pos), "start": start, "end": end,
            "length": data.get("L"), "confidence": confidence,
            "source": source,
            "span_source": DORIC_SPAN_SOURCE if span else POINT_SPAN_SOURCE,
            "replicon": data.get("replicon") or "",
            "doric": data.get("doric"),
            "region": span.as_dict() if span else None,
            "point": {"pos": int(pos), "kind": chosen.get("cfg", "ORCA"),
                      "confidence": confidence},
            "alternatives": [],
            "discrepancy": {},
        }
    return None
