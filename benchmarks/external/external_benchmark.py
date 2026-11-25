#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark externo CPU y reproducible contra Bakta, Prokka y DFAST.

Diseño
------
- Lee un genoma FASTA y una anotación de referencia en GFF3 o GenBank (GBFF).
- Detecta ``bakta``, ``prokka`` y ``dfast`` con :func:`shutil.which`.
- Solo ejecuta herramientas cuando se pide explícitamente ``run=True``.
- Si una herramienta no está instalada, su resultado queda con estado
  ``not_available`` y **sin métricas inventadas**.
- Compara CDS, tRNA y rRNA por tipo, contig, hebra y solapamiento recíproco
  1:1; reporta TP/FP/FN, precisión, recall, F1, coincidencia exacta,
  diferencias de bordes y diferencias concretas.
- Permite comparar GFF/GBFF precalculados con ``query_gffs`` aunque el binario
  no exista; en ese caso la fuente se marca ``precomputed``.

El módulo no descarga bases de datos, imágenes ni datos de red. La receta CPU
está en ``benchmarks/external/Dockerfile`` y el uso en
``docs/EXTERNAL_BENCHMARK_RUN.md``.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from urllib.parse import unquote

SCHEMA = "external-benchmark-run/v1"

CATEGORIES: tuple[str, ...] = ("CDS", "tRNA", "rRNA")
TOOLS: tuple[str, ...] = ("bakta", "prokka", "dfast")

TOOL_SPECS: dict[str, dict[str, Any]] = {
    "bakta": {
        "label": "Bakta",
        "binary": "bakta",
        "version_args": ["--version"],
        "requires_database": True,
    },
    "prokka": {
        "label": "Prokka",
        "binary": "prokka",
        "version_args": ["--version"],
        "requires_database": False,
    },
    "dfast": {
        "label": "DFAST",
        "binary": "dfast",
        "version_args": ["--version"],
        "requires_database": False,
    },
}

TYPE_ALIASES = {
    "cds": "CDS",
    "trna": "tRNA",
    "rrna": "rRNA",
}

STATUS_LABELS: dict[str, str] = {
    "ok": "ok",
    "not_available": "no disponible (not_available)",
    "not_run": "no ejecutado (not_run)",
    "failed": "falló (failed)",
    "timeout": "timeout",
    "not_comparable": "no comparable (not_comparable)",
    "invalid_input": "entrada inválida",
}

HYPOTHETICAL_RE = re.compile(
    r"\b(?:hypothetical|uncharacterized|unknown\s+protein|"
    r"protein\s+of\s+unknown\s+function|duf\d+|upf\d+)\b",
    re.IGNORECASE,
)

_ACCESSION_RE = re.compile(r"^[A-Za-z]{1,8}_?\d+(?:\.\d+)?$")
_COMPOUND_WARNING = (
    "Las ubicaciones compuestas de GenBank se representan por su intervalo "
    "envolvente (start..end); no se separan intrones ni segmentos."
)

LIMITATIONS = [
    "La referencia es el gold de esta comparación; sus criterios y errores se heredan.",
    "El emparejamiento es 1:1 voraz por solapamiento recíproco y no un óptimo bipartito.",
    "La coincidencia literal de 'product' no mide función; los nombres dependen de cada base.",
    "No se comparan ncRNA, tmRNA, pseudogenes (si se excluyen), dominios ni rutas.",
    "Una herramienta ausente no es un cero: queda como not_available sin métricas.",
    "No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.",
    "Las ubicaciones compuestas de GenBank se comparan por su intervalo envolvente.",
]


# ---------------------------------------------------------------------------
# Errores y utilidades generales
# ---------------------------------------------------------------------------
class BenchmarkError(Exception):
    """Error controlado del benchmark."""


class InputError(BenchmarkError):
    """Entrada ausente, ilegible o sin features evaluables."""


class DependencyError(BenchmarkError):
    """Falta una dependencia opcional, por ejemplo Biopython para GBFF."""


def utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def display_path(path: str | Path | None) -> str | None:
    """Ruta relativa al cwd cuando sea posible; absoluta si no."""
    if path is None:
        return None
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return str(resolved)


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: str | Path):
    """Abre texto plano o gzip detectando magic bytes, no solo la extensión."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"no existe el archivo: {path}")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _round(value: float | None, ndigits: int = 4) -> float | None:
    return None if value is None else round(float(value), ndigits)


def _ratio(numerator: float, denominator: float) -> float | None:
    if not denominator:
        return None
    return _round(numerator / denominator)


def _mean(values: Iterable[float]) -> float | None:
    values = list(values)
    if not values:
        return None
    return _round(sum(values) / len(values))


def _pct(numerator: int, denominator: int) -> float | None:
    if not denominator:
        return None
    return round(100.0 * numerator / denominator, 2)


def _tail(value: Any, limit: int = 4000) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    return str(value)[-limit:]


def _cap(items: Sequence[Any], limit: int) -> tuple[list[Any], bool]:
    if limit < 0:
        return list(items), False
    return list(items[:limit]), len(items) > limit


def _parse_attributes(text: str) -> dict[str, str]:
    """Parsea GFF3 (``clave=valor``) y tolera GTF (``clave "valor"``)."""
    attrs: dict[str, str] = {}
    for part in (text or "").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            match = re.match(r'^(\S+)\s+"(.*)"$', part)
            if not match:
                continue
            key, value = match.group(1), match.group(2)
        key = unquote(key.strip())
        value = unquote(value.strip().strip('"'))
        if key:
            attrs[key] = value
    return attrs


def _first_qualifier(qualifiers: Mapping[str, Any], key: str) -> str:
    value = qualifiers.get(key)
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return " ".join(str(item) for item in value if item is not None).strip()
    return str(value).strip()


def _is_pseudo(attrs: Mapping[str, str]) -> bool:
    for key in ("pseudo", "pseudogene", "is_pseudogene"):
        value = str(attrs.get(key) or "").strip().lower()
        if value in {"", "true", "1", "yes"} and key in attrs:
            return True
    return False


def normalize_type(raw: str) -> str:
    text = (raw or "").strip()
    return TYPE_ALIASES.get(text.lower(), text)


def normalize_text(value: str) -> str:
    text = (value or "").lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def product_tokens(value: str) -> set[str]:
    return set(normalize_text(value).split())


def product_token_jaccard(a: str, b: str) -> float | None:
    left, right = product_tokens(a), product_tokens(b)
    union = left | right
    if not union:
        return None
    return len(left & right) / len(union)


def is_hypothetical(product: str) -> bool:
    """Clasificador conservador; ``putative`` a secas no cuenta como hipotética."""
    text = (product or "").strip()
    return bool(text and HYPOTHETICAL_RE.search(text))


# ---------------------------------------------------------------------------
# Modelo de datos
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Feature:
    seqid: str
    type: str
    start: int
    end: int
    strand: str
    attrs: Mapping[str, str] = field(default_factory=dict)
    source: str = "."

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def attr(self, key: str) -> str:
        return str(self.attrs.get(key) or "").strip()

    def with_seqid(self, seqid: str) -> "Feature":
        return Feature(
            seqid=seqid,
            type=self.type,
            start=self.start,
            end=self.end,
            strand=self.strand,
            attrs=self.attrs,
            source=self.source,
        )

    def summary(self) -> dict[str, Any]:
        return {
            "seqid": self.seqid,
            "type": self.type,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "product": self.attr("product"),
            "gene": self.attr("gene"),
            "locus_tag": self.attr("locus_tag"),
        }


@dataclass
class AnnotationSet:
    path: str
    format: str
    features: list[Feature]
    excluded_pseudo: int = 0
    skipped_lines: int = 0
    warnings: list[str] = field(default_factory=list)
    sequence_regions: dict[str, int] = field(default_factory=dict)

    @property
    def counts(self) -> dict[str, int]:
        counter = Counter(feature.type for feature in self.features)
        return {category: int(counter.get(category, 0)) for category in CATEGORIES}

    @property
    def contigs(self) -> list[str]:
        return sorted({feature.seqid for feature in self.features})

    def summary(self, sha256: str | None = None) -> dict[str, Any]:
        return {
            "path": self.path,
            "sha256": sha256,
            "format": self.format,
            "features": len(self.features),
            "counts": self.counts,
            "contigs": len(self.contigs),
            "contig_ids": self.contigs[:200],
            "contig_ids_truncated": len(self.contigs) > 200,
            "excluded_pseudo": self.excluded_pseudo,
            "skipped_lines": self.skipped_lines,
            "sequence_regions": len(self.sequence_regions),
            "warnings": list(self.warnings),
        }


@dataclass
class GenomeInfo:
    path: str
    sha256: str | None
    contigs: list[dict[str, Any]]
    warnings: list[str] = field(default_factory=list)

    @property
    def total_length(self) -> int:
        return sum(int(item["length"]) for item in self.contigs)

    def summary(self) -> dict[str, Any]:
        ids = [str(item["id"]) for item in self.contigs]
        return {
            "path": self.path,
            "sha256": self.sha256,
            "contigs": len(self.contigs),
            "total_length": self.total_length,
            "contig_ids": ids[:200],
            "contig_ids_truncated": len(ids) > 200,
            "warnings": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# Parseo GFF3 y GenBank
# ---------------------------------------------------------------------------
def _merge_gff_segments(features: Sequence[Feature]) -> tuple[list[Feature], int]:
    """Fusiona líneas GFF con el mismo ``ID`` (p. ej. CDS con frameshift).

    GFF3 puede representar un CDS en varias líneas; GenBank lo representa en un
    solo feature. Sin esta fusión, un frameshift aparecería como dos CDS en la
    consulta y como uno en la referencia.
    """
    order: list[tuple[str, str, str, str]] = []
    groups: dict[tuple[str, str, str, str], list[Feature]] = {}
    for feature in features:
        identifier = feature.attr("ID")
        if not identifier:
            key = ("", feature.type, feature.seqid, f"{feature.start}-{feature.end}-{feature.strand}")
        else:
            key = (identifier, feature.type, feature.seqid, feature.strand)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(feature)

    merged: list[Feature] = []
    merged_groups = 0
    for key in order:
        group = groups[key]
        if len(group) == 1:
            merged.append(group[0])
            continue
        merged_groups += 1
        first = group[0]
        merged.append(
            Feature(
                seqid=first.seqid,
                type=first.type,
                start=min(item.start for item in group),
                end=max(item.end for item in group),
                strand=first.strand,
                attrs=first.attrs,
                source=first.source,
            )
        )
    return merged, merged_groups


def parse_gff(
    path: str | Path,
    *,
    exclude_pseudo: bool = True,
    categories: Sequence[str] = CATEGORIES,
) -> AnnotationSet:
    """Lee GFF3 1-based inclusivo. Tolera gzip y líneas inválidas."""
    path = Path(path)
    wanted = {normalize_type(category) for category in categories}
    features: list[Feature] = []
    sequence_regions: dict[str, int] = {}
    excluded_pseudo = 0
    skipped_lines = 0
    warnings: list[str] = []

    with open_text(path) as fh:
        for line_number, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            sequence_regions[parts[1]] = int(parts[3])
                        except ValueError:
                            warnings.append(
                                f"línea {line_number}: sequence-region inválida"
                            )
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                skipped_lines += 1
                continue
            feature_type = normalize_type(columns[2])
            if feature_type not in wanted:
                continue
            try:
                start, end = int(columns[3]), int(columns[4])
            except ValueError:
                skipped_lines += 1
                continue
            if start > end:
                start, end = end, start
                warnings.append(
                    f"línea {line_number}: start > end, coordenadas intercambiadas"
                )
            attrs = _parse_attributes(columns[8])
            if exclude_pseudo and feature_type == "CDS" and _is_pseudo(attrs):
                excluded_pseudo += 1
                continue
            strand = columns[6] if columns[6] in {"+", "-", "."} else "."
            features.append(
                Feature(
                    seqid=columns[0],
                    type=feature_type,
                    start=start,
                    end=end,
                    strand=strand,
                    attrs=attrs,
                    source=columns[1],
                )
            )

    features, merged_groups = _merge_gff_segments(features)
    if merged_groups:
        warnings.append(
            f"{merged_groups} feature(s) GFF multi-segmento fusionados por ID; "
            "se comparan por su intervalo envolvente."
        )
    if skipped_lines:
        warnings.append(f"{skipped_lines} línea(s) GFF inválidas omitidas")
    return AnnotationSet(
        path=display_path(path) or str(path),
        format="gff3",
        features=features,
        excluded_pseudo=excluded_pseudo,
        skipped_lines=skipped_lines,
        warnings=warnings,
        sequence_regions=sequence_regions,
    )


def parse_genbank(
    path: str | Path,
    *,
    exclude_pseudo: bool = True,
    categories: Sequence[str] = CATEGORIES,
) -> AnnotationSet:
    """Lee GenBank/GBFF con Biopython. No se importa Biopython al cargar el módulo."""
    try:
        from Bio import SeqIO  # noqa: PLC0415 - dependencia opcional
    except ImportError as exc:  # pragma: no cover - depende del entorno
        raise DependencyError(
            "Biopython es necesario para leer referencias GenBank/GBFF. "
            "Instalar con: pip install biopython"
        ) from exc

    path = Path(path)
    wanted = {normalize_type(category) for category in categories}
    features: list[Feature] = []
    sequence_regions: dict[str, int] = {}
    excluded_pseudo = 0
    compound = 0
    skipped_records = 0
    warnings: list[str] = []

    with open_text(path) as fh:
        try:
            records = SeqIO.parse(fh, "genbank")
            for record in records:
                seqid = str(record.id or record.name or "").strip()
                if not seqid:
                    skipped_records += 1
                    continue
                try:
                    sequence_regions[seqid] = len(record.seq)
                except (TypeError, AttributeError):
                    pass
                for raw_feature in getattr(record, "features", []):
                    feature_type = normalize_type(str(raw_feature.type))
                    if feature_type not in wanted:
                        continue
                    location = raw_feature.location
                    try:
                        start = int(location.start) + 1
                        end = int(location.end)
                    except (TypeError, ValueError):
                        skipped_records += 1
                        continue
                    if start > end:
                        start, end = end, start
                    if len(getattr(location, "parts", [])) > 1:
                        compound += 1
                    strand_number = getattr(location, "strand", None)
                    strand = "+" if strand_number == 1 else "-" if strand_number == -1 else "."
                    qualifiers = getattr(raw_feature, "qualifiers", {}) or {}
                    attrs = {
                        "product": _first_qualifier(qualifiers, "product"),
                        "gene": _first_qualifier(qualifiers, "gene"),
                        "locus_tag": _first_qualifier(qualifiers, "locus_tag"),
                    }
                    pseudo = any(
                        key in qualifiers for key in ("pseudo", "pseudogene")
                    ) or _first_qualifier(qualifiers, "pseudo").lower() in {
                        "true",
                        "1",
                        "yes",
                    }
                    if pseudo:
                        attrs["pseudo"] = "true"
                    if exclude_pseudo and feature_type == "CDS" and pseudo:
                        excluded_pseudo += 1
                        continue
                    features.append(
                        Feature(
                            seqid=seqid,
                            type=feature_type,
                            start=max(1, start),
                            end=end,
                            strand=strand,
                            attrs=attrs,
                            source="genbank",
                        )
                    )
        except (ValueError, OSError) as exc:
            raise InputError(f"no se pudo parsear GenBank {path}: {exc}") from exc

    if compound:
        warnings.append(
            f"{compound} feature(s) con ubicación compuesta: {_COMPOUND_WARNING}"
        )
    if skipped_records:
        warnings.append(f"{skipped_records} registro/feature(s) omitidos")
    return AnnotationSet(
        path=display_path(path) or str(path),
        format="genbank",
        features=features,
        excluded_pseudo=excluded_pseudo,
        skipped_lines=skipped_records,
        warnings=warnings,
        sequence_regions=sequence_regions,
    )


def detect_format(path: str | Path, requested: str = "auto") -> str:
    requested = (requested or "auto").lower()
    if requested in {"gff", "gff3"}:
        return "gff3"
    if requested in {"genbank", "gbff", "gbk", "gb"}:
        return "genbank"
    if requested != "auto":
        raise InputError(f"formato no soportado: {requested}")

    path = Path(path)
    name = path.name.lower()
    if name.endswith((".gbff", ".gbk", ".gb", ".genbank")) or name.endswith(
        (".gbff.gz", ".gbk.gz", ".gb.gz", ".genbank.gz")
    ):
        return "genbank"
    if name.endswith((".gff", ".gff3", ".gtf")) or name.endswith(
        (".gff.gz", ".gff3.gz", ".gtf.gz")
    ):
        return "gff3"

    with open_text(path) as fh:
        for _ in range(40):
            line = fh.readline()
            if not line:
                break
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("LOCUS") or stripped == "ORIGIN":
                return "genbank"
            if stripped.startswith("##gff-version") or "\t" in stripped:
                return "gff3"
    raise InputError(
        f"no se pudo inferir el formato de {path}; usar format='gff3' o 'genbank'"
    )


def load_annotations(
    path: str | Path,
    *,
    format: str = "auto",
    exclude_pseudo: bool = True,
) -> AnnotationSet:
    """Carga GFF3 o GenBank/GBFF según ``format`` (auto por defecto)."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"no existe la anotación: {path}")
    resolved = detect_format(path, format)
    if resolved == "gff3":
        return parse_gff(path, exclude_pseudo=exclude_pseudo)
    return parse_genbank(path, exclude_pseudo=exclude_pseudo)


def read_fasta(path: str | Path, *, compute_sha256: bool = True) -> GenomeInfo:
    """Lee un FASTA (plano o gzip) y devuelve contigs con longitud."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"no existe el genoma FASTA: {path}")
    contigs: list[dict[str, Any]] = []
    warnings: list[str] = []
    current: dict[str, Any] | None = None
    seen: set[str] = set()

    with open_text(path) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    warnings.append(f"línea {line_number}: cabecera FASTA vacía")
                    current = {"id": f"contig_{len(contigs) + 1}", "description": "", "length": 0}
                else:
                    identifier = header.split()[0]
                    if identifier in seen:
                        warnings.append(f"ID FASTA duplicado: {identifier}")
                    seen.add(identifier)
                    current = {
                        "id": identifier,
                        "description": header,
                        "length": 0,
                    }
                contigs.append(current)
                continue
            if not line.strip():
                continue
            if current is None:
                raise InputError(
                    f"FASTA inválido: secuencia antes de la primera cabecera (línea {line_number})"
                )
            current["length"] += sum(1 for char in line if not char.isspace())

    if not contigs:
        raise InputError(f"el FASTA no contiene contigs: {path}")
    empty = [item["id"] for item in contigs if item["length"] == 0]
    if empty:
        warnings.append(f"contigs sin longitud: {', '.join(empty[:10])}")

    return GenomeInfo(
        path=display_path(path) or str(path),
        sha256=sha256_file(path) if compute_sha256 else None,
        contigs=contigs,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Normalización de seqid (NCBI vs cabecera FASTA vs salida de herramienta)
# ---------------------------------------------------------------------------
def _seqid_aliases(value: str) -> set[str]:
    value = (value or "").strip()
    if not value:
        return set()
    aliases = {value}
    first = value.split()[0]
    aliases.add(first)
    for token in re.split(r"[|\s]+", value):
        if token:
            aliases.add(token)
    for token in list(aliases):
        if _ACCESSION_RE.match(token):
            aliases.add(re.sub(r"\.\d+$", "", token))
    return aliases


@dataclass
class SeqidResolver:
    mapping: dict[str, str]
    ambiguous: list[str] = field(default_factory=list)

    def _lookup(self, seqid: str) -> str | None:
        for alias in sorted(_seqid_aliases(seqid), key=len, reverse=True):
            key = alias.lower()
            if key in self.mapping:
                return self.mapping[key]
            root = re.sub(r"\.\d+$", "", key)
            if root in self.mapping:
                return self.mapping[root]
        return None

    def has(self, seqid: str) -> bool:
        return self._lookup(seqid) is not None

    def resolve(self, seqid: str) -> str:
        found = self._lookup(seqid)
        return found if found is not None else (seqid or "").strip()

    def resolve_features(
        self, features: Sequence[Feature]
    ) -> tuple[list[Feature], int, list[str]]:
        remapped: list[Feature] = []
        changed = 0
        unknown: set[str] = set()
        for feature in features:
            if not self.has(feature.seqid):
                unknown.add(feature.seqid)
            target = self.resolve(feature.seqid)
            if target != feature.seqid:
                changed += 1
            remapped.append(feature.with_seqid(target))
        return remapped, changed, sorted(unknown)


def build_seqid_resolver(genome: GenomeInfo) -> SeqidResolver:
    mapping: dict[str, str] = {}
    ambiguous: list[str] = []
    for contig in genome.contigs:
        canonical = str(contig["id"])
        for text in (str(contig["id"]), str(contig["description"])):
            for alias in _seqid_aliases(text):
                key = alias.lower()
                if not key:
                    continue
                previous = mapping.get(key)
                if previous is not None and previous != canonical:
                    if alias not in ambiguous:
                        ambiguous.append(alias)
                    continue
                mapping[key] = canonical
    return SeqidResolver(mapping=mapping, ambiguous=ambiguous)


# ---------------------------------------------------------------------------
# Solapamiento y emparejamiento
# ---------------------------------------------------------------------------
def interval_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int
) -> int:
    """Longitud de intersección con coordenadas 1-based inclusivas."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def reciprocal_overlap(a: Feature, b: Feature) -> float:
    overlap = interval_overlap(a.start, a.end, b.start, b.end)
    if overlap <= 0:
        return 0.0
    return min(overlap / a.length, overlap / b.length)


def match_features(
    reference_features: Sequence[Feature],
    query_features: Sequence[Feature],
    *,
    threshold: float = 0.5,
) -> tuple[list[dict[str, Any]], list[Feature], list[Feature]]:
    """Emparejamiento 1:1 voraz por tipo, contig, hebra y solapamiento recíproco.

    El orden es determinista: mayor solapamiento, luego coordenadas. No es un
    alineamiento global ni un óptimo bipartito; se documenta como limitación.
    """
    if not 0.0 < threshold <= 1.0:
        raise ValueError("threshold debe estar en (0, 1]")

    ref_by_scope: dict[tuple[str, str, str], list[tuple[int, Feature]]] = defaultdict(list)
    query_by_scope: dict[tuple[str, str, str], list[tuple[int, Feature]]] = defaultdict(list)
    for index, feature in enumerate(reference_features):
        ref_by_scope[(feature.type, feature.seqid, feature.strand)].append((index, feature))
    for index, feature in enumerate(query_features):
        query_by_scope[(feature.type, feature.seqid, feature.strand)].append((index, feature))

    matches: list[dict[str, Any]] = []
    used_ref: set[int] = set()
    used_query: set[int] = set()

    for scope in sorted(set(ref_by_scope) | set(query_by_scope)):
        refs = ref_by_scope.get(scope, [])
        queries = query_by_scope.get(scope, [])
        if not refs or not queries:
            continue
        candidates: list[tuple[float, int, int, int, int, int, int]] = []
        for ref_index, ref in refs:
            for query_index, query in queries:
                score = reciprocal_overlap(ref, query)
                if score + 1e-12 < threshold or score <= 0.0:
                    continue
                candidates.append(
                    (
                        -score,
                        ref.start,
                        ref.end,
                        query.start,
                        query.end,
                        ref_index,
                        query_index,
                    )
                )
        candidates.sort()
        local_refs: set[int] = set()
        local_queries: set[int] = set()
        for neg_score, _, _, _, _, ref_index, query_index in candidates:
            if ref_index in local_refs or query_index in local_queries:
                continue
            local_refs.add(ref_index)
            local_queries.add(query_index)
            used_ref.add(ref_index)
            used_query.add(query_index)
            ref = reference_features[ref_index]
            query = query_features[query_index]
            matches.append(
                {
                    "reference": ref,
                    "query": query,
                    "reciprocal_overlap": _round(-neg_score),
                    "exact": ref.start == query.start and ref.end == query.end,
                    "start_delta": query.start - ref.start,
                    "end_delta": query.end - ref.end,
                }
            )

    unmatched_reference = [
        feature for index, feature in enumerate(reference_features) if index not in used_ref
    ]
    unmatched_query = [
        feature for index, feature in enumerate(query_features) if index not in used_query
    ]
    matches.sort(
        key=lambda item: (
            item["reference"].type,
            item["reference"].seqid,
            item["reference"].start,
            item["reference"].end,
            item["query"].start,
        )
    )
    return matches, unmatched_reference, unmatched_query


def _metrics(tp: int, fp: int, fn: int) -> dict[str, Any]:
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    denominator = 2 * tp + fp + fn
    f1 = _round((2 * tp) / denominator) if denominator else None
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def hypothetical_stats(features: Sequence[Feature]) -> dict[str, Any]:
    cds = [feature for feature in features if feature.type == "CDS"]
    with_product = [feature for feature in cds if feature.attr("product")]
    hypothetical = [feature for feature in with_product if is_hypothetical(feature.attr("product"))]
    return {
        "cds": len(cds),
        "with_product": len(with_product),
        "without_product": len(cds) - len(with_product),
        "hypothetical": len(hypothetical),
        "pct_hypothetical": _pct(len(hypothetical), len(with_product)),
    }


def _attribute_summary(
    matches: Sequence[Mapping[str, Any]],
    reference_features: Sequence[Feature],
    query_features: Sequence[Feature],
) -> dict[str, Any]:
    pairs = len(matches)
    both = Counter()
    literal = Counter()
    token_scores: list[float] = []

    for match in matches:
        ref = match["reference"]
        query = match["query"]
        for key in ("product", "gene", "locus_tag"):
            left, right = ref.attr(key), query.attr(key)
            if left and right:
                both[key] += 1
                if key == "product":
                    if normalize_text(left) == normalize_text(right):
                        literal[key] += 1
                    score = product_token_jaccard(left, right)
                    if score is not None:
                        token_scores.append(score)
                elif key == "locus_tag":
                    if left.lower() == right.lower():
                        literal[key] += 1
                elif normalize_text(left) == normalize_text(right):
                    literal[key] += 1

    def coverage(features: Sequence[Feature], key: str) -> dict[str, Any]:
        n = sum(1 for feature in features if feature.attr(key))
        return {"n": n, "pct": _pct(n, len(features))}

    return {
        "matched_pairs": pairs,
        "reference_coverage": {
            key: coverage(reference_features, key)
            for key in ("product", "gene", "locus_tag")
        },
        "query_coverage": {
            key: coverage(query_features, key)
            for key in ("product", "gene", "locus_tag")
        },
        "product": {
            "both_present": both["product"],
            "literal_matches": literal["product"],
            "literal_agreement": _ratio(literal["product"], both["product"]),
            "token_jaccard_mean": _mean(token_scores),
        },
        "gene": {
            "both_present": both["gene"],
            "literal_matches": literal["gene"],
            "literal_agreement": _ratio(literal["gene"], both["gene"]),
        },
        "locus_tag": {
            "both_present": both["locus_tag"],
            "literal_matches": literal["locus_tag"],
            "literal_agreement": _ratio(literal["locus_tag"], both["locus_tag"]),
        },
    }


def _strand_mismatches(
    reference_features: Sequence[Feature],
    query_features: Sequence[Feature],
    *,
    threshold: float,
) -> list[dict[str, Any]]:
    """Pares con solapamiento suficiente pero hebra distinta (diagnóstico)."""
    by_scope: dict[tuple[str, str], tuple[list[Feature], list[Feature]]] = {}
    refs_by_type: dict[tuple[str, str], list[Feature]] = defaultdict(list)
    queries_by_type: dict[tuple[str, str], list[Feature]] = defaultdict(list)
    for feature in reference_features:
        refs_by_type[(feature.type, feature.seqid)].append(feature)
    for feature in query_features:
        queries_by_type[(feature.type, feature.seqid)].append(feature)

    for key in set(refs_by_type) & set(queries_by_type):
        by_scope[key] = (refs_by_type[key], queries_by_type[key])

    candidates: list[tuple[float, str, str, Feature, Feature]] = []
    for key, (refs, queries) in by_scope.items():
        for ref in refs:
            for query in queries:
                if ref.strand == query.strand:
                    continue
                score = reciprocal_overlap(ref, query)
                if score + 1e-12 >= threshold:
                    candidates.append((-score, key[0], key[1], ref, query))
    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3].start, item[4].start))

    used_ref: set[tuple[str, str, str, int, int]] = set()
    used_query: set[tuple[str, str, str, int, int]] = set()
    output: list[dict[str, Any]] = []
    for neg_score, _, _, ref, query in candidates:
        ref_key = (ref.type, ref.seqid, ref.strand, ref.start, ref.end)
        query_key = (query.type, query.seqid, query.strand, query.start, query.end)
        if ref_key in used_ref or query_key in used_query:
            continue
        used_ref.add(ref_key)
        used_query.add(query_key)
        output.append(
            {
                "reference": ref.summary(),
                "query": query.summary(),
                "reciprocal_overlap": _round(-neg_score),
            }
        )
    return output


def _not_comparable(
    reference: AnnotationSet,
    query: AnnotationSet,
    reason: str,
    *,
    threshold: float,
    sequence_resolution: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "status": "not_comparable",
        "reason": reason,
        "threshold": threshold,
        "reference_counts": reference.counts,
        "query_counts": query.counts,
        "sequence_resolution": dict(sequence_resolution),
        "categories": {},
        "aggregate": {},
        "attributes": {},
        "hypothetical": {},
        "differences": {
            "counts": {
                "false_negatives": 0,
                "false_positives": 0,
                "strand_mismatches": 0,
                "boundary_shifts": 0,
                "exact_matches": 0,
            },
            "false_negatives": [],
            "false_positives": [],
            "strand_mismatches": [],
            "boundary_shifts": [],
            "truncated": False,
            "max_items": 0,
        },
        "warnings": [reason],
        "limitations": list(LIMITATIONS),
    }


def compare_annotations(
    reference: AnnotationSet,
    query: AnnotationSet,
    *,
    threshold: float = 0.5,
    resolver: SeqidResolver | None = None,
    max_differences: int = 100,
) -> dict[str, Any]:
    """Compara dos anotaciones cargadas y devuelve métricas y diferencias."""
    if not 0.0 < threshold <= 1.0:
        raise ValueError("threshold debe estar en (0, 1]")
    if max_differences < 0:
        raise ValueError("max_differences no puede ser negativo")

    warnings = list(reference.warnings) + list(query.warnings)
    reference_features = list(reference.features)
    query_features = list(query.features)
    ref_changed = query_changed = 0
    ref_unknown: list[str] = []
    query_unknown: list[str] = []

    if resolver is not None:
        reference_features, ref_changed, ref_unknown = resolver.resolve_features(reference_features)
        query_features, query_changed, query_unknown = resolver.resolve_features(query_features)
        if ref_unknown:
            warnings.append(
                "Seqids de referencia fuera del FASTA: " + ", ".join(ref_unknown[:20])
            )
        if query_unknown:
            warnings.append(
                "Seqids de la consulta fuera del FASTA: " + ", ".join(query_unknown[:20])
            )

    common_contigs = sorted(
        {feature.seqid for feature in reference_features}
        & {feature.seqid for feature in query_features}
    )
    sequence_resolution = {
        "enabled": resolver is not None,
        "reference_features_remapped": ref_changed,
        "query_features_remapped": query_changed,
        "reference_unmapped_seqids": ref_unknown[:100],
        "query_unmapped_seqids": query_unknown[:100],
        "reference_contigs": len({feature.seqid for feature in reference_features}),
        "query_contigs": len({feature.seqid for feature in query_features}),
        "common_contigs": len(common_contigs),
        "common_contig_ids": common_contigs[:100],
    }

    if not reference_features:
        return _not_comparable(
            reference,
            query,
            "La referencia no tiene features CDS/tRNA/rRNA evaluables.",
            threshold=threshold,
            sequence_resolution=sequence_resolution,
        )
    if not query_features:
        return _not_comparable(
            reference,
            query,
            "La anotación de la herramienta no tiene features CDS/tRNA/rRNA.",
            threshold=threshold,
            sequence_resolution=sequence_resolution,
        )
    if not common_contigs:
        return _not_comparable(
            reference,
            query,
            "No hay contigs en común entre referencia y herramienta tras normalizar seqid.",
            threshold=threshold,
            sequence_resolution=sequence_resolution,
        )

    matches, unmatched_reference, unmatched_query = match_features(
        reference_features,
        query_features,
        threshold=threshold,
    )

    categories: dict[str, Any] = {}
    for category in CATEGORIES:
        category_matches = [
            match for match in matches if match["reference"].type == category
        ]
        category_ref = [feature for feature in reference_features if feature.type == category]
        category_query = [feature for feature in query_features if feature.type == category]
        tp = len(category_matches)
        fp = sum(1 for feature in unmatched_query if feature.type == category)
        fn = sum(1 for feature in unmatched_reference if feature.type == category)
        deltas_start = [abs(match["start_delta"]) for match in category_matches]
        deltas_end = [abs(match["end_delta"]) for match in category_matches]
        categories[category] = {
            "reference_n": len(category_ref),
            "query_n": len(category_query),
            **_metrics(tp, fp, fn),
            "exact_matches": sum(1 for match in category_matches if match["exact"]),
            "boundary_shifted": sum(
                1 for match in category_matches if not match["exact"]
            ),
            "mean_reciprocal_overlap": _mean(
                match["reciprocal_overlap"] for match in category_matches
            ),
            "mean_abs_start_delta": _mean(deltas_start),
            "mean_abs_end_delta": _mean(deltas_end),
            "attributes": _attribute_summary(
                category_matches, category_ref, category_query
            ),
            "hypothetical": {
                "reference": hypothetical_stats(category_ref),
                "query": hypothetical_stats(category_query),
            },
        }

    tp_total = len(matches)
    fp_total = len(unmatched_query)
    fn_total = len(unmatched_reference)
    deltas_start_all = [abs(match["start_delta"]) for match in matches]
    deltas_end_all = [abs(match["end_delta"]) for match in matches]
    f1_values = [
        block["f1"] for block in categories.values() if block["f1"] is not None
    ]
    aggregate = {
        "reference_n": len(reference_features),
        "query_n": len(query_features),
        **_metrics(tp_total, fp_total, fn_total),
        "exact_matches": sum(1 for match in matches if match["exact"]),
        "boundary_shifted": sum(1 for match in matches if not match["exact"]),
        "mean_reciprocal_overlap": _mean(
            match["reciprocal_overlap"] for match in matches
        ),
        "mean_abs_start_delta": _mean(deltas_start_all),
        "mean_abs_end_delta": _mean(deltas_end_all),
        "macro_f1": _mean(f1_values),
        "hypothetical_reference": hypothetical_stats(reference_features),
        "hypothetical_query": hypothetical_stats(query_features),
    }

    strand_mismatches_full = _strand_mismatches(
        unmatched_reference,
        unmatched_query,
        threshold=threshold,
    )
    boundary_shifts_full = [
        {
            "type": match["reference"].type,
            "seqid": match["reference"].seqid,
            "strand": match["reference"].strand,
            "reference": match["reference"].summary(),
            "query": match["query"].summary(),
            "start_delta": match["start_delta"],
            "end_delta": match["end_delta"],
            "reciprocal_overlap": match["reciprocal_overlap"],
        }
        for match in matches
        if not match["exact"]
    ]
    false_negatives = [
        feature.summary()
        for feature in sorted(
            unmatched_reference, key=lambda f: (f.type, f.seqid, f.start, f.end)
        )
    ]
    false_positives = [
        feature.summary()
        for feature in sorted(
            unmatched_query, key=lambda f: (f.type, f.seqid, f.start, f.end)
        )
    ]

    fn_capped, fn_truncated = _cap(false_negatives, max_differences)
    fp_capped, fp_truncated = _cap(false_positives, max_differences)
    strand_capped, strand_truncated = _cap(strand_mismatches_full, max_differences)
    shifts_capped, shifts_truncated = _cap(boundary_shifts_full, max_differences)

    if tp_total == 0:
        warnings.append(
            "Ningún feature supera el umbral de solapamiento; revisar contigs, "
            "hebras o el umbral antes de interpretar los ceros."
        )
    if aggregate["query_n"] and categories["CDS"]["attributes"]["query_coverage"]["product"]["n"] == 0:
        warnings.append(
            "La anotación de la herramienta no trae 'product' en CDS; "
            "el % hypothetical de la consulta no es evaluable."
        )

    differences = {
        "counts": {
            "false_negatives": len(false_negatives),
            "false_positives": len(false_positives),
            "strand_mismatches": len(strand_mismatches_full),
            "boundary_shifts": len(boundary_shifts_full),
            "exact_matches": aggregate["exact_matches"],
        },
        "false_negatives": fn_capped,
        "false_positives": fp_capped,
        "strand_mismatches": strand_capped,
        "boundary_shifts": shifts_capped,
        "truncated": bool(
            fn_truncated or fp_truncated or strand_truncated or shifts_truncated
        ),
        "max_items": max_differences,
    }

    return {
        "status": "ok",
        "reason": "",
        "threshold": threshold,
        "reference_counts": reference.counts,
        "query_counts": query.counts,
        "sequence_resolution": sequence_resolution,
        "categories": categories,
        "aggregate": aggregate,
        "attributes": _attribute_summary(matches, reference_features, query_features),
        "hypothetical": {
            "reference": hypothetical_stats(reference_features),
            "query": hypothetical_stats(query_features),
        },
        "differences": differences,
        "warnings": warnings,
        "limitations": list(LIMITATIONS),
    }


def compare_gffs(
    reference_path: str | Path,
    query_path: str | Path,
    *,
    reference_format: str = "auto",
    query_format: str = "auto",
    threshold: float = 0.5,
    exclude_pseudo: bool = True,
    resolver: SeqidResolver | None = None,
    max_differences: int = 100,
) -> dict[str, Any]:
    """Atajo: carga dos archivos (GFF/GBFF) y los compara."""
    reference = load_annotations(
        reference_path, format=reference_format, exclude_pseudo=exclude_pseudo
    )
    query = load_annotations(
        query_path, format=query_format, exclude_pseudo=exclude_pseudo
    )
    return compare_annotations(
        reference,
        query,
        threshold=threshold,
        resolver=resolver,
        max_differences=max_differences,
    )


# ---------------------------------------------------------------------------
# Detección, versión y ejecución de herramientas
# ---------------------------------------------------------------------------
def detect_tool(name: str, which=None) -> dict[str, Any]:
    if name not in TOOL_SPECS:
        raise ValueError(f"herramienta desconocida: {name}")
    which = which or shutil.which
    spec = TOOL_SPECS[name]
    path = which(spec["binary"])
    return {
        "tool": name,
        "label": spec["label"],
        "binary": spec["binary"],
        "available": bool(path),
        "path": path,
        "reason": (
            ""
            if path
            else f"binario '{spec['binary']}' no encontrado en PATH (shutil.which)"
        ),
        "requires_database": bool(spec["requires_database"]),
    }


def detect_tools(which=None) -> dict[str, dict[str, Any]]:
    return {name: detect_tool(name, which=which) for name in TOOLS}


def probe_version(
    tool: str,
    binary_path: str,
    *,
    timeout: float = 30.0,
) -> dict[str, Any]:
    """Intenta registrar la versión de una herramienta; nunca falla el benchmark."""
    if tool not in TOOL_SPECS:
        raise ValueError(f"herramienta desconocida: {tool}")
    command = [str(binary_path), *TOOL_SPECS[tool]["version_args"]]
    try:
        proc = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        text = "\n".join(
            part.strip() for part in (proc.stdout, proc.stderr) if part and part.strip()
        )
        return {
            "status": "ok" if proc.returncode == 0 else "failed",
            "command": command,
            "returncode": proc.returncode,
            "text": _tail(text, 500),
        }
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "status": "failed",
            "command": command,
            "returncode": None,
            "text": _tail(str(exc), 500),
        }


def build_command(
    tool: str,
    binary_path: str,
    genome: str | Path,
    outdir: str | Path,
    *,
    prefix: str = "benchmark",
    threads: int = 1,
    bakta_db: str | Path | None = None,
    extra_args: Sequence[str] = (),
) -> list[str]:
    """Comando determinista y monocPU para cada herramienta."""
    if tool not in TOOL_SPECS:
        raise ValueError(f"herramienta desconocida: {tool}")
    if threads < 1:
        raise ValueError("threads debe ser >= 1")
    genome_s = str(Path(genome).resolve())
    out_s = str(Path(outdir).resolve())
    extras = [str(item) for item in extra_args]

    if tool == "bakta":
        command = [
            str(binary_path),
            "--output",
            out_s,
            "--prefix",
            prefix,
            "--threads",
            str(threads),
            "--force",
            "--skip-plot",
        ]
        if bakta_db:
            command.extend(["--db", str(Path(bakta_db).resolve())])
        command.extend(extras)
        command.append(genome_s)
        return command
    if tool == "prokka":
        return [
            str(binary_path),
            "--outdir",
            out_s,
            "--prefix",
            prefix,
            "--cpus",
            str(threads),
            "--force",
            "--quiet",
            *extras,
            genome_s,
        ]
    return [
        str(binary_path),
        "--genome",
        genome_s,
        "--out",
        out_s,
        "--cpu",
        str(threads),
        *extras,
    ]


def _max_rss_kb() -> int | None:
    try:
        import resource  # noqa: PLC0415 - opcional por plataforma

        return int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    except (ImportError, ValueError, OSError):
        return None


def execute_command(
    command: Sequence[str],
    *,
    timeout: float = 3600.0,
    cwd: str | Path | None = None,
) -> dict[str, Any]:
    """Ejecuta un comando externo y mide pared y RSS máximo de hijos."""
    start = time.perf_counter()
    rss_before = _max_rss_kb()
    try:
        completed = subprocess.run(
            list(command),
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        elapsed = time.perf_counter() - start
        rss_after = _max_rss_kb()
        return {
            "status": "ok" if completed.returncode == 0 else "failed",
            "command": [str(item) for item in command],
            "returncode": completed.returncode,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": rss_after,
            "max_rss_mib": round(rss_after / 1024.0, 2) if rss_after is not None else None,
            "max_rss_delta_kb": (
                max(0, rss_after - rss_before)
                if rss_after is not None and rss_before is not None
                else None
            ),
            "stdout_tail": _tail(completed.stdout),
            "stderr_tail": _tail(completed.stderr),
            "timeout_seconds": timeout,
            "rss_note": (
                "ru_maxrss de RUSAGE_CHILDREN (KiB en Linux); es el máximo "
                "acumulado de hijos, no un perfil instantáneo."
            ),
        }
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - start
        return {
            "status": "timeout",
            "command": [str(item) for item in command],
            "returncode": None,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": _max_rss_kb(),
            "max_rss_mib": None,
            "max_rss_delta_kb": None,
            "stdout_tail": _tail(exc.stdout),
            "stderr_tail": _tail(exc.stderr),
            "timeout_seconds": timeout,
            "rss_note": "Timeout durante la ejecución; RSS no es interpretable.",
        }
    except OSError as exc:
        elapsed = time.perf_counter() - start
        return {
            "status": "failed",
            "command": [str(item) for item in command],
            "returncode": None,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": _max_rss_kb(),
            "max_rss_mib": None,
            "max_rss_delta_kb": None,
            "stdout_tail": "",
            "stderr_tail": _tail(str(exc)),
            "timeout_seconds": timeout,
            "rss_note": "No se pudo iniciar el proceso externo.",
        }


def find_output_annotation(
    outdir: str | Path,
) -> tuple[Path | None, list[str]]:
    """Busca la anotación de salida, prefiriendo GFF3/GFF sobre GenBank."""
    outdir = Path(outdir)
    if not outdir.is_dir():
        return None, []
    suffixes = {".gff3": 0, ".gff": 1, ".gbff": 2, ".gbk": 3, ".gb": 4}
    candidates: list[tuple[int, str, Path]] = []
    for path in outdir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix not in suffixes:
            continue
        candidates.append((suffixes[suffix], path.as_posix(), path))
    candidates.sort(key=lambda item: (item[0], item[1]))
    names = [item[2].name for item in candidates]
    return (candidates[0][2] if candidates else None), names


def _copy_evidence(tool: str, save_dir: Path, outdir: Path, run: Mapping[str, Any]) -> None:
    destination = Path(save_dir) / tool
    destination.mkdir(parents=True, exist_ok=True)
    if outdir.is_dir():
        for path in outdir.rglob("*"):
            if path.is_file() and path.suffix.lower() in {
                ".gff3",
                ".gff",
                ".gbff",
                ".gbk",
                ".gb",
            }:
                target = destination / path.name
                shutil.copy2(path, target)
    (destination / "command.txt").write_text(
        " ".join(str(item) for item in run.get("command") or []) + "\n",
        encoding="utf-8",
    )
    (destination / "stdout.txt").write_text(
        str(run.get("stdout_tail") or ""), encoding="utf-8"
    )
    (destination / "stderr.txt").write_text(
        str(run.get("stderr_tail") or ""), encoding="utf-8"
    )


def evaluate_tool(
    tool: str,
    *,
    genome: GenomeInfo,
    reference: AnnotationSet,
    resolver: SeqidResolver | None,
    run: bool = False,
    query_gff: str | Path | None = None,
    which=None,
    threshold: float = 0.5,
    timeout: float = 3600.0,
    threads: int = 1,
    bakta_db: str | Path | None = None,
    prefix: str = "benchmark",
    exclude_pseudo: bool = True,
    probe_versions: bool = True,
    save_gffs: str | Path | None = None,
    max_differences: int = 100,
    extra_args: Sequence[str] = (),
) -> dict[str, Any]:
    """Evalúa una herramienta: precalculada, no ejecutada o ejecutada."""
    if tool not in TOOL_SPECS:
        raise ValueError(f"herramienta desconocida: {tool}")
    detection = detect_tool(tool, which=which)
    result: dict[str, Any] = {
        "tool": tool,
        "label": TOOL_SPECS[tool]["label"],
        "binary": TOOL_SPECS[tool]["binary"],
        "available": detection["available"],
        "binary_path": detection["path"],
        "status": "",
        "reason": "",
        "source": None,
        "version": None,
        "execution": None,
        "query_gff": None,
        "comparison": None,
    }

    # 1) Anotación precalculada: no requiere binario ni ejecución.
    if query_gff is not None:
        query_path = Path(query_gff)
        if not query_path.is_file():
            result.update(
                status="invalid_input",
                reason=f"no existe la anotación precalculada: {query_path}",
            )
            return result
        try:
            query_set = load_annotations(
                query_path, exclude_pseudo=exclude_pseudo
            )
        except (InputError, DependencyError) as exc:
            result.update(status="invalid_input", reason=str(exc))
            return result
        comparison = compare_annotations(
            reference,
            query_set,
            threshold=threshold,
            resolver=resolver,
            max_differences=max_differences,
        )
        result.update(
            status=comparison["status"],
            reason=comparison.get("reason", ""),
            source="precomputed",
            query_gff={
                **query_set.summary(sha256=sha256_file(query_path)),
            },
            comparison=comparison,
        )
        return result

    # 2) Sin binario: no se ejecuta nada y no se inventan métricas.
    if not detection["available"]:
        result.update(status="not_available", reason=detection["reason"])
        return result

    # 3) Binario disponible pero no se pidió ejecutar.
    if not run:
        result.update(
            status="not_run",
            reason="binario disponible, pero no se pidió --run",
        )
        return result

    if probe_versions:
        result["version"] = probe_version(
            tool,
            str(detection["path"]),
            timeout=min(30.0, max(5.0, float(timeout))),
        )

    with tempfile.TemporaryDirectory(prefix=f"external_benchmark_{tool}_") as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        outdir.mkdir(parents=True, exist_ok=True)
        command = build_command(
            tool,
            str(detection["path"]),
            genome.path,
            outdir,
            prefix=prefix,
            threads=threads,
            bakta_db=bakta_db,
            extra_args=extra_args,
        )
        execution = execute_command(command, timeout=timeout, cwd=tmpdir)
        result["execution"] = execution
        if save_gffs is not None:
            _copy_evidence(tool, Path(save_gffs), outdir, execution)

        if execution["status"] != "ok":
            result.update(
                status=execution["status"],
                reason=(
                    "la ejecución externa hizo timeout"
                    if execution["status"] == "timeout"
                    else f"la ejecución externa terminó con returncode={execution.get('returncode')}"
                ),
            )
            return result

        annotation_path, candidates = find_output_annotation(outdir)
        execution["annotation_candidates"] = candidates
        if annotation_path is None:
            result.update(
                status="failed",
                reason="la herramienta terminó con código 0, pero no se encontró GFF/GBFF de salida",
            )
            return result

        try:
            query_set = load_annotations(
                annotation_path, exclude_pseudo=exclude_pseudo
            )
        except (InputError, DependencyError) as exc:
            result.update(status="failed", reason=f"no se pudo leer la salida: {exc}")
            return result

        saved_path: str | None = None
        if save_gffs is not None:
            saved = Path(save_gffs) / tool / annotation_path.name
            if saved.is_file():
                saved_path = display_path(saved)
        result["query_gff"] = {
            **query_set.summary(sha256=sha256_file(annotation_path)),
            "path": saved_path or f"<temporal>/{annotation_path.name}",
        }
        comparison = compare_annotations(
            reference,
            query_set,
            threshold=threshold,
            resolver=resolver,
            max_differences=max_differences,
        )
        result.update(
            status=comparison["status"],
            reason=comparison.get("reason", ""),
            source="executed",
            comparison=comparison,
        )
        return result


# ---------------------------------------------------------------------------
# Orquestación y reporte
# ---------------------------------------------------------------------------
def _normalize_tools(tools: Sequence[str] | str) -> tuple[str, ...]:
    if isinstance(tools, str):
        if tools.strip().lower() == "all":
            return TOOLS
        selected = [item.strip() for item in tools.split(",") if item.strip()]
    else:
        selected = [str(item).strip() for item in tools if str(item).strip()]
    if not selected:
        raise ValueError("no se seleccionó ninguna herramienta")
    unknown = [item for item in selected if item not in TOOLS]
    if unknown:
        raise ValueError(f"herramientas desconocidas: {', '.join(unknown)}")
    return tuple(dict.fromkeys(selected))


def run_benchmark(
    genome_path: str | Path,
    reference_path: str | Path,
    *,
    tools: Sequence[str] | str = TOOLS,
    run: bool = False,
    query_gffs: Mapping[str, str | Path] | None = None,
    reference_format: str = "auto",
    exclude_pseudo: bool = True,
    threshold: float = 0.5,
    timeout: float = 3600.0,
    threads: int = 1,
    bakta_db: str | Path | None = None,
    prefix: str = "benchmark",
    seqid_aliasing: bool = True,
    probe_versions: bool = True,
    save_gffs: str | Path | None = None,
    max_differences: int = 100,
    which=None,
    extra_args: Mapping[str, Sequence[str]] | None = None,
    command: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Orquesta el benchmark completo y devuelve un dict serializable a JSON."""
    selected = _normalize_tools(tools)
    query_gffs = dict(query_gffs or {})
    unknown = [name for name in query_gffs if name not in TOOLS]
    if unknown:
        raise ValueError(f"query_gffs para herramientas desconocidas: {', '.join(unknown)}")
    for name, path in query_gffs.items():
        if not Path(path).is_file():
            raise InputError(f"query_gff de {name} no existe: {path}")

    genome = read_fasta(genome_path)
    reference = load_annotations(
        reference_path,
        format=reference_format,
        exclude_pseudo=exclude_pseudo,
    )
    if not reference.features:
        raise InputError(
            "la referencia no contiene CDS/tRNA/rRNA evaluables; revisar formato y archivo"
        )
    resolver = build_seqid_resolver(genome) if seqid_aliasing else None

    detection = detect_tools(which=which)
    results: dict[str, Any] = {}
    for name in selected:
        results[name] = evaluate_tool(
            name,
            genome=genome,
            reference=reference,
            resolver=resolver,
            run=run,
            query_gff=query_gffs.get(name),
            which=which,
            threshold=threshold,
            timeout=timeout,
            threads=threads,
            bakta_db=bakta_db,
            prefix=prefix,
            exclude_pseudo=exclude_pseudo,
            probe_versions=probe_versions,
            save_gffs=save_gffs,
            max_differences=max_differences,
            extra_args=(extra_args or {}).get(name, ()),
        )

    environment_tools: dict[str, Any] = {}
    for name in TOOLS:
        item = dict(detection[name])
        item["version"] = None
        if name in results and results[name].get("version"):
            version = results[name]["version"]
            item["version"] = version.get("text") if version.get("status") == "ok" else None
        environment_tools[name] = item

    not_available = [
        name for name in selected if results[name]["status"] == "not_available"
    ]
    warnings = list(genome.warnings) + list(reference.warnings)
    if resolver is not None and resolver.ambiguous:
        warnings.append(
            "Alias de seqid ambiguos (se conserva el primero): "
            + ", ".join(resolver.ambiguous[:20])
        )
    for name, result in results.items():
        if result["status"] == "not_available":
            warnings.append(f"{name}: {result['reason']}")
        if result["status"] == "failed":
            warnings.append(f"{name}: {result['reason']}")

    reference_summary = reference.summary(sha256=sha256_file(reference_path))
    return {
        "schema": SCHEMA,
        "generated_at": utc_now(),
        "scope": (
            "Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST "
            "contra una referencia GFF3/GBFF. CPU, sin descargas de datos."
        ),
        "command": [str(item) for item in command] if command else None,
        "parameters": {
            "genome": display_path(genome_path),
            "reference": display_path(reference_path),
            "reference_format": detect_format(reference_path, reference_format),
            "tools": list(selected),
            "run": bool(run),
            "query_gffs": {
                name: display_path(path) for name, path in query_gffs.items()
            },
            "overlap_threshold": threshold,
            "exclude_pseudo": bool(exclude_pseudo),
            "seqid_aliasing": bool(seqid_aliasing),
            "timeout_seconds": timeout,
            "threads": threads,
            "bakta_db": display_path(bakta_db),
            "probe_versions": bool(probe_versions),
            "max_differences": max_differences,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": f"{platform.system()} {platform.release()} {platform.machine()}",
            "cpu_count": os.cpu_count(),
            "tools": environment_tools,
            "available_tools": [
                name for name in TOOLS if environment_tools[name]["available"]
            ],
            "missing_tools": [
                name for name in TOOLS if not environment_tools[name]["available"]
            ],
        },
        "inputs": {
            "genome": genome.summary(),
            "reference": reference_summary,
        },
        "results": results,
        "not_available": not_available,
        "warnings": warnings,
        "limitations": list(LIMITATIONS),
        "reproduction": {
            "command": [str(item) for item in command] if command else None,
            "notes": [
                "Fijar las versiones de Bakta/Prokka/DFAST y registrar la digestión de la imagen.",
                "Conservar SHA256 de FASTA, referencia y GFF/GBFF de cada herramienta.",
                "Ejecutar con 1 CPU/hilo para comparar tiempos de forma consistente.",
                "Una comparación válida necesita un gold independiente y ejecuciones pareadas.",
            ],
        },
    }


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "sí" if value else "no"
    if isinstance(value, int):
        return str(value)
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def _cell(value: Any) -> str:
    text = str(value if value is not None else "—")
    return text.replace("|", "\\|").replace("\n", " ")


def _table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(_cell(header) for header in headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    rows = list(rows)
    if not rows:
        rows = [["—"] * len(headers)]
    for row in rows:
        lines.append("| " + " | ".join(_cell(cell) for cell in row) + " |")
    return "\n".join(lines)


def _status_label(status: str | None) -> str:
    return STATUS_LABELS.get(status or "", status or "—")


def _result_summary_row(name: str, result: Mapping[str, Any]) -> list[Any]:
    comparison = result.get("comparison") or {}
    categories = comparison.get("categories") or {}
    aggregate = comparison.get("aggregate") or {}
    execution = result.get("execution") or {}
    return [
        result.get("label", name),
        result.get("source") or "—",
        _status_label(result.get("status")),
        _fmt((categories.get("CDS") or {}).get("f1")),
        _fmt((categories.get("tRNA") or {}).get("f1")),
        _fmt((categories.get("rRNA") or {}).get("f1")),
        _fmt(aggregate.get("precision")),
        _fmt(aggregate.get("recall")),
        _fmt(aggregate.get("f1")),
        _fmt(execution.get("runtime_seconds"), 1),
    ]


def _category_table(comparison: Mapping[str, Any]) -> str:
    rows = []
    for category, block in (comparison.get("categories") or {}).items():
        rows.append(
            [
                category,
                block.get("reference_n"),
                block.get("query_n"),
                block.get("tp"),
                block.get("fp"),
                block.get("fn"),
                _fmt(block.get("precision")),
                _fmt(block.get("recall")),
                _fmt(block.get("f1")),
                block.get("exact_matches"),
                block.get("boundary_shifted"),
                _fmt(block.get("mean_reciprocal_overlap")),
            ]
        )
    return _table(
        [
            "tipo",
            "ref",
            "query",
            "TP",
            "FP",
            "FN",
            "precisión",
            "recall",
            "F1",
            "exactos",
            "borde distinto",
            "solape medio",
        ],
        rows,
    )


def _attribute_table(comparison: Mapping[str, Any]) -> str:
    rows = []
    for category, block in (comparison.get("categories") or {}).items():
        attrs = block.get("attributes") or {}
        rows.append(
            [
                category,
                attrs.get("matched_pairs"),
                (attrs.get("product") or {}).get("both_present"),
                _fmt((attrs.get("product") or {}).get("literal_agreement")),
                _fmt((attrs.get("product") or {}).get("token_jaccard_mean")),
                _fmt((attrs.get("gene") or {}).get("literal_agreement")),
                _fmt((attrs.get("locus_tag") or {}).get("literal_agreement")),
            ]
        )
    return _table(
        [
            "tipo",
            "pares",
            "product en ambos",
            "product literal",
            "product tokens",
            "gene literal",
            "locus literal",
        ],
        rows,
    )


def _hypothetical_table(comparison: Mapping[str, Any]) -> str:
    rows = []
    for category, block in (comparison.get("categories") or {}).items():
        hyp = block.get("hypothetical") or {}
        ref = hyp.get("reference") or {}
        query = hyp.get("query") or {}
        rows.append(
            [
                category,
                ref.get("cds"),
                ref.get("with_product"),
                ref.get("hypothetical"),
                _fmt(ref.get("pct_hypothetical"), 1),
                query.get("cds"),
                query.get("with_product"),
                query.get("hypothetical"),
                _fmt(query.get("pct_hypothetical"), 1),
            ]
        )
    return _table(
        [
            "tipo",
            "ref CDS",
            "ref con product",
            "ref hyp",
            "% ref",
            "query CDS",
            "query con product",
            "query hyp",
            "% query",
        ],
        rows,
    )


def _differences_markdown(comparison: Mapping[str, Any]) -> list[str]:
    differences = comparison.get("differences") or {}
    counts = differences.get("counts") or {}
    lines = [
        "#### Diferencias",
        "",
        "- Falsos negativos (referencia sin pareja): **{}**.".format(
            counts.get("false_negatives", 0)
        ),
        "- Falsos positivos (herramienta sin pareja): **{}**.".format(
            counts.get("false_positives", 0)
        ),
        "- Pares solo con hebra distinta: **{}**.".format(
            counts.get("strand_mismatches", 0)
        ),
        "- Pares emparejados con algún borde distinto: **{}**.".format(
            counts.get("boundary_shifts", 0)
        ),
        "",
    ]
    if differences.get("truncated"):
        lines.append(
            "Las listas de diferencias están truncadas a "
            f"{differences.get('max_items')} elementos por lista; los conteos son completos."
        )
        lines.append("")

    fn_rows = []
    for feature in differences.get("false_negatives") or []:
        fn_rows.append(
            [
                "FN",
                feature.get("type"),
                feature.get("seqid"),
                f"{feature.get('start')}..{feature.get('end')}",
                feature.get("strand"),
                feature.get("product") or feature.get("locus_tag") or "—",
            ]
        )
    fp_rows = []
    for feature in differences.get("false_positives") or []:
        fp_rows.append(
            [
                "FP",
                feature.get("type"),
                feature.get("seqid"),
                f"{feature.get('start')}..{feature.get('end')}",
                feature.get("strand"),
                feature.get("product") or feature.get("locus_tag") or "—",
            ]
        )
    shift_rows = []
    for item in differences.get("boundary_shifts") or []:
        ref = item.get("reference") or {}
        shift_rows.append(
            [
                item.get("type"),
                item.get("seqid"),
                f"{ref.get('start')}..{ref.get('end')}",
                f"{item.get('start_delta'):+d}" if item.get("start_delta") is not None else "—",
                f"{item.get('end_delta'):+d}" if item.get("end_delta") is not None else "—",
                _fmt(item.get("reciprocal_overlap")),
            ]
        )
    if fn_rows or fp_rows:
        lines += ["> Primeras diferencias no emparejadas", ""]
        lines += [
            _table(
                ["tipo", "feature", "seqid", "coordenadas", "hebra", "product/locus"],
                fn_rows + fp_rows,
            ),
            "",
        ]
    if shift_rows:
        lines += ["> Primeros pares con coordenadas no exactas", ""]
        lines += [
            _table(
                ["tipo", "seqid", "ref", "Δstart", "Δend", "solape"],
                shift_rows,
            ),
            "",
        ]
    return lines


def render_markdown(data: Mapping[str, Any]) -> str:
    parameters = data.get("parameters") or {}
    environment = data.get("environment") or {}
    inputs = data.get("inputs") or {}
    results = data.get("results") or {}
    genome = inputs.get("genome") or {}
    reference = inputs.get("reference") or {}

    lines: list[str] = [
        "# Benchmark externo CPU: Bakta, Prokka y DFAST",
        "",
        f"**Generado:** {data.get('generated_at', '—')}",
        "",
        f"**Esquema:** `{data.get('schema', SCHEMA)}`",
        "",
        f"**Alcance:** {data.get('scope', '—')}",
        "",
    ]
    if data.get("command"):
        lines += [
            "**Comando reproducible:**",
            "",
            "```bash",
            " ".join(str(item) for item in data["command"]),
            "```",
            "",
        ]
    lines += [
        "Este informe no declara superioridad entre herramientas. Las métricas "
        "comparan features contra la referencia local; no evalúan función, "
        "dominios ni rutas.",
        "",
        "## Entorno",
        "",
        f"- Python: `{environment.get('python', '—')}`",
        f"- Plataforma: `{environment.get('platform', '—')}`",
        f"- CPU lógicas: `{environment.get('cpu_count', '—')}`",
        f"- Umbral de solapamiento recíproco: **{_fmt(parameters.get('overlap_threshold'), 2)}**",
        f"- `--run`: **{_fmt(parameters.get('run'))}**",
        f"- Pseudogenes excluidos: **{_fmt(parameters.get('exclude_pseudo'))}**",
        f"- Normalización de seqid: **{_fmt(parameters.get('seqid_aliasing'))}**",
        "",
    ]

    tool_rows = []
    for name, item in (environment.get("tools") or {}).items():
        tool_rows.append(
            [
                item.get("label", name),
                f"`{item.get('binary')}`",
                _fmt(item.get("available")),
                item.get("version") or "—",
                item.get("path") or item.get("reason") or "—",
            ]
        )
    lines += [
        _table(["herramienta", "binario", "disponible", "versión", "PATH / razón"], tool_rows),
        "",
        f"- Herramientas disponibles: **{', '.join(environment.get('available_tools') or []) or 'ninguna'}**.",
        f"- Herramientas ausentes: **{', '.join(environment.get('missing_tools') or []) or 'ninguna'}**.",
        "",
        "## Entradas",
        "",
        f"- FASTA: `{genome.get('path', '—')}` "
        f"(SHA256 `{str(genome.get('sha256') or '—')[:16]}…`, "
        f"{genome.get('contigs', '—')} "
        f"{'contig' if genome.get('contigs') == 1 else 'contigs'}, "
        f"{genome.get('total_length', '—')} bp).",
        f"- Referencia: `{reference.get('path', '—')}` "
        f"(formato `{reference.get('format', '—')}`, "
        f"SHA256 `{str(reference.get('sha256') or '—')[:16]}…`).",
        f"- Conteos de referencia: CDS **{reference.get('counts', {}).get('CDS', '—')}**, "
        f"tRNA **{reference.get('counts', {}).get('tRNA', '—')}**, "
        f"rRNA **{reference.get('counts', {}).get('rRNA', '—')}**.",
        f"- CDS pseudo excluidos en referencia: **{reference.get('excluded_pseudo', 0)}**.",
        "",
        "## Resultados",
        "",
    ]
    summary_rows = [
        _result_summary_row(name, result) for name, result in results.items()
    ]
    lines += [
        _table(
            [
                "herramienta",
                "fuente",
                "estado",
                "F1 CDS",
                "F1 tRNA",
                "F1 rRNA",
                "precisión global",
                "recall global",
                "F1 global",
                "segundos",
            ],
            summary_rows,
        ),
        "",
        "La precisión/recall/F1 se calculan contra la referencia local. Si una "
        "herramienta no está instalada, su estado es `not_available` y no hay "
        "métricas que mostrar.",
        "",
    ]

    for name, result in results.items():
        lines += [f"## {result.get('label', name)}", ""]
        lines.append(f"- Estado: **{_status_label(result.get('status'))}**.")
        if result.get("source"):
            lines.append(f"- Fuente: `{result['source']}`.")
        if result.get("reason"):
            lines.append(f"- Razón: {result['reason']}.")
        if result.get("version"):
            version = result["version"]
            lines.append(
                f"- Versión: `{version.get('text') or '—'}` "
                f"(estado `{version.get('status')}`)."
            )
        execution = result.get("execution")
        if execution:
            lines.append(
                f"- Ejecución: **{_status_label(execution.get('status'))}**, "
                f"{_fmt(execution.get('runtime_seconds'), 2)} s, "
                f"RSS máx. {_fmt(execution.get('max_rss_mib'), 1)} MiB."
            )
            if execution.get("command"):
                lines.append(f"- Comando: `{' '.join(execution['command'])}`.")
            if execution.get("stderr_tail"):
                lines.append(
                    f"- stderr (cola): `{execution['stderr_tail'].strip()[-600:]}`"
                )
        query = result.get("query_gff") or {}
        if query:
            lines.append(
                f"- Anotación comparada: `{query.get('path', '—')}` "
                f"(formato `{query.get('format', '—')}`, "
                f"SHA256 `{str(query.get('sha256') or '—')[:16]}…`)."
            )
        comparison = result.get("comparison")
        if comparison:
            lines += ["", f"### Features ({result.get('label', name)})", ""]
            lines.append(_category_table(comparison))
            lines += ["", f"### Atributos ({result.get('label', name)})", ""]
            lines.append(_attribute_table(comparison))
            lines += ["", f"### Hipotéticas ({result.get('label', name)})", ""]
            lines.append(_hypothetical_table(comparison))
            lines += [""] + _differences_markdown(comparison)
            for warning in comparison.get("warnings") or []:
                lines.append(f"- Aviso: {warning}")
            lines.append("")
        else:
            lines += [
                "",
                "Sin comparación: no se ejecutó la herramienta o no produjo una "
                "anotación legible. No hay métricas para esta fila.",
                "",
            ]

    lines += ["## Herramientas no disponibles", ""]
    not_available = data.get("not_available") or []
    if not_available:
        for name in not_available:
            result = results.get(name) or {}
            lines.append(
                f"- **{result.get('label', name)}**: `not_available`; "
                f"{result.get('reason', 'sin binario')}. No se ejecutó y no se "
                "reportan métricas para esta herramienta."
            )
    else:
        lines.append("- Ninguna herramienta seleccionada quedó `not_available`.")

    lines += ["", "## Limitaciones", ""]
    for item in data.get("limitations") or []:
        lines.append(f"- {item}")

    lines += [
        "",
        "## Reproducción",
        "",
        "```bash",
        "python benchmarks/external/external_benchmark.py \\",
        f"  --genome {parameters.get('genome') or '<genoma.fna>'} \\",
        f"  --reference {parameters.get('reference') or '<referencia.gff3|gbff>'} \\",
        f"  --tools {','.join(parameters.get('tools') or list(TOOLS))} \\",
        "  --overlap-threshold "
        f"{_fmt(parameters.get('overlap_threshold'), 2)} \\",
        "  --json reports/external_benchmark_run/benchmark.json \\",
        "  --md reports/external_benchmark_run/benchmark.md",
        "```",
        "",
        "Añadir `--run` para ejecutar las herramientas instaladas. Para comparar "
        "salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin "
        "`--run`; el origen queda marcado como `precomputed`.",
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


def write_reports(
    data: Mapping[str, Any],
    *,
    json_path: str | Path,
    md_path: str | Path,
) -> tuple[Path, Path]:
    json_path = Path(json_path)
    md_path = Path(md_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    md_path.write_text(render_markdown(data), encoding="utf-8")
    return json_path, md_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_tool_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("usar TOOL=RUTA, por ejemplo prokka=out.gff")
    tool, path = value.split("=", 1)
    tool = tool.strip().lower()
    if tool not in TOOLS:
        raise argparse.ArgumentTypeError(f"herramienta desconocida: {tool}")
    if not path.strip():
        raise argparse.ArgumentTypeError("ruta de --query-gff vacía")
    return tool, Path(path)


def _parse_extra_arg(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("usar TOOL=ARG, por ejemplo dfast=--config")
    tool, argument = value.split("=", 1)
    tool = tool.strip().lower()
    if tool not in TOOLS:
        raise argparse.ArgumentTypeError(f"herramienta desconocida: {tool}")
    if not argument.strip():
        raise argparse.ArgumentTypeError("argumento extra vacío")
    return tool, argument


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark externo CPU contra Bakta, Prokka y DFAST. "
            "Sin --run solo detecta; las herramientas ausentes quedan not_available."
        )
    )
    parser.add_argument("--genome", type=Path, required=True, help="FASTA del genoma")
    parser.add_argument(
        "--reference",
        "--reference-gff",
        dest="reference",
        type=Path,
        required=True,
        help="Anotación de referencia GFF3 o GenBank/GBFF",
    )
    parser.add_argument(
        "--reference-format",
        choices=("auto", "gff3", "genbank"),
        default="auto",
        help="Formato de la referencia (auto por defecto)",
    )
    parser.add_argument(
        "--tools",
        "--tool",
        dest="tools",
        default="all",
        help="all o lista separada por comas: bakta,prokka,dfast",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="ejecuta cada herramienta disponible (sin descargar bases)",
    )
    parser.add_argument(
        "--query-gff",
        action="append",
        type=_parse_tool_path,
        default=[],
        metavar="TOOL=RUTA",
        help="compara un GFF/GBFF precalculado en vez de ejecutar la herramienta",
    )
    parser.add_argument("--bakta-db", type=Path, default=None, help="base Bakta local")
    parser.add_argument("--threads", type=int, default=1, help="CPUs/hilos (1 por defecto)")
    parser.add_argument(
        "--timeout", type=float, default=3600.0, help="timeout por herramienta en segundos"
    )
    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=0.5,
        help="solapamiento recíproco mínimo en (0, 1]",
    )
    parser.add_argument(
        "--keep-pseudo",
        action="store_true",
        help="no excluir CDS marcados como pseudo",
    )
    parser.add_argument(
        "--no-seqid-aliasing",
        action="store_true",
        help="no normalizar seqids contra las cabeceras del FASTA",
    )
    parser.add_argument(
        "--no-probe-versions",
        action="store_true",
        help="no ejecutar <tool> --version al usar --run",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        type=_parse_extra_arg,
        default=[],
        metavar="TOOL=ARG",
        help="argumento extra para una herramienta (repetible)",
    )
    parser.add_argument(
        "--save-gffs",
        type=Path,
        default=None,
        help="copia GFF/GBFF de salida y logs a este directorio",
    )
    parser.add_argument(
        "--max-differences",
        type=int,
        default=100,
        help="máximo de diferencias listadas por categoría en JSON/Markdown",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("reports/external_benchmark_run/benchmark.json"),
        help="ruta del JSON",
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=Path("reports/external_benchmark_run/benchmark.md"),
        help="ruta del Markdown",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="salir con 3 si alguna herramienta seleccionada no queda en estado ok",
    )
    parser.add_argument("--quiet", action="store_true", help="no imprimir resumen")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    raw_argv = list(argv) if argv is not None else sys.argv[1:]

    try:
        tools = _normalize_tools(args.tools)
    except ValueError as exc:
        parser.error(str(exc))

    if not 0.0 < args.overlap_threshold <= 1.0:
        parser.error("--overlap-threshold debe estar en (0, 1]")
    if args.threads < 1:
        parser.error("--threads debe ser >= 1")
    if args.timeout <= 0:
        parser.error("--timeout debe ser > 0")
    if args.max_differences < 0:
        parser.error("--max-differences no puede ser negativo")

    query_gffs: dict[str, Path] = {}
    for tool, path in args.query_gff:
        if tool in query_gffs:
            parser.error(f"--query-gff duplicado para {tool}")
        query_gffs[tool] = path

    extra_args: dict[str, list[str]] = defaultdict(list)
    for tool, argument in args.extra_arg:
        extra_args[tool].append(argument)

    bakta_db = args.bakta_db
    if bakta_db is None and os.environ.get("BAKTA_DB"):
        bakta_db = Path(os.environ["BAKTA_DB"])

    try:
        data = run_benchmark(
            args.genome,
            args.reference,
            tools=tools,
            run=args.run,
            query_gffs=query_gffs,
            reference_format=args.reference_format,
            exclude_pseudo=not args.keep_pseudo,
            threshold=args.overlap_threshold,
            timeout=args.timeout,
            threads=args.threads,
            bakta_db=bakta_db,
            seqid_aliasing=not args.no_seqid_aliasing,
            probe_versions=not args.no_probe_versions,
            save_gffs=args.save_gffs,
            max_differences=args.max_differences,
            extra_args=extra_args,
            command=[
                sys.executable,
                "benchmarks/external/external_benchmark.py",
                *raw_argv,
            ],
        )
    except BenchmarkError as exc:
        print(f"[external-benchmark] error: {exc}", file=sys.stderr)
        return 2

    write_reports(data, json_path=args.json, md_path=args.md)
    if not args.quiet:
        print(
            "[external-benchmark] herramientas disponibles: "
            f"{data['environment']['available_tools'] or 'ninguna'}"
        )
        print(
            "[external-benchmark] herramientas ausentes: "
            f"{data['not_available'] or 'ninguna'}"
        )
        for name, result in data["results"].items():
            comparison = result.get("comparison") or {}
            aggregate = comparison.get("aggregate") or {}
            print(
                f"[external-benchmark] {name}: {result['status']} "
                f"F1={aggregate.get('f1')} fuente={result.get('source')}"
            )
        print(f"[external-benchmark] escrito {args.json}")
        print(f"[external-benchmark] escrito {args.md}")

    if args.strict:
        bad = [
            name
            for name, result in data["results"].items()
            if result.get("status") != "ok"
        ]
        if bad:
            print(
                "[external-benchmark] --strict: estados no ok: " + ", ".join(bad),
                file=sys.stderr,
            )
            return 3
    return 0


__all__ = [
    "SCHEMA",
    "CATEGORIES",
    "TOOLS",
    "TOOL_SPECS",
    "AnnotationSet",
    "BenchmarkError",
    "DependencyError",
    "Feature",
    "GenomeInfo",
    "InputError",
    "SeqidResolver",
    "build_command",
    "build_seqid_resolver",
    "compare_annotations",
    "compare_gffs",
    "detect_format",
    "detect_tool",
    "detect_tools",
    "execute_command",
    "evaluate_tool",
    "find_output_annotation",
    "hypothetical_stats",
    "interval_overlap",
    "is_hypothetical",
    "load_annotations",
    "main",
    "match_features",
    "parse_genbank",
    "parse_gff",
    "probe_version",
    "product_token_jaccard",
    "read_fasta",
    "reciprocal_overlap",
    "render_markdown",
    "run_benchmark",
    "sha256_file",
    "write_reports",
]


if __name__ == "__main__":
    raise SystemExit(main())
