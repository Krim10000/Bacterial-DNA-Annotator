"""Colocación tentativa por vecinos de k-mers (CPU, sin red ni GPU).

Este módulo responde a una pregunta limitada: *¿a qué genoma de referencia
local se parece más una secuencia nueva según un MinHash de k-mers?* No
responde a "¿qué especie es?".

API pública
-----------

``sketch(sequence, k=21, size=512)``
    Devuelve un ``set[int]`` con los ``size`` hashes de k-mers canónicos más
    pequeños de una secuencia (ADN; se aceptan minúsculas y U como T).

``sketch_fasta(path, k=21, size=512)``
    Igual que ``sketch`` pero lee un FASTA (``.fna``, ``.fa``, ``.fasta``,
    ``.fas``, ``.ffn``), con o sin compresión gzip. Los k-mers no cruzan los
    límites entre registros.

``nearest_genomes(new_fasta_or_sequences, data_root, top_n=5)``
    Compara la consulta con todos los FASTA de ``<data_root>/data/raw`` y
    ``<data_root>/data/raw_excluded`` (también acepta ``data_root`` apuntando
    a ``data`` o a una de esas carpetas). Devuelve una lista de diccionarios
    con ``genome``, ``jaccard`` (estimado), ``lineage`` (solo si hay metadatos
    taxonómicos locales) y metadatos de la comparación.

``placement(new_sequences, data_root)``
    Envuelve ``nearest_genomes`` en un contrato serializable para la webapp:
    ``status`` (``ok``/``not_available``), ``nearest``, ``lineage``,
    ``confidence="low"``, ``method="kmer_nearest_not_ani"`` y ``limitation``.

Límites científicos (deliberados)
---------------------------------
- **No consulta NCBI, GTDB ni ninguna API.** Solo lee FASTA y JSON/TSV
  locales ya precalculados.
- **No calcula ANI real.** La similitud es un Jaccard estimado sobre hashes
  de k-mers; no mide identidad nucleotídica, cobertura ni sintenía.
- **No inventa especie ni linaje.** El linaje se copia de
  ``<data_root>/taxonomy/full/<genoma>.json`` (o de un TSV local de
  respaldo); un rango ausente queda ausente. Un vecino por k-mers **no** es
  una identificación taxonómica.
- El resultado siempre tiene ``confidence="low"`` y el método declarado
  ``"kmer_nearest_not_ani"``.

El backend vectorizado usa ``numpy`` solo si ya está instalado; si no está,
funciona igual con una implementación pura de la biblioteca estándar. Ambos
caminos comparten la misma función de mezcla de 64 bits y producen el mismo
sketch.

Rendimiento: el barrido completo de los 400 genomas tarda minutos en CPU
(el sketch de cada referencia se cachea en memoria, hasta 512 entradas, para
consultas posteriores dentro del mismo proceso). No se escribe caché en
disco ni se descarga nada.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import heapq
import json
import re
import zlib
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable, Iterator

try:  # numpy es una dependencia existente del proyecto, no una nueva.
    import numpy as _np
    from numpy.lib.stride_tricks import (
        sliding_window_view as _sliding_window_view,
    )
except Exception:  # pragma: no cover - entorno sin numpy o numpy roto
    _np = None
    _sliding_window_view = None  # type: ignore[assignment]

DEFAULT_K = 21
DEFAULT_SKETCH_SIZE = 512
DEFAULT_TOP_N = 5

METHOD = "kmer_nearest_not_ani"
CONFIDENCE = "low"
STATUS_OK = "ok"
STATUS_NOT_AVAILABLE = "not_available"

REFERENCE_DIRS = ("raw", "raw_excluded")
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fas", ".ffn", ".fsa")

LIMITATION = (
    "Similitud de k-mers (MinHash), no ANI: el Jaccard estimado no mide "
    "identidad nucleotídica ni cobertura. Un genoma vecino no es una "
    "identificación taxonómica. El linaje solo se copia de metadatos locales "
    "precalculados; no se consulta NCBI/GTDB y no se infiere especie nueva."
)

CONFIDENCE_REASON = (
    "El método no valida especie: propone el genoma de referencia local más "
    "parecido por k-mers y deja la decisión taxonómica a revisión humana."
)

_RANKS = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)
_PREFIX_TO_RANK = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
    "t": "strain",
}
_PLACEHOLDERS = {
    "",
    "unknown",
    "unclassified",
    "unranked",
    "na",
    "n/a",
    "none",
    "null",
    "not available",
    "not_available",
}

_MASK64 = (1 << 64) - 1
_C1 = 0x9E3779B97F4A7C15
_C2 = 0xBF58476D1CE4E5B9
_C3 = 0x94D049BB133111EB

# Traduce ADN a mayúsculas y U->T. Los bytes desconocidos se conservan y
# cortan la ventana de k-mers (no se adivinan bases).
_TRANS = bytes.maketrans(b"acgtuACGTU", b"ACGTTACGTT")
_WHITESPACE = b" \t\r\n\v\f"
_ACGT_SET = frozenset(b"ACGT")
_COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")

_ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")
_GENOMIC_SUFFIX_RE = re.compile(r"_genomic$", re.IGNORECASE)

_SKETCH_CACHE: "OrderedDict[tuple[Any, ...], set[int]]" = OrderedDict()
_SKETCH_CACHE_LIMIT = 512
_TABLE_CACHE: "OrderedDict[tuple[Any, ...], tuple[dict[str, dict[str, str]], str]]" = (
    OrderedDict()
)
_TABLE_CACHE_LIMIT = 8

_NUMPY_MIN_BASES = 512
_NUMPY_BLOCK_WINDOWS = 200_000


def _clean(value: Any) -> str:
    """Texto limpio; los valores que significan "sin dato" se vuelven ``""``."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in _PLACEHOLDERS:
        return ""
    return text


def _validate_params(k: int, size: int) -> tuple[int, int]:
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k debe ser un entero >= 1")
    if isinstance(size, bool) or not isinstance(size, int) or size < 1:
        raise ValueError("size debe ser un entero >= 1")
    return k, size


def _normalize_sequence(sequence: str | bytes | bytearray | memoryview) -> bytes:
    """Convierte una secuencia a bytes ACGT (U->T) sin espacios."""
    if isinstance(sequence, str):
        raw = sequence.encode("ascii", "replace")
    elif isinstance(sequence, (bytes, bytearray, memoryview)):
        raw = bytes(sequence)
    else:
        raise TypeError(
            "la secuencia debe ser str o bytes, no "
            f"{type(sequence).__name__}"
        )
    raw = raw.translate(_TRANS)
    return raw.translate(None, _WHITESPACE)


def _mix64(value: int) -> int:
    """Mezcla determinista de 64 bits (estilo splitmix64), sin aleatoriedad."""
    value = (value + _C1) & _MASK64
    value = ((value ^ (value >> 30)) * _C2) & _MASK64
    value = ((value ^ (value >> 27)) * _C3) & _MASK64
    return value ^ (value >> 31)


def _mix64_numpy(values: Any) -> Any:
    """Versión vectorizada de ``_mix64``; debe dar exactamente lo mismo."""
    values = values.copy()
    values += _np.uint64(_C1)
    values ^= values >> _np.uint64(30)
    values *= _np.uint64(_C2)
    values ^= values >> _np.uint64(27)
    values *= _np.uint64(_C3)
    values ^= values >> _np.uint64(31)
    return values


def _iter_kmer_hashes(sequence: bytes, k: int) -> Iterator[int]:
    """Genera hashes de k-mers canónicos de ``sequence`` (k <= 32)."""
    if k > 32:
        yield from _iter_long_kmer_hashes(sequence, k)
        return

    length = len(sequence)
    if length < k:
        return

    shift = 2 * (k - 1)
    mask = (1 << (2 * k)) - 1
    forward = 0
    reverse = 0
    valid = 0
    for base in sequence:
        code = _BASE_TO_CODE.get(base)
        if code is None:
            forward = 0
            reverse = 0
            valid = 0
            continue
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << shift)
        valid += 1
        if valid >= k:
            canonical = forward if forward <= reverse else reverse
            yield _mix64(canonical)


def _iter_long_kmer_hashes(sequence: bytes, k: int) -> Iterator[int]:
    """Camino lento para k > 32: blake2b del k-mer canónico como bytes."""
    for start in range(len(sequence) - k + 1):
        window = sequence[start:start + k]
        if not _ACGT_SET.issuperset(window):
            continue
        reverse = window.translate(_COMPLEMENT)[::-1]
        canonical = window if window <= reverse else reverse
        digest = hashlib.blake2b(canonical, digest_size=8).digest()
        yield int.from_bytes(digest, "big")


_BASE_TO_CODE = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}


def _sketch_sequence_python(sequence: bytes, k: int, size: int) -> set[int]:
    """Sketch bottom-k de una secuencia usando solo la biblioteca estándar."""
    heap: list[int] = []
    present: set[int] = set()
    for value in _iter_kmer_hashes(sequence, k):
        if len(heap) < size:
            if value not in present:
                present.add(value)
                heapq.heappush(heap, -value)
        elif value < -heap[0] and value not in present:
            removed = -heapq.heapreplace(heap, -value)
            present.discard(removed)
            present.add(value)
    return present


def _sketch_sequence_numpy(
    sequence: bytes, k: int, size: int
) -> set[int] | None:
    """Sketch bottom-k vectorizado. ``None`` si numpy no puede aplicarse."""
    if _np is None or _sliding_window_view is None or k > 32:
        return None
    length = len(sequence)
    if length < k:
        return set()

    codes = _BASE_CODE_LUT[_np.frombuffer(sequence, dtype=_np.uint8)]
    invalid = _np.flatnonzero(codes == 255)
    runs: list[tuple[int, int]] = []
    start = 0
    for index in invalid:
        index = int(index)
        if index - start >= k:
            runs.append((start, index))
        start = index + 1
    if length - start >= k:
        runs.append((start, length))
    if not runs:
        return set()

    candidates: list[Any] = []
    for run_start, run_end in runs:
        run_codes = codes[run_start:run_end]
        total_windows = run_end - run_start - k + 1
        for block_start in range(0, total_windows, _NUMPY_BLOCK_WINDOWS):
            block_end = min(
                block_start + _NUMPY_BLOCK_WINDOWS, total_windows
            )
            windows = _sliding_window_view(
                run_codes[block_start:block_end + k - 1], k
            )
            n_windows = block_end - block_start
            forward = _np.zeros(n_windows, dtype=_np.uint64)
            reverse = _np.zeros(n_windows, dtype=_np.uint64)
            for offset in range(k):
                forward <<= _np.uint64(2)
                forward |= windows[:, offset]
                reverse <<= _np.uint64(2)
                reverse |= 3 - windows[:, k - 1 - offset]
            canonical = _np.minimum(forward, reverse)
            hashes = _mix64_numpy(canonical)
            unique = _np.unique(hashes)
            candidates.append(unique[:size])

    if not candidates:
        return set()
    merged = _np.unique(_np.concatenate(candidates))
    return {int(value) for value in merged[:size].tolist()}


if _np is not None:
    _BASE_CODE_LUT = _np.full(256, 255, dtype=_np.uint8)
    for _index, _base in enumerate(b"ACGT"):
        _BASE_CODE_LUT[_base] = _index
else:  # pragma: no cover - solo para mantener el nombre definido
    _BASE_CODE_LUT = None  # type: ignore[assignment]


def _sketch_sequences(
    sequences: Iterable[bytes], k: int, size: int
) -> set[int]:
    """Une los sketches de varios registros tomando el bottom-k global."""
    present: set[int] = set()
    for sequence in sequences:
        if not sequence:
            continue
        if (
            _np is not None
            and k <= 32
            and len(sequence) >= _NUMPY_MIN_BASES
        ):
            candidate = _sketch_sequence_numpy(sequence, k, size)
            if candidate is None:  # pragma: no cover - defensivo
                candidate = _sketch_sequence_python(sequence, k, size)
        else:
            candidate = _sketch_sequence_python(sequence, k, size)
        if not candidate:
            continue
        if not present:
            present = candidate
        elif len(present) < size or min(candidate) < max(present):
            merged = present | candidate
            if len(merged) > size:
                merged = set(heapq.nsmallest(size, merged))
            present = merged
    return present


def sketch(
    sequence: str | bytes | bytearray | memoryview,
    k: int = DEFAULT_K,
    size: int = DEFAULT_SKETCH_SIZE,
) -> set[int]:
    """MinHash de k-mers canónicos de una secuencia.

    Devuelve como máximo ``size`` enteros de 64 bits: los hashes más pequeños
    de los k-mers distintos. Secuencias más cortas que ``k`` o sin bases ACGT
    suficientes producen un conjunto vacío. El resultado es determinista para
    la misma versión del módulo y no depende de ``PYTHONHASHSEED``.
    """
    k, size = _validate_params(k, size)
    normalized = _normalize_sequence(sequence)
    return _sketch_sequences([normalized], k, size)


def _open_fasta(path: Path) -> Any:
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return path.open("rb")


def _iter_fasta_records(path: Path) -> Iterator[bytes]:
    """Itera secuencias normalizadas de un FASTA (plano o gzip)."""
    with _open_fasta(path) as handle:
        parts: list[bytes] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if parts:
                    yield _normalize_sequence(b"".join(parts))
                    parts = []
            else:
                parts.append(line)
        if parts:
            yield _normalize_sequence(b"".join(parts))


def sketch_fasta(
    path: str | Path,
    k: int = DEFAULT_K,
    size: int = DEFAULT_SKETCH_SIZE,
) -> set[int]:
    """MinHash de un FASTA completo (varios registros, plano o ``.gz``)."""
    k, size = _validate_params(k, size)
    fasta_path = Path(path).expanduser()
    if not fasta_path.is_file():
        raise FileNotFoundError(f"no existe el FASTA: {fasta_path}")
    return _sketch_sequences(_iter_fasta_records(fasta_path), k, size)


def _is_fasta_name(name: str) -> bool:
    lowered = name.lower()
    if lowered.endswith(".gz"):
        lowered = lowered[:-3]
    return lowered.endswith(FASTA_SUFFIXES)


def _genome_id_from_path(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return _GENOMIC_SUFFIX_RE.sub("", name)


def _accession(genome: str) -> str:
    match = _ACCESSION_RE.match(genome)
    return match.group(1) if match else ""


def _resolve_layout(data_root: Path | str) -> tuple[Path, list[Path]]:
    """Acepta raíz del proyecto, ``data`` o una de las carpetas de FASTA."""
    root = Path(data_root).expanduser()
    try:
        root = root.resolve()
    except OSError:  # pragma: no cover - rutas raras
        pass

    if root.name in REFERENCE_DIRS:
        base = root.parent
        names = [root.name] + [
            name for name in REFERENCE_DIRS if name != root.name
        ]
        directories = [
            base / name for name in names if (base / name).is_dir()
        ]
        return base, directories

    candidates = [root]
    if (root / "data").is_dir():
        candidates.append(root / "data")
    for base in candidates:
        directories = [
            base / name for name in REFERENCE_DIRS if (base / name).is_dir()
        ]
        if directories:
            return base, directories
    return root, []


def _reference_files(
    data_root: Path | str,
) -> tuple[list[dict[str, Any]], Path, list[str]]:
    """Descubre FASTA en ``raw`` y ``raw_excluded`` (sin abrirlos aún)."""
    base, directories = _resolve_layout(data_root)
    if not directories and base.is_dir():
        # Comodidad para tests y para raíces que son directamente un panel.
        try:
            has_fasta = any(
                item.is_file() and _is_fasta_name(item.name)
                for item in base.iterdir()
            )
        except OSError:
            has_fasta = False
        if has_fasta:
            directories = [base]

    references: list[dict[str, Any]] = []
    labels: list[str] = []
    for directory in directories:
        label = directory.name
        if label not in labels:
            labels.append(label)
        try:
            entries = sorted(directory.iterdir(), key=lambda item: item.name)
        except OSError:  # pragma: no cover - permisos/carreras
            continue
        for path in entries:
            if not path.is_file() or not _is_fasta_name(path.name):
                continue
            genome = _genome_id_from_path(path)
            if not genome:
                continue
            try:
                relative = str(path.relative_to(base))
            except ValueError:  # pragma: no cover - otra raíz
                relative = str(path)
            references.append(
                {
                    "genome": genome,
                    "accession": _accession(genome),
                    "path": path,
                    "relative": relative,
                    "source": label,
                }
            )
    return references, base, labels


def _cached_sketch(path: Path, k: int, size: int) -> set[int]:
    try:
        stat = path.stat()
    except OSError:
        return sketch_fasta(path, k, size)
    key = (str(path), stat.st_mtime_ns, stat.st_size, k, size)
    cached = _SKETCH_CACHE.get(key)
    if cached is not None:
        _SKETCH_CACHE.move_to_end(key)
        return cached
    value = sketch_fasta(path, k, size)
    _SKETCH_CACHE[key] = value
    if len(_SKETCH_CACHE) > _SKETCH_CACHE_LIMIT:
        _SKETCH_CACHE.popitem(last=False)
    return value


def _jaccard_details(
    left: set[int], right: set[int], size: int
) -> tuple[float, int, int]:
    """Jaccard bottom-k: exacto si la unión cabe; si no, sobre la muestra."""
    if not left or not right:
        return 0.0, 0, 0
    union = left | right
    if len(union) <= size:
        sample = union
    else:
        sample = set(sorted(union)[:size])
    shared = len(left & right & sample)
    return shared / len(sample), shared, len(sample)


def _parse_lineage(text: Any) -> dict[str, str]:
    """``d__Bacteria;p__...`` -> ``{"domain": "Bacteria", ...}``."""
    cleaned = _clean(text)
    if not cleaned:
        return {}
    found: dict[str, str] = {}
    for token in cleaned.split(";"):
        token = token.strip()
        if len(token) < 4 or token[1:3] != "__":
            continue
        rank = _PREFIX_TO_RANK.get(token[0].lower())
        if rank is None:
            continue
        taxon = _clean(token[3:])
        if taxon and rank not in found:
            found[rank] = taxon
    return {rank: found[rank] for rank in _RANKS if rank in found}


def _ordered_ranks(ranks: dict[str, str]) -> dict[str, str]:
    return {rank: ranks[rank] for rank in _RANKS if ranks.get(rank)}


def _ranks_from_evidence(data: dict[str, Any]) -> dict[str, str]:
    """Respaldo mínimo para JSON antiguos con ``evidence.ncbi.ranks``."""
    evidence = data.get("evidence")
    if not isinstance(evidence, dict):
        return {}
    ncbi = evidence.get("ncbi")
    if not isinstance(ncbi, dict):
        return {}
    ranks = ncbi.get("ranks")
    if not isinstance(ranks, dict):
        return {}
    found: dict[str, str] = {}
    for rank in _RANKS:
        key = "species_ncbi" if rank == "species" else rank
        taxon = _clean(ranks.get(key))
        if taxon:
            found[rank] = taxon
    return found


def _lineage_from_json(
    path: Path, genome: str, data_dir: Path
) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    if not isinstance(data, dict):
        return None

    evidence = data.get("evidence")
    if not isinstance(evidence, dict):
        evidence = {}
    gtdb_evidence = evidence.get("gtdb")
    if not isinstance(gtdb_evidence, dict):
        gtdb_evidence = {}
    ncbi_evidence = evidence.get("ncbi")
    if not isinstance(ncbi_evidence, dict):
        ncbi_evidence = {}

    lineage_text = ""
    method = ""
    ranks: dict[str, str] = {}
    for candidate, label in (
        (data.get("gtdb_lineage"), "gtdb"),
        (gtdb_evidence.get("lineage"), "gtdb"),
        (data.get("ncbi_lineage"), "ncbi"),
        (ncbi_evidence.get("lineage"), "ncbi"),
    ):
        parsed = _parse_lineage(candidate)
        if parsed:
            ranks = parsed
            lineage_text = _clean(candidate)
            method = label
            break
    if not ranks:
        ranks = _ranks_from_evidence(data)
        if ranks:
            method = "ncbi"

    species = ranks.get("species", "")
    species_source = method if species else ""
    if not species:
        if method == "gtdb":
            candidates = (data.get("species_gtdb"), data.get("species"))
        elif method == "ncbi":
            candidates = (data.get("species_ncbi"), data.get("species"))
        else:
            candidates = (
                data.get("species"),
                data.get("species_ncbi"),
                data.get("species_gtdb"),
            )
        for candidate in candidates:
            cleaned = _clean(candidate)
            if cleaned:
                species = cleaned
                species_source = "gtdb" if method == "gtdb" else (
                    "ncbi" if method == "ncbi" else "species_json"
                )
                ranks["species"] = cleaned
                break

    if not ranks and not species:
        return None

    declared = _clean(data.get("confidence"))
    try:
        relative = str(path.relative_to(data_dir))
    except ValueError:  # pragma: no cover - ruta no estándar
        relative = str(path)
    return {
        "genome": genome,
        "lineage": lineage_text or None,
        "ranks": _ordered_ranks(ranks),
        "species": species or None,
        "species_source": species_source or None,
        "lineage_method": method or None,
        "declared_confidence": declared or None,
        "source": relative,
        "source_kind": "taxonomy_json",
        "partial": not bool(lineage_text),
    }


def _table_index(data_dir: Path) -> tuple[dict[str, dict[str, str]], str]:
    """Índice de metadatos locales de respaldo (``selected.tsv``, etc.)."""
    for relative in ("selected.tsv", "species_phylum.tsv"):
        path = data_dir / relative
        if not path.is_file():
            continue
        try:
            stat = path.stat()
        except OSError:  # pragma: no cover
            continue
        key = (str(path), stat.st_mtime_ns, stat.st_size)
        cached = _TABLE_CACHE.get(key)
        if cached is not None:
            _TABLE_CACHE.move_to_end(key)
            return cached
        index: dict[str, dict[str, str]] = {}
        try:
            with path.open(encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    accession = _clean(
                        row.get("assembly_accession") or row.get("accession")
                    )
                    if accession:
                        index[accession] = row
        except (OSError, UnicodeDecodeError, csv.Error):
            continue
        value = (index, relative)
        _TABLE_CACHE[key] = value
        if len(_TABLE_CACHE) > _TABLE_CACHE_LIMIT:
            _TABLE_CACHE.popitem(last=False)
        return value
    return {}, ""


def _lineage_from_table(
    data_dir: Path, genome: str, accession: str
) -> dict[str, Any] | None:
    index, relative = _table_index(data_dir)
    if not index:
        return None
    row = index.get(accession) or index.get(genome)
    if row is None:
        return None
    phylum = _clean(row.get("phylum"))
    species = _clean(row.get("species"))
    organism = _clean(row.get("organism_name"))
    if not phylum and not species:
        return None
    ranks: dict[str, str] = {}
    if phylum:
        ranks["phylum"] = phylum
    if species:
        ranks["species"] = species
    return {
        "genome": genome,
        "lineage": None,
        "ranks": _ordered_ranks(ranks),
        "species": species or None,
        "organism": organism or None,
        "species_source": "table" if species else None,
        "lineage_method": "table",
        "declared_confidence": None,
        "source": relative,
        "source_kind": "metadata_table",
        "partial": True,
    }


def _lineage_for(
    data_dir: Path, genome: str, accession: str
) -> dict[str, Any] | None:
    names = [genome]
    if accession and accession != genome:
        names.append(accession)
    for name in names:
        for folder in (data_dir / "taxonomy" / "full",
                       data_dir / "taxonomy"):
            path = folder / f"{name}.json"
            if path.is_file():
                lineage = _lineage_from_json(path, genome, data_dir)
                if lineage is not None:
                    return lineage
    return _lineage_from_table(data_dir, genome, accession)


def _is_probable_path_text(value: str) -> bool:
    if "/" in value or "\\" in value:
        return True
    lowered = value.strip().lower()
    if lowered.endswith(".gz"):
        lowered = lowered[:-3]
    return lowered.endswith(FASTA_SUFFIXES)


def _parse_fasta_text(text: str) -> list[bytes]:
    records: list[bytes] = []
    parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if parts:
                records.append(_normalize_sequence("".join(parts)))
                parts = []
        else:
            parts.append(line)
    if parts:
        records.append(_normalize_sequence("".join(parts)))
    return records


def _load_query(value: Any) -> list[bytes]:
    """Normaliza la consulta: ruta FASTA, texto FASTA o secuencia(s)."""
    if isinstance(value, Path):
        if not value.is_file():
            raise FileNotFoundError(f"no existe el FASTA: {value}")
        return list(_iter_fasta_records(value))
    if isinstance(value, (bytes, bytearray, memoryview)):
        return [_normalize_sequence(value)]
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith(">"):
            return _parse_fasta_text(value)
        path = Path(value).expanduser()
        try:
            is_file = path.is_file()
        except OSError:  # pragma: no cover
            is_file = False
        if is_file:
            return list(_iter_fasta_records(path))
        if _is_probable_path_text(value):
            raise FileNotFoundError(f"no existe el FASTA: {value}")
        return [_normalize_sequence(value)]
    if isinstance(value, Iterable):
        sequences: list[bytes] = []
        for item in value:
            if isinstance(item, Path):
                if not item.is_file():
                    raise FileNotFoundError(f"no existe el FASTA: {item}")
                sequences.extend(_iter_fasta_records(item))
            elif isinstance(item, (str, bytes, bytearray, memoryview)):
                sequences.append(_normalize_sequence(item))
            else:
                raise TypeError(
                    "cada elemento de la consulta debe ser str o bytes, no "
                    f"{type(item).__name__}"
                )
        return sequences
    raise TypeError(
        "la consulta debe ser una ruta FASTA, una secuencia o una lista de "
        f"secuencias, no {type(value).__name__}"
    )


def _nearest_impl(
    new_fasta_or_sequences: Any,
    data_root: Path | str,
    top_n: int,
    k: int,
    size: int,
) -> dict[str, Any]:
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n debe ser un entero >= 1")
    k, size = _validate_params(k, size)

    sequences = _load_query(new_fasta_or_sequences)
    query = _sketch_sequences(sequences, k, size)
    references, base, labels = _reference_files(data_root)

    diagnostics: dict[str, Any] = {
        "nearest": [],
        "query_hashes": len(query),
        "query_sequences": len(sequences),
        "n_references": len(references),
        "n_sketched": 0,
        "n_usable": 0,
        "skipped": [],
        "data_root": str(base),
        "directories": labels,
        "best_jaccard": 0.0,
    }
    if not query or not references:
        return diagnostics

    results: list[dict[str, Any]] = []
    for reference in references:
        try:
            reference_sketch = _cached_sketch(reference["path"], k, size)
        except (OSError, EOFError, zlib.error) as exc:
            diagnostics["skipped"].append(
                {"file": reference["relative"], "error": str(exc)}
            )
            continue
        diagnostics["n_sketched"] += 1
        if not reference_sketch:
            continue
        diagnostics["n_usable"] += 1
        jaccard, shared, sample_size = _jaccard_details(
            query, reference_sketch, size
        )
        results.append(
            {
                "genome": reference["genome"],
                "accession": reference["accession"],
                "jaccard": round(jaccard, 6),
                "shared_hashes": shared,
                "sketch_size": sample_size,
                "source": reference["source"],
                "file": reference["relative"],
            }
        )

    results.sort(key=lambda item: (-item["jaccard"], item["genome"]))
    top = results[:top_n]
    for item in top:
        lineage = _lineage_for(base, item["genome"], item["accession"])
        item["lineage"] = lineage
        item["lineage_text"] = lineage.get("lineage") if lineage else None
    diagnostics["nearest"] = top
    diagnostics["best_jaccard"] = top[0]["jaccard"] if top else 0.0
    return diagnostics


def nearest_genomes(
    new_fasta_or_sequences: Any,
    data_root: Path | str,
    top_n: int = DEFAULT_TOP_N,
    *,
    k: int = DEFAULT_K,
    size: int = DEFAULT_SKETCH_SIZE,
) -> list[dict[str, Any]]:
    """Genomas de referencia locales más parecidos a la consulta.

    Acepta una ruta a FASTA (plano o ``.gz``), una secuencia, texto FASTA o
    una lista de secuencias. Compara contra los FASTA de
    ``<data_root>/raw`` y ``<data_root>/raw_excluded`` (o de ``data/raw`` y
    ``data/raw_excluded`` si ``data_root`` es la raíz del proyecto).

    Devuelve una lista ordenada de diccionarios con ``genome``, ``jaccard``,
    ``lineage`` (diccionario local o ``None``), ``shared_hashes`` y
    ``sketch_size``. Si no hay referencias legibles o la consulta no produce
    k-mers, devuelve ``[]``; no lanza una excepción por ausencia de datos.
    """
    diagnostics = _nearest_impl(
        new_fasta_or_sequences, data_root, top_n, k, size
    )
    return diagnostics["nearest"]


def _not_available(
    reason: str,
    diagnostics: dict[str, Any] | None,
    top_n: int,
    k: int,
    size: int,
) -> dict[str, Any]:
    diagnostics = diagnostics or {}
    return {
        "status": STATUS_NOT_AVAILABLE,
        "nearest": [],
        "lineage": None,
        "lineage_text": None,
        "lineage_from": None,
        "confidence": CONFIDENCE,
        "confidence_reason": CONFIDENCE_REASON,
        "method": METHOD,
        "limitation": LIMITATION,
        "reason": reason,
        "top_n": top_n,
        "query": {
            "n_sequences": diagnostics.get("query_sequences", 0),
            "n_hashes": diagnostics.get("query_hashes", 0),
            "k": k,
            "sketch_size": size,
        },
        "references": {
            "root": diagnostics.get("data_root"),
            "directories": diagnostics.get("directories", []),
            "n_fastas": diagnostics.get("n_references", 0),
            "n_sketched": diagnostics.get("n_sketched", 0),
            "n_usable": diagnostics.get("n_usable", 0),
            "skipped": diagnostics.get("skipped", []),
        },
    }


def placement(
    new_sequences: Any,
    data_root: Path | str,
    top_n: int = DEFAULT_TOP_N,
    *,
    k: int = DEFAULT_K,
    size: int = DEFAULT_SKETCH_SIZE,
) -> dict[str, Any]:
    """Resultado de colocación tentativa listo para la webapp.

    ``status`` es ``ok`` solo si hay consulta, referencias legibles y al menos
    un k-mer compartido. En cualquier otro caso es ``not_available`` con una
    ``reason`` explícita. Nunca devuelve una especie inventada: ``lineage``
    es ``None`` o copia exacta de los metadatos locales del mejor vecino.
    """
    try:
        diagnostics = _nearest_impl(
            new_sequences, data_root, top_n, k, size
        )
    except (OSError, ValueError, TypeError) as exc:
        return _not_available(
            f"No se pudo procesar la consulta: {exc}",
            None,
            top_n,
            k,
            size,
        )

    if diagnostics["query_hashes"] == 0:
        return _not_available(
            "La consulta no produjo k-mers válidos "
            f"(k={k}; se requieren al menos {k} bases ACGT consecutivas).",
            diagnostics,
            top_n,
            k,
            size,
        )
    if diagnostics["n_references"] == 0:
        return _not_available(
            "No hay FASTA de referencia en data/raw ni data/raw_excluded "
            f"bajo {diagnostics['data_root']}.",
            diagnostics,
            top_n,
            k,
            size,
        )
    if diagnostics["n_usable"] == 0:
        return _not_available(
            "No se pudo leer ningún FASTA de referencia con k-mers válidos "
            f"(revisados {diagnostics['n_references']}).",
            diagnostics,
            top_n,
            k,
            size,
        )

    nearest = diagnostics["nearest"]
    if not nearest or nearest[0]["jaccard"] <= 0:
        return _not_available(
            "No hay k-mers compartidos con ninguna referencia local; no se "
            "propone un vecino.",
            diagnostics,
            top_n,
            k,
            size,
        )

    top = nearest[0]
    lineage = top.get("lineage")
    return {
        "status": STATUS_OK,
        "nearest": nearest,
        "lineage": lineage,
        "lineage_text": lineage.get("lineage") if lineage else None,
        "lineage_from": top["genome"],
        "confidence": CONFIDENCE,
        "confidence_reason": CONFIDENCE_REASON,
        "method": METHOD,
        "limitation": LIMITATION,
        "reason": "",
        "top_n": top_n,
        "best_jaccard": top["jaccard"],
        "query": {
            "n_sequences": diagnostics["query_sequences"],
            "n_hashes": diagnostics["query_hashes"],
            "k": k,
            "sketch_size": size,
        },
        "references": {
            "root": diagnostics["data_root"],
            "directories": diagnostics["directories"],
            "n_fastas": diagnostics["n_references"],
            "n_sketched": diagnostics["n_sketched"],
            "n_usable": diagnostics["n_usable"],
            "skipped": diagnostics["skipped"],
        },
    }


__all__ = [
    "CONFIDENCE",
    "DEFAULT_K",
    "DEFAULT_SKETCH_SIZE",
    "DEFAULT_TOP_N",
    "LIMITATION",
    "METHOD",
    "STATUS_NOT_AVAILABLE",
    "STATUS_OK",
    "nearest_genomes",
    "placement",
    "sketch",
    "sketch_fasta",
]
