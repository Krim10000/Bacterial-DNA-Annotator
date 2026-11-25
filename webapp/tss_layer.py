"""TSS/TTS precalculados para genomas nuevos (CPU, sin red ni GPU).

Esta capa **no predice** inicios ni finales de transcripción. Solo localiza
tablas ya calculadas, las normaliza y las deja listas para el navegador o para
un gen. Si no hay datos para el genoma pedido, devuelve ``not_available``: no
se inventa ninguna coordenada.

Fuentes locales reconocidas
---------------------------
- ``data/dbtbs/*_tss_tts.tsv`` y ``*_tss_tts.json`` (tabla combinada).
- ``data/regulation/<genoma>/`` con ``tss_tts.*``, ``tss.*``, ``tts.*`` o
  ``terminators.*`` (también un nivel de subdirectorio).
- Respaldo para las especies ya curadas de este repositorio:
  ``data/dbtbs/gt_tss.tsv`` + ``gt_terminators.tsv`` (*B. subtilis*, BSGatlas)
  y ``data/regulondb/tss_k12_union.tsv`` + ``TerminatorSet.tsv`` (*E. coli*
  K-12). Solo se usan si no hay una tabla específica del genoma.

Normalización
-------------
``strand`` se reduce a ``+``, ``-`` o ``.`` (desconocida); ``pos`` debe ser un
entero 1-based; ``source`` y ``class`` se mapean a un vocabulario corto y
documentado. ``class`` para TSS es el factor sigma (``sigmaA``, ``sigma70``…);
para TTS es la fuerza estructural (``strong``/``medium``/``weak``) cuando la
fuente la trae. La ubicación anómala (``intergenic``, ``antisense``…) se guarda
aparte en ``location``.

Los duplicados exactos (mismo tipo, contig, posición, hebra y fuente) se
eliminan y se cuentan. Dos fuentes distintas en la misma coordenada **no** se
borran: son corroboración, no duplicado.

API
---
``load(genome, data_root)``
    Devuelve ``{status, tss, tts, sources, ...}``. ``status`` es ``available``,
    ``partial`` (solo TSS o solo TTS) o ``not_available``.
``attach(features, tss, tts)``
    Añade a cada feature el TSS más cercano aguas arriba y el TTS más cercano
    aguas abajo (misma hebra y contig, dentro de ``max_distance``), además de
    contar los que caen dentro del feature. No crea coordenadas: solo enlaza
    registros ya cargados. Muta y devuelve ``features``.
``summary(loaded)`` o ``summary(tss, tts)``
    Conteos por tipo, clase, ubicación, fuente, hebra y cromosoma.

Nada de esto demuestra que un TSS sea un promotor funcional ni que un TTS
termine el transcrito de un gen concreto; son marcas de transcripción
precalculadas, con la procedencia que declare la fuente.
"""
from __future__ import annotations

import csv
import json
import re
import unicodedata
from bisect import bisect_left, bisect_right
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = ROOT / "data"

TSS = "TSS"
TTS = "TTS"

STATUS_AVAILABLE = "available"
STATUS_PARTIAL = "partial"
STATUS_NOT_AVAILABLE = "not_available"

#: Ventana por defecto para enlazar un TSS/TTS con un gen (pb). No es un
#: umbral biológico: es un límite conservador para no enlazar marcas lejanas.
MAX_ATTACH_DISTANCE = 600
MAX_DETAIL = 300

#: Procedencias canónicas. El valor de ``source`` de una fila se mapea aquí.
SOURCE_BSGATLAS = "BSGatlas"
SOURCE_DBTBS = "DBTBS"
SOURCE_REGULONDB = "RegulonDB"
SOURCE_HT = "RegulonDB alto rendimiento"
SOURCE_TRANSFER = "transferido K-12"
SOURCE_RULE = "rule_based"
SOURCE_BATTER = "BATTER"
SOURCE_PRECOMPUTED = "precalculado"
SOURCE_UNKNOWN = "desconocido"

LIMITATION = (
    "TSS/TTS precalculados y normalizados, no predichos por esta capa. Un TSS "
    "no es necesariamente un promotor funcional y un TTS no siempre termina el "
    "transcrito del gen vecino. La procedencia y la clase se conservan; sin "
    "datos locales el estado es not_available y no se inventan posiciones."
)

#: Alias de genoma → tokens normalizados. Sirve para que ``subtilis``, ``bsu``
#: o ``Bacillus subtilis 168`` encuentren la misma tabla precalculada.
GENOME_ALIASES: dict[str, frozenset[str]] = {
    "subtilis": frozenset({
        "subtilis", "bacillus", "bacillus_subtilis", "b_subtilis", "bsu",
        "bsub", "bacillus_subtilis_168", "bacillus_subtilis_subsp_subtilis",
        "168",
    }),
    "ecoli_k12": frozenset({
        "ecoli", "e_coli", "escherichia_coli", "ecoli_k12", "e_coli_k12",
        "k12", "k_12", "mg1655", "ecoli_mg1655",
    }),
}

_SOURCE_ALIASES: dict[str, str] = {
    "bsgatlas": SOURCE_BSGATLAS,
    "bsg": SOURCE_BSGATLAS,
    "gt": SOURCE_BSGATLAS,
    "gt_tss": SOURCE_BSGATLAS,
    "gt_terminators": SOURCE_BSGATLAS,
    "dbtbs": SOURCE_DBTBS,
    "regulondb": SOURCE_REGULONDB,
    "regulon_db": SOURCE_REGULONDB,
    "rdb": SOURCE_REGULONDB,
    "curado": SOURCE_REGULONDB,
    "curated": SOURCE_REGULONDB,
    "alto_rendimiento": SOURCE_HT,
    "high_throughput": SOURCE_HT,
    "ht": SOURCE_HT,
    "regulondb_ht": SOURCE_HT,
    "regulondb_alto_rendimiento": SOURCE_HT,
    "transferido_k12": SOURCE_TRANSFER,
    "transferido": SOURCE_TRANSFER,
    "transferido_desde_k12": SOURCE_TRANSFER,
    "transferred_k12": SOURCE_TRANSFER,
    "transfer_k12": SOURCE_TRANSFER,
    "rule_based": SOURCE_RULE,
    "reglas": SOURCE_RULE,
    "promoter_rule_based": SOURCE_RULE,
    "batter": SOURCE_BATTER,
    "precalculado": SOURCE_PRECOMPUTED,
    "precomputed": SOURCE_PRECOMPUTED,
    "local": SOURCE_PRECOMPUTED,
}

_LOCATION_ALIASES: dict[str, str] = {
    "intergenic": "intergenic",
    "intergenica": "intergenic",
    "inter_gene": "intergenic",
    "between_genes": "intergenic",
    "antisense": "antisense",
    "antisentido": "antisense",
    "antisense_transcript": "antisense",
    "inside_gene": "inside_gene",
    "inside": "inside_gene",
    "internal": "inside_gene",
    "intragenic": "inside_gene",
    "intragenica": "inside_gene",
    "within_gene": "inside_gene",
    "dentro_de_gen": "inside_gene",
    "gene_edge": "gene_edge",
    "edge_gene": "gene_edge",
    "gene_boundary": "gene_edge",
    "borde_de_gen": "gene_edge",
    "contig_edge": "contig_edge",
    "edge_contig": "contig_edge",
    "contig_boundary": "contig_edge",
    "borde_de_contig": "contig_edge",
}

_STRENGTH_ALIASES: dict[str, str] = {
    "a_fuerte": "strong",
    "a_fuerte_fuerte": "strong",
    "fuerte": "strong",
    "strong": "strong",
    "a": "strong",
    "b_media": "medium",
    "b_media_media": "medium",
    "media": "medium",
    "medium": "medium",
    "b": "medium",
    "c_debil": "weak",
    "c_debil_debil": "weak",
    "debil": "weak",
    "weak": "weak",
    "c": "weak",
    "batter": "batter",
}

_MECHANISM_ALIASES: dict[str, str] = {
    "rho_independent": "rho_independent",
    "rho_independiente": "rho_independent",
    "intrinsic": "rho_independent",
    "intrinseco": "rho_independent",
    "rho_dependent": "rho_dependent",
    "rho_dependiente": "rho_dependent",
}

_POS_FIELDS = (
    "pos", "position", "position_1based", "tss", "start", "leftEndPos",
    "2)leftEndPos", "coordinate", "coord", "middle", "center", "site",
)
_END_FIELDS = (
    "end", "stop", "rightEndPos", "3)rightEndPos", "fin", "end_pos",
)
_CONTIG_FIELDS = (
    "contig", "chrom", "chromosome", "chromosome_name", "seqid", "seq_id",
    "replicon", "refseq", "accession", "scaffold", "locus",
)
_STRAND_FIELDS = (
    "strand", "hebra", "orientation", "orientacion", "direccion", "4)strand",
)
_CLASS_FIELDS = (
    "class", "clase", "sigma", "strength", "fuerza", "class_type",
)
_LOCATION_FIELDS = (
    "location", "ubicacion", "anomaly", "anomalia", "position_class",
    "region_class",
)
_MECHANISM_FIELDS = (
    "mechanism", "mecanismo", "terminator_type", "7)type", "termination_type",
)
_SOURCE_FIELDS = (
    "source", "fuente", "provenance", "procedencia", "origin", "origen",
    "evidence_source",
)
_EVIDENCE_FIELDS = (
    "evidence", "evidencia", "comment", "comentario", "detail", "detalle",
    "description", "descripcion", "note", "nota",
)

_TYPE_ALIASES_TSS = frozenset({
    "tss", "transcription_start", "transcription_start_site",
    "inicio_transcripcion", "inicio_de_transcripcion", "start_site",
})
_TYPE_ALIASES_TTS = frozenset({
    "tts", "transcription_end", "transcription_end_site", "terminator",
    "terminador", "termination", "termination_site", "end_site",
})


# --------------------------------------------------------------------------
# Normalización de valores
# --------------------------------------------------------------------------
def _text(value: Any) -> str:
    """Texto limpio; ``None``/NaN/None-like se convierten en ``""``."""
    if value is None:
        return ""
    s = str(value).strip()
    if s.lower() in {"nan", "none", "null", "na", "n/a"}:
        return ""
    return s


def _ascii_key(value: Any) -> str:
    """Clave comparable: sin acentos, minúsculas y separadores ``_``."""
    s = _text(value).lower().replace("σ", "sigma").replace("ς", "sigma")
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


def normalize_strand(value: Any) -> str:
    """Normaliza hebra a ``+``, ``-`` o ``.`` (desconocida).

    Acepta ``+``/``-``, ``1``/``-1``, ``plus``/``minus``,
    ``forward``/``reverse`` y sus variantes en español.
    """
    raw = _text(value).lower()
    key = _ascii_key(raw)
    plus = {
        "+", "1", "1.0", "+1", "plus", "positive", "forward", "fwd", "for",
        "f", "direct", "directo", "watson", "mas", "sentido_directo",
    }
    minus = {
        "-", "-1", "-1.0", "minus", "negative", "reverse", "rev", "r",
        "complement", "complementario", "antisentido", "menos",
        "sentido_inverso",
    }
    if raw in plus or key in {
        "plus", "positive", "forward", "fwd", "for", "direct", "directo",
        "watson", "mas", "sentido_directo",
    }:
        return "+"
    if raw in minus or key in {
        "minus", "negative", "reverse", "rev", "complement",
        "complementario", "antisentido", "menos", "sentido_inverso",
    }:
        return "-"
    return "."


def normalize_pos(value: Any) -> int | None:
    """Entero 1-based o ``None``.

    Rechaza 0 y negativos a propósito: las tablas de este proyecto son
    1-based; convertir un 0 de BED sumando 1 sería inventar una coordenada.
    """
    if isinstance(value, bool):
        return None
    try:
        if isinstance(value, int):
            number = value
        elif isinstance(value, float):
            if not value.is_integer():
                return None
            number = int(value)
        else:
            text = _text(value).replace(",", "").replace(" ", "")
            if not text:
                return None
            parsed = float(text)
            if not parsed.is_integer():
                return None
            number = int(parsed)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if number >= 1 else None


def normalize_location(value: Any) -> str:
    """Ubicación anómala canónica o ``""``."""
    return _LOCATION_ALIASES.get(_ascii_key(value), "")


def normalize_mechanism(value: Any) -> str:
    """Mecanismo de terminación canónico o ``""``."""
    return _MECHANISM_ALIASES.get(_ascii_key(value), "")


def normalize_source(value: Any, fallback: str = "") -> str:
    """Procedencia canónica; una etiqueta desconocida se conserva tal cual."""
    raw = _text(value)
    key = _ascii_key(raw)
    if key in _SOURCE_ALIASES:
        return _SOURCE_ALIASES[key]
    if raw:
        return raw
    return fallback or SOURCE_UNKNOWN


def normalize_class(value: Any, kind: str) -> str:
    """Clase canónica según el tipo.

    - TSS: ``sigmaA``, ``sigma70``… (factor sigma), o ``unknown``.
    - TTS: ``strong``/``medium``/``weak``/``batter``, o ``unknown``.

    Una clase que en realidad sea una ubicación (``intergenic``…) devuelve
    ``unknown`` y se conserva en ``location``.
    """
    raw = _text(value)
    key = _ascii_key(raw)
    if not raw or key in {
        "unknown", "desconocido", "unresolved", "unresuelto", "sin_dato",
        "none", "na",
    }:
        return "unknown"
    if normalize_location(raw):
        return "unknown"
    if kind == TSS:
        if re.match(r"^rpo", raw, flags=re.IGNORECASE):
            # RpoD/RpoS son genes, no etiquetas de factor sigma; mapearlos a
            # "sigmaD"/"sigmaS" sería una interpretación falsa.
            return "unknown"
        label = re.sub(
            r"^(?:sigma|sig|σ)[\s:_-]*", "", raw, flags=re.IGNORECASE
        ).strip()
        if not label or _ascii_key(label) in {"unknown", "unresolved", ""}:
            return "unknown"
        if not re.fullmatch(r"[A-Za-z0-9]{1,4}", label):
            return "unknown"
        if len(label) == 1:
            label = label.upper()
        return "sigma" + label
    if kind == TTS:
        strength = _STRENGTH_ALIASES.get(key)
        if strength:
            return strength
        first = key.split("_")[0] if key else ""
        return _STRENGTH_ALIASES.get(first, "unknown")
    return "unknown"


# --------------------------------------------------------------------------
# Lectura de campos
# --------------------------------------------------------------------------
def _first(raw: Mapping[str, Any], keys: Iterable[str]) -> str:
    for key in keys:
        if key in raw:
            value = _text(raw[key])
            if value:
                return value
    return ""


def _clean_key(key: Any) -> str:
    return str(key).lstrip("\ufeff").strip()


def _record_type(raw: Mapping[str, Any], hint: str = "") -> str:
    kind = _ascii_key(_first(raw, ("type", "kind", "feature_type", "tipo")))
    if kind in _TYPE_ALIASES_TSS or kind.startswith("tss"):
        return TSS
    if kind in _TYPE_ALIASES_TTS or kind.startswith("tts") \
            or kind.startswith("terminator"):
        return TTS
    return hint if hint in (TSS, TTS) else ""


def _default_source(path: Path, genome_group: str = "") -> str:
    stem = _ascii_key(path.stem)
    if "bsgatlas" in stem or stem.startswith("gt_"):
        return SOURCE_BSGATLAS
    if "dbtbs" in stem:
        return SOURCE_DBTBS
    if ("regulondb" in stem or stem.startswith("terminatorset")
            or "tss_k12" in stem or "tts_k12" in stem):
        return SOURCE_REGULONDB
    if "transfer" in stem:
        return SOURCE_TRANSFER
    if "batter" in stem:
        return SOURCE_BATTER
    if "rule" in stem:
        return SOURCE_RULE
    if genome_group == "subtilis":
        # La tabla combinada de B. subtilis de este repositorio viene de
        # BSGatlas; si no trae columna source, esta es la procedencia correcta.
        return SOURCE_BSGATLAS
    return SOURCE_PRECOMPUTED


def _row_source(raw: Mapping[str, Any], default: str) -> str:
    explicit = _first(raw, _SOURCE_FIELDS)
    if explicit:
        return normalize_source(explicit, default)
    sigma = _first(raw, ("sigma", "sigma_factor", "factor_sigma"))
    if default == SOURCE_REGULONDB and _ascii_key(sigma).startswith("ht"):
        return SOURCE_HT
    return normalize_source(default, SOURCE_UNKNOWN)


def _detail(raw: Mapping[str, Any], *, class_raw: str, class_norm: str,
            kind: str) -> str:
    parts: list[str] = []
    if class_raw and class_norm == "unknown" and not normalize_location(class_raw):
        parts.append(f"class_raw={class_raw}")
    if kind == TSS:
        sigma = _first(raw, ("sigma", "sigma_factor", "factor_sigma"))
        if sigma and not class_norm.startswith("sigma"):
            parts.append(f"sigma={sigma}")
    labels = (
        ("evidence", _EVIDENCE_FIELDS),
        ("operon", ("operon", "operon_id")),
        ("gene", ("gene", "gene_name", "feature")),
        ("dist_gene", ("dist_gene", "distance_to_gene")),
        ("energy", ("energy", "energia", "dg", "delta_g")),
        ("parent", ("parent", "terminator_parent")),
        ("stem", ("stem",)),
        ("polyU", ("polyU", "poly_u", "polyu")),
        ("gc", ("gc",)),
        ("confidence", (
            "confidence", "confianza", "confidenceLevel", "11)confidenceLevel",
        )),
    )
    for label, fields in labels:
        value = _first(raw, fields)
        if value:
            parts.append(f"{label}={value}")
    text = "; ".join(parts)
    if len(text) > MAX_DETAIL:
        text = text[:MAX_DETAIL - 3] + "..."
    return text


def _normalize_record(raw: Mapping[str, Any], *, kind_hint: str,
                      genome: str, source_default: str
                      ) -> tuple[dict[str, Any] | None, str]:
    """Convierte una fila cruda en un registro canónico.

    Devuelve ``(registro, "")`` o ``(None, motivo)``.
    """
    kind = _record_type(raw, kind_hint)
    if not kind:
        return None, "sin tipo TSS/TTS"
    pos = normalize_pos(_first(raw, _POS_FIELDS))
    if pos is None:
        return None, "posicion invalida o no 1-based"
    end = normalize_pos(_first(raw, _END_FIELDS))
    if end is not None and end < pos:
        pos, end = end, pos

    contig = _first(raw, _CONTIG_FIELDS)
    strand = normalize_strand(_first(raw, _STRAND_FIELDS))
    class_raw = _first(raw, _CLASS_FIELDS)
    location_raw = _first(raw, _LOCATION_FIELDS)
    class_location = normalize_location(class_raw)
    location = normalize_location(location_raw) or class_location
    class_norm = normalize_class(class_raw, kind)
    if class_location and not normalize_location(location_raw):
        # La columna `class` era en realidad una ubicación.
        class_norm = "unknown"
    mechanism = normalize_mechanism(_first(raw, _MECHANISM_FIELDS))
    if not mechanism and kind == TTS:
        mechanism = normalize_mechanism(class_raw)
    source = _row_source(raw, source_default)
    return {
        "type": kind,
        "genome": genome,
        "contig": contig,
        "pos": pos,
        "end": end,
        "strand": strand,
        "source": source,
        "class": class_norm,
        "location": location,
        "mechanism": mechanism,
        "detail": _detail(raw, class_raw=class_raw, class_norm=class_norm,
                          kind=kind),
    }, ""


# --------------------------------------------------------------------------
# Lectura de archivos
# --------------------------------------------------------------------------
def _read_tsv(path: Path, *, genome: str, kind_hint: str,
              source_default: str) -> tuple[list[dict[str, Any]], int,
                                            list[str]]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as fh:
            lines = [line for line in fh
                     if line.strip() and not line.lstrip().startswith("#")]
    except OSError as exc:  # pragma: no cover - permiso/IO poco común
        return [], 0, [f"no se pudo leer {path.name}: {exc}"]
    if not lines:
        return [], 0, [f"{path.name} está vacío"]
    reader = csv.DictReader(lines, delimiter="\t")
    records: list[dict[str, Any]] = []
    rejected = 0
    for raw in reader:
        row = {_clean_key(k): v for k, v in raw.items() if k is not None}
        record, _why = _normalize_record(
            row, kind_hint=kind_hint, genome=genome,
            source_default=source_default,
        )
        if record is None:
            rejected += 1
            continue
        records.append(record)
    return records, rejected, []


def _json_items(data: Any) -> tuple[list[tuple[str, Any]], list[str]]:
    """Extrae ``(tipo_sugerido, item)`` de un JSON local.

    Acepta una lista de registros, ``{"tss": [...], "tts": [...]}`` o
    ``{"rows": [...]}``. Un JSON de resumen (solo conteos) no aporta ítems.
    """
    items: list[tuple[str, Any]] = []
    notes: list[str] = []

    def add_section(value: Any, hint: str) -> None:
        if isinstance(value, Mapping):
            for key in ("rows", "records", "positions", "data", "sites"):
                if isinstance(value.get(key), list):
                    value = value[key]
                    break
            else:
                if any(k in value for k in ("pos", "position")):
                    value = [value]
                else:
                    return
        if isinstance(value, list):
            items.extend((hint, item) for item in value)

    if isinstance(data, list):
        items.extend(("", item) for item in data)
    elif isinstance(data, Mapping):
        for key in ("tss", "TSS", "Tss", "tss_sites", "TSS_sites"):
            if key in data:
                add_section(data[key], TSS)
        for key in ("tts", "TTS", "Tts", "tts_sites", "TTS_sites"):
            if key in data:
                add_section(data[key], TTS)
        if not items:
            for key in ("rows", "records", "data", "features", "sites"):
                if isinstance(data.get(key), list):
                    items.extend(("", item) for item in data[key])
                    break
            if not items and any(k in data for k in ("pos", "position")):
                items.append(("", data))
        if not items:
            notes.append("JSON sin listas de posiciones (probablemente resumen)")
    else:
        notes.append("JSON no es una lista ni un objeto")
    return items, notes


def _read_json(path: Path, *, genome: str, kind_hint: str,
               source_default: str) -> tuple[list[dict[str, Any]], int,
                                             list[str]]:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        return [], 0, [f"no se pudo leer {path.name}: {exc}"]
    items, notes = _json_items(data)
    records: list[dict[str, Any]] = []
    rejected = 0
    for hint, item in items:
        if isinstance(item, Mapping):
            raw = {_clean_key(k): v for k, v in item.items()}
        elif isinstance(item, (int, float, str)) and not isinstance(item, bool):
            raw = {"pos": item}
        else:
            rejected += 1
            continue
        record, _why = _normalize_record(
            raw, kind_hint=hint or kind_hint, genome=genome,
            source_default=source_default,
        )
        if record is None:
            rejected += 1
            continue
        records.append(record)
    return records, rejected, notes


def _parse_file(path: Path, *, genome: str, kind_hint: str,
                source_default: str) -> tuple[list[dict[str, Any]], int,
                                              list[str]]:
    if path.suffix.lower() == ".json":
        return _read_json(path, genome=genome, kind_hint=kind_hint,
                          source_default=source_default)
    return _read_tsv(path, genome=genome, kind_hint=kind_hint,
                     source_default=source_default)


# --------------------------------------------------------------------------
# Descubrimiento de archivos
# --------------------------------------------------------------------------
def _resolve_data_root(data_root: str | Path | None) -> Path:
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    if (root / "dbtbs").is_dir() or (root / "regulation").is_dir():
        return root
    nested = root / "data"
    if (nested / "dbtbs").is_dir() or (nested / "regulation").is_dir():
        return nested
    return root


def _genome_context(genome: str) -> dict[str, Any]:
    """Claves y tokens del nombre de genoma, expandidos por alias."""
    raw = Path(str(genome)).name
    stem = raw
    pattern = re.compile(
        r"\.(?:fna|fa|fasta|fas|gbff|gbk|gff3?|gz|zip|bz2|xz)$", re.IGNORECASE)
    for _ in range(4):
        new = pattern.sub("", stem)
        if new == stem:
            break
        stem = new
    key = _ascii_key(stem)
    keys = {key} if key else set()
    for suffix in ("_genomic", "_genome", "_chromosome", "_dna", "_sequence"):
        if key.endswith(suffix):
            keys.add(key[:-len(suffix)].strip("_"))
    tokens: set[str] = set()
    for candidate in keys:
        tokens.update(part for part in candidate.split("_") if part)
    group = ""
    for name, aliases in GENOME_ALIASES.items():
        if keys & aliases or tokens & aliases:
            group = name
            keys |= set(aliases)
            for alias in aliases:
                tokens.update(part for part in alias.split("_") if part)
            break
    return {"keys": keys, "tokens": tokens, "group": group, "stem": key}


def _prefix_score(prefix: str, context: Mapping[str, Any],
                  allow_empty: bool = False) -> int | None:
    if not prefix:
        return 0 if allow_empty else None
    if prefix in context["keys"]:
        return 0
    tokens = set(prefix.split("_"))
    if tokens & context["tokens"]:
        return 2
    return None


def _file_kind(path: Path) -> str:
    """``TSS``, ``TTS``, ``both`` o ``""`` a partir del nombre."""
    stem = _ascii_key(path.stem)
    tokens = set(stem.split("_")) if stem else set()
    has_tss = ("tss" in tokens or "transcription_start" in stem
               or "start_site" in stem)
    has_tts = ("tts" in tokens or "terminator" in stem
               or "termination" in stem or "transcription_end" in stem)
    if has_tss and has_tts:
        return "both"
    if "transcription" in tokens:
        return "both"
    if has_tss:
        return TSS
    if has_tts:
        return TTS
    return ""


def _candidates(directory: Path) -> list[Path]:
    try:
        entries = sorted(directory.iterdir(), key=lambda p: p.name.lower())
    except OSError:  # pragma: no cover - directorio ilegible
        return []
    return [
        path for path in entries
        if path.is_file() and not path.name.startswith(".")
        and path.suffix.lower() in {".tsv", ".tab", ".txt", ".json"}
    ]


def _is_meta(path: Path) -> bool:
    tokens = set(_ascii_key(path.stem).split("_"))
    return bool(tokens & {"meta", "metadata", "summary", "resumen"})


def _pick_best(paths: Sequence[Path]) -> Path:
    priority = {
        "tss_tts": 0, "tss": 1, "tts": 1, "terminators": 1,
        "terminatorset": 1, "tss_tts_meta": 4,
    }

    def key(path: Path) -> tuple[int, int, int, str]:
        tsv = 0 if path.suffix.lower() in {".tsv", ".tab", ".txt"} else 1
        stem = _ascii_key(path.stem)
        return (tsv, priority.get(stem, 3), len(path.name), path.name.lower())

    return min(paths, key=key)


def _combined_dbtbs(directory: Path,
                    context: Mapping[str, Any]) -> Path | None:
    grouped: dict[str, list[tuple[int, Path]]] = {}
    for path in _candidates(directory):
        if _file_kind(path) != "both":
            continue
        core = _ascii_key(path.stem)
        prefix = core
        for suffix in ("_tss_tts", "_tts_tss"):
            if prefix.endswith(suffix):
                prefix = prefix[:-len(suffix)]
                break
        if core in {"tss_tts", "tts_tss"}:
            prefix = ""
        score = _prefix_score(prefix, context)
        if score is None:
            continue
        grouped.setdefault(prefix, []).append((score, path))
    if not grouped:
        return None
    best_prefix = min(
        grouped,
        key=lambda prefix: (min(score for score, _ in grouped[prefix]), prefix),
    )
    return _pick_best([path for _, path in grouped[best_prefix]])


def _find_regulation_dir(root: Path,
                         context: Mapping[str, Any]) -> Path | None:
    regulation = root / "regulation"
    if not regulation.is_dir():
        return None
    for candidate in sorted(context["keys"]):
        path = regulation / candidate
        if path.is_dir():
            return path
    for path in sorted(regulation.iterdir(), key=lambda p: p.name.lower()):
        if path.is_dir() and not path.name.startswith(".") \
                and _ascii_key(path.name) in context["keys"]:
            return path
    return None


def _regulation_selection(directory: Path,
                          genome_group: str = "") -> list[tuple[Path, str, str]]:
    """Elige una tabla combinada o una de TSS y otra de TTS del genoma."""
    files = [path for path in _candidates(directory) if _file_kind(path)]
    for child in sorted(directory.iterdir(), key=lambda p: p.name.lower()):
        if child.is_dir() and not child.name.startswith("."):
            files.extend(path for path in _candidates(child) if _file_kind(path))
    if not files:
        return []
    combined = [path for path in files if _file_kind(path) == "both"]
    real_combined = [path for path in combined if not _is_meta(path)]
    if real_combined:
        chosen = _pick_best(real_combined)
        return [(chosen, "", _default_source(chosen, genome_group))]
    selection: list[tuple[Path, str, str]] = []
    for kind in (TSS, TTS):
        matches = [
            path for path in files
            if _file_kind(path) == kind and not _is_meta(path)
        ]
        if matches:
            chosen = _pick_best(matches)
            selection.append(
                (chosen, kind, _default_source(chosen, genome_group)))
    return selection


def _legacy_selection(root: Path,
                      context: Mapping[str, Any]) -> list[tuple[Path, str, str]]:
    """Respaldo explícito para las dos especies curadas del repositorio."""
    selection: list[tuple[Path, str, str]] = []
    if context["group"] == "subtilis":
        tss = root / "dbtbs" / "gt_tss.tsv"
        tts = root / "dbtbs" / "gt_terminators.tsv"
        if tss.is_file():
            selection.append((tss, TSS, SOURCE_BSGATLAS))
        if tts.is_file():
            selection.append((tts, TTS, SOURCE_BSGATLAS))
    elif context["group"] == "ecoli_k12":
        tss = root / "regulondb" / "tss_k12_union.tsv"
        tts = root / "regulondb" / "TerminatorSet.tsv"
        if tss.is_file():
            selection.append((tss, TSS, SOURCE_REGULONDB))
        if tts.is_file():
            selection.append((tts, TTS, SOURCE_REGULONDB))
    return selection


# --------------------------------------------------------------------------
# Deduplicación y resumen
# --------------------------------------------------------------------------
def _sort_key(record: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        str(record.get("contig", "")), int(record.get("pos", 0) or 0),
        str(record.get("strand", ".")), str(record.get("source", "")),
        int(record.get("end", 0) or 0), str(record.get("class", "")),
    )


def _dedupe(records: Sequence[Mapping[str, Any]]
            ) -> tuple[list[dict[str, Any]], int]:
    seen: set[tuple[Any, ...]] = set()
    kept: list[dict[str, Any]] = []
    duplicates = 0
    for record in sorted(records, key=_sort_key):
        key = (
            record.get("type"), str(record.get("contig", "")),
            int(record.get("pos", 0) or 0), str(record.get("strand", ".")),
            str(record.get("source", "")),
        )
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)
        kept.append(dict(record))
    return kept, duplicates


def _shared_positions(tss: Sequence[Mapping[str, Any]],
                      tts: Sequence[Mapping[str, Any]]) -> int:
    shared = 0
    for records in (tss, tts):
        by_position: dict[tuple[Any, ...], set[str]] = {}
        for record in records:
            key = (record.get("type"), record.get("contig"),
                   record.get("pos"), record.get("strand"))
            by_position.setdefault(key, set()).add(
                str(record.get("source", "")))
        shared += sum(1 for sources in by_position.values() if len(sources) > 1)
    return shared


def _counts(records: Sequence[Mapping[str, Any]],
            field: str) -> dict[str, int]:
    counter = Counter(
        str(record.get(field) or "") for record in records
        if str(record.get(field) or "")
    )
    return {key: counter[key] for key in sorted(counter)}


def _positions(records: Sequence[Mapping[str, Any]]) -> list[int]:
    values: list[int] = []
    for record in records:
        pos = normalize_pos(record.get("pos"))
        if pos is not None:
            values.append(pos)
    return values


# --------------------------------------------------------------------------
# API pública
# --------------------------------------------------------------------------
def load(genome: str, data_root: str | Path | None = None) -> dict[str, Any]:
    """Carga TSS/TTS precalculados para ``genome``.

    Busca primero en ``data/regulation/<genoma>/`` (lo más específico), luego
    ``data/dbtbs/*_tss_tts.tsv|json`` y por último los respaldos curados de
    *B. subtilis* / *E. coli* K-12. Nunca descarga ni predice.

    Devuelve un dict con ``status`` (``available``/``partial``/
    ``not_available``), ``tss``, ``tts``, ``sources``, conteos y ``reason``.
    """
    root = _resolve_data_root(data_root)
    context = _genome_context(genome)
    genome_name = str(genome)
    searched: list[str] = []
    selection: list[tuple[Path, str, str]] = []

    regulation_dir = _find_regulation_dir(root, context)
    if regulation_dir is not None:
        searched.append(str(regulation_dir))
        selection = _regulation_selection(
            regulation_dir, context["group"])

    if not selection:
        dbtbs = root / "dbtbs"
        if dbtbs.is_dir():
            searched.append(str(dbtbs / "*_tss_tts.tsv|json"))
            combined = _combined_dbtbs(dbtbs, context)
            if combined is not None:
                selection = [(combined, "",
                              _default_source(combined, context["group"]))]

    if not selection:
        selection = _legacy_selection(root, context)

    tss: list[dict[str, Any]] = []
    tts: list[dict[str, Any]] = []
    sources: list[dict[str, Any]] = []
    notes: list[str] = []
    rejected = 0
    for path, hint, source_default in selection:
        records, file_rejected, file_notes = _parse_file(
            path, genome=genome_name, kind_hint=hint,
            source_default=source_default,
        )
        rejected += file_rejected
        notes.extend(file_notes)
        n_tss = sum(1 for record in records if record["type"] == TSS)
        n_tts = sum(1 for record in records if record["type"] == TTS)
        tss.extend(record for record in records if record["type"] == TSS)
        tts.extend(record for record in records if record["type"] == TTS)
        sources.append({
            "path": str(path),
            "format": path.suffix.lower().lstrip("."),
            "hint": hint or "both",
            "source": source_default,
            "n_tss": n_tss,
            "n_tts": n_tts,
            "n_rejected": file_rejected,
            "notes": list(file_notes),
        })

    tss, tss_duplicates = _dedupe(tss)
    tts, tts_duplicates = _dedupe(tts)
    duplicates = tss_duplicates + tts_duplicates
    shared = _shared_positions(tss, tts)

    if tss and tts:
        status = STATUS_AVAILABLE
        reason = ""
    elif tss or tts:
        status = STATUS_PARTIAL
        reason = ("solo hay TSS precalculados para este genoma; no hay TTS"
                  if tss else
                  "solo hay TTS precalculados para este genoma; no hay TSS")
    else:
        status = STATUS_NOT_AVAILABLE
        reason = (
            f"sin TSS/TTS precalculados para '{genome_name}'; no se predicen "
            "posiciones"
        )
        unique_notes = list(dict.fromkeys(note for note in notes if note))
        if unique_notes:
            reason += " (" + "; ".join(unique_notes) + ")"

    return {
        "genome": genome_name,
        "data_root": str(root),
        "status": status,
        "available": status != STATUS_NOT_AVAILABLE,
        "reason": reason,
        "tss": tss,
        "tts": tts,
        "n_tss": len(tss),
        "n_tts": len(tts),
        "n_duplicates": duplicates,
        "n_rejected": rejected,
        "n_shared_positions": shared,
        "sources": sources,
        "searched": searched,
        "limitation": LIMITATION,
    }


def _index_landmarks(records: Sequence[Mapping[str, Any]] | None
                     ) -> dict[tuple[str, str], tuple[list[int],
                                                      list[dict[str, Any]]]]:
    grouped: dict[tuple[str, str], list[tuple[int, dict[str, Any]]]] = {}
    for record in records or []:
        if not isinstance(record, Mapping):
            continue
        pos = normalize_pos(record.get("pos"))
        strand = normalize_strand(record.get("strand"))
        if pos is None or strand not in {"+", "-"}:
            continue
        key = (str(record.get("contig") or ""), strand)
        grouped.setdefault(key, []).append((pos, dict(record)))
    index: dict[tuple[str, str], tuple[list[int], list[dict[str, Any]]]] = {}
    for key, values in grouped.items():
        values.sort(key=lambda pair: pair[0])
        index[key] = ([pos for pos, _ in values],
                      [record for _, record in values])
    return index


def _upstream(index: Mapping[tuple[str, str], tuple[list[int],
                                                    list[dict[str, Any]]]],
              contig: str, strand: str, start: int, end: int,
              max_distance: int) -> tuple[dict[str, Any], int] | None:
    values = index.get((contig, strand))
    if not values:
        return None
    positions, records = values
    if strand == "+":
        i = bisect_left(positions, start) - 1
        if i < 0:
            return None
        distance = start - positions[i]
    else:
        i = bisect_right(positions, end)
        if i >= len(positions):
            return None
        distance = positions[i] - end
    if distance > max_distance:
        return None
    return records[i], int(distance)


def _downstream(index: Mapping[tuple[str, str], tuple[list[int],
                                                      list[dict[str, Any]]]],
                contig: str, strand: str, start: int, end: int,
                max_distance: int) -> tuple[dict[str, Any], int] | None:
    values = index.get((contig, strand))
    if not values:
        return None
    positions, records = values
    if strand == "+":
        i = bisect_right(positions, end)
        if i >= len(positions):
            return None
        distance = positions[i] - end
    else:
        i = bisect_left(positions, start) - 1
        if i < 0:
            return None
        distance = start - positions[i]
    if distance > max_distance:
        return None
    return records[i], int(distance)


def _inside(index: Mapping[tuple[str, str], tuple[list[int],
                                                  list[dict[str, Any]]]],
            contig: str, strand: str, start: int, end: int) -> int:
    values = index.get((contig, strand))
    if not values:
        return 0
    positions, _records = values
    return max(0, bisect_right(positions, end) - bisect_left(positions, start))


def attach(features: list[dict[str, Any]],
           tss: Sequence[Mapping[str, Any]] | None,
           tts: Sequence[Mapping[str, Any]] | None,
           *, max_distance: int = MAX_ATTACH_DISTANCE) -> list[dict[str, Any]]:
    """Enlaza el TSS aguas arriba y el TTS aguas abajo más cercanos.

    Solo considera marcas del mismo ``contig`` y hebra, a una distancia
    ``<= max_distance``. Para hebra ``+`` el TSS está antes de ``start`` y el
    TTS después de ``end``; para hebra ``-`` al revés. También deja
    ``n_tss_inside``/``n_tts_inside`` con las marcas dentro del feature.

    No crea posiciones: si no hay una marca cercana, deja ``None``. Muta las
    features (dicts) y devuelve la misma lista.
    """
    try:
        limit = max(0, int(max_distance))
    except (TypeError, ValueError):
        limit = MAX_ATTACH_DISTANCE
    tss_index = _index_landmarks(tss)
    tts_index = _index_landmarks(tts)
    for feature in features or []:
        if not isinstance(feature, dict):
            continue
        strand = normalize_strand(feature.get("strand"))
        start = normalize_pos(feature.get("start"))
        end = normalize_pos(feature.get("end"))
        if strand not in {"+", "-"} or start is None or end is None:
            continue
        start, end = sorted((start, end))
        contig = str(feature.get("contig") or "")
        upstream = _upstream(tss_index, contig, strand, start, end, limit)
        downstream = _downstream(tts_index, contig, strand, start, end, limit)
        feature["tss_upstream"] = upstream[0] if upstream else None
        feature["tss_upstream_distance"] = upstream[1] if upstream else None
        feature["tts_downstream"] = downstream[0] if downstream else None
        feature["tts_downstream_distance"] = (
            downstream[1] if downstream else None)
        feature["n_tss_inside"] = _inside(
            tss_index, contig, strand, start, end)
        feature["n_tts_inside"] = _inside(
            tts_index, contig, strand, start, end)
    return features


def summary(data_or_tss: Any, tts: Sequence[Mapping[str, Any]] | None = None
            ) -> dict[str, Any]:
    """Resume una carga de :func:`load` o pares ``(tss, tts)``.

    Acepta también directamente listas: ``summary(tss, tts)``. No modifica
    las filas.
    """
    if isinstance(data_or_tss, Mapping):
        payload = data_or_tss
        tss_rows = list(payload.get("tss") or [])
        tts_rows = list(payload.get("tts") or [])
        genome = str(payload.get("genome", ""))
        duplicates = int(payload.get("n_duplicates", 0) or 0)
        rejected = int(payload.get("n_rejected", 0) or 0)
        shared = int(payload.get("n_shared_positions", 0) or 0)
        sources = list(payload.get("sources") or [])
        reason = str(payload.get("reason", ""))
        status = str(payload.get("status", ""))
    else:
        tss_rows = list(data_or_tss or [])
        tts_rows = list(tts or [])
        genome = ""
        duplicates = rejected = shared = 0
        sources = []
        reason = ""
        status = ""

    n_tss = len(tss_rows)
    n_tts = len(tts_rows)
    if not status:
        if n_tss and n_tts:
            status = STATUS_AVAILABLE
        elif n_tss or n_tts:
            status = STATUS_PARTIAL
        else:
            status = STATUS_NOT_AVAILABLE

    combined = list(tss_rows) + list(tts_rows)
    positions = _positions(combined)
    strand_counts = {"+": 0, "-": 0, ".": 0}
    for record in combined:
        strand = normalize_strand(record.get("strand"))
        strand_counts[strand] += 1

    return {
        "genome": genome,
        "status": status,
        "available": status != STATUS_NOT_AVAILABLE,
        "n_tss": n_tss,
        "n_tts": n_tts,
        "n_total": n_tss + n_tts,
        "n_tss_unique_positions": len({
            (record.get("contig"), record.get("pos"), record.get("strand"))
            for record in tss_rows
        }),
        "n_tts_unique_positions": len({
            (record.get("contig"), record.get("pos"), record.get("strand"))
            for record in tts_rows
        }),
        "n_duplicates": duplicates,
        "n_rejected": rejected,
        "n_shared_positions": shared,
        "by_strand": strand_counts,
        "tss_classes": _counts(tss_rows, "class"),
        "tts_classes": _counts(tts_rows, "class"),
        "tss_locations": _counts(tss_rows, "location"),
        "tts_locations": _counts(tts_rows, "location"),
        "tss_sources": _counts(tss_rows, "source"),
        "tts_sources": _counts(tts_rows, "source"),
        "mechanisms": _counts(tts_rows, "mechanism"),
        "contigs": sorted({
            str(record.get("contig")) for record in combined
            if record.get("contig")
        }),
        "pos_min": min(positions) if positions else None,
        "pos_max": max(positions) if positions else None,
        "sources": sources,
        "reason": reason,
        "limitation": LIMITATION,
    }


__all__ = [
    "TSS",
    "TTS",
    "STATUS_AVAILABLE",
    "STATUS_PARTIAL",
    "STATUS_NOT_AVAILABLE",
    "MAX_ATTACH_DISTANCE",
    "LIMITATION",
    "load",
    "attach",
    "summary",
    "normalize_strand",
    "normalize_pos",
    "normalize_source",
    "normalize_class",
    "normalize_location",
    "normalize_mechanism",
]
