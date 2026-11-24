#!/usr/bin/env python3
"""Paridad de anotación base contra un GFF3 de referencia local.

Compara nuestra anotación (``data/annotations_ab_initio/<genoma>/annotation.gff3``
o un GFF de ``out_fl``) con una referencia local (NCBI/RefSeq, Bakta, Prokka,
etc.) y describe diferencias accionables: genes ausentes/extra, fragmentación,
fusiones, cambios de hebra, desplazamientos de borde, ORFs pequeños, ncRNA y
productos/hipotéticas.

No usa red, GPU ni binarios externos. Si no hay referencia, el estado es
``not_evaluable``: no se inventan precision/recall ni diferencias.

Uso típico::

    .venv/bin/python tools/annotation_parity.py --genome GCF_000006685.1_ASM668v1

El GFF de referencia se autodescubre en ``data/annotations/<genoma>.gff.gz``.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shlex
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOME = "GCF_000006685.1_ASM668v1"
DEFAULT_MD = ROOT / "reports" / "annotation_parity.md"
DEFAULT_JSON = ROOT / "reports" / "annotation_parity.json"

SCHEMA_VERSION = 1
CLASSES = ("CDS", "tRNA", "rRNA", "ncRNA")
SMALL_ORF_MAX_AA = 150

# Tipos GFF3 que cuentan como ncRNA. Se acepta además cualquier tipo que
# termine en ``_RNA``/``RNA`` y no sea tRNA/rRNA (p. ej. tmRNA, ncRNA).
NC_RNA_TYPES = {
    "ncrna",
    "ncrna_gene",
    "misc_rna",
    "rnase_p_rna",
    "rnase_mrp_rna",
    "srp_rna",
    "tmrna",
    "snrna",
    "snorna",
    "lncrna",
    "mirna",
    "scrna",
    "precursor_rna",
    "antisense_rna",
    "vault_rna",
    "y_rna",
    "guide_rna",
}

HYPOTHETICAL_RE = re.compile(
    r"\b(?:hypothetical|uncharacterized|unknown\s+protein|"
    r"protein\s+of\s+unknown\s+function|duf\d+|upf\d+)\b",
    re.IGNORECASE,
)

SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2}

STATUS_LABELS = {
    "evaluated": "evaluado",
    "not_evaluable": "no evaluable",
    "not_present": "no presente",
    "reference_absent": "ausente en la referencia",
    "ours_absent": "ausente en nuestra anotación",
    "not_evaluable_contigs": "no evaluable (sin contigs compartidos)",
}


class ParityError(RuntimeError):
    """Error de uso: falta nuestra anotación o un insumo es ilegible."""


# ---------------------------------------------------------------------------
# Utilidades generales
# ---------------------------------------------------------------------------
def _open_text(path: Path):
    """Abre un GFF/TSV normal o gzip sin depender de la extensión."""
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _rel(root: Path, path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return Path(path).as_posix()


def _pct(numerator: int | float, denominator: int | float) -> float | None:
    if not denominator:
        return None
    return round(100.0 * float(numerator) / float(denominator), 2)


def _ratio(numerator: int | float, denominator: int | float) -> float | None:
    if not denominator:
        return None
    return round(float(numerator) / float(denominator), 4)


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None or (precision + recall) == 0:
        return None
    return round(2 * precision * recall / (precision + recall), 4)


def _mean(values: Iterable[float]) -> float | None:
    vals = list(values)
    if not vals:
        return None
    return round(sum(vals) / len(vals), 4)


def _median(values: Iterable[float]) -> float | None:
    vals = list(values)
    if not vals:
        return None
    return round(statistics.median(vals), 4)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def _fmt_int(value: Any) -> str:
    if value is None:
        return "—"
    return f"{int(value):,}".replace(",", ".")


def _fmt_pct(value: Any) -> str:
    if value is None:
        return "—"
    return f"{float(value):.1f} %".replace(".", ",")


def _fmt_ratio(value: Any) -> str:
    if value is None:
        return "—"
    return f"{float(value):.3f}".replace(".", ",")


def _fmt_nt(value: Any) -> str:
    if value is None:
        return "—"
    return f"{int(value)} nt"


def _escape_md(text: Any) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def _truncate(text: str, limit: int = 180) -> str:
    text = (text or "").replace("\n", " ").strip()
    if len(text) <= limit:
        return text
    return text[: limit - 1].rstrip() + "…"


# ---------------------------------------------------------------------------
# Atributos y features GFF3
# ---------------------------------------------------------------------------
def parse_attributes(text: str) -> dict[str, str]:
    """Parsea ``clave=valor;...`` con unquote; conserva el último valor."""
    out: dict[str, str] = {}
    for part in (text or "").split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = unquote(key.strip())
        value = unquote(value.strip())
        if key:
            out[key] = value
    return out


def _is_pseudo(attrs: Mapping[str, str]) -> bool:
    for key in ("pseudo", "pseudogene", "is_pseudogene"):
        value = (attrs.get(key) or "").strip().lower()
        if value in {"true", "1", "yes"}:
            return True
    return False


def feature_class(raw_type: str) -> str | None:
    """Mapea un tipo GFF3 a una de las cuatro clases comparables."""
    raw = (raw_type or "").strip()
    low = raw.lower()
    if low == "cds":
        return "CDS"
    if low == "trna":
        return "tRNA"
    if low == "rrna":
        return "rRNA"
    if low in NC_RNA_TYPES:
        return "ncRNA"
    # Tipos bacterianos habituales: tmRNA, RNase_P_RNA, SRP_RNA, ncRNA...
    if low.endswith("rna") and "trna" not in low and "rrna" not in low:
        return "ncRNA"
    return None


def _approx_aa(feature: Mapping[str, Any]) -> int:
    """Estimación de aminoácidos a partir del ADN; no traduce la secuencia.

    Se resta un codón de parada salvo que el atributo ``partial`` indique que
    falta el extremo 3'. Para CDS multiparte se suma el largo de los segmentos.
    """
    nt = int(feature.get("span_nt") or 0)
    partial = str((feature.get("attrs") or {}).get("partial", "") or "")
    no_stop = (len(partial) >= 2 and partial[1] == "1") or partial.lower() == "true"
    aa = nt // 3 - (0 if no_stop else 1)
    return max(0, aa)


def _is_small_orf(feature: Mapping[str, Any], max_aa: int = SMALL_ORF_MAX_AA) -> bool:
    return feature["class"] == "CDS" and _approx_aa(feature) < max_aa


def _first_product(attrs: Mapping[str, str]) -> str:
    product = (attrs.get("product") or "").strip()
    if product:
        return product
    return ""


def _make_feature(
    *,
    seqid: str,
    raw_type: str,
    cls: str,
    start: int,
    end: int,
    strand: str,
    source: str,
    attrs: dict[str, str],
    line: int,
    segments: Sequence[tuple[int, int]] | None = None,
) -> dict[str, Any]:
    segs = sorted(segments or [(start, end)])
    span_nt = sum(e - s + 1 for s, e in segs)
    return {
        "class": cls,
        "type": raw_type,
        "seqid": seqid,
        "start": start,
        "end": end,
        "strand": strand,
        "source": source,
        "attrs": attrs,
        "line": line,
        "id": (attrs.get("ID") or "").strip(),
        "name": (attrs.get("Name") or "").strip(),
        "gene": (attrs.get("gene") or "").strip(),
        "locus_tag": (attrs.get("locus_tag") or "").strip(),
        "product": _first_product(attrs),
        "is_pseudo": _is_pseudo(attrs),
        "segments": segs,
        "span_nt": span_nt,
        "length_nt": end - start + 1,
        "multipart": len(segs) > 1,
        "approx_aa": 0,  # se rellena al final del parseo
    }


def _merge_cds_group(group: list[dict[str, Any]]) -> dict[str, Any]:
    """Fusiona líneas CDS con el mismo ID (CDS multiparte) en una feature."""
    group = sorted(group, key=lambda f: (f["start"], f["end"], f["line"]))
    attrs = dict(group[0]["attrs"])
    for feature in group[1:]:
        for key, value in feature["attrs"].items():
            attrs.setdefault(key, value)
    segments = sorted({(f["start"], f["end"]) for f in group})
    start = min(f["start"] for f in group)
    end = max(f["end"] for f in group)
    merged = _make_feature(
        seqid=group[0]["seqid"],
        raw_type="CDS",
        cls="CDS",
        start=start,
        end=end,
        strand=group[0]["strand"],
        source=group[0]["source"],
        attrs=attrs,
        line=group[0]["line"],
        segments=segments,
    )
    merged["is_pseudo"] = any(f["is_pseudo"] for f in group)
    merged["multipart"] = True
    return merged


def parse_gff(path: Path, exclude_pseudo: bool = True) -> dict[str, Any]:
    """Parsea un GFF3 y agrupa CDS multiparte por ``ID``.

    Las coordenadas se conservan 1-based inclusive. ``span_nt`` es la suma de
    segmentos (útil para aminoácidos); ``length_nt`` es el sobre mínimo..máximo
    (útil para solapamiento).
    """
    path = Path(path)
    features: list[dict[str, Any]] = []
    raw_type_counts: Counter[str] = Counter()
    sequence_regions: dict[str, int] = {}
    excluded_pseudo = 0
    invalid_lines = 0
    cds_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    cds_no_id: list[dict[str, Any]] = []
    rna_features: list[dict[str, Any]] = []
    pseudo_features: list[dict[str, Any]] = []

    with _open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            sequence_regions[parts[1]] = int(parts[3])
                        except ValueError:
                            pass
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                invalid_lines += 1
                continue
            raw_type = cols[2].strip()
            raw_type_counts[raw_type] += 1
            cls = feature_class(raw_type)
            if cls is None:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                invalid_lines += 1
                continue
            if start > end:
                start, end = end, start
            attrs = parse_attributes(cols[8])
            strand = cols[6] if cols[6] in {"+", "-"} else "."
            feature = _make_feature(
                seqid=cols[0],
                raw_type=raw_type,
                cls=cls,
                start=start,
                end=end,
                strand=strand,
                source=cols[1],
                attrs=attrs,
                line=line_no,
            )
            if cls == "CDS":
                key = (feature["seqid"], feature["id"])
                if feature["id"]:
                    cds_groups[key].append(feature)
                else:
                    cds_no_id.append(feature)
            else:
                rna_features.append(feature)

    # Grupos con hebras distintas no se fusionan: se conservan como líneas.
    cds_features: list[dict[str, Any]] = []
    for group in cds_groups.values():
        strands = {f["strand"] for f in group}
        if len(group) == 1 or len(strands) > 1:
            cds_features.extend(group)
        else:
            cds_features.append(_merge_cds_group(group))
    cds_features.extend(cds_no_id)

    # El filtro de pseudogenes se aplica a la feature ya fusionada, no a una
    # línea suelta de un CDS multiparte.
    for feature in cds_features:
        if exclude_pseudo and feature["is_pseudo"]:
            excluded_pseudo += 1
            pseudo_features.append(feature)
            continue
        features.append(feature)
    features.extend(rna_features)

    for feature in features + pseudo_features:
        feature["approx_aa"] = _approx_aa(feature)

    features.sort(key=lambda f: (f["seqid"], f["start"], f["end"], f["class"]))
    counts = Counter(f["class"] for f in features)
    contigs = sorted({f["seqid"] for f in features} | set(sequence_regions))
    return {
        "path": str(path),
        "features": features,
        "feature_counts": {cls: int(counts.get(cls, 0)) for cls in CLASSES},
        "raw_type_counts": dict(sorted(raw_type_counts.items())),
        "contigs": contigs,
        "sequence_regions": sequence_regions,
        "excluded_pseudo": excluded_pseudo,
        "pseudo_features": pseudo_features,
        "invalid_lines": invalid_lines,
        "multipart_cds": sum(1 for f in features if f["class"] == "CDS" and f["multipart"]),
        "total_features": len(features),
    }


# ---------------------------------------------------------------------------
# Comparación de features
# ---------------------------------------------------------------------------
def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Longitud de intersección 1-based inclusive; 0 si no se tocan."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def reciprocal_overlap(a: Mapping[str, Any], b: Mapping[str, Any]) -> float:
    """Mínimo de las dos fracciones de solapamiento (0.0 a 1.0)."""
    ov = interval_overlap(a["start"], a["end"], b["start"], b["end"])
    if ov <= 0:
        return 0.0
    return min(ov / max(1, a["length_nt"]), ov / max(1, b["length_nt"]))


def _segments_match(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    return list(a.get("segments") or []) == list(b.get("segments") or [])


def match_features(
    reference_features: Sequence[Mapping[str, Any]],
    query_features: Sequence[Mapping[str, Any]],
    threshold: float = 0.5,
) -> list[dict[str, Any]]:
    """Emparejamiento 1:1 voraz por solapamiento recíproco.

    El alcance es ``(clase, contig, hebra)``. Es determinista, pero no es un
    emparejamiento bipartito óptimo; se documenta para que sea reproducible.
    """
    ref_by_scope: dict[tuple[str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    query_by_scope: dict[tuple[str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for feature in reference_features:
        ref_by_scope[(feature["class"], feature["seqid"], feature["strand"])].append(feature)
    for feature in query_features:
        query_by_scope[(feature["class"], feature["seqid"], feature["strand"])].append(feature)

    matches: list[dict[str, Any]] = []
    for scope in sorted(set(ref_by_scope) | set(query_by_scope)):
        refs = ref_by_scope.get(scope, [])
        queries = query_by_scope.get(scope, [])
        if not refs or not queries:
            continue
        pairs: list[tuple[float, int, int]] = []
        for i, ref in enumerate(refs):
            for j, query in enumerate(queries):
                score = reciprocal_overlap(ref, query)
                if score >= threshold and score > 0.0:
                    pairs.append((score, i, j))
        pairs.sort(
            key=lambda p: (
                -p[0],
                refs[p[1]]["start"],
                refs[p[1]]["end"],
                queries[p[2]]["start"],
                queries[p[2]]["end"],
            )
        )
        used_ref: set[int] = set()
        used_query: set[int] = set()
        for score, i, j in pairs:
            if i in used_ref or j in used_query:
                continue
            used_ref.add(i)
            used_query.add(j)
            ref, query = refs[i], queries[j]
            matches.append(
                {
                    "ref": ref,
                    "query": query,
                    "reciprocal_overlap": round(score, 4),
                    "exact": bool(
                        ref["start"] == query["start"]
                        and ref["end"] == query["end"]
                        and ref["strand"] == query["strand"]
                    ),
                    "segments_match": _segments_match(ref, query),
                }
            )
    return matches


def _boundary_deltas(ref: Mapping[str, Any], query: Mapping[str, Any]) -> dict[str, int]:
    """Deltas genómicos y 5'/3' (valores absolutos en nt)."""
    start_delta = int(query["start"]) - int(ref["start"])
    end_delta = int(query["end"]) - int(ref["end"])
    if ref["strand"] == "-":
        five = abs(end_delta)
        three = abs(start_delta)
    else:
        five = abs(start_delta)
        three = abs(end_delta)
    return {
        "start_delta": start_delta,
        "end_delta": end_delta,
        "five_prime_delta_nt": five,
        "three_prime_delta_nt": three,
        "max_boundary_delta_nt": max(five, three),
    }


def _is_multipart_structure_mismatch(m: Mapping[str, Any]) -> bool:
    """Solo compara estructura si al menos un lado es multiparte y el sobre coincide."""
    ref, query = m["ref"], m["query"]
    if not (ref.get("multipart") or query.get("multipart")):
        return False
    return bool(
        ref["start"] == query["start"]
        and ref["end"] == query["end"]
        and not m["segments_match"]
    )


def _length_ratio(ref: Mapping[str, Any], query: Mapping[str, Any]) -> float | None:
    if not ref["span_nt"]:
        return None
    return round(float(query["span_nt"]) / float(ref["span_nt"]), 4)


def _is_length_discordant(
    ref: Mapping[str, Any], query: Mapping[str, Any], min_ratio: float = 0.8,
    max_ratio: float = 1.25, min_diff_nt: int = 50,
) -> bool:
    ratio = _length_ratio(ref, query)
    if ratio is None:
        return False
    diff = abs(int(query["span_nt"]) - int(ref["span_nt"]))
    return diff >= min_diff_nt and (ratio < min_ratio or ratio > max_ratio)


def _partition_matches(
    reference_features: Sequence[dict],
    query_features: Sequence[dict],
    matches: Sequence[Mapping[str, Any]],
) -> tuple[list[dict], list[dict]]:
    matched_ref = {id(m["ref"]) for m in matches}
    matched_query = {id(m["query"]) for m in matches}
    unmatched_ref = [f for f in reference_features if id(f) not in matched_ref]
    unmatched_query = [f for f in query_features if id(f) not in matched_query]
    return unmatched_ref, unmatched_query


# ---------------------------------------------------------------------------
# Productos e hipotéticas
# ---------------------------------------------------------------------------
def normalize_product(text: str) -> str:
    text = (text or "").lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def product_tokens(text: str) -> set[str]:
    return {tok for tok in normalize_product(text).split() if tok}


def product_token_jaccard(a: str, b: str) -> float | None:
    ta, tb = product_tokens(a), product_tokens(b)
    if not ta and not tb:
        return None
    union = ta | tb
    if not union:
        return None
    return len(ta & tb) / len(union)


def is_hypothetical_product(product: str) -> bool:
    """Clasificador conservador; no usa 'putative' a secas."""
    text = (product or "").strip()
    if not text:
        return False
    return bool(HYPOTHETICAL_RE.search(text))


def _clean_product_text(text: str, from_functional_tsv: bool = False) -> str:
    """Limpia descripciones Swiss-Prot/TrEMBL sin cambiar el significado."""
    value = (text or "").strip()
    if not value:
        return ""
    if from_functional_tsv and (value.startswith("sp|") or value.startswith("tr|")):
        # "sp|Q9PLU4|HEM2_CHLMU Delta-aminolevulinic acid dehydratase OS=..."
        parts = value.split(" ", 1)
        value = parts[1] if len(parts) > 1 else ""
    # Quita metadatos Swiss-Prot que no son parte del nombre del producto.
    for marker in (" OS=", " OX=", " GN=", " PE=", " SV=", " EC="):
        if marker in value:
            value = value.split(marker, 1)[0]
    return value.strip()


def _product_stats(features: Sequence[Mapping[str, Any]], source: str) -> dict[str, Any]:
    cds = [f for f in features if f["class"] == "CDS"]
    with_product = [f for f in cds if (f.get("product") or "").strip()]
    hypothetical = [f for f in with_product if is_hypothetical_product(f["product"])]
    status = "evaluated" if with_product else "not_available"
    return {
        "source": source,
        "status": status,
        "cds": len(cds),
        "with_product": len(with_product),
        "without_product": len(cds) - len(with_product),
        "product_coverage_pct": _pct(len(with_product), len(cds)),
        "hypothetical": len(hypothetical),
        "pct_hypothetical": _pct(len(hypothetical), len(with_product)),
    }


def load_products_tsv(path: Path) -> dict[str, str]:
    """Carga un TSV ``id<TAB>...<TAB>producto`` (functional.tsv o 2 columnas)."""
    mapping: dict[str, str] = {}
    with _open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2:
                continue
            key = (row[0] or "").strip()
            if not key or key.startswith("#"):
                continue
            raw = (row[-1] or "").strip()
            if not raw:
                continue
            product = _clean_product_text(raw, from_functional_tsv=True)
            if product:
                mapping[key] = product
    return mapping


def _candidate_product_files(ours_path: Path) -> list[Path]:
    parent = ours_path.parent
    stem = ours_path.name
    for suffix in (".gff3.gz", ".gff.gz", ".gff3", ".gff"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    candidates = [
        parent / "functional.tsv",
        parent / "products.tsv",
        parent / f"{stem}.functional.tsv",
        parent / f"{stem}.products.tsv",
    ]
    if parent.name == "annotation":
        candidates.extend(
            [
                parent.parent / f"{stem}.functional.tsv",
                parent.parent / "functional.tsv",
            ]
        )
    return candidates


def resolve_products(
    ours_path: Path,
    ours_parsed: Mapping[str, Any],
    requested: str,
) -> tuple[str, str, Path | None]:
    """Devuelve ``(fuente, etiqueta, path)`` para los productos de nuestra GFF.

    Prioridad: atributo ``product`` del GFF; si no hay, ``--ours-products`` o
    ``functional.tsv`` vecino. No se mezclan fuentes.
    """
    cds = [f for f in ours_parsed["features"] if f["class"] == "CDS"]
    gff_products = sum(1 for f in cds if (f.get("product") or "").strip())
    if gff_products:
        return "gff", "atributo product del GFF", None
    if requested == "none":
        return "none", "sin productos", None
    if requested != "auto":
        path = Path(requested)
        if not path.is_file():
            raise ParityError(f"no existe --ours-products: {path}")
        return "tsv", f"TSV explícito: {path.name}", path
    for path in _candidate_product_files(ours_path):
        if path.is_file():
            return "tsv", f"TSV vecino: {path.name}", path
    return "none", "sin productos", None


def apply_products(
    ours_parsed: dict[str, Any], products: Mapping[str, str], source: str
) -> None:
    if source not in {"tsv", "gff"}:
        return
    for feature in ours_parsed["features"]:
        if feature["class"] != "CDS":
            continue
        if source == "gff" and feature.get("product"):
            continue
        for key in (feature.get("id"), feature.get("name"), feature.get("locus_tag")):
            if key and key in products:
                feature["product"] = products[key]
                break


# ---------------------------------------------------------------------------
# Descubrimiento de entradas
# ---------------------------------------------------------------------------
def _infer_genome(path: Path) -> str:
    name = path.name
    for suffix in (".gff3.gz", ".gff.gz", ".gff3", ".gff"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name == "annotation":
        return path.parent.name
    if name.endswith("_genomic"):
        name = name[: -len("_genomic")]
    return name or path.parent.name


def _existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path
    return None


def discover_ours(root: Path, genome: str) -> Path | None:
    base = root / "data" / "annotations_ab_initio" / genome
    candidates = [
        base / "annotation.gff3",
        base / "annotation.gff",
        base / f"{genome}.gff3",
        root / "out_fl" / "annotation" / f"{genome}.gff3",
        root / "out_fl" / f"{genome}.gff3",
        root / "out_fl" / genome / "annotation.gff3",
        root / "out_fl" / genome / f"{genome}.gff3",
        root / "out_fl" / f"{genome}.annotation.gff3",
    ]
    found = _existing(candidates)
    if found:
        return found
    if base.is_dir():
        for pattern in ("*.gff3", "*.gff"):
            hits = sorted(base.glob(pattern))
            if hits:
                return hits[0]
    return None


def discover_reference(root: Path, genome: str) -> Path | None:
    candidates: list[Path] = []
    for suffix in (".gff.gz", ".gff3.gz", ".gff", ".gff3"):
        candidates.extend(
            [
                root / "data" / "annotations" / f"{genome}{suffix}",
                root / "data" / "annotations" / f"{genome}_genomic{suffix}",
                root / "data" / "references" / f"{genome}{suffix}",
                root / "data" / "references" / f"{genome}_genomic{suffix}",
                root / "data" / "annotations" / genome / f"annotation{suffix}",
            ]
        )
    return _existing(candidates)


def _resolve_file(value: str | Path | None, *, directory_patterns: Sequence[str]) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_dir():
        for pattern in directory_patterns:
            hits = sorted(path.glob(pattern))
            if hits:
                return hits[0]
        return None
    return path


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
def _feature_ref(feature: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if feature is None:
        return None
    return {
        "id": feature.get("id") or None,
        "type": feature.get("type"),
        "class": feature.get("class"),
        "seqid": feature.get("seqid"),
        "start": int(feature.get("start", 0)),
        "end": int(feature.get("end", 0)),
        "strand": feature.get("strand"),
        "length_nt": int(feature.get("span_nt", 0)),
        "envelope_nt": int(feature.get("length_nt", 0)),
        "approx_aa": int(feature.get("approx_aa", 0)),
        "multipart": bool(feature.get("multipart")),
        "segments": [list(s) for s in feature.get("segments", [])],
        "product": feature.get("product") or "",
        "gene": feature.get("gene") or "",
        "locus_tag": feature.get("locus_tag") or "",
        "source": feature.get("source") or "",
    }


def _feature_refs(features: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [f for f in (_feature_ref(x) for x in features) if f is not None]


def _feature_label(feature: Mapping[str, Any] | None) -> str:
    if not feature:
        return "—"
    return f"{feature['seqid']}:{feature['start']}-{feature['end']}({feature['strand']})"


def _features_label(features: Sequence[Mapping[str, Any]] | None) -> str:
    if not features:
        return "—"
    labels = [_feature_label(f) for f in features]
    if len(labels) <= 3:
        return "; ".join(labels)
    return "; ".join(labels[:3]) + f"; +{len(labels) - 3}"


def _class_metrics(
    cls: str,
    ref_all: Sequence[Mapping[str, Any]],
    ours_all: Sequence[Mapping[str, Any]],
    ref_scope: Sequence[Mapping[str, Any]],
    ours_scope: Sequence[Mapping[str, Any]],
    matches: Sequence[Mapping[str, Any]],
    unmatched_ref: Sequence[Mapping[str, Any]],
    unmatched_ours: Sequence[Mapping[str, Any]],
    *,
    shared_contigs: bool,
    threshold: float,
    boundary_alert: int,
    coord_tolerance: int,
) -> dict[str, Any]:
    ref_scope_n = len(ref_scope)
    ours_scope_n = len(ours_scope)
    matched = len(matches)
    status = "evaluated"
    if not shared_contigs:
        status = "not_evaluable_contigs"
    elif ref_scope_n == 0 and ours_scope_n == 0:
        status = "not_present"
    elif ref_scope_n == 0:
        status = "reference_absent" if not ref_all else "not_evaluable_contigs"
    elif ours_scope_n == 0:
        status = "ours_absent" if not ours_all else "not_evaluable_contigs"

    evaluate = status == "evaluated"
    exact = sum(1 for m in matches if m["exact"])
    within = sum(
        1
        for m in matches
        if abs(int(m["query"]["start"]) - int(m["ref"]["start"])) <= coord_tolerance
        and abs(int(m["query"]["end"]) - int(m["ref"]["end"])) <= coord_tolerance
        and m["query"]["strand"] == m["ref"]["strand"]
    )
    overlaps = [float(m["reciprocal_overlap"]) for m in matches]
    boundary_alert_n = 0
    length_discordant = 0
    multipart_mismatch = 0
    for m in matches:
        deltas = _boundary_deltas(m["ref"], m["query"])
        if deltas["max_boundary_delta_nt"] >= boundary_alert:
            boundary_alert_n += 1
        if _is_length_discordant(m["ref"], m["query"]):
            length_discordant += 1
        if _is_multipart_structure_mismatch(m):
            multipart_mismatch += 1

    recall = _ratio(matched, ref_scope_n) if evaluate and ref_scope_n else None
    precision = _ratio(matched, ours_scope_n) if evaluate and ours_scope_n else None
    return {
        "class": cls,
        "status": status,
        "reference_total": len(ref_all),
        "ours_total": len(ours_all),
        "delta_ours_minus_reference": len(ours_all) - len(ref_all),
        "reference_in_scope": ref_scope_n,
        "ours_in_scope": ours_scope_n,
        "matched": matched if evaluate else None,
        "recall": recall,
        "precision": precision,
        "f1": _f1(precision, recall) if evaluate else None,
        "exact": exact if evaluate else None,
        "exact_rate": _ratio(exact, matched) if evaluate and matched else None,
        "within_tolerance": within if evaluate else None,
        "coordinate_tolerance_nt": coord_tolerance,
        "within_tolerance_rate": (
            _ratio(within, matched) if evaluate and matched else None
        ),
        "mean_reciprocal_overlap": _mean(overlaps) if evaluate else None,
        "median_reciprocal_overlap": _median(overlaps) if evaluate else None,
        "boundary_shift_alert": boundary_alert_n if evaluate else None,
        "boundary_alert_nt": boundary_alert,
        "length_discordant": length_discordant if evaluate else None,
        "multipart_structure_mismatch": multipart_mismatch if evaluate else None,
        "unmatched_reference": len(unmatched_ref) if evaluate else None,
        "unmatched_ours": len(unmatched_ours) if evaluate else None,
        "reciprocal_overlap_threshold": threshold,
    }


def _small_orf_stats(
    ref: Sequence[Mapping[str, Any]],
    ours: Sequence[Mapping[str, Any]],
    matches: Sequence[Mapping[str, Any]],
    actionable_items: Sequence[Mapping[str, Any]],
    max_aa: int,
) -> dict[str, Any]:
    ref_cds = [f for f in ref if f["class"] == "CDS"]
    ours_cds = [f for f in ours if f["class"] == "CDS"]
    ref_small = [f for f in ref_cds if _approx_aa(f) < max_aa]
    ours_small = [f for f in ours_cds if _approx_aa(f) < max_aa]
    matched_small = sum(
        1
        for m in matches
        if m["ref"]["class"] == "CDS"
        and m["query"]["class"] == "CDS"
        and _approx_aa(m["ref"]) < max_aa
        and _approx_aa(m["query"]) < max_aa
    )
    missing_small = sum(
        1
        for item in actionable_items
        if str(item.get("category", "")).startswith("missing_") and item.get("small_orf")
    )
    extra_small = sum(
        1
        for item in actionable_items
        if str(item.get("category", "")).startswith("extra_") and item.get("small_orf")
    )
    return {
        "threshold_aa": max_aa,
        "ours": {
            "cds": len(ours_cds),
            "small": len(ours_small),
            "pct": _pct(len(ours_small), len(ours_cds)),
        },
        "reference": {
            "cds": len(ref_cds),
            "small": len(ref_small),
            "pct": _pct(len(ref_small), len(ref_cds)),
        },
        "matched_small": matched_small,
        "missing_small": missing_small,
        "extra_small": extra_small,
    }


def _ncrna_stats(
    ref: Sequence[Mapping[str, Any]],
    ours: Sequence[Mapping[str, Any]],
    matches: Sequence[Mapping[str, Any]],
    unmatched_ref: Sequence[Mapping[str, Any]],
    unmatched_ours: Sequence[Mapping[str, Any]],
    *,
    evaluated: bool,
) -> dict[str, Any]:
    ref_nc = [f for f in ref if f["class"] == "ncRNA"]
    ours_nc = [f for f in ours if f["class"] == "ncRNA"]
    if not ref_nc and not ours_nc:
        status = "not_present"
    elif not ref_nc:
        status = "reference_absent"
    elif not ours_nc:
        status = "ours_absent"
    elif not evaluated:
        status = "not_evaluable_contigs"
    else:
        status = "evaluated"

    if status == "evaluated":
        matched = len(matches)
        missing = len(unmatched_ref)
        extra = len(unmatched_ours)
    elif status == "ours_absent":
        matched = 0
        missing = len(ref_nc)
        extra = 0
    else:
        matched = None
        missing = None
        extra = None
    type_mismatch = [
        {
            "reference_type": m["ref"]["type"],
            "ours_type": m["query"]["type"],
            "reference": _feature_ref(m["ref"]),
            "ours": _feature_ref(m["query"]),
        }
        for m in matches
        if m["ref"]["type"] != m["query"]["type"]
    ]
    return {
        "status": status,
        "reference_total": len(ref_nc),
        "ours_total": len(ours_nc),
        "reference_types": dict(sorted(Counter(f["type"] for f in ref_nc).items())),
        "ours_types": dict(sorted(Counter(f["type"] for f in ours_nc).items())),
        "matched": matched,
        "missing": missing,
        "extra": extra,
        "type_mismatches": type_mismatch,
    }


def _is_product_mismatch(reference_product: str, ours_product: str) -> bool:
    """Dos productos anotados con poco solapamiento de tokens (no hipotéticas)."""
    if not reference_product or not ours_product:
        return False
    if is_hypothetical_product(reference_product) or is_hypothetical_product(ours_product):
        return False
    score = product_token_jaccard(reference_product, ours_product)
    if score is None:
        return False
    return score < 0.34


def _paired_products(
    matches: Sequence[Mapping[str, Any]],
    ours_source: str,
) -> dict[str, Any]:
    cds_matches = [m for m in matches if m["ref"]["class"] == "CDS"]
    both = [
        m
        for m in cds_matches
        if (m["ref"].get("product") or "").strip()
        and (m["query"].get("product") or "").strip()
    ]
    literal = sum(
        1
        for m in both
        if normalize_product(m["ref"]["product"]) == normalize_product(m["query"]["product"])
    )
    jaccards = [
        product_token_jaccard(m["ref"]["product"], m["query"]["product"]) for m in both
    ]
    jaccards = [j for j in jaccards if j is not None]
    hypothetical_disagreements = [
        m
        for m in both
        if is_hypothetical_product(m["ref"]["product"])
        != is_hypothetical_product(m["query"]["product"])
    ]
    product_mismatches = [
        m
        for m in both
        if _is_product_mismatch(m["ref"]["product"], m["query"]["product"])
    ]
    return {
        "status": "evaluated" if both else "not_evaluable",
        "matched_cds": len(cds_matches),
        "both_with_product": len(both),
        "literal_matches": literal,
        "literal_agreement": _ratio(literal, len(both)),
        "token_jaccard_mean": _mean(jaccards),
        "token_jaccard_median": _median(jaccards),
        "hypothetical_disagreements": len(hypothetical_disagreements),
        "product_mismatches": len(product_mismatches),
        "our_product_source": ours_source,
        "examples_hypothetical_disagreements": [
            {
                "reference": _feature_ref(m["ref"]),
                "ours": _feature_ref(m["query"]),
                "reference_product": m["ref"]["product"],
                "ours_product": m["query"]["product"],
            }
            for m in hypothetical_disagreements[:10]
        ],
        "examples_product_mismatches": [
            {
                "reference": _feature_ref(m["ref"]),
                "ours": _feature_ref(m["query"]),
                "reference_product": m["ref"]["product"],
                "ours_product": m["query"]["product"],
                "token_jaccard": round(
                    product_token_jaccard(m["ref"]["product"], m["query"]["product"]) or 0.0,
                    4,
                ),
            }
            for m in product_mismatches[:10]
        ],
    }


# ---------------------------------------------------------------------------
# Construcción de diferencias accionables
# ---------------------------------------------------------------------------
def _item(
    category: str,
    severity: str,
    cls: str,
    seqid: str | None,
    message: str,
    *,
    reference: Mapping[str, Any] | None = None,
    ours: Mapping[str, Any] | None = None,
    reference_features: Sequence[Mapping[str, Any]] | None = None,
    ours_features: Sequence[Mapping[str, Any]] | None = None,
    metrics: Mapping[str, Any] | None = None,
    tags: Sequence[str] | None = None,
    small_orf: bool = False,
) -> dict[str, Any]:
    return {
        "category": category,
        "severity": severity,
        "feature_class": cls,
        "seqid": seqid,
        "reference": dict(reference) if reference else None,
        "ours": dict(ours) if ours else None,
        "reference_features": [dict(f) for f in (reference_features or [])] or None,
        "ours_features": [dict(f) for f in (ours_features or [])] or None,
        "metrics": dict(metrics or {}),
        "tags": list(tags or []),
        "small_orf": bool(small_orf),
        "message": message,
    }


def _sort_items(items: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    def key(item: Mapping[str, Any]) -> tuple:
        ref = item.get("reference") or {}
        ours = item.get("ours") or {}
        length = max(
            int(ref.get("length_nt") or 0),
            int(ours.get("length_nt") or 0),
            sum(int(f.get("length_nt") or 0) for f in (item.get("reference_features") or [])),
            sum(int(f.get("length_nt") or 0) for f in (item.get("ours_features") or [])),
        )
        return (
            SEVERITY_ORDER.get(item["severity"], 9),
            -length,
            str(item.get("seqid") or ""),
            int(ref.get("start") or ours.get("start") or 0),
            item["category"],
        )

    return sorted((dict(item) for item in items), key=key)


def _best_overlap(
    feature: Mapping[str, Any],
    candidates: Sequence[Mapping[str, Any]],
    threshold: float,
) -> Mapping[str, Any] | None:
    hits = [
        c
        for c in candidates
        if c["seqid"] == feature["seqid"]
        and c["strand"] == feature["strand"]
        and reciprocal_overlap(feature, c) >= threshold
    ]
    if not hits:
        return None
    return max(hits, key=lambda c: reciprocal_overlap(feature, c))


def _categorize_unmatched(
    *,
    cls: str,
    unmatched_ref: Sequence[Mapping[str, Any]],
    unmatched_ours: Sequence[Mapping[str, Any]],
    pseudo_ref: Sequence[Mapping[str, Any]],
    pseudo_ours: Sequence[Mapping[str, Any]],
    actionable: list[dict[str, Any]],
    review: list[dict[str, Any]],
    threshold: float,
) -> None:
    """Convierte no emparejados en categorías explicativas, sin duplicar."""
    consumed_ours: set[int] = set()
    consumed_ref: set[int] = set()

    # 1) Un gen de referencia partido en varias features nuestras.
    for ref in unmatched_ref:
        pieces = [
            q
            for q in unmatched_ours
            if id(q) not in consumed_ours
            and q["seqid"] == ref["seqid"]
            and q["strand"] == ref["strand"]
            and interval_overlap(
                ref["start"], ref["end"], q["start"], q["end"]
            ) > 0
        ]
        if len(pieces) >= 2:
            covered = sum(
                interval_overlap(ref["start"], ref["end"], q["start"], q["end"])
                for q in pieces
            )
            if covered >= 0.5 * ref["length_nt"]:
                for q in pieces:
                    consumed_ours.add(id(q))
                consumed_ref.add(id(ref))
                actionable.append(
                    _item(
                        "fragmented_in_ours",
                        "high",
                        cls,
                        ref["seqid"],
                        (
                            f"La referencia tiene 1 {cls} de {ref['span_nt']} nt "
                            f"({_feature_label(_feature_ref(ref))}) y nuestra anotación "
                            f"lo cubre con {len(pieces)} features separadas."
                        ),
                        reference=_feature_ref(ref),
                        ours_features=_feature_refs(pieces),
                        metrics={
                            "reference_length_nt": ref["span_nt"],
                            "pieces": len(pieces),
                            "covered_nt": covered,
                        },
                        tags=["fragmentation"],
                    )
                )

    # 2) Varios genes de referencia fusionados en una feature nuestra.
    for ours in unmatched_ours:
        if id(ours) in consumed_ours:
            continue
        refs = [
            r
            for r in unmatched_ref
            if id(r) not in consumed_ref
            and r["seqid"] == ours["seqid"]
            and r["strand"] == ours["strand"]
            and interval_overlap(
                r["start"], r["end"], ours["start"], ours["end"]
            ) > 0
        ]
        if len(refs) >= 2:
            covered = sum(
                interval_overlap(r["start"], r["end"], ours["start"], ours["end"])
                for r in refs
            )
            if covered >= 0.5 * ours["length_nt"]:
                for r in refs:
                    consumed_ref.add(id(r))
                consumed_ours.add(id(ours))
                actionable.append(
                    _item(
                        "merged_in_ours",
                        "high",
                        cls,
                        ours["seqid"],
                        (
                            f"Nuestra anotación tiene 1 {cls} de {ours['span_nt']} nt "
                            f"({_feature_label(_feature_ref(ours))}) que solapa "
                            f"{len(refs)} genes separados de la referencia."
                        ),
                        ours=_feature_ref(ours),
                        reference_features=_feature_refs(refs),
                        metrics={
                            "ours_length_nt": ours["span_nt"],
                            "reference_features": len(refs),
                            "covered_nt": covered,
                        },
                        tags=["merge"],
                    )
                )

    # 3) Conflictos de hebra y solapamientos débiles.
    for ref in unmatched_ref:
        if id(ref) in consumed_ref:
            continue
        opposite = [
            q
            for q in unmatched_ours
            if id(q) not in consumed_ours
            and q["seqid"] == ref["seqid"]
            and q["strand"] != ref["strand"]
            and reciprocal_overlap(ref, q) >= threshold
        ]
        if opposite:
            best = max(opposite, key=lambda q: reciprocal_overlap(ref, q))
            consumed_ref.add(id(ref))
            consumed_ours.add(id(best))
            actionable.append(
                _item(
                    "strand_conflict",
                    "high",
                    cls,
                    ref["seqid"],
                    (
                        f"El mismo locus está anotado en hebra {ref['strand']} en la "
                        f"referencia y en hebra {best['strand']} en nuestra anotación."
                    ),
                    reference=_feature_ref(ref),
                    ours=_feature_ref(best),
                    metrics={"reciprocal_overlap": reciprocal_overlap(ref, best)},
                    tags=["strand"],
                )
            )
            continue
        weak = [
            (reciprocal_overlap(ref, q), q)
            for q in unmatched_ours
            if id(q) not in consumed_ours
            and q["seqid"] == ref["seqid"]
            and q["strand"] == ref["strand"]
            and 0.1 <= reciprocal_overlap(ref, q) < threshold
        ]
        if weak:
            score, best = max(weak, key=lambda pair: pair[0])
            consumed_ref.add(id(ref))
            consumed_ours.add(id(best))
            actionable.append(
                _item(
                    "weak_overlap",
                    "low",
                    cls,
                    ref["seqid"],
                    (
                        f"Solapamiento recíproco {score:.2f} (< umbral {threshold:.2f}); "
                        "puede ser un fragmento, un borde distinto o una feature no homóloga."
                    ),
                    reference=_feature_ref(ref),
                    ours=_feature_ref(best),
                    metrics={"reciprocal_overlap": round(score, 4)},
                    tags=["weak_overlap"],
                )
            )

    # 4) Genes de referencia ausentes en nuestra anotación.
    prefix = "ncrna" if cls == "ncRNA" else cls.lower()
    for ref in unmatched_ref:
        if id(ref) in consumed_ref:
            continue
        pseudo = _best_overlap(ref, pseudo_ours, threshold)
        if pseudo is not None:
            review.append(
                _item(
                    "pseudo_in_ours",
                    "low",
                    cls,
                    ref["seqid"],
                    (
                        "Nuestra anotación marca este locus como pseudogén "
                        f"(pseudo=true); la referencia lo anota como {cls}. "
                        "No se cuenta como gen ausente."
                    ),
                    reference=_feature_ref(ref),
                    ours=_feature_ref(pseudo),
                    metrics={"reciprocal_overlap": reciprocal_overlap(ref, pseudo)},
                    tags=["pseudo"],
                )
            )
            continue
        small = _is_small_orf(ref)
        severity = "high"
        if cls == "CDS":
            severity = "high" if ref["approx_aa"] >= SMALL_ORF_MAX_AA else "medium"
        product = (ref.get("product") or "").strip()
        detail = (
            f"La referencia anota un {cls} de {ref['span_nt']} nt "
            f"(~{ref['approx_aa']} aa) que no aparece en nuestra anotación."
        )
        if product:
            detail += f" product={product}."
        actionable.append(
            _item(
                f"missing_{prefix}",
                severity,
                cls,
                ref["seqid"],
                detail,
                reference=_feature_ref(ref),
                metrics={"approx_aa": ref["approx_aa"], "small_orf": small},
                tags=["missing"] + (["small_orf"] if small else []),
                small_orf=small,
            )
        )

    # 5) Features nuestras sin contraparte en la referencia.
    for ours in unmatched_ours:
        if id(ours) in consumed_ours:
            continue
        pseudo = _best_overlap(ours, pseudo_ref, threshold)
        if pseudo is not None:
            review.append(
                _item(
                    "pseudo_in_reference",
                    "low",
                    cls,
                    ours["seqid"],
                    (
                        "La referencia anota este locus como pseudogén "
                        f"(pseudo=true); nuestro {cls} no es un extra real, "
                        "pero conviene revisar marco/stop."
                    ),
                    reference=_feature_ref(pseudo),
                    ours=_feature_ref(ours),
                    metrics={"reciprocal_overlap": reciprocal_overlap(pseudo, ours)},
                    tags=["pseudo"],
                )
            )
            continue
        small = _is_small_orf(ours)
        severity = "medium"
        if cls == "CDS" and small:
            severity = "low"
        product = (ours.get("product") or "").strip()
        detail = (
            f"Nuestra anotación tiene un {cls} de {ours['span_nt']} nt "
            f"(~{ours['approx_aa']} aa) sin contraparte en la referencia."
        )
        if product:
            detail += f" product={product}."
        if small:
            detail += " Es un ORF pequeño (<150 aa): posible sobrellamada o gen real no anotado en la referencia."
        actionable.append(
            _item(
                f"extra_{prefix}",
                severity,
                cls,
                ours["seqid"],
                detail,
                ours=_feature_ref(ours),
                metrics={"approx_aa": ours["approx_aa"], "small_orf": small},
                tags=["extra"] + (["small_orf"] if small else []),
                small_orf=small,
            )
        )


def _class_absent_items(
    cls: str,
    ref_total: int,
    ours_total: int,
    actionable: list[dict[str, Any]],
    review: list[dict[str, Any]],
) -> None:
    if ref_total == 0 and ours_total > 0:
        review.append(
            _item(
                "class_absent_in_reference",
                "low",
                cls,
                None,
                (
                    f"La referencia no anota ningún {cls} ({ours_total} en nuestra "
                    "anotación). No se evalúa paridad para esta clase."
                ),
                metrics={"ours_total": ours_total, "reference_total": 0},
                tags=["reference_limitation"],
            )
        )
    elif ours_total == 0 and ref_total > 0:
        actionable.append(
            _item(
                f"missing_{'ncrna' if cls == 'ncRNA' else cls.lower()}",
                "high",
                cls,
                None,
                (
                    f"La referencia anota {ref_total} {cls} y nuestra anotación "
                    "no tiene ninguno."
                ),
                metrics={"ours_total": 0, "reference_total": ref_total},
                tags=["class_absent"],
            )
        )


def build_differences(
    *,
    parsed_ref: Mapping[str, Any] | None,
    parsed_ours: Mapping[str, Any],
    matches_by_class: Mapping[str, Sequence[Mapping[str, Any]]],
    unmatched_ref_by_class: Mapping[str, Sequence[Mapping[str, Any]]],
    unmatched_ours_by_class: Mapping[str, Sequence[Mapping[str, Any]]],
    scope_by_class: Mapping[str, tuple[Sequence[dict], Sequence[dict]]],
    shared_contigs: Sequence[str],
    threshold: float,
    boundary_alert: int,
    max_details: int,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    actionable: list[dict[str, Any]] = []
    review: list[dict[str, Any]] = []

    ref_contigs = set(parsed_ref["contigs"]) if parsed_ref else set()
    ours_contigs = set(parsed_ours["contigs"])
    if parsed_ref is not None:
        if not shared_contigs:
            actionable.append(
                _item(
                    "contig_mismatch",
                    "high",
                    "genome",
                    None,
                    (
                        "La referencia y nuestra anotación no comparten ningún "
                        "nombre de contig; no se pueden comparar coordenadas. "
                        f"Ref={sorted(ref_contigs)[:5]}; ours={sorted(ours_contigs)[:5]}."
                    ),
                    tags=["contig"],
                )
            )
        else:
            for contig in sorted(ref_contigs - set(shared_contigs)):
                actionable.append(
                    _item(
                        "contig_only_reference",
                        "medium",
                        "genome",
                        contig,
                        (
                            f"El contig {contig} está en la referencia pero no en "
                            "nuestra anotación; sus features quedan fuera del alcance."
                        ),
                        tags=["contig"],
                    )
                )
            for contig in sorted(ours_contigs - set(shared_contigs)):
                actionable.append(
                    _item(
                        "contig_only_ours",
                        "low",
                        "genome",
                        contig,
                        (
                            f"El contig {contig} está en nuestra anotación pero no "
                            "en la referencia; sus features quedan fuera del alcance."
                        ),
                        tags=["contig"],
                    )
                )

    for cls in CLASSES:
        ref_all = (
            [f for f in parsed_ref["features"] if f["class"] == cls]
            if parsed_ref
            else []
        )
        ours_all = [f for f in parsed_ours["features"] if f["class"] == cls]
        if parsed_ref is None:
            continue
        ref_scope, ours_scope = scope_by_class.get(cls, ([], []))
        if shared_contigs:
            _class_absent_items(cls, len(ref_all), len(ours_all), actionable, review)
        if not shared_contigs:
            continue
        # Si una clase falta por completo en un lado, el ítem de clase ya la
        # describe; no se enumeran decenas de features como ausentes/extra.
        if not ref_scope or not ours_scope:
            continue
        matches = list(matches_by_class.get(cls, []))
        unmatched_ref = list(unmatched_ref_by_class.get(cls, []))
        unmatched_ours = list(unmatched_ours_by_class.get(cls, []))

        # Pares emparejados: bordes, largo, estructura multiparte.
        for m in matches:
            ref, ours = m["ref"], m["query"]
            deltas = _boundary_deltas(ref, ours)
            if deltas["max_boundary_delta_nt"] >= boundary_alert:
                severity = "high" if deltas["max_boundary_delta_nt"] > 90 else "medium"
                actionable.append(
                    _item(
                        "boundary_shift",
                        severity,
                        cls,
                        ref["seqid"],
                        (
                            f"Bordes desplazados: 5' {deltas['five_prime_delta_nt']} nt, "
                            f"3' {deltas['three_prime_delta_nt']} nt respecto de la "
                            f"referencia (umbral {boundary_alert} nt)."
                        ),
                        reference=_feature_ref(ref),
                        ours=_feature_ref(ours),
                        metrics={**deltas, "reciprocal_overlap": m["reciprocal_overlap"]},
                        tags=["coordinates"],
                    )
                )
            if _is_length_discordant(ref, ours):
                ratio = _length_ratio(ref, ours)
                severity = "high" if ratio is not None and (ratio < 0.5 or ratio > 2.0) else "medium"
                actionable.append(
                    _item(
                        "length_discordance",
                        severity,
                        cls,
                        ref["seqid"],
                        (
                            f"Largo discordante: referencia {ref['span_nt']} nt vs "
                            f"nuestra {ours['span_nt']} nt (ratio {ratio})."
                        ),
                        reference=_feature_ref(ref),
                        ours=_feature_ref(ours),
                        metrics={
                            "reference_length_nt": ref["span_nt"],
                            "ours_length_nt": ours["span_nt"],
                            "length_ratio": ratio,
                        },
                        tags=["coordinates"],
                    )
                )
            if _is_multipart_structure_mismatch(m):
                actionable.append(
                    _item(
                        "multipart_structure_mismatch",
                        "medium",
                        cls,
                        ref["seqid"],
                        (
                            "El sobre de coordenadas coincide, pero los segmentos "
                            "multiparte no son idénticos."
                        ),
                        reference=_feature_ref(ref),
                        ours=_feature_ref(ours),
                        metrics={
                            "reference_segments": len(ref.get("segments") or []),
                            "ours_segments": len(ours.get("segments") or []),
                        },
                        tags=["multipart"],
                    )
                )

        _categorize_unmatched(
            cls=cls,
            unmatched_ref=unmatched_ref,
            unmatched_ours=unmatched_ours,
            pseudo_ref=[
                f
                for f in (parsed_ref.get("pseudo_features") or [])
                if f["class"] == cls
            ],
            pseudo_ours=[
                f
                for f in (parsed_ours.get("pseudo_features") or [])
                if f["class"] == cls
            ],
            actionable=actionable,
            review=review,
            threshold=threshold,
        )

        if cls == "ncRNA":
            for m in matches:
                if m["ref"]["type"] != m["query"]["type"]:
                    review.append(
                        _item(
                            "ncrna_type_mismatch",
                            "low",
                            cls,
                            m["ref"]["seqid"],
                            (
                                f"Subtipo ncRNA distinto: referencia "
                                f"{m['ref']['type']} vs nuestra {m['query']['type']}."
                            ),
                            reference=_feature_ref(m["ref"]),
                            ours=_feature_ref(m["query"]),
                            tags=["ncrna"],
                        )
                    )

    # Productos: desacuerdos de hipotética son accionables; nombres distintos
    # entre dos productos anotados son solo revisión (convención de nombres).
    if parsed_ref is not None:
        for m in matches_by_class.get("CDS", []):
            ref, ours = m["ref"], m["query"]
            ref_product = (ref.get("product") or "").strip()
            ours_product = (ours.get("product") or "").strip()
            if not ref_product or not ours_product:
                continue
            ref_hyp = is_hypothetical_product(ref_product)
            ours_hyp = is_hypothetical_product(ours_product)
            if ref_hyp != ours_hyp:
                actionable.append(
                    _item(
                        "hypothetical_disagreement",
                        "medium",
                        "CDS",
                        ref["seqid"],
                        (
                            "La referencia y nuestra anotación discrepan en si el "
                            f"producto es hipotético: ref={ref_product}; "
                            f"ours={ours_product}."
                        ),
                        reference=_feature_ref(ref),
                        ours=_feature_ref(ours),
                        metrics={"reference_hypothetical": ref_hyp, "ours_hypothetical": ours_hyp},
                        tags=["product"],
                    )
                )
                continue
            if ref_hyp and ours_hyp:
                continue
            if not _is_product_mismatch(ref_product, ours_product):
                continue
            score = product_token_jaccard(ref_product, ours_product) or 0.0
            review.append(
                    _item(
                        "product_mismatch",
                        "low",
                        "CDS",
                        ref["seqid"],
                        (
                            f"Productos con poco solapamiento de tokens "
                            f"(Jaccard {score:.2f}): ref={ref_product}; "
                            f"ours={ours_product}."
                        ),
                        reference=_feature_ref(ref),
                        ours=_feature_ref(ours),
                        metrics={"token_jaccard": round(score, 4)},
                        tags=["product_naming"],
                    )
                )

    actionable_sorted = _sort_items(actionable)
    review_sorted = _sort_items(review)
    actionable_block = {
        "total": len(actionable_sorted),
        "by_category": dict(
            sorted(Counter(i["category"] for i in actionable_sorted).items())
        ),
        "by_severity": dict(
            sorted(Counter(i["severity"] for i in actionable_sorted).items())
        ),
        "items": actionable_sorted[:max_details],
        "truncated": len(actionable_sorted) > max_details,
        "max_details": max_details,
    }
    review_block = {
        "total": len(review_sorted),
        "by_category": dict(sorted(Counter(i["category"] for i in review_sorted).items())),
        "items": review_sorted[:max_details],
        "truncated": len(review_sorted) > max_details,
        "max_details": max_details,
    }
    return actionable_block, review_block, actionable_sorted


# ---------------------------------------------------------------------------
# Orquestación
# ---------------------------------------------------------------------------
def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path | None, str, str | None]:
    root = Path(args.root).resolve()
    genome = args.genome
    ours_path = _resolve_file(args.ours, directory_patterns=("annotation.gff3", "*.gff3", "*.gff"))
    ref_path = _resolve_file(
        args.reference, directory_patterns=("*.gff3.gz", "*.gff.gz", "*.gff3", "*.gff")
    )

    if genome is None:
        if ours_path is not None:
            genome = _infer_genome(ours_path)
        elif ref_path is not None:
            genome = _infer_genome(ref_path)
        else:
            genome = DEFAULT_GENOME

    if ours_path is None:
        ours_path = discover_ours(root, genome)
    if ours_path is None:
        raise ParityError(
            "no se encontró nuestra anotación; use --ours o revise "
            f"data/annotations_ab_initio/{genome}/annotation.gff3"
        )
    if not ours_path.is_file():
        raise ParityError(f"no existe nuestra anotación: {ours_path}")

    if ref_path is None:
        ref_path = discover_reference(root, genome)

    ref_note: str | None = None
    if args.reference is not None and (ref_path is None or not ref_path.is_file()):
        ref_note = f"la ruta de referencia indicada no existe: {args.reference}"
        ref_path = None
    elif ref_path is not None and not ref_path.is_file():
        ref_note = f"la referencia autodescubierta no existe: {ref_path}"
        ref_path = None
    return ours_path, ref_path, genome, ref_note


def run_parity(args: argparse.Namespace) -> dict[str, Any]:
    """Ejecuta la comparación y devuelve el informe serializable."""
    ours_path, ref_path, genome, ref_note = _resolve_inputs(args)
    root = Path(args.root).resolve()

    try:
        ours_parsed = parse_gff(ours_path)
    except OSError as exc:
        raise ParityError(f"no se pudo leer nuestra anotación {ours_path}: {exc}") from exc

    product_source, product_label, product_path = resolve_products(
        ours_path, ours_parsed, args.ours_products
    )
    if product_path is not None:
        try:
            products = load_products_tsv(product_path)
        except OSError as exc:
            raise ParityError(f"no se pudo leer {product_path}: {exc}") from exc
        apply_products(ours_parsed, products, "tsv")
    elif product_source == "gff":
        apply_products(ours_parsed, {}, "gff")
    ours_product_stats = _product_stats(ours_parsed["features"], product_label)

    ref_parsed: dict[str, Any] | None = None
    reference_unusable_reason: str | None = None
    if ref_path is not None:
        try:
            ref_parsed = parse_gff(ref_path)
        except OSError as exc:
            reference_unusable_reason = f"no se pudo leer la referencia {ref_path}: {exc}"
            ref_parsed = None
        else:
            if ref_parsed["total_features"] == 0:
                reference_unusable_reason = (
                    "el GFF de referencia no contiene CDS/tRNA/rRNA/ncRNA comparables"
                )
                ref_parsed = None

    usable_reference = ref_parsed is not None
    if not usable_reference:
        reason = ref_note or reference_unusable_reason or (
            "no se encontró un GFF de referencia local para este genoma"
        )
        return _build_not_evaluable(
            args=args,
            root=root,
            genome=genome,
            ours_path=ours_path,
            ref_path=ref_path,
            ours_parsed=ours_parsed,
            ours_product_stats=ours_product_stats,
            product_label=product_label,
            product_source=product_source,
            reason=reason,
        )

    assert ref_parsed is not None
    shared_contigs = sorted(set(ref_parsed["contigs"]) & set(ours_parsed["contigs"]))
    shared_set = set(shared_contigs)
    threshold = float(args.overlap_threshold)
    boundary_alert = int(args.boundary_alert)
    coord_tolerance = int(args.coord_tolerance)

    summary: dict[str, Any] = {}
    matches_by_class: dict[str, list[dict[str, Any]]] = {}
    unmatched_ref_by_class: dict[str, list[dict[str, Any]]] = {}
    unmatched_ours_by_class: dict[str, list[dict[str, Any]]] = {}
    scope_by_class: dict[str, tuple[list[dict], list[dict]]] = {}
    all_matches: list[dict[str, Any]] = []

    for cls in CLASSES:
        ref_all = [f for f in ref_parsed["features"] if f["class"] == cls]
        ours_all = [f for f in ours_parsed["features"] if f["class"] == cls]
        ref_scope = [f for f in ref_all if f["seqid"] in shared_set]
        ours_scope = [f for f in ours_all if f["seqid"] in shared_set]
        scope_by_class[cls] = (ref_scope, ours_scope)
        matches = match_features(ref_scope, ours_scope, threshold) if shared_contigs else []
        unmatched_ref, unmatched_ours = _partition_matches(ref_scope, ours_scope, matches)
        matches_by_class[cls] = matches
        unmatched_ref_by_class[cls] = unmatched_ref
        unmatched_ours_by_class[cls] = unmatched_ours
        all_matches.extend(matches)
        summary[cls] = _class_metrics(
            cls,
            ref_all,
            ours_all,
            ref_scope,
            ours_scope,
            matches,
            unmatched_ref,
            unmatched_ours,
            shared_contigs=bool(shared_contigs),
            threshold=threshold,
            boundary_alert=boundary_alert,
            coord_tolerance=coord_tolerance,
        )

    actionable_block, review_block, all_actionable = build_differences(
        parsed_ref=ref_parsed,
        parsed_ours=ours_parsed,
        matches_by_class=matches_by_class,
        unmatched_ref_by_class=unmatched_ref_by_class,
        unmatched_ours_by_class=unmatched_ours_by_class,
        scope_by_class=scope_by_class,
        shared_contigs=shared_contigs,
        threshold=threshold,
        boundary_alert=boundary_alert,
        max_details=int(args.max_details),
    )

    ref_has_products = any(
        (f.get("product") or "").strip()
        for f in ref_parsed["features"]
        if f["class"] == "CDS"
    )
    ref_product_label = (
        "atributo product del GFF" if ref_has_products else "sin productos"
    )
    ref_product_stats = _product_stats(ref_parsed["features"], ref_product_label)
    paired_products = _paired_products(all_matches, product_label)
    ref_small = _small_orf_stats(
        ref_parsed["features"],
        ours_parsed["features"],
        all_matches,
        all_actionable,
        int(args.small_orf_max_aa),
    )
    nc_matches = matches_by_class.get("ncRNA", [])
    nc_stats = _ncrna_stats(
        ref_parsed["features"],
        ours_parsed["features"],
        nc_matches,
        unmatched_ref_by_class.get("ncRNA", []),
        unmatched_ours_by_class.get("ncRNA", []),
        evaluated=bool(shared_contigs)
        and bool(scope_by_class.get("ncRNA", ([], []))[0])
        and bool(scope_by_class.get("ncRNA", ([], []))[1]),
    )

    exact_total = sum(int(s["exact"] or 0) for s in summary.values() if s["status"] == "evaluated")
    matched_total = sum(int(s["matched"] or 0) for s in summary.values() if s["status"] == "evaluated")
    ref_scope_total = sum(
        int(s["reference_in_scope"]) for s in summary.values() if s["status"] == "evaluated"
    )
    ours_scope_total = sum(
        int(s["ours_in_scope"]) for s in summary.values() if s["status"] == "evaluated"
    )
    micro_precision = _ratio(matched_total, ours_scope_total)
    micro_recall = _ratio(matched_total, ref_scope_total)

    coordinates = {
        "status": "evaluated" if shared_contigs else "not_evaluable",
        "reason": None
        if shared_contigs
        else "la referencia y nuestra anotación no comparten contigs",
        "shared_contigs": shared_contigs,
        "only_reference": sorted(set(ref_parsed["contigs"]) - shared_set),
        "only_ours": sorted(set(ours_parsed["contigs"]) - shared_set),
        "matched_total": matched_total,
        "exact_total": exact_total,
        "exact_rate": _ratio(exact_total, matched_total) if matched_total else None,
        "micro_scope": "clases evaluadas (ambos lados y contigs compartidos)",
        "micro_precision": micro_precision,
        "micro_recall": micro_recall,
        "micro_f1": _f1(micro_precision, micro_recall),
        "boundary_shift_alert": sum(
            int(s["boundary_shift_alert"] or 0) for s in summary.values()
        ),
        "length_discordant": sum(
            int(s["length_discordant"] or 0) for s in summary.values()
        ),
        "multipart_structure_mismatch": sum(
            int(s["multipart_structure_mismatch"] or 0) for s in summary.values()
        ),
        "scope_excluded_reference": sum(
            len([f for f in ref_parsed["features"] if f["class"] == cls and f["seqid"] not in shared_set])
            for cls in CLASSES
        )
        if shared_contigs
        else None,
        "scope_excluded_ours": sum(
            len([f for f in ours_parsed["features"] if f["class"] == cls and f["seqid"] not in shared_set])
            for cls in CLASSES
        )
        if shared_contigs
        else None,
    }

    limitations = _limitations(
        has_reference=True,
        has_products=bool(ref_product_stats["with_product"]),
        our_product_source=product_source,
        shared_contigs=bool(shared_contigs),
    )
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "tool": "annotation_parity",
        "generated_at": _utc_now(),
        "status": "evaluated",
        "reason": None,
        "genome": genome,
        "inputs": {
            "ours": _input_metadata(root, ours_path, ours_parsed, product_label),
            "reference": _input_metadata(root, ref_path, ref_parsed, ref_product_label),
        },
        "parameters": _parameters(args, threshold, boundary_alert, coord_tolerance),
        "contigs": {
            "reference": ref_parsed["contigs"],
            "ours": ours_parsed["contigs"],
            "shared": shared_contigs,
            "only_reference": coordinates["only_reference"],
            "only_ours": coordinates["only_ours"],
            "coordinate_comparison": coordinates["status"],
        },
        "summary": summary,
        "coordinates": coordinates,
        "products": {
            "ours": ours_product_stats,
            "reference": ref_product_stats,
            "paired": paired_products,
        },
        "small_orfs": ref_small,
        "ncrna": nc_stats,
        "differences": {
            "actionable": actionable_block,
            "review": review_block,
        },
        "limitations": limitations,
        "command": _command(args),
    }
    return result


def _parameters(
    args: argparse.Namespace, threshold: float, boundary_alert: int, coord_tolerance: int
) -> dict[str, Any]:
    return {
        "overlap_threshold": threshold,
        "boundary_alert_nt": boundary_alert,
        "coordinate_tolerance_nt": coord_tolerance,
        "small_orf_max_aa": int(args.small_orf_max_aa),
        "max_details": int(args.max_details),
        "exclude_pseudo_cds": bool(args.exclude_pseudo),
        "ours_products": args.ours_products,
    }


def _input_metadata(
    root: Path,
    path: Path | None,
    parsed: Mapping[str, Any] | None,
    product_label: str,
) -> dict[str, Any] | None:
    if path is None or parsed is None:
        return None
    return {
        "path": _rel(root, path),
        "sha256": _sha256(path),
        "feature_counts": parsed["feature_counts"],
        "raw_type_counts": parsed["raw_type_counts"],
        "contigs": parsed["contigs"],
        "excluded_pseudo": parsed["excluded_pseudo"],
        "pseudo_features": len(parsed.get("pseudo_features") or []),
        "invalid_lines": parsed["invalid_lines"],
        "multipart_cds": parsed["multipart_cds"],
        "product_source": product_label,
    }


def _limitations(
    *,
    has_reference: bool,
    has_products: bool,
    our_product_source: str,
    shared_contigs: bool,
) -> list[str]:
    limits = [
        "La comparación es de anotación (coordenadas, clases y etiquetas), no de función, homología ni identidad de secuencia.",
        "El emparejamiento es 1:1 voraz por solapamiento recíproco; no es un óptimo bipartito ni un alineamiento global.",
        "Las coordenadas GFF3 son 1-based inclusive. No se hace liftover ni se alinean secuencias.",
        "El largo en aminoácidos es una estimación desde el ADN (span/3 menos un codón de parada), no una traducción.",
        "Las diferencias de nombre de producto entre bases de datos no son errores de anotación; se separan como 'revisión'.",
        "No se descargan referencias ni se ejecutan Bakta/Prokka/DFAST; solo se leen GFF locales.",
        "Un GFF de referencia puede omitir tRNA/rRNA/ncRNA o incluir pseudogenes; 'reference_absent' no indica un error nuestro.",
        "La ausencia de referencia no es un resultado: el estado es 'not_evaluable' y no se calculan precision/recall.",
    ]
    if not has_reference:
        limits.append(
            "Sin referencia local no se pueden identificar genes ausentes/extra ni medir paridad de coordenadas."
        )
    if has_reference and not has_products:
        limits.append(
            "La referencia no trae atributo 'product': no se calcula % hypothetical de referencia ni comparación pareada."
        )
    if our_product_source == "tsv":
        limits.append(
            "Nuestros productos vienen de un TSV funcional (p. ej. Swiss-Prot), no del GFF; "
            "la comparación de nombres con la referencia es orientativa."
        )
    elif our_product_source == "none":
        limits.append(
            "Nuestra anotación no trae atributo 'product' ni un TSV vecino: "
            "no se calcula % hypothetical propio ni comparación pareada de productos."
        )
    if has_reference and not shared_contigs:
        limits.append(
            "No hay contigs compartidos: los conteos son descriptivos, pero la paridad de coordenadas no es evaluable."
        )
    return limits


def _build_not_evaluable(
    *,
    args: argparse.Namespace,
    root: Path,
    genome: str,
    ours_path: Path,
    ref_path: Path | None,
    ours_parsed: Mapping[str, Any],
    ours_product_stats: Mapping[str, Any],
    product_label: str,
    product_source: str,
    reason: str,
) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for cls in CLASSES:
        ours_all = [f for f in ours_parsed["features"] if f["class"] == cls]
        summary[cls] = {
            "class": cls,
            "status": "not_evaluable",
            "reference_total": None,
            "ours_total": len(ours_all),
            "delta_ours_minus_reference": None,
            "reference_in_scope": None,
            "ours_in_scope": None,
            "matched": None,
            "recall": None,
            "precision": None,
            "f1": None,
            "exact": None,
            "exact_rate": None,
            "within_tolerance": None,
            "coordinate_tolerance_nt": int(args.coord_tolerance),
            "within_tolerance_rate": None,
            "mean_reciprocal_overlap": None,
            "median_reciprocal_overlap": None,
            "boundary_shift_alert": None,
            "boundary_alert_nt": int(args.boundary_alert),
            "length_discordant": None,
            "multipart_structure_mismatch": None,
            "unmatched_reference": None,
            "unmatched_ours": None,
            "reciprocal_overlap_threshold": float(args.overlap_threshold),
        }
    small = _small_orf_stats(
        [], ours_parsed["features"], [], [], int(args.small_orf_max_aa)
    )
    small["reference"] = None
    small["matched_small"] = None
    small["missing_small"] = None
    small["extra_small"] = None
    return {
        "schema_version": SCHEMA_VERSION,
        "tool": "annotation_parity",
        "generated_at": _utc_now(),
        "status": "not_evaluable",
        "reason": reason,
        "genome": genome,
        "inputs": {
            "ours": _input_metadata(root, ours_path, ours_parsed, product_label),
            "reference": None,
        },
        "parameters": _parameters(
            args,
            float(args.overlap_threshold),
            int(args.boundary_alert),
            int(args.coord_tolerance),
        ),
        "contigs": {
            "reference": None,
            "ours": ours_parsed["contigs"],
            "shared": [],
            "only_reference": None,
            "only_ours": ours_parsed["contigs"],
            "coordinate_comparison": "not_evaluable",
        },
        "summary": summary,
        "coordinates": {
            "status": "not_evaluable",
            "reason": "no hay referencia local utilizable",
            "shared_contigs": [],
            "only_reference": None,
            "only_ours": ours_parsed["contigs"],
            "matched_total": None,
            "exact_total": None,
            "exact_rate": None,
            "micro_scope": "no evaluable",
            "micro_precision": None,
            "micro_recall": None,
            "micro_f1": None,
            "boundary_shift_alert": None,
            "length_discordant": None,
            "multipart_structure_mismatch": None,
            "scope_excluded_reference": None,
            "scope_excluded_ours": None,
        },
        "products": {
            "ours": ours_product_stats,
            "reference": None,
            "paired": {
                "status": "not_evaluable",
                "reason": "sin referencia local",
                "matched_cds": None,
                "both_with_product": None,
                "literal_matches": None,
                "literal_agreement": None,
                "token_jaccard_mean": None,
                "token_jaccard_median": None,
                "hypothetical_disagreements": None,
                "product_mismatches": None,
                "our_product_source": product_label,
            },
        },
        "small_orfs": small,
        "ncrna": {
            "status": "not_evaluable",
            "reference_total": None,
            "ours_total": sum(
                1 for f in ours_parsed["features"] if f["class"] == "ncRNA"
            ),
            "reference_types": None,
            "ours_types": dict(
                sorted(
                    Counter(
                        f["type"]
                        for f in ours_parsed["features"]
                        if f["class"] == "ncRNA"
                    ).items()
                )
            ),
            "matched": None,
            "missing": None,
            "extra": None,
            "type_mismatches": [],
        },
        "differences": {
            "actionable": {
                "total": None,
                "by_category": {},
                "by_severity": {},
                "items": [],
                "truncated": False,
                "max_details": int(args.max_details),
                "note": "No evaluable sin referencia: no se inventan diferencias.",
            },
            "review": {
                "total": None,
                "by_category": {},
                "items": [],
                "truncated": False,
                "max_details": int(args.max_details),
            },
        },
        "limitations": _limitations(
            has_reference=False,
            has_products=False,
            our_product_source=product_source,
            shared_contigs=False,
        ),
        "command": _command(args),
    }


def _command(args: argparse.Namespace) -> str:
    argv = getattr(args, "_argv", None) or []
    rendered = " ".join(shlex.quote(str(a)) for a in argv)
    return f"python tools/annotation_parity.py {rendered}".strip()


# ---------------------------------------------------------------------------
# Render Markdown
# ---------------------------------------------------------------------------
def _table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    for row in rows:
        cells = [_escape_md(c) for c in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _status_label(status: str | None) -> str:
    return STATUS_LABELS.get(status or "", status or "—")


def render_markdown(data: Mapping[str, Any]) -> str:
    out: list[str] = []
    out.append("# Paridad de anotación base")
    out.append("")
    out.append(f"**Generado:** `{data['generated_at']}`  ")
    out.append(f"**Genoma:** `{data['genome']}`  ")
    status = data["status"]
    out.append(f"**Estado:** **{_status_label(status)}**")
    out.append("")
    if status == "not_evaluable":
        out.append("> **No evaluable.** " + str(data.get("reason") or ""))
        out.append("> No se calculan precision/recall, genes ausentes/extra ni")
        out.append("> diferencias accionables. El informe solo describe nuestra anotación.")
        out.append("")

    out.append("## Alcance")
    out.append("")
    out.append(
        "Comparación de anotación base contra un GFF3 de referencia local. "
        "Mide CDS/tRNA/rRNA/ncRNA, coordenadas, productos/hipotéticas, ORFs "
        "pequeños (<150 aa) y diferencias explicables. No mide identidad de "
        "secuencia ni función."
    )
    out.append("")
    out.append("```bash")
    out.append(data.get("command") or "python tools/annotation_parity.py")
    out.append("```")
    out.append("")

    out.append("## Entradas")
    out.append("")
    rows = []
    ours = (data.get("inputs") or {}).get("ours") or {}
    ref = (data.get("inputs") or {}).get("reference") or {}
    rows.append(
        [
            "nuestra anotación",
            f"`{ours.get('path', '—')}`",
            (ours.get("sha256") or "—")[:12],
            _fmt_int((ours.get("feature_counts") or {}).get("CDS")),
            ours.get("product_source") or "—",
        ]
    )
    if ref:
        rows.append(
            [
                "referencia",
                f"`{ref.get('path', '—')}`",
                (ref.get("sha256") or "—")[:12],
                _fmt_int((ref.get("feature_counts") or {}).get("CDS")),
                ref.get("product_source") or "—",
            ]
        )
    else:
        rows.append(["referencia", "—", "—", "—", "no disponible"])
    out.append(_table(["rol", "ruta", "sha256", "CDS", "productos"], rows))
    out.append("")

    params = data.get("parameters") or {}
    out.append(
        "Parámetros: umbral de solapamiento recíproco "
        f"**{_fmt_ratio(params.get('overlap_threshold'))}**, alerta de borde "
        f"**{_fmt_int(params.get('boundary_alert_nt'))} nt**, ORF pequeño "
        f"**<{_fmt_int(params.get('small_orf_max_aa'))} aa**, "
        f"pseudo-CDS excluidos: **{'sí' if params.get('exclude_pseudo_cds') else 'no'}**."
    )
    out.append("")

    if status == "not_evaluable":
        out.append("## Nuestra anotación (descriptivo, sin referencia)")
        out.append("")
        products = data.get("products") or {}
        ours_products = products.get("ours") or {}
        out.append(
            _table(
                ["CDS", "con product", "sin product", "hipotéticas", "% hypothetical"],
                [
                    [
                        _fmt_int(ours_products.get("cds")),
                        _fmt_int(ours_products.get("with_product")),
                        _fmt_int(ours_products.get("without_product")),
                        _fmt_int(ours_products.get("hypothetical")),
                        _fmt_pct(ours_products.get("pct_hypothetical")),
                    ]
                ],
            )
        )
        out.append("")
        small = (data.get("small_orfs") or {}).get("ours") or {}
        out.append(
            _table(
                ["CDS", "<150 aa", "% del total"],
                [
                    [
                        _fmt_int(small.get("cds")),
                        _fmt_int(small.get("small")),
                        _fmt_pct(small.get("pct")),
                    ]
                ],
            )
        )
        out.append("")
        ncrna = data.get("ncrna") or {}
        ours_types = ", ".join(
            f"{k}×{v}" for k, v in (ncrna.get("ours_types") or {}).items()
        )
        out.append(f"ncRNA propios: {ours_types or 'ninguno'}.")
        out.append("")

    if status != "not_evaluable":
        out.append("## Contigs")
        out.append("")
        contigs = data.get("contigs") or {}
        rows = [
            ["referencia", _fmt_int(len(contigs.get("reference") or [])),
             ", ".join(contigs.get("reference") or [])[:160] or "—"],
            ["nuestra", _fmt_int(len(contigs.get("ours") or [])),
             ", ".join(contigs.get("ours") or [])[:160] or "—"],
            ["compartidos", _fmt_int(len(contigs.get("shared") or [])),
             ", ".join(contigs.get("shared") or [])[:160] or "—"],
            ["solo referencia", _fmt_int(len(contigs.get("only_reference") or [])),
             ", ".join(contigs.get("only_reference") or [])[:160] or "—"],
            ["solo nuestra", _fmt_int(len(contigs.get("only_ours") or [])),
             ", ".join(contigs.get("only_ours") or [])[:160] or "—"],
        ]
        out.append(_table(["lado", "n", "contigs"], rows))
        out.append("")

    out.append("## Conteos por clase")
    out.append("")
    summary = data.get("summary") or {}
    rows = []
    for cls in CLASSES:
        block = summary.get(cls) or {}
        if status == "not_evaluable":
            rows.append(
                [
                    cls,
                    _fmt_int(block.get("ours_total")),
                    "—",
                    "—",
                    _status_label(block.get("status")),
                ]
            )
        else:
            rows.append(
                [
                    cls,
                    _fmt_int(block.get("ours_total")),
                    _fmt_int(block.get("reference_total")),
                    _fmt_int(block.get("delta_ours_minus_reference")),
                    f"{_fmt_int(block.get('ours_in_scope'))} / {_fmt_int(block.get('reference_in_scope'))}",
                    _status_label(block.get("status")),
                ]
            )
    if status == "not_evaluable":
        out.append(_table(["clase", "nuestra", "referencia", "delta", "estado"], rows))
    else:
        out.append(
            _table(
                ["clase", "nuestra", "referencia", "delta", "en alcance (nuestra/ref)", "estado"],
                rows,
            )
        )
    out.append("")

    if status != "not_evaluable":
        out.append("## Paridad de coordenadas")
        out.append("")
        rows = []
        for cls in CLASSES:
            block = summary.get(cls) or {}
            rows.append(
                [
                    cls,
                    _fmt_int(block.get("matched")),
                    _fmt_ratio(block.get("recall")),
                    _fmt_ratio(block.get("precision")),
                    _fmt_ratio(block.get("f1")),
                    _fmt_int(block.get("exact")),
                    _fmt_pct(
                        None
                        if block.get("exact_rate") is None
                        else 100 * float(block["exact_rate"])
                    ),
                    _fmt_int(block.get("boundary_shift_alert")),
                    _fmt_int(block.get("length_discordant")),
                ]
            )
        out.append(
            _table(
                [
                    "clase",
                    "emparejados",
                    "recall",
                    "precision",
                    "F1",
                    "exactos",
                    "% exactos",
                    "alertas de borde",
                    "largo discordante",
                ],
                rows,
            )
        )
        out.append("")
        coords = data.get("coordinates") or {}
        out.append(
            f"Micro recall **{_fmt_ratio(coords.get('micro_recall'))}**, "
            f"micro precision **{_fmt_ratio(coords.get('micro_precision'))}**, "
            f"micro F1 **{_fmt_ratio(coords.get('micro_f1'))}** "
            "(solo clases evaluadas); "
            f"tasa de coordenadas exactas **{_fmt_pct(None if coords.get('exact_rate') is None else 100 * float(coords['exact_rate']))}** "
            f"sobre {_fmt_int(coords.get('matched_total'))} emparejamientos."
        )
        out.append("")

        out.append("## Productos e hipotéticas")
        out.append("")
        products = data.get("products") or {}
        rows = []
        for label, key in (("nuestra", "ours"), ("referencia", "reference")):
            block = products.get(key)
            if not block:
                rows.append([label, "—", "—", "—", "—", "—", "—"])
                continue
            rows.append(
                [
                    label,
                    block.get("source") or "—",
                    _fmt_int(block.get("cds")),
                    _fmt_int(block.get("with_product")),
                    _fmt_int(block.get("without_product")),
                    _fmt_int(block.get("hypothetical")),
                    _fmt_pct(block.get("pct_hypothetical")),
                ]
            )
        out.append(
            _table(
                ["lado", "fuente", "CDS", "con product", "sin product", "hipotéticas", "% hypothetical"],
                rows,
            )
        )
        out.append("")
        paired = products.get("paired") or {}
        if paired.get("status") == "evaluated":
            out.append(
                f"Pares CDS con producto en ambos lados: **{_fmt_int(paired.get('both_with_product'))}**; "
                f"coincidencia literal **{_fmt_pct(None if paired.get('literal_agreement') is None else 100 * float(paired['literal_agreement']))}**; "
                f"Jaccard de tokens medio **{_fmt_ratio(paired.get('token_jaccard_mean'))}**, "
                f"mediana **{_fmt_ratio(paired.get('token_jaccard_median'))}**."
            )
            out.append(
                f"Desacuerdos hipotética/no-hipotética: **{_fmt_int(paired.get('hypothetical_disagreements'))}**; "
                f"posibles nombres distintos: **{_fmt_int(paired.get('product_mismatches'))}**."
            )
        else:
            out.append(
                "Comparación pareada de productos: **no evaluable** "
                "(falta producto en uno de los lados)."
            )
        out.append("")

        out.append("## ORFs pequeños (<150 aa)")
        out.append("")
        small = data.get("small_orfs") or {}
        rows = []
        for label, key in (("nuestra", "ours"), ("referencia", "reference")):
            block = small.get(key)
            if not block:
                rows.append([label, "—", "—", "—"])
                continue
            rows.append(
                [
                    label,
                    _fmt_int(block.get("cds")),
                    _fmt_int(block.get("small")),
                    _fmt_pct(block.get("pct")),
                ]
            )
        out.append(_table(["lado", "CDS", "<150 aa", "% del total"], rows))
        out.append("")
        out.append(
            f"Emparejados pequeños: **{_fmt_int(small.get('matched_small'))}**; "
            f"ausentes en nuestra anotación: **{_fmt_int(small.get('missing_small'))}**; "
            f"extra en nuestra anotación: **{_fmt_int(small.get('extra_small'))}**."
        )
        out.append("")

        out.append("## ncRNA")
        out.append("")
        ncrna = data.get("ncrna") or {}
        rows = [
            ["nuestra", _fmt_int(ncrna.get("ours_total")),
             ", ".join(f"{k}×{v}" for k, v in (ncrna.get("ours_types") or {}).items()) or "—"],
            ["referencia", _fmt_int(ncrna.get("reference_total")),
             ", ".join(f"{k}×{v}" for k, v in (ncrna.get("reference_types") or {}).items()) or "—"],
            ["estado", _status_label(ncrna.get("status")), "—"],
        ]
        out.append(_table(["lado", "total", "tipos"], rows))
        out.append("")
        if ncrna.get("type_mismatches"):
            out.append(
                f"Subtipos ncRNA distintos en pares emparejados: "
                f"**{len(ncrna['type_mismatches'])}** (ver JSON)."
            )
            out.append("")

        out.append("## Diferencias accionables")
        out.append("")
        differences = data.get("differences") or {}
        actionable = differences.get("actionable") or {}
        total = actionable.get("total")
        items = actionable.get("items") or []
        shown = min(len(items), 40)
        json_note = (
            f"el JSON guarda los primeros {_fmt_int(len(items))} (truncado)"
            if actionable.get("truncated")
            else f"el JSON incluye los {_fmt_int(len(items))}"
        )
        out.append(
            f"**Total:** **{_fmt_int(total)}**. "
            f"En este Markdown se muestran {_fmt_int(shown)}; {json_note}."
        )
        out.append("")
        by_cat = actionable.get("by_category") or {}
        if by_cat:
            out.append(_table(["categoría", "n"], sorted(by_cat.items(), key=lambda kv: (-kv[1], kv[0]))))
            out.append("")
        if items:
            rows = []
            for item in items[:40]:
                rows.append(
                    [
                        item.get("severity"),
                        item.get("category"),
                        item.get("feature_class"),
                        _feature_label(item.get("reference")) if item.get("reference") else _features_label(item.get("reference_features")),
                        _feature_label(item.get("ours")) if item.get("ours") else _features_label(item.get("ours_features")),
                        _truncate(item.get("message") or "", 220),
                    ]
                )
            out.append(_table(["severidad", "categoría", "clase", "referencia", "nuestra", "detalle"], rows))
            out.append("")
        else:
            out.append("No se detectaron diferencias accionables con los umbrales actuales.")
            out.append("")

        review = differences.get("review") or {}
        out.append("## Revisión (productos y pseudogenes)")
        out.append("")
        review_items = review.get("items") or []
        out.append(
            f"**Total:** **{_fmt_int(review.get('total'))}**. "
            "Son diferencias de nomenclatura o limitaciones de la referencia; "
            "no se cuentan como errores de anotación. "
            f"En este Markdown se muestran {_fmt_int(min(len(review_items), 20))}; "
            f"el JSON guarda {_fmt_int(len(review_items))}."
        )
        out.append("")
        if review_items:
            rows = []
            for item in review_items[:20]:
                rows.append(
                    [
                        item.get("category"),
                        _feature_label(item.get("reference")),
                        _feature_label(item.get("ours")),
                        _truncate(item.get("message") or "", 220),
                    ]
                )
            out.append(_table(["categoría", "referencia", "nuestra", "detalle"], rows))
            out.append("")

    out.append("## Límites")
    out.append("")
    for limit in data.get("limitations") or []:
        out.append(f"- {limit}")
    out.append("")
    out.append("## Reproducir")
    out.append("")
    out.append("```bash")
    out.append(data.get("command") or "python tools/annotation_parity.py")
    out.append(".venv/bin/python -m pytest tests/test_annotation_parity.py -q")
    out.append("```")
    out.append("")
    out.append(
        "Este informe no declara superioridad entre anotadores. Sin una verdad "
        "independiente y ejecuciones pareadas, solo describe diferencias."
    )
    out.append("")
    return "\n".join(out)


def write_reports(data: Mapping[str, Any], md_path: Path, json_path: Path) -> None:
    md_path = Path(md_path)
    json_path = Path(json_path)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    md_path.write_text(render_markdown(data) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compara anotación base contra un GFF3 de referencia local y "
            "reporta diferencias accionables. Sin referencia: not_evaluable."
        )
    )
    parser.add_argument(
        "--genome",
        default=None,
        help=(
            "Nombre del genoma para autodescubrir entradas y referencia. "
            f"Sin argumentos se usa {DEFAULT_GENOME}."
        ),
    )
    parser.add_argument(
        "--ours",
        default=None,
        help="GFF3 propio o directorio (por defecto se autodescubre).",
    )
    parser.add_argument(
        "--reference",
        default=None,
        help=(
            "GFF3 de referencia local (NCBI/Bakta/Prokka) o directorio. "
            "Si se omite, se autodescubre en data/annotations/."
        ),
    )
    parser.add_argument(
        "--ours-products",
        default="auto",
        help=(
            "Productos de nuestra anotación: 'auto' (GFF, luego functional.tsv "
            "vecino), 'none' o ruta a un TSV id<TAB>...<TAB>producto."
        ),
    )
    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=0.5,
        help="Solapamiento recíproco mínimo para emparejar (0-1; 0.5 por defecto).",
    )
    parser.add_argument(
        "--boundary-alert",
        type=int,
        default=30,
        help="Diferencia de borde (nt) que genera alerta de coordenadas.",
    )
    parser.add_argument(
        "--coord-tolerance",
        type=int,
        default=0,
        help="Tolerancia (nt) para considerar coordenadas equivalentes.",
    )
    parser.add_argument(
        "--small-orf-max-aa",
        type=int,
        default=150,
        help="Umbral de ORF pequeño en aminoácidos estimados (estricto <).",
    )
    parser.add_argument(
        "--max-details",
        type=int,
        default=200,
        help="Máximo de diferencias detalladas por bloque en el JSON.",
    )
    parser.add_argument(
        "--exclude-pseudo",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Excluir CDS con pseudo=true/pseudogene=true (por defecto sí).",
    )
    parser.add_argument(
        "--root",
        default=str(ROOT),
        help="Raíz del proyecto para autodescubrir datos (por defecto el repo).",
    )
    parser.add_argument("--md", type=Path, default=None, help="Salida Markdown.")
    parser.add_argument("--json", type=Path, default=None, help="Salida JSON.")
    parser.add_argument("--quiet", action="store_true", help="No imprimir resumen.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = argv  # type: ignore[attr-defined]
    if args.genome is None and args.ours is None and args.reference is None:
        args.genome = DEFAULT_GENOME
    root = Path(args.root).resolve()
    md_path = args.md or (root / "reports" / "annotation_parity.md")
    json_path = args.json or (root / "reports" / "annotation_parity.json")
    try:
        data = run_parity(args)
    except ParityError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    write_reports(data, md_path, json_path)
    if not args.quiet:
        if data["status"] == "not_evaluable":
            print(f"estado: not_evaluable — {data.get('reason')}")
        else:
            actionable = (data.get("differences") or {}).get("actionable") or {}
            print(
                f"estado: evaluated | CDS ours={data['summary']['CDS']['ours_total']} "
                f"ref={data['summary']['CDS']['reference_total']} "
                f"matched={data['summary']['CDS']['matched']} "
                f"actionable={actionable.get('total')}"
            )
        print(f"md: {md_path}")
        print(f"json: {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
