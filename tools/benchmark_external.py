#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark EXTERNO contra anotadores bacterianos, sin red ni descargas.

Qué hace
--------
- Detecta ``bakta``, ``prokka`` y ``dfast`` con :func:`shutil.which`.
- Si el binario no existe, la capa se declara ``not_evaluable`` con la razón.
- Compara un GFF de referencia (gold local) con el GFF producido por la
  herramienta: CDS, tRNA y rRNA por coordenadas/hebra, más atributos y
  porcentaje de proteínas hipotéticas si el GFF trae ``product``.
- Solo ejecuta una herramienta si el usuario pasa ``--run`` y el binario
  existe. La ejecución ocurre en un directorio temporal.
- Mide tiempo de pared y RSS máximo del proceso hijo cuando hay ejecución.
- No descarga bases de datos ni consulta servicios web. No modifica pipeline,
  GUI ni app.

Uso
---
    python tools/benchmark_external.py \
        --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
        --reference-gff data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
        --tool all

    # Ejecución externa explícita (solo si el binario está instalado):
    python tools/benchmark_external.py ... --tool prokka --run

Salidas por defecto:
    reports/external_benchmark.md
    reports/external_benchmark.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MD = ROOT / "reports" / "external_benchmark.md"
DEFAULT_JSON = ROOT / "reports" / "external_benchmark.json"

CATEGORIES = ("CDS", "tRNA", "rRNA")
TOOL_SPECS: dict[str, dict[str, str]] = {
    "bakta": {"binary": "bakta", "label": "Bakta"},
    "prokka": {"binary": "prokka", "label": "Prokka"},
    "dfast": {"binary": "dfast", "label": "DFAST"},
}
ALL_TOOLS = tuple(TOOL_SPECS)
TYPE_MAP = {"cds": "CDS", "trna": "tRNA", "rrna": "rRNA"}

HYPOTHETICAL_RE = re.compile(
    r"\b(?:hypothetical|uncharacterized|unknown\s+protein|"
    r"protein\s+of\s+unknown\s+function|duf\d+|upf\d+)\b",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Utilidades
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


def _mean(values: Iterable[float]) -> float | None:
    vals = list(values)
    if not vals:
        return None
    return round(sum(vals) / len(vals), 4)


def _tail(value: Any, limit: int = 2000) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    text = str(value)
    return text[-limit:]


def _max_rss_kb() -> int | None:
    """RSS máximo acumulado de procesos hijos (KiB en Linux)."""
    try:
        import resource  # noqa: PLC0415 — opcional según plataforma

        return int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    except (ImportError, ValueError, OSError):
        return None


def parse_attributes(text: str) -> dict[str, str]:
    """Parsea atributos GFF3 ``clave=valor;...`` con unquote y sin romper valores."""
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


def _type_name(raw: str) -> str:
    raw = (raw or "").strip()
    return TYPE_MAP.get(raw.lower(), raw)


# ---------------------------------------------------------------------------
# Parseo de GFF3
# ---------------------------------------------------------------------------
def parse_gff(
    path: Path,
    categories: Sequence[str] = CATEGORIES,
    exclude_pseudo: bool = True,
) -> dict[str, Any]:
    """Devuelve features y conteos de un GFF3.

    Las coordenadas se conservan 1-based inclusive, tal como vienen en GFF3.
    """
    path = Path(path)
    wanted = set(categories)
    features: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    sequence_regions: dict[str, int] = {}
    excluded_pseudo = 0
    invalid = 0

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
                invalid += 1
                continue
            ftype = _type_name(cols[2])
            if ftype not in wanted:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                invalid += 1
                continue
            if start > end:
                start, end = end, start
            attrs = parse_attributes(cols[8])
            if exclude_pseudo and ftype == "CDS" and _is_pseudo(attrs):
                excluded_pseudo += 1
                continue
            strand = cols[6] if cols[6] in {"+", "-"} else "."
            features.append({
                "seqid": cols[0],
                "type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "source": cols[1],
                "attrs": attrs,
                "line": line_no,
            })
            counts[ftype] += 1

    contigs = sorted({f["seqid"] for f in features})
    return {
        "path": str(path),
        "features": features,
        "counts": {cat: int(counts.get(cat, 0)) for cat in categories},
        "contigs": contigs,
        "sequence_regions": sequence_regions,
        "excluded_pseudo": excluded_pseudo,
        "invalid_lines": invalid,
    }


# ---------------------------------------------------------------------------
# Comparación de features
# ---------------------------------------------------------------------------
def _length(feature: Mapping[str, Any]) -> int:
    return int(feature["end"]) - int(feature["start"]) + 1


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Longitud de intersección 1-based inclusive; 0 si no se tocan."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def reciprocal_overlap(a: Mapping[str, Any], b: Mapping[str, Any]) -> float:
    """Mínimo de las dos fracciones de solapamiento (0.0 a 1.0)."""
    ov = interval_overlap(a["start"], a["end"], b["start"], b["end"])
    if ov <= 0:
        return 0.0
    return min(ov / _length(a), ov / _length(b))


def match_features(
    reference_features: Sequence[Mapping[str, Any]],
    query_features: Sequence[Mapping[str, Any]],
    threshold: float = 0.5,
) -> list[dict[str, Any]]:
    """Emparejamiento 1:1 voraz por solapamiento recíproco, por tipo/contig/hebra.

    Es determinista: ordena por ``(-overlap, start_ref, end_ref, start_query,
    end_query)``. No es un alineamiento global ni un óptimo bipartito; se
    documenta para que las cifras sean reproducibles.
    """
    ref_by_scope: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    query_by_scope: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for feature in reference_features:
        ref_by_scope[(feature["type"], feature["seqid"], feature["strand"])].append(feature)
    for feature in query_features:
        query_by_scope[(feature["type"], feature["seqid"], feature["strand"])].append(feature)

    matches: list[dict[str, Any]] = []
    scopes = sorted(set(ref_by_scope) | set(query_by_scope))
    for scope in scopes:
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
        pairs.sort(key=lambda p: (
            -p[0],
            refs[p[1]]["start"], refs[p[1]]["end"],
            queries[p[2]]["start"], queries[p[2]]["end"],
        ))
        used_ref: set[int] = set()
        used_query: set[int] = set()
        for score, i, j in pairs:
            if i in used_ref or j in used_query:
                continue
            used_ref.add(i)
            used_query.add(j)
            ref, query = refs[i], queries[j]
            matches.append({
                "ref": ref,
                "query": query,
                "reciprocal_overlap": round(score, 4),
                "exact": bool(
                    ref["start"] == query["start"]
                    and ref["end"] == query["end"]
                    and ref["strand"] == query["strand"]
                ),
            })
    return matches


def _prf(tp: int, fp: int, fn: int) -> dict[str, Any]:
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    if precision is None or recall is None or (precision + recall) == 0:
        f1 = None
    else:
        f1 = round(2 * precision * recall / (precision + recall), 4)
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


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
    """Clasificador conservador de ``product``; no usa 'putative' a secas."""
    text = (product or "").strip()
    if not text:
        return False
    return bool(HYPOTHETICAL_RE.search(text))


def hypothetical_stats(features: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    cds = [f for f in features if f["type"] == "CDS"]
    with_product = [
        f for f in cds
        if (f["attrs"].get("product") or "").strip()
    ]
    hypothetical = [
        f for f in with_product
        if is_hypothetical_product(f["attrs"].get("product", ""))
    ]
    return {
        "cds": len(cds),
        "with_product": len(with_product),
        "without_product": len(cds) - len(with_product),
        "hypothetical": len(hypothetical),
        "pct_hypothetical": _pct(len(hypothetical), len(with_product)),
    }


def _attributes_for_pair(ref: Mapping[str, Any], query: Mapping[str, Any]) -> dict[str, Any]:
    ra = ref.get("attrs") or {}
    qa = query.get("attrs") or {}
    return {
        "product_ref": (ra.get("product") or "").strip(),
        "product_query": (qa.get("product") or "").strip(),
        "gene_ref": (ra.get("gene") or "").strip(),
        "gene_query": (qa.get("gene") or "").strip(),
        "locus_ref": (ra.get("locus_tag") or "").strip(),
        "locus_query": (qa.get("locus_tag") or "").strip(),
    }


def attribute_stats(
    matches: Sequence[Mapping[str, Any]],
    query_features: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    pairs = 0
    ref_product = 0
    query_product = 0
    literal_product = 0
    token_scores: list[float] = []
    ref_gene = 0
    query_gene = 0
    literal_gene = 0
    ref_locus = 0
    query_locus = 0
    literal_locus = 0

    for match in matches:
        vals = _attributes_for_pair(match["ref"], match["query"])
        pairs += 1
        if vals["product_ref"]:
            ref_product += 1
        if vals["product_query"]:
            query_product += 1
        if vals["product_ref"] and vals["product_query"]:
            if normalize_product(vals["product_ref"]) == normalize_product(vals["product_query"]):
                literal_product += 1
            score = product_token_jaccard(vals["product_ref"], vals["product_query"])
            if score is not None:
                token_scores.append(score)
        if vals["gene_ref"]:
            ref_gene += 1
        if vals["gene_query"]:
            query_gene += 1
        if vals["gene_ref"] and vals["gene_query"]:
            if normalize_product(vals["gene_ref"]) == normalize_product(vals["gene_query"]):
                literal_gene += 1
        if vals["locus_ref"]:
            ref_locus += 1
        if vals["locus_query"]:
            query_locus += 1
        if vals["locus_ref"] and vals["locus_query"]:
            if vals["locus_ref"] == vals["locus_query"]:
                literal_locus += 1

    query_cov = {}
    for key in ("product", "gene", "locus_tag"):
        n = sum(1 for f in query_features if (f["attrs"].get(key) or "").strip())
        query_cov[key] = {"n": n, "pct": _pct(n, len(query_features))}
    ref_cov = {}
    refs = [m["ref"] for m in matches]
    for key in ("product", "gene", "locus_tag"):
        n = sum(1 for f in refs if (f["attrs"].get(key) or "").strip())
        ref_cov[key] = {"n": n, "pct": _pct(n, len(refs))}

    both_products = sum(
        1 for m in matches
        if (m["ref"]["attrs"].get("product") or "").strip()
        and (m["query"]["attrs"].get("product") or "").strip()
    )
    both_genes = sum(
        1 for m in matches
        if (m["ref"]["attrs"].get("gene") or "").strip()
        and (m["query"]["attrs"].get("gene") or "").strip()
    )
    both_locus = sum(
        1 for m in matches
        if (m["ref"]["attrs"].get("locus_tag") or "").strip()
        and (m["query"]["attrs"].get("locus_tag") or "").strip()
    )

    return {
        "matched_pairs": pairs,
        "reference_attribute_coverage": ref_cov,
        "query_attribute_coverage": query_cov,
        "product": {
            "both_present": both_products,
            "literal_matches": literal_product,
            "literal_agreement": _ratio(literal_product, both_products),
            "token_jaccard_mean": _mean(token_scores),
        },
        "gene": {
            "both_present": both_genes,
            "literal_matches": literal_gene,
            "literal_agreement": _ratio(literal_gene, both_genes),
        },
        "locus_tag": {
            "both_present": both_locus,
            "literal_matches": literal_locus,
            "literal_agreement": _ratio(literal_locus, both_locus),
        },
    }


def compare_gffs(
    reference_gff: Path,
    query_gff: Path,
    threshold: float = 0.5,
    exclude_pseudo: bool = True,
    query_label: str | None = None,
) -> dict[str, Any]:
    """Compara un GFF de referencia (gold) con el GFF de una herramienta."""
    if not 0.0 < threshold <= 1.0:
        raise ValueError("threshold debe estar en (0, 1]")
    try:
        ref = parse_gff(reference_gff, exclude_pseudo=exclude_pseudo)
        query = parse_gff(query_gff, exclude_pseudo=exclude_pseudo)
    except OSError as exc:
        return _not_evaluable_comparison(
            reference_gff, query_gff, f"No se pudo leer un GFF: {exc}"
        )

    if not any(ref["counts"].values()):
        return _not_evaluable_comparison(
            reference_gff, query_gff,
            "El GFF de referencia no tiene features CDS/tRNA/rRNA evaluables.",
        )
    if not any(query["counts"].values()):
        return _not_evaluable_comparison(
            reference_gff, query_gff,
            "El GFF de la herramienta no tiene features CDS/tRNA/rRNA.",
        )
    common_contigs = sorted(set(ref["contigs"]) & set(query["contigs"]))
    if not common_contigs:
        return _not_evaluable_comparison(
            reference_gff, query_gff,
            "No hay contigs en común entre referencia y herramienta.",
        )

    matches = match_features(ref["features"], query["features"], threshold=threshold)
    categories: dict[str, Any] = {}
    for cat in CATEGORIES:
        cat_matches = [m for m in matches if m["ref"]["type"] == cat]
        ref_n = ref["counts"].get(cat, 0)
        query_n = query["counts"].get(cat, 0)
        tp = len(cat_matches)
        metrics = _prf(tp, query_n - tp, ref_n - tp)
        overlaps = [m["reciprocal_overlap"] for m in cat_matches]
        cat_query_features = [f for f in query["features"] if f["type"] == cat]
        categories[cat] = {
            "reference_n": ref_n,
            "query_n": query_n,
            **metrics,
            "exact_matches": sum(1 for m in cat_matches if m["exact"]),
            "mean_reciprocal_overlap": _mean(overlaps),
            "attributes": attribute_stats(cat_matches, cat_query_features),
            "hypothetical": {
                "reference": hypothetical_stats(
                    [f for f in ref["features"] if f["type"] == cat]
                ),
                "query": hypothetical_stats(
                    [f for f in query["features"] if f["type"] == cat]
                ),
            },
        }

    ref_total = sum(ref["counts"].values())
    query_total = sum(query["counts"].values())
    tp_total = len(matches)
    aggregate = {
        "reference_n": ref_total,
        "query_n": query_total,
        **_prf(tp_total, query_total - tp_total, ref_total - tp_total),
        "exact_matches": sum(1 for m in matches if m["exact"]),
        "mean_reciprocal_overlap": _mean(m["reciprocal_overlap"] for m in matches),
        "hypothetical_reference": categories["CDS"]["hypothetical"]["reference"],
        "hypothetical_query": categories["CDS"]["hypothetical"]["query"],
    }

    not_evaluable: list[dict[str, str]] = []
    if categories["CDS"]["attributes"]["query_attribute_coverage"]["product"]["n"] == 0:
        not_evaluable.append({
            "metric": "product y % hypothetical",
            "reason": (
                "El GFF de la herramienta no tiene atributo 'product' en ningún "
                "CDS; no se puede calcular el porcentaje de hipotéticas."
            ),
        })
    if aggregate["tp"] == 0:
        not_evaluable.append({
            "metric": "coincidencia de features",
            "reason": (
                "Ningún feature supera el umbral de solapamiento recíproco; "
                "revisar contigs/hebra o el umbral antes de interpretar."
            ),
        })

    comparison = {
        "status": "partial" if not_evaluable else "ok",
        "reference_gff": _rel(ROOT, reference_gff),
        "query_gff": query_label or _rel(ROOT, query_gff),
        "threshold": threshold,
        "exclude_pseudo": exclude_pseudo,
        "reference_counts": ref["counts"],
        "query_counts": query["counts"],
        "reference_contigs": len(ref["contigs"]),
        "query_contigs": len(query["contigs"]),
        "common_contigs": len(common_contigs),
        "excluded_pseudo_reference": ref["excluded_pseudo"],
        "excluded_pseudo_query": query["excluded_pseudo"],
        "categories": categories,
        "aggregate": aggregate,
        "not_evaluable": not_evaluable,
        "limitations": [
            "La referencia es el gold de esta comparación; no se evalúa función, dominios ni rutas.",
            "El emparejamiento es 1:1 voraz por solapamiento recíproco; no es un alineamiento global.",
            "La coincidencia literal de 'product' no es accuracy: los nombres dependen de la base y la convención de cada herramienta.",
            "CDS, tRNA y rRNA se comparan por separado; no se comparan ncRNA, tmRNA ni pseudogenes (si se excluyen).",
            "Las coordenadas se asumen 1-based inclusive y se comparan por contig y hebra.",
        ],
    }
    return comparison


def _not_evaluable_comparison(
    reference_gff: Path,
    query_gff: Path,
    reason: str,
) -> dict[str, Any]:
    return {
        "status": "not_evaluable",
        "reference_gff": _rel(ROOT, reference_gff),
        "query_gff": _rel(ROOT, query_gff),
        "categories": {},
        "aggregate": {},
        "not_evaluable": [{"metric": "comparación de GFF", "reason": reason}],
        "limitations": [reason],
    }


# ---------------------------------------------------------------------------
# Detección y ejecución de herramientas
# ---------------------------------------------------------------------------
def detect_tool(name: str, which=shutil.which) -> dict[str, Any]:
    spec = TOOL_SPECS[name]
    path = which(spec["binary"])
    return {
        "tool": name,
        "label": spec["label"],
        "binary": spec["binary"],
        "available": bool(path),
        "path": path,
        "reason": (
            "" if path
            else f"binario '{spec['binary']}' no encontrado en PATH (shutil.which)"
        ),
    }


def detect_tools(which=None) -> dict[str, Any]:
    which = which or shutil.which
    return {name: detect_tool(name, which) for name in ALL_TOOLS}


def build_command(
    tool: str,
    binary_path: str,
    genome: Path,
    outdir: Path,
    prefix: str = "result",
    bakta_db: Path | None = None,
) -> list[str]:
    genome_s = str(Path(genome).resolve())
    out_s = str(Path(outdir).resolve())
    if tool == "bakta":
        cmd = [binary_path, "--output", out_s, "--prefix", prefix,
               "--threads", "1", "--force"]
        if bakta_db:
            cmd.extend(["--db", str(Path(bakta_db).resolve())])
        cmd.append(genome_s)
        return cmd
    if tool == "prokka":
        return [binary_path, "--outdir", out_s, "--prefix", prefix,
                "--cpus", "1", "--force", genome_s]
    if tool == "dfast":
        return [binary_path, "--genome", genome_s, "--out", out_s, "--cpu", "1"]
    raise ValueError(f"herramienta desconocida: {tool}")


def execute_command(cmd: Sequence[str], timeout: float, cwd: Path | None = None) -> dict[str, Any]:
    """Ejecuta un comando y mide tiempo de pared y RSS máximo (aproximado)."""
    start = time.perf_counter()
    rss_before = _max_rss_kb()
    try:
        proc = subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        elapsed = time.perf_counter() - start
        rss_after = _max_rss_kb()
        return {
            "status": "ok" if proc.returncode == 0 else "failed",
            "command": list(cmd),
            "returncode": proc.returncode,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": rss_after,
            "max_rss_mib": round(rss_after / 1024.0, 2) if rss_after is not None else None,
            "max_rss_delta_kb": (
                max(0, rss_after - rss_before)
                if rss_after is not None and rss_before is not None else None
            ),
            "rss_note": (
                "ru_maxrss de RUSAGE_CHILDREN (KiB en Linux); es el máximo "
                "acumulado de procesos hijos, no un perfil por instante."
            ),
            "stdout_tail": _tail(proc.stdout),
            "stderr_tail": _tail(proc.stderr),
        }
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - start
        rss_after = _max_rss_kb()
        return {
            "status": "timeout",
            "command": list(cmd),
            "returncode": None,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": rss_after,
            "max_rss_mib": round(rss_after / 1024.0, 2) if rss_after is not None else None,
            "max_rss_delta_kb": None,
            "rss_note": "Timeout durante la ejecución externa.",
            "stdout_tail": _tail(exc.stdout),
            "stderr_tail": _tail(exc.stderr),
        }
    except OSError as exc:
        elapsed = time.perf_counter() - start
        return {
            "status": "failed",
            "command": list(cmd),
            "returncode": None,
            "runtime_seconds": round(elapsed, 3),
            "max_rss_kb": _max_rss_kb(),
            "max_rss_mib": None,
            "max_rss_delta_kb": None,
            "rss_note": "No se pudo iniciar el proceso externo.",
            "stdout_tail": "",
            "stderr_tail": _tail(str(exc)),
        }


def find_output_gff(outdir: Path) -> tuple[Path | None, list[str]]:
    candidates = sorted(
        {p for p in Path(outdir).rglob("*.gff3") if p.is_file()}
        | {p for p in Path(outdir).rglob("*.gff") if p.is_file()}
    )
    names = [p.name for p in candidates]
    return (candidates[0] if candidates else None), names


def benchmark_tool(
    tool: str,
    genome: Path | None,
    reference_gff: Path | None,
    run: bool = False,
    which=None,
    threshold: float = 0.5,
    timeout: float = 1800.0,
    bakta_db: Path | None = None,
) -> dict[str, Any]:
    """Detecta y, solo si ``run``, ejecuta una herramienta y compara su GFF."""
    which = which or shutil.which
    spec = TOOL_SPECS[tool]
    binary_path = which(spec["binary"])
    result: dict[str, Any] = {
        "tool": tool,
        "label": spec["label"],
        "binary": spec["binary"],
        "available": bool(binary_path),
        "binary_path": binary_path,
        "run": None,
        "comparison": None,
        "not_evaluable": [],
        "limitations": [],
    }
    if not binary_path:
        reason = f"binario '{spec['binary']}' no encontrado en PATH (shutil.which)"
        result.update({"status": "not_evaluable", "reason": reason})
        result["not_evaluable"].append({"metric": "ejecución externa", "reason": reason})
        result["limitations"].append(
            "No se ejecutó ninguna herramienta: el binario no está instalado."
        )
        return result

    if not run:
        reason = "binario disponible, pero no se pidió --run"
        result.update({"status": "not_run", "reason": reason})
        result["run"] = {"status": "not_run", "reason": reason}
        result["not_evaluable"].append({"metric": "ejecución externa", "reason": reason})
        return result

    if genome is None or not Path(genome).is_file():
        reason = "falta --genome o el archivo no existe; no se ejecuta la herramienta"
        result.update({"status": "not_evaluable", "reason": reason})
        result["not_evaluable"].append({"metric": "ejecución externa", "reason": reason})
        return result

    with tempfile.TemporaryDirectory(prefix=f"extbench_{tool}_") as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = build_command(tool, binary_path, genome, outdir, bakta_db=bakta_db)
        run_info = execute_command(cmd, timeout=timeout, cwd=tmpdir)
        gff_path, gff_names = find_output_gff(outdir)
        run_info["gff_candidates"] = gff_names
        if gff_path is not None:
            run_info["gff"] = str(gff_path.relative_to(tmpdir))
        result["run"] = run_info

        if run_info["status"] != "ok":
            result["status"] = "failed"
            result["reason"] = f"la ejecución externa terminó con estado {run_info['status']}"
            result["not_evaluable"].append({
                "metric": "comparación de GFF",
                "reason": result["reason"],
            })
            return result
        if gff_path is None:
            result["status"] = "run_ok_no_gff"
            result["reason"] = "la herramienta terminó pero no se encontró un GFF de salida"
            result["not_evaluable"].append({
                "metric": "comparación de GFF",
                "reason": result["reason"],
            })
            return result
        if reference_gff is None or not Path(reference_gff).is_file():
            result["status"] = "run_ok_no_reference"
            result["reason"] = "no hay --reference-gff local para comparar"
            result["not_evaluable"].append({
                "metric": "comparación de GFF",
                "reason": result["reason"],
            })
            return result

        comparison = compare_gffs(
            reference_gff,
            gff_path,
            threshold=threshold,
            query_label=f"<temp>/{gff_path.name}",
        )
        result["comparison"] = comparison
        result["status"] = comparison["status"]
        result["not_evaluable"].extend(comparison.get("not_evaluable", []))
        result["limitations"].extend(comparison.get("limitations", []))
        return result


# ---------------------------------------------------------------------------
# Informe
# ---------------------------------------------------------------------------
def summarize_reference(path: Path | None, exclude_pseudo: bool = True) -> dict[str, Any]:
    if path is None:
        return {"status": "not_evaluable", "reason": "no se pasó --reference-gff"}
    if not Path(path).is_file():
        return {"status": "not_evaluable", "reason": f"no existe {path}"}
    try:
        parsed = parse_gff(path, exclude_pseudo=exclude_pseudo)
    except OSError as exc:
        return {"status": "not_evaluable", "reason": f"no se pudo leer {path}: {exc}"}
    if not any(parsed["counts"].values()):
        return {
            "status": "not_evaluable",
            "reason": "el GFF de referencia no tiene CDS/tRNA/rRNA",
        }
    return {
        "status": "ok",
        "path": _rel(ROOT, path),
        "counts": parsed["counts"],
        "contigs": len(parsed["contigs"]),
        "sequence_regions": len(parsed["sequence_regions"]),
        "excluded_pseudo": parsed["excluded_pseudo"],
        "invalid_lines": parsed["invalid_lines"],
        "hypothetical": hypothetical_stats(parsed["features"]),
    }


def build_data(
    genome: Path | None = None,
    reference_gff: Path | None = None,
    tool: str = "all",
    run: bool = False,
    which=None,
    threshold: float = 0.5,
    timeout: float = 1800.0,
    bakta_db: Path | None = None,
) -> dict[str, Any]:
    if tool not in {"all", *ALL_TOOLS}:
        raise ValueError(f"herramienta desconocida: {tool}")
    selected = ALL_TOOLS if tool == "all" else (tool,)
    detection = detect_tools(which=which)
    reference = summarize_reference(reference_gff)

    results: dict[str, Any] = {}
    for name in selected:
        results[name] = benchmark_tool(
            name,
            genome=genome,
            reference_gff=reference_gff,
            run=run,
            which=which,
            threshold=threshold,
            timeout=timeout,
            bakta_db=bakta_db,
        )

    global_ne: list[dict[str, str]] = []
    for name, result in results.items():
        for item in result.get("not_evaluable", []):
            global_ne.append({
                "tool": name,
                "metric": item.get("metric", ""),
                "reason": item.get("reason", ""),
            })
    if reference.get("status") == "not_evaluable":
        global_ne.append({
            "tool": "referencia",
            "metric": "gold de comparación",
            "reason": reference.get("reason", ""),
        })

    limitations = [
        "No se descargó nada: el script solo lee datos locales y, con --run, ejecuta un binario ya instalado.",
        "Una herramienta ausente no es un cero ni un fracaso del método: queda como not_evaluable.",
        "La referencia GFF es el gold de esta comparación; no se evalúa función, rutas ni dominios.",
        "Los porcentajes de hipotéticas dependen de cómo cada herramienta nombra 'product'.",
        "No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.",
    ]

    return {
        "schema": "external-benchmark/v1",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "project_root": str(ROOT),
        "scope": (
            "Benchmark externo con datos locales, sin red ni descargas. "
            "La ejecución de terceros solo ocurre con --run."
        ),
        "parameters": {
            "tool": tool,
            "run": bool(run),
            "genome": _rel(ROOT, genome),
            "reference_gff": _rel(ROOT, reference_gff),
            "overlap_threshold": threshold,
            "timeout_seconds": timeout,
            "exclude_pseudo": True,
            "bakta_db": _rel(ROOT, bakta_db),
        },
        "environment": {
            "tools": detection,
            "available_tools": [n for n, d in detection.items() if d["available"]],
            "missing_tools": [n for n, d in detection.items() if not d["available"]],
        },
        "reference": reference,
        "results": results,
        "not_evaluable": global_ne,
        "limitations": limitations,
    }


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _fmt_int(value: Any) -> str:
    if value is None:
        return "—"
    try:
        return f"{int(value):,}".replace(",", ".")
    except (TypeError, ValueError):
        return str(value)


def _fmt_num(value: Any, ndigits: int = 2) -> str:
    if value is None:
        return "—"
    try:
        text = f"{float(value):,.{ndigits}f}"
    except (TypeError, ValueError):
        return str(value)
    return text.replace(",", "X").replace(".", ",").replace("X", ".")


def _fmt_pct(value: Any) -> str:
    return "—" if value is None else f"{_fmt_num(value, 1)} %"


def _fmt_ratio(value: Any) -> str:
    if value is None:
        return "—"
    return _fmt_num(value, 3)


def _yes_no(value: Any) -> str:
    return "sí" if value else "no"


def _table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(str(h) for h in headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return "\n".join(lines)


def _status_label(status: str | None) -> str:
    return {
        "ok": "ok",
        "partial": "parcial",
        "not_evaluable": "no evaluable",
        "not_run": "no ejecutado",
        "failed": "falló",
        "run_ok_no_gff": "ejecutó sin GFF",
        "run_ok_no_reference": "ejecutó sin referencia",
    }.get(status or "", status or "—")


def _category_rows(comparison: Mapping[str, Any]) -> list[list[str]]:
    rows = []
    for cat, block in (comparison.get("categories") or {}).items():
        rows.append([
            cat,
            _fmt_int(block.get("reference_n")),
            _fmt_int(block.get("query_n")),
            _fmt_int(block.get("tp")),
            _fmt_int(block.get("fp")),
            _fmt_int(block.get("fn")),
            _fmt_ratio(block.get("precision")),
            _fmt_ratio(block.get("recall")),
            _fmt_ratio(block.get("f1")),
            _fmt_int(block.get("exact_matches")),
            _fmt_ratio(block.get("mean_reciprocal_overlap")),
        ])
    return rows


def _attributes_rows(comparison: Mapping[str, Any]) -> list[list[str]]:
    rows = []
    for cat, block in (comparison.get("categories") or {}).items():
        attrs = block.get("attributes") or {}
        rows.append([
            cat,
            _fmt_int(attrs.get("matched_pairs")),
            _fmt_int((attrs.get("product") or {}).get("both_present")),
            _fmt_ratio((attrs.get("product") or {}).get("literal_agreement")),
            _fmt_ratio((attrs.get("product") or {}).get("token_jaccard_mean")),
            _fmt_ratio((attrs.get("gene") or {}).get("literal_agreement")),
            _fmt_ratio((attrs.get("locus_tag") or {}).get("literal_agreement")),
        ])
    return rows


def _hypothetical_rows(comparison: Mapping[str, Any]) -> list[list[str]]:
    rows = []
    for cat, block in (comparison.get("categories") or {}).items():
        hyp = block.get("hypothetical") or {}
        ref = hyp.get("reference") or {}
        query = hyp.get("query") or {}
        rows.append([
            cat,
            _fmt_int(ref.get("cds")),
            _fmt_int(ref.get("with_product")),
            _fmt_int(ref.get("hypothetical")),
            _fmt_pct(ref.get("pct_hypothetical")),
            _fmt_int(query.get("cds")),
            _fmt_int(query.get("with_product")),
            _fmt_int(query.get("hypothetical")),
            _fmt_pct(query.get("pct_hypothetical")),
        ])
    return rows


def _run_line(run: Mapping[str, Any] | None) -> str:
    if not run:
        return "—"
    status = _status_label(run.get("status"))
    if run.get("status") == "ok" or run.get("status") == "failed":
        return (
            f"{status}; {_fmt_num(run.get('runtime_seconds'), 3)} s; "
            f"{_fmt_num(run.get('max_rss_mib'), 1)} MiB RSS"
        )
    if run.get("status") == "timeout":
        return f"timeout; {_fmt_num(run.get('runtime_seconds'), 3)} s"
    return status


def render_markdown(data: Mapping[str, Any]) -> str:
    params = data.get("parameters") or {}
    env = data.get("environment") or {}
    reference = data.get("reference") or {}
    results = data.get("results") or {}
    lines: list[str] = [
        "# Benchmark externo contra anotadores",
        "",
        f"**Generado:** {data.get('generated_at', '—')}",
        "",
        f"**Alcance:** {data.get('scope', '—')}",
        "",
        "**Comando de ejemplo:**",
        "",
        "```bash",
        "python tools/benchmark_external.py \\",
        f"  --genome {params.get('genome') or '<genoma.fna>'} \\",
        f"  --reference-gff {params.get('reference_gff') or '<referencia.gff.gz>'} \\",
        f"  --tool {params.get('tool', 'all')}",
        "# añadir --run para ejecutar un binario ya instalado (sin red)",
        "```",
        "",
        "Este informe no declara superioridad entre herramientas: sin gold "
        "independiente y ejecuciones pareadas solo se describen diferencias.",
        "",
        "## Entorno",
        "",
    ]
    tool_rows = []
    for name, item in (env.get("tools") or {}).items():
        tool_rows.append([
            item.get("label", name),
            f"`{item.get('binary')}`",
            "sí" if item.get("available") else "no",
            f"`{item.get('path')}`" if item.get("path") else item.get("reason", "—"),
        ])
    lines += [
        _table(["herramienta", "binario", "disponible", "PATH / razón"], tool_rows),
        "",
        f"- Herramientas disponibles: **{', '.join(env.get('available_tools') or []) or 'ninguna'}**.",
        f"- Herramientas ausentes: **{', '.join(env.get('missing_tools') or []) or 'ninguna'}**.",
        f"- `--run` solicitado: **{_yes_no(params.get('run'))}**.",
        f"- Umbral de solapamiento recíproco: **{_fmt_num(params.get('overlap_threshold'), 2)}**.",
        "",
        "## Referencia",
        "",
    ]
    if reference.get("status") == "ok":
        counts = reference.get("counts") or {}
        hyp = reference.get("hypothetical") or {}
        lines += [
            f"- GFF: `{reference.get('path')}`.",
            f"- Conteos: CDS **{_fmt_int(counts.get('CDS'))}**, tRNA **{_fmt_int(counts.get('tRNA'))}**, "
            f"rRNA **{_fmt_int(counts.get('rRNA'))}**.",
            f"- Contigs: **{_fmt_int(reference.get('contigs'))}**; "
            f"CDS pseudo excluidos: **{_fmt_int(reference.get('excluded_pseudo'))}**.",
            f"- CDS con product: **{_fmt_int(hyp.get('with_product'))}**; "
            f"hipotéticas: **{_fmt_int(hyp.get('hypothetical'))}** "
            f"({_fmt_pct(hyp.get('pct_hypothetical'))}).",
        ]
    else:
        lines.append(f"- **No evaluable:** {reference.get('reason', '—')}.")
    lines += ["", "## Resultados", ""]
    result_rows = []
    for name, result in results.items():
        comp = result.get("comparison") or {}
        agg = comp.get("aggregate") or {}
        result_rows.append([
            result.get("label", name),
            "sí" if result.get("available") else "no",
            _status_label((result.get("run") or {}).get("status")) if result.get("available") else "no",
            _status_label(result.get("status")),
            _fmt_ratio(agg.get("precision")),
            _fmt_ratio(agg.get("recall")),
            _fmt_ratio(agg.get("f1")),
            _fmt_pct((agg.get("hypothetical_query") or {}).get("pct_hypothetical")),
        ])
    lines += [
        _table(
            ["herramienta", "binario", "ejecución", "comparación", "precisión",
             "recall", "F1", "% hypothetical (query)"],
            result_rows,
        ),
        "",
        "La precisión/recall/F1 se calculan contra el GFF de referencia, no contra "
        "una verdad funcional. No son una nota global del anotador.",
        "",
    ]

    for name, result in results.items():
        lines += [f"### {result.get('label', name)}", ""]
        lines.append(f"- Estado: **{_status_label(result.get('status'))}**.")
        if result.get("reason"):
            lines.append(f"- Razón: {result.get('reason')}.")
        if result.get("run"):
            run = result["run"]
            lines.append(f"- Ejecución: **{_run_line(run)}**.")
            if run.get("command"):
                lines.append(f"- Comando: `{' '.join(run['command'])}`.")
            if run.get("gff"):
                lines.append(f"- GFF temporal detectado: `{run['gff']}`.")
            if run.get("stderr_tail"):
                lines.append(f"- stderr (cola): `{run['stderr_tail'].strip()[-500:]}`")
        comp = result.get("comparison")
        if comp:
            lines.append(f"- Comparación: **{_status_label(comp.get('status'))}** "
                         f"contra `{comp.get('reference_gff')}`.")
            lines.append("")
            lines.append("#### Features")
            lines.append("")
            lines.append(_table(
                ["tipo", "ref", "query", "TP", "FP", "FN", "precision", "recall",
                 "F1", "exactos", "solape medio"],
                _category_rows(comp),
            ))
            lines.append("")
            lines.append("#### Atributos en pares emparejados")
            lines.append("")
            lines.append(_table(
                ["tipo", "pares", "product en ambos", "product literal",
                 "product tokens", "gene literal", "locus literal"],
                _attributes_rows(comp),
            ))
            lines.append("")
            lines.append("#### Hipotéticas")
            lines.append("")
            lines.append(_table(
                ["tipo", "ref CDS", "ref con product", "ref hyp", "% ref",
                 "query CDS", "query con product", "query hyp", "% query"],
                _hypothetical_rows(comp),
            ))
            lines.append("")
            for item in comp.get("not_evaluable") or []:
                lines.append(f"- **No evaluable:** {item.get('metric')}: {item.get('reason')}")
            lines.append("")
        else:
            lines.append("- Comparación: **no evaluable** (ver razones abajo).")
            lines.append("")

    lines += ["## No evaluable", ""]
    items = data.get("not_evaluable") or []
    if items:
        for item in items:
            tool = item.get("tool", "—")
            lines.append(f"- **{tool} / {item.get('metric', '?')}**: {item.get('reason', '')}")
    else:
        lines.append("- Ninguna métrica quedó no evaluable.")
    lines += ["", "## Limitaciones", ""]
    for item in data.get("limitations") or []:
        lines.append(f"- {item}")
    lines += [
        "",
        "## Reproducción",
        "",
        "```bash",
        "python tools/benchmark_external.py \\",
        f"  --genome {params.get('genome') or '<genoma.fna>'} \\",
        f"  --reference-gff {params.get('reference_gff') or '<referencia.gff.gz>'} \\",
        f"  --tool {params.get('tool', 'all')} \\",
        "  --md reports/external_benchmark.md \\",
        "  --json reports/external_benchmark.json",
        "```",
        "",
        "Para ejecutar de verdad una herramienta instalada, añadir `--run`. "
        "La ejecución usa un directorio temporal; el script no descarga bases "
        "ni consulta servicios web.",
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark externo contra bakta/prokka/dfast con datos locales.",
    )
    parser.add_argument("--genome", type=Path, default=None,
                        help="FASTA del genoma (solo necesario con --run)")
    parser.add_argument("--reference-gff", type=Path, default=None,
                        help="GFF3 de referencia (gold local) para comparar")
    parser.add_argument("--tool", choices=("all", *ALL_TOOLS), default="all")
    parser.add_argument("--run", action="store_true",
                        help="ejecuta la herramienta si el binario existe (sin red)")
    parser.add_argument("--bakta-db", type=Path, default=None,
                        help="ruta a la base Bakta local (opcional)")
    parser.add_argument("--overlap-threshold", type=float, default=0.5)
    parser.add_argument("--timeout", type=float, default=1800.0,
                        help="timeout por herramienta en segundos")
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    data = build_data(
        genome=args.genome,
        reference_gff=args.reference_gff,
        tool=args.tool,
        run=args.run,
        threshold=args.overlap_threshold,
        timeout=args.timeout,
        bakta_db=args.bakta_db,
    )
    raw_args = list(argv) if argv is not None else sys.argv[1:]
    data["command"] = [sys.executable, "tools/benchmark_external.py", *raw_args]
    md = render_markdown(data)
    args.md.parent.mkdir(parents=True, exist_ok=True)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.md.write_text(md, encoding="utf-8")
    args.json.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not args.quiet:
        print(f"[external] herramientas disponibles: {data['environment']['available_tools'] or 'ninguna'}")
        for name, result in data["results"].items():
            print(f"[external] {name}: {result.get('status')}")
        print(f"[external] escrito {args.md}")
        print(f"[external] escrito {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
