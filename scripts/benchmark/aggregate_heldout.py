#!/usr/bin/env python3
"""Agrega los JSON de un benchmark held-out por genoma y herramienta.

Lee los documentos producidos por ``benchmarks/external/external_benchmark.py``
(esquema ``external-benchmark-run/v1``) y escribe un resumen agregado en JSON y
Markdown. No ejecuta herramientas ni modifica los JSON de entrada.

Métricas agregadas por herramienta:

- TP/FP/FN sumados (micro) y precisión/recall/F1 micro;
- F1 medio, mínimo y máximo por genoma (macro de genomas);
- coincidencias exactas de bordes y pares con algún borde distinto, con su
  tasa sobre los TP;
- tiempo de pared sumado, medio y máximo;
- agregado micro por categoría (CDS/tRNA/rRNA).

El resumen no declara superioridad. Si la diferencia de F1 entre herramientas
es <= 1 punto porcentual, la etiqueta explícita es "sin ganador".
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA_INPUT = "external-benchmark-run/v1"
SCHEMA_OUTPUT = "external-benchmark-heldout-aggregate/v1"
TOOLS = ("bakta", "prokka", "dfast")
CATEGORIES = ("CDS", "tRNA", "rRNA")
DEFAULT_INPUT_DIR = Path("reports/external_benchmark/heldout")
DEFAULT_TUNING_GENOMES = ("cglutamicum", "mtb_h37rv", "paeruginosa", "pputida")
SMALL_DIFFERENCE_F1 = 0.01


class AggregateError(RuntimeError):
    """Error de datos que impide construir un resumen fiable."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _safe_div(numerator: float | int | None, denominator: float | int | None) -> float | None:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def _ratio(numerator: float | int | None, denominator: float | int | None) -> float | None:
    value = _safe_div(numerator, denominator)
    return None if value is None else round(value, 6)


def _f1(tp: float | int | None, fp: float | int | None, fn: float | int | None) -> float | None:
    if tp is None or fp is None or fn is None:
        return None
    denominator = 2 * tp + fp + fn
    if denominator == 0:
        return None
    return round(2 * tp / denominator, 6)


def _number(value: Any) -> float | int | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def _round(value: float | None, digits: int = 6) -> float | None:
    return None if value is None else round(value, digits)


def _fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _pct(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "—"
    return f"{100 * value:.{digits}f}%"


def load_document(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AggregateError(f"no se pudo leer {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise AggregateError(f"{path}: el documento raíz no es un objeto JSON")
    if document.get("schema") != SCHEMA_INPUT:
        raise AggregateError(
            f"{path}: schema inesperado {document.get('schema')!r}; "
            f"se esperaba {SCHEMA_INPUT!r}"
        )
    return document


def assembly_metadata(genome_id: str, repo_root: Path) -> dict[str, str]:
    path = repo_root / "data" / "reference_genomes" / genome_id / "assembly.json"
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(data, dict):
        return {}
    return {
        key: str(data[key])
        for key in ("accession", "assembly_name", "species", "status")
        if key in data
    }


def image_id_from_log(log_path: Path) -> str | None:
    if not log_path.is_file():
        return None
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    matches = re.findall(r"image_id=(sha256:[0-9a-f]{64})", text)
    return matches[-1] if matches else None


def summarize_metrics(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    for row in ok_rows:
        for key in ("tp", "fp", "fn"):
            if row.get(key) is None:
                raise AggregateError(
                    f"{row['genome']}/{row['tool']}: status ok sin {key} en aggregate"
                )

    tp = sum(int(row["tp"]) for row in ok_rows)
    fp = sum(int(row["fp"]) for row in ok_rows)
    fn = sum(int(row["fn"]) for row in ok_rows)
    exact = sum(int(row.get("exact_matches") or 0) for row in ok_rows)
    shifted = sum(int(row.get("boundary_shifted") or 0) for row in ok_rows)
    runtime_values = [
        float(row["runtime_seconds"])
        for row in ok_rows
        if row.get("runtime_seconds") is not None
    ]
    f1_values = [float(row["f1"]) for row in ok_rows if row.get("f1") is not None]
    precision_values = [
        float(row["precision"]) for row in ok_rows if row.get("precision") is not None
    ]
    recall_values = [
        float(row["recall"]) for row in ok_rows if row.get("recall") is not None
    ]

    per_genome = []
    for row in rows:
        per_genome.append(
            {
                "genome": row["genome"],
                "status": row.get("status"),
                "reason": row.get("reason") or None,
                "source": row.get("source"),
                "tp": row.get("tp"),
                "fp": row.get("fp"),
                "fn": row.get("fn"),
                "precision": row.get("precision"),
                "recall": row.get("recall"),
                "f1": row.get("f1"),
                "exact_matches": row.get("exact_matches"),
                "boundary_shifted": row.get("boundary_shifted"),
                "exact_match_rate": _ratio(row.get("exact_matches"), row.get("tp")),
                "boundary_shift_rate": _ratio(row.get("boundary_shifted"), row.get("tp")),
                "runtime_seconds": row.get("runtime_seconds"),
                "version": row.get("version"),
            }
        )

    return {
        "genomes_total": len(rows),
        "genomes_ok": len(ok_rows),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": _ratio(tp, tp + fp),
        "recall": _ratio(tp, tp + fn),
        "f1": _f1(tp, fp, fn),
        "mean_genome_f1": _round(statistics.fmean(f1_values)) if f1_values else None,
        "min_genome_f1": _round(min(f1_values)) if f1_values else None,
        "max_genome_f1": _round(max(f1_values)) if f1_values else None,
        "mean_genome_precision": _round(statistics.fmean(precision_values))
        if precision_values
        else None,
        "mean_genome_recall": _round(statistics.fmean(recall_values))
        if recall_values
        else None,
        "exact_matches": exact,
        "exact_match_rate": _ratio(exact, tp),
        "boundary_shifted": shifted,
        "boundary_shift_rate": _ratio(shifted, tp),
        "total_runtime_seconds": _round(sum(runtime_values)) if runtime_values else None,
        "mean_runtime_seconds": _round(statistics.fmean(runtime_values))
        if runtime_values
        else None,
        "max_runtime_seconds": _round(max(runtime_values)) if runtime_values else None,
        "per_genome": per_genome,
    }


def summarize_categories(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    for category in CATEGORIES:
        category_rows = [
            row["categories"][category]
            for row in ok_rows
            if isinstance(row.get("categories"), Mapping)
            and isinstance(row["categories"].get(category), Mapping)
        ]
        tp = sum(int(item.get("tp") or 0) for item in category_rows)
        fp = sum(int(item.get("fp") or 0) for item in category_rows)
        fn = sum(int(item.get("fn") or 0) for item in category_rows)
        summary[category] = {
            "genomes_with_category": len(category_rows),
            "reference_n": sum(int(item.get("reference_n") or 0) for item in category_rows),
            "query_n": sum(int(item.get("query_n") or 0) for item in category_rows),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": _ratio(tp, tp + fp),
            "recall": _ratio(tp, tp + fn),
            "f1": _f1(tp, fp, fn),
            "exact_matches": sum(int(item.get("exact_matches") or 0) for item in category_rows),
            "boundary_shifted": sum(
                int(item.get("boundary_shifted") or 0) for item in category_rows
            ),
        }
    return summary


def _unique(items: Iterable[str]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for item in items:
        text = str(item).strip()
        if text and text not in seen:
            seen.add(text)
            output.append(text)
    return output


def summarize_excluded(
    documents: Sequence[tuple[Path, Mapping[str, Any]]],
    *,
    reason: str,
) -> list[dict[str, Any]]:
    excluded: list[dict[str, Any]] = []
    for path, document in documents:
        genome_id = path.stem
        results = document.get("results") or {}
        tools: dict[str, Any] = {}
        for tool in TOOLS:
            result = results.get(tool) or {}
            aggregate = (result.get("comparison") or {}).get("aggregate") or {}
            tools[tool] = {
                "status": result.get("status"),
                "reason": result.get("reason") or None,
                "f1": _number(aggregate.get("f1")),
                "tp": _number(aggregate.get("tp")),
                "fp": _number(aggregate.get("fp")),
                "fn": _number(aggregate.get("fn")),
            }
        excluded.append(
            {
                "genome": genome_id,
                "json": str(path),
                "sha256": sha256_file(path),
                "reason": reason,
                "tools": tools,
            }
        )
    return excluded


def build_summary(
    documents: Sequence[tuple[Path, Mapping[str, Any]]],
    *,
    input_dir: Path,
    repo_root: Path,
    tuning_genomes: Sequence[str],
    log_path: Path | None,
    excluded_documents: Sequence[tuple[Path, Mapping[str, Any]]] = (),
    exclude_reason: str = "",
) -> dict[str, Any]:
    if not documents:
        raise AggregateError(f"no se encontraron JSON de benchmark en {input_dir}")

    genomes: list[dict[str, Any]] = []
    rows_by_tool: dict[str, list[dict[str, Any]]] = {tool: [] for tool in TOOLS}
    limitations: list[str] = []
    warnings: list[str] = []
    commands: dict[str, list[str]] = {}
    threads_seen: set[int] = set()
    timeouts_seen: set[float] = set()
    thresholds_seen: set[float] = set()
    tools_seen: set[str] = set()
    versions: dict[str, str] = {}

    for path, document in documents:
        genome_id = path.stem
        if not genome_id:
            raise AggregateError(f"{path}: no se pudo derivar el id de genoma")
        results = document.get("results")
        if not isinstance(results, Mapping):
            raise AggregateError(f"{path}: falta 'results'")
        parameters = document.get("parameters") or {}
        if parameters.get("run") is not True:
            raise AggregateError(f"{path}: el JSON no corresponde a una corrida --run")
        if parameters.get("threads") is not None:
            threads_seen.add(int(parameters["threads"]))
        if parameters.get("timeout_seconds") is not None:
            timeouts_seen.add(float(parameters["timeout_seconds"]))
        if parameters.get("overlap_threshold") is not None:
            thresholds_seen.add(float(parameters["overlap_threshold"]))
        selected_tools = parameters.get("tools")
        if isinstance(selected_tools, list):
            tools_seen.update(str(item) for item in selected_tools)

        inputs = document.get("inputs") or {}
        genome_input = inputs.get("genome") or {}
        reference_input = inputs.get("reference") or {}
        genome_entry = {
            "id": genome_id,
            "species": None,
            "accession": None,
            "assembly_name": None,
            "assembly_status": None,
            "fasta": genome_input.get("path"),
            "fasta_sha256": genome_input.get("sha256"),
            "total_length": genome_input.get("total_length"),
            "reference": reference_input.get("path"),
            "reference_sha256": reference_input.get("sha256"),
            "reference_counts": reference_input.get("counts"),
            "tools": {},
        }
        genome_entry.update(assembly_metadata(genome_id, repo_root))
        genomes.append(genome_entry)

        for tool in TOOLS:
            result = results.get(tool)
            if result is None:
                row = {
                    "genome": genome_id,
                    "tool": tool,
                    "status": "missing",
                    "reason": "la herramienta no aparece en el JSON",
                }
            else:
                comparison = result.get("comparison") or {}
                aggregate = comparison.get("aggregate") or {}
                execution = result.get("execution") or {}
                version_field = result.get("version")
                if isinstance(version_field, Mapping):
                    version_text = version_field.get("version") or version_field.get("text")
                else:
                    version_text = version_field
                row = {
                    "genome": genome_id,
                    "tool": tool,
                    "status": result.get("status"),
                    "reason": result.get("reason") or None,
                    "source": result.get("source"),
                    "version": str(version_text) if version_text else None,
                    "tp": _number(aggregate.get("tp")),
                    "fp": _number(aggregate.get("fp")),
                    "fn": _number(aggregate.get("fn")),
                    "precision": _number(aggregate.get("precision")),
                    "recall": _number(aggregate.get("recall")),
                    "f1": _number(aggregate.get("f1")),
                    "exact_matches": _number(aggregate.get("exact_matches")),
                    "boundary_shifted": _number(aggregate.get("boundary_shifted")),
                    "runtime_seconds": _number(execution.get("runtime_seconds")),
                    "categories": comparison.get("categories") or {},
                }
                if version_text:
                    versions[tool] = str(version_text)
                command = execution.get("command")
                if isinstance(command, list):
                    commands.setdefault(genome_id, []).append(" ".join(str(part) for part in command))
                for item in result.get("limitations") or comparison.get("limitations") or []:
                    limitations.append(str(item))
            rows_by_tool[tool].append(row)
            genome_entry["tools"][tool] = row["status"]

        for item in document.get("limitations") or []:
            limitations.append(str(item))
        for item in document.get("warnings") or []:
            warnings.append(f"{genome_id}: {item}")

    tools_summary: dict[str, Any] = {}
    for tool in TOOLS:
        rows = rows_by_tool[tool]
        metrics = summarize_metrics(rows)
        metrics["categories"] = summarize_categories(rows)
        metrics["label"] = {
            "bakta": "Bakta",
            "prokka": "Prokka",
            "dfast": "DFAST",
        }[tool]
        metrics["version"] = versions.get(tool)
        tools_summary[tool] = metrics

    incomplete = {
        tool: {
            "genomes_ok": data["genomes_ok"],
            "genomes_total": data["genomes_total"],
            "failed_or_missing": [
                f"{row['genome']}:{row['status']}"
                for row in data["per_genome"]
                if row["status"] != "ok"
            ],
        }
        for tool, data in tools_summary.items()
        if data["genomes_ok"] != data["genomes_total"]
    }

    f1_values = {
        tool: data["f1"] for tool, data in tools_summary.items() if data["f1"] is not None
    }
    if len(f1_values) >= 2:
        spread = max(f1_values.values()) - min(f1_values.values())
    else:
        spread = None
    if incomplete:
        verdict = "incomplete"
        verdict_text = (
            "No se puede comparar: al menos una herramienta no quedó 'ok' en "
            "todos los genomas."
        )
    elif spread is not None and spread <= SMALL_DIFFERENCE_F1:
        verdict = "no_winner_small_differences"
        verdict_text = (
            f"Sin ganador: la dispersión de F1 micro es {spread * 100:.2f} puntos "
            f"porcentuales (≤ {SMALL_DIFFERENCE_F1 * 100:.0f} pp). Las diferencias "
            "no permiten declarar superioridad."
        )
    else:
        winner = max(f1_values, key=f1_values.get) if f1_values else None
        verdict = "highest_observed_f1_no_superiority_claim"
        verdict_text = (
            f"Mayor F1 observado: {winner}. No se declara ganador: el diseño no "
            "incluye réplicas ni gold curado a mano y la referencia hereda sus errores."
        )

    overlap = sorted(set(genome_id for genome_id in (g["id"] for g in genomes)) & set(tuning_genomes))
    if overlap:
        warnings.append(
            "Genomas de la selección que también están en el set de tuning: "
            + ", ".join(overlap)
        )

    standard_limitations = [
        "La referencia es el gold de esta comparación; sus criterios y errores se heredan.",
        "El emparejamiento es 1:1 voraz por solapamiento recíproco al umbral indicado; "
        "no es un óptimo bipartito y otro umbral cambiaría TP/FP/FN.",
        "La coincidencia de 'product' no mide función; los nombres dependen de cada base.",
        "No se comparan ncRNA, tmRNA, pseudogenes (si se excluyen), dominios ni rutas.",
        "Los tiempos son de una sola ejecución secuencial por genoma en el mismo host; "
        "no hay réplicas ni intervalo de confianza.",
        "El F1 agregado es micro (TP/FP/FN sumados); el F1 medio es la media aritmética "
        "por genoma.",
        "El resumen no declara superioridad entre herramientas.",
    ]

    return {
        "schema": SCHEMA_OUTPUT,
        "generated_at": utc_now(),
        "scope": (
            "Resumen agregado de un benchmark held-out por genoma para Bakta, Prokka "
            "y DFAST contra anotaciones de referencia locales."
        ),
        "verdict": {
            "winner_declared": False,
            "label": verdict,
            "text": verdict_text,
            "f1_spread": _round(spread, 6),
            "f1_spread_percentage_points": _round(spread * 100, 4) if spread is not None else None,
            "small_difference_threshold": SMALL_DIFFERENCE_F1,
        },
        "selection": {
            "genomes": [genome["id"] for genome in genomes],
            "excluded_from_aggregate": [path.stem for path, _doc in excluded_documents],
            "tuning_genomes": list(tuning_genomes),
            "tuning_overlap": overlap,
            "heldout_note": (
                "Selección fijada antes de ejecutar; no incluye los genomas de tuning "
                "ni el gold de 60 kb."
            ),
        },
        "parameters": {
            "threads": sorted(threads_seen),
            "timeout_seconds": sorted(timeouts_seen),
            "overlap_threshold": sorted(thresholds_seen),
            "tools_selected": sorted(tools_seen),
        },
        "provenance": {
            "input_dir": str(input_dir),
            "input_files": [
                {
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
                for path, _document in [*documents, *excluded_documents]
            ],
            "run_log": str(log_path) if log_path and log_path.is_file() else None,
            "image_id": image_id_from_log(log_path) if log_path else None,
            "tool_versions": versions,
            "commands": commands,
        },
        "genomes": genomes,
        "excluded_genomes": summarize_excluded(excluded_documents, reason=exclude_reason),
        "tools": tools_summary,
        "incomplete_tools": incomplete,
        "warnings": _unique(warnings),
        "limitations": _unique([*limitations, *standard_limitations]),
    }


def render_markdown(summary: Mapping[str, Any]) -> str:
    tools = summary["tools"]
    lines: list[str] = []
    lines.append("# Benchmark held-out: Bakta, Prokka y DFAST")
    lines.append("")
    lines.append(f"**Generado:** {summary['generated_at']}  ")
    lines.append(f"**Alcance:** {summary['scope']}")
    lines.append("")
    selection = summary["selection"]
    lines.append(
        "**Genomas held-out:** " + ", ".join(f"`{g}`" for g in selection["genomes"]) + "  "
    )
    lines.append(
        "**Set de tuning excluido:** "
        + ", ".join(f"`{g}`" for g in selection["tuning_genomes"])
    )
    lines.append("")
    verdict = summary["verdict"]
    lines.append(f"**Veredicto:** {verdict['text']}")
    lines.append("")

    lines.append("## Resumen agregado por herramienta")
    lines.append("")
    lines.append(
        "| Herramienta | Versión | Genomas ok | TP | FP | FN | Precisión | Recall | "
        "F1 micro | F1 medio | Exactos | Tasa exacta | Bordes distintos | "
        "Tiempo total (s) | Tiempo medio (s) |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for tool in TOOLS:
        data = tools[tool]
        lines.append(
            "| {label} | {version} | {ok}/{total} | {tp} | {fp} | {fn} | {p} | {r} | "
            "{f1} | {mean_f1} | {exact} | {exact_rate} | {shifted} | {total_time} | "
            "{mean_time} |".format(
                label=data["label"],
                version=data["version"] or "—",
                ok=data["genomes_ok"],
                total=data["genomes_total"],
                tp=data["tp"],
                fp=data["fp"],
                fn=data["fn"],
                p=_fmt(data["precision"], 4),
                r=_fmt(data["recall"], 4),
                f1=_fmt(data["f1"], 4),
                mean_f1=_fmt(data["mean_genome_f1"], 4),
                exact=data["exact_matches"],
                exact_rate=_pct(data["exact_match_rate"]),
                shifted=data["boundary_shifted"],
                total_time=_fmt(data["total_runtime_seconds"], 1),
                mean_time=_fmt(data["mean_runtime_seconds"], 1),
            )
        )
    lines.append("")
    lines.append(
        "> TP/FP/FN y F1 micro se calculan sumando los conteos de todos los genomas. "
        "La tasa exacta es `exact_matches / TP`; un par puede ser TP con bordes "
        "distintos. El tiempo es la suma de `runtime_seconds` de cada herramienta."
    )
    lines.append("")

    lines.append("## Resultado por genoma")
    lines.append("")
    lines.append(
        "| Genoma | Herramienta | Estado | TP | FP | FN | P | R | F1 | Exactos | "
        "Bordes distintos | Tiempo (s) |"
    )
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for genome in summary["genomes"]:
        for tool in TOOLS:
            row = next(
                (item for item in tools[tool]["per_genome"] if item["genome"] == genome["id"]),
                None,
            )
            if row is None:
                continue
            lines.append(
                "| {genome} | {label} | `{status}` | {tp} | {fp} | {fn} | {p} | {r} | "
                "{f1} | {exact} | {shifted} | {time} |".format(
                    genome=genome["id"],
                    label=tools[tool]["label"],
                    status=row["status"],
                    tp=_fmt(row["tp"]),
                    fp=_fmt(row["fp"]),
                    fn=_fmt(row["fn"]),
                    p=_fmt(row["precision"], 4),
                    r=_fmt(row["recall"], 4),
                    f1=_fmt(row["f1"], 4),
                    exact=_fmt(row["exact_matches"]),
                    shifted=_fmt(row["boundary_shifted"]),
                    time=_fmt(row["runtime_seconds"], 1),
                )
            )
            if row["status"] != "ok" and row.get("reason"):
                lines.append(
                    f"| {genome['id']} | {tools[tool]['label']} | razón | — | — | — | — | — | — | — | — | — |"
                )
    lines.append("")

    lines.append("## Agregado por categoría")
    lines.append("")
    lines.append("| Herramienta | Categoría | TP | FP | FN | P | R | F1 | Exactos |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for tool in TOOLS:
        for category in CATEGORIES:
            data = tools[tool]["categories"][category]
            lines.append(
                "| {label} | {category} | {tp} | {fp} | {fn} | {p} | {r} | {f1} | {exact} |".format(
                    label=tools[tool]["label"],
                    category=category,
                    tp=data["tp"],
                    fp=data["fp"],
                    fn=data["fn"],
                    p=_fmt(data["precision"], 4),
                    r=_fmt(data["recall"], 4),
                    f1=_fmt(data["f1"], 4),
                    exact=data["exact_matches"],
                )
            )
    lines.append("")

    lines.append("## TP / FP / FN")
    lines.append("")
    for tool in TOOLS:
        data = tools[tool]
        lines.append(
            f"- **{data['label']}**: TP={data['tp']}, FP={data['fp']}, FN={data['fn']} "
            f"(precisión {_fmt(data['precision'], 4)}, recall {_fmt(data['recall'], 4)}, "
            f"F1 {_fmt(data['f1'], 4)})."
        )
    lines.append("")

    incomplete = summary.get("incomplete_tools") or {}
    if incomplete:
        lines.append("## Herramientas incompletas")
        lines.append("")
        for tool, data in incomplete.items():
            lines.append(
                f"- **{tools[tool]['label']}**: {data['genomes_ok']}/{data['genomes_total']} "
                f"genomas ok; no ok: {', '.join(data['failed_or_missing']) or '—'}."
            )
        lines.append("")

    excluded = summary.get("excluded_genomes") or []
    if excluded:
        lines.append("## Genomas excluidos del agregado")
        lines.append("")
        for item in excluded:
            statuses = ", ".join(
                f"{tools[name]['label']}={data['status']}" for name, data in item["tools"].items()
            )
            lines.append(
                f"- `{item['genome']}`: {item['reason']} Estado por herramienta: {statuses}."
            )
            for name, data in item["tools"].items():
                if data.get("status") == "ok" and data.get("f1") is not None:
                    lines.append(
                        f"  - {tools[name]['label']}: F1 {_fmt(data['f1'], 4)}, "
                        f"TP/FP/FN {data['tp']}/{data['fp']}/{data['fn']}."
                    )
        lines.append("")

    lines.append("## Entradas y procedencia")
    lines.append("")
    provenance = summary["provenance"]
    lines.append(f"- Imagen: `{provenance.get('image_id') or 'no registrada'}`.")
    lines.append(
        f"- Hilos: {summary['parameters']['threads'] or '—'}; "
        f"timeout: {summary['parameters']['timeout_seconds'] or '—'} s; "
        f"umbral de solapamiento: {summary['parameters']['overlap_threshold'] or '—'}."
    )
    for genome in summary["genomes"]:
        species = genome.get("species") or genome["id"]
        accession = genome.get("accession") or "—"
        lines.append(
            f"- `{genome['id']}` — {species} ({accession}); "
            f"FASTA `{genome['fasta_sha256'] or '—'}`; "
            f"referencia `{genome['reference_sha256'] or '—'}`."
        )
    lines.append("")

    lines.append("## Limitaciones")
    lines.append("")
    for item in summary["limitations"]:
        lines.append(f"- {item}")
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Agrega JSON de benchmark held-out (esquema external-benchmark-run/v1) "
            "en un resumen JSON/Markdown por herramienta."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"directorio con los JSON por genoma (por defecto: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        help="ruta del resumen JSON (por defecto: <input-dir>/summary.json)",
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=None,
        help="ruta del resumen Markdown (por defecto: <input-dir>/summary.md)",
    )
    parser.add_argument(
        "--run-log",
        type=Path,
        default=None,
        help="log opcional de la corrida para registrar la imagen Docker",
    )
    parser.add_argument(
        "--tuning-genomes",
        default=",".join(DEFAULT_TUNING_GENOMES),
        help="ids separados por coma que no deben aparecer en el set held-out",
    )
    parser.add_argument(
        "--exclude-genome",
        action="append",
        default=[],
        metavar="ID",
        help="id de genoma a excluir del agregado (se reporta en una sección aparte); repetible",
    )
    parser.add_argument(
        "--exclude-reason",
        default="excluido del agregado por decisión explícita",
        help="motivo que se registra para los genomas excluidos",
    )
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="no falla si alguna herramienta no está 'ok' en todos los genomas",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="no imprime el resumen por stdout",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    input_dir = args.input_dir
    if not input_dir.is_dir():
        print(f"error: no existe el directorio {input_dir}", file=sys.stderr)
        return 2

    output_json = args.json or input_dir / "summary.json"
    output_md = args.md or input_dir / "summary.md"
    excluded_names = {output_json.name, output_md.name, "summary.json", "summary.md"}
    paths = sorted(
        path
        for path in input_dir.glob("*.json")
        if path.is_file() and path.name not in excluded_names and not path.name.startswith(".")
    )
    if not paths:
        print(f"error: no hay JSON de benchmark en {input_dir}", file=sys.stderr)
        return 2

    repo_root = Path(__file__).resolve().parents[2]
    documents_all: list[tuple[Path, Mapping[str, Any]]] = []
    try:
        for path in paths:
            documents_all.append((path, load_document(path)))
        tuning = [item.strip() for item in args.tuning_genomes.split(",") if item.strip()]
        exclude = {item.strip() for item in args.exclude_genome if item.strip()}
        documents = [(path, doc) for path, doc in documents_all if path.stem not in exclude]
        excluded_documents = [
            (path, doc) for path, doc in documents_all if path.stem in exclude
        ]
        missing_exclusions = sorted(exclude - {path.stem for path, _doc in excluded_documents})
        if missing_exclusions:
            print(
                "error: --exclude-genome no encontrado en los JSON: "
                + ", ".join(missing_exclusions),
                file=sys.stderr,
            )
            return 2
        summary = build_summary(
            documents,
            input_dir=input_dir,
            repo_root=repo_root,
            tuning_genomes=tuning,
            log_path=args.run_log,
            excluded_documents=excluded_documents,
            exclude_reason=args.exclude_reason,
        )
    except AggregateError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    atomic_write_text(output_json, json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    atomic_write_text(output_md, render_markdown(summary))

    if not args.quiet:
        print(render_markdown(summary), end="")

    incomplete = summary.get("incomplete_tools") or {}
    if incomplete and not args.allow_incomplete:
        print(
            "error: hay herramientas sin status=ok en todos los genomas "
            f"({', '.join(incomplete)}); revisar el resumen",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
