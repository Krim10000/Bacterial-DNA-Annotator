#!/usr/bin/env python3
"""Validador offline del golden set TSS-sigma.

No usa red. Lee ``schema.json`` y los TSV de ``data/sigma/gold/`` y comprueba:

- cabecera y columnas contra el esquema;
- tipos, enums y patrones;
- coordenadas 1-based mayores que cero;
- hebra ``+``/``-``;
- sigma permitido;
- coherencia con el genoma declarado (especie, accession, contig, fuente);
- duplicados por clave primaria;
- formato de PMIDs;
- resumen por genoma, sigma, fuente y confianza.

Uso::

    python scripts/sigma/40_validate_gold.py
    python scripts/sigma/40_validate_gold.py --json
    python scripts/sigma/40_validate_gold.py --gold-dir data/sigma/gold

El validador no reconstruye el gold ni consulta las fuentes primarias. Para
eso estan los tests de procedencia en ``tests/test_sigma_gold.py``.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GOLD_DIR = ROOT / "data" / "sigma" / "gold"
DEFAULT_SCHEMA = DEFAULT_GOLD_DIR / "schema.json"
PMID_SPLIT = re.compile(r"[;,]\s*")


# ---------------------------------------------------------------------------
# Esquema
# ---------------------------------------------------------------------------
def load_schema(path: Path | str) -> dict[str, Any]:
    """Carga y valida estructuralmente el esquema JSON."""
    path = Path(path)
    with path.open(encoding="utf-8") as fh:
        schema = json.load(fh)
    errors = validate_schema(schema)
    if errors:
        raise ValueError("; ".join(errors))
    return schema


def validate_schema(schema: Any) -> list[str]:
    """Comprueba que el esquema tenga las piezas que el validador necesita."""
    errors: list[str] = []
    if not isinstance(schema, dict):
        return ["schema: debe ser un objeto JSON"]

    for key in ("columns", "primary_key", "sigma_allowlist", "strand_values",
                "confidence_values", "genomes"):
        if key not in schema:
            errors.append(f"schema: falta la clave {key!r}")

    columns = schema.get("columns")
    if not isinstance(columns, list) or not columns:
        errors.append("schema: 'columns' debe ser una lista no vacia")
        columns = []

    names: list[str] = []
    by_name: dict[str, dict[str, Any]] = {}
    for i, col in enumerate(columns):
        if not isinstance(col, dict):
            errors.append(f"schema: columns[{i}] no es un objeto")
            continue
        name = col.get("name")
        if not name or not isinstance(name, str):
            errors.append(f"schema: columns[{i}] no tiene 'name' valido")
            continue
        if name in by_name:
            errors.append(f"schema: columna duplicada {name!r}")
        names.append(name)
        by_name[name] = col
        if not isinstance(col.get("nullable"), bool):
            errors.append(f"schema: {name}.nullable debe ser booleano")
        if "enum" in col and not isinstance(col["enum"], list):
            errors.append(f"schema: {name}.enum debe ser una lista")
        if "pattern" in col:
            try:
                re.compile(col["pattern"])
            except (re.error, TypeError) as exc:
                errors.append(f"schema: {name}.pattern invalido: {exc}")

    primary_key = schema.get("primary_key")
    if not isinstance(primary_key, list) or not primary_key:
        errors.append("schema: 'primary_key' debe ser una lista no vacia")
    else:
        if len(set(primary_key)) != len(primary_key):
            errors.append("schema: 'primary_key' tiene columnas repetidas")
        for key in primary_key:
            if key not in by_name:
                errors.append(f"schema: primary_key {key!r} no es una columna")

    for list_name in ("sigma_allowlist", "strand_values", "confidence_values"):
        value = schema.get(list_name)
        if not isinstance(value, list) or not value:
            errors.append(f"schema: {list_name} debe ser una lista no vacia")

    # El enum de las columnas clave no deberia separarse de las listas globales.
    for col_name, list_name in (
        ("sigma", "sigma_allowlist"),
        ("strand", "strand_values"),
        ("confidence", "confidence_values"),
    ):
        col = by_name.get(col_name)
        if col and isinstance(col.get("enum"), list):
            expected = schema.get(list_name)
            if isinstance(expected, list) and set(col["enum"]) != set(expected):
                errors.append(
                    f"schema: enum de {col_name} no coincide con {list_name}"
                )

    genomes = schema.get("genomes")
    if not isinstance(genomes, dict) or not genomes:
        errors.append("schema: 'genomes' debe ser un objeto no vacio")
    else:
        for genome_id, meta in genomes.items():
            if not isinstance(meta, dict):
                errors.append(f"schema: genomes.{genome_id} no es un objeto")
                continue
            for field in ("file", "species", "accession", "contig", "source"):
                if not meta.get(field):
                    errors.append(f"schema: genomes.{genome_id}.{field} falta")

    return errors


# ---------------------------------------------------------------------------
# Lectura y validacion
# ---------------------------------------------------------------------------
def read_tsv(path: Path | str) -> tuple[list[str], list[tuple[int, list[str]]]]:
    """Devuelve cabecera y filas con su numero de linea (1-based, cabecera=1)."""
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return [], []
        rows: list[tuple[int, list[str]]] = []
        for line, row in enumerate(reader, start=2):
            rows.append((line, row))
        return header, rows


def check_column(name: str, col: dict[str, Any], value: str) -> list[str]:
    """Valida una celda contra la definicion de su columna."""
    errors: list[str] = []
    nullable = bool(col.get("nullable", False))
    if value == "":
        if not nullable:
            errors.append(f"{name}: campo obligatorio vacio")
        return errors

    col_type = col.get("type", "string")
    if col_type == "integer":
        try:
            number = int(value)
        except (TypeError, ValueError):
            errors.append(f"{name}: se esperaba un entero, no {value!r}")
        else:
            minimum = col.get("minimum")
            if minimum is not None and number < minimum:
                errors.append(f"{name}: {number} es menor que el minimo {minimum}")

    enum = col.get("enum")
    if isinstance(enum, list) and value not in enum:
        errors.append(
            f"{name}: valor {value!r} no permitido "
            f"({', '.join(str(item) for item in enum)})"
        )

    pattern = col.get("pattern")
    if pattern and not re.fullmatch(pattern, value):
        errors.append(f"{name}: {value!r} no cumple el patron {pattern}")
    return errors


def check_pmid(value: str) -> tuple[list[str], list[str]]:
    """Valida PMIDs separados por coma o punto y coma."""
    errors: list[str] = []
    warnings: list[str] = []
    if not value:
        return errors, warnings
    parts = [part.strip() for part in PMID_SPLIT.split(value) if part.strip()]
    if not parts:
        errors.append("pmid: no contiene ningun identificador")
        return errors, warnings
    for part in parts:
        if not part.isdigit():
            errors.append(f"pmid: {part!r} no es numerico")
        elif int(part) <= 0:
            errors.append(f"pmid: {part!r} no es un PMID valido")
    if len(parts) != len(set(parts)):
        warnings.append("pmid: hay identificadores repetidos dentro de la fila")
    return errors, warnings


def check_genome(values: dict[str, str], genome_key: str,
                 meta: dict[str, Any]) -> list[str]:
    """Comprueba que la fila pertenezca al genoma que declara el esquema."""
    errors: list[str] = []
    if values.get("genome_id", "") != genome_key:
        errors.append(
            f"genome_id: {values.get('genome_id', '')!r} != {genome_key!r} "
            "para este archivo"
        )
    for field in ("species", "accession", "contig", "source"):
        expected = meta.get(field)
        observed = values.get(field, "")
        if expected and observed != expected:
            errors.append(
                f"{field}: {observed!r} != {expected!r} para {genome_key}"
            )
    return errors


def validate_file(path: Path | str, schema: dict[str, Any],
                  genome_key: str) -> dict[str, Any]:
    """Valida un TSV completo y devuelve errores, avisos y contadores."""
    path = Path(path)
    result: dict[str, Any] = {
        "file": str(path),
        "genome": genome_key,
        "rows": 0,
        "duplicates": 0,
        "errors": [],
        "warnings": [],
        "by_sigma": Counter(),
        "by_confidence": Counter(),
        "by_source": Counter(),
        "with_pmid": 0,
        "with_gene": 0,
        "keys": {},
    }

    if not path.is_file():
        result["errors"].append(f"{path.name}: no existe")
        return result

    header, raw_rows = read_tsv(path)
    if not header:
        result["errors"].append(f"{path.name}: archivo vacio o sin cabecera")
        return result
    if len(set(header)) != len(header):
        result["errors"].append(f"{path.name}: cabecera con nombres duplicados")

    expected_header = [col["name"] for col in schema["columns"]]
    if header != expected_header:
        missing = [name for name in expected_header if name not in header]
        extra = [name for name in header if name not in expected_header]
        detail = []
        if missing:
            detail.append(f"faltan={missing}")
        if extra:
            detail.append(f"sobran={extra}")
        if not missing and not extra:
            detail.append("orden distinto")
        result["errors"].append(
            f"{path.name}: cabecera no coincide con schema.json ({'; '.join(detail)})"
        )
        return result

    try:
        meta = schema["genomes"][genome_key]
    except KeyError:
        result["errors"].append(f"schema: no define el genoma {genome_key!r}")
        return result

    by_name = {col["name"]: col for col in schema["columns"]}
    seen: dict[tuple[str, ...], int] = {}

    for line, row in raw_rows:
        if not row or (len(row) == 1 and row[0] == ""):
            result["errors"].append(f"{path.name}:{line}: linea vacia")
            continue
        if len(row) != len(header):
            result["errors"].append(
                f"{path.name}:{line}: {len(row)} columnas, "
                f"se esperaban {len(header)}"
            )
            continue

        values = {
            name: (raw.strip() if isinstance(raw, str) else raw)
            for name, raw in zip(header, row)
        }
        row_errors: list[str] = []
        row_warnings: list[str] = []

        for name, col in by_name.items():
            value = values.get(name, "")
            if name == "pmid":
                if value == "":
                    if not col.get("nullable", False):
                        row_errors.append("pmid: campo obligatorio vacio")
                else:
                    pmid_errors, pmid_warnings = check_pmid(value)
                    row_errors.extend(pmid_errors)
                    row_warnings.extend(pmid_warnings)
                continue
            row_errors.extend(check_column(name, col, value))

        row_errors.extend(check_genome(values, genome_key, meta))

        for error in row_errors:
            result["errors"].append(f"{path.name}:{line}: {error}")
        for warning in row_warnings:
            result["warnings"].append(f"{path.name}:{line}: {warning}")

        key = tuple(values.get(name, "") for name in schema["primary_key"])
        if key in seen:
            result["duplicates"] += 1
            result["errors"].append(
                f"{path.name}:{line}: clave duplicada {key} "
                f"(primera aparicion en linea {seen[key]})"
            )
        else:
            seen[key] = line
            result["keys"][key] = line

        result["rows"] += 1
        result["by_sigma"][values.get("sigma", "")] += 1
        result["by_confidence"][values.get("confidence", "")] += 1
        result["by_source"][values.get("source", "")] += 1
        if values.get("pmid", ""):
            result["with_pmid"] += 1
        if values.get("gene", ""):
            result["with_gene"] += 1

    return result


def validate_gold(schema: dict[str, Any], gold_dir: Path | str) -> dict[str, Any]:
    """Valida todos los archivos declarados en ``schema['genomes']``."""
    gold_dir = Path(gold_dir)
    schema_errors = validate_schema(schema)
    if schema_errors:
        return {
            "ok": False,
            "schema_path": None,
            "schema_version": schema.get("version") if isinstance(schema, dict) else None,
            "gold_dir": str(gold_dir),
            "files": [],
            "summary": {},
            "errors": [f"schema: {error}" for error in schema_errors],
            "warnings": [],
        }

    errors: list[str] = []
    warnings: list[str] = []
    file_results: list[dict[str, Any]] = []
    file_summaries: list[dict[str, Any]] = []
    global_keys: dict[tuple[str, ...], tuple[str, int]] = {}

    for genome_key, meta in schema["genomes"].items():
        result = validate_file(gold_dir / meta["file"], schema, genome_key)
        file_results.append(result)
        errors.extend(result["errors"])
        warnings.extend(result["warnings"])

        for key, line in result["keys"].items():
            if key in global_keys:
                other_file, other_line = global_keys[key]
                errors.append(
                    f"clave duplicada entre archivos: {key} "
                    f"({Path(result['file']).name}:{line} y "
                    f"{Path(other_file).name}:{other_line})"
                )
            else:
                global_keys[key] = (result["file"], line)

        file_summaries.append({
            "file": result["file"],
            "genome": result["genome"],
            "rows": result["rows"],
            "duplicates": result["duplicates"],
            "errors": list(result["errors"]),
            "warnings": list(result["warnings"]),
            "by_sigma": dict(sorted(result["by_sigma"].items())),
            "by_confidence": dict(sorted(result["by_confidence"].items())),
            "by_source": dict(sorted(result["by_source"].items())),
            "with_pmid": result["with_pmid"],
            "with_gene": result["with_gene"],
        })

    by_genome: Counter[str] = Counter()
    by_sigma: Counter[str] = Counter()
    by_confidence: Counter[str] = Counter()
    by_source: Counter[str] = Counter()
    with_pmid = 0
    with_gene = 0
    duplicates = 0
    for result in file_results:
        by_genome[result["genome"]] += result["rows"]
        by_sigma.update(result["by_sigma"])
        by_confidence.update(result["by_confidence"])
        by_source.update(result["by_source"])
        with_pmid += result["with_pmid"]
        with_gene += result["with_gene"]
        duplicates += result["duplicates"]

    total_rows = sum(by_genome.values())
    summary = {
        "rows": total_rows,
        "files": len(file_results),
        "genomes": dict(sorted(by_genome.items())),
        "sigma": dict(sorted(by_sigma.items())),
        "confidence": dict(sorted(by_confidence.items())),
        "source": dict(sorted(by_source.items())),
        "pmid": {"with": with_pmid, "total": total_rows},
        "gene": {"with": with_gene, "total": total_rows},
        "duplicates": duplicates,
    }

    return {
        "ok": not errors,
        "schema_path": None,
        "schema_version": schema.get("version"),
        "gold_dir": str(gold_dir),
        "files": file_summaries,
        "summary": summary,
        "errors": errors,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Salida
# ---------------------------------------------------------------------------
def _relative(path: str) -> str:
    try:
        return str(Path(path).relative_to(ROOT))
    except ValueError:
        return str(path)


def _counts_text(label: str, values: dict[str, int]) -> str:
    if not values:
        return f"  {label}: -"
    text = ", ".join(f"{key or '<vacio>'}={value}" for key, value in values.items())
    return f"  {label}: {text}"


def format_text(report: dict[str, Any], max_errors: int = 50) -> str:
    """Formatea el reporte para humanos."""
    lines: list[str] = []
    lines.append("Validador sigma gold (offline)")
    schema_path = report.get("schema_path")
    if schema_path:
        lines.append(f"schema: {_relative(str(schema_path))} "
                     f"(v{report.get('schema_version', '?')})")
    lines.append(f"gold:   {_relative(str(report.get('gold_dir', '')))}")
    lines.append("")

    for file in report.get("files", []):
        status = "OK" if not file["errors"] else "FAIL"
        suffix = f" ({len(file['errors'])} errores)" if file["errors"] else ""
        lines.append(
            f"[{status}] {file['genome']}: "
            f"{_relative(file['file'])} -> {file['rows']} filas{suffix}"
        )

    summary = report.get("summary") or {}
    lines.append("")
    lines.append("Resumen")
    lines.append(f"  filas: {summary.get('rows', 0)}")
    lines.append(_counts_text("por genoma", summary.get("genomes", {})))
    lines.append(_counts_text("por sigma", summary.get("sigma", {})))
    lines.append(_counts_text("por fuente", summary.get("source", {})))
    lines.append(_counts_text("por confianza", summary.get("confidence", {})))

    pmid = summary.get("pmid", {})
    total = pmid.get("total", 0)
    with_pmid = pmid.get("with", 0)
    pct = (100.0 * with_pmid / total) if total else 0.0
    lines.append(f"  PMID: {with_pmid}/{total} filas ({pct:.1f}%)")
    gene = summary.get("gene", {})
    total_gene = gene.get("total", 0)
    with_gene = gene.get("with", 0)
    gene_pct = (100.0 * with_gene / total_gene) if total_gene else 0.0
    lines.append(f"  gen: {with_gene}/{total_gene} filas ({gene_pct:.1f}%)")
    lines.append(f"  duplicados: {summary.get('duplicates', 0)}")
    lines.append(f"  errores: {len(report.get('errors', []))}")
    lines.append(f"  avisos: {len(report.get('warnings', []))}")

    errors = report.get("errors", [])
    warnings = report.get("warnings", [])
    if errors:
        lines.append("")
        lines.append(f"Errores ({len(errors)})")
        for error in errors[:max_errors]:
            lines.append(f"  - {error}")
        if len(errors) > max_errors:
            lines.append(f"  ... ({len(errors) - max_errors} mas)")
    if warnings:
        lines.append("")
        lines.append(f"Avisos ({len(warnings)})")
        for warning in warnings[:max_errors]:
            lines.append(f"  - {warning}")
        if len(warnings) > max_errors:
            lines.append(f"  ... ({len(warnings) - max_errors} mas)")

    lines.append("")
    if report.get("ok"):
        lines.append("OK: gold valido")
    else:
        lines.append("FALLO: gold invalido")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Valida data/sigma/gold/ contra schema.json (sin red)."
    )
    parser.add_argument(
        "--gold-dir",
        type=Path,
        default=DEFAULT_GOLD_DIR,
        help="Directorio del golden set (default: data/sigma/gold).",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=None,
        help="Ruta al schema.json (default: <gold-dir>/schema.json).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="Emite el reporte completo en JSON.",
    )
    parser.add_argument(
        "--max-errors",
        type=int,
        default=50,
        help="Maximo de errores/avisos a imprimir en modo texto (default: 50).",
    )
    args = parser.parse_args(argv)

    gold_dir = args.gold_dir
    schema_path = args.schema or (gold_dir / "schema.json")
    try:
        schema = load_schema(schema_path)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"ERROR: no se pudo cargar el schema {schema_path}: {exc}",
              file=sys.stderr)
        return 2

    report = validate_gold(schema, gold_dir)
    report["schema_path"] = str(schema_path)

    if args.as_json:
        print(json.dumps(report, ensure_ascii=False, indent=2, default=str))
    else:
        print(format_text(report, max_errors=max(0, args.max_errors)))

    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
