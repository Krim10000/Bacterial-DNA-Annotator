#!/usr/bin/env python3
"""Paridad de conteos: anotación del sistema vs GFF de genomas de referencia.

Compara CDS/tRNA/rRNA del resumen del sistema (``results/<id>/<id>.summary.json``)
contra los features del GFF de referencia declarado en
``data/reference_genomes.json``.

Alcance de esta versión:

- solo conteos por clase (CDS, tRNA, rRNA);
- no evalúa función, productos, coordenadas ni identidad de secuencia;
- no modifica ``results/``, los GFF ni el manifiesto;
- no usa red ni binarios externos.

Salidas: ``reports/reference_parity.tsv`` y ``reports/reference_parity.md``.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
ELIGIBLE_STATUSES = ("annotated", "available")
CLASSES = ("CDS", "tRNA", "rRNA")
_ID_RE = re.compile(r"(?:^|;)\s*ID=([^;]+)")


# ---------------------------------------------------------------------------
# Entradas
# ---------------------------------------------------------------------------
def open_text(path: Path):
    """Abre texto plano o gzip según la extensión."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def load_manifest(path: Path | str) -> list[dict[str, Any]]:
    """Lee ``reference_genomes.json`` y devuelve sus registros."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no existe el manifiesto: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"manifiesto ilegible ({path}): {exc}") from exc
    if not isinstance(data, list):
        raise ValueError(f"el manifiesto debe ser una lista: {path}")
    records = [rec for rec in data if isinstance(rec, dict)]
    if not records:
        raise ValueError(f"el manifiesto no contiene registros: {path}")
    return records


def _resolve(root: Path, raw: Any) -> Path | None:
    if raw is None or str(raw).strip() == "":
        return None
    path = Path(str(raw))
    return path if path.is_absolute() else root / path


def _feature_id(attributes: str, fallback: str) -> str:
    match = _ID_RE.search(attributes)
    if match:
        value = match.group(1).strip()
        if value:
            return value
    return fallback


def count_gff_features(path: Path | str) -> dict[str, Any]:
    """Cuenta CDS/tRNA/rRNA de un GFF3 agrupando filas con el mismo ID.

    Un CDS multiparte (varias filas con el mismo ``ID``) cuenta como un
    feature. Las filas sin ``ID`` cuentan una a una. No interpreta
    ``pseudo`` ni ningún otro atributo: solo mira la columna tipo.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no existe el GFF: {path}")

    raw_counts: Counter[str] = Counter()
    feature_ids: dict[str, set[tuple[str, str]]] = {kind: set() for kind in CLASSES}
    invalid_lines = 0
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                invalid_lines += 1
                continue
            seqid, ftype, attributes = fields[0], fields[2], fields[8]
            raw_counts[ftype] += 1
            if ftype in feature_ids:
                fid = _feature_id(attributes, f"__row_{lineno}")
                feature_ids[ftype].add((seqid, fid))

    counts = {kind: len(feature_ids[kind]) for kind in CLASSES}
    return {
        "counts": counts,
        "raw_counts": dict(raw_counts),
        "invalid_lines": invalid_lines,
    }


def read_system_counts(path: Path | str) -> tuple[dict[str, int], dict[str, Any]]:
    """Extrae ``features.CDS/tRNA/rRNA`` del resumen JSON del sistema."""
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"resumen ilegible ({path}): {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"el resumen no es un objeto JSON: {path}")
    features = data.get("features")
    if not isinstance(features, dict):
        raise ValueError(f"el resumen no tiene el bloque 'features': {path}")

    counts: dict[str, int] = {}
    for kind in CLASSES:
        raw = features.get(kind, 0)
        try:
            counts[kind] = int(raw or 0)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"conteo inválido para {kind} en {path}: {raw!r}"
            ) from exc
    return counts, data


# ---------------------------------------------------------------------------
# Comparación
# ---------------------------------------------------------------------------
def _percent(system: int | None, reference: int | None) -> float | None:
    if system is None or reference is None or reference <= 0:
        return None
    return round(100.0 * system / reference, 2)


def compare_record(record: dict[str, Any], root: Path | str) -> dict[str, Any]:
    """Compara un registro del manifiesto contra su resumen y GFF."""
    root = Path(root)
    genome = str(record.get("id") or "")
    summary_raw = str(record.get("system_summary") or "").strip()
    annotation_raw = str(record.get("annotation_path") or "").strip()
    summary_path = _resolve(root, summary_raw)
    annotation_path = _resolve(root, annotation_raw)

    row: dict[str, Any] = {
        "genome": genome,
        "species": str(record.get("species") or ""),
        "strain": str(record.get("strain") or ""),
        "manifest_status": str(record.get("status") or ""),
        "state": "not_available",
        "reason": "",
        # Se muestran las rutas tal como están en el manifiesto (portabilidad).
        "system_summary": summary_raw,
        "annotation_path": annotation_raw,
        "classes": {
            kind: {"system": None, "reference": None, "delta": None, "pct": None}
            for kind in CLASSES
        },
    }

    if summary_path is None:
        row["reason"] = "el manifiesto no declara system_summary"
        return row
    if not summary_path.is_file():
        row["reason"] = f"falta el resumen del sistema: {summary_raw}"
        return row
    try:
        system_counts, _ = read_system_counts(summary_path)
    except ValueError as exc:
        row["reason"] = str(exc)
        return row
    for kind in CLASSES:
        row["classes"][kind]["system"] = system_counts[kind]

    if annotation_path is None:
        row["reason"] = "el manifiesto no declara annotation_path"
        return row
    if not annotation_path.is_file():
        row["reason"] = f"falta el GFF de referencia: {annotation_raw}"
        return row
    try:
        gff = count_gff_features(annotation_path)
    except (OSError, EOFError, ValueError) as exc:
        row["reason"] = f"GFF ilegible ({annotation_raw}): {exc}"
        return row

    reference_counts = gff["counts"]
    if sum(reference_counts.values()) == 0:
        row["reason"] = f"el GFF no contiene CDS/tRNA/rRNA: {annotation_raw}"
        return row

    missing: list[str] = []
    differences: list[str] = []
    for kind in CLASSES:
        system = system_counts[kind]
        reference = reference_counts[kind]
        delta = system - reference
        row["classes"][kind].update(
            {
                "reference": reference,
                "delta": delta,
                "pct": _percent(system, reference),
            }
        )
        if reference == 0:
            missing.append(kind)
        elif delta != 0:
            differences.append(f"{kind}={delta:+d}")

    if missing:
        row["state"] = "partial"
        parts = [f"la referencia no tiene {', '.join(missing)}"]
        if differences:
            parts.append("delta " + ", ".join(differences))
        row["reason"] = "; ".join(parts)
    elif differences:
        row["state"] = "partial"
        row["reason"] = "delta " + ", ".join(differences)
    else:
        row["state"] = "ok"
        row["reason"] = "conteos idénticos en CDS, tRNA y rRNA"
    return row


def _aggregate(rows: Sequence[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    aggregate: dict[str, dict[str, Any]] = {}
    for kind in CLASSES:
        comparable = [
            row
            for row in rows
            if row["classes"][kind]["system"] is not None
            and row["classes"][kind]["reference"] not in (None, 0)
        ]
        system = sum(row["classes"][kind]["system"] for row in comparable)
        reference = sum(row["classes"][kind]["reference"] for row in comparable)
        aggregate[kind] = {
            "system": system,
            "reference": reference,
            "delta": system - reference,
            "pct": _percent(system, reference),
        }
    return aggregate


def build_report(
    manifest_path: Path | str,
    root: Path | str,
    generated_at: str | None = None,
) -> dict[str, Any]:
    """Construye el informe completo a partir del manifiesto y las entradas."""
    manifest_path = Path(manifest_path)
    root = Path(root)
    records = load_manifest(manifest_path)

    selected: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for record in records:
        status = str(record.get("status") or "").strip().lower()
        if status not in ELIGIBLE_STATUSES:
            skipped.append(
                {
                    "genome": str(record.get("id") or ""),
                    "status": status,
                    "reason": f"status '{status}' fuera de annotated/available",
                }
            )
        elif not record.get("system_summary"):
            skipped.append(
                {
                    "genome": str(record.get("id") or ""),
                    "status": status,
                    "reason": "sin system_summary en el manifiesto",
                }
            )
        else:
            selected.append(record)

    rows = [compare_record(record, root) for record in selected]
    state_counts = Counter(row["state"] for row in rows)
    eligible_skipped = [
        item for item in skipped if item["status"] in ELIGIBLE_STATUSES
    ]
    return {
        "generated_at": generated_at
        or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifest": str(manifest_path),
        "root": str(root),
        "rows": rows,
        "skipped": skipped,
        "eligible_skipped": eligible_skipped,
        "summary": {
            "candidates": len(selected) + len(eligible_skipped),
            "selected": len(rows),
            "ok": state_counts.get("ok", 0),
            "partial": state_counts.get("partial", 0),
            "not_available": state_counts.get("not_available", 0),
            "skipped": len(skipped),
        },
        "aggregate": _aggregate(rows),
    }


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------
TSV_COLUMNS = (
    "genome",
    "species",
    "strain",
    "manifest_status",
    "state",
    "cds_system",
    "cds_reference",
    "cds_delta",
    "cds_pct",
    "trna_system",
    "trna_reference",
    "trna_delta",
    "trna_pct",
    "rrna_system",
    "rrna_reference",
    "rrna_delta",
    "rrna_pct",
    "annotation_path",
    "system_summary",
    "reason",
)


def _class_cell(value: Any) -> Any:
    return "" if value is None else value


def _pct_cell(value: float | None) -> str:
    return "" if value is None else f"{value:.2f}"


def write_tsv(rows: Iterable[dict[str, Any]], path: Path | str) -> None:
    """Escribe la tabla TSV máquina-legible."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for row in rows:
            values: list[Any] = [
                row["genome"],
                row["species"],
                row["strain"],
                row["manifest_status"],
                row["state"],
            ]
            for kind in CLASSES:
                cell = row["classes"][kind]
                values.extend(
                    [
                        _class_cell(cell["system"]),
                        _class_cell(cell["reference"]),
                        _class_cell(cell["delta"]),
                        _pct_cell(cell["pct"]),
                    ]
                )
            values.extend(
                [row["annotation_path"], row["system_summary"], row["reason"]]
            )
            writer.writerow(values)


def _md(value: Any) -> str:
    if value is None or value == "":
        return "—"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _md_pct(value: float | None) -> str:
    return "—" if value is None else f"{value:.2f}"


def _md_delta(value: int | None) -> str:
    if value is None:
        return "—"
    return "0" if value == 0 else f"{value:+d}"


def render_markdown(result: dict[str, Any]) -> str:
    """Renderiza el informe legible en Markdown."""
    summary = result["summary"]
    lines: list[str] = [
        "# Paridad de anotación: sistema vs GFF de referencia",
        "",
        f"**Generado:** `{result['generated_at']}`  ",
        f"**Manifiesto:** `{result['manifest']}`  ",
        f"**Raíz de datos:** `{result['root']}`",
        "",
        "## Alcance",
        "",
        "Compara los conteos de **CDS, tRNA y rRNA** del resumen del sistema "
        "(`results/<id>/<id>.summary.json`) contra los features del GFF de "
        "referencia (`annotation_path`) para cada genoma con `status` "
        "`annotated` o `available` y `system_summary` en el manifiesto.",
        "",
        "Esta versión **no evalúa función, productos, coordenadas ni identidad "
        "de secuencia**. Tampoco modifica `results/`, los GFF ni el manifiesto.",
        "",
        "## Cobertura",
        "",
        f"- Registros con status `annotated`/`available`: "
        f"**{summary['candidates']}**",
        f"- Evaluados (con `system_summary`): **{summary['selected']}** "
        f"(ok **{summary['ok']}**, partial **{summary['partial']}**, "
        f"not_available **{summary['not_available']}**)",
        f"- Omitidos: **{summary['skipped']}** "
        f"({len(result['eligible_skipped'])} elegibles sin `system_summary`; "
        f"{summary['skipped'] - len(result['eligible_skipped'])} con otro status)",
        "",
    ]

    if result["eligible_skipped"]:
        lines.extend(
            [
                "### Omitidos con status elegible",
                "",
                "| genoma | status | motivo |",
                "|---|---|---|",
            ]
        )
        for item in result["eligible_skipped"]:
            lines.append(
                f"| `{_md(item['genome'])}` | {_md(item['status'])} | "
                f"{_md(item['reason'])} |"
            )
        lines.append("")

    lines.extend(
        [
            "## Resultados por genoma",
            "",
            "| genoma | estado | CDS sis | CDS ref | Δ CDS | % CDS | "
            "tRNA sis | tRNA ref | Δ tRNA | % tRNA | rRNA sis | rRNA ref | "
            "Δ rRNA | % rRNA |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["rows"]:
        cells = [f"`{_md(row['genome'])}`", _md(row["state"])]
        for kind in CLASSES:
            c = row["classes"][kind]
            cells.extend(
                [
                    _md(c["system"]),
                    _md(c["reference"]),
                    _md_delta(c["delta"]),
                    _md_pct(c["pct"]),
                ]
            )
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    lines.extend(
        [
            "## Estado y motivos",
            "",
            "| genoma | estado | motivo |",
            "|---|---|---|",
        ]
    )
    for row in result["rows"]:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['state'])} | "
            f"{_md(row['reason'])} |"
        )
    lines.append("")

    lines.extend(
        [
            "## Entradas y metadatos",
            "",
            "| genoma | especie | cepa | status manifiesto | GFF de referencia | "
            "resumen del sistema |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in result["rows"]:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['species'])} | "
            f"{_md(row['strain'])} | {_md(row['manifest_status'])} | "
            f"`{_md(row['annotation_path'])}` | "
            f"`{_md(row['system_summary'])}` |"
        )
    lines.append("")

    aggregate = result["aggregate"]
    lines.extend(
        [
            "## Agregado (micro)",
            "",
            "| clase | sistema | referencia | delta | % sistema/referencia |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for kind in CLASSES:
        agg = aggregate[kind]
        lines.append(
            f"| {kind} | {agg['system']} | {agg['reference']} | "
            f"{agg['delta']:+d} | {_md_pct(agg['pct'])} |"
        )
    lines.extend(
        [
            "",
            "El agregado suma solo genomas con conteo de sistema y referencia "
            "disponibles para esa clase.",
            "",
            "## Criterio de estado",
            "",
            "- **ok**: los tres conteos (CDS, tRNA y rRNA) coinciden exactamente "
            "con la referencia y la referencia tiene las tres clases.",
            "- **partial**: la comparación fue posible pero al menos un conteo "
            "difiere, o la referencia carece de alguna de las tres clases.",
            "- **not_available**: falta o es ilegible el resumen, el GFF, o el GFF "
            "no contiene ninguna de las tres clases.",
            "",
            "En los recuentos de referencia, un CDS multiparte (varias filas con "
            "el mismo `ID`) cuenta como un solo feature. No se filtran "
            "pseudogenes por atributo: solo se usa la columna tipo del GFF.",
            "",
            "## Límites",
            "",
            "- Solo mide presencia y número de features; no mide exactitud de "
            "coordenadas, función, productos ni identidad de secuencia.",
            "- El gold hereda los criterios y errores de la fuente RefSeq/NCBI; "
            "no es una verdad absoluta.",
            "- Las diferencias de rRNA pueden reflejar que el pipeline no llama "
            "todas las copias o tipos (p. ej. 23S/5S); este informe no lo "
            "atribuye a una causa sin evidencia adicional.",
            "- No se comparan ncRNA ni otros tipos, aunque estén presentes en el "
            "GFF de referencia.",
            "- Los genomas sin `system_summary` en el manifiesto quedan fuera del "
            "cálculo; se listan arriba cuando su status era elegible.",
            "",
            "## Reproducir",
            "",
            "```bash",
            ".venv/bin/python scripts/reference/reference_parity.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def write_markdown(result: dict[str, Any], path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_markdown(result), encoding="utf-8")


def write_reports(
    result: dict[str, Any], tsv_path: Path | str, md_path: Path | str
) -> tuple[Path, Path]:
    """Escribe ambas salidas y devuelve sus rutas."""
    tsv_path, md_path = Path(tsv_path), Path(md_path)
    write_tsv(result["rows"], tsv_path)
    write_markdown(result, md_path)
    return tsv_path, md_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compara conteos CDS/tRNA/rRNA del sistema contra los GFF de "
            "reference_genomes.json."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="raíz del proyecto (default: la del script)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="manifiesto (default: <root>/data/reference_genomes.json)",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=None,
        help="salida TSV (default: <root>/reports/reference_parity.tsv)",
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=None,
        help="salida Markdown (default: <root>/reports/reference_parity.md)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = (args.root or ROOT).resolve()
    manifest = (args.manifest or root / "data" / "reference_genomes.json").resolve()
    tsv_path = (args.tsv or root / "reports" / "reference_parity.tsv").resolve()
    md_path = (args.md or root / "reports" / "reference_parity.md").resolve()

    try:
        result = build_report(manifest, root)
    except (OSError, EOFError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    write_reports(result, tsv_path, md_path)
    summary = result["summary"]
    print(
        f"{summary['selected']} genomas evaluados: "
        f"{summary['ok']} ok, {summary['partial']} partial, "
        f"{summary['not_available']} not_available; "
        f"{summary['skipped']} omitidos"
    )
    print(f"TSV: {tsv_path}")
    print(f"MD:  {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
