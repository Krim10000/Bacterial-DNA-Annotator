#!/usr/bin/env python3
"""Audita la cobertura de la transferencia K-12 en los genomas del panel.

La transferencia por ortología (``transfer_<genoma>.tsv``) y el modelo
genómico K-12 (``model_<genoma>.complexes.tsv``) son **inferencias**. No son
catálogos ComplexPortal independientes y no se tratan como validación.

El script lee, sin red y sin descargar ComplexPortal:

- ``data/complexes/transfer_<genoma>.tsv``;
- ``data/complexes/model_<genoma>.complexes.tsv``;
- ``data/complexes/gold_k12.tsv`` (catálogo K-12, solo para mapear b-numbers a
  nombres de gen y comparar con las firmas curadas);
- ``data/taxonomy/full/<genoma>.json`` (opcional, solo para aceptar catálogos
  independientes nombrados por taxid);
- ``data/complexes/catalog_<genoma|taxid>.tsv`` o
  ``data/complexes/gold_<genoma|taxid>.tsv`` (catálogo independiente, si existe).

No modifica ``transfer_*.tsv``, ``model_*.tsv`` ni ningún catálogo.

Salidas por defecto:

- ``reports/complex_transfer_coverage.tsv``: una fila por genoma.
- ``reports/complex_transfer_coverage.md``: resumen, estados, cobertura y
  tabla completa.

Estados y convenciones
----------------------
``not_available``
    No hay catálogo ComplexPortal independiente para ese genoma. No significa
    que ComplexPortal no tenga una entrada remota: el script no consulta la
    red. Tampoco se rellena con ``transfer_*.tsv`` ni ``model_*.tsv``.
``transfer_role=inferred_transfer``
    La transferencia K-12 es una hipótesis por ortología, no una validación.
``transfer_suspected``
    Complejo transferido con ``frac_mapped >= 0.8`` (la capa que el pipeline
    llama ``complex_suspected``).
``transfer_candidate``
    Complejo transferido con ``frac_mapped < 0.8``.
``signature_*``
    Coincidencias con las firmas curadas de ``webapp.complexes.SIGNATURES``
    usando los nombres de gen K-12 del catálogo. ``signature_complete`` exige
    todos los genes esenciales de la firma; no confirma el complejo.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
NOT_AVAILABLE = "not_available"
CATALOG_AVAILABLE = "available"
TRANSFER_ROLE = "inferred_transfer"
SIGNATURE_EVALUATED = "evaluated"

TSV_COLUMNS = (
    "genome",
    "independent_catalog",
    "catalog_status",
    "catalog_complexes",
    "transfer_role",
    "transfer_complexes",
    "transfer_suspected",
    "transfer_candidate",
    "transfer_stoich_complete",
    "transfer_stoich_incomplete",
    "transfer_stoich_unknown",
    "transfer_full_subunits",
    "transfer_with_ec",
    "subunits_expected",
    "subunits_mapped",
    "subunit_coverage",
    "mean_frac_mapped",
    "mean_identity",
    "signature_status",
    "signature_matches",
    "signature_rows",
    "signature_complete",
    "signature_incomplete",
    "signature_ids",
    "model_complexes",
    "model_candidates",
    "model_complete",
    "model_incomplete",
    "model_other",
    "model_subunits",
    "limitations",
)

STATUS_MEANING = {
    "transfer_suspected": "complejo transferido con frac_mapped >= 0.8; sospecha fuerte",
    "transfer_candidate": "complejo transferido con frac_mapped < 0.8; candidato",
    "transfer_stoich_complete": "observed_n >= expected_n en el TSV de transferencia",
    "transfer_stoich_incomplete": "expected_n presente y observed_n < expected_n",
    "transfer_stoich_unknown": "sin expected_n; estequiometría no evaluable",
    "transfer_full_subunits": "n_mapped == n_subunits del catálogo K-12",
    "signature_complete": "todos los genes esenciales de la firma están en la transferencia",
    "signature_incomplete": "firma con >=2 genes pero sin todos los esenciales",
    "model_candidates": "filas del modelo K-12 con status complex_candidate",
    "model_complete": "filas del modelo K-12 con status complex_complete",
    "model_incomplete": "filas del modelo K-12 con status complex_incomplete",
    "model_other": "filas del modelo K-12 con otro status",
}


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _int(value: object, default: int = 0) -> int:
    try:
        return int(float(_text(value)))
    except (TypeError, ValueError):
        return default


def _float(value: object, default: float = 0.0) -> float:
    try:
        return float(_text(value))
    except (TypeError, ValueError):
        return default


def _ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def _num(value: int) -> str:
    """Entero con separador de miles estilo es-CL (19.085)."""
    return f"{value:,}".replace(",", ".")


def _fmt_ratio(value: float) -> str:
    return f"{value:.4f}"


def _dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = _text(value).replace("\t", " ").replace("\n", " ")
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _taxid_from_taxonomy(root: Path, genome: str) -> int | None:
    """Taxid del ensamblado, si el JSON local existe. No consulta la red."""
    for relative in ("data/taxonomy/full", "data/taxonomy"):
        path = root / relative / f"{genome}.json"
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if not isinstance(data, dict):
            return None
        candidates: list[object] = [data.get("ncbi_taxid"), data.get("taxid")]
        evidence = data.get("evidence")
        if isinstance(evidence, dict):
            ncbi = evidence.get("ncbi")
            if isinstance(ncbi, dict):
                candidates.append(ncbi.get("taxid"))
        for value in candidates:
            number = _int(value, default=0)
            if number > 0:
                return number
        return None
    return None


def _catalog_rows(path: Path) -> int:
    """Cuenta filas con algún contenido, sin contar la cabecera."""
    try:
        with path.open(encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader, None)
            return sum(1 for row in reader if any(_text(cell) for cell in row))
    except OSError:
        return 0


def find_independent_catalog(root: Path, genome: str) -> Path | None:
    """Busca un catálogo independiente nombrado por genoma o taxid.

    Nunca devuelve ``gold_k12.tsv`` para un genoma del panel: el catálogo
    K-12 es el prior de la transferencia y no es una validación independiente.
    """
    complexes = root / "data" / "complexes"
    suffixes = [genome]
    taxid = _taxid_from_taxonomy(root, genome)
    if taxid is not None:
        suffixes.append(str(taxid))
    for suffix in suffixes:
        for directory in (complexes, complexes / "catalogs"):
            for name in (f"catalog_{suffix}.tsv", f"gold_{suffix}.tsv"):
                path = directory / name
                if path.is_file() and _catalog_rows(path) > 0:
                    return path
    return None


def load_signatures() -> list[dict[str, list[str] | str]]:
    """Firmas curadas de ``webapp.complexes``, normalizadas a minúsculas."""
    try:
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        from webapp.complexes import SIGNATURES  # noqa: PLC0415
    except Exception:  # pragma: no cover - solo si el paquete no está instalado
        return []
    out: list[dict[str, list[str] | str]] = []
    for signature in SIGNATURES:
        genes = [_text(gene).lower() for gene in signature.get("genes", [])]
        essential = [
            _text(gene).lower() for gene in signature.get("essential", [])
        ]
        genes = [gene for gene in genes if gene]
        essential = [gene for gene in essential if gene]
        signature_id = _text(signature.get("id"))
        if signature_id and genes and essential:
            out.append({
                "id": signature_id,
                "genes": genes,
                "essential": essential,
            })
    return out


def read_k12_gene_names(path: Path) -> dict[str, str]:
    """b-number -> nombre de gen, desde ``gold_k12.tsv``."""
    mapping: dict[str, str] = {}
    if not path.is_file():
        return mapping
    try:
        with path.open(encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                tags = [_text(tag) for tag in
                        _text(row.get("locus_tags")).split(",")]
                names = [_text(name) for name in
                         _text(row.get("gene_names")).split(",")]
                for tag, name in zip(tags, names):
                    if tag and name:
                        mapping[tag.lower()] = name.lower()
    except (OSError, csv.Error):
        return {}
    return mapping


def signature_hits(
    genes_k12: str,
    gene_by_b: dict[str, str],
    signatures: list[dict[str, list[str] | str]],
) -> list[tuple[str, bool]]:
    """Devuelve (signature_id, completo) para una fila de transferencia."""
    present: set[str] = set()
    for tag in _text(genes_k12).split(","):
        name = gene_by_b.get(tag.strip().lower(), "")
        if name:
            present.add(name)
    hits: list[tuple[str, bool]] = []
    for signature in signatures:
        genes = set(signature["genes"])  # type: ignore[arg-type]
        essential = list(signature["essential"])  # type: ignore[arg-type]
        found = present & genes
        if len(found) < 2:
            continue
        complete = all(str(gene) in present for gene in essential)
        hits.append((str(signature["id"]), complete))
    return hits


def _model_subunits(value: object) -> int:
    """Cuenta subunidades únicas de la columna ``subunits`` del modelo."""
    ids: set[str] = set()
    for part in _text(value).split(";"):
        part = part.strip()
        if not part:
            continue
        ids.add(part.split(":", 1)[1].strip() if ":" in part else part)
    return len(ids)


def discover_genomes(complexes: Path) -> list[str]:
    """Genomas con transferencia, modelo o ambos; ordenados."""
    genomes: set[str] = set()
    if not complexes.is_dir():
        return []
    for path in complexes.iterdir():
        if not path.is_file() or path.name.startswith("."):
            continue
        match = re.match(r"^transfer_(.+)\.tsv$", path.name)
        if match:
            genomes.add(match.group(1))
            continue
        match = re.match(r"^model_(.+)\.complexes\.tsv$", path.name)
        if match:
            genomes.add(match.group(1))
    return sorted(genomes)


def _transfer_stats(
    path: Path,
    *,
    gene_by_b: dict[str, str],
    signatures: list[dict[str, list[str] | str]],
) -> dict[str, Any]:
    stats: Counter[str] = Counter()
    expected = mapped = 0
    fractions: list[float] = []
    identities: list[float] = []
    signature_ids: set[str] = set()
    limitations: list[str] = []
    if not path.is_file():
        return {
            "stats": stats,
            "subunits_expected": 0,
            "subunits_mapped": 0,
            "fractions": fractions,
            "identities": identities,
            "signature_ids": signature_ids,
            "limitations": limitations,
        }
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stats["transfer_complexes"] += 1
            n_subunits = max(0, _int(row.get("n_subunits")))
            n_mapped = max(0, _int(row.get("n_mapped")))
            expected += n_subunits
            mapped += n_mapped
            fraction = _float(row.get("frac_mapped"))
            if not _text(row.get("frac_mapped")):
                fraction = _ratio(n_mapped, n_subunits)
            fractions.append(fraction)
            if fraction >= 0.8:
                stats["transfer_suspected"] += 1
            else:
                stats["transfer_candidate"] += 1
            stoich = _text(row.get("stoich_status")).lower()
            if stoich not in {"complete", "incomplete", "unknown"}:
                stoich = "unknown"
            stats[f"transfer_stoich_{stoich}"] += 1
            if n_subunits and n_mapped >= n_subunits:
                stats["transfer_full_subunits"] += 1
            if _text(row.get("ecs")):
                stats["transfer_with_ec"] += 1
            if _text(row.get("identity_mean")):
                identities.append(_float(row.get("identity_mean")))
            limitation = _text(row.get("limitation"))
            if limitation:
                limitations.append(limitation)
            hits = signature_hits(row.get("genes_k12"), gene_by_b, signatures)
            if hits:
                stats["signature_rows"] += 1
            for signature_id, complete in hits:
                stats["signature_matches"] += 1
                signature_ids.add(signature_id)
                if complete:
                    stats["signature_complete"] += 1
                else:
                    stats["signature_incomplete"] += 1
    return {
        "stats": stats,
        "subunits_expected": expected,
        "subunits_mapped": mapped,
        "fractions": fractions,
        "identities": identities,
        "signature_ids": signature_ids,
        "limitations": limitations,
    }


def _model_stats(path: Path) -> dict[str, Any]:
    stats: Counter[str] = Counter()
    subunits = 0
    limitations: list[str] = []
    if not path.is_file():
        return {"stats": stats, "subunits": 0, "limitations": limitations}
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stats["model_complexes"] += 1
            status = _text(row.get("status")).lower()
            if status == "complex_candidate":
                stats["model_candidates"] += 1
            elif status == "complex_complete":
                stats["model_complete"] += 1
            elif status == "complex_incomplete":
                stats["model_incomplete"] += 1
            else:
                stats["model_other"] += 1
            subunits += _model_subunits(row.get("subunits"))
            limitation = _text(row.get("limitation"))
            if limitation:
                limitations.append(limitation)
    return {"stats": stats, "subunits": subunits, "limitations": limitations}


def audit_genome(
    root: Path,
    genome: str,
    *,
    gene_by_b: dict[str, str],
    signatures: list[dict[str, list[str] | str]],
) -> dict[str, Any]:
    complexes = root / "data" / "complexes"
    transfer = _transfer_stats(
        complexes / f"transfer_{genome}.tsv",
        gene_by_b=gene_by_b,
        signatures=signatures,
    )
    model = _model_stats(complexes / f"model_{genome}.complexes.tsv")
    catalog = find_independent_catalog(root, genome)
    catalog_rows = _catalog_rows(catalog) if catalog is not None else 0
    if catalog is not None and catalog_rows == 0:
        catalog = None

    tstats: Counter[str] = transfer["stats"]
    mstats: Counter[str] = model["stats"]
    fractions: list[float] = transfer["fractions"]
    identities: list[float] = transfer["identities"]
    coverage = _ratio(transfer["subunits_mapped"], transfer["subunits_expected"])
    mean_fraction = sum(fractions) / len(fractions) if fractions else 0.0
    mean_identity = sum(identities) / len(identities) if identities else 0.0

    limitations: list[str] = []
    if catalog is None:
        limitations.append(
            "sin catálogo ComplexPortal independiente; not_available"
        )
    else:
        limitations.append(
            "catálogo independiente disponible; no se usa como validación"
        )
    limitations.append("transferencia K-12 = inferencia; no es validación")
    if tstats["transfer_complexes"] == 0:
        limitations.append("sin complejos transferidos")
    if mstats["model_complexes"] == 0:
        limitations.append("sin candidatos del modelo K-12")
    if not gene_by_b or not signatures:
        limitations.append("firmas no evaluables: falta gold_k12.tsv o firmas")
    if tstats["transfer_stoich_unknown"]:
        limitations.append(
            "estequiometría no evaluable en "
            f"{tstats['transfer_stoich_unknown']} complejos"
        )
    limitations.extend(transfer["limitations"])
    limitations.extend(model["limitations"])

    try:
        catalog_label = (
            catalog.relative_to(root).as_posix()
            if catalog is not None else NOT_AVAILABLE
        )
    except ValueError:
        catalog_label = catalog.as_posix() if catalog is not None else NOT_AVAILABLE

    return {
        "genome": genome,
        "independent_catalog": catalog_label,
        "catalog_status": (
            CATALOG_AVAILABLE if catalog is not None else NOT_AVAILABLE
        ),
        "catalog_complexes": catalog_rows,
        "transfer_role": TRANSFER_ROLE,
        "transfer_complexes": tstats["transfer_complexes"],
        "transfer_suspected": tstats["transfer_suspected"],
        "transfer_candidate": tstats["transfer_candidate"],
        "transfer_stoich_complete": tstats["transfer_stoich_complete"],
        "transfer_stoich_incomplete": tstats["transfer_stoich_incomplete"],
        "transfer_stoich_unknown": tstats["transfer_stoich_unknown"],
        "transfer_full_subunits": tstats["transfer_full_subunits"],
        "transfer_with_ec": tstats["transfer_with_ec"],
        "subunits_expected": transfer["subunits_expected"],
        "subunits_mapped": transfer["subunits_mapped"],
        "subunit_coverage": coverage,
        "mean_frac_mapped": mean_fraction,
        "mean_identity": mean_identity,
        "signature_status": (
            SIGNATURE_EVALUATED if gene_by_b and signatures else NOT_AVAILABLE
        ),
        "signature_matches": tstats["signature_matches"],
        "signature_rows": tstats["signature_rows"],
        "signature_complete": tstats["signature_complete"],
        "signature_incomplete": tstats["signature_incomplete"],
        "signature_ids": ";".join(sorted(transfer["signature_ids"])),
        "model_complexes": mstats["model_complexes"],
        "model_candidates": mstats["model_candidates"],
        "model_complete": mstats["model_complete"],
        "model_incomplete": mstats["model_incomplete"],
        "model_other": mstats["model_other"],
        "model_subunits": model["subunits"],
        "limitations": "; ".join(_dedupe(limitations)),
    }


def build_audit(root: Path = ROOT) -> list[dict[str, Any]]:
    """Construye una fila por genoma. No escribe ni descarga nada."""
    root = Path(root)
    complexes = root / "data" / "complexes"
    gene_by_b = read_k12_gene_names(complexes / "gold_k12.tsv")
    signatures = load_signatures()
    return [
        audit_genome(
            root, genome, gene_by_b=gene_by_b, signatures=signatures,
        )
        for genome in discover_genomes(complexes)
    ]


def _format_value(column: str, value: object) -> str:
    if column in {"subunit_coverage", "mean_frac_mapped"}:
        return _fmt_ratio(float(value))
    if column == "mean_identity":
        return f"{float(value):.1f}"
    return str(value)


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                column: _format_value(column, row.get(column, ""))
                for column in TSV_COLUMNS
            })


def _md(value: object) -> str:
    return _text(value).replace("|", "\\|").replace("\n", " ")


def _sum(rows: list[dict[str, Any]], key: str) -> int:
    return sum(int(row.get(key, 0)) for row in rows)


def render_markdown(rows: list[dict[str, Any]], generated_at: str) -> str:
    total_transfer = _sum(rows, "transfer_complexes")
    total_suspected = _sum(rows, "transfer_suspected")
    total_candidate = _sum(rows, "transfer_candidate")
    total_expected = _sum(rows, "subunits_expected")
    total_mapped = _sum(rows, "subunits_mapped")
    coverage = _ratio(total_mapped, total_expected)
    catalog_available = sum(
        1 for row in rows if row["catalog_status"] == CATALOG_AVAILABLE
    )
    catalog_missing = len(rows) - catalog_available
    signature_genomes = sum(
        1 for row in rows if int(row.get("signature_matches", 0)) > 0
    )
    without_transfer = sum(
        1 for row in rows if int(row.get("transfer_complexes", 0)) == 0
    )
    without_model = sum(
        1 for row in rows if int(row.get("model_complexes", 0)) == 0
    )
    total_model = _sum(rows, "model_complexes")
    total_signature_matches = _sum(rows, "signature_matches")
    total_signature_complete = _sum(rows, "signature_complete")
    total_signature_incomplete = _sum(rows, "signature_incomplete")

    lines = [
        "# Cobertura de la transferencia de complejos K-12",
        "",
        f"Generado: `{generated_at}`",
        "",
        "Alcance: genomas con `data/complexes/transfer_*.tsv` o "
        "`data/complexes/model_*.complexes.tsv`. Solo lectura local: no se "
        "descargó ComplexPortal, no se consultó la red y no se modificó ningún "
        "complejo.",
        "",
        "`not_available` significa que **no hay catálogo ComplexPortal "
        "independiente local** para ese genoma. No significa que ComplexPortal "
        "no tenga una entrada remota. La transferencia K-12 y el modelo K-12 "
        "son inferencias: no son validación.",
        "",
        "## Resumen",
        "",
        "| métrica | valor |",
        "|---|---:|",
        f"| genomas auditados | {_num(len(rows))} |",
        f"| con catálogo ComplexPortal independiente | {_num(catalog_available)} |",
        f"| sin catálogo independiente (`not_available`) | {_num(catalog_missing)} |",
        f"| genomas con transferencia | {_num(len(rows) - without_transfer)} |",
        f"| genomas sin transferencia | {_num(without_transfer)} |",
        f"| genomas con modelo K-12 | {_num(len(rows) - without_model)} |",
        f"| genomas sin modelo K-12 | {_num(without_model)} |",
        f"| complejos transferidos | {_num(total_transfer)} |",
        f"| transferidos con `frac_mapped >= 0.8` (sospechosos) | {_num(total_suspected)} |",
        f"| transferidos con `frac_mapped < 0.8` (candidatos) | {_num(total_candidate)} |",
        f"| subunidades esperadas (suma `n_subunits`) | {_num(total_expected)} |",
        f"| subunidades mapeadas (suma `n_mapped`) | {_num(total_mapped)} |",
        f"| cobertura global de subunidades | {_fmt_ratio(coverage)} |",
        f"| candidatos del modelo K-12 | {_num(total_model)} |",
        f"| coincidencias con firmas curadas | {_num(total_signature_matches)} |",
        f"| firmas completas / incompletas | {_num(total_signature_complete)} / {_num(total_signature_incomplete)} |",
        f"| genomas con alguna firma | {_num(signature_genomes)} |",
        "",
        "## Estados",
        "",
        "| estado | n | significado |",
        "|---|---:|---|",
    ]
    for key in (
        "transfer_suspected",
        "transfer_candidate",
        "transfer_stoich_complete",
        "transfer_stoich_incomplete",
        "transfer_stoich_unknown",
        "transfer_full_subunits",
        "signature_complete",
        "signature_incomplete",
        "model_candidates",
        "model_complete",
        "model_incomplete",
        "model_other",
    ):
        lines.append(
            f"| `{key}` | {_num(_sum(rows, key))} | {STATUS_MEANING[key]} |"
        )

    lines.extend([
        "",
        "`transfer_full_subunits` cuenta filas con `n_mapped == n_subunits`; "
        "`transfer_stoich_complete` cuenta la estequiometría declarada "
        "(`observed_n >= expected_n`). No son lo mismo: el catálogo K-12 puede "
        "tener más subunidades que el `expected_n` de `properties`.",
        "",
        "## Cobertura de subunidades",
        "",
        "| métrica | valor |",
        "|---|---:|",
        f"| subunidades esperadas | {_num(total_expected)} |",
        f"| subunidades mapeadas | {_num(total_mapped)} |",
        f"| cobertura global | {_fmt_ratio(coverage)} |",
        f"| media no ponderada de las medias de `frac_mapped` por genoma | {_fmt_ratio(sum(float(r.get('mean_frac_mapped', 0.0)) for r in rows) / len(rows) if rows else 0.0)} |",
        f"| complejos con todas las subunidades del catálogo | {_num(_sum(rows, 'transfer_full_subunits'))} |",
        "",
        "## Genomas con más complejos transferidos",
        "",
        "| genome | transferidos | sospechosos | cobertura | firmas | modelo K-12 |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    top = sorted(
        rows,
        key=lambda row: (
            int(row.get("transfer_complexes", 0)),
            float(row.get("subunit_coverage", 0.0)),
        ),
        reverse=True,
    )[:10]
    for row in top:
        lines.append(
            f"| `{_md(row['genome'])}` | {_num(int(row['transfer_complexes']))} | "
            f"{_num(int(row['transfer_suspected']))} | "
            f"{_fmt_ratio(float(row['subunit_coverage']))} | "
            f"{_num(int(row['signature_matches']))} | "
            f"{_num(int(row['model_complexes']))} |"
        )

    lines.extend([
        "",
        "## Genomas con menor cobertura de subunidades",
        "",
        "| genome | transferidos | cobertura | subunidades mapeadas / esperadas |",
        "|---|---:|---:|---:|",
    ])
    low = sorted(
        (row for row in rows if int(row.get("transfer_complexes", 0)) > 0),
        key=lambda row: float(row.get("subunit_coverage", 0.0)),
    )[:10]
    for row in low:
        lines.append(
            f"| `{_md(row['genome'])}` | {_num(int(row['transfer_complexes']))} | "
            f"{_fmt_ratio(float(row['subunit_coverage']))} | "
            f"{_num(int(row['subunits_mapped']))} / "
            f"{_num(int(row['subunits_expected']))} |"
        )

    lines.extend([
        "",
        "## Limitaciones",
        "",
        "- **Sin catálogo independiente.** Los genomas del panel no tienen un "
        "catálogo ComplexPortal propio; el único catálogo local es el de K-12. "
        "Por eso la columna `independent_catalog` queda en `not_available` y "
        "no se usa como validación.",
        "- **La transferencia no es validación.** Es RBH K-12, con identidad "
        ">= 30 % y evalue <= 1e-5, pero sin cobertura de alineamiento en los "
        "TSV actuales. `frac_mapped` bajo puede significar complejo ausente, "
        "anotación incompleta u ortología no detectada.",
        "- **Estequiometría.** `expected_n` solo está presente en una parte de "
        "las filas; `stoich_status=unknown` no es un fallo del complejo, es un "
        "dato no evaluable. No hay ninguna fila `complete` en el panel actual.",
        "- **Firmas curadas.** Se comparan solo nombres de gen K-12; no son un "
        "catálogo exhaustivo ni prueban interacción física. "
        "`signature_complete` significa que están los genes esenciales de la "
        "firma en la transferencia, no que el complejo exista.",
        "- **Modelo K-12.** Cada fila es `complex_candidate`: una sospecha por "
        "vecindad/operón, sin estructura ni estequiometría. No cierra brechas.",
        "- **Sin red.** No se descargó ComplexPortal ni se consultó el listado "
        "remoto. Un `not_available` es una ausencia local, no una afirmación "
        "sobre la base remota.",
        "",
        "## Detalle por genoma",
        "",
        "| genome | catálogo independiente | transferidos | sospechosos | cobertura | firmas (compl./inc.) | modelo K-12 | estoq. (C/I/U) |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            f"| `{_md(row['genome'])}` | `{_md(row['independent_catalog'])}` | "
            f"{_num(int(row['transfer_complexes']))} | "
            f"{_num(int(row['transfer_suspected']))} | "
            f"{_fmt_ratio(float(row['subunit_coverage']))} | "
            f"{_num(int(row['signature_complete']))} / "
            f"{_num(int(row['signature_incomplete']))} | "
            f"{_num(int(row['model_complexes']))} | "
            f"{_num(int(row['transfer_stoich_complete']))} / "
            f"{_num(int(row['transfer_stoich_incomplete']))} / "
            f"{_num(int(row['transfer_stoich_unknown']))} |"
        )

    lines.extend([
        "",
        "## Regenerar",
        "",
        "```bash",
        "python scripts/complexes/transfer_coverage_audit.py",
        "```",
        "",
        "El detalle completo de limitaciones y firmas está en "
        "`reports/complex_transfer_coverage.tsv`.",
        "",
    ])
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audita la cobertura de la transferencia K-12 (sin red)."
    )
    parser.add_argument("--root", type=Path, default=ROOT,
                        help="raíz del proyecto (default: la del script)")
    parser.add_argument("--out-tsv", type=Path, default=None,
                        help="default: <root>/reports/complex_transfer_coverage.tsv")
    parser.add_argument("--out-md", type=Path, default=None,
                        help="default: <root>/reports/complex_transfer_coverage.md")
    parser.add_argument("--generated-at", default=None,
                        help="timestamp fijo (para reproducibilidad/tests)")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    rows = build_audit(root)
    out_tsv = args.out_tsv or (
        root / "reports" / "complex_transfer_coverage.tsv"
    )
    out_md = args.out_md or (
        root / "reports" / "complex_transfer_coverage.md"
    )
    generated_at = args.generated_at or datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )

    write_tsv(rows, out_tsv)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(render_markdown(rows, generated_at), encoding="utf-8")

    catalog_available = sum(
        1 for row in rows if row["catalog_status"] == CATALOG_AVAILABLE
    )
    print(f"genomas: {len(rows)}")
    print(
        f"catálogo independiente: {catalog_available} "
        f"(not_available: {len(rows) - catalog_available})"
    )
    print(f"complejos transferidos: {_sum(rows, 'transfer_complexes')}")
    print(
        "sospechosos/candidatos: "
        f"{_sum(rows, 'transfer_suspected')}/"
        f"{_sum(rows, 'transfer_candidate')}"
    )
    print(
        "subunidades mapeadas/esperadas: "
        f"{_sum(rows, 'subunits_mapped')}/{_sum(rows, 'subunits_expected')}"
    )
    print(f"modelo K-12: {_sum(rows, 'model_complexes')}")
    print(f"firmas: {_sum(rows, 'signature_matches')}")
    print(f"escrito: {out_tsv}")
    print(f"escrito: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
