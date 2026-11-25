#!/usr/bin/env python3
"""Auditoría de estado del catálogo AFDB/PAE (solo lectura).

- Lee el catálogo, el universo de accesiones y la caché local.
- Clasifica cada fila en ``ok`` / ``rate_limited`` / ``not_found`` / ``error``.
- Trata 429 como reintentable y 404 como permanente.
- Calcula accesiones pendientes y cobertura sobre el universo (111.569).
- Genera ``reports/afdb_status.tsv`` y ``reports/afdb_status.md``.
- **No** hace red, **no** descarga y **no** reescribe el catálogo ni la caché.

Ejemplos:
  # Informe y artefactos, sin red ni escrituras en la caché
  python scripts/structure/afdb_status.py

  # Solo imprimir el plan y el siguiente lote; no escribe reportes
  python scripts/structure/afdb_status.py --dry-run

  # Usar un catálogo sintético (pruebas)
  python scripts/structure/afdb_status.py --catalog /tmp/cat.tsv \\
      --targets /tmp/targets --reports-dir /tmp/reports --dry-run
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
# Permite ``python scripts/structure/afdb_status.py`` desde cualquier cwd.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.structure.fetch_afdb_cache import (  # noqa: E402
    is_permanent as fetcher_is_permanent,
    load_catalog,
)

STATUS_OK = "ok"
STATUS_RATE_LIMITED = "rate_limited"
STATUS_NOT_FOUND = "not_found"
STATUS_ERROR = "error"
STATUS_LOCAL_ONLY = "local_only"

# Los reintentables se muestran primero en el TSV.
STATUS_ORDER = {
    STATUS_RATE_LIMITED: 0,
    STATUS_ERROR: 1,
    STATUS_NOT_FOUND: 2,
    STATUS_LOCAL_ONLY: 3,
    STATUS_OK: 4,
}

TSV_COLS = [
    "accession", "status", "catalog_status", "local_pdb", "local_pae",
    "pending", "retryable", "permanent", "in_catalog", "in_universe",
    "http_code", "pdb_path", "pae_path", "error",
]

_HTTP_CODE_RE = re.compile(r"HTTP(?:\s+Error)?[:\s]+(\d{3})", re.IGNORECASE)
_RATE_LIMIT_RE = re.compile(r"\b429\b|too many requests|rate[ -]?limit",
                            re.IGNORECASE)
_NOT_FOUND_RE = re.compile(r"\b404\b|not found", re.IGNORECASE)

# Supuestos del ETA; se documentan en el informe.
LATENCY_PER_REQUEST = 0.5
REQUESTS_PER_ACCESSION = 2  # metadatos + PDB (sin PAE)
BATCH_PAUSE = 30.0
FULL_DEBT_PAUSE = 30.0

__all__ = [
    "AccessionStatus",
    "Audit",
    "build_audit",
    "build_plan",
    "classify_catalog_row",
    "estimate_full_pending",
    "estimate_seconds",
    "extract_http_code",
    "load_universe",
    "main",
    "next_batch",
    "recommended_limit",
    "render_markdown",
    "scan_local",
    "summarize",
    "write_reports",
]


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def _rel(path: str | os.PathLike[str]) -> str:
    """Ruta relativa a la raíz del proyecto cuando sea posible."""
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except (OSError, ValueError):
        return str(path)


def extract_http_code(text: Any) -> int | None:
    """Extrae el código de un mensaje tipo ``HTTPError: HTTP 429: ...``."""
    match = _HTTP_CODE_RE.search(str(text or ""))
    return int(match.group(1)) if match else None


def classify_catalog_row(row: Mapping[str, Any],
                         has_local_pdb: bool) -> tuple[str, int | None]:
    """Devuelve ``(estado_efectivo, código_http)`` para una fila de catálogo.

    Reglas:
    - ``ok``/``already_present`` con archivo local -> ``ok``.
    - ``ok``/``already_present`` sin archivo -> ``error`` (integridad).
    - ``not_found``, HTTP 404 o texto "not found" -> ``not_found`` (permanente).
    - ``rate_limited``, HTTP 429 o texto de rate limit -> ``rate_limited``.
    - cualquier otro ``error`` -> ``error`` (reintentable).
    """
    raw = str(row.get("status") or "").strip().lower()
    error = str(row.get("error") or "")
    code = extract_http_code(error)
    if raw in {STATUS_OK, "already_present"}:
        return (STATUS_OK if has_local_pdb else STATUS_ERROR), code
    if raw == STATUS_NOT_FOUND or code == 404 or _NOT_FOUND_RE.search(error):
        return STATUS_NOT_FOUND, code
    if raw == STATUS_RATE_LIMITED or code == 429 or _RATE_LIMIT_RE.search(error):
        return STATUS_RATE_LIMITED, code
    return STATUS_ERROR, code


def _as_bool(value: bool) -> str:
    return "true" if value else "false"


# ---------------------------------------------------------------------------
# Modelo
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AccessionStatus:
    """Estado efectivo de una accesión con evidencia local o de catálogo."""

    accession: str
    status: str
    catalog_status: str = ""
    error: str = ""
    http_code: int | None = None
    local_pdb: str = ""
    local_pae: str = ""
    in_catalog: bool = False
    in_universe: bool = False

    @property
    def retryable(self) -> bool:
        return self.status in {STATUS_RATE_LIMITED, STATUS_ERROR}

    @property
    def permanent(self) -> bool:
        return self.status == STATUS_NOT_FOUND

    @property
    def pending(self) -> bool:
        return (self.in_universe and not self.local_pdb and not self.permanent)


@dataclass
class Audit:
    """Resultado de leer catálogo + caché + universo, sin modificar nada."""

    universe: list[str]
    catalog_path: Path
    targets_path: Path
    catalog_rows: dict[str, dict[str, str]]
    local_pdb: dict[str, Path]
    local_pae: dict[str, Path]
    items: dict[str, AccessionStatus]
    pending: list[str]

    @property
    def n_universe(self) -> int:
        return len(self.universe)


# ---------------------------------------------------------------------------
# Lectura local
# ---------------------------------------------------------------------------
def load_universe(data_root: str | os.PathLike[str]) -> list[str]:
    """Accesiones únicas de ``data/regulation/*/genes.tsv`` (solo lectura).

    Es una variante de ``fetch_afdb_cache.accessions`` que evita construir un
    diccionario por fila: el universo tiene ~1,37 M de filas y la auditoría
    debe ser rápida. Busca la columna ``acc`` en la cabecera y tolera archivos
    sin esa columna.
    """
    seen: set[str] = set()
    for path in sorted((Path(data_root) / "regulation").glob("*/genes.tsv")):
        try:
            with path.open(encoding="utf-8", errors="replace", newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                try:
                    header = next(reader)
                except StopIteration:
                    continue
                try:
                    index = header.index("acc")
                except ValueError:
                    continue
                for row in reader:
                    if len(row) > index:
                        accession = row[index].strip()
                        if accession:
                            seen.add(accession)
        except OSError:
            continue
    return sorted(seen)


def scan_local(targets: str | os.PathLike[str]
               ) -> tuple[dict[str, Path], dict[str, Path]]:
    """Indexa los PDB/PAE locales (planos o ``.gz``) sin modificar nada.

    Prefiere el archivo sin comprimir si existen ambos, igual que
    ``cached_path`` en el descargador.
    """
    pdb: dict[str, Path] = {}
    pae: dict[str, Path] = {}
    base = Path(targets)
    if not base.is_dir():
        return pdb, pae
    suffixes = (
        (".pdb.gz", "pdb"),
        (".pdb", "pdb"),
        (".pae.json.gz", "pae"),
        (".pae.json", "pae"),
    )
    for entry in base.iterdir():
        try:
            if not entry.is_file():
                continue
        except OSError:
            continue
        name = entry.name
        kind = ""
        accession = ""
        for suffix, suffix_kind in suffixes:
            if name.endswith(suffix):
                accession = name[:-len(suffix)]
                kind = suffix_kind
                break
        if not accession or not kind:
            continue
        bucket = pdb if kind == "pdb" else pae
        current = bucket.get(accession)
        if current is None or (current.suffix == ".gz" and entry.suffix != ".gz"):
            bucket[accession] = entry
    return pdb, pae


def build_audit(universe: list[str], catalog: Mapping[str, Mapping[str, Any]],
                targets: str | os.PathLike[str],
                catalog_path: str | os.PathLike[str] | None = None) -> Audit:
    """Cruza universo, catálogo y caché local. Todo es lectura."""
    universe_sorted = sorted(set(universe))
    universe_set = set(universe_sorted)
    pdb, pae = scan_local(targets)

    names = sorted(set(catalog) | set(pdb) | set(pae))
    items: dict[str, AccessionStatus] = {}
    for accession in names:
        row = catalog.get(accession)
        has_pdb = accession in pdb
        if row is not None:
            status, code = classify_catalog_row(row, has_pdb)
            raw_status = str(row.get("status") or "")
            error = str(row.get("error") or "")
        else:
            status, code = STATUS_LOCAL_ONLY, None
            raw_status, error = "", ""
        items[accession] = AccessionStatus(
            accession=accession,
            status=status,
            catalog_status=raw_status,
            error=error,
            http_code=code,
            local_pdb=str(pdb[accession]) if accession in pdb else "",
            local_pae=str(pae[accession]) if accession in pae else "",
            in_catalog=row is not None,
            in_universe=accession in universe_set,
        )

    # Las accesiones del universo sin fila de catálogo ni archivos también
    # están pendientes y deben contarse; no se materializan una por una.
    known = set(items)
    pending = [accession for accession, item in items.items() if item.pending]
    pending.extend(universe_set - known)
    pending.sort()
    return Audit(
        universe=universe_sorted,
        catalog_path=Path(catalog_path) if catalog_path is not None else Path(""),
        targets_path=Path(targets),
        catalog_rows=dict(catalog),
        local_pdb=pdb,
        local_pae=pae,
        items=items,
        pending=pending,
    )


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
def summarize(audit: Audit) -> dict[str, Any]:
    """Cuenta estados, cobertura y pendientes sobre el universo."""
    catalog_items = [item for item in audit.items.values() if item.in_catalog]
    catalog_counts = Counter(item.status for item in catalog_items)
    n_universe = audit.n_universe
    denominator = n_universe or 1

    local_pdb = len(audit.local_pdb)
    local_pae = len(audit.local_pae)
    local_pdb_set = set(audit.local_pdb)
    local_pae_set = set(audit.local_pae)
    covered_in_universe = len(local_pdb_set & set(audit.universe))
    pae_in_local = len(local_pdb_set & local_pae_set)

    catalog_ok_with_file = sum(
        1 for item in catalog_items
        if item.status == STATUS_OK and item.local_pdb and item.in_universe
    )
    catalog_ok_missing_file = sum(
        1 for item in catalog_items
        if item.catalog_status.strip().lower() in {STATUS_OK, "already_present"}
        and not item.local_pdb
    )
    permanent_uncovered = sum(
        1 for item in catalog_items
        if item.status == STATUS_NOT_FOUND and not item.local_pdb and item.in_universe
    )
    # El descargador solo reconoce "HTTP 404", no "HTTP Error 404". Esas filas
    # son permanentes para la auditoría pero siguen en su cola de pendientes.
    fetcher_gap_items = sorted(
        item.accession for item in catalog_items
        if item.status == STATUS_NOT_FOUND and not item.local_pdb
        and item.in_universe
        and not fetcher_is_permanent(audit.catalog_rows.get(item.accession, {}))
    )
    retryable_pending = sum(
        1 for item in catalog_items
        if item.retryable and not item.local_pdb and item.in_universe
    )
    local_only = sum(
        1 for item in audit.items.values()
        if not item.in_catalog and item.local_pdb
    )
    raw_error = sum(
        1 for row in audit.catalog_rows.values()
        if str(row.get("status") or "").strip().lower() == STATUS_ERROR
    )
    raw_error_items = [
        item for item in catalog_items
        if item.catalog_status.strip().lower() == STATUS_ERROR
    ]
    raw_error_rate_limited = sum(
        1 for item in raw_error_items if item.status == STATUS_RATE_LIMITED)
    raw_error_not_found = sum(
        1 for item in raw_error_items if item.status == STATUS_NOT_FOUND)
    raw_error_other = (len(raw_error_items) - raw_error_rate_limited
                       - raw_error_not_found)
    pdb_bytes = sum(path.stat().st_size for path in audit.local_pdb.values())
    pae_bytes = sum(path.stat().st_size for path in audit.local_pae.values())

    pending_total = len(audit.pending)
    never_attempted = sum(1 for acc in audit.pending
                          if acc not in audit.catalog_rows)
    accounted = covered_in_universe + permanent_uncovered

    return {
        "universe_total": n_universe,
        "catalog_rows": len(audit.catalog_rows),
        "catalog_ok": catalog_counts[STATUS_OK],
        "catalog_rate_limited": catalog_counts[STATUS_RATE_LIMITED],
        "catalog_not_found": catalog_counts[STATUS_NOT_FOUND],
        "catalog_error": catalog_counts[STATUS_ERROR],
        "catalog_ok_with_file": catalog_ok_with_file,
        "catalog_ok_missing_file": catalog_ok_missing_file,
        "catalog_raw_error": raw_error,
        "catalog_raw_error_rate_limited": raw_error_rate_limited,
        "catalog_raw_error_not_found": raw_error_not_found,
        "catalog_raw_error_other": raw_error_other,
        "local_pdb": local_pdb,
        "local_pae": local_pae,
        "local_pdb_bytes": pdb_bytes,
        "local_pae_bytes": pae_bytes,
        "local_only": local_only,
        "pae_in_local": pae_in_local,
        "pae_missing_for_local": len(local_pdb_set - local_pae_set),
        "orphan_pae": len(local_pae_set - local_pdb_set),
        "covered_in_universe": covered_in_universe,
        "permanent_uncovered": permanent_uncovered,
        "fetcher_gap": len(fetcher_gap_items),
        "fetcher_gap_items": fetcher_gap_items,
        "fetcher_pending_total": pending_total + len(fetcher_gap_items),
        "pending_total": pending_total,
        "never_attempted": never_attempted,
        "retryable_pending": retryable_pending,
        "catalog_coverage": catalog_ok_with_file / denominator,
        "disk_coverage": covered_in_universe / denominator,
        "accounted_coverage": accounted / denominator,
        "pae_coverage": pae_in_local / max(1, local_pdb),
        "observed_429_rate": (catalog_counts[STATUS_RATE_LIMITED]
                              / max(1, len(audit.catalog_rows))),
    }


# ---------------------------------------------------------------------------
# Plan de lotes y ETA
# ---------------------------------------------------------------------------
def estimate_seconds(limit: int, workers: int, delay: float,
                     requests_per_accession: int = REQUESTS_PER_ACCESSION,
                     latency: float = LATENCY_PER_REQUEST,
                     retry_factor: float = 1.0) -> float:
    """ETA aproximada de un lote, con ``workers`` hilos y ``delay`` por request."""
    if limit <= 0 or workers <= 0:
        return 0.0
    waves = math.ceil(limit / workers)
    return (waves * requests_per_accession
            * (max(0.0, delay) + max(0.0, latency))
            * max(1.0, retry_factor))


def estimate_full_pending(pending: int, limit: int, workers: int, delay: float,
                          retry_factor: float = 1.0,
                          pause: float = FULL_DEBT_PAUSE) -> float:
    """ETA de cerrar ``pending`` con lotes de ``limit`` (sin PAE)."""
    if pending <= 0 or limit <= 0 or workers <= 0:
        return 0.0
    batch_size = min(limit, pending)
    batches = math.ceil(pending / batch_size)
    return batches * (estimate_seconds(batch_size, workers, delay,
                                       retry_factor=retry_factor) + pause)


def build_plan(summary: Mapping[str, Any],
               probe_limit: int = 10) -> dict[str, Any]:
    """Plan adaptativo de lotes; no descarga nada.

    Con 429 observados se empieza con una sonda de ``probe_limit`` accesiones
    (10 por defecto) a 1 worker y se sube solo si el lote previo queda limpio.
    Sin 429 se puede sostener el ritmo de producción desde el inicio.
    """
    pending = int(summary.get("pending_total", 0))
    # El descargador puede arrastrar 404 que aún no reconoce como permanentes.
    queue = int(summary.get("fetcher_pending_total", pending))
    rate_limited = int(summary.get("catalog_rate_limited", 0))
    retry_factor = 1.0 + min(0.5, float(summary.get("observed_429_rate", 0.0)))
    probe_limit = max(1, int(probe_limit))

    if pending <= 0:
        return {
            "steps": [], "retry_factor": retry_factor,
            "sustained_seconds": 0.0, "accelerated_seconds": 0.0,
            "slow_seconds": 0.0, "sustained_step": None,
        }

    if rate_limited > 0:
        specs = [
            ("sonda", 1, 2.0, probe_limit, 3,
             "Siempre. Mide si la API se recuperó de los 429."),
            ("rampa", 1, 1.5, 25, 3,
             "Solo si la sonda deja <= 1/10 rate_limited."),
            ("estable", 2, 1.0, 50, 3,
             "Ritmo sostenido si el 429 baja a <= 5%."),
            ("produccion", 2, 0.75, 100, 3,
             "Solo si el 429 es <= 2% en el paso anterior."),
        ]
    else:
        specs = [
            ("produccion", 2, 0.5, 100, 3,
             "Sin 429 observados en el catálogo."),
        ]

    steps: list[dict[str, Any]] = []
    for name, workers, delay, limit, retries, condition in specs:
        effective_limit = min(limit, queue)
        steps.append({
            "name": name,
            "workers": workers,
            "delay": delay,
            "limit": effective_limit,
            "retries": retries,
            "pae": False,
            "condition": condition,
            "eta_seconds": estimate_seconds(effective_limit, workers, delay,
                                            retry_factor=retry_factor),
            "command": (
                "python scripts/structure/fetch_afdb_cache.py"
                f" --limit {effective_limit} --workers {workers}"
                f" --delay {delay:g} --retries {retries}"
            ),
        })

    by_name = {step["name"]: step for step in steps}
    sustained = by_name.get("estable") or steps[-1]
    accelerated = by_name.get("produccion") or sustained
    slow = steps[0]

    return {
        "steps": steps,
        "retry_factor": retry_factor,
        "sustained_step": sustained,
        "sustained_seconds": estimate_full_pending(
            queue, sustained["limit"], sustained["workers"],
            sustained["delay"], retry_factor),
        "accelerated_seconds": estimate_full_pending(
            queue, accelerated["limit"], accelerated["workers"],
            accelerated["delay"], retry_factor),
        "slow_seconds": estimate_full_pending(
            queue, slow["limit"], slow["workers"], slow["delay"],
            retry_factor),
    }


def next_batch(audit: Audit, limit: int) -> list[str]:
    """Primeras accesiones pendientes en el orden del descargador."""
    if limit <= 0:
        return []
    return audit.pending[:limit]


def recommended_limit(audit: Audit, summary: Mapping[str, Any],
                      desired: int) -> int:
    """``--limit`` real para que el descargador incluya ``desired`` pendientes.

    El descargador puede tener por delante filas 404 que aún no reconoce como
    permanentes (``fetcher_gap``). Se suman las que caen antes del último
    pendiente deseado para que el lote haga los ``desired`` reintentables.
    """
    if desired <= 0 or not audit.pending:
        return 0
    desired = min(desired, len(audit.pending))
    last = audit.pending[desired - 1]
    gap_prefix = sum(1 for accession in summary.get("fetcher_gap_items", [])
                     if accession <= last)
    return desired + gap_prefix


# ---------------------------------------------------------------------------
# Artefactos
# ---------------------------------------------------------------------------
def _sorted_items(audit: Audit) -> list[AccessionStatus]:
    return sorted(audit.items.values(),
                  key=lambda item: (STATUS_ORDER.get(item.status, 99),
                                    item.accession))


def render_tsv(audit: Audit) -> str:
    """TSV por accesión (catálogo + caché local), reintentables primero."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=TSV_COLS, delimiter="\t",
                            extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for item in _sorted_items(audit):
        writer.writerow({
            "accession": item.accession,
            "status": item.status,
            "catalog_status": item.catalog_status,
            "local_pdb": _as_bool(bool(item.local_pdb)),
            "local_pae": _as_bool(bool(item.local_pae)),
            "pending": _as_bool(item.pending),
            "retryable": _as_bool(item.retryable),
            "permanent": _as_bool(item.permanent),
            "in_catalog": _as_bool(item.in_catalog),
            "in_universe": _as_bool(item.in_universe),
            "http_code": "" if item.http_code is None else item.http_code,
            "pdb_path": _rel(item.local_pdb) if item.local_pdb else "",
            "pae_path": _rel(item.local_pae) if item.local_pae else "",
            "error": item.error,
        })
    return buffer.getvalue()


def fmt_int(value: int) -> str:
    """Entero con separador de miles (punto, convención del proyecto)."""
    return f"{int(value):,}".replace(",", ".")


def fmt_pct(ratio: float, digits: int = 4) -> str:
    return f"{ratio * 100:.{digits}f}%"


def fmt_duration(seconds: float) -> str:
    if seconds <= 0:
        return "0 s"
    if seconds < 90:
        return f"{seconds:.0f} s"
    minutes = seconds / 60.0
    if minutes < 90:
        return f"{minutes:.1f} min"
    hours = minutes / 60.0
    if hours < 48:
        return f"{hours:.1f} h"
    return f"{hours / 24.0:.1f} días"


def _plan_rows(plan: Mapping[str, Any]) -> list[str]:
    lines = [
        "| paso | workers | delay (s) | limit | PAE | ETA del lote | condición |",
        "|---|---:|---:|---:|:---:|---:|---|",
    ]
    for step in plan["steps"]:
        lines.append(
            f"| `{step['name']}` | {step['workers']} | {step['delay']:g} | "
            f"{step['limit']} | {'sí' if step['pae'] else 'no'} | "
            f"{fmt_duration(step['eta_seconds'])} | {step['condition']} |"
        )
    return lines


def render_markdown(audit: Audit, summary: Mapping[str, Any],
                    plan: Mapping[str, Any], next_accs: list[str],
                    generated: str, tsv_path: Path | None = None) -> str:
    """Informe legible con métricas, cobertura, plan y siguiente lote."""
    n = summary["universe_total"]
    catalog = summary["catalog_rows"]
    catalog_denominator = max(1, catalog)
    catalog_total_pct = "100.00%" if catalog else "0.00%"
    lines: list[str] = []
    add = lines.append

    add("# Estado del catálogo AFDB/PAE")
    add("")
    add(f"Generado: `{generated}`")
    add("")
    add("Script: `scripts/structure/afdb_status.py` · modo solo lectura.")
    add(f"Catálogo: `{_rel(audit.catalog_path)}` · "
        f"caché: `{_rel(audit.targets_path)}`.")
    add("")
    add("Este informe **no descarga**, **no consulta la red** y **no modifica** "
        "el catálogo ni los archivos de la caché. Solo lee rutas, tamaños y "
        "estados locales.")
    add("")
    add("## Resumen")
    add("")
    add("| métrica | valor |")
    add("|---|---:|")
    add(f"| universo de accesiones | {fmt_int(n)} |")
    add(f"| filas del catálogo | {fmt_int(catalog)} |")
    add(f"| filas `ok` | {fmt_int(summary['catalog_ok'])} |")
    add(f"| filas `rate_limited` (429) | "
        f"{fmt_int(summary['catalog_rate_limited'])} |")
    add(f"| filas `not_found` (404) | "
        f"{fmt_int(summary['catalog_not_found'])} |")
    add(f"| filas `error` (otros) | {fmt_int(summary['catalog_error'])} |")
    add(f"| PDB locales | {fmt_int(summary['local_pdb'])} |")
    add(f"| PAE locales | {fmt_int(summary['local_pae'])} |")
    add(f"| PDB locales sin fila de catálogo | "
        f"{fmt_int(summary['local_only'])} |")
    add(f"| PDB locales sin PAE | "
        f"{fmt_int(summary['pae_missing_for_local'])} |")
    add(f"| pendientes de PDB | {fmt_int(summary['pending_total'])} |")
    if summary["fetcher_gap"]:
        add(f"| pendientes que ve el descargador | "
            f"{fmt_int(summary['fetcher_pending_total'])} |")
        add(f"| 404 que el descargador aún reintenta | "
            f"{fmt_int(summary['fetcher_gap'])} |")
    add(f"| cobertura PDB sobre el universo | "
        f"{fmt_pct(summary['disk_coverage'])} |")
    add(f"| cobertura contable (PDB + 404) | "
        f"{fmt_pct(summary['accounted_coverage'])} |")
    add("")
    add("## Inventario del catálogo")
    add("")
    add("| estado | filas | % del catálogo | ¿reintentable? | "
        "¿permanente? | criterio |")
    add("|---|---:|---:|:---:|:---:|---|")
    add(f"| `ok` | {fmt_int(summary['catalog_ok'])} | "
        f"{fmt_pct(summary['catalog_ok'] / catalog_denominator, 2)} | no | no | "
        "fila `ok`/`already_present` con archivo local |")
    add(f"| `rate_limited` | {fmt_int(summary['catalog_rate_limited'])} | "
        f"{fmt_pct(summary['catalog_rate_limited'] / catalog_denominator, 2)} | "
        "sí | no | 429 o `Retry-After` agotado |")
    add(f"| `not_found` | {fmt_int(summary['catalog_not_found'])} | "
        f"{fmt_pct(summary['catalog_not_found'] / catalog_denominator, 2)} | "
        "no | sí | 404; no se vuelve a pedir |")
    add(f"| `error` | {fmt_int(summary['catalog_error'])} | "
        f"{fmt_pct(summary['catalog_error'] / catalog_denominator, 2)} | sí | no | "
        "5xx, timeout, red u otro 4xx; se reintenta en otro lote |")
    add(f"| **total** | **{fmt_int(catalog)}** | **{catalog_total_pct}** | | | |")
    add("")
    if summary["catalog_raw_error"]:
        detail = (
            f"**{fmt_int(summary['catalog_raw_error_rate_limited'])}** son 429 "
            "(reintentables) y "
            f"**{fmt_int(summary['catalog_raw_error_not_found'])}** son 404 "
            "(permanentes)"
        )
        if summary["catalog_raw_error_other"]:
            detail += (f"; **{fmt_int(summary['catalog_raw_error_other'])}** "
                       "quedan como error reintentable")
        add(f"> El catálogo guarda **{fmt_int(summary['catalog_raw_error'])}** "
            "filas con `status=error` en bruto. La auditoría las reclasifica "
            f"por el código HTTP del campo `error`: {detail}. El catálogo "
            "**no** se reescribe; la clasificación es efectiva para el plan.")
        add("")
    if summary["catalog_ok_missing_file"]:
        add(f"> **Integridad:** {fmt_int(summary['catalog_ok_missing_file'])} "
            "fila(s) `ok` no tienen archivo local; se cuentan como `error` "
            "reintentable.")
        add("")
    add("## Cobertura sobre el universo")
    add("")
    add("| medida | n | % del universo |")
    add("|---|---:|---:|")
    add(f"| modelos PDB presentes en disco | "
        f"{fmt_int(summary['covered_in_universe'])} | "
        f"{fmt_pct(summary['disk_coverage'])} |")
    add(f"| filas de catálogo `ok` con archivo | "
        f"{fmt_int(summary['catalog_ok_with_file'])} | "
        f"{fmt_pct(summary['catalog_coverage'])} |")
    add(f"| PDB locales + 404 permanentes | "
        f"{fmt_int(summary['covered_in_universe'] + summary['permanent_uncovered'])}"
        f" | {fmt_pct(summary['accounted_coverage'])} |")
    add(f"| pendientes de PDB | {fmt_int(summary['pending_total'])} | "
        f"{fmt_pct(summary['pending_total'] / max(1, n))} |")
    add("")
    add("## Pendientes")
    add("")
    add("| concepto | n |")
    add("|---|---:|")
    add(f"| nunca intentadas (sin fila de catálogo y sin PDB) | "
        f"{fmt_int(summary['never_attempted'])} |")
    add(f"| reintentables con fila (429/error, sin PDB) | "
        f"{fmt_int(summary['retryable_pending'])} |")
    add(f"| permanentes 404 (no se reintentan) | "
        f"{fmt_int(summary['permanent_uncovered'])} |")
    if summary["fetcher_gap"]:
        add(f"| 404 que el descargador aún cuenta como pendientes | "
            f"{fmt_int(summary['fetcher_gap'])} |")
    add(f"| **total pendiente de PDB** | "
        f"**{fmt_int(summary['pending_total'])}** |")
    add("")
    if summary["fetcher_gap"]:
        add("## Discrepancia con el descargador")
        add("")
        add(f"`fetch_afdb_cache.is_permanent()` reconoce `\"HTTP 404\"` pero no "
            f"`\"HTTP Error 404\"`. Por eso **{fmt_int(summary['fetcher_gap'])}** "
            "filas que esta auditoría marca `not_found` siguen en la cola del "
            f"descargador (él ve {fmt_int(summary['fetcher_pending_total'])} "
            f"pendientes, no {fmt_int(summary['pending_total'])}). Son:")
        add("")
        add("| accesión |")
        add("|---|")
        for accession in summary["fetcher_gap_items"]:
            add(f"| `{accession}` |")
        add("")
        add("Se cerrarán como `not_found` cuando el descargador las pida (una "
            "petición cada una). Para no re-pedirlas, la corrección de fondo es "
            "ampliar `is_permanent()` a `\"HTTP Error 404\"`; este cambio no "
            "modifica el descargador.")
        add("")
    add("## Caché local")
    add("")
    add("| archivo | n | bytes | cobertura |")
    add("|---|---:|---:|---:|")
    add(f"| PDB | {fmt_int(summary['local_pdb'])} | "
        f"{fmt_int(summary['local_pdb_bytes'])} | "
        f"{fmt_pct(summary['disk_coverage'])} del universo |")
    add(f"| PAE | {fmt_int(summary['local_pae'])} | "
        f"{fmt_int(summary['local_pae_bytes'])} | "
        f"{fmt_pct(summary['pae_coverage'])} de los PDB locales |")
    add(f"| PAE huérfanos (sin PDB) | {fmt_int(summary['orphan_pae'])} | | |")
    add("")
    add("## Plan de lotes recomendado")
    add("")
    add("Sin descargar nada: solo se calcula el próximo ritmo. Supuestos del ETA: "
        f"{REQUESTS_PER_ACCESSION} solicitudes por accesión nueva (metadatos + "
        f"PDB), {LATENCY_PER_REQUEST:g} s de latencia media por solicitud, "
        f"factor de reintento {plan['retry_factor']:.4f} por el "
        f"{fmt_pct(summary['observed_429_rate'], 2)} de 429 observado, y "
        f"{FULL_DEBT_PAUSE:g} s de pausa entre lotes. No incluye esperas largas "
        "de `Retry-After`: el ETA es una estimación optimista.")
    add("")
    lines.extend(_plan_rows(plan))
    add("")
    if plan["steps"]:
        sustained = plan["sustained_step"]
        add(f"- Ritmo sostenido (`{sustained['name']}`, "
            f"{sustained['workers']} workers, delay {sustained['delay']:g}, "
            f"limit {sustained['limit']}): "
            f"**{fmt_duration(plan['sustained_seconds'])}** para los "
            f"{fmt_int(summary['pending_total'])} PDB pendientes "
            f"(~{math.ceil(summary['pending_total'] / max(1, sustained['limit']))} "
            "lotes).")
        if plan["accelerated_seconds"] < plan["sustained_seconds"]:
            add(f"- Si el 429 cae a <= 2%, el paso `produccion` baja el total a "
                f"**{fmt_duration(plan['accelerated_seconds'])}**.")
        if plan["slow_seconds"] > plan["sustained_seconds"]:
            add(f"- Si el 429 persiste, quedarse en `{plan['steps'][0]['name']}`: "
                f"**{fmt_duration(plan['slow_seconds'])}**. A ese ritmo la API "
                "no admite cerrar la deuda completa; hay que esperar.")
        add(f"- PAE de los {fmt_int(summary['pae_missing_for_local'])} PDB "
            "locales que no lo tienen: hacerlo **después** del PDB y solo si "
            "esos modelos se van a evaluar; con el mismo ritmo son ~"
            f"{fmt_duration(estimate_full_pending(summary['pae_missing_for_local'], 50, 2, 1.0, plan['retry_factor'], pause=BATCH_PAUSE))}.")
    add("")
    add("## Siguiente lote recomendado")
    add("")
    if next_accs:
        first = plan["steps"][0] if plan["steps"] else None
        if first:
            add("```bash")
            add(first["command"])
            add("```")
            add("")
        if summary["fetcher_gap"]:
            gap_prefix = [accession for accession in summary["fetcher_gap_items"]
                          if accession <= next_accs[-1]]
            if gap_prefix:
                names = ", ".join(f"`{accession}`" for accession in gap_prefix)
                add(f"> El comando usa `--limit {first['limit']}`: además de "
                    f"los {len(next_accs)} reintentables de la tabla, el "
                    f"descargador procesará {len(gap_prefix)} fila(s) 404 que "
                    f"aún no reconoce como permanentes ({names}). Se cerrarán "
                    "como `not_found` en esta pasada.")
                add("")
        add("Accesiones en el orden natural del descargador (alfabético):")
        add("")
        add("| # | accesión |")
        add("|---:|---|")
        for i, accession in enumerate(next_accs, start=1):
            add(f"| {i} | `{accession}` |")
        add("")
        add("Hoy las primeras pendientes son las que quedaron en 429; el "
            "descargador las reintenta porque solo los 404 son permanentes.")
    else:
        add("No hay accesiones pendientes.")
    add("")
    add("## Criterio de avance")
    add("")
    add("- Después del lote, ejecutar "
        "`python scripts/structure/afdb_status.py --dry-run` y contar los "
        "nuevos `ok`/`rate_limited`.")
    add("- Si >= 9/10 quedan `ok`, pasar al paso `rampa`.")
    add("- Si <= 8/10, mantener 1 worker / 2,0 s y repetir 10; **no acelerar**.")
    add("- No usar `--pae` durante la recuperación de PDB: duplica solicitudes.")
    add("")
    add("## Artefacto TSV")
    add("")
    tsv_label = _rel(tsv_path) if tsv_path is not None else "reports/afdb_status.tsv"
    add(f"`{tsv_label}` tiene una fila por accesión con evidencia local o de "
        f"catálogo ({fmt_int(len(audit.items))} filas), ordenada por estado "
        "(reintentables primero). Columnas: "
        + ", ".join(f"`{column}`" for column in TSV_COLS) + ".")
    add("")
    add("Las accesiones nunca intentadas no se listan una por una: se cuentan "
        "en el resumen y el siguiente lote se toma del orden determinista del "
        "descargador.")
    add("")
    add("## Reproducir")
    add("")
    add("```bash")
    add(".venv/bin/python scripts/structure/afdb_status.py --dry-run")
    add(".venv/bin/python scripts/structure/afdb_status.py")
    add(".venv/bin/python -m pytest tests/test_afdb_status.py -q")
    add("```")
    add("")
    return "\n".join(lines)


def _write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def write_reports(audit: Audit, summary: Mapping[str, Any],
                  plan: Mapping[str, Any],
                  next_accs: list[str],
                  reports_dir: str | os.PathLike[str],
                  generated: str | None = None) -> tuple[Path, Path]:
    """Escribe el TSV y el MD; no toca catálogo ni caché."""
    generated = generated or datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")
    reports = Path(reports_dir)
    tsv_path = reports / "afdb_status.tsv"
    md_path = reports / "afdb_status.md"
    _write_text_atomic(tsv_path, render_tsv(audit))
    _write_text_atomic(md_path, render_markdown(audit, summary, plan,
                                                next_accs, generated,
                                                tsv_path=tsv_path))
    return tsv_path, md_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audita el catálogo AFDB/PAE sin red y sin modificar la caché.")
    parser.add_argument("--data-root", type=Path, default=ROOT / "data")
    parser.add_argument("--catalog", type=Path,
                        default=ROOT / "data" / "metabolism" / "targets"
                        / "afdb_catalog.tsv")
    parser.add_argument("--targets", type=Path,
                        default=ROOT / "data" / "metabolism" / "targets"
                        / "afdb25k")
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    parser.add_argument("--limit-next", type=int, default=10,
                        help="accesiones del siguiente lote (def. 10)")
    parser.add_argument("--dry-run", action="store_true",
                        help="imprimir el plan sin escribir reportes")
    args = parser.parse_args(argv)

    if args.limit_next < 0:
        parser.error("--limit-next debe ser >= 0")
    if not args.data_root.is_dir():
        print(f"error: no existe --data-root {args.data_root}", file=sys.stderr)
        return 2
    if not args.catalog.is_file():
        print(f"error: no existe --catalog {args.catalog}", file=sys.stderr)
        return 2

    universe = load_universe(args.data_root)
    if not universe:
        print(f"error: no hay accesiones en {args.data_root}/regulation/*/genes.tsv",
              file=sys.stderr)
        return 2

    catalog = load_catalog(args.catalog)
    audit = build_audit(universe, catalog, args.targets,
                        catalog_path=args.catalog)
    summary = summarize(audit)
    next_accs = next_batch(audit, args.limit_next)
    probe_limit = recommended_limit(audit, summary, args.limit_next)
    plan = build_plan(summary, probe_limit=probe_limit or args.limit_next)

    print(f"universo={summary['universe_total']} "
          f"catálogo={summary['catalog_rows']} "
          f"ok={summary['catalog_ok']} "
          f"rate_limited={summary['catalog_rate_limited']} "
          f"not_found={summary['catalog_not_found']} "
          f"error={summary['catalog_error']}")
    print(f"pdb_local={summary['local_pdb']} pae_local={summary['local_pae']} "
          f"pendientes={summary['pending_total']} "
          f"cobertura_pdb={fmt_pct(summary['disk_coverage'])}")
    if next_accs:
        print("siguiente lote: " + " ".join(next_accs))
    else:
        print("siguiente lote: (ninguno)")

    if args.dry_run:
        print("plan recomendado:")
        for step in plan["steps"]:
            print(f"  {step['name']}: workers={step['workers']} "
                  f"delay={step['delay']:g} limit={step['limit']} "
                  f"eta={fmt_duration(step['eta_seconds'])}")
        print(f"eta_deuda={fmt_duration(plan['sustained_seconds'])}")
        print("dry-run: no se escriben reportes ni se toca la caché")
        return 0

    tsv_path, md_path = write_reports(audit, summary, plan, next_accs,
                                      args.reports_dir)
    print(f"reportes: {_rel(tsv_path)} {_rel(md_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
