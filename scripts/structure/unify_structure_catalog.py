#!/usr/bin/env python3
"""Unifica los catálogos AFDB y ESMFold en un único catálogo con procedencia.

Entradas (por defecto):

- ``data/metabolism/targets/afdb_catalog.tsv``
- ``data/metabolism/targets/esmfold_catalog.tsv``
- ``data/metabolism/targets/afdb25k/*.pdb``
- ``data/metabolism/targets/esmfold_kaggle/*.pdb``

Salidas (por defecto):

- ``data/metabolism/targets/structure_catalog.tsv``
- ``reports/structure_catalog_report.tsv``
- ``reports/structure_catalog_report.md``

Reglas:

- una fila por accesión;
- si una accesión existe en AFDB y ESMFold, se prefiere AFDB cuando su
  archivo es utilizable; si AFDB falta o está en error, se conserva ESMFold;
- si dos accesiones distintas comparten el mismo SHA-256, se conserva una
  como canónica (AFDB > ESMFold, luego mayor pLDDT) y la otra queda con
  ``duplicate_of``;
- nunca se borran PDB ni filas de los catálogos de origen;
- ``--dry-run`` calcula todo pero no escribe;
- sin red, ordenado y determinista (idempotente).

Uso:
  python scripts/structure/unify_structure_catalog.py --dry-run
  python scripts/structure/unify_structure_catalog.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

try:  # ejecución como módulo (pytest, ``python -m``)
    from scripts.structure.ingest_kaggle_zips import (
        read_tsv_rows, write_tsv_rows,
    )
except ImportError:  # ejecución directa del script
    from ingest_kaggle_zips import read_tsv_rows, write_tsv_rows

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "data" / "metabolism" / "targets"

CATALOG_COLUMNS = [
    "accession", "source", "plddt", "ptm", "length", "path", "sha256",
    "bytes", "status", "duplicate_of", "updated",
]
REPORT_COLUMNS = [
    "kind", "accession", "source", "path", "sha256", "bytes", "related",
    "detail",
]
SOURCE_PRIORITY = {"AFDB": 0, "ESMFold": 1}
STATUS_PRIORITY = {"ok": 0, "missing": 1, "not_found": 2, "error": 3}


@dataclass
class Entry:
    accession: str
    source: str
    plddt: str = ""
    ptm: str = ""
    length: str = ""
    path: str = ""
    sha256: str = ""
    nbytes: int = 0
    status: str = ""
    duplicate_of: str = ""
    updated: str = ""
    abs_path: Path | None = None
    raw: dict = field(default_factory=dict)

    def usable(self) -> bool:
        return (self.status == "ok" and self.abs_path is not None
                and self.abs_path.exists())

    def as_row(self) -> dict:
        return {
            "accession": self.accession,
            "source": self.source,
            "plddt": self.plddt,
            "ptm": self.ptm,
            "length": self.length,
            "path": self.path,
            "sha256": self.sha256,
            "bytes": str(self.nbytes),
            "status": self.status,
            "duplicate_of": self.duplicate_of,
            "updated": self.updated,
        }


# ---------------------------------------------------------------------------
# Normalización y lectura de archivos
# ---------------------------------------------------------------------------


def _num(value) -> float | None:
    try:
        num = float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


def _as_int(value, default: int = 0) -> int:
    num = _num(value)
    return int(num) if num is not None else default


def _norm_plddt(value) -> str:
    num = _num(value)
    if num is None:
        return ""
    if num > 1.5:  # AFDB suele guardar pLDDT 0-100 en el B-factor
        num /= 100.0
    return f"{num:.3f}"


def _norm_ptm(value) -> str:
    num = _num(value)
    return "" if num is None else f"{num:.4f}"


def _norm_length(value) -> str:
    if value is None or not str(value).strip():
        return ""
    num = _num(value)
    return "" if num is None else str(int(num))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _pdb_stats(path: Path) -> tuple[str, str]:
    """(pLDDT medio de CA, nº de residuos) leídos del PDB, sin red.

    AFDB guarda pLDDT 0-100 en el B-factor; se normaliza a 0-1. Solo se
    considera el primer MODEL. Ante cualquier problema devuelve ("", "").
    """
    plddts: list[float] = []
    residues: set[tuple[str, str, str]] = set()
    model = 0
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                record = line[:6].strip()
                if record == "MODEL":
                    model += 1
                    continue
                if record in ("ENDMDL", "END"):
                    continue
                if record not in ("ATOM", "HETATM"):
                    continue
                if model > 1:
                    continue
                altloc = line[16:17]
                if altloc not in (" ", "", "A"):
                    continue
                chain = line[21:22].strip()
                resseq = line[22:26].strip()
                icode = line[26:27].strip()
                residues.add((chain, resseq, icode))
                if line[12:16].strip() != "CA":
                    continue
                try:
                    plddts.append(float(line[60:66]))
                except ValueError:
                    continue
    except OSError:
        return "", ""
    length = str(len(residues)) if residues else ""
    if not plddts:
        return "", length
    mean = sum(plddts) / len(plddts)
    if mean > 1.5:
        mean /= 100.0
    return f"{mean:.3f}", length


def _resolve_path(raw: str, catalog_path: Path) -> Path | None:
    text = (raw or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if path.is_absolute():
        return path
    # Una ruta relativa se interpreta desde el catálogo; si no existe, se
    # prueba desde la raíz del repo (catálogos antiguos con rutas tipo
    # ``data/metabolism/targets/...``).
    candidate = catalog_path.parent / path
    if candidate.exists():
        return candidate
    root_candidate = ROOT / path
    if root_candidate.exists():
        return root_candidate
    return candidate


def _display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except (OSError, ValueError):
        return str(path.resolve())


def _status(raw: str, path: Path | None) -> str:
    raw = (raw or "").strip().lower()
    exists = bool(path and path.exists())
    if raw in ("ok", "already_present", "done"):
        return "ok" if exists else "missing"
    if raw in ("not_found", "not found", "404"):
        return "not_found"
    if raw == "missing":
        return "missing"
    if not raw:
        return "ok" if exists else "error"
    return "error"


def _acc(row: dict) -> str:
    return (row.get("accession") or row.get("uniprot") or "").strip()


# ---------------------------------------------------------------------------
# Construcción de entradas por fuente
# ---------------------------------------------------------------------------


def _afdb_entries(rows: list[dict], catalog_path: Path) -> list[Entry]:
    entries: list[Entry] = []
    for row in rows:
        acc = _acc(row)
        if not acc:
            continue
        raw_path = (row.get("pdb") or "").strip()
        abs_path = _resolve_path(raw_path, catalog_path)
        status = _status(row.get("status", ""), abs_path)
        plddt = _norm_plddt(row.get("plddt"))
        ptm = _norm_ptm(row.get("ptm"))
        length = _norm_length(row.get("length"))
        calc_plddt = calc_length = ""
        if abs_path is not None and abs_path.exists():
            calc_plddt, calc_length = _pdb_stats(abs_path)
        plddt = plddt or calc_plddt
        length = length or calc_length
        sha = (row.get("sha256") or "").strip()
        nbytes = _as_int(row.get("pdb_bytes") or row.get("bytes"), 0)
        if abs_path is not None and abs_path.exists():
            if not nbytes:
                nbytes = abs_path.stat().st_size
            if not sha:
                sha = _sha256(abs_path)
        updated = (row.get("updated") or row.get("ingested_at")
                   or "").strip()
        entries.append(Entry(
            accession=acc, source="AFDB", plddt=plddt, ptm=ptm,
            length=length,
            path=_display_path(abs_path) if abs_path else "",
            sha256=sha, nbytes=nbytes, status=status,
            duplicate_of="", updated=updated, abs_path=abs_path, raw=row,
        ))
    return entries


def _esm_entries(rows: list[dict], esmfold_dir: Path) -> list[Entry]:
    entries: list[Entry] = []
    for row in rows:
        acc = _acc(row)
        if not acc:
            continue
        abs_path = esmfold_dir / f"{acc}.pdb"
        status = _status(row.get("status", ""), abs_path)
        plddt = _norm_plddt(row.get("plddt"))
        ptm = _norm_ptm(row.get("ptm"))
        length = _norm_length(row.get("length"))
        sha = (row.get("sha256") or "").strip()
        nbytes = _as_int(row.get("bytes") or row.get("pdb_bytes"), 0)
        if abs_path.exists():
            if not nbytes:
                nbytes = abs_path.stat().st_size
            if not sha:
                sha = _sha256(abs_path)
            if not length:
                _, length = _pdb_stats(abs_path)
        updated = (row.get("ingested_at") or row.get("updated")
                   or "").strip()
        entries.append(Entry(
            accession=acc, source="ESMFold", plddt=plddt, ptm=ptm,
            length=length, path=_display_path(abs_path), sha256=sha,
            nbytes=nbytes, status=status, duplicate_of="", updated=updated,
            abs_path=abs_path, raw=row,
        ))
    return entries


# ---------------------------------------------------------------------------
# Deduplicación
# ---------------------------------------------------------------------------


def _plddt_num(entry: Entry) -> float:
    num = _num(entry.plddt)
    return -1.0 if num is None else num


def _best(entries: list[Entry]) -> Entry:
    """Mejor fila dentro de una misma fuente (determinista)."""
    return max(entries, key=lambda e: (
        1 if e.usable() else 0,
        -STATUS_PRIORITY.get(e.status, 9),
        _plddt_num(e),
        e.updated,
        e.sha256,
        e.path,
    ))


def _canonical_hash_key(entry: Entry) -> tuple:
    return (
        0 if entry.usable() else 1,
        STATUS_PRIORITY.get(entry.status, 9),
        SOURCE_PRIORITY.get(entry.source, 9),
        -_plddt_num(entry),
        entry.accession,
        entry.path,
    )


def _entry_report(kind: str, entry: Entry, related: str = "",
                  detail: str = "") -> dict:
    return {
        "kind": kind,
        "accession": entry.accession,
        "source": entry.source,
        "path": entry.path,
        "sha256": entry.sha256,
        "bytes": str(entry.nbytes),
        "related": related,
        "detail": detail,
    }


def _file_report(kind: str, accession: str, source: str, path: str,
                 sha256: str, nbytes: int, detail: str) -> dict:
    return {
        "kind": kind,
        "accession": accession,
        "source": source,
        "path": path,
        "sha256": sha256,
        "bytes": str(nbytes),
        "related": "",
        "detail": detail,
    }


def _group(entries: list[Entry]) -> dict[str, list[Entry]]:
    grouped: dict[str, list[Entry]] = defaultdict(list)
    for entry in entries:
        grouped[entry.accession].append(entry)
    return dict(grouped)


# ---------------------------------------------------------------------------
# Huérfanos
# ---------------------------------------------------------------------------


def _scan_orphans(directory: Path, source: str,
                  by_acc: dict[str, list[Entry]],
                  hash_files: bool = True) -> list[dict]:
    """PDB en disco que el catálogo de esa fuente no registra (o registra
    con otra ruta). No modifica nada."""
    if not directory or not Path(directory).is_dir():
        return []
    directory = Path(directory)
    entries = [e for group in by_acc.values() for e in group]
    registered = {e.abs_path.resolve() for e in entries if e.abs_path}
    registered_paths: dict[str, list[Path]] = defaultdict(list)
    for entry in entries:
        if entry.abs_path:
            registered_paths[entry.accession].append(entry.abs_path)
    out: list[dict] = []
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() != ".pdb":
            continue
        acc = path.stem
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path
        if acc not in by_acc:
            kind = "orphan_pdb"
            detail = f"accesión {acc} no registrada en el catálogo {source}"
        elif not registered_paths.get(acc):
            kind = "orphan_pdb"
            estados = ", ".join(sorted({e.status for e in by_acc[acc]}))
            detail = (f"accesión {acc} registrada sin ruta en {source} "
                      f"(estado: {estados or 'desconocido'})")
        elif resolved not in registered:
            kind = "path_mismatch"
            expected = ", ".join(_display_path(p) for p in
                                 registered_paths.get(acc, []))
            detail = f"el catálogo {source} apunta a {expected}"
        else:
            continue
        try:
            sha = _sha256(path) if hash_files else ""
            size = path.stat().st_size
        except OSError:
            sha, size = "", 0
        out.append(_file_report(
            kind, acc, source, _display_path(path), sha,
            size, detail,
        ))
    return out


# ---------------------------------------------------------------------------
# Reporte Markdown
# ---------------------------------------------------------------------------


def _md_cell(value) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _markdown(summary: dict, report_rows: list[dict], limit: int) -> str:
    lines = [
        "# Catálogo unificado de estructuras AFDB + ESMFold",
        "",
        "Generado por `scripts/structure/unify_structure_catalog.py`.",
        "Este reporte no es evidencia experimental: los modelos son "
        "predicciones computacionales.",
        "",
        "## Resumen",
        "",
        "| métrica | valor |",
        "|---|---:|",
    ]
    labels = [
        ("Filas AFDB de entrada", "afdb_rows"),
        ("Filas ESMFold de entrada", "esmfold_rows"),
        ("Accesiones unificadas", "unified"),
        ("  · fuente AFDB", "source_AFDB"),
        ("  · fuente ESMFold", "source_ESMFold"),
        ("  · estado ok", "status_ok"),
        ("  · estado missing", "status_missing"),
        ("  · estado error", "status_error"),
        ("  · estado not_found", "status_not_found"),
        ("Duplicados por accesión", "duplicate_accessions"),
        ("Duplicados por hash", "hash_duplicates"),
        ("PDB huérfanos", "orphans"),
        ("Rutas no coincidentes", "path_mismatches"),
        ("Archivos faltantes", "missing_files"),
    ]
    for label, key in labels:
        lines.append(f"| {label} | {summary.get(key, 0)} |")
    lines += ["", "## Reglas aplicadas", "",
              "- Una fila por accesión.",
              "- Misma accesión en AFDB y ESMFold: se prefiere AFDB si su "
              "archivo es utilizable; si no, se conserva ESMFold.",
              "- Mismo SHA-256 en accesiones distintas: canónica = archivo "
              "utilizable > estado > AFDB > mayor pLDDT > accesión menor.",
              "- Los duplicados no se borran: quedan en el catálogo con "
              "`duplicate_of` y en este reporte.",
              ""]
    sections = [
        ("Duplicados por accesión", "duplicate_accession"),
        ("Duplicados por hash", "duplicate_hash"),
        ("PDB huérfanos", "orphan_pdb"),
        ("Rutas no coincidentes", "path_mismatch"),
        ("Archivos faltantes", "missing_file"),
    ]
    for title, kind in sections:
        subset = sorted(
            (r for r in report_rows if r["kind"] == kind),
            key=lambda r: (r["accession"], r["source"], r["path"]),
        )
        lines += [f"## {title}", ""]
        if not subset:
            lines += ["Sin hallazgos.", ""]
            continue
        lines += ["| accession | source | path | bytes | related | detail |",
                  "|---|---|---|---:|---|---|"]
        for row in subset[:limit]:
            lines.append(
                "| {accession} | {source} | {path} | {bytes} | "
                "{related} | {detail} |".format(
                    accession=_md_cell(row["accession"]),
                    source=_md_cell(row["source"]),
                    path=_md_cell(row["path"]),
                    bytes=_md_cell(row["bytes"]),
                    related=_md_cell(row["related"]),
                    detail=_md_cell(row["detail"]),
                )
            )
        if len(subset) > limit:
            lines += ["", f"_{len(subset) - limit} filas adicionales "
                          f"en el TSV._"]
        lines.append("")
    lines += [
        "## Límites",
        "",
        "- `plddt` se normaliza a 0-1. En AFDB se calcula como la media de "
        "B-factor de los CA del primer modelo; es una aproximación, no el "
        "valor oficial por residuo.",
        "- `ptm` solo existe en el catálogo ESMFold; en AFDB queda vacío.",
        "- `length` en AFDB se cuenta desde los residuos del PDB, no desde "
        "la secuencia UniProt.",
        "- La deduplicación es por accesión y por SHA-256 del archivo PDB; "
        "no detecta la misma proteína con coordenadas distintas ni "
        "duplicados por secuencia.",
        "- Solo se inspeccionan `*.pdb` en los directorios configurados; no "
        "hay mmCIF, ensamblados biológicos ni validación geométrica.",
        "- No se borra ni se mueve ningún archivo, tampoco los huérfanos.",
        "- El script no usa red.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Núcleo
# ---------------------------------------------------------------------------


def unify(afdb_catalog: Path, esmfold_catalog: Path, *,
          afdb_dir: Path | None = None, esmfold_dir: Path | None = None,
          out: Path | None = None, report_tsv: Path | None = None,
          report_md: Path | None = None, dry_run: bool = False,
          scan_orphans: bool = True, hash_orphans: bool = True,
          report_limit: int = 200,
          allow_empty: bool = False) -> dict:
    """Construye el catálogo unificado y, salvo ``dry_run``, lo escribe."""
    afdb_catalog = Path(afdb_catalog)
    esmfold_catalog = Path(esmfold_catalog)
    afdb_dir = Path(afdb_dir) if afdb_dir else (
        afdb_catalog.parent / "afdb25k")
    esmfold_dir = Path(esmfold_dir) if esmfold_dir else (
        esmfold_catalog.parent / "esmfold_kaggle")

    afdb_all = _afdb_entries(read_tsv_rows(afdb_catalog), afdb_catalog)
    esm_all = _esm_entries(read_tsv_rows(esmfold_catalog), esmfold_dir)
    afdb_by = _group(afdb_all)
    esm_by = _group(esm_all)

    report_rows: list[dict] = []

    # Duplicados dentro de una misma fuente.
    for source, grouped in (("AFDB", afdb_by), ("ESMFold", esm_by)):
        for acc in sorted(grouped):
            entries = grouped[acc]
            if len(entries) < 2:
                continue
            best = _best(entries)
            for entry in entries:
                if entry is best:
                    continue
                report_rows.append(_entry_report(
                    "duplicate_accession", entry,
                    related=f"{best.source}:{best.path}",
                    detail="fila duplicada en el mismo catálogo; se "
                           "conserva la de mayor pLDDT",
                ))

    # Una fila por accesión: AFDB si es utilizable; si no, ESMFold.
    chosen: dict[str, Entry] = {}
    for acc in sorted(set(afdb_by) | set(esm_by)):
        afdb = _best(afdb_by[acc]) if acc in afdb_by else None
        esm = _best(esm_by[acc]) if acc in esm_by else None
        if afdb is not None and esm is not None:
            if afdb.usable():
                winner, loser = afdb, esm
                detail = "AFDB preferido para la misma accesión"
            elif esm.usable():
                winner, loser = esm, afdb
                detail = "AFDB no utilizable; se conserva ESMFold"
            else:
                winner, loser = afdb, esm
                detail = "Ninguna copia utilizable; se conserva AFDB"
            report_rows.append(_entry_report(
                "duplicate_accession", loser,
                related=f"{winner.source}:{winner.path}", detail=detail,
            ))
            chosen[acc] = winner
        else:
            chosen[acc] = afdb if afdb is not None else esm

    # Duplicados por hash entre accesiones distintas.
    by_hash: dict[str, list[Entry]] = defaultdict(list)
    for entry in chosen.values():
        if entry.sha256:
            by_hash[entry.sha256].append(entry)
    hash_duplicates = 0
    for sha in sorted(by_hash):
        group = by_hash[sha]
        if len(group) < 2:
            continue
        canonical = min(group, key=_canonical_hash_key)
        for entry in group:
            if entry is canonical:
                continue
            entry.duplicate_of = canonical.accession
            hash_duplicates += 1
            report_rows.append(_entry_report(
                "duplicate_hash", entry, related=canonical.accession,
                detail=f"sha256 idéntico a {canonical.accession} "
                       f"({canonical.source})",
            ))

    # Archivos que el catálogo declara pero no existen en disco.
    missing_seen: set[tuple[str, str, str]] = set()
    missing_files = 0
    for entries in (afdb_all, esm_all):
        for entry in entries:
            if entry.abs_path is None or entry.abs_path.exists():
                continue
            key = (entry.accession, entry.source, entry.path)
            if key in missing_seen:
                continue
            missing_seen.add(key)
            missing_files += 1
            report_rows.append(_entry_report(
                "missing_file", entry,
                detail=f"estado del catálogo: {entry.status}",
            ))

    # PDB huérfanos o con ruta distinta a la registrada.
    orphans = 0
    path_mismatches = 0
    if scan_orphans:
        for directory, source, grouped in (
                (afdb_dir, "AFDB", afdb_by),
                (esmfold_dir, "ESMFold", esm_by)):
            for row in _scan_orphans(directory, source, grouped,
                                     hash_files=hash_orphans):
                if row["kind"] == "orphan_pdb":
                    orphans += 1
                else:
                    path_mismatches += 1
                report_rows.append(row)

    report_rows.sort(key=lambda r: (r["kind"], r["accession"],
                                    r["source"], r["path"], r["detail"]))
    source_counts = Counter(e.source for e in chosen.values())
    status_counts = Counter(e.status for e in chosen.values())
    duplicate_accessions = sum(
        1 for r in report_rows if r["kind"] == "duplicate_accession")
    summary = {
        "afdb_rows": len(afdb_all),
        "esmfold_rows": len(esm_all),
        "unified": len(chosen),
        "source_AFDB": source_counts.get("AFDB", 0),
        "source_ESMFold": source_counts.get("ESMFold", 0),
        "status_ok": status_counts.get("ok", 0),
        "status_missing": status_counts.get("missing", 0),
        "status_error": status_counts.get("error", 0),
        "status_not_found": status_counts.get("not_found", 0),
        "duplicate_accessions": duplicate_accessions,
        "hash_duplicates": hash_duplicates,
        "orphans": orphans,
        "path_mismatches": path_mismatches,
        "missing_files": missing_files,
        "report_rows": len(report_rows),
        "written": False,
    }

    catalog_rows = [chosen[acc].as_row() for acc in sorted(chosen)]
    wrote = False
    if not dry_run and (catalog_rows or allow_empty):
        if out is not None:
            write_tsv_rows(out, CATALOG_COLUMNS, catalog_rows)
            wrote = True
        if report_tsv is not None:
            write_tsv_rows(report_tsv, REPORT_COLUMNS, report_rows)
            wrote = True
        if report_md is not None:
            report_md = Path(report_md)
            report_md.parent.mkdir(parents=True, exist_ok=True)
            tmp = report_md.with_suffix(report_md.suffix + ".tmp")
            try:
                tmp.write_text(_markdown(summary, report_rows, report_limit),
                               encoding="utf-8")
                os.replace(tmp, report_md)
            finally:
                if tmp.exists():
                    tmp.unlink()
            wrote = True
    summary["written"] = wrote
    return {
        "summary": summary,
        "rows": {row["accession"]: row for row in catalog_rows},
        "catalog_rows": catalog_rows,
        "report_rows": report_rows,
        "entries": chosen,
    }


def _print_summary(summary: dict, dry_run: bool, out: Path,
                   report_tsv: Path, report_md: Path) -> None:
    prefix = "[dry-run] " if dry_run else ""
    print(
        f"{prefix}AFDB={summary['afdb_rows']} "
        f"ESMFold={summary['esmfold_rows']} "
        f"unificado={summary['unified']} "
        f"(AFDB={summary['source_AFDB']}, "
        f"ESMFold={summary['source_ESMFold']}) "
        f"ok={summary['status_ok']} missing={summary['status_missing']} "
        f"error={summary['status_error']} "
        f"not_found={summary['status_not_found']} "
        f"dups_acc={summary['duplicate_accessions']} "
        f"dups_hash={summary['hash_duplicates']} "
        f"huérfanos={summary['orphans']} "
        f"path_mismatch={summary['path_mismatches']} "
        f"faltantes={summary['missing_files']}"
    )
    if not dry_run:
        print(f"catálogo -> {out}")
        print(f"reporte TSV -> {report_tsv}")
        print(f"reporte MD -> {report_md}")
    else:
        print("dry-run: no se escribió ningún archivo")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--afdb-catalog", type=Path,
                    default=TARGETS / "afdb_catalog.tsv")
    ap.add_argument("--esmfold-catalog", type=Path,
                    default=TARGETS / "esmfold_catalog.tsv")
    ap.add_argument("--afdb-dir", type=Path,
                    default=TARGETS / "afdb25k")
    ap.add_argument("--esmfold-dir", type=Path,
                    default=TARGETS / "esmfold_kaggle")
    ap.add_argument("--out", type=Path,
                    default=TARGETS / "structure_catalog.tsv")
    ap.add_argument("--report-tsv", type=Path,
                    default=ROOT / "reports" / "structure_catalog_report.tsv")
    ap.add_argument("--report-md", type=Path,
                    default=ROOT / "reports" / "structure_catalog_report.md")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-orphans", action="store_true",
                    help="no escanear directorios de PDB")
    ap.add_argument("--no-hash-orphans", action="store_true",
                    help="no calcular SHA-256 de PDB huérfanos")
    ap.add_argument("--report-limit", type=int, default=200,
                    help="filas por sección en el reporte Markdown")
    ap.add_argument("--allow-empty", action="store_true",
                    help="escribir catálogo vacío si no hay entradas")
    ap.add_argument("--json", action="store_true",
                    help="imprimir el resumen como JSON")
    args = ap.parse_args(argv)
    if not args.afdb_catalog.exists() and not args.esmfold_catalog.exists():
        ap.error("no existe ninguno de los dos catálogos de entrada")
    result = unify(
        args.afdb_catalog, args.esmfold_catalog,
        afdb_dir=args.afdb_dir, esmfold_dir=args.esmfold_dir,
        out=args.out, report_tsv=args.report_tsv,
        report_md=args.report_md, dry_run=args.dry_run,
        scan_orphans=not args.no_orphans,
        hash_orphans=not args.no_hash_orphans,
        report_limit=args.report_limit, allow_empty=args.allow_empty,
    )
    summary = result["summary"]
    if args.json:
        print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    else:
        _print_summary(summary, args.dry_run, args.out, args.report_tsv,
                       args.report_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
