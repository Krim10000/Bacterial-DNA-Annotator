#!/usr/bin/env python3
"""Ingiere ZIP de ESMFold de data/Kaggler sin duplicar accesiones.

Cada ZIP contiene ``<lote>/metrics.tsv`` y ``<lote>/pdb/<acc>.pdb``.
El catálogo conserva pLDDT/pTM/longitud, origen y hash. La copia local queda
en ``data/metabolism/targets/esmfold_kaggle/<acc>.pdb``. Si una accesión
aparece en dos lotes, se conserva el PDB con mayor pLDDT y se reporta el
duplicado; nunca se sobrescribe con una copia de peor calidad.

Una carpeta de lote incompleta (sin ``metrics.tsv``, PDB truncado o PDB sin
fila de métricas) se reporta y **no se ingiere**. La CLI es estricta por
defecto; ``ingest(..., require_complete=True)`` expone lo mismo como API.

``--report`` audita, sin borrar nada:

- ZIP que aún no tienen filas en el catálogo;
- PDB en el directorio de salida que no tienen fila en el catálogo;
- accesiones repetidas entre lotes (se conserva la de mayor pLDDT);
- PDB truncados, ZIP incompletos y hashes repetidos.

``--prune-missing`` añade al reporte las filas del catálogo cuyo PDB local
falta. No elimina filas ni archivos: solo marca ``status=missing`` (salvo en
``--dry-run``) y conserva pLDDT/pTM/longitud/sha256/source_zip.

El lector/escritor TSV se exporta para que
``unify_structure_catalog.py`` reutilice exactamente el mismo formato.

Uso:
  python scripts/structure/ingest_kaggle_zips.py --dry-run
  python scripts/structure/ingest_kaggle_zips.py --dry-run --report
  python scripts/structure/ingest_kaggle_zips.py --dry-run --report --prune-missing
  python scripts/structure/ingest_kaggle_zips.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
import os
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COLS = ["accession", "status", "plddt", "ptm", "length", "source_zip",
        "source_member", "sha256", "bytes", "duplicate_of", "error",
        "ingested_at"]
CATALOG_COLUMNS = COLS  # alias explícito para otros módulos

REPORT_COLUMNS = [
    "kind", "accession", "source_zip", "path", "sha256", "bytes",
    "related", "detail",
]
REPORT_KINDS = [
    "new_zip", "incomplete_zip", "corrupt_pdb", "missing_metrics",
    "metrics_status", "duplicate_accession", "duplicate_hash",
    "orphan_pdb", "missing_file",
]
OK_METRIC_STATUSES = {"", "done", "ok", "already_present"}


def read_tsv_rows(path: Path) -> list[dict]:
    """Lee un TSV como lista de diccionarios (tolera BOM y filas vacías)."""
    path = Path(path)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            clean = {}
            for key, value in row.items():
                if key is None:
                    continue
                key = str(key).strip()
                if isinstance(value, str):
                    value = value.strip()
                clean[key] = value
            if any(v for v in clean.values()):
                rows.append(clean)
    return rows


def write_tsv_rows(path: Path, columns: list[str], rows) -> None:
    """Escribe un TSV de forma atómica; no toca el archivo si falla."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    try:
        with tmp.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({col: row.get(col, "") for col in columns})
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def _report(kind: str, *, accession: str = "", source_zip: str = "",
            path: str = "", sha256: str = "", nbytes: int | str = "",
            related: str = "", detail: str = "") -> dict:
    return {
        "kind": kind,
        "accession": accession,
        "source_zip": source_zip,
        "path": path,
        "sha256": sha256,
        "bytes": str(nbytes),
        "related": related,
        "detail": detail,
    }


def _catalog(path: Path) -> dict[str, dict]:
    return {r["accession"]: r for r in read_tsv_rows(path)
            if r.get("accession")}


def _write(path: Path, rows: dict[str, dict]) -> None:
    write_tsv_rows(path, COLS, [rows[acc] for acc in sorted(rows)])


def _num(value) -> float | None:
    try:
        num = float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


def _metrics(zf: zipfile.ZipFile, member: str) -> dict[str, dict]:
    try:
        text = zf.read(member).decode("utf-8-sig", errors="replace")
    except (KeyError, zipfile.BadZipFile, OSError, RuntimeError):
        return {}
    out = {}
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        acc = (r.get("uniprot") or r.get("accession") or r.get("acc")
               or r.get("protein_id") or "").strip()
        if acc:
            out[acc] = r
    return out


def _pdb_problem(data: bytes) -> str:
    """Devuelve el motivo por el que un PDB está truncado ('' si está bien)."""
    if not data:
        return "PDB vacío"
    if b"\x00" in data:
        return "PDB contiene bytes NUL"
    text = data.decode("utf-8", errors="replace")
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    if not lines:
        return "PDB sin contenido"
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    if not atoms:
        return "PDB sin registros ATOM/HETATM"
    incomplete = [line for line in atoms if len(line) < 54]
    if incomplete:
        return f"PDB truncado: {len(incomplete)} línea(s) ATOM incompletas"
    if not lines[-1].startswith("END"):
        return "PDB truncado: falta registro END final"
    return ""


def _member_base(member: str) -> str:
    if "/pdb/" in member:
        return member.rsplit("/pdb/", 1)[0]
    return member.rsplit("/", 1)[0] if "/" in member else ""


def _inspect_zip(zp: Path, require_complete: bool) -> tuple[
        list[dict], list[dict], dict]:
    """Lee y valida un ZIP completo. No escribe nada.

    Devuelve ``(candidatos, hallazgos, estado)``. En modo estricto un ZIP
    con cualquier problema de completitud se descarta por completo.
    """
    findings: list[dict] = []
    candidates: list[dict] = []
    info: dict = {
        "zip": zp.name, "metrics": False, "pdbs": 0, "skipped": 0,
        "incomplete": False, "error": "",
    }
    try:
        zf = zipfile.ZipFile(zp)
    except (zipfile.BadZipFile, OSError) as exc:
        info["error"] = str(exc)
        info["incomplete"] = True
        findings.append(_report("corrupt_zip", source_zip=zp.name,
                                path=str(zp),
                                detail=f"ZIP ilegible: {exc}"))
        return candidates, findings, info

    with zf:
        names = [n for n in zf.namelist()
                 if not n.startswith("__MACOSX/") and not n.endswith("/")]
        metric_files = [n for n in names
                        if n.endswith("/metrics.tsv") or n == "metrics.tsv"]
        info["metrics"] = bool(metric_files)
        metrics_by_base: dict[str, dict[str, dict]] = {}
        for mpath in metric_files:
            base = mpath.rsplit("/", 1)[0] if "/" in mpath else ""
            metrics_by_base[base] = _metrics(zf, mpath)

        pdb_members = [n for n in names if n.lower().endswith(".pdb")]
        info["pdbs"] = len(pdb_members)
        if not pdb_members:
            info["incomplete"] = True
            findings.append(_report(
                "incomplete_zip", source_zip=zp.name,
                detail="ZIP sin PDB; no se ingiere",
            ))
            return candidates, findings, info

        if not metric_files:
            info["incomplete"] = True
            findings.append(_report(
                "incomplete_zip", source_zip=zp.name,
                detail="ZIP sin metrics.tsv; no se ingiere" if
                require_complete else
                "ZIP sin metrics.tsv; se ingieren PDB por compatibilidad",
            ))

        for member in pdb_members:
            acc = Path(member).stem.strip()
            if not acc:
                continue
            base = _member_base(member)
            row = metrics_by_base.get(base, {}).get(acc)
            missing_metrics = not row
            if missing_metrics:
                detail = ("ZIP sin metrics.tsv" if not metric_files else
                          f"sin fila en metrics.tsv para {acc}")
                findings.append(_report(
                    "missing_metrics", accession=acc, source_zip=zp.name,
                    path=member, detail=detail,
                ))
                if require_complete:
                    info["incomplete"] = True
                    info["skipped"] += 1
                    continue
            status = ((row or {}).get("status") or "").strip().lower()
            if status not in OK_METRIC_STATUSES:
                info["skipped"] += 1
                findings.append(_report(
                    "metrics_status", accession=acc, source_zip=zp.name,
                    path=member, detail=f"status={status or '?'}; no se ingiere",
                ))
                continue
            try:
                data = zf.read(member)
            except (KeyError, zipfile.BadZipFile, OSError) as exc:
                info["incomplete"] = True
                info["skipped"] += 1
                findings.append(_report(
                    "corrupt_pdb", accession=acc, source_zip=zp.name,
                    path=member,
                    detail=f"no se pudo leer el miembro del ZIP: {exc}",
                ))
                continue
            problem = _pdb_problem(data)
            if problem:
                info["incomplete"] = True
                info["skipped"] += 1
                findings.append(_report(
                    "corrupt_pdb", accession=acc, source_zip=zp.name,
                    path=member,
                    sha256=hashlib.sha256(data).hexdigest(),
                    nbytes=len(data), detail=problem,
                ))
                continue
            candidates.append({
                "accession": acc,
                "member": member,
                "data": data,
                "sha256": hashlib.sha256(data).hexdigest(),
                "row": row or {},
            })
        if require_complete and info["incomplete"] and not any(
                f["kind"] == "incomplete_zip" for f in findings):
            findings.append(_report(
                "incomplete_zip", source_zip=zp.name,
                detail="ZIP incompleto (PDB truncado o métricas faltantes); "
                       "no se ingiere",
            ))
    return candidates, findings, info


def _scan_orphans(outdir: Path, cat: dict[str, dict]) -> list[dict]:
    """PDB en disco sin fila de catálogo. Solo reporta; no borra."""
    outdir = Path(outdir)
    if not outdir.is_dir():
        return []
    findings: list[dict] = []
    for path in sorted(outdir.glob("*.pdb")):
        acc = path.stem
        if acc in cat:
            continue
        try:
            data = path.read_bytes()
            sha, size = hashlib.sha256(data).hexdigest(), len(data)
        except OSError:
            sha, size = "", 0
        findings.append(_report(
            "orphan_pdb", accession=acc, source_zip="", path=str(path),
            sha256=sha, nbytes=size,
            detail="PDB presente en disco sin fila de catálogo; no se borra",
        ))
    return findings


def _prune_missing(outdir: Path, cat: dict[str, dict], dry_run: bool,
                   report_rows: list[dict],
                   skip: set[str] | None = None) -> int:
    """Reporta PDB locales ausentes y marca ``status=missing`` sin borrar.

    ``skip`` omite accesiones que esta corrida acaba de escribir (en
    ``dry_run`` todavía no existen, pero no son un faltante real).
    """
    marked = 0
    outdir = Path(outdir)
    skip = skip or set()
    for acc in sorted(cat):
        row = cat[acc]
        if acc in skip:
            continue
        target = outdir / f"{acc}.pdb"
        if target.exists():
            continue
        report_rows.append(_report(
            "missing_file", accession=acc,
            source_zip=row.get("source_zip", ""), path=str(target),
            sha256=row.get("sha256", ""), nbytes=row.get("bytes", ""),
            detail="PDB local ausente; fila conservada, no se borra nada",
        ))
        if not dry_run and row.get("status") == "ok":
            row["status"] = "missing"
            row["error"] = "pdb local ausente; fila conservada"
            marked += 1
    return marked


def _print_report(report_rows: list[dict], limit: int = 20) -> None:
    print("reporte de reingesta (solo lectura; no se borra nada):")
    for kind in REPORT_KINDS:
        rows = [r for r in report_rows if r["kind"] == kind]
        if not rows:
            continue
        print(f"  {kind}: {len(rows)}")
        for row in rows[:limit]:
            bits = [b for b in (row["accession"], row["source_zip"],
                                row["related"], row["detail"]) if b]
            where = f" [{row['path']}]" if row["path"] else ""
            print(f"    - {' | '.join(bits)}{where}")
        if len(rows) > limit:
            print(f"    ... ({len(rows) - limit} más)")


def ingest(zips: list[Path], outdir: Path, catalog_path: Path,
           dry_run: bool = False, *, require_complete: bool = False,
           report: bool = False, report_tsv: Path | None = None,
           prune_missing: bool = False, report_limit: int = 20) -> dict:
    """Ingiere ZIPs y devuelve conteos. ``dry_run`` no escribe catálogo/PDB.

    ``require_complete`` descarta por completo un ZIP sin metrics.tsv, con
    PDB truncado o con PDB sin fila de métricas. La CLI lo activa por
    defecto; el valor ``False`` conserva la compatibilidad histórica de los
    llamados que sí ingerían lotes sin metrics.

    ``report``/``prune_missing`` activan la auditoría. ``report_tsv``, si se
    indica explícitamente, se escribe incluso en ``dry_run`` porque es un
    diagnóstico solicitado; el catálogo y los PDB nunca se tocan en seco.
    """
    outdir = Path(outdir)
    catalog_path = Path(catalog_path)
    catalog_existed = catalog_path.exists()
    cat = _catalog(catalog_path)
    known_zips = {r.get("source_zip", "") for r in cat.values()
                  if r.get("source_zip")}
    report_rows: list[dict] = []
    audit = report or prune_missing or report_tsv is not None

    # Deduplicar rutas de entrada preservando el orden.
    unique_zips: list[Path] = []
    seen_paths: set[Path] = set()
    for zp in zips:
        zp = Path(zp)
        try:
            key = zp.resolve()
        except OSError:
            key = zp
        if key in seen_paths:
            continue
        seen_paths.add(key)
        unique_zips.append(zp)

    added = dup = improved = errors = written = skipped = 0
    incomplete_zips = 0
    zips_seen = 0
    batches: list[dict] = []
    touched: set[str] = set()
    seen_hashes: dict[str, str] = {
        row.get("sha256", ""): acc for acc, row in cat.items()
        if row.get("sha256") and row.get("status") == "ok"
    }
    hash_dups: set[tuple[str, str, str]] = set()
    # Origen(es) de cada accesión: catálogo previo + PDB aceptados de esta
    # corrida. Sirve para detectar duplicados entre lotes, no relecturas.
    occurrences: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for acc, row in cat.items():
        if row.get("source_zip"):
            occurrences[acc].add(
                (row.get("source_zip", ""), row.get("source_member", "")))

    for zp in unique_zips:
        zips_seen += 1
        zip_added = zip_dup = zip_improved = zip_skipped = 0
        if zp.name not in known_zips and audit:
            report_rows.append(_report(
                "new_zip", source_zip=zp.name, path=str(zp),
                detail="ZIP sin filas previas en el catálogo",
            ))
        candidates, findings, info = _inspect_zip(zp, require_complete)
        report_rows.extend(findings)
        if info["error"]:
            errors += 1
        strict_skip = bool(require_complete and info["incomplete"])
        if strict_skip:
            incomplete_zips += 1
            zip_skipped = info["pdbs"]
            skipped += zip_skipped
            batches.append({
                "zip": zp.name, "added": 0, "duplicate": 0,
                "improved": 0, "skipped": zip_skipped,
                "status": "incompleto",
            })
            print(f"{zp.name}: OMITIDO (ZIP incompleto); "
                  f"omitidos={zip_skipped}")
            continue
        zip_skipped = info["skipped"]
        skipped += zip_skipped

        for cand in candidates:
            acc = cand["accession"]
            member = cand["member"]
            data = cand["data"]
            sha = cand["sha256"]
            row = cand["row"]
            plddt = _num(row.get("plddt"))
            plddt = 0.0 if plddt is None else plddt
            occurrences[acc].add((zp.name, member))
            old = cat.get(acc)
            old_ok = bool(old and old.get("status") == "ok")
            is_new = not old_ok
            if old_ok:
                old_plddt = _num(old.get("plddt")) or 0.0
                if old_plddt >= plddt:
                    dup += 1
                    zip_dup += 1
                    # Trazabilidad del lote desplazado, sin tocar la copia
                    # conservada ni repetir el mismo ZIP.
                    if (zp.name != old.get("source_zip")
                            and not old.get("duplicate_of")):
                        old["duplicate_of"] = zp.name
                    continue
                improved += 1
                zip_improved += 1
                old_sha = old.get("sha256", "")
                if old_sha and seen_hashes.get(old_sha) == acc:
                    seen_hashes.pop(old_sha, None)
            target = outdir / f"{acc}.pdb"
            if not dry_run:
                outdir.mkdir(parents=True, exist_ok=True)
                tmp = target.with_suffix(".pdb.tmp")
                tmp.write_bytes(data)
                os.replace(tmp, target)
                written += 1
            cat[acc] = {
                "accession": acc,
                "status": "ok",
                "plddt": f"{plddt:.3f}",
                "ptm": row.get("ptm", ""),
                "length": row.get("length", ""),
                "source_zip": zp.name,
                "source_member": member,
                "sha256": sha,
                "bytes": len(data),
                "duplicate_of": old.get("source_zip", "") if old else "",
                "error": "",
                "ingested_at": "",  # se rellena al final
            }
            touched.add(acc)
            prev = seen_hashes.get(sha)
            if prev and prev != acc:
                key = (sha, prev, acc)
                if key not in hash_dups:
                    report_rows.append(_report(
                        "duplicate_hash", accession=acc,
                        source_zip=zp.name, path=member, sha256=sha,
                        nbytes=len(data), related=prev,
                        detail=f"sha256 idéntico a {prev}",
                    ))
                hash_dups.add(key)
            seen_hashes.setdefault(sha, acc)
            if is_new:
                added += 1
                zip_added += 1

        batches.append({
            "zip": zp.name, "added": zip_added, "duplicate": zip_dup,
            "improved": zip_improved, "skipped": zip_skipped,
            "status": "incompleto" if info["incomplete"] else "ok",
        })
        print(f"{zp.name}: nuevos={zip_added} duplicados={zip_dup} "
              f"mejorados={zip_improved} omitidos={zip_skipped}")

    # Duplicados entre lotes: misma accesión en ZIP distintos (o repetida
    # dentro de un ZIP). Una relectura del mismo lote no cuenta como duplicado.
    for acc in sorted(occurrences):
        sources = occurrences[acc]
        if len(sources) <= 1:
            continue
        kept = cat.get(acc, {})
        kept_key = (kept.get("source_zip", ""), kept.get("source_member", ""))
        for src_zip, src_member in sorted(sources):
            if (src_zip, src_member) == kept_key:
                continue
            report_rows.append(_report(
                "duplicate_accession", accession=acc, source_zip=src_zip,
                path=src_member, related=kept.get("source_zip", ""),
                detail="accesión presente en más de un lote; se conserva "
                       "la copia de mayor pLDDT",
            ))

    pruned_missing = 0
    if prune_missing:
        pruned_missing = _prune_missing(outdir, cat, dry_run, report_rows,
                                        skip=touched)
    orphan_pdbs = 0
    if audit:
        orphan_rows = _scan_orphans(outdir, cat)
        report_rows.extend(orphan_rows)
        orphan_pdbs = len(orphan_rows)

    if not dry_run and (cat or catalog_existed):
        now = time.strftime("%Y-%m-%dT%H:%M:%S")
        for r in cat.values():
            if not r.get("ingested_at"):
                r["ingested_at"] = now
        _write(catalog_path, cat)

    report_rows.sort(key=lambda r: (
        REPORT_KINDS.index(r["kind"])
        if r["kind"] in REPORT_KINDS else len(REPORT_KINDS),
        r["accession"], r["source_zip"], r["path"], r["related"],
        r["detail"],
    ))
    counts = Counter(r["kind"] for r in report_rows)
    if report_tsv is not None:
        write_tsv_rows(Path(report_tsv), REPORT_COLUMNS, report_rows)
    if report or prune_missing:
        _print_report(report_rows, report_limit)

    print(f"total catálogo={len(cat)} nuevos={added} duplicados={dup} "
          f"mejorados={improved} omitidos={skipped} errores={errors}")
    if counts.get("new_zip") or counts.get("orphan_pdb") or \
            counts.get("incomplete_zip"):
        print(f"reporte: zips_nuevos={counts.get('new_zip', 0)} "
              f"pdb_sin_catalogo={counts.get('orphan_pdb', 0)} "
              f"zip_incompletos={counts.get('incomplete_zip', 0)} "
              f"accesiones_duplicadas={counts.get('duplicate_accession', 0)} "
              f"pdb_corruptos={counts.get('corrupt_pdb', 0)}")
    if prune_missing:
        print(f"prune-missing: pdb_ausentes={counts.get('missing_file', 0)} "
              f"marcados={pruned_missing} (no se borró nada)")
    if hash_dups:
        print(f"AVISO hashes duplicados entre accesiones: {len(hash_dups)}")
    return {
        "total": len(cat),
        "added": added,
        "duplicate": dup,
        "improved": improved,
        "error": errors,
        "errors": errors,
        "written": written,
        "skipped": skipped,
        "zips": zips_seen,
        "incomplete_zips": incomplete_zips,
        "new_zips": counts.get("new_zip", 0),
        "orphan_pdbs": orphan_pdbs,
        "duplicate_accessions": counts.get("duplicate_accession", 0),
        "duplicate_hashes": len(hash_dups),
        "hash_duplicates": len(hash_dups),
        "hash_duplicate_pairs": sorted(f"{a}<->{b}" for _, a, b in hash_dups),
        "missing_files": counts.get("missing_file", 0),
        "pruned_missing": pruned_missing,
        "report_rows": report_rows,
        "report_counts": dict(counts),
        "batches": batches,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zips", type=Path, default=ROOT / "data" / "Kaggler",
                    help="carpeta con ZIP o un ZIP concreto")
    ap.add_argument("--outdir", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets"
                    / "esmfold_kaggle")
    ap.add_argument("--catalog", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets"
                    / "esmfold_catalog.tsv")
    ap.add_argument("--dry-run", action="store_true",
                    help="calcula y reporta, sin escribir catálogo ni PDB")
    ap.add_argument("--report", action="store_true",
                    help="audita ZIP nuevos, PDB sin catálogo y duplicados")
    ap.add_argument("--prune-missing", action="store_true",
                    help="reporta y marca PDB locales ausentes; no borra nada")
    ap.add_argument("--report-tsv", type=Path, default=None,
                    help="escribe el reporte de auditoría en este TSV")
    ap.add_argument("--report-limit", type=int, default=20,
                    help="filas impresas por tipo de hallazgo (default 20)")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="compatibilidad: ingestar ZIP sin metrics.tsv")
    a = ap.parse_args(argv)
    zips = (sorted(a.zips.glob("*.zip")) if a.zips.is_dir()
            else [a.zips])
    res = ingest(
        zips, a.outdir, a.catalog, a.dry_run,
        require_complete=not a.allow_incomplete,
        report=a.report or a.prune_missing,
        report_tsv=a.report_tsv,
        prune_missing=a.prune_missing,
        report_limit=max(0, a.report_limit),
    )
    return 0 if res["error"] == 0 and res["incomplete_zips"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
