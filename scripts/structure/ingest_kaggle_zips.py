#!/usr/bin/env python3
"""Ingiere ZIP de ESMFold de data/Kaggler sin duplicar accesiones.

Cada ZIP contiene ``<lote>/metrics.tsv`` y ``<lote>/pdb/<acc>.pdb``.
El catálogo conserva pLDDT/pTM/longitud, origen y hash. La copia local queda
en ``data/metabolism/targets/esmfold_kaggle/<acc>.pdb``. Si una accesión
aparece en dos lotes, se conserva el PDB con mayor pLDDT y se registra el
duplicado; nunca se sobrescribe con una copia de peor calidad.

El lector/escritor TSV se exporta para que
``unify_structure_catalog.py`` reutilice exactamente el mismo formato.

Uso:
  python scripts/structure/ingest_kaggle_zips.py --dry-run
  python scripts/structure/ingest_kaggle_zips.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COLS = ["accession", "status", "plddt", "ptm", "length", "source_zip",
        "source_member", "sha256", "bytes", "duplicate_of", "error",
        "ingested_at"]
CATALOG_COLUMNS = COLS  # alias explícito para otros módulos


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


def _catalog(path: Path) -> dict[str, dict]:
    return {r["accession"]: r for r in read_tsv_rows(path)
            if r.get("accession")}


def _write(path: Path, rows: dict[str, dict]) -> None:
    write_tsv_rows(path, COLS, [rows[acc] for acc in sorted(rows)])


def _num(value) -> float | None:
    try:
        return float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return None


def _metrics(zf: zipfile.ZipFile, member: str) -> dict[str, dict]:
    try:
        text = zf.read(member).decode("utf-8-sig", errors="replace")
    except KeyError:
        return {}
    out = {}
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        acc = (r.get("uniprot") or r.get("accession") or r.get("acc")
               or r.get("protein_id") or "").strip()
        if acc:
            out[acc] = r
    return out


def ingest(zips: list[Path], outdir: Path, catalog_path: Path,
           dry_run: bool = False) -> dict:
    """Ingiere ZIPs y devuelve conteos. ``dry_run`` no escribe nada."""
    outdir = Path(outdir)
    catalog_path = Path(catalog_path)
    cat = _catalog(catalog_path)
    added = dup = improved = errors = written = 0
    seen_hashes: dict[str, str] = {
        row.get("sha256", ""): acc for acc, row in cat.items()
        if row.get("sha256") and row.get("status") == "ok"
    }
    hash_dups: set[tuple[str, str, str]] = set()

    for zp in zips:
        zp = Path(zp)
        zip_added = zip_dup = zip_improved = 0
        try:
            zf = zipfile.ZipFile(zp)
        except zipfile.BadZipFile as exc:
            print(f"ERROR {zp.name}: {exc}")
            errors += 1
            continue
        with zf:
            names = [n for n in zf.namelist()
                     if not n.startswith("__MACOSX/") and not n.endswith("/")]
            metrics_by_base: dict[str, dict[str, dict]] = {}
            metric_files = [n for n in names if n.endswith("/metrics.tsv")
                            or n == "metrics.tsv"]
            for mpath in metric_files:
                base = mpath.rsplit("/", 1)[0] if "/" in mpath else ""
                metrics_by_base[base] = _metrics(zf, mpath)
            if not metric_files:
                print(f"AVISO {zp.name}: sin metrics.tsv; se ingieren PDB igual")
            pdb_members = [n for n in names if n.lower().endswith(".pdb")]
            if not pdb_members:
                print(f"AVISO {zp.name}: sin PDB")
                continue
            for member in pdb_members:
                if "/pdb/" in member:
                    base = member.rsplit("/pdb/", 1)[0]
                else:
                    base = member.rsplit("/", 1)[0] if "/" in member else ""
                row = metrics_by_base.get(base, {}).get(Path(member).stem, {})
                acc = Path(member).stem.strip()
                if not acc:
                    continue
                if (row.get("status") or "").strip().lower() not in (
                        "", "done", "ok", "already_present"):
                    continue
                try:
                    data = zf.read(member)
                except KeyError:
                    continue
                if not data:
                    continue
                sha = hashlib.sha256(data).hexdigest()
                plddt = _num(row.get("plddt"))
                plddt = 0.0 if plddt is None else plddt
                old = cat.get(acc)
                old_ok = bool(old and old.get("status") == "ok")
                if old_ok:
                    old_plddt = _num(old.get("plddt")) or 0.0
                    if old_plddt >= plddt:
                        dup += 1
                        zip_dup += 1
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
                prev = seen_hashes.get(sha)
                if prev and prev != acc:
                    hash_dups.add((sha, prev, acc))
                seen_hashes.setdefault(sha, acc)
                added += 1
                zip_added += 1
        print(f"{zp.name}: nuevos={zip_added} duplicados={zip_dup} "
              f"mejorados={zip_improved}")

    if not dry_run:
        now = time.strftime("%Y-%m-%dT%H:%M:%S")
        for r in cat.values():
            if not r.get("ingested_at"):
                r["ingested_at"] = now
        _write(catalog_path, cat)
    print(f"total catálogo={len(cat)} nuevos={added} duplicados={dup} "
          f"mejorados={improved} errores={errors}")
    if hash_dups:
        print(f"AVISO hashes duplicados entre accesiones: {len(hash_dups)}")
    return {
        "total": len(cat),
        "added": added,
        "duplicate": dup,
        "improved": improved,
        "error": errors,
        "errors": errors,
        "hash_duplicates": len(hash_dups),
        "hash_duplicate_pairs": sorted(f"{a}<->{b}" for _, a, b in hash_dups),
        "written": written,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zips", type=Path, default=ROOT / "data" / "Kaggler")
    ap.add_argument("--outdir", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets" / "esmfold_kaggle")
    ap.add_argument("--catalog", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets" / "esmfold_catalog.tsv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    zips = sorted(a.zips.glob("*.zip")) if a.zips.is_dir() else [a.zips]
    res = ingest(zips, a.outdir, a.catalog, a.dry_run)
    return 0 if res["error"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
