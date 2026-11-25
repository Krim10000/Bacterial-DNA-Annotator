#!/usr/bin/env python3
"""Ingiere ZIP de ESMFold de data/Kaggler sin duplicar accesiones.

Cada ZIP contiene ``<lote>/metrics.tsv`` y ``<lote>/pdb/<acc>.pdb``.
El catálogo conserva pLDDT/pTM/longitud, origen y hash. La copia local queda
en ``data/metabolism/targets/esmfold_kaggle/<acc>.pdb``. Si una accesión
aparece en dos lotes, se conserva el PDB con mayor pLDDT y se registra el
duplicado; nunca se sobrescribe con una copia de peor calidad.

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
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COLS = ["accession", "status", "plddt", "ptm", "length", "source_zip",
        "source_member", "sha256", "bytes", "duplicate_of", "error",
        "ingested_at"]


def _catalog(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("accession")}


def _write(path: Path, rows: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for acc in sorted(rows):
            w.writerow(rows[acc])
    os.replace(tmp, path)


def _metrics(zf: zipfile.ZipFile, member: str) -> dict[str, dict]:
    try:
        text = zf.read(member).decode("utf-8", errors="replace")
    except KeyError:
        return {}
    out = {}
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        acc = (r.get("uniprot") or r.get("accession") or "").strip()
        if acc:
            out[acc] = r
    return out


def ingest(zips: list[Path], outdir: Path, catalog_path: Path,
           dry_run: bool = False) -> dict:
    cat = _catalog(catalog_path)
    added = dup = skipped = errors = 0
    for zp in zips:
        try:
            zf = zipfile.ZipFile(zp)
        except zipfile.BadZipFile as exc:
            print(f"ERROR {zp.name}: {exc}")
            errors += 1
            continue
        with zf:
            names = zf.namelist()
            mets = [n for n in names if n.endswith("/metrics.tsv")]
            if not mets:
                print(f"AVISO {zp.name}: sin metrics.tsv")
                continue
            for mpath in mets:
                base = mpath.rsplit("/", 1)[0]
                metrics = _metrics(zf, mpath)
                for member in names:
                    if not member.startswith(base + "/pdb/") or not member.endswith(".pdb"):
                        continue
                    acc = Path(member).stem
                    row = metrics.get(acc, {})
                    if row.get("status") not in (None, "", "done"):
                        continue
                    try:
                        data = zf.read(member)
                    except KeyError:
                        continue
                    sha = hashlib.sha256(data).hexdigest()
                    try:
                        plddt = float(row.get("plddt") or 0)
                    except ValueError:
                        plddt = 0.0
                    old = cat.get(acc)
                    target = outdir / f"{acc}.pdb"
                    if old and old.get("status") == "ok":
                        if float(old.get("plddt") or 0) >= plddt:
                            dup += 1
                            continue
                        skipped += 1
                    if not dry_run:
                        outdir.mkdir(parents=True, exist_ok=True)
                        tmp = target.with_suffix(".pdb.tmp")
                        tmp.write_bytes(data)
                        os.replace(tmp, target)
                    cat[acc] = {
                        "accession": acc, "status": "ok", "plddt": f"{plddt:.3f}",
                        "ptm": row.get("ptm", ""), "length": row.get("length", ""),
                        "source_zip": zp.name, "source_member": member,
                        "sha256": sha, "bytes": len(data),
                        "duplicate_of": old.get("source_zip", "") if old else "",
                        "error": "", "ingested_at": "",  # se rellena al final
                    }
                    added += 1
        print(f"{zp.name}: PDB válidos={added} duplicados={dup}")
    if not dry_run:
        import time
        for r in cat.values():
            if not r.get("ingested_at"):
                r["ingested_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _write(catalog_path, cat)
    print(f"total catálogo={len(cat)} nuevos={added} duplicados={dup} "
          f"mejorados={skipped} errores={errors}")
    return {"total": len(cat), "added": added, "duplicate": dup,
            "improved": skipped, "error": errors}


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
