"""Contratos del catálogo unificado AFDB + ESMFold, sin red.

Usa catálogos y PDB sintéticos en ``tmp_path``; nunca toca
``data/metabolism/targets``.
"""
from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from pathlib import Path

from scripts.structure.ingest_kaggle_zips import ingest
from scripts.structure.unify_structure_catalog import (
    CATALOG_COLUMNS,
    REPORT_COLUMNS,
    main,
    unify,
)

AFDB_COLS = [
    "accession", "status", "pdb", "pae", "sha256", "pdb_bytes",
    "pae_bytes", "pdb_url", "pae_url", "error", "updated",
]
ESM_COLS = [
    "accession", "status", "plddt", "ptm", "length", "source_zip",
    "source_member", "sha256", "bytes", "duplicate_of", "error",
    "ingested_at",
]


# ---------------------------------------------------------------------------
# Utilidades sintéticas
# ---------------------------------------------------------------------------


def _atom(serial: int, resseq: int, bfactor: float, chain: str = "A") -> str:
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain:1s}{resseq:>4d}    "
        f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
        f"{1.0:>6.2f}{bfactor:>6.2f}           C\n"
    )


def _pdb_text(bfactor: float = 80.0, residues: int = 3,
              chain: str = "A") -> str:
    lines = [_atom(i, i, bfactor, chain) for i in range(1, residues + 1)]
    return "".join(lines) + "END\n"


def _write_pdb(path: Path, bfactor: float = 80.0, residues: int = 3,
               chain: str = "A") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_pdb_text(bfactor, residues, chain), encoding="utf-8")
    return path


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _afdb_row(acc: str, pdb: Path | str = "", *, status: str = "ok",
              updated: str = "2026-01-01T00:00:00") -> dict:
    pdb_text = str(pdb) if pdb else ""
    pdb_path = Path(pdb) if pdb else None
    exists = bool(pdb_path and pdb_path.exists())
    sha = _sha(pdb_path) if exists else ""
    size = pdb_path.stat().st_size if exists else 0
    return {
        "accession": acc, "status": status, "pdb": pdb_text, "pae": "",
        "sha256": sha, "pdb_bytes": str(size), "pae_bytes": "0",
        "pdb_url": "", "pae_url": "", "error": "", "updated": updated,
    }


def _esm_row(acc: str, pdb: Path, *, plddt: float = 0.800,
             ptm: float = 0.7000, length: int = 3, status: str = "ok",
             updated: str = "2026-01-02T00:00:00") -> dict:
    return {
        "accession": acc, "status": status, "plddt": f"{plddt:.3f}",
        "ptm": f"{ptm:.4f}", "length": str(length), "source_zip": "lote.zip",
        "source_member": f"lote/pdb/{acc}.pdb", "sha256": _sha(pdb),
        "bytes": str(pdb.stat().st_size), "duplicate_of": "", "error": "",
        "ingested_at": updated,
    }


def _paths(tmp_path: Path) -> dict[str, Path]:
    return {
        "afdb_dir": tmp_path / "afdb25k",
        "esm_dir": tmp_path / "esmfold_kaggle",
        "afdb_cat": tmp_path / "afdb_catalog.tsv",
        "esm_cat": tmp_path / "esmfold_catalog.tsv",
        "out": tmp_path / "out" / "structure_catalog.tsv",
        "report_tsv": tmp_path / "out" / "report.tsv",
        "report_md": tmp_path / "out" / "report.md",
    }


def _run(tmp_path: Path, afdb_rows: list[dict], esm_rows: list[dict],
         **kwargs):
    p = _paths(tmp_path)
    p["afdb_dir"].mkdir(exist_ok=True)
    p["esm_dir"].mkdir(exist_ok=True)
    _write_tsv(p["afdb_cat"], AFDB_COLS, afdb_rows)
    _write_tsv(p["esm_cat"], ESM_COLS, esm_rows)
    kwargs.setdefault("afdb_dir", p["afdb_dir"])
    kwargs.setdefault("esmfold_dir", p["esm_dir"])
    kwargs.setdefault("out", p["out"])
    kwargs.setdefault("report_tsv", p["report_tsv"])
    kwargs.setdefault("report_md", p["report_md"])
    result = unify(p["afdb_cat"], p["esm_cat"], **kwargs)
    result["_paths"] = p
    return result


def _read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------


def test_schema_and_prefer_afdb(tmp_path):
    p = _paths(tmp_path)
    afdb_file = _write_pdb(p["afdb_dir"] / "P1.pdb", bfactor=80.0, residues=3)
    esm_file = _write_pdb(p["esm_dir"] / "P1.pdb", bfactor=90.0, residues=4)
    esm_only = _write_pdb(p["esm_dir"] / "P4.pdb", bfactor=85.0, residues=5)
    result = _run(
        tmp_path,
        [_afdb_row("P1", afdb_file)],
        [_esm_row("P1", esm_file, plddt=0.950),
         _esm_row("P4", esm_only, plddt=0.900, ptm=0.8000, length=5)],
    )
    s = result["summary"]
    assert s["afdb_rows"] == 1
    assert s["esmfold_rows"] == 2
    assert s["unified"] == 2
    assert s["source_AFDB"] == 1
    assert s["source_ESMFold"] == 1
    assert s["duplicate_accessions"] == 1
    assert s["written"] is True

    rows = _read_tsv(result["_paths"]["out"])
    assert list(rows[0]) == CATALOG_COLUMNS
    by_acc = {r["accession"]: r for r in rows}
    assert by_acc["P1"]["source"] == "AFDB"
    assert by_acc["P1"]["plddt"] == "0.800"
    assert by_acc["P1"]["ptm"] == ""
    assert by_acc["P1"]["length"] == "3"
    assert by_acc["P1"]["status"] == "ok"
    assert by_acc["P4"]["source"] == "ESMFold"
    assert by_acc["P4"]["plddt"] == "0.900"
    assert by_acc["P4"]["ptm"] == "0.8000"

    report = _read_tsv(result["_paths"]["report_tsv"])
    assert list(report[0]) == REPORT_COLUMNS
    assert any(r["kind"] == "duplicate_accession" and r["source"] == "ESMFold"
               and r["accession"] == "P1" for r in report)
    # No se borran archivos: ambos P1 siguen en disco.
    assert afdb_file.exists() and esm_file.exists()
    md = result["_paths"]["report_md"].read_text(encoding="utf-8")
    assert "# Catálogo unificado de estructuras" in md
    assert "## Límites" in md


def test_esmfold_kept_when_afdb_unusable(tmp_path):
    p = _paths(tmp_path)
    esm_file = _write_pdb(p["esm_dir"] / "P2.pdb", bfactor=70.0, residues=2)
    missing = p["afdb_dir"] / "missing" / "P3.pdb"
    result = _run(
        tmp_path,
        [_afdb_row("P2", "", status="error"),
         _afdb_row("P3", missing)],
        [_esm_row("P2", esm_file, plddt=0.710, length=2)],
        dry_run=True,
    )
    rows = result["rows"]
    assert rows["P2"]["source"] == "ESMFold"
    assert rows["P2"]["status"] == "ok"
    assert rows["P2"]["plddt"] == "0.710"
    assert rows["P3"]["source"] == "AFDB"
    assert rows["P3"]["status"] == "missing"
    assert result["summary"]["missing_files"] == 1
    assert result["summary"]["written"] is False
    report_kinds = {r["kind"] for r in result["report_rows"]}
    assert "duplicate_accession" in report_kinds
    assert "missing_file" in report_kinds


def test_hash_duplicates_marked_not_deleted(tmp_path):
    p = _paths(tmp_path)
    a1 = _write_pdb(p["afdb_dir"] / "A1.pdb", bfactor=70.0, residues=3,
                    chain="A")
    a2 = _write_pdb(p["esm_dir"] / "A2.pdb", bfactor=70.0, residues=3,
                    chain="A")
    a3 = _write_pdb(p["esm_dir"] / "A3.pdb", bfactor=70.0, residues=3,
                    chain="A")
    b1 = _write_pdb(p["esm_dir"] / "B1.pdb", bfactor=71.0, residues=3,
                    chain="B")
    b2 = _write_pdb(p["esm_dir"] / "B2.pdb", bfactor=71.0, residues=3,
                    chain="B")
    result = _run(
        tmp_path,
        [_afdb_row("A1", a1)],
        [_esm_row("A2", a2, plddt=0.900),
         _esm_row("A3", a3, plddt=0.990),
         _esm_row("B1", b1, plddt=0.600),
         _esm_row("B2", b2, plddt=0.950)],
    )
    rows = result["rows"]
    assert rows["A1"]["duplicate_of"] == ""
    assert rows["A2"]["duplicate_of"] == "A1"
    assert rows["A3"]["duplicate_of"] == "A1"
    assert rows["B2"]["duplicate_of"] == ""
    assert rows["B1"]["duplicate_of"] == "B2"
    assert result["summary"]["hash_duplicates"] == 3
    for path in (a1, a2, a3, b1, b2):
        assert path.exists()
    kinds = [r["kind"] for r in result["report_rows"]]
    assert kinds.count("duplicate_hash") == 3


def test_orphans_path_mismatch_and_duplicate_rows(tmp_path):
    p = _paths(tmp_path)
    registered = _write_pdb(p["afdb_dir"] / "P1.pdb", bfactor=80.0)
    orphan = _write_pdb(p["afdb_dir"] / "ORPHAN.pdb", bfactor=60.0)
    moved = _write_pdb(p["afdb_dir"] / "P3.pdb", bfactor=75.0)
    other = p["afdb_dir"] / "alt" / "P3.pdb"  # ruta declarada, no existe
    no_route = _write_pdb(p["afdb_dir"] / "P4.pdb", bfactor=65.0)
    result = _run(
        tmp_path,
        [
            _afdb_row("P1", registered),
            # P3 se registra en alt/, pero el archivo huérfano está en raíz.
            _afdb_row("P3", other),
            # P4 está en disco pero su fila no declara ruta.
            _afdb_row("P4", "", status="error"),
        ],
        [],
        dry_run=True,
    )
    s = result["summary"]
    assert s["orphans"] == 2  # ORPHAN y P4 (registrada sin ruta)
    assert s["path_mismatches"] == 1  # P3 en raíz, el catálogo apunta a alt/
    assert s["missing_files"] == 1  # alt/P3.pdb
    kinds = {r["kind"] for r in result["report_rows"]}
    assert {"orphan_pdb", "path_mismatch", "missing_file"} <= kinds
    assert orphan.exists() and moved.exists() and no_route.exists()


def test_dry_run_writes_nothing_and_rerun_is_idempotent(tmp_path):
    p = _paths(tmp_path)
    afdb_file = _write_pdb(p["afdb_dir"] / "P1.pdb", bfactor=80.0)
    esm_file = _write_pdb(p["esm_dir"] / "P2.pdb", bfactor=82.0)
    afdb_rows = [_afdb_row("P1", afdb_file)]
    esm_rows = [_esm_row("P2", esm_file, plddt=0.820)]

    dry = _run(tmp_path, afdb_rows, esm_rows, dry_run=True)
    assert dry["summary"]["unified"] == 2
    assert dry["summary"]["written"] is False
    assert not p["out"].exists()
    assert not p["report_tsv"].exists()
    assert not p["report_md"].exists()

    _run(tmp_path, afdb_rows, esm_rows)
    first = p["out"].read_bytes()
    first_tsv = p["report_tsv"].read_bytes()
    first_md = p["report_md"].read_bytes()
    _run(tmp_path, afdb_rows, esm_rows)
    assert p["out"].read_bytes() == first
    assert p["report_tsv"].read_bytes() == first_tsv
    assert p["report_md"].read_bytes() == first_md
    assert p["out"].read_text(encoding="utf-8").count("\n") == 3


def test_duplicate_accession_within_source_keeps_best(tmp_path):
    p = _paths(tmp_path)
    low = _write_pdb(p["afdb_dir"] / "P1_low.pdb", bfactor=60.0)
    high = _write_pdb(p["afdb_dir"] / "P1_high.pdb", bfactor=90.0)
    result = _run(
        tmp_path,
        [_afdb_row("P1", low), _afdb_row("P1", high)],
        [],
        dry_run=True,
    )
    assert result["summary"]["duplicate_accessions"] == 1
    assert result["rows"]["P1"]["plddt"] == "0.900"
    assert result["rows"]["P1"]["path"].endswith("P1_high.pdb")
    assert any(r["kind"] == "duplicate_accession"
               and r["accession"] == "P1"
               and r["path"].endswith("P1_low.pdb")
               for r in result["report_rows"])


def test_ingest_synthetic_zip_then_unify(tmp_path):
    p = _paths(tmp_path)
    p["afdb_dir"].mkdir(parents=True, exist_ok=True)
    p["esm_dir"].mkdir(parents=True, exist_ok=True)
    # Sin catálogo AFDB previo: solo se escribe el encabezado.
    _write_tsv(p["afdb_cat"], AFDB_COLS, [])
    zip_path = tmp_path / "lote.zip"
    pdb = _pdb_text(bfactor=88.0, residues=3)
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(
            "lote/metrics.tsv",
            "uniprot\tstatus\tlength\tplddt\tptm\n"
            "Z1\tdone\t3\t0.880\t0.7000\n",
        )
        zf.writestr("lote/pdb/Z1.pdb", pdb)
    # Segundo ZIP sin metrics.tsv: el PDB debe ingerirse igualmente.
    z2 = tmp_path / "sin_metrics.zip"
    pdb2 = _pdb_text(bfactor=66.0, residues=2)
    with zipfile.ZipFile(z2, "w") as zf:
        zf.writestr("otro/pdb/Z2.pdb", pdb2)
    res = ingest([zip_path, z2], p["esm_dir"], p["esm_cat"])
    assert res["added"] == 2
    assert (p["esm_dir"] / "Z1.pdb").exists()
    assert (p["esm_dir"] / "Z2.pdb").exists()

    result = unify(
        p["afdb_cat"], p["esm_cat"], afdb_dir=p["afdb_dir"],
        esmfold_dir=p["esm_dir"], out=p["out"],
        report_tsv=p["report_tsv"], report_md=p["report_md"],
    )
    rows = result["rows"]
    assert rows["Z1"]["source"] == "ESMFold"
    assert rows["Z1"]["plddt"] == "0.880"
    assert rows["Z1"]["ptm"] == "0.7000"
    assert rows["Z2"]["source"] == "ESMFold"
    assert rows["Z2"]["status"] == "ok"
    assert result["summary"]["source_ESMFold"] == 2


def test_relative_afdb_paths_resolve_from_catalog(tmp_path):
    p = _paths(tmp_path)
    p["afdb_dir"].mkdir(parents=True, exist_ok=True)
    p["esm_dir"].mkdir(parents=True, exist_ok=True)
    _write_pdb(p["afdb_dir"] / "P1.pdb", bfactor=80.0)
    row = _afdb_row("P1", p["afdb_dir"] / "P1.pdb")
    row["pdb"] = "afdb25k/P1.pdb"  # relativa al catálogo
    missing = _afdb_row("P9", "")
    missing["pdb"] = "afdb25k/P9.pdb"  # relativa, no existe
    _write_tsv(p["afdb_cat"], AFDB_COLS, [row, missing])
    _write_tsv(p["esm_cat"], ESM_COLS, [])
    result = unify(
        p["afdb_cat"], p["esm_cat"], afdb_dir=p["afdb_dir"],
        esmfold_dir=p["esm_dir"], out=p["out"],
        report_tsv=p["report_tsv"], report_md=p["report_md"],
        dry_run=True,
    )
    assert result["rows"]["P1"]["status"] == "ok"
    assert result["rows"]["P1"]["path"].endswith("afdb25k/P1.pdb")
    assert result["rows"]["P9"]["status"] == "missing"
    assert result["rows"]["P9"]["path"].endswith("afdb25k/P9.pdb")


def test_cli_dry_run_json(tmp_path, capsys):
    p = _paths(tmp_path)
    afdb_file = _write_pdb(p["afdb_dir"] / "P1.pdb", bfactor=80.0)
    _run(tmp_path, [_afdb_row("P1", afdb_file)], [], dry_run=True)
    rc = main([
        "--afdb-catalog", str(p["afdb_cat"]),
        "--esmfold-catalog", str(p["esm_cat"]),
        "--afdb-dir", str(p["afdb_dir"]),
        "--esmfold-dir", str(p["esm_dir"]),
        "--out", str(p["out"]),
        "--report-tsv", str(p["report_tsv"]),
        "--report-md", str(p["report_md"]),
        "--dry-run", "--json",
    ])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["unified"] == 1
    assert payload["written"] is False
    assert not p["out"].exists()
