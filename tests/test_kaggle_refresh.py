"""Reingesta incremental de ZIP ESMFold de ``data/Kaggler``.

Cubre el refresco real: ZIP nuevos, accesiones duplicadas entre lotes y lotes
corruptos o incompletos. Todo ocurre en ``tmp_path``; no toca
``data/metabolism/targets``.
"""
from __future__ import annotations

import csv
import sys
import zipfile
from pathlib import Path

from scripts.structure.ingest_kaggle_zips import (
    COLS,
    REPORT_COLUMNS,
    ingest,
    main,
)


# ---------------------------------------------------------------------------
# Utilidades sintéticas
# ---------------------------------------------------------------------------


def _atom(serial: int, resseq: int, bfactor: float = 80.0,
          chain: str = "A") -> str:
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain:1s}{resseq:>4d}    "
        f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
        f"{1.0:>6.2f}{bfactor:>6.2f}           C\n"
    )


def _pdb_text(bfactor: float = 80.0, residues: int = 3,
              terminate: bool = True) -> str:
    body = "".join(_atom(i, i, bfactor) for i in range(1, residues + 1))
    return body + ("END\n" if terminate else "")


def _make_zip(path: Path, batch: str, pdbs: dict[str, str],
              metrics: dict[str, dict] | None = None) -> Path:
    """``metrics=None`` genera un ZIP sin ``metrics.tsv`` (incompleto)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as zf:
        if metrics is not None:
            lines = ["uniprot\tstatus\tlength\tplddt\tptm"]
            for acc, row in metrics.items():
                lines.append(
                    f"{acc}\t{row.get('status', 'done')}\t"
                    f"{row.get('length', 3)}\t{row.get('plddt', 0.8)}\t"
                    f"{row.get('ptm', 0.7)}"
                )
            zf.writestr(f"{batch}/metrics.tsv", "\n".join(lines) + "\n")
        for acc, text in pdbs.items():
            zf.writestr(f"{batch}/pdb/{acc}.pdb", text)
    return path


def _read_catalog(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# ZIP nuevo, campos del catálogo e idempotencia
# ---------------------------------------------------------------------------


def test_new_zip_catalogs_fields_and_rerun_is_idempotent(tmp_path):
    zp = _make_zip(
        tmp_path / "kaggle_nuevo.zip", "lote_nuevo",
        {"P1": _pdb_text(80.0), "P2": _pdb_text(70.0, residues=4)},
        {"P1": {"plddt": 0.910, "ptm": 0.8100, "length": 3},
         "P2": {"plddt": 0.720, "ptm": 0.6500, "length": 4}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"

    res = ingest([zp], out, cat, require_complete=True, report=True)
    assert res["added"] == 2
    assert res["duplicate"] == 0
    assert res["improved"] == 0
    assert res["skipped"] == 0
    assert res["incomplete_zips"] == 0
    assert res["new_zips"] == 1
    assert res["error"] == 0

    rows = _read_catalog(cat)
    assert list(rows[0]) == COLS
    by_acc = {r["accession"]: r for r in rows}
    assert by_acc["P1"]["source_zip"] == "kaggle_nuevo.zip"
    assert by_acc["P1"]["plddt"] == "0.910"
    assert by_acc["P1"]["ptm"] == "0.81"
    assert by_acc["P1"]["length"] == "3"
    assert len(by_acc["P1"]["sha256"]) == 64
    assert by_acc["P1"]["duplicate_of"] == ""
    assert (out / "P1.pdb").exists()
    assert (out / "P2.pdb").exists()

    # Reejecutar con las mismas entradas no cambia ni un byte ni duplica filas.
    before = cat.read_bytes()
    res2 = ingest([zp], out, cat, require_complete=True)
    assert res2["added"] == 0
    assert res2["duplicate"] == 2
    assert res2["improved"] == 0
    assert res2["written"] == 0
    assert cat.read_bytes() == before
    assert len(_read_catalog(cat)) == 2


def test_dry_run_counts_but_writes_nothing(tmp_path):
    zp = _make_zip(
        tmp_path / "nuevo.zip", "lote",
        {"P1": _pdb_text()}, {"P1": {"plddt": 0.9}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"
    res = ingest([zp], out, cat, dry_run=True, require_complete=True,
                 report=True)
    assert res["added"] == 1
    assert res["written"] == 0
    assert not cat.exists()
    assert not out.exists()


# ---------------------------------------------------------------------------
# Accesión duplicada entre lotes: se conserva la de mayor pLDDT
# ---------------------------------------------------------------------------


def test_duplicate_accession_keeps_best_and_reports_it(tmp_path):
    low = _make_zip(
        tmp_path / "a_low.zip", "a",
        {"P1": _pdb_text(50.0)},
        {"P1": {"plddt": 0.550, "ptm": 0.5000, "length": 3}},
    )
    high = _make_zip(
        tmp_path / "b_high.zip", "b",
        {"P1": _pdb_text(95.0)},
        {"P1": {"plddt": 0.950, "ptm": 0.9000, "length": 3}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"

    res = ingest([low, high], out, cat, require_complete=True, report=True)
    assert res["added"] == 1
    assert res["improved"] == 1
    assert res["duplicate"] == 0
    assert res["duplicate_accessions"] >= 1

    row = _read_catalog(cat)[0]
    assert row["plddt"] == "0.950"
    assert row["source_zip"] == "b_high.zip"
    # El lote desplazado queda registrado; no se borra nada.
    assert row["duplicate_of"] == "a_low.zip"

    dup = [r for r in res["report_rows"]
           if r["kind"] == "duplicate_accession"]
    assert len(dup) == 1
    assert dup[0]["accession"] == "P1"
    assert dup[0]["related"] == "b_high.zip"
    assert "P1.pdb" in [p.name for p in out.glob("*.pdb")]

    # La copia de peor calidad nunca sobrescribe a la mejor.
    res3 = ingest([low], out, cat, require_complete=True, report=True)
    assert res3["duplicate"] == 1
    assert _read_catalog(cat)[0]["plddt"] == "0.950"


def test_lower_duplicate_after_best_is_traced_in_duplicate_of(tmp_path):
    high = _make_zip(
        tmp_path / "a_high.zip", "a",
        {"P1": _pdb_text(90.0)},
        {"P1": {"plddt": 0.900}},
    )
    low = _make_zip(
        tmp_path / "b_low.zip", "b",
        {"P1": _pdb_text(40.0)},
        {"P1": {"plddt": 0.400}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"
    res = ingest([high, low], out, cat, require_complete=True, report=True)
    assert res["added"] == 1
    assert res["improved"] == 0
    assert res["duplicate"] == 1
    row = _read_catalog(cat)[0]
    assert row["plddt"] == "0.900"
    assert row["source_zip"] == "a_high.zip"
    assert row["duplicate_of"] == "b_low.zip"


def test_duplicate_hash_between_accessions_is_reported(tmp_path):
    same_pdb = _pdb_text(90.0)
    zp = _make_zip(
        tmp_path / "hash.zip", "h",
        {"PX": same_pdb, "PY": same_pdb},
        {"PX": {"plddt": 0.9}, "PY": {"plddt": 0.8}},
    )
    res = ingest([zp], tmp_path / "out", tmp_path / "cat.tsv",
                 require_complete=True, report=True)
    assert res["added"] == 2
    assert res["duplicate_hashes"] == 1
    assert any(r["kind"] == "duplicate_hash" and r["accession"] == "PY"
               for r in res["report_rows"])


# ---------------------------------------------------------------------------
# ZIP incompleto/corrupto: se reporta y no se ingiere
# ---------------------------------------------------------------------------


def test_truncated_pdb_skips_whole_zip(tmp_path):
    zp = _make_zip(
        tmp_path / "truncado.zip", "t",
        {"P1": _pdb_text(80.0), "P2": _pdb_text(70.0, terminate=False)},
        {"P1": {"plddt": 0.9}, "P2": {"plddt": 0.8}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "cat.tsv"
    res = ingest([zp], out, cat, require_complete=True, report=True)
    assert res["added"] == 0
    assert res["incomplete_zips"] == 1
    assert res["skipped"] == 2
    assert not cat.exists()
    assert not out.exists()
    assert any(r["kind"] == "corrupt_pdb" and r["accession"] == "P2"
               for r in res["report_rows"])
    assert any(r["kind"] == "incomplete_zip" for r in res["report_rows"])


def test_zip_without_metrics_is_not_ingested(tmp_path):
    zp = _make_zip(tmp_path / "sin_metrics.zip", "s",
                   {"P1": _pdb_text(80.0)}, metrics=None)
    res = ingest([zp], tmp_path / "out", tmp_path / "cat.tsv",
                 require_complete=True, report=True)
    assert res["added"] == 0
    assert res["incomplete_zips"] == 1
    assert res["skipped"] == 1
    kinds = {r["kind"] for r in res["report_rows"]}
    assert "incomplete_zip" in kinds
    assert "missing_metrics" in kinds
    # Sin metrics.tsv tampoco se escribe catálogo.
    assert not (tmp_path / "cat.tsv").exists()


def test_pdb_without_metrics_row_is_not_ingested(tmp_path):
    zp = _make_zip(
        tmp_path / "sin_fila.zip", "f",
        {"P1": _pdb_text(80.0), "P2": _pdb_text(70.0)},
        {"P1": {"plddt": 0.9}},  # metrics solo para P1
    )
    res = ingest([zp], tmp_path / "out", tmp_path / "cat.tsv",
                 require_complete=True, report=True)
    assert res["added"] == 0
    assert res["incomplete_zips"] == 1
    assert any(r["kind"] == "missing_metrics" and r["accession"] == "P2"
               for r in res["report_rows"])


def test_metrics_error_status_is_not_ingested(tmp_path):
    zp = _make_zip(
        tmp_path / "fallo.zip", "f",
        {"P1": _pdb_text(80.0), "P2": _pdb_text(70.0)},
        {"P1": {"plddt": 0.9},
         "P2": {"plddt": 0.0, "status": "error"}},
    )
    res = ingest([zp], tmp_path / "out", tmp_path / "cat.tsv",
                 require_complete=True, report=True)
    assert res["added"] == 1
    assert res["skipped"] == 1
    assert res["incomplete_zips"] == 0
    assert any(r["kind"] == "metrics_status" and r["accession"] == "P2"
               for r in res["report_rows"])


def test_incomplete_zip_can_be_ingested_with_legacy_opt_in(tmp_path):
    """Compatibilidad: ``require_complete=False`` mantiene la ingesta vieja."""
    zp = _make_zip(tmp_path / "legacy.zip", "l",
                   {"P1": _pdb_text(80.0)}, metrics=None)
    res = ingest([zp], tmp_path / "out", tmp_path / "cat.tsv",
                 require_complete=False)
    assert res["added"] == 1
    assert res["incomplete_zips"] == 0


# ---------------------------------------------------------------------------
# Auditoría: ZIP nuevos, PDB huérfanos y PDB ausentes (--prune-missing)
# ---------------------------------------------------------------------------


def test_report_and_prune_missing_detect_without_deleting(tmp_path):
    base = _make_zip(
        tmp_path / "base.zip", "b",
        {"P1": _pdb_text(80.0), "P2": _pdb_text(70.0)},
        {"P1": {"plddt": 0.9}, "P2": {"plddt": 0.8}},
    )
    out = tmp_path / "out"
    cat = tmp_path / "cat.tsv"
    ingest([base], out, cat, require_complete=True)

    # Un PDB desaparece y aparece un PDB que ningún catálogo registra.
    (out / "P2.pdb").unlink()
    (out / "ORPH.pdb").write_text(_pdb_text(60.0), encoding="utf-8")
    new = _make_zip(
        tmp_path / "z_nuevo.zip", "z",
        {"P3": _pdb_text(85.0)}, {"P3": {"plddt": 0.85}},
    )

    res = ingest([new], out, cat, require_complete=True, report=True,
                 prune_missing=True)
    assert res["added"] == 1
    assert res["new_zips"] == 1
    assert res["orphan_pdbs"] == 1
    assert res["missing_files"] == 1
    assert res["pruned_missing"] == 1
    assert any(r["kind"] == "new_zip" and r["source_zip"] == "z_nuevo.zip"
               for r in res["report_rows"])
    assert any(r["kind"] == "orphan_pdb" and r["accession"] == "ORPH"
               for r in res["report_rows"])

    rows = {r["accession"]: r for r in _read_catalog(cat)}
    assert set(rows) == {"P1", "P2", "P3"}       # no se borró la fila
    assert rows["P2"]["status"] == "missing"
    assert rows["P2"]["plddt"] == "0.800"        # metadatos conservados
    assert rows["P2"]["source_zip"] == "base.zip"
    assert not (out / "P2.pdb").exists()         # no se inventó el archivo
    assert (out / "ORPH.pdb").exists()           # no se borra el huérfano

    # Segunda pasada: mismo diagnóstico, sin cambios de catálogo.
    before = cat.read_bytes()
    res2 = ingest([new], out, cat, require_complete=True, report=True,
                  prune_missing=True)
    assert res2["pruned_missing"] == 0
    assert res2["missing_files"] == 1
    assert res2["orphan_pdbs"] == 1
    assert cat.read_bytes() == before
    assert len(_read_catalog(cat)) == 3


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_report_tsv_is_written_even_in_dry_run(tmp_path):
    zp = _make_zip(tmp_path / "r.zip", "r",
                   {"P1": _pdb_text(80.0)}, {"P1": {"plddt": 0.9}})
    cat = tmp_path / "cat.tsv"
    tsv = tmp_path / "report.tsv"
    res = ingest([zp], tmp_path / "out", cat, dry_run=True,
                 require_complete=True, report_tsv=tsv)
    assert res["added"] == 1
    assert not cat.exists()
    with tsv.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert list(rows[0]) == REPORT_COLUMNS
    assert rows[0]["kind"] == "new_zip"
    assert rows[0]["source_zip"] == "r.zip"


def test_cli_dry_run_report_writes_nothing(tmp_path, monkeypatch, capsys):
    zp = _make_zip(tmp_path / "cli.zip", "c",
                   {"P1": _pdb_text(80.0)}, {"P1": {"plddt": 0.9}})
    out = tmp_path / "out"
    cat = tmp_path / "cat.tsv"
    monkeypatch.setattr(sys, "argv", [
        "ingest_kaggle_zips.py", "--zips", str(zp), "--outdir", str(out),
        "--catalog", str(cat), "--dry-run", "--report", "--prune-missing",
    ])
    assert main() == 0
    assert not cat.exists() and not out.exists()
    captured = capsys.readouterr().out
    assert "reporte de reingesta" in captured
    assert "total catálogo" in captured


def test_cli_returns_nonzero_for_incomplete_zip(tmp_path, monkeypatch):
    zp = _make_zip(tmp_path / "corrupto.zip", "x",
                   {"P1": _pdb_text(80.0)}, metrics=None)
    monkeypatch.setattr(sys, "argv", [
        "ingest_kaggle_zips.py", "--zips", str(zp),
        "--outdir", str(tmp_path / "out"),
        "--catalog", str(tmp_path / "cat.tsv"), "--dry-run",
    ])
    assert main() == 1
