"""Contratos de ``structural_b.evaluate``: filtros B sin red por defecto.

Las pruebas no ejecutan Foldseek, EnzyMM, P2Rank ni AlphaFill reales. La rama
sin base local se verifica con dobles que fallan si alguien intenta usarlos.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest
import requests

import webapp.foldseek_local as fsl
import webapp.structural_b as sb
import webapp.structural_filters as sf

AFDB_REL = Path("metabolism") / "targets" / "afdb25k"


def _pdb(tmp_path: Path, acc: str = "P1") -> Path:
    afdb = tmp_path / "data" / AFDB_REL
    afdb.mkdir(parents=True, exist_ok=True)
    pdb = afdb / f"{acc}.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n",
        encoding="utf-8",
    )
    return pdb


def _rows(*accs: str) -> list[dict]:
    return [{"protein_id": f"p{i}", "path": "B", "uniprot": acc}
            for i, acc in enumerate(accs, 1)]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _forbid(monkeypatch, module, name: str, calls: list[str]):
    def forbidden(*args, **kwargs):
        calls.append(f"{module.__name__}.{name}")
        raise AssertionError(f"no debía llamarse {module.__name__}.{name}")

    monkeypatch.setattr(module, name, forbidden)


@pytest.mark.parametrize("db_arg", [None, "missing"])
def test_evaluate_without_local_db_does_not_call_filters_or_network(
        tmp_path: Path, monkeypatch, db_arg):
    _pdb(tmp_path, "P1")
    rows = _rows("P1", "P2")
    rows.append({"protein_id": "p3", "path": "C", "uniprot": "P3"})
    outdir = tmp_path / "out"
    data = tmp_path / "data"

    # Conserva el manifiesto/status existente de run().
    manifest = sb.run(rows, outdir, "g", data)
    assert Path(manifest["manifest"]).exists()
    assert rows[0]["structural_status"] == "afdb_model_available"

    calls: list[str] = []
    _forbid(monkeypatch, fsl, "batch_search", calls)
    _forbid(monkeypatch, sf, "run_filters", calls)
    _forbid(monkeypatch, requests, "get", calls)
    _forbid(monkeypatch, requests, "post", calls)

    db = None if db_arg is None else tmp_path / "missing_db"
    result = sb.evaluate(rows, outdir, "g", data, foldseek_db=db, cpu=2,
                         max_models=25)

    assert calls == []
    assert result["network"] == "not_used"
    assert result["n_evaluated"] == 0
    assert result["n_no_local_db"] == 1
    assert result["n_no_model"] == 1

    # La fila B con modelo no se evalúa, pero tampoco se marca como negativa.
    assert rows[0]["structural_status"] == "afdb_model_available"
    assert rows[0]["filter_status"] == "not_evaluated_no_local_db"
    for field in sb.BOOL_FILTER_FIELDS:
        assert rows[0][field] == ""
    assert rows[0]["filter_verdict"] == "not_evaluated"
    assert "no se puede afirmar" in rows[0]["filter_limitation"].lower()
    assert "no se ejecutaron" in rows[0]["filter_note"].lower()
    assert "no es evidencia negativa" in rows[0]["filter_note"].lower()

    # La fila B sin modelo tiene estado propio y tampoco inventa F1-F5.
    assert rows[1]["filter_status"] == "no_afdb_model"
    for field in sb.BOOL_FILTER_FIELDS:
        assert rows[1][field] == ""
    assert rows[1]["filter_verdict"] == "not_evaluated"

    # La fila C no se toca.
    assert "filter_status" not in rows[2]

    header, records = _read_tsv(Path(result["filters"]))
    assert header == list(sb.FILTER_COLUMNS)
    assert [r["protein_id"] for r in records] == ["p1", "p2"]
    assert records[0]["filter_status"] == "not_evaluated_no_local_db"
    assert records[1]["filter_status"] == "no_afdb_model"
    assert records[0]["filter_f1_foldseek"] == ""
    assert records[0]["filter_verdict"] == "not_evaluated"
    assert records[0]["filter_note"]


def test_evaluate_with_local_db_uses_local_hits_and_maps_filters(
        tmp_path: Path, monkeypatch):
    _pdb(tmp_path, "P1")
    _pdb(tmp_path, "P2")
    rows = _rows("P1", "P2")
    outdir = tmp_path / "out"
    data = tmp_path / "data"
    db = tmp_path / "local_db"
    db.write_text("", encoding="utf-8")

    calls: dict[str, object] = {}

    def fake_batch(pdbs, target_db, workdir, threads=8, max_seqs=5):
        calls["batch"] = {
            "pdbs": [Path(p).stem for p in pdbs],
            "target_db": str(target_db),
            "threads": threads,
            "max_seqs": max_seqs,
        }
        return {
            "hits": {"P1": [{"pdb_id": "1abc", "bits": "123"}]},
            "error": "",
        }

    def fake_filters(pdb, **kwargs):
        calls.setdefault("filters", []).append(
            {"pdb": Path(pdb).stem, **kwargs}
        )
        return {
            "f1_foldseek": True,
            "f2_enzymm": False,
            "f3_pocket": True,
            "f4_alphafill": False,
            "f5_active_site": True,
            "n_pass": 3,
            "priority": "review",
            "gating": "evaluado",
            "quality_gate": "unknown",
            "foldseek_top": "1abc",
            "source": "local_foldseek_swissprot",
            "limitation": "Filtros computacionales; no confirman función.",
            "error": "",
        }

    monkeypatch.setattr(fsl, "batch_search", fake_batch)
    monkeypatch.setattr(sf, "run_filters", fake_filters)

    result = sb.evaluate(rows, outdir, "g", data, foldseek_db=db, cpu=3)

    batch = calls["batch"]
    assert batch["pdbs"] == ["P1", "P2"]
    assert batch["threads"] == 3
    assert batch["max_seqs"] == 5
    filters = calls["filters"]
    assert [c["pdb"] for c in filters] == ["P1", "P2"]
    assert filters[0]["genome"] == "g" and filters[0]["uniprot"] == "P1"
    assert filters[0]["foldseek_hits"][0]["pdb_id"] == "1abc"
    assert filters[1]["foldseek_hits"] == []

    assert result["n_evaluated"] == 2
    for row in rows:
        assert row["filter_status"] == "evaluated_local_foldseek"
        assert row["filter_verdict"] == "review"
        assert isinstance(row["filter_f1_foldseek"], bool)
        assert "no confirman" in row["filter_limitation"]
    assert rows[0]["filter_f1_foldseek"] is True
    assert rows[0]["filter_f4_alphafill"] is False

    header, records = _read_tsv(Path(result["filters"]))
    assert header == list(sb.FILTER_COLUMNS)
    assert records[0]["filter_f1_foldseek"] == "true"
    assert records[0]["filter_f2_enzymm"] == "false"
    assert records[0]["filter_f5_active_site"] == "true"


def test_evaluate_limits_models_and_marks_overflow_not_evaluated(
        tmp_path: Path, monkeypatch):
    _pdb(tmp_path, "P1")
    _pdb(tmp_path, "P2")
    _pdb(tmp_path, "P3")
    rows = _rows("P1", "P2", "P3")
    db = tmp_path / "db"
    db.write_text("", encoding="utf-8")
    calls: list[str] = []

    def fake_batch(pdbs, target_db, workdir, threads=8, max_seqs=5):
        calls.extend(Path(p).stem for p in pdbs)
        return {"hits": {}, "error": ""}

    def fake_filters(pdb, **kwargs):
        return {
            "f1_foldseek": False, "f2_enzymm": False, "f3_pocket": False,
            "f4_alphafill": False, "f5_active_site": False,
            "n_pass": 0, "priority": "low", "gating": "no_structural_signal",
            "quality_gate": "unknown", "foldseek_top": "",
            "source": "local_foldseek_swissprot",
            "limitation": "Filtros computacionales.", "error": "",
        }

    monkeypatch.setattr(fsl, "batch_search", fake_batch)
    monkeypatch.setattr(sf, "run_filters", fake_filters)

    result = sb.evaluate(rows, tmp_path / "out", "g", tmp_path / "data",
                         foldseek_db=db, max_models=2)

    assert calls == ["P1", "P2"]
    assert result["n_evaluated"] == 2
    assert result["n_skipped_max_models"] == 1
    assert rows[2]["filter_status"] == "not_evaluated_max_models"
    assert rows[2]["filter_verdict"] == "not_evaluated"
    for field in sb.BOOL_FILTER_FIELDS:
        assert rows[2][field] == ""
    assert "max_models" in rows[2]["filter_note"]


def test_evaluate_foldseek_error_does_not_invent_filter_results(
        tmp_path: Path, monkeypatch):
    _pdb(tmp_path, "P1")
    rows = _rows("P1")
    db = tmp_path / "db"
    db.write_text("", encoding="utf-8")
    calls: list[str] = []

    monkeypatch.setattr(
        fsl, "batch_search",
        lambda *a, **k: {"hits": {}, "error": "foldseek no instalado"},
    )
    _forbid(monkeypatch, sf, "run_filters", calls)

    result = sb.evaluate(rows, tmp_path / "out", "g", tmp_path / "data",
                         foldseek_db=db)

    assert calls == []
    assert result["n_evaluated"] == 0
    assert result["n_errors"] == 1
    assert rows[0]["filter_status"] == "not_evaluated_foldseek_error"
    for field in sb.BOOL_FILTER_FIELDS:
        assert rows[0][field] == ""
    assert "foldseek no instalado" in rows[0]["filter_note"]
