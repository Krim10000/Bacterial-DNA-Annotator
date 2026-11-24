"""Contratos de ``structural_b.evaluate``: filtros B sin red por defecto.

Las pruebas no ejecutan Foldseek, EnzyMM, P2Rank ni AlphaFill reales. La rama
sin base local se verifica con dobles que fallan si alguien intenta usarlos.
"""
from __future__ import annotations

import csv
import subprocess
import sys
import textwrap
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

    result = sb.evaluate(rows, outdir, "g", data, foldseek_db=db, cpu=3,
                         external=True)

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


def test_evaluate_external_false_uses_only_local_sources(
        tmp_path: Path, monkeypatch):
    _pdb(tmp_path, "P1")
    rows = _rows("P1")
    outdir = tmp_path / "out"
    data = tmp_path / "data"
    db = tmp_path / "local_db"
    db.write_text("", encoding="utf-8")
    calls: list[str] = []

    def fake_batch(pdbs, target_db, workdir, threads=8, max_seqs=5):
        calls.append("batch")
        return {
            "hits": {"P1": [{"pdb_id": "1abc", "bits": "123"}]},
            "error": "",
        }

    monkeypatch.setattr(fsl, "batch_search", fake_batch)
    monkeypatch.setattr(
        sb, "_enzymm_local",
        lambda pdb, workdir, cpu: (True, "", 2, "9xyz"),
    )
    monkeypatch.setattr(
        sb, "_active_site_local",
        lambda data_root, genome, uniprot: (True, "", [42]),
    )
    # Ninguna de estas rutas debe ejecutarse en modo external=False.
    _forbid(monkeypatch, sb, "_run_filters", calls)
    _forbid(monkeypatch, sf, "run_filters", calls)
    _forbid(monkeypatch, sf, "p2rank", calls)
    _forbid(monkeypatch, sf, "alphafill", calls)
    _forbid(monkeypatch, sf, "foldseek_search", calls)
    _forbid(monkeypatch, requests, "get", calls)
    _forbid(monkeypatch, requests, "post", calls)

    result = sb.evaluate(rows, outdir, "g", data, foldseek_db=db,
                         external=False)

    assert calls == ["batch"]
    assert result["external"] is False
    assert result["network"] == "local_only"
    assert result["n_evaluated"] == 1
    assert result["n_evaluated_local_only"] == 1
    assert result["n_evaluated_delegated"] == 0

    row = rows[0]
    assert row["filter_status"] == "evaluated_local_no_external"
    assert row["filter_f1_foldseek"] is True
    assert row["filter_f2_enzymm"] is True
    assert row["filter_f3_pocket"] == ""
    assert row["filter_f4_alphafill"] == ""
    assert row["filter_f5_active_site"] is True
    assert row["filter_verdict"] == "local_only_partial"
    assert "sin red" in row["filter_limitation"].lower()
    assert "P2Rank/AlphaFill/APIs=no_consultados" in row["filter_note"]
    assert "active_positions=42" in row["filter_note"]

    header, records = _read_tsv(Path(result["filters"]))
    assert header == list(sb.FILTER_COLUMNS)
    assert records[0]["filter_status"] == "evaluated_local_no_external"
    assert records[0]["filter_f3_pocket"] == ""
    assert records[0]["filter_f4_alphafill"] == ""
    assert records[0]["filter_f5_active_site"] == "true"


def test_evaluate_external_false_leaves_unavailable_f2_f5_empty(
        tmp_path: Path, monkeypatch):
    _pdb(tmp_path, "P1")
    rows = _rows("P1")
    db = tmp_path / "local_db"
    db.write_text("", encoding="utf-8")
    calls: list[str] = []

    monkeypatch.setattr(
        fsl, "batch_search",
        lambda *a, **k: {"hits": {}, "error": ""},
    )
    monkeypatch.setattr(
        sb, "_enzymm_local",
        lambda pdb, workdir, cpu: (None, "EnzyMM no instalado", 0, ""),
    )
    monkeypatch.setattr(
        sb, "_active_site_local",
        lambda data_root, genome, uniprot: (None, "mapa ausente", []),
    )
    _forbid(monkeypatch, sb, "_run_filters", calls)
    _forbid(monkeypatch, sf, "p2rank", calls)
    _forbid(monkeypatch, sf, "alphafill", calls)

    result = sb.evaluate(rows, tmp_path / "out", "g", tmp_path / "data",
                         foldseek_db=db, external=False)

    assert calls == []
    assert result["n_evaluated_local_only"] == 1
    row = rows[0]
    assert row["filter_status"] == "evaluated_local_no_external"
    assert row["filter_f1_foldseek"] is False
    assert row["filter_f2_enzymm"] == ""
    assert row["filter_f3_pocket"] == ""
    assert row["filter_f4_alphafill"] == ""
    assert row["filter_f5_active_site"] == ""
    assert "F2 no disponible" in row["filter_limitation"]
    assert "F5 no disponible" in row["filter_limitation"]
    assert "f2_no_evaluado=EnzyMM no instalado" in row["filter_note"]
    assert "f5_no_evaluado=mapa ausente" in row["filter_note"]


def test_active_site_local_reads_data_root(tmp_path: Path):
    base = tmp_path / "data" / "metabolism" / "targets"
    base.mkdir(parents=True)
    (base / "ACTIVE_SITE_STRUCTURAL_MAP.tsv").write_text(
        "genoma\tacc\tfeature\tdescription\tstart\tend\tref_pos\tref_aa\t"
        "qry_pos\tqry_aa\tconserved\n"
        "g\tP1\tActive site\tdesc\t10\t10\t10\tA\t42\tA\tTrue\n"
        "g\tP1\tBinding site\tdesc\t11\t11\t11\tA\t43\tA\tTrue\n",
        encoding="utf-8",
    )

    available, note, positions = sb._active_site_local(
        tmp_path / "data", "g", "P1")

    assert available is True, note
    assert positions == [42]


def test_external_false_does_not_import_api_modules(tmp_path: Path):
    """Verifica en un intérprete limpio que no se importan módulos de API."""
    _pdb(tmp_path, "P1")
    (tmp_path / "db").write_text("", encoding="utf-8")
    root = Path(__file__).resolve().parents[1]
    script = textwrap.dedent("""
        import sys
        from pathlib import Path

        import webapp.foldseek_local as fsl
        import webapp.structural_b as sb

        BLOCKED = {"webapp.structural_filters", "structure_evidence", "requests"}

        class BlockAPI:
            def find_spec(self, fullname, path=None, target=None):
                if fullname in BLOCKED:
                    raise AssertionError(f"import prohibido: {fullname}")
                return None

        sys.meta_path.insert(0, BlockAPI())
        root = Path(sys.argv[1])

        def batch(pdbs, target_db, workdir, threads=8, max_seqs=5):
            return {"hits": {"P1": [{"pdb_id": "1abc"}]}, "error": ""}

        fsl.batch_search = batch
        sb._enzymm_local = lambda *a, **k: (False, "", 0, "")
        sb._active_site_local = lambda *a, **k: (False, "", [])
        rows = [{"protein_id": "p1", "path": "B", "uniprot": "P1"}]
        sb.evaluate(rows, root / "out", "g", root / "data",
                    foldseek_db=root / "db", external=False)
        assert rows[0]["filter_status"] == "evaluated_local_no_external"
        assert rows[0]["filter_f3_pocket"] == ""
        assert rows[0]["filter_f4_alphafill"] == ""
        assert "webapp.structural_filters" not in sys.modules
        assert "structure_evidence" not in sys.modules
        assert "requests" not in sys.modules
        print("ok")
    """)
    proc = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path)],
        cwd=root, capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
    assert "ok" in proc.stdout


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
                         foldseek_db=db, max_models=2, external=True)

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
