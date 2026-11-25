"""Pruebas sintéticas de la auditoría AFDB/PAE.

No usan red ni tocan la caché real: el catálogo, el universo y los archivos
locales se crean en ``tmp_path``. El doble de ``urlopen`` falla si algo intenta
conectarse.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from scripts.structure import afdb_status as status
from scripts.structure import fetch_afdb_cache as cache

PDB_TEXT = (
    "HEADER    TEST PROTEIN\n"
    "SEQRES   1 A    4  MET ALA LEU GLY\n"
    "ATOM      1  CA  MET A     1       0.000   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)


def _write_universe(data_root: Path, accessions: list[str]) -> None:
    folder = data_root / "regulation" / "g1"
    folder.mkdir(parents=True)
    with (folder / "genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome", "gene_id", "acc"])
        for accession in accessions:
            writer.writerow(["g1", f"gene_{accession}", accession])


def _write_catalog(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cache.CAT_COLS, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "")
                             for column in cache.CAT_COLS})


def _case(tmp_path: Path) -> dict:
    """Caso sintético: 20 accesiones, 9 filas de catálogo y 4 PDB locales.

    - A0A00001: ok con PDB + PAE.
    - A0A00002: error 429 sin archivo (reintentable).
    - A0A00003: error 429 con PDB ``.gz`` (cubierta, no pendiente).
    - A0A00004: error 404 (permanente).
    - A0A00005: not_found (permanente).
    - A0A00006: timeout (error reintentable).
    - A0A00007: rate_limited sin archivo.
    - A0A00008: PDB local sin fila de catálogo.
    - A0A00009: ok sin archivo (integridad: se vuelve error).
    - A0A00010: already_present con PDB local.
    - A0A00011..A0A00020: nunca intentadas.
    """
    universe = [f"A0A000{index:02d}" for index in range(1, 21)]
    data_root = tmp_path / "data"
    targets = tmp_path / "targets"
    targets.mkdir()
    _write_universe(data_root, universe)

    (targets / "A0A00001.pdb").write_text(PDB_TEXT, encoding="utf-8")
    (targets / "A0A00001.pae.json").write_text("{}", encoding="utf-8")
    with gzip.open(targets / "A0A00003.pdb.gz", "wb") as fh:
        fh.write(PDB_TEXT.encode("ascii"))
    (targets / "A0A00008.pdb").write_text(PDB_TEXT, encoding="utf-8")
    (targets / "A0A00010.pdb").write_text(PDB_TEXT, encoding="utf-8")

    catalog = tmp_path / "afdb_catalog.tsv"
    rows = [
        {"accession": "A0A00001", "status": "ok",
         "pdb": str(targets / "A0A00001.pdb"),
         "pae": str(targets / "A0A00001.pae.json")},
        {"accession": "A0A00002", "status": "error",
         "error": "HTTPError: HTTP Error 429: Too Many Requests"},
        {"accession": "A0A00003", "status": "error",
         "error": "HTTPError: HTTP Error 429: Too Many Requests"},
        {"accession": "A0A00004", "status": "error",
         "error": "HTTPError: HTTP 404: Not Found"},
        {"accession": "A0A00005", "status": "not_found",
         "error": "HTTP 404: Not Found"},
        {"accession": "A0A00006", "status": "error",
         "error": "TimeoutError: timed out"},
        {"accession": "A0A00007", "status": "rate_limited", "error": "429"},
        {"accession": "A0A00009", "status": "ok"},
        {"accession": "A0A00010", "status": "already_present",
         "pdb": str(targets / "A0A00010.pdb")},
    ]
    _write_catalog(catalog, rows)
    return {
        "universe": universe,
        "data_root": data_root,
        "targets": targets,
        "catalog": catalog,
        "rows": rows,
    }


def _audit(case: dict) -> status.Audit:
    return status.build_audit(
        case["universe"],
        cache.load_catalog(case["catalog"]),
        case["targets"],
        catalog_path=case["catalog"],
    )


# ---------------------------------------------------------------------------
# Clasificación y métricas
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("text,code", [
    ("HTTPError: HTTP Error 429: Too Many Requests", 429),
    ("HTTPError: HTTP 404: Not Found", 404),
    ("HTTP 500: boom", 500),
    ("TimeoutError: timed out", None),
    ("", None),
])
def test_extract_http_code(text: str, code: int | None) -> None:
    assert status.extract_http_code(text) == code


def test_classifies_429_as_retryable_and_404_as_permanent(tmp_path: Path) -> None:
    audit = _audit(_case(tmp_path))
    summary = status.summarize(audit)

    assert summary["universe_total"] == 20
    assert summary["catalog_rows"] == 9
    assert (summary["catalog_ok"],
            summary["catalog_rate_limited"],
            summary["catalog_not_found"],
            summary["catalog_error"]) == (2, 3, 2, 2)

    # La fila bruta ``status=error`` se reparte por código HTTP.
    assert summary["catalog_raw_error"] == 4
    assert summary["catalog_raw_error_rate_limited"] == 2
    assert summary["catalog_raw_error_not_found"] == 1
    assert summary["catalog_raw_error_other"] == 1

    assert audit.items["A0A00002"].retryable is True
    assert audit.items["A0A00002"].permanent is False
    assert audit.items["A0A00003"].retryable is True  # 429, aunque tenga PDB
    assert audit.items["A0A00004"].permanent is True
    assert audit.items["A0A00004"].retryable is False
    assert audit.items["A0A00005"].permanent is True
    assert audit.items["A0A00006"].retryable is True
    assert audit.items["A0A00007"].retryable is True
    assert audit.items["A0A00009"].status == status.STATUS_ERROR  # ok sin archivo
    assert audit.items["A0A00003"].local_pdb.endswith(".pdb.gz")


def test_coverage_pending_and_local_cache(tmp_path: Path) -> None:
    audit = _audit(_case(tmp_path))
    summary = status.summarize(audit)

    assert summary["local_pdb"] == 4
    assert summary["local_pae"] == 1
    assert summary["local_only"] == 1
    assert summary["catalog_ok_missing_file"] == 1
    assert summary["pae_missing_for_local"] == 3
    assert summary["orphan_pae"] == 0

    assert summary["covered_in_universe"] == 4
    assert summary["permanent_uncovered"] == 2
    assert summary["retryable_pending"] == 4   # A2, A6, A7, A9
    assert summary["never_attempted"] == 10    # A11..A20
    assert summary["pending_total"] == 14      # 10 + 4

    assert summary["catalog_coverage"] == pytest.approx(2 / 20)
    assert summary["disk_coverage"] == pytest.approx(4 / 20)
    assert summary["accounted_coverage"] == pytest.approx(6 / 20)
    assert summary["pae_coverage"] == pytest.approx(1 / 4)
    assert summary["observed_429_rate"] == pytest.approx(3 / 9)
    # El caso sintético usa "HTTP 404" (reconocido por el descargador).
    assert summary["fetcher_gap"] == 0
    assert summary["fetcher_pending_total"] == summary["pending_total"]


def test_next_batch_is_deterministic_and_respects_order(tmp_path: Path) -> None:
    audit = _audit(_case(tmp_path))

    assert status.next_batch(audit, 3) == ["A0A00002", "A0A00006", "A0A00007"]
    assert status.next_batch(audit, 0) == []
    assert status.next_batch(audit, 100) == audit.pending
    assert audit.pending[0] == "A0A00002"
    assert "A0A00003" not in audit.pending  # tiene PDB local
    assert "A0A00005" not in audit.pending  # 404 permanente


def test_scan_local_detects_gz_and_prefers_plain(tmp_path: Path) -> None:
    targets = tmp_path / "targets"
    targets.mkdir()
    (targets / "X1.pdb").write_text("plain", encoding="utf-8")
    (targets / "X1.pdb.gz").write_bytes(gzip.compress(b"gz"))
    (targets / "X2.pae.json.gz").write_bytes(gzip.compress(b"{}"))

    pdb, pae = status.scan_local(targets)

    assert pdb["X1"].name == "X1.pdb"
    assert pdb["X1"].read_text(encoding="utf-8") == "plain"
    assert pae["X2"].name == "X2.pae.json.gz"
    assert status.scan_local(tmp_path / "missing") == ({}, {})


def test_load_universe_matches_downloader_and_deduplicates(tmp_path: Path) -> None:
    data_root = tmp_path / "data"
    _write_universe(data_root, ["A0A00002", "A0A00001", "A0A00002"])
    # Carpeta sin columna ``acc``: se ignora sin romper la auditoría.
    bad = data_root / "regulation" / "g2"
    bad.mkdir()
    (bad / "genes.tsv").write_text("genome\tgene_id\n", encoding="utf-8")

    assert status.load_universe(data_root) == ["A0A00001", "A0A00002"]
    assert status.load_universe(data_root) == cache.accessions(data_root)
    assert status.load_universe(tmp_path / "missing") == []


# ---------------------------------------------------------------------------
# Plan y ETA
# ---------------------------------------------------------------------------
def test_fetcher_gap_is_surfaced_and_limit_is_adjusted(tmp_path: Path) -> None:
    universe = ["A1", "A2", "A3"]
    data_root = tmp_path / "data"
    _write_universe(data_root, universe)
    catalog = tmp_path / "afdb_catalog.tsv"
    _write_catalog(catalog, [
        {"accession": "A1", "status": "error",
         "error": "HTTPError: HTTP Error 429: Too Many Requests"},
        {"accession": "A2", "status": "error",
         "error": "HTTPError: HTTP Error 404: Not Found"},
    ])
    audit = status.build_audit(universe, cache.load_catalog(catalog),
                               tmp_path / "targets", catalog_path=catalog)
    summary = status.summarize(audit)

    # A2 es permanente y el descargador ya la reconoce: no hay brecha.
    assert summary["catalog_not_found"] == 1
    assert summary["fetcher_gap"] == 0
    assert summary["fetcher_gap_items"] == []
    assert summary["pending_total"] == 2          # A1 y A3
    assert summary["fetcher_pending_total"] == 2
    assert status.next_batch(audit, 2) == ["A1", "A3"]
    assert status.recommended_limit(audit, summary, 2) == 2
    plan = status.build_plan(summary, probe_limit=2)
    assert plan["steps"][0]["limit"] == 2
    assert "--limit 2" in plan["steps"][0]["command"]


def test_plan_starts_conservative_when_429_exist(tmp_path: Path) -> None:
    summary = status.summarize(_audit(_case(tmp_path)))
    plan = status.build_plan(summary)

    assert [step["name"] for step in plan["steps"]] == [
        "sonda", "rampa", "estable", "produccion",
    ]
    probe = plan["steps"][0]
    assert probe["workers"] == 1
    assert probe["delay"] == 2.0
    assert probe["limit"] == 10
    assert probe["pae"] is False
    assert "fetch_afdb_cache.py" in probe["command"]
    assert "--pae" not in probe["command"]
    assert probe["eta_seconds"] > 0

    assert plan["sustained_seconds"] > 0
    assert plan["accelerated_seconds"] <= plan["sustained_seconds"]
    assert plan["slow_seconds"] >= plan["sustained_seconds"]
    assert plan["retry_factor"] == pytest.approx(1 + 3 / 9)


def test_plan_skips_probe_without_429() -> None:
    summary = {
        "pending_total": 1000,
        "catalog_rate_limited": 0,
        "observed_429_rate": 0.0,
    }
    plan = status.build_plan(summary)

    assert len(plan["steps"]) == 1
    step = plan["steps"][0]
    assert (step["workers"], step["delay"], step["limit"]) == (2, 0.5, 100)
    assert plan["sustained_seconds"] == plan["accelerated_seconds"]
    assert plan["sustained_seconds"] == plan["slow_seconds"]


def test_plan_with_nothing_pending() -> None:
    plan = status.build_plan({"pending_total": 0})
    assert plan["steps"] == []
    assert plan["sustained_seconds"] == 0


# ---------------------------------------------------------------------------
# CLI / artefactos / solo lectura
# ---------------------------------------------------------------------------
def test_main_dry_run_writes_nothing_and_real_run_writes_reports(
        tmp_path: Path, monkeypatch, capsys) -> None:
    case = _case(tmp_path)
    reports = tmp_path / "reports"
    args = [
        "--data-root", str(case["data_root"]),
        "--catalog", str(case["catalog"]),
        "--targets", str(case["targets"]),
        "--reports-dir", str(reports),
        "--limit-next", "3",
    ]

    def _no_network(*_args, **_kwargs):
        raise AssertionError("la auditoría no debe hacer red")

    monkeypatch.setattr(cache.urllib.request, "urlopen", _no_network)
    catalog_before = case["catalog"].read_bytes()
    targets_before = {path.name: path.read_bytes()
                      for path in case["targets"].iterdir()}

    assert status.main(args + ["--dry-run"]) == 0
    assert not reports.exists()
    output = capsys.readouterr().out
    assert "universo=20" in output
    assert "rate_limited=3" in output
    assert "A0A00002" in output
    assert "dry-run" in output

    assert status.main(args) == 0
    tsv = reports / "afdb_status.tsv"
    md = reports / "afdb_status.md"
    assert tsv.is_file() and md.is_file()

    with tsv.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == status.TSV_COLS
    assert len(rows) == 10  # 9 filas de catálogo + 1 PDB solo local
    assert rows[0]["accession"] == "A0A00002"
    assert rows[0]["status"] == "rate_limited"
    by_acc = {row["accession"]: row for row in rows}
    assert by_acc["A0A00003"]["local_pdb"] == "true"
    assert by_acc["A0A00003"]["pending"] == "false"
    assert by_acc["A0A00004"]["permanent"] == "true"
    assert by_acc["A0A00004"]["retryable"] == "false"
    assert by_acc["A0A00008"]["status"] == "local_only"
    assert by_acc["A0A00009"]["status"] == "error"

    md_text = md.read_text(encoding="utf-8")
    assert "# Estado del catálogo AFDB/PAE" in md_text
    assert "modo solo lectura" in md_text
    assert "A0A00002" in md_text
    assert "python scripts/structure/fetch_afdb_cache.py" in md_text
    assert "no modifica" in md_text

    # Ni el catálogo ni la caché cambian: solo se escriben los reportes.
    assert case["catalog"].read_bytes() == catalog_before
    assert {path.name: path.read_bytes()
            for path in case["targets"].iterdir()} == targets_before
    assert sorted(path.name for path in reports.iterdir()) == [
        "afdb_status.md", "afdb_status.tsv",
    ]


def test_main_returns_2_when_catalog_is_missing(tmp_path: Path) -> None:
    data_root = tmp_path / "data"
    _write_universe(data_root, ["A0A00001"])
    reports = tmp_path / "reports"
    rc = status.main([
        "--data-root", str(data_root),
        "--catalog", str(tmp_path / "missing.tsv"),
        "--targets", str(tmp_path / "targets"),
        "--reports-dir", str(reports),
    ])
    assert rc == 2
    assert not reports.exists()
