"""Contrato del centro de descargas.

Cubre lo que no necesita Streamlit: forma del catálogo, rutas que existen,
filtrado de archivos ausentes, selección de reportes por módulo y tamaños.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import downloads_ui as dl

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_KEYS = {"title", "path", "description", "license", "evidence", "modules"}

#: Los ocho informes de auditoría que deben estar en el catálogo si existen
#: en disco, en cualquiera de sus dos formatos.
REQUESTED_REPORTS = (
    "annotation_parity",
    "mobile_coverage",
    "gap_closure_audit",
    "layer_benchmark",
    "external_benchmark",
    "complex_holdout",
    "model_clusters_calibration",
    "operon_benchmark",
)


def test_catalog_contract_and_no_duplicates():
    assert dl.CATALOG
    paths = [entry["path"] for entry in dl.CATALOG]
    titles = [entry["title"] for entry in dl.CATALOG]
    assert len(paths) == len(set(paths)), "rutas duplicadas en CATALOG"
    assert len(titles) == len(set(titles)), "títulos duplicados en CATALOG"

    for entry in dl.CATALOG:
        assert REQUIRED_KEYS <= set(entry), entry
        assert isinstance(entry["title"], str) and entry["title"].strip()
        assert isinstance(entry["description"], str) and entry["description"].strip()
        assert isinstance(entry["license"], str) and entry["license"].strip()
        assert entry["evidence"] in dl.EVIDENCE_LEVELS
        assert isinstance(entry["modules"], tuple)
        assert all(isinstance(module, str) and module for module in entry["modules"])
        path = Path(entry["path"])
        assert not path.is_absolute()
        assert ".." not in path.parts


def test_every_catalog_path_exists_and_is_a_file():
    # Un checkout público puede no incluir data/; se valida lo que sí existe.
    for entry in dl.catalog_entries(ROOT):
        assert Path(entry["abs_path"]).is_file(), f"archivo inventado o ausente: {entry['path']}"


def test_requested_reports_listed_in_available_formats():
    """Cada base pedida lista exactamente los formatos .md/.json que existen."""
    catalog_paths = {entry["path"] for entry in dl.CATALOG}
    for base in REQUESTED_REPORTS:
        on_disk = {
            f"reports/{base}{suffix}"
            for suffix in (".md", ".json")
            if (ROOT / "reports" / f"{base}{suffix}").is_file()
        }
        assert on_disk, f"no hay ningún informe para {base}"
        listed = {
            path for path in catalog_paths
            if Path(path).stem == base
            and Path(path).suffix in {".md", ".json"}
        }
        assert listed == on_disk, f"{base}: catálogo {listed} != disco {on_disk}"


def test_catalog_entries_resolve_existing_files_and_sizes():
    entries = dl.catalog_entries(ROOT)
    assert len(entries) <= len(dl.CATALOG)
    for entry in entries:
        path = ROOT / entry["path"]
        assert entry["exists"] is True
        assert entry["size_bytes"] == path.stat().st_size > 0
        assert entry["size_human"] == dl._human(path.stat().st_size)
        assert entry["abs_path"] == str(path)


def test_missing_entry_is_filtered_by_default(monkeypatch):
    fake = dict(dl.CATALOG[0])
    fake.update({
        "title": "Reporte ausente de prueba",
        "path": "reports/__reporte_que_no_existe__.json",
        "description": "No debe aparecer en el catálogo disponible.",
        "license": "Prueba",
        "evidence": "not_evaluable",
        "modules": ("annotation",),
    })
    monkeypatch.setattr(dl, "CATALOG", [*dl.CATALOG, fake])

    available = {entry["path"] for entry in dl.catalog_entries(ROOT)}
    assert fake["path"] not in available

    all_entries = dl.catalog_entries(ROOT, include_missing=True)
    missing = [entry for entry in all_entries if entry["path"] == fake["path"]]
    assert len(missing) == 1
    assert missing[0]["exists"] is False
    assert missing[0]["size_bytes"] is None
    assert missing[0]["size_human"] == "—"
    assert missing[0]["rows"] is None


@pytest.mark.parametrize(
    ("module", "expected"),
    [
        ("annotation", {
            "reports/annotation_parity.md",
            "reports/annotation_parity.json",
        }),
        ("mobile", {
            "reports/mobile_coverage.md",
            "reports/mobile_coverage.json",
        }),
        ("complex", {
            "reports/complex_holdout.md",
            "reports/complex_holdout.json",
            "reports/model_clusters_calibration.md",
            "reports/model_clusters_calibration.json",
        }),
        ("gap", {
            "reports/gap_closure_audit.md",
            "reports/gap_closure_audit.json",
        }),
        ("operon", {
            "reports/operon_benchmark.md",
            "reports/operon_benchmark.json",
        }),
        ("layers", {
            "reports/layer_benchmark.md",
            "reports/layer_benchmark.json",
        }),
        ("membrane", {
            "reports/layer_benchmark.md",
            "reports/layer_benchmark.json",
        }),
    ],
)
def test_module_reports_known_modules(module, expected):
    reports = dl.module_reports(module, ROOT)
    assert reports
    paths = {report["path"] for report in reports}
    assert expected <= paths
    canonical = dl.MODULE_ALIASES[module]
    assert all(report["module"] == canonical for report in reports)
    assert all(canonical in report["modules"] for report in reports)
    assert all(report["exists"] is True for report in reports)


def test_module_reports_accepts_aliases_and_case():
    assert dl.module_reports("MOBILE", ROOT) == dl.module_reports("mobile", ROOT)
    assert dl.module_reports("complejos", ROOT) == dl.module_reports("complex", ROOT)
    assert dl.module_reports("gap_closure", ROOT) == dl.module_reports("metabolism", ROOT)


def test_module_reports_unknown_module_is_empty():
    for module in ("no-existe", "", "annotationn"):
        assert dl.module_reports(module, ROOT) == []


def test_module_reports_known_module_without_reports_is_empty():
    # Son ids reales de modules/ que no tienen informe propio en reports/.
    for module in ("taxonomy", "ghost", "targets"):
        assert module in dl.MODULE_ALIASES
        assert dl.module_reports(module, ROOT) == []


def test_module_reports_filters_missing_under_other_root(tmp_path):
    assert dl.module_reports("annotation", tmp_path) == []
    missing = dl.module_reports("annotation", tmp_path, include_missing=True)
    assert missing
    assert all(report["exists"] is False for report in missing)
    assert all(report["size_bytes"] is None for report in missing)
    assert all(report["size_human"] == "—" for report in missing)


def test_module_report_sizes_match_disk():
    for module in ("annotation", "mobile", "complex"):
        for report in dl.module_reports(module, ROOT):
            path = ROOT / report["path"]
            assert report["size_bytes"] == path.stat().st_size > 0
            assert report["size_human"] == dl._human(path.stat().st_size)


@pytest.mark.parametrize(
    ("size", "expected"),
    [
        (0, "0 B"),
        (1, "1 B"),
        (1023, "1023 B"),
        (1024, "1.0 KB"),
        (1536, "1.5 KB"),
        (1024 ** 2, "1.0 MB"),
        (1024 ** 3, "1.0 GB"),
        (5 * 1024 ** 3, "5.0 GB"),
    ],
)
def test_human_sizes(size, expected):
    assert dl._human(size) == expected


def test_rows_only_for_tabular_files():
    tabular = next((entry for entry in dl.catalog_entries(ROOT)
                    if entry["path"].endswith(".tsv")), None)
    if tabular is None:
        return
    report = dl.module_reports("annotation", ROOT)[0]
    assert isinstance(tabular["rows"], int) and tabular["rows"] >= 0
    assert report["rows"] is None
