"""La demo incluida debe existir y cargar sin datos locales."""
from __future__ import annotations

from pathlib import Path

from browser import DEMO_DIR, load_tables


def test_demo_files_exist():
    for name in ("genes.tsv", "operons.tsv", "promoters.tsv", "hairpins.tsv"):
        assert (DEMO_DIR / name).exists(), name


def test_demo_loads_without_project_data(tmp_path: Path):
    tables = load_tables(tmp_path, "demo")
    assert len(tables["genes"]) > 50
    assert len(tables["operons"]) > 1
    assert "locus_tag" in tables["genes"].columns
    assert "operon_id" in tables["operons"].columns


def test_demo_operons_are_marked_low_confidence():
    tables = load_tables(Path("/nonexistent"), "demo")
    assert set(tables["operons"]["confidence"]) == {"low"}
    assert all("demo" == s for s in tables["operons"]["status"])
