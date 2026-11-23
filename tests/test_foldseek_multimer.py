"""Contrato de Foldseek-Multimer sin base de ensamblajes."""
from pathlib import Path

from scripts.complexes.foldseek_multimer import available


def test_requires_database(tmp_path: Path):
    ok, why = available(tmp_path / "no-db")
    assert not ok
    assert "db" in why.lower() or "foldseek" in why.lower()
