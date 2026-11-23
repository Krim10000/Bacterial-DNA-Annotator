"""Membrana en pipeline nuevo: contrato de disponibilidad."""
from __future__ import annotations

from pathlib import Path

from webapp.membrane_layer import available


def test_missing_panel_is_unavailable(tmp_path: Path):
    ok, why = available(tmp_path)
    assert not ok
    assert "pfam" in why.lower() or "hmmscan" in why.lower()
