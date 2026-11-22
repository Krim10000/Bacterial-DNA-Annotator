"""Prueba de humo de la anotación: un FASTA de demostración debe producir
features, cobertura y un GFF no vacío, sin red ni GPU."""
from __future__ import annotations

from pathlib import Path

import pytest

from annotator import annotate

SAMPLE = Path(__file__).resolve().parents[1] / "webapp" / "samples" / "demo.fna"


@pytest.mark.skipif(not SAMPLE.exists(), reason="falta el FASTA de demostración")
def test_annotate_demo_produces_features():
    pytest.importorskip("pyrodigal")
    res = annotate("demo.fna", SAMPLE.read_bytes(),
                   {"cds": True, "ghost": False, "trna": False, "rrna": False,
                    "esm": False, "deserts": True, "min_region_len": 150,
                    "min_aa": 60})
    assert res["bp"] > 0
    assert res["stats"].get("cds", 0) > 10
    assert res["coverage"] > 0.5
    assert res["gff"].strip()
    assert "##gff-version" in res["gff"] or res["gff"].startswith("##")
    assert res["report_md"]
