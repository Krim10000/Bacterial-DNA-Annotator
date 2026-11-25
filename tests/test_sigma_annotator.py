"""Pruebas del módulo sigma: scanner, confianza, esquema y contrato GUI."""
from __future__ import annotations

from pathlib import Path

import pytest

from scripts.sigma import sigma_annotator as S

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def test_selftest():
    assert "OK" in S.selftest()


def test_pwm_library_has_architectures():
    pwms = S.load_pwm_library(DATA)
    assert "sigma70" in pwms and "sigma54" in pwms
    assert pwms["sigma54"]["arch"] == "-24/-12"
    assert pwms["sigma70"]["arch"] == "-35/-10"
    assert pwms["sigma70"]["up"]["width"] == 6


def test_confidence_tiers():
    high = S._confidence(0.001, 4.0, 2.5, True, True, "not_available",
                         0.05, 2.0, 1.0)
    assert high[0] == "high" and high[1] == "sospecha_fuerte"
    medium = S._confidence(0.04, 2.2, 1.6, True, False, "not_available",
                           0.05, 2.0, 1.0)
    assert medium[0] == "medium" and medium[1] == "sospecha_media"
    low = S._confidence(0.2, 0.5, 0.2, False, False, "not_available",
                        0.05, 2.0, 1.0)
    assert low[0] == "low" and low[1] == "sin_sospecha"


def test_module_contract():
    from modules import get

    mod = get("sigma")
    assert mod is not None
    assert mod.input_kind == "genome"
    assert mod.run_fn is not None
    assert "OK" in mod.selftest()
    assert "window" in mod.defaults()


@pytest.mark.skipif(not (DATA / "annotations" / "GCF_000006685.1_ASM668v1_genomic.gbff.gz").exists(),
                    reason="sin GBFF de respaldo")
def test_sequence_fallback_for_genome_without_seq_source():
    seqs = S.load_sequence(DATA, "GCF_000006685.1_ASM668v1")
    assert seqs
    assert any(len(v) > 1000 for v in seqs.values())


@pytest.mark.skipif(not (DATA / "regulation" / "ecoli_k12" / "genes.tsv").exists(),
                    reason="sin datos de K-12")
@pytest.mark.skipif(not (DATA / "regulation" / "ecoli_k12" / "genes.tsv").exists(),
                    reason="sin datos de K-12")
def test_module_run_has_preview_and_outputs(tmp_path: Path):
    from modules import get
    from modules.base import RunContext

    mod = get("sigma")
    cfg = mod.defaults()
    cfg.update({"window": 120, "null_windows": 20})
    res = mod.run(RunContext(data_root=DATA, out_root=tmp_path,
                             config=cfg), "ecoli_k12")
    assert res.ok, res.error
    assert "preview" in res.outputs and res.outputs["preview"].exists()
    assert "promoters" in res.outputs and res.outputs["promoters"].exists()
    assert "factors" in res.outputs and res.outputs["factors"].exists()
    assert res.metrics["n_operons"] > 0


def test_k12_schema_and_evidence(tmp_path: Path):
    res = S.annotate_genome(DATA, "ecoli_k12", tmp_path / "sigma",
                            window=120, null_windows=40, seed=3)
    assert res["promoters"].exists()
    assert res["factors"].exists()
    assert res["report"].exists()
    import csv

    with open(res["promoters"], encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    for col in ("sigma", "sigma_status", "confidence_label", "evidence",
                "detected_by", "pipeline_step", "alternatives", "provenance"):
        assert col in rows[0]
    with open(res["factors"], encoding="utf-8") as fh:
        factors = list(csv.DictReader(fh, delimiter="\t"))
    assert factors
    assert factors[0]["sigma_type"]
    assert factors[0]["locus_tag"]
    assert "countable" in factors[0]
    assert "quality_flag" in factors[0]
