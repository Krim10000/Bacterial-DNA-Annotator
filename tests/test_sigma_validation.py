"""Pruebas del benchmark TSS-level de sigma (``30_validate_sigma.py``).

Cubren mapeo de etiquetas, ventana orientada, métricas/confusión/calibración y
un smoke test con datos reales. Los tests de datos se saltan si faltan archivos,
para que la suite no dependa de un checkout completo.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma import sigma_annotator as S  # noqa: E402

V = importlib.import_module("scripts.sigma.30_validate_sigma")
DATA = ROOT / "data"


# --------------------------------------------------------------------------
# Mapeo de etiquetas
# --------------------------------------------------------------------------
def test_map_k12_labels():
    assert V.map_k12_sigma("sigma70")["mapping_status"] == "evaluable"
    assert V.map_k12_sigma("sigma70")["sigma_truth"] == "sigma70"
    assert V.map_k12_sigma("sigma54")["sigma_truth"] == "sigma54"
    assert V.map_k12_sigma("sigma19")["mapping_status"] == "no_pwm"
    assert V.map_k12_sigma("sigma19")["mapping_reason"]
    assert V.map_k12_sigma("")["mapping_status"] == "unlabeled"
    assert V.map_k12_sigma("  ")["mapping_status"] == "unlabeled"


def test_map_bsu_labels():
    # Etiquetas de letra y formato SigX.
    assert V.map_bsu_sigma("A")["sigma_truth"] == "sigma70"
    assert V.map_bsu_sigma("SigA")["sigma_truth"] == "sigma70"
    assert V.map_bsu_sigma("B")["sigma_truth"] == "sigma38"
    assert V.map_bsu_sigma("SigB")["sigma_truth"] == "sigma38"
    assert V.map_bsu_sigma("D")["sigma_truth"] == "sigma28"
    assert V.map_bsu_sigma("I")["sigma_truth"] == "sigma32"
    assert V.map_bsu_sigma("L")["sigma_truth"] == "sigma54"
    # Sigma sin PWM equivalente: no cuentan como fallo.
    for label in ("E", "F", "G", "K", "H", "W", "X", "M", "YlaC"):
        item = V.map_bsu_sigma(label)
        assert item["mapping_status"] == "no_pwm", label
        assert item["sigma_truth"] == ""
        assert item["mapping_reason"]
    assert V.map_bsu_sigma("?")["mapping_status"] == "unlabeled"
    assert V.map_bsu_sigma("")["mapping_status"] == "unlabeled"


# --------------------------------------------------------------------------
# Ventana orientada
# --------------------------------------------------------------------------
def test_tss_window_plus_and_minus():
    seq = "ACGT" * 60  # 240 nt
    pos = 121         # índice 0-based 120
    plus = V.tss_window(seq, pos, "+")
    minus = V.tss_window(seq, pos, "-")
    assert len(plus) == V.TSS_WINDOW == 101
    assert len(minus) == V.TSS_WINDOW
    assert plus[V.TSS_UPSTREAM] == seq[120]
    assert plus[:V.TSS_UPSTREAM] == seq[40:120]
    assert plus[V.TSS_UPSTREAM + 1:] == seq[121:141]
    assert minus == S.revcomp(seq[100:201])
    assert minus[V.TSS_UPSTREAM] == S.revcomp(seq[120:121])
    assert V.tss_window(seq, 0, "+") == ""
    assert V.tss_window(seq, 5, "+") == ""
    assert V.tss_window(seq, 235, "-") == ""
    assert V.tss_window(seq, None, "+") == ""


# --------------------------------------------------------------------------
# Métricas sintéticas
# --------------------------------------------------------------------------
def _rows() -> list[dict]:
    return [
        {"label_status": "evaluable", "mapping_status": "evaluable",
         "sigma_truth": "sigma70", "pred_sigma": "sigma70",
         "correct_top1": 1, "correct_top3": 1, "pred_rank": 1, "pred_z": 4.0,
         "confidence_label": "sospecha_fuerte", "confidence_score": 0.9},
        {"label_status": "evaluable", "mapping_status": "evaluable",
         "sigma_truth": "sigma70", "pred_sigma": "sigma38",
         "correct_top1": 0, "correct_top3": 0, "pred_rank": 5, "pred_z": 1.0,
         "confidence_label": "sin_sospecha", "confidence_score": 0.2},
        {"label_status": "evaluable", "mapping_status": "evaluable",
         "sigma_truth": "sigma38", "pred_sigma": "sigma70",
         "correct_top1": 0, "correct_top3": 1, "pred_rank": 2, "pred_z": 3.0,
         "confidence_label": "sospecha_media", "confidence_score": 0.5},
        {"label_status": "no_pwm", "mapping_status": "no_pwm",
         "sigma_truth": "", "pred_sigma": "", "correct_top1": "",
         "correct_top3": "", "pred_rank": "", "pred_z": "",
         "confidence_label": "", "confidence_score": ""},
    ]


def test_metrics_basic():
    metrics = V.compute_metrics(_rows(), "synthetic")
    assert metrics["status"] == "evaluated"
    assert metrics["counts"]["n_evaluable"] == 3
    assert metrics["counts"]["n_no_pwm"] == 1
    assert metrics["counts"]["n_labeled_valid"] == 4
    assert metrics["coverage"]["evaluable_over_labeled_valid"] == pytest.approx(0.75)
    assert metrics["accuracy"]["top1"] == pytest.approx(1 / 3)
    assert metrics["accuracy"]["top3"] == pytest.approx(2 / 3)
    assert metrics["accuracy"]["majority_class"] == "sigma70"
    assert metrics["confusion_matrix"][0][0] == 1
    labels = metrics["confusion_labels"]
    assert metrics["confusion_matrix"][0][labels.index("sigma38")] == 1
    assert sum(sum(row) for row in metrics["confusion_matrix"]) == 3
    assert metrics["per_sigma"]["sigma70"]["tp"] == 1
    assert metrics["per_sigma"]["sigma70"]["fn"] == 1
    assert metrics["per_sigma"]["sigma38"]["tp"] == 0
    assert metrics["per_sigma"]["sigma38"]["fp"] == 1
    assert metrics["no_pwm_counts"] == {}
    assert metrics["macro_f1_supported"] is not None
    cal = metrics["calibration"]
    assert cal["status"] == "available" and cal["n"] == 3
    assert cal["expected_calibration_error"] is not None
    assert cal["brier_score"] is not None
    assert set(cal["by_confidence_label"]) == {"sospecha_fuerte", "sospecha_media", "sin_sospecha"}


def test_metrics_not_available_without_evaluable_rows():
    rows = [{"label_status": "no_pwm", "mapping_status": "no_pwm",
             "sigma_truth": "", "pred_sigma": ""}]
    metrics = V.compute_metrics(rows, "synthetic")
    assert metrics["status"] == "not_available"
    assert metrics["accuracy"]["top1"] is None
    assert metrics["calibration"]["status"] == "not_available"


# --------------------------------------------------------------------------
# Selftest / CLI
# --------------------------------------------------------------------------
def test_selftest():
    assert "OK" in V.selftest()


def test_cli_selftest(capsys):
    assert V.main(["--selftest"]) == 0
    assert "OK" in capsys.readouterr().out


# --------------------------------------------------------------------------
# Smoke tests con datos reales (se saltan si faltan)
# --------------------------------------------------------------------------
K12_GT = DATA / "regulondb" / "gt" / "gt_promoters.tsv"
K12_SEQ = DATA / "allostery" / "work" / "ecoli_k12.fna"
K12_EXISTING = DATA / "sigma" / "ecoli_k12" / "sigma_promoters.tsv"
BSU_GT = DATA / "dbtbs" / "gt_tss.tsv"
PWMS = DATA / "regulondb" / "gt" / "promoter_pwms.tsv"


@pytest.mark.skipif(not (K12_GT.exists() and K12_SEQ.exists() and PWMS.exists()),
                    reason="sin datos de K-12")
def test_real_k12_tss_scoring_smoke():
    pwms = S.load_pwm_library(DATA)
    rows = [row for row in V.load_k12_tss(DATA, pwms)
            if row["label_status"] == "evaluable"][:3]
    assert rows
    seqs, _source = V.load_genome_sequence(DATA, "ecoli_k12")
    assert seqs
    null = S.build_null(seqs, pwms, V.TSS_WINDOW, 20, 17)
    scored = V.score_tss_rows(rows, seqs, pwms, null)
    assert len(scored) == len(rows)
    assert all(row["contig"] for row in scored)
    assert all(row["pred_sigma"] in V.SIGMA_LABELS for row in scored)
    assert all(row["z_sigma70"] != "" for row in scored)


@pytest.mark.skipif(not (K12_GT.exists() and K12_EXISTING.exists() and PWMS.exists()),
                    reason="sin salida existente de K-12")
def test_k12_operon_comparison_smoke():
    pwms = S.load_pwm_library(DATA)
    rows = V.load_k12_tss(DATA, pwms)
    comparison = V.compare_k12_existing(DATA, rows)
    assert comparison["status"] == "evaluated"
    assert comparison["n_existing_rows"] > 0
    assert comparison["n_gt_genes_with_mappable"] > 0
    assert comparison["n_existing_matched"] > 0
    assert 0.0 <= comparison["exact_rate"] <= 1.0


@pytest.mark.skipif(not (BSU_GT.exists() and PWMS.exists()),
                    reason="sin datos de B. subtilis")
def test_bsu_mapping_smoke():
    pwms = S.load_pwm_library(DATA)
    rows, secondary = V.load_bsu_tss(DATA, pwms)
    assert rows
    assert any(row["sigma_truth"] == "sigma70" for row in rows)
    assert any(row["mapping_status"] == "no_pwm" and row["sigma_raw"] == "E"
               for row in rows)
    assert any(row["mapping_status"] == "unlabeled" for row in rows)
    assert secondary["n_rows"] >= 0
