"""Pruebas del comparador de ventanas/fondos de la capa 1 sigma.

No corren el barrido completo: validan las piezas que pueden romper la
comparación (ventanas, shuffles, flancos, estadística y scorer vectorizado).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SPEC = importlib.util.spec_from_file_location(
    "sigma_local_background",
    ROOT / "scripts" / "sigma" / "70_local_background.py")
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)

PWMS = DATA / "regulondb" / "gt" / "promoter_pwms.tsv"


def test_selftest():
    assert "OK" in MOD.selftest()


def test_tss_window_matches_current_evaluator():
    seq = "ACGTTGCA" * 200
    for tss in (1, 2, 40, 300, 1599, 1600):
        for strand in ("+", "-"):
            for up, down in ((60, 20), (80, 20), (100, 20), (120, 20)):
                got = MOD.tss_window(seq, tss, strand, up, down)
                expected = MOD._EVAL._window(seq, tss, strand, up, down)
                assert got == expected, (tss, strand, up, down)


def test_flanks_are_local_and_do_not_include_the_core_window():
    seq = "ACGT" * 1000
    tss = 1000
    up, down = MOD.flank_segments(seq, tss, "+")
    assert up == seq[499:899]
    assert down == seq[1099:1499]
    core = MOD.tss_window(seq, tss, "+", 60, 20)
    assert core == seq[939:1020]
    # Intervalos disjuntos: el flanco -500/-100 no toca la ventana -60/+20.
    assert 899 <= 939 and 1020 <= 1099
    # En hebra '-' los flancos se revcomp y siguen siendo locales.
    up_r, down_r = MOD.flank_segments(seq, tss, "-")
    assert up_r == MOD.revcomp(seq[1100:1500])
    assert down_r == MOD.revcomp(seq[500:900])


def test_local_shuffles_preserve_composition_and_are_reproducible():
    window = "AACCGGTTAACCGTTAA"
    a = MOD.local_shuffles(window, 25, 123)
    b = MOD.local_shuffles(window, 25, 123)
    assert a == b
    assert len(a) == 25
    assert all(len(x) == len(window) for x in a)
    assert all(sorted(x) == sorted(window) for x in a)
    # Semillas distintas deben producir órdenes distintos con alta probabilidad.
    c = MOD.local_shuffles(window, 25, 124)
    assert a != c


@pytest.mark.skipif(not PWMS.exists(), reason="sin PWMs en esta instalación")
def test_batch_scorer_matches_annotator_exactly():
    pwms = MOD.load_pwm_library(DATA)
    base = ("ACGTTGCAAGCTTAGCGATCG" * 10)
    windows = [base[i:i + 80] for i in range(0, 40, 10)]
    windows.append(MOD.revcomp(base)[:80])
    batched = MOD.score_batch(windows, pwms)
    for sigma, rec in pwms.items():
        expected = np.array([MOD.score_sigma(w, rec, len(w))["score"]
                             for w in windows])
        assert np.allclose(batched[sigma], expected, atol=1e-9), sigma


def test_null_stats_match_evaluator_definitions():
    null = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    stats = MOD.null_stats(6.0, null)
    assert stats["n"] == 5
    assert stats["pvalue"] == pytest.approx(1 / 6)
    assert stats["pct"] == pytest.approx(1.0)
    assert stats["z"] == pytest.approx(
        (6.0 - null.mean()) / null.std(ddof=1))
    # Empate: `pct` cuenta estrictamente menores, `pvalue` incluye el empate.
    stats = MOD.null_stats(3.0, null)
    assert stats["pct"] == pytest.approx(2 / 5)
    assert stats["pvalue"] == pytest.approx(4 / 6)


def test_summarize_stats_and_support_rule():
    rows = [
        {"pvalue": 0.01, "pct": 0.99, "z": 2.5},
        {"pvalue": 0.04, "pct": 0.96, "z": 1.2},
        {"pvalue": 0.20, "pct": 0.90, "z": 0.5},
    ]
    summary = MOD._summarize_stats(rows, 0.05, 0.95, 1.0)
    assert summary["n"] == 3
    assert summary["frac_p05"] == pytest.approx(2 / 3)
    assert summary["frac_pct95"] == pytest.approx(2 / 3)
    assert summary["frac_supported"] == pytest.approx(2 / 3)
    assert summary["median_z"] == pytest.approx(1.2)


def test_recommendation_uses_guardrail_and_tie_break():
    def result(genome, up, bg, support, p05=0.5, top1=0.5, top3=0.88):
        return {
            "genome": genome,
            "window": [-up, 20],
            "window_up": up,
            "window_down": 20,
            "background": bg,
            "global_summary": {"frac_supported": support, "frac_p05": p05},
            "assignment": {"top1_accuracy": top1, "top3_any": top3},
        }

    results = [
        # Referencia -80/genomic.
        result("ecoli_k12", 80, "genomic", 0.45, top1=0.50, top3=0.88),
        result("subtilis", 80, "genomic", 0.40, top1=0.48, top3=0.86),
        # Alto soporte pero asignación degradada: debe quedar fuera.
        result("ecoli_k12", 60, "flank", 0.70, top1=0.40, top3=0.70),
        result("subtilis", 60, "flank", 0.60, top1=0.38, top3=0.68),
        # Mejor soporte sin degradar asignación.
        result("ecoli_k12", 60, "genomic", 0.55, p05=0.53, top1=0.52, top3=0.89),
        result("subtilis", 60, "genomic", 0.45, p05=0.47, top1=0.49, top3=0.87),
        # Empata en soporte con -60/genomic; gana por p05 y ventana más corta.
        result("ecoli_k12", 70, "genomic", 0.55, p05=0.40, top1=0.52, top3=0.89),
        result("subtilis", 70, "genomic", 0.45, p05=0.40, top1=0.49, top3=0.87),
    ]
    rec = MOD._recommend(results, ["genomic", "flank"])
    assert rec["chosen"]["window"] == [-60, 20]
    assert rec["chosen"]["background"] == "genomic"
    # La combinación de mayor soporte no pasó el guardarraíl.
    assert rec["n_accepted"] == 3
    assert rec["chosen"]["allowed_by_guardrail"] is True
