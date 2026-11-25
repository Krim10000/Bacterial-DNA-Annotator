"""Pruebas de la calibración de la capa 2 de motivo para sigma54 y sigma38."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "sigma_motif_calibration",
    ROOT / "scripts" / "sigma" / "72_motif_calibration.py")
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def _prob_pwm(probs: list[list[float]]) -> np.ndarray:
    return np.log2(np.asarray(probs, dtype=float) / 0.25)


def _rec() -> dict:
    up = _prob_pwm([
        [0.97, 0.01, 0.01, 0.01],
        [0.01, 0.97, 0.01, 0.01],
        [0.40, 0.30, 0.20, 0.10],
        [0.10, 0.20, 0.30, 0.40],
    ])
    down = _prob_pwm([
        [0.01, 0.01, 0.01, 0.97],
        [0.40, 0.30, 0.20, 0.10],
        [0.10, 0.20, 0.30, 0.40],
    ])
    return {
        "sigma": "sigma54",
        "class": "alternative",
        "arch": "-24/-12",
        "up": {"width": 4, "pwm": up},
        "down": {"width": 3, "pwm": down},
    }


def _rows(n_true: int = 5, n_false: int = 5) -> list[dict[str, str]]:
    rows = []
    for i in range(n_true):
        rows.append({"true": "yes", "is_top1": "yes", "z": str(4 - 0.5 * i),
                     "pvalue": "0.001", "margin": "2",
                     "dist_to_start": "26"})
    for i in range(n_false):
        rows.append({"true": "no", "is_top1": "no", "z": str(0 - 0.5 * i),
                     "pvalue": "0.5", "margin": "-1",
                     "dist_to_start": "60"})
    return rows


def test_selftest():
    assert "OK" in MOD.selftest()


def test_position_information_identifies_conserved_columns():
    rec = _rec()
    info_up = MOD.position_information(rec["up"]["pwm"])
    info_down = MOD.position_information(rec["down"]["pwm"])
    assert info_up.shape == (4,)
    assert info_down.shape == (3,)
    assert info_up[0] > info_up[2]
    assert info_up[1] > info_up[2]
    assert info_down[0] > info_down[1]


def test_build_core_preserves_widths_and_masks_columns():
    rec = _rec()
    original_up = rec["up"]["pwm"].copy()
    core, kept = MOD.build_core_pwm(rec, 2, "box")
    assert core["up"]["width"] == 4
    assert core["down"]["width"] == 3
    assert len(kept["up"]) == 2
    assert len(kept["down"]) == 2
    # Las columnas descartadas quedan todas a cero; las conservadas no cambian.
    assert int(np.sum(~np.any(core["up"]["pwm"], axis=1))) == 2
    assert int(np.sum(~np.any(core["down"]["pwm"], axis=1))) == 1
    for pos in range(4):
        if pos + 1 in kept["up"]:
            assert np.array_equal(core["up"]["pwm"][pos], original_up[pos])
        else:
            assert not np.any(core["up"]["pwm"][pos])
    # El modelo original no se muta.
    assert np.array_equal(rec["up"]["pwm"], original_up)


def test_build_core_joint_keeps_exactly_k_positions():
    rec = _rec()
    core, kept = MOD.build_core_pwm(rec, 3, "joint")
    assert len(kept["up"]) + len(kept["down"]) == 3
    zero_up = int(np.sum(~np.any(core["up"]["pwm"], axis=1)))
    zero_down = int(np.sum(~np.any(core["down"]["pwm"], axis=1)))
    assert zero_up == 4 - len(kept["up"])
    assert zero_down == 3 - len(kept["down"])


def test_auc_score_perfect_reversed_and_ties():
    assert MOD.auc_score([True, False, True, False], [3, 1, 2, 0]) == 1.0
    assert MOD.auc_score([False, True], [3, 1]) == 0.0
    assert MOD.auc_score([True, False], [1, 1]) == 0.5
    assert MOD.auc_score([True, True], [1, 2]) is None


def test_binary_metrics_counts():
    metrics = MOD.binary_metrics_from_mask(
        [True, True, False, False], [True, False, True, False])
    assert metrics["tp"] == 1
    assert metrics["fp"] == 1
    assert metrics["fn"] == 1
    assert metrics["n_pred"] == 2
    assert metrics["precision"] == pytest.approx(0.5)
    assert metrics["recall"] == pytest.approx(0.5)


def test_calibrate_single_prefers_strict_threshold_on_ties():
    arrays = MOD.prepare_arrays(_rows())
    chosen = MOD.calibrate_single(arrays, "z", "f1")
    assert chosen is not None
    # z >= 2.0 y z >= 2.5 rinden igual en K-12; se prefiere el más exigente.
    assert chosen["threshold"] == pytest.approx(2.0)
    assert chosen["metrics"]["f1"] == pytest.approx(1.0)


def test_calibration_rule_transfers_without_retuning():
    calibration = MOD.prepare_arrays(_rows())
    chosen = MOD.calibrate_single(calibration, "z", "f1")
    transfer_rows = [
        {"true": "yes", "is_top1": "no", "z": "3.0", "pvalue": "0.01",
         "margin": "1", "dist_to_start": "25"},
        {"true": "yes", "is_top1": "no", "z": "2.2", "pvalue": "0.02",
         "margin": "1", "dist_to_start": "28"},
        {"true": "no", "is_top1": "yes", "z": "1.0", "pvalue": "0.2",
         "margin": "0", "dist_to_start": "50"},
        {"true": "no", "is_top1": "no", "z": "-1.0", "pvalue": "0.8",
         "margin": "-2", "dist_to_start": "70"},
    ]
    metrics = MOD.evaluate_rule(chosen["rule"], MOD.prepare_arrays(transfer_rows))
    assert metrics["tp"] == 2
    assert metrics["fp"] == 0
    assert metrics["recall"] == pytest.approx(1.0)


def test_two_feature_rule_has_two_features_and_ignores_gold_metadata():
    rows = []
    for i in range(20):
        true = i < 8
        rows.append({
            "true": "yes" if true else "no",
            "is_top1": "yes" if true else "no",
            "z": str(3.0 if true else 0.0),
            "pvalue": "0.01" if true else "0.5",
            "margin": str(1.5 if true else -0.5),
            "dist_to_start": str(28 if true else 60),
        })
    base = MOD.prepare_arrays(rows)
    chosen = MOD.select_two_feature(MOD.two_feature_candidates(base), "f1")
    assert chosen is not None
    assert len(chosen["rule"]["conditions"]) == 2
    features = {c["feature"] for c in chosen["rule"]["conditions"]}
    assert features <= {"z", "p", "margin", "dist"}
    assert not features & set(MOD.FORBIDDEN_COLUMNS)
    # Cambiar metadatos del gold no puede cambiar la regla elegida.
    for row in rows:
        row["n_true"] = "3"
        row["gold_confidence"] = "S"
        row["source"] = "RegulonDB"
    chosen2 = MOD.select_two_feature(
        MOD.two_feature_candidates(MOD.prepare_arrays(rows)), "f1")
    assert chosen2 is not None
    assert chosen2["rule"]["label"] == chosen["rule"]["label"]


def test_distance_and_margin_are_candidate_features():
    assert MOD.FEATURE_OPS["dist"] == "<="
    assert MOD.FEATURE_OPS["margin"] == ">="
    assert "spacer" not in MOD.FEATURE_GRID
    assert "inventory" not in MOD.FEATURE_GRID


def test_limit_rows_keeps_all_sigmas_for_selected_tss():
    rows = [
        {"contig": "c", "tss": "10", "strand": "+", "sigma": s}
        for s in ("sigma54", "sigma38")
    ] + [
        {"contig": "c", "tss": "20", "strand": "+", "sigma": s}
        for s in ("sigma54", "sigma38")
    ]
    limited = MOD._limit_rows(rows, 1)
    assert len(limited) == 2
    assert {r["tss"] for r in limited} == {"10"}


def test_generated_outputs_match_contract_if_present():
    out = ROOT / "data" / "sigma" / "motif_calibration"
    if not (out / "summary.json").exists():
        pytest.skip("todavía no se generaron las salidas")
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["version"] == MOD.VERSION
    assert set(summary["baselines_top1"]) == {"ecoli_k12", "subtilis"}
    core = summary["core"]
    assert "ecoli_k12/sigma54" in core
    assert "subtilis/sigma38" in core
    rows = summary.get("calibration", {})
    assert set(rows) == {"sigma54", "sigma38"}
    two = summary.get("two_feature", {})
    assert set(two) == {"sigma54", "sigma38"}
    for info in two.values():
        assert len({c for c in ("z", "p", "margin", "dist")
                    if c in info["rule"].split()}) <= 2
