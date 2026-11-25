"""Pruebas del minero de combinaciones de features sigma."""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "sigma_feature_combinations", ROOT / "scripts" / "sigma" / "60_feature_combinations.py")
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def test_selftest():
    assert "OK" in MOD.selftest()


def test_rule_metrics():
    rows = [
        {"sigma": "sigma54", "true": "yes", "is_top1": "yes", "z": "4",
         "pvalue": "0.001", "margin": "3", "spacer_dev": "0",
         "dist_to_start": "30", "inventory": "yes"},
        {"sigma": "sigma54", "true": "no", "is_top1": "no", "z": "0.5",
         "pvalue": "0.5", "margin": "0", "spacer_dev": "0",
         "dist_to_start": "120", "inventory": "yes"},
    ]
    res = MOD.evaluate_sigma(rows, "sigma54", min_support=1, top=5)
    assert res["baseline"]["tp"] == 1
    assert res["rules"]
    assert all(r["tp"] >= 1 for r in res["rules"])


def test_condition_helpers():
    assert MOD.CONDITIONS["z>=2"][1]({"z": "2.5"})
    assert not MOD.CONDITIONS["z>=3"][1]({"z": "2.5"})
    assert MOD.CONDITIONS["spacer_ok"][1]({"spacer_dev": "0"})
    assert not MOD.CONDITIONS["inventory"][1]({"inventory": "no"})
