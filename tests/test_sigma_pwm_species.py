import json
from pathlib import Path

from scripts.sigma.sigma_annotator import load_pwm_library, score_sigma


def test_one_box_pwm_loads_and_scores(tmp_path):
    pwm = tmp_path / "pwm.tsv"
    pwm.write_text(
        "sigma\tbox\tn_sites\twidth\tconsensus\tinfo_bits\tpwm\n"
        "sigma70\t-10\t10\t4\tTATA\t1.0\t"
        + json.dumps([[1, -1, -1, -1], [-1, 1, -1, -1],
                      [-1, -1, 1, -1], [-1, -1, -1, 1]]) + "\n",
        encoding="utf-8",
    )
    lib = load_pwm_library(tmp_path, pwm)
    assert "sigma70" in lib and "up" not in lib["sigma70"]
    hit = score_sigma("N" * 5 + "TATA" + "N" * 5, lib["sigma70"], 14)
    assert hit is not None
    assert hit["arch"] == "single-box"
    assert hit["down_box"] == "TATA"


def test_bsu_pwm_has_species_matrices():
    import pytest
    path = Path("data/dbtbs/promoter_pwms_subtilis_all.tsv")
    if not path.exists():
        pytest.skip("sin PWM de B. subtilis")
    lib = load_pwm_library(Path("data"), path)
    assert {"sigma70", "sigma38", "sigma28", "sigma32", "sigma54"} <= set(lib)
    assert all("up" not in rec for rec in lib.values())
