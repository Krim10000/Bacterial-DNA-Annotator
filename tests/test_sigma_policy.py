import pytest
from webapp.sigma_policy import argmax_order, two_stage_order, hybrid_order


def test_argmax():
    scores = {"sigma70": 2.0, "sigma38": 3.0}
    assert argmax_order(scores)[0] == "sigma38"


def test_two_stage_keeps_sigma70_when_margin_small():
    scores = {"sigma70": 2.0, "sigma38": 2.2}
    assert two_stage_order(scores, delta=1.0)[0] == "sigma70"


def test_two_stage_switches_when_margin_large():
    scores = {"sigma70": 2.0, "sigma38": 3.5}
    assert two_stage_order(scores, delta=1.0)[0] == "sigma38"


def test_hybrid_falls_back_when_external_not_conclusive():
    scores = {"sigma70": 2.0, "sigma38": 2.2}
    external = {"sigma38": 0.4, "sigma70": 0.6}
    assert hybrid_order(scores, external, threshold=0.75)[0] == "sigma70"


def test_hybrid_uses_external_when_conclusive():
    scores = {"sigma70": 2.0, "sigma38": 2.2}
    external = {"sigma38": 0.9, "sigma70": 0.1}
    assert hybrid_order(scores, external, threshold=0.75)[0] == "sigma38"


def test_select_best_two_stage_uses_sigma70_fallback():
    from scripts.sigma.sigma_annotator import select_best
    ranked = [{"sigma": "sigma38", "z": 2.2}, {"sigma": "sigma70", "z": 2.0}]
    best, margin = select_best(ranked, policy="two_stage", delta=1.0)
    assert best["sigma"] == "sigma70"
    assert margin == pytest.approx(-0.2)


def test_select_best_rejects_unknown_policy():
    from scripts.sigma.sigma_annotator import select_best
    with pytest.raises(ValueError):
        select_best([{"sigma": "sigma70", "z": 1}], policy="magic")
