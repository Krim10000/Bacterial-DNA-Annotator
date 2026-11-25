#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pruebas rápidas de scripts/sigma/74_multilabel.py.

No usan datos reales ni red: arman features, gold y promoters sintéticos en
``tmp_path`` para fijar el contrato de columnas, ranking, evaluación
multietiqueta y reglas de refinamiento.
"""
from __future__ import annotations

import csv
import importlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ml = importlib.import_module("scripts.sigma.74_multilabel")

FEATURE_COLUMNS = [
    "genome", "contig", "tss", "strand", "sigma", "true", "is_top1",
    "z", "pvalue", "score", "spacer", "spacer_dev", "dist_to_start",
    "margin", "top_z", "inventory", "n_true", "gold_confidence", "source",
]


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def _feature(contig: str, tss: int, strand: str, sigma: str, z: float,
             pvalue: float, margin: float, score: float | None = None,
             spacer: int = 17, dist: int = 30, inventory: str = "yes",
             true: str = "no", is_top1: str = "no", n_true: int = 1,
             gold_confidence: str = "S", source: str = "") -> dict[str, object]:
    return {
        "genome": "g",
        "contig": contig,
        "tss": tss,
        "strand": strand,
        "sigma": sigma,
        "true": true,
        "is_top1": is_top1,
        "z": z,
        "pvalue": pvalue,
        "score": score if score is not None else z + 5,
        "spacer": spacer,
        "spacer_dev": 0,
        "dist_to_start": dist,
        "margin": margin,
        "top_z": z,
        "inventory": inventory,
        "n_true": n_true,
        "gold_confidence": gold_confidence,
        "source": source,
    }


def _write_gold(path: Path, rows: list[tuple[str, int, str, str]]) -> None:
    _write_tsv(path, ["contig", "pos", "strand", "sigma", "gene", "confidence"],
               [{"contig": c, "pos": t, "strand": s, "sigma": sig,
                 "gene": "", "confidence": "S"}
                for c, t, s, sig in rows])


def _six_sigmas(contig: str, tss: int, strand: str,
                ranking: list[tuple[str, float]]) -> list[dict[str, object]]:
    """Devuelve seis filas; ``ranking`` fija el orden z y el resto baja."""
    known = {sigma: z for sigma, z in ranking}
    rows = []
    for idx, sigma in enumerate(ml.SIGMA_ORDER):
        z = known.get(sigma, -1.0 - idx * 0.1)
        margin = z - max([v for s, v in known.items() if s != sigma] or [z])
        rows.append(_feature(contig, tss, strand, sigma, z,
                             pvalue=min(0.9, 0.01 + idx * 0.05),
                             margin=margin, dist=30))
    return rows


def test_build_ranks_and_does_not_use_gold_columns(tmp_path: Path) -> None:
    features = tmp_path / "g" / "layer_sigma_features.tsv"
    rows = [
        _feature("chr", 100, "+", "sigma70", 5.0, 0.01, 3.0, true="no"),
        _feature("chr", 100, "+", "sigma32", 4.0, 0.02, 2.0, true="yes"),
        _feature("chr", 100, "+", "sigma54", 3.0, 0.03, 1.0, true="yes"),
        _feature("chr", 200, "-", "sigma38", 2.0, 0.04, 1.0, true="yes"),
        _feature("chr", 200, "-", "sigma24", 1.0, 0.05, 0.0, true="no"),
    ]
    _write_tsv(features, FEATURE_COLUMNS, rows)

    result = ml.build_candidates_from_layers(features, "g", top_k=2)
    assert [(r["tss"], r["rank"], r["sigma"]) for r in result.candidates] == [
        (100, 1, "sigma70"), (100, 2, "sigma32"),
        (200, 1, "sigma38"), (200, 2, "sigma24"),
    ]
    assert all("true" not in row for row in result.candidates)
    assert all("is_top1" not in row for row in result.candidates)
    assert result.stats["n_candidates_all"] == 5
    assert result.stats["n_candidates_topk"] == 4
    assert result.stats["n_input_loci"] == 2

    # Cambiar por completo las columnas del gold no debe cambiar el ranking.
    flipped = [dict(r) for r in rows]
    for row in flipped:
        row["true"] = "yes" if row["true"] == "no" else "no"
        row["is_top1"] = "yes" if row["is_top1"] == "no" else "no"
        row["n_true"] = 3
        row["gold_confidence"] = "W"
        row["source"] = "otra"
    features2 = tmp_path / "g" / "layer_sigma_features_2.tsv"
    _write_tsv(features2, FEATURE_COLUMNS, flipped)
    result2 = ml.build_candidates_from_layers(features2, "g", top_k=2)
    key = lambda r: (r["tss"], r["rank"], r["sigma"], r["z"], r["margin"])
    assert [key(r) for r in result.candidates] == [key(r) for r in result2.candidates]


def test_thresholds_filter_and_recompute_ranks(tmp_path: Path) -> None:
    features = tmp_path / "g" / "layer_sigma_features.tsv"
    _write_tsv(features, FEATURE_COLUMNS, [
        _feature("chr", 1, "+", "sigma70", 3.0, 0.01, 2.0),
        _feature("chr", 1, "+", "sigma32", 1.5, 0.02, 1.0),
        _feature("chr", 1, "+", "sigma54", 0.5, 0.30, -1.0),
    ])
    result = ml.build_candidates_from_layers(features, "g", top_k=3,
                                              min_z=1.0)
    assert [r["sigma"] for r in result.candidates] == ["sigma70", "sigma32"]
    assert [r["rank"] for r in result.candidates] == [1, 2]
    assert result.stats["n_filtered_rows"] == 1
    assert result.stats["n_candidates_topk"] == 2

    empty = ml.build_candidates_from_layers(features, "g", top_k=3, min_z=99.0)
    assert empty.candidates == []
    assert empty.stats["n_tss_without_candidates"] == 1

    by_p = ml.build_candidates_from_layers(features, "g", top_k=3,
                                            max_pvalue=0.05)
    assert [r["sigma"] for r in by_p.candidates] == ["sigma70", "sigma32"]


def _multilabel_fixture(tmp_path: Path) -> tuple[ml.BuildResult, dict, dict]:
    features = tmp_path / "g" / "layer_sigma_features.tsv"
    rows: list[dict[str, object]] = []
    rows += _six_sigmas("chr", 100, "+", [
        ("sigma70", 5.0), ("sigma32", 4.0), ("sigma54", 3.0)])
    rows += _six_sigmas("chr", 200, "+", [
        ("sigma32", 5.0), ("sigma70", 4.0), ("sigma54", 3.0)])
    rows += _six_sigmas("chr", 300, "-", [
        ("sigma38", 5.0), ("sigma24", 4.0), ("sigma70", 3.0)])
    rows += _six_sigmas("chr", 400, "-", [
        ("sigma28", 5.0), ("sigma70", 4.0), ("sigma24", 3.0)])
    _write_tsv(features, FEATURE_COLUMNS, rows)
    result = ml.build_candidates_from_layers(features, "g", top_k=3,
                                             refine=False)

    gold_path = tmp_path / "gold" / "g" / "tss_sigma.tsv"
    _write_gold(gold_path, [
        ("chr", 100, "+", "sigma70"),
        ("chr", 200, "+", "sigma32"),
        ("chr", 200, "+", "sigma54"),
        ("chr", 300, "-", "sigma38"),
        ("chr", 300, "-", "sigma24"),
        ("chr", 300, "-", "sigma28"),
        ("chr", 400, "-", "sigma24"),
    ])
    gold, gold_stats = ml.load_gold(gold_path)
    return result, gold, gold_stats


def test_multilabel_metrics_single_and_multi(tmp_path: Path) -> None:
    result, gold, gold_stats = _multilabel_fixture(tmp_path)
    evaluation = ml.evaluate_candidates(result, gold)

    assert gold_stats["n_tss"] == 4
    assert gold_stats["n_single"] == 2
    assert gold_stats["n_multi"] == 2
    ov = evaluation["overall"]
    assert (ov["top1_any"], ov["top1_all"]) == (3, 1)
    assert (ov["top3_any"], ov["top3_all"]) == (4, 3)
    assert ov["jaccard_top1"] == pytest.approx(11.0 / 24.0)
    assert ov["jaccard_top3"] == pytest.approx(11.0 / 24.0)
    assert ov["micro"]["precision"] == pytest.approx(0.5)
    assert ov["micro"]["recall"] == pytest.approx(6.0 / 7.0)
    assert ov["micro"]["f1"] == pytest.approx(2 * 0.5 * (6 / 7) / (0.5 + 6 / 7))
    assert ov["macro"]["n_labels"] == 6
    assert ov["macro"]["precision"] == pytest.approx(3.25 / 6.0)
    assert ov["macro"]["recall"] == pytest.approx(5.0 / 6.0)
    assert ov["macro"]["f1"] == pytest.approx((0.4 + 2 * (2 / 3) + 1 + 1 + 0) / 6.0)

    assert evaluation["single"]["n_tss"] == 2
    assert evaluation["single"]["top3_all"] == 2
    assert evaluation["multi"]["n_tss"] == 2
    assert evaluation["multi"]["top1_all"] == 0
    assert evaluation["multi"]["top3_any"] == 2
    assert evaluation["multi"]["top3_all"] == 1

    by_sigma = {r["sigma"]: r for r in evaluation["per_sigma"]}
    assert by_sigma["sigma70"]["support"] == 1
    assert by_sigma["sigma70"]["predicted"] == 4
    assert by_sigma["sigma70"]["tp"] == 1
    assert by_sigma["sigma38"]["f1"] == pytest.approx(1.0)
    assert by_sigma["sigma28"]["recall"] == 0.0
    assert by_sigma["sigma28"]["f1"] == 0.0

    coverage = evaluation["multi_sigma_coverage"]
    assert coverage["n_multi_tss"] == 2
    assert coverage["multi_any"] == 2
    assert coverage["multi_all"] == 1
    assert coverage["multi_not_fully_covered"] == 1
    assert coverage["multi_missed_entirely"] == 0
    assert coverage["gold_labels_multi"] == 5
    assert coverage["gold_labels_covered"] == 4
    assert coverage["gold_labels_outside"] == 1

    k3 = {c["k"]: c for c in evaluation["topk_sensitivity"]}[3]
    assert k3["multi_all"] == 1
    assert k3["gold_labels_outside"] == 1
    assert evaluation["multi_examples"][0]["missing"] == ["sigma28"]


def test_refinement_changes_confidence_not_sigma(tmp_path: Path) -> None:
    features = tmp_path / "g" / "layer_sigma_features.tsv"
    rows = [
        _feature("chr", 1, "+", "sigma54", 4.0, 0.01, 2.0, dist=30),
        _feature("chr", 2, "+", "sigma54", 4.0, 0.01, 2.0, dist=90),
        _feature("chr", 3, "+", "sigma38", 3.0, 0.20, 1.5, dist=40),
        _feature("chr", 4, "+", "sigma38", 3.0, 0.20, 0.5, dist=40),
        _feature("chr", 1, "+", "sigma70", 2.0, 0.02, 0.0),
        _feature("chr", 2, "+", "sigma70", 2.0, 0.02, 0.0),
        _feature("chr", 3, "+", "sigma70", 2.0, 0.02, 0.0),
        _feature("chr", 4, "+", "sigma70", 2.0, 0.02, 0.0),
    ]
    _write_tsv(features, FEATURE_COLUMNS, rows)
    result = ml.build_candidates_from_layers(features, "g", top_k=3,
                                             refine=True)
    by_key = {(r["tss"], r["sigma"]): r for r in result.candidates}
    assert by_key[(1, "sigma54")]["confidence_label"] == "alta"
    assert by_key[(1, "sigma54")]["refined_confidence"] == "alta"
    assert by_key[(2, "sigma54")]["confidence_label"] == "alta"
    assert by_key[(2, "sigma54")]["refined_confidence"] == "media"
    assert by_key[(3, "sigma38")]["confidence_label"] == "media"
    assert by_key[(3, "sigma38")]["refined_confidence"] == "alta"
    assert by_key[(4, "sigma38")]["refined_confidence"] == "baja"
    assert by_key[(1, "sigma70")]["refined_confidence"] == \
        by_key[(1, "sigma70")]["confidence_label"]
    assert {r["sigma"] for r in result.candidates} == {"sigma54", "sigma38", "sigma70"}

    gold_path = tmp_path / "gold" / "g" / "tss_sigma.tsv"
    _write_gold(gold_path, [
        ("chr", 1, "+", "sigma54"), ("chr", 2, "+", "sigma54"),
        ("chr", 3, "+", "sigma38"), ("chr", 4, "+", "sigma38"),
    ])
    gold, _ = ml.load_gold(gold_path)
    refinement = ml.evaluate_refinement(result, gold)
    by_sigma = {r["sigma"]: r for r in refinement["per_sigma"]}
    assert by_sigma["sigma54"]["n_pass"] == 1
    assert by_sigma["sigma54"]["tp"] == 1
    assert by_sigma["sigma54"]["fn"] == 1
    assert by_sigma["sigma38"]["n_pass"] == 1
    assert by_sigma["sigma38"]["tp"] == 1
    assert by_sigma["sigma38"]["fn"] == 1

    no_refine = ml.build_candidates_from_layers(features, "g", top_k=3,
                                                refine=False)
    assert all(r["refined_confidence"] == r["confidence_label"]
               for r in no_refine.candidates)
    assert ml.evaluate_refinement(no_refine, gold)["enabled"] is False


def test_promoter_fallback_uses_operon_anchor(tmp_path: Path) -> None:
    promoters = tmp_path / "sigma" / "g" / "sigma_promoters.tsv"
    alternatives = [
        {"sigma": "sigma54", "z": 3.0, "pvalue": 0.01, "score": 9.0, "spacer": 11},
        {"sigma": "sigma70", "z": 2.0, "pvalue": 0.02, "score": 8.0, "spacer": 17},
        {"sigma": "sigma38", "z": 1.0, "pvalue": 0.10, "score": 7.0, "spacer": 17},
    ]
    _write_tsv(promoters,
               ["operon_id", "contig", "strand", "sigma", "zscore", "pvalue",
                "score", "score_margin", "dist_to_start", "evidence",
                "alternatives"],
               [{
                   "operon_id": "g:chr1:500:+",
                   "contig": "chr1",
                   "strand": "1",
                   "sigma": "sigma54",
                   "zscore": 3.0,
                   "pvalue": 0.01,
                   "score": 9.0,
                   "score_margin": 1.0,
                   "dist_to_start": 20,
                   "evidence": "evidencia del operón",
                   "alternatives": json.dumps(alternatives),
               }])
    result = ml.build_candidates_from_promoters(promoters, "g", top_k=3,
                                                refine=True)
    assert result.source == "promoters"
    assert [r["sigma"] for r in result.candidates] == [
        "sigma54", "sigma70", "sigma38"]
    assert all(r["tss"] == "" for r in result.candidates)
    assert all(r["coord_type"] == "operon_anchor" for r in result.candidates)
    assert all(r["operon_anchor"] == 500 for r in result.candidates)
    assert all("NO TSS" in r["evidence"] for r in result.candidates)
    assert result.candidates[0]["spacer_dev"] == 0
    assert result.candidates[0]["margin"] == pytest.approx(1.0)

    source, path = ml.resolve_source("g", "auto", tmp_path / "layers",
                                     tmp_path / "sigma")
    assert source == "promoters"
    assert path == promoters

    gold_path = tmp_path / "gold" / "g" / "tss_sigma.tsv"
    _write_gold(gold_path, [("chr1", 500, "+", "sigma54")])
    gold, gold_stats = ml.load_gold(gold_path)
    evaluation = ml.evaluate_candidates(result, gold)
    refinement = ml.evaluate_refinement(result, gold)
    assert evaluation["status"] == "not_applicable"
    assert refinement["status"] == "not_applicable"
    report = ml.build_report("g", result, gold_stats, evaluation, refinement)
    assert "No se calculan métricas de TSS" in report


def test_cli_end_to_end(tmp_path: Path) -> None:
    features = tmp_path / "layers" / "g" / "layer_sigma_features.tsv"
    _write_tsv(features, FEATURE_COLUMNS, [
        _feature("chr", 10, "+", "sigma70", 5.0, 0.01, 2.0),
        _feature("chr", 10, "+", "sigma32", 4.0, 0.02, 1.0),
        _feature("chr", 10, "+", "sigma54", 3.0, 0.03, 0.5),
        _feature("chr", 10, "+", "sigma38", 2.0, 0.04, 0.0),
        _feature("chr", 10, "+", "sigma28", 1.0, 0.05, -0.5),
        _feature("chr", 10, "+", "sigma24", 0.5, 0.06, -1.0),
        _feature("chr", 20, "-", "sigma32", 5.0, 0.01, 2.0),
        _feature("chr", 20, "-", "sigma70", 4.0, 0.02, 1.0),
        _feature("chr", 20, "-", "sigma24", 3.0, 0.03, 0.5),
        _feature("chr", 20, "-", "sigma54", 2.0, 0.04, 0.0),
        _feature("chr", 20, "-", "sigma38", 1.0, 0.05, -0.5),
        _feature("chr", 20, "-", "sigma28", 0.5, 0.06, -1.0),
    ])
    gold = tmp_path / "gold" / "g" / "tss_sigma.tsv"
    _write_gold(gold, [
        ("chr", 10, "+", "sigma70"),
        ("chr", 20, "-", "sigma32"),
        ("chr", 20, "-", "sigma24"),
    ])
    out = tmp_path / "out"

    code = ml.main([
        "--genome", "g",
        "--features-root", str(tmp_path / "layers"),
        "--gold-root", str(tmp_path / "gold"),
        "--out-root", str(out),
        "--top-k", "3",
    ])
    assert code == 0
    candidates = out / "g" / "sigma_candidates.tsv"
    metrics = out / "g" / "metrics.json"
    report = out / "g" / "report.md"
    assert candidates.exists() and metrics.exists() and report.exists()

    with candidates.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        header = list(reader.fieldnames or [])
    assert header[:10] == [
        "genome", "contig", "tss", "strand", "locus_id", "coord_type",
        "operon_anchor", "rank", "sigma", "z"]
    assert len(rows) == 6
    assert {r["coord_type"] for r in rows} == {"tss"}
    assert all(r["tss"] for r in rows)

    payload = json.loads(metrics.read_text(encoding="utf-8"))
    assert payload["genome"] == "g"
    assert payload["evaluation"]["status"] == "ok"
    assert payload["evaluation"]["overall"]["top3_any"] == 2
    assert payload["config"]["refine"] is True
    assert payload["refinement"]["enabled"] is True
    assert "refinement" in payload
    text = report.read_text(encoding="utf-8")
    assert "top3_any" in text
    assert "## 9. Límites" in text
    assert "sigma_candidates.tsv" in text


def test_top_k_under_three_keeps_top3_undefined(tmp_path: Path) -> None:
    result, gold, _ = _multilabel_fixture(tmp_path)
    result.top_k = 2
    evaluation = ml.evaluate_candidates(result, gold)
    assert evaluation["overall"]["top3_any"] is None
    assert evaluation["overall"]["top3_all"] is None
    assert evaluation["overall"]["topk_any"] is not None
