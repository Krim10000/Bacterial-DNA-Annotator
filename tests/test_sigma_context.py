"""Pruebas de la capa de contexto sigma.

Cubren:
  * índice genómico y de operones (hebras + y -);
  * bins de contexto y la regla anti-sobreajuste de `fam`;
  * log-odds con `none` neutral, `min_support` y pesos exclusivos (ablación);
  * combinación, guardarraíl de top-1 y geometría sigma54;
  * calibración (ECE/Brier/NLL) y bootstrap;
  * el flag `--pwm` de `50_eval_layers.py` (bug: se parseaba y se ignoraba);
  * una corrida sintética de extremo a extremo de `51_eval_context.py`.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

from scripts.sigma import sigma_context as C

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------
# Índices
# --------------------------------------------------------------------------
def test_selftest():
    assert "OK" in C.selftest()


def test_gene_index_plus_minus():
    idx = C.GeneIndex([
        {"contig": "c", "start": 100, "end": 200, "strand": "+", "gene_name": "a"},
        {"contig": "c", "start": 300, "end": 400, "strand": "+", "gene_name": "b"},
        {"contig": "c", "start": 500, "end": 600, "strand": "-", "gene_name": "c"},
    ])
    assert idx.downstream("c", 50, "+")["gene_name"] == "a"
    assert idx.downstream("c", 250, "+")["gene_name"] == "b"
    assert idx.downstream("c", 650, "-")["gene_name"] == "c"
    assert idx.inside("c", 150, "+")["gene_name"] == "a"
    assert idx.inside("c", 250, "+") is None
    assert idx.inside("c", 550, "-")["gene_name"] == "c"


def test_operon_index_containing_and_downstream():
    idx = C.OperonIndex([
        {"contig": "c", "start": 90, "end": 200, "strand": "+", "n_genes": "2"},
        {"contig": "c", "start": 290, "end": 400, "strand": "+", "n_genes": "1"},
    ])
    assert idx.containing_or_downstream("c", 50, "+")["n_genes"] == "2"
    assert idx.containing_or_downstream("c", 150, "+")["n_genes"] == "2"
    assert idx.containing_or_downstream("c", 250, "+")["n_genes"] == "1"


def test_context_bins():
    assert C._bin_leader(None) == "none"
    assert C._bin_leader(10) == "0-50"
    assert C._bin_leader(100) == "50-150"
    assert C._bin_leader(300) == "150-500"
    assert C._bin_leader(900) == ">500"
    assert C._bin_operon_size(1) == "1"
    assert C._bin_operon_size(4) == "3-5"
    assert C._bin_operon_size(9) == ">5"
    assert C._bin_conservation(3) == "low"
    assert C._bin_conservation(30) == "mid"
    assert C._bin_conservation(300) == "high"
    assert C._ec_class("2.7.1.1") == "2"
    assert C._ec_class("") == "none"


# --------------------------------------------------------------------------
# Construcción de contexto
# --------------------------------------------------------------------------
def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_build_context_genomic(tmp_path: Path):
    data = tmp_path / "data"
    _write_tsv(data / "regulation" / "toy" / "genes.tsv",
               ["contig", "start", "end", "strand", "gene_id", "gene_name", "product", "ec"],
               [["c", "100", "200", "+", "g1", "a", "enzyme", "2.7.1.1"],
                ["c", "300", "400", "+", "g2", "b", "other", ""]])
    _write_tsv(data / "regulation" / "toy" / "operons.tsv",
               ["operon_id", "contig", "start", "end", "strand", "n_genes"],
               [["op1", "c", "90", "200", "+", "2"]])
    sites = [("c", 50, "+"), ("c", 150, "+"), ("c", 250, "+")]
    ctx = C.build_context(data, "toy", sites)
    assert ctx.available["operon_position"] is True
    assert ctx.available["pathway_ec"] is True
    assert ctx.available["operon_size"] is True
    row = ctx.row(("c", 50, "+"))
    assert row["operon_position"] == "gene_edge"
    assert row["leader"] == "0-50"
    assert row["pathway_ec"] == "2"
    assert row["operon_size"] == "2"
    assert ctx.row(("c", 150, "+"))["operon_position"] == "inside_gene"
    assert ctx.coverage["operon_position"] == 1.0


def test_fam_is_not_used_as_domain(tmp_path: Path):
    data = tmp_path / "data"
    _write_tsv(data / "regulation" / "toy" / "genes.tsv",
               ["contig", "start", "end", "strand", "gene_id", "gene_name", "fam"],
               [["c", "100", "200", "+", "g1", "a", "INH10"]])
    ctx = C.build_context(data, "toy", [("c", 50, "+")])
    assert ctx.available["domain"] is False
    assert ctx.row(("c", 50, "+"))["domain"] == "none"


def test_active_weights_respect_coverage():
    ctx = C.ContextTable(genome="toy")
    ctx.available = {name: True for name in C.EVIDENCE_SPECS}
    ctx.coverage = {name: 0.0 for name in C.EVIDENCE_SPECS}
    ctx.coverage["operon_position"] = 1.0
    weights = C.active_weights(ctx, min_coverage=0.05)
    assert weights["operon_position"] == 1.0
    assert weights["pathway_ec"] == 0.0
    assert weights["expression"] == 0.0


# --------------------------------------------------------------------------
# Log-odds
# --------------------------------------------------------------------------
def _toy_table() -> tuple[C.ContextTable, dict]:
    ctx = C.ContextTable(genome="toy")
    rows = [
        {"operon_position": "gene_edge", "pathway_ec": "1"},
        {"operon_position": "gene_edge", "pathway_ec": "1"},
        {"operon_position": "inside_gene", "pathway_ec": "2"},
        {"operon_position": "inside_gene", "pathway_ec": "2"},
    ]
    keys = [("c", i, "+") for i in range(4)]
    for key, row in zip(keys, rows):
        ctx.sites[key] = {**row, "leader": "0-50", "domain": "none",
                          "conservation": "none", "expression": "none",
                          "operon_size": "1", "sigma54_geometry": "none"}
    labels = {keys[0]: {"sigma70"}, keys[1]: {"sigma70"},
              keys[2]: {"sigma32"}, keys[3]: {"sigma32"}}
    return ctx, labels


def test_logodds_sign_and_missing_neutral():
    ctx, labels = _toy_table()
    model = C.fit_logodds(ctx, labels, ["operon_position", "pathway_ec"],
                          min_support=1)
    edge = ctx.row(("c", 0, "+"))
    inside = ctx.row(("c", 2, "+"))
    assert model.llr(edge, "sigma70") > 0
    assert model.llr(inside, "sigma70") < 0
    assert model.llr({"operon_position": "none", "pathway_ec": "none"},
                     "sigma70") == 0.0


def test_llr_weights_are_exclusive():
    """Un dict explícito no debe volver a sumar las demás evidencias."""
    ctx, labels = _toy_table()
    model = C.fit_logodds(ctx, labels, ["operon_position", "pathway_ec"],
                          min_support=1)
    row = ctx.row(("c", 0, "+"))
    only_pos = model.llr(row, "sigma70", {"operon_position": 1.0})
    only_ec = model.llr(row, "sigma70", {"pathway_ec": 1.0})
    both = model.llr(row, "sigma70", {"operon_position": 1.0, "pathway_ec": 1.0})
    assert only_pos != 0.0
    assert only_ec != 0.0
    assert abs(both - (only_pos + only_ec)) < 1e-9


def test_min_support_makes_rare_category_neutral():
    ctx, labels = _toy_table()
    model = C.fit_logodds(ctx, labels, ["operon_position"], min_support=5)
    assert model.llr(ctx.row(("c", 0, "+")), "sigma70") == 0.0


# --------------------------------------------------------------------------
# Combinación y guardarraíl
# --------------------------------------------------------------------------
def test_combine_scale_zero_is_motif():
    z = {"sigma70": 3.0, "sigma32": 2.0}
    model = C.LogOddsModel()
    ctx_row = {"operon_position": "gene_edge"}
    assert C.combine(z, model, ctx_row, scale=0.0) == z


def test_select_top1_guardrail():
    z = {"sigma70": 3.0, "sigma32": 2.0}
    final = {"sigma70": 3.0, "sigma32": 5.0}
    assert C.select_top1(final, z, gap=0.0) == "sigma70"
    assert C.select_top1(final, z, gap=None) == "sigma32"
    assert C.select_top1(final, z, gap=1.5) == "sigma32"


def test_sigma54_geometry():
    assert C.sigma54_geometry_category({"arch": "-24/-12", "spacer": 11}, "x") == "canonical"
    assert C.sigma54_geometry_category({"arch": "-24/-12", "spacer": 17}, "x") == "shifted"
    assert C.sigma54_geometry_category({"arch": "single-box", "dist_to_start": 12}, "x") == "canonical"
    assert C.sigma54_geometry_category({"arch": "single-box", "dist_to_start": 60}, "x") == "shifted"
    assert C.sigma54_geometry_category(None, "x") == "none"


# --------------------------------------------------------------------------
# Calibración
# --------------------------------------------------------------------------
def test_calibration_metrics_perfect():
    records = [
        {"probs": {"sigma70": 1.0, "sigma32": 0.0}, "true": {"sigma70"}},
        {"probs": {"sigma70": 0.0, "sigma32": 1.0}, "true": {"sigma32"}},
    ]
    m = C.calibration_metrics(records, sigmas=("sigma70", "sigma32"))
    assert m["ece"] is not None and m["ece"] < 1e-6
    assert m["brier"] == 0.0
    assert m["top1_accuracy"] == 1.0


def test_bootstrap_delta_identical():
    a = [1, 0, 1, 1, 0]
    out = C.bootstrap_delta(a, a, n_boot=100, seed=17)
    assert out["delta"] == 0.0
    assert out["ci_low"] == 0.0 and out["ci_high"] == 0.0


# --------------------------------------------------------------------------
# Resolución de PWM por especie
# --------------------------------------------------------------------------
def test_resolve_pwm_per_species():
    ev = _load(ROOT / "scripts" / "sigma" / "51_eval_context.py", "eval_ctx_pwm")
    assert ev.resolve_pwm(DATA, "subtilis", None).name == "promoter_pwms_subtilis_all.tsv"
    assert ev.resolve_pwm(DATA, "ecoli_k12", None).name == "promoter_pwms.tsv"
    assert ev.is_ecoli_transfer("subtilis", DATA / "regulondb" / "gt" / "promoter_pwms.tsv")
    assert not ev.is_ecoli_transfer("subtilis",
                                    DATA / "dbtbs" / "promoter_pwms_subtilis_all.tsv")


# --------------------------------------------------------------------------
# Bug --pwm de 50_eval_layers.py
# --------------------------------------------------------------------------
def test_eval_layers_honours_pwm_flag(tmp_path: Path, monkeypatch):
    mod = _load(ROOT / "scripts" / "sigma" / "50_eval_layers.py", "eval_layers_pwm")
    gold = tmp_path / "gold" / "toy" / "tss_sigma.tsv"
    _write_tsv(gold, ["contig", "pos", "strand", "sigma"],
               [["c", "10", "+", "sigma70"]])
    pwm = tmp_path / "pwm.tsv"
    pwm.write_text("sigma\tbox\tn_sites\twidth\tconsensus\tinfo_bits\tpwm\n",
                   encoding="utf-8")
    captured: list[Path | None] = []
    fake = {
        "genome": "toy", "source": "gold",
        "config": {"window": [-80, 20], "decoys": 200, "seed": 17,
                   "zmin": 2.0, "pmax": 0.05},
        "layer0_data": {"n_gold_tss": 1, "n_with_sequence": 1, "n_unlabeled": 0,
                        "n_no_pwm": 0, "n_missing_sequence": 0, "coverage_sequence": 1.0},
        "layer1_localization": {"n": 1, "n_supported": 1, "frac_supported": 1.0},
        "layer2_motif": {},
        "layer3_assignment": {
            "n_evaluable": 1, "n_single": 1, "n_multi": 0,
            "top1_accuracy": 1.0, "top3_any": 1.0, "top3_all": 1.0,
            "single_top1": 1.0, "multi_top1_any": None, "multi_exact_top1": None,
            "majority_class": "sigma70", "majority_baseline": 1.0, "per_sigma": {}},
        "layer4_context": {"sigma_promoters_available": False,
                           "n_with_context_gene": 0, "frac_with_context_gene": 0.0},
        "classes": {}, "rows": [], "sigma_features": [],
    }

    def fake_evaluate(*args, **kwargs):
        captured.append(mod.PWMS_PATH)
        return fake

    monkeypatch.setattr(mod, "evaluate", fake_evaluate)
    code = mod.main(["--genome", "toy", "--source", "gold",
                     "--gold-root", str(tmp_path / "gold"),
                     "--out", str(tmp_path / "out"), "--pwm", str(pwm)])
    assert code == 0
    assert captured and Path(captured[0]) == pwm


# --------------------------------------------------------------------------
# Integración sintética de 51_eval_context.py
# --------------------------------------------------------------------------
def _write_fasta(path: Path, name: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{seq}\n", encoding="utf-8")


def test_eval_context_synthetic_end_to_end(tmp_path: Path):
    ev = _load(ROOT / "scripts" / "sigma" / "51_eval_context.py", "eval_ctx_e2e")
    data = tmp_path / "data"
    genome = "toy"
    chrom = "c"
    # 400 pb: motivo ACGT en posiciones conocidas.
    seq = ("A" * 400)
    _write_fasta(data / "raw" / f"{genome}_genomic.fna", chrom, seq)
    (data / "regulation" / genome).mkdir(parents=True, exist_ok=True)
    (data / "regulation" / genome / "genes_meta.json").write_text(json.dumps({
        "seq_source": {"kind": "fna", "path": str(data / "raw" / f"{genome}_genomic.fna")}
    }), encoding="utf-8")
    _write_tsv(data / "regulation" / genome / "genes.tsv",
               ["contig", "start", "end", "strand", "gene_id", "gene_name",
                "product", "ec"],
               [["c", "100", "200", "+", "g1", "a", "enzyme", "2.7.1.1"],
                ["c", "250", "350", "+", "g2", "b", "other", ""]])
    _write_tsv(data / "regulation" / genome / "operons.tsv",
               ["operon_id", "contig", "start", "end", "strand", "n_genes"],
               [["op1", "c", "90", "200", "+", "2"]])
    _write_tsv(data / "sigma" / "gold" / genome / "tss_sigma.tsv",
               ["species", "genome_id", "contig", "pos", "strand", "sigma", "gene"],
               [["Toy", genome, "c", "90", "+", "sigma70", "a"],
                ["Toy", genome, "c", "95", "+", "sigma70", "a"],
                ["Toy", genome, "c", "150", "+", "sigma32", "a"],
                ["Toy", genome, "c", "260", "+", "sigma32", "b"]])
    pwm = data / "regulondb" / "gt" / "promoter_pwms.tsv"
    pwm.parent.mkdir(parents=True, exist_ok=True)
    matrix = json.dumps([[1.0, -1.0, -1.0, -1.0], [-1.0, 1.0, -1.0, -1.0],
                         [-1.0, -1.0, 1.0, -1.0], [-1.0, -1.0, -1.0, 1.0]])
    pwm.write_text(
        "sigma\tbox\tn_sites\twidth\tconsensus\tinfo_bits\tpwm\n"
        f"sigma70\t-10\t10\t4\tACGT\t2.0\t{matrix}\n"
        f"sigma32\t-10\t10\t4\tACGT\t2.0\t{matrix}\n",
        encoding="utf-8")
    result = ev.evaluate_genome(
        data, genome, tmp_path / "out", decoys=5, folds=2, scales=[0.0, 1.0],
        min_support=1, min_gain=0.0)
    assert result["data"]["n_evaluable"] == 4
    out = tmp_path / "out" / genome
    for name in ("context_metrics.json", "context_sites.tsv",
                 "context_predictions.tsv", "context_evidence.tsv",
                 "context_evidence_coverage.tsv", "context_report.md"):
        assert (out / name).exists(), name
    metrics = json.loads((out / "context_metrics.json").read_text())
    assert metrics["baseline"]["top1"] is not None
