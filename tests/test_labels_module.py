from pathlib import Path
from modules.labels.module import MODULE
from modules.base import RunContext


def test_labels_module_is_registered_after_regulation_before_sigma():
    assert MODULE.id == "labels"
    assert MODULE.order == 32
    assert MODULE.run_fn is not None


def test_labels_selftest():
    assert "OK" in MODULE.selftest()


def test_labels_run_k12_if_available():
    if not Path("data/regulation/ecoli_k12/genes.tsv").exists():
        import pytest; pytest.skip("sin K-12 local")
    ctx = RunContext(data_root=Path("data"), out_root=Path(".scratch/test_labels"),
                     config=MODULE.defaults())
    res = MODULE.run(ctx, "ecoli_k12")
    assert res.ok
    assert res.metrics["n_operons"] > 0
