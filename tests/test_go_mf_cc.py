from pathlib import Path
from modules.go_mf_cc.module import MODULE
from modules.base import RunContext


def test_module_registered_before_labels():
    assert MODULE.id == "go_mf_cc"
    assert MODULE.order == 31
    assert MODULE.run_fn is not None


def test_selftest():
    assert "OK" in MODULE.selftest()


def test_k12_has_mf_cc_if_index():
    if not Path("data/labels/go_mf_cc/goa.sqlite").exists():
        import pytest; pytest.skip("sin índice GO")
    ctx = RunContext(data_root=Path("data"), out_root=Path(".scratch/test_go"),
                     config=MODULE.defaults())
    res = MODULE.run(ctx, "ecoli_k12")
    assert res.ok
    assert res.metrics["n_mf"] > 0
    assert res.metrics["n_cc"] > 0
