from pathlib import Path
from modules.go_bp.module import MODULE
from modules.base import RunContext


def test_bp_module_order_and_source():
    assert MODULE.id == "go_bp"
    assert MODULE.order == 33
    assert MODULE.run_fn is not None


def test_selftest():
    assert "OK" in MODULE.selftest()


def test_k12_bp_if_data():
    if not Path("data/labels/external/go/go.obo").exists():
        import pytest; pytest.skip("sin go.obo")
    ctx = RunContext(data_root=Path("data"), out_root=Path(".scratch/test_bp"),
                     config=MODULE.defaults())
    res = MODULE.run(ctx, "ecoli_k12")
    assert res.ok
    assert res.metrics["n_bp_rows"] > 0
    assert res.metrics["source"].startswith("go.obo")
