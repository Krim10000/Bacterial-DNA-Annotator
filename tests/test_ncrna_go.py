from pathlib import Path
from modules.ncrna_go.module import MODULE
from modules.base import RunContext


def test_module_is_independent_and_before_labels():
    assert MODULE.id == "ncrna_go"
    assert MODULE.order == 29
    assert MODULE.run_fn is not None


def test_selftest():
    assert "OK" in MODULE.selftest()


def test_k12_rrna_trna_if_data():
    if not Path("data/annotations").exists():
        import pytest; pytest.skip("sin anotaciones")
    ctx = RunContext(data_root=Path("data"), out_root=Path(".scratch/test_ncrna"),
                     config=MODULE.defaults())
    res = MODULE.run(ctx, "ecoli_k12")
    assert res.ok
    assert res.metrics["n_rows"] > 0
    assert {"rRNA", "tRNA"} <= set(res.metrics["by_type"])
