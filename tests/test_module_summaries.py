"""Ejecuciones de resumen de módulos que leen datos ya calculados."""
from __future__ import annotations

from pathlib import Path

import pytest

from modules import get
from modules.base import RunContext


def test_regulation_summarizes_demo(tmp_path: Path):
    mod = get("regulation")
    assert mod is not None and mod.input_kind == "genome"
    ctx = RunContext(data_root=tmp_path, out_root=tmp_path / "out",
                     config=mod.defaults())
    res = mod.run(ctx, "demo")
    assert res.ok, res.error
    assert res.metrics["genes"] > 100
    assert res.metrics["operons"] > 10
    assert any(p.exists() for p in res.outputs.values())


def test_validation_summary(tmp_path: Path):
    mod = get("validation")
    assert mod is not None
    gt = Path(__file__).resolve().parents[1] / "data" / "regulondb" / "gt" / "gt_operons.tsv"
    if not gt.exists():
        pytest.skip("sin verdad de referencia local")
    ctx = RunContext(data_root=gt.parents[2], out_root=tmp_path / "out",
                     config=mod.defaults())
    res = mod.run(ctx)
    assert res.ok, res.error
    assert res.metrics["rows"]["operons"] > 0
