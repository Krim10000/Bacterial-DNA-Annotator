"""Contrato y funcionamiento de los módulos autodescubiertos."""
from __future__ import annotations

from pathlib import Path

import pytest

from modules import discover, get
from modules.base import RunContext


def test_registry_discovers_modules():
    mods = discover()
    ids = [m.id for m in mods]
    assert len(mods) >= 10
    assert len(ids) == len(set(ids))
    assert {"annotation", "metabolism", "regulation", "structure",
            "mobile", "ghost", "taxonomy", "membrane", "targets",
            "validation"} <= set(ids)


def test_every_module_declares_unique_options_and_defaults():
    for mod in discover():
        keys = [o.key for o in mod.options]
        assert len(keys) == len(set(keys)), mod.id
        cfg = mod.config(mod.defaults())
        assert set(cfg) == set(keys)
        assert isinstance(mod.name, str) and mod.name
        assert isinstance(mod.summary, str) and mod.summary


def test_option_validation_rejects_out_of_range():
    ann = get("annotation")
    assert ann is not None
    with pytest.raises(ValueError):
        ann.config({"min_aa": 5})       # por debajo del mínimo
    with pytest.raises(ValueError):
        ann.config({"no_existe": 1})


def test_annotation_module_runs_on_demo(tmp_path: Path):
    mod = get("annotation")
    assert mod is not None
    demo = Path(__file__).resolve().parents[1] / "webapp" / "samples" / "demo.fna"
    ctx = RunContext(data_root=tmp_path, out_root=tmp_path / "out",
                     config=mod.config({"cds": True, "trna": False,
                                        "rrna": False, "ghost": False,
                                        "esm": False, "deserts": False}))
    res = mod.run(ctx, [("demo.fna", demo.read_bytes())])
    assert res.ok, res.error
    assert res.metrics["demo"]["cds"] > 50
    assert res.metrics["demo"]["coverage"] > 0.5
    assert any(p.exists() for p in res.outputs.values())


def test_selftest_is_fast_and_informative():
    ann = get("annotation")
    assert ann is not None
    msg = ann.selftest()
    assert "CDS" in msg
