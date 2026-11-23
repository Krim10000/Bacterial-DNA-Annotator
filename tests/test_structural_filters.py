"""Gating barato de filtros estructurales (sin GPU ni binarios externos)."""
from __future__ import annotations

from pathlib import Path

import webapp.structural_filters as sf


def _fake_cache(monkeypatch):
    store: dict = {}

    def cache(key, value=None):
        if value is None:
            return store.get(key)
        store[key] = value
        return value

    monkeypatch.setattr(sf, "_cache", cache)


def _pdb(tmp_path: Path) -> Path:
    p = tmp_path / "model.pdb"
    p.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n")
    return p


def test_low_quality_skips_expensive_filters(tmp_path, monkeypatch):
    _fake_cache(monkeypatch)
    calls: list[str] = []
    monkeypatch.setattr(sf, "foldseek_search", lambda p: {"hits": [{"pdb_id": "1abc"}]})
    monkeypatch.setattr(sf, "scan_pdbs", lambda *a, **k: {"hits": []})
    monkeypatch.setattr(sf, "summarize_enzymm", lambda h: [])
    monkeypatch.setattr(sf, "p2rank", lambda p: calls.append("p2rank") or {"pockets": []})
    monkeypatch.setattr(sf, "alphafill", lambda u: calls.append("alphafill") or {"hits": []})
    out = sf.run_filters(_pdb(tmp_path), uniprot="P00001", plddt=0.30, ptm=0.20)
    assert out["quality_gate"] == "low"
    assert out["gating"] == "low_quality"
    assert out["f3_pocket"] is False and out["f4_alphafill"] is False
    assert calls == []


def test_no_signal_skips_expensive_filters(tmp_path, monkeypatch):
    _fake_cache(monkeypatch)
    calls: list[str] = []
    monkeypatch.setattr(sf, "foldseek_search", lambda p: {"hits": []})
    monkeypatch.setattr(sf, "scan_pdbs", lambda *a, **k: {"hits": []})
    monkeypatch.setattr(sf, "summarize_enzymm", lambda h: [])
    monkeypatch.setattr(sf, "p2rank", lambda p: calls.append("p2rank") or {"pockets": []})
    out = sf.run_filters(_pdb(tmp_path), plddt=0.90, ptm=0.80)
    assert out["gating"] == "no_structural_signal"
    assert calls == []


def test_high_quality_with_signal_runs_f3(tmp_path, monkeypatch):
    _fake_cache(monkeypatch)
    calls: list[str] = []
    monkeypatch.setattr(sf, "foldseek_search", lambda p: {"hits": [{"pdb_id": "1abc"}]})
    monkeypatch.setattr(sf, "scan_pdbs", lambda *a, **k: {"hits": []})
    monkeypatch.setattr(sf, "summarize_enzymm", lambda h: [])
    monkeypatch.setattr(sf, "p2rank", lambda p: calls.append("p2rank") or {"pockets": []})
    out = sf.run_filters(_pdb(tmp_path), plddt=0.90, ptm=0.80)
    assert out["quality_gate"] == "ok"
    assert out["gating"] == "evaluado"
    assert calls == ["p2rank"]
