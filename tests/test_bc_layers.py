"""Contratos de B estructural y C dominios."""
from __future__ import annotations

from pathlib import Path

from webapp.domain_layer import available as domains_available
from webapp.structural_b import run as run_b


def test_domains_require_panel(tmp_path: Path):
    ok, why = domains_available(tmp_path)
    assert not ok
    assert "pfam" in why.lower() or "hmmscan" in why.lower()


def test_b_structural_manifest(tmp_path: Path):
    data = tmp_path / "data"
    afdb = data / "metabolism" / "targets" / "afdb25k"
    afdb.mkdir(parents=True)
    (afdb / "P1.pdb").write_text("ATOM\n", encoding="utf-8")
    rows = [
        {"protein_id": "p1", "path": "B", "uniprot": "P1"},
        {"protein_id": "p2", "path": "B", "uniprot": "P2"},
        {"protein_id": "p3", "path": "C"},
    ]
    out = run_b(rows, tmp_path / "out", "g", data)
    assert out["ok"]
    assert out["n_available"] == 1
    assert rows[0]["structural_status"] == "afdb_model_available"
    assert rows[1]["structural_status"] == "no_afdb_model"
    assert "p1" in Path(out["manifest"]).read_text()
