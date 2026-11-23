"""Construcción de GenBank y estados de oriC para genomas nuevos."""
from __future__ import annotations

from pathlib import Path

import webapp.oric as oric


def test_gbk_uses_accession_version(tmp_path: Path):
    seq = "ATG" * 300
    out = tmp_path / "x.gbk"
    name = oric._write_gbk(
        {"contig_1": seq},
        [{"contig": "contig_1", "type": "CDS", "start": 1, "end": 300,
          "strand": "+", "id": "cds_1"}],
        {"cds_1": {"gene": "dnaA", "product": "DnaA"}}, out)
    text = out.read_text()
    assert name == "contig_1"
    assert "contig_1.1" in text
    assert "dnaA" in text


def test_run_without_model_is_unavailable(tmp_path: Path):
    res = {"genome": "x", "sequences": {"c": "ATG" * 500},
           "features": [], "functional": {"hits": {}}}
    out = oric.run(res, tmp_path)  # sin modelo en tmp_path
    assert out["status"] == "unavailable"
    assert "modelo" in out["reason"].lower() or "orcapy" in out["reason"].lower()


def test_doric_lookup(tmp_path: Path):
    d = tmp_path / "doric"
    d.mkdir()
    (d / "x.csv").write_text(
        'refseq,nc,orics\nGCF_000001.1,NC_000001.1,"[{""oric_start"": 100, '
        '""oric_end"": 300}]"\n', encoding="utf-8")
    assert oric.doric_position(tmp_path, "GCF_000001.1") == 200
    assert oric.doric_position(tmp_path, "NC_999") is None
