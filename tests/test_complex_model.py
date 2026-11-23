"""Puente al modelo K-12 de scripts/complexes."""
from __future__ import annotations

from pathlib import Path

from webapp.complex_model import _write_genes, _write_operons, available


def test_unavailable_without_model(tmp_path: Path):
    ok, why = available(tmp_path)
    assert not ok
    assert "model" in why.lower() or "script" in why.lower()


def test_write_genes_and_operons(tmp_path: Path):
    feats = [
        {"id": "a", "type": "CDS", "contig": "c", "start": 100, "end": 500,
         "strand": "+"},
        {"id": "b", "type": "CDS", "contig": "c", "start": 700, "end": 900,
         "strand": "+"},
        {"id": "r", "type": "tRNA", "contig": "c", "start": 950, "end": 1000,
         "strand": "+"},
    ]
    hits = {"a": {"gene": "atpA", "ec": "", "uniprot": "P1"},
            "b": {"gene": "atpB", "ec": "", "uniprot": "P2"}}
    genes = tmp_path / "genes.tsv"
    ops = tmp_path / "operons.tsv"
    _write_genes(feats, hits, genes)
    _write_operons(feats, ops)
    text = genes.read_text()
    assert "atpA" in text and "P1" in text
    assert "cx_operon_" in ops.read_text()
