"""Contrato del descargador ComplexPortal sin red."""
from __future__ import annotations

from scripts.complexes.fetch_complexportal import gene_names, uniprot_map


def test_gene_names_and_uniprot(tmp_path):
    g = tmp_path / "genes.tsv"
    g.write_text("locus_tag\tgene_name\tacc\nb1\tatpA\tP1\n", encoding="utf-8")
    f = tmp_path / "functional.tsv"
    f.write_text("q1\tsp|P1|ATP1\t99\t100\t0\t0\t1\t100\t1\t100\t1e-9\t200\n",
                 encoding="utf-8")
    assert gene_names(g)["b1"] == "atpA"
    assert uniprot_map(f)["P1"] == "q1"
