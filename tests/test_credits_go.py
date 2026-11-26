from webapp import credits


def test_go_and_goa_are_credited():
    names = {c["name"]: c for c in credits.CREDITS}
    assert "Gene Ontology (GO)" in names
    assert "GOA / UniProt-GO" in names
    assert "CC BY 4.0" in names["Gene Ontology (GO)"]["license"]
    assert "2026-07-26" in names["Gene Ontology (GO)"]["cite"]
    assert "gene_go.tsv" in names["Gene Ontology (GO)"]["use"]
