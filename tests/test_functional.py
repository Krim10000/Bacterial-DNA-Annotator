"""Anotación funcional y tres caminos de proteínas."""
from __future__ import annotations

from pathlib import Path

from webapp.functional import gene_positions, parse_hit, parse_m8


def test_paths_a_b_c():
    a = parse_hit("p1", "sp|P03004|DNAA_ECOLI",
                  "Chromosomal replication initiator protein DnaA "
                  "OS=Escherichia coli OX=83333 GN=dnaA PE=1 SV=4",
                  99.0, 400, 0.0, 800.0)
    b = parse_hit("p2", "sp|Q12345|YXXX_ECOLI",
                  "Uncharacterized protein Yxxx OS=Escherichia coli GN=yxxx",
                  40.0, 200, 1e-20, 100.0)
    c = parse_hit("p3", "sp|Q99999|DUF1234_ECOLI",
                  "DUF1234 domain-containing protein OS=Escherichia coli",
                  30.0, 150, 1e-10, 50.0)
    assert a.path == "A" and a.gene == "dnaA" and a.uniprot == "P03004"
    assert b.path == "B"
    assert c.path == "B"


def test_parse_m8_and_gene_positions(tmp_path: Path):
    p = tmp_path / "hits.tsv"
    p.write_text(
        "cds_1\tsp|P03004|DNAA_ECOLI\t100\t400\t0.0\t800\t"
        "DnaA OS=Escherichia coli GN=dnaA\n", encoding="utf-8")
    hits = parse_m8(p)
    assert hits["cds_1"].gene == "dnaA"
    feats = [{"id": "cds_1", "start": 1398}, {"id": "cds_2", "start": 2000}]
    pos = gene_positions(feats, {k: {"gene": v.gene} for k, v in hits.items()})
    assert pos == {"dnaa": 1398}


def test_subunit_hint_is_only_a_name():
    from webapp.functional import _row_for
    h = parse_hit("p1", "sp|P0ABB0|ATPA_ECOLI",
                  "ATP synthase subunit alpha OS=Escherichia coli GN=atpA",
                  99.0, 400, 0.0, 800.0)
    row = _row_for({"id": "p1", "aa": "M" * 400}, h)
    assert row["subunit_hint"] == "named_subunit"
    assert row["subunit_evidence"] == "product_name_only"


def test_interpretable_labels():
    from webapp.functional import _row_for
    a = parse_hit("p1", "sp|P03004|DNAA_ECOLI",
                  "DnaA OS=Escherichia coli GN=dnaA", 99.0, 400, 0.0, 800.0)
    row = _row_for({"id": "p1", "aa": "M" * 400}, a)
    assert row["label"] == "A_named_function"
    assert row["evidence_strength"] == "strong"
    assert row["evidence_toward"] == "named_function"
    assert row["subunit_hint"] == "none"
    b = parse_hit("p2", "sp|Q1|DUF1_ECOLI",
                  "DUF1 domain-containing protein", 40.0, 100, 1e-5, 50.0)
    row = _row_for({"id": "p2", "aa": "M" * 400}, b)
    assert row["label"] == "B_uncharacterized_homolog"
    assert row["evidence_strength"] == "weak"
    c = _row_for({"id": "p3", "aa": "M" * 100}, None)
    assert c["label"] == "C_no_swissprot_hit"
    assert c["evidence_strength"] == "none"
