"""Tests de la auditoría de CDS solapados (``webapp.overlap_analysis``).

Features y secuencias sintéticas: no se usan genomas reales, herramientas
externas ni red. Se comprueban los cinco casos pedidos (mismo marco, marco
distinto, hebras opuestas, sin solapamiento y sin secuencia), el contrato TSV
y el resumen. La auditoría no debe mutar ni eliminar genes.
"""
from __future__ import annotations

import csv
import io
from copy import deepcopy

import pytest

from webapp.overlap_analysis import (
    FIELDS,
    FRAME_STATUSES,
    analyze,
    summary,
    to_tsv,
)

# Contig sintético suficientemente largo para todas las coordenadas usadas.
SEQUENCES = {
    "chr1": "ACGT" * 400,
    "chr2": "ACGT" * 400,
    "chr3": "ACGT" * 400,
}


def _cds(gene_id: str, start: int, end: int, strand: str = "+",
         contig: str = "chr1", **extra) -> dict:
    return {
        "contig": contig,
        "type": "CDS",
        "start": start,
        "end": end,
        "strand": strand,
        "id": gene_id,
        "product": "hypothetical protein",
        **extra,
    }


def test_same_frame_plus_strand():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 160, 500),  # inicio +60 pb -> mod 3 = 0
    ]
    snapshot = deepcopy(feats)

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    row = rows[0]
    assert row["contig"] == "chr1"
    assert (row["gene_a"], row["gene_b"]) == ("a", "b")
    assert (row["start_a"], row["end_a"], row["strand_a"]) == (100, 300, "+")
    assert (row["start_b"], row["end_b"], row["strand_b"]) == (160, 500, "+")
    assert row["strand_relation"] == "same"
    assert (row["overlap_start"], row["overlap_end"]) == (160, 300)
    assert row["overlap_bp"] == 141
    assert row["frame_shift"] == 0
    assert row["frame_status"] == "same_frame"
    assert (row["frame_reference"], row["frame_target"]) == ("a", "b")
    assert row["severity"] == "warning"
    assert "mismo marco" in row["warning"]
    assert "No se descarta" in row["warning"]
    assert feats == snapshot, "analyze no debe mutar las features"


def test_different_frame_plus_strand():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 161, 400),  # inicio +61 pb -> mod 3 = 1
    ]

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    row = rows[0]
    assert row["overlap_bp"] == 140
    assert row["frame_shift"] == 1
    assert row["frame_status"] == "different_frame"
    assert row["severity"] == "warning"
    assert "marcos distintos" in row["frame_reason"]


def test_minus_strand_frame_uses_five_prime_end():
    feats = [
        _cds("a", 100, 300, strand="-"),
        # En hebra - el 5' es ``end``: 450 - 300 = 150 -> mismo marco.
        # Un cálculo ingenuo con ``start`` daría 61 % 3 = 1.
        _cds("b", 161, 450, strand="-"),
    ]

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    row = rows[0]
    assert row["strand_relation"] == "same"
    assert row["overlap_bp"] == 140
    assert row["frame_reference"] == "b"
    assert row["frame_target"] == "a"
    assert row["frame_shift"] == 0
    assert row["frame_status"] == "same_frame"
    assert "hebra -" in row["frame_reason"]


@pytest.mark.parametrize("strand_a,strand_b", [("+", "-"), ("-", "+")])
def test_opposite_strands_are_antisense_not_same_frame(strand_a, strand_b):
    feats = [
        _cds("a", 100, 300, strand=strand_a),
        _cds("b", 200, 400, strand=strand_b),
    ]

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    row = rows[0]
    assert row["strand_relation"] == "opposite"
    assert row["overlap_bp"] == 101
    assert row["frame_shift"] == "unknown"
    assert row["frame_status"] == "opposite_strand"
    assert (row["frame_reference"], row["frame_target"]) == ("", "")
    assert row["severity"] == "info"
    assert "antisentido" in row["warning"]
    assert "opposite_strand" in FRAME_STATUSES


def test_no_overlap_no_rows():
    feats = [
        _cds("a", 100, 200),
        _cds("b", 201, 300),  # contiguo, no solapado
        _cds("other", 100, 300, contig="chr2"),  # otro contig
        _cds("minus", 400, 500, strand="-", contig="chr2"),  # no solapa
    ]

    assert analyze(feats, SEQUENCES) == []
    assert analyze(feats) == []


def test_without_sequences_frame_is_unknown():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 160, 500),
    ]

    rows = analyze(feats, None)

    assert len(rows) == 1
    row = rows[0]
    assert row["frame_shift"] == "unknown"
    assert row["frame_status"] == "unknown"
    assert row["severity"] == "warning"
    assert "sin secuencias" in row["frame_reason"]
    # Las coordenadas del solapamiento sí se auditan sin genoma.
    assert row["overlap_bp"] == 141


def test_sequences_do_not_cover_contig_or_coordinates():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 160, 500),
    ]

    missing = analyze(feats, {"other": "ACGT" * 400})
    assert missing[0]["frame_shift"] == "unknown"
    assert missing[0]["frame_status"] == "unknown"
    assert "ausente" in missing[0]["frame_reason"]

    short = analyze(feats, {"chr1": "A" * 250})
    assert short[0]["frame_shift"] == "unknown"
    assert "fuera de la secuencia" in short[0]["frame_reason"]


def test_to_tsv_contract():
    rows = analyze([_cds("a", 100, 300), _cds("b", 160, 500)], SEQUENCES)

    text = to_tsv(rows)
    parsed = list(csv.DictReader(io.StringIO(text), delimiter="\t"))

    assert tuple(parsed[0].keys()) == FIELDS
    assert len(parsed) == 1
    assert parsed[0]["contig"] == "chr1"
    assert parsed[0]["overlap_bp"] == "141"
    assert parsed[0]["frame_shift"] == "0"
    assert parsed[0]["frame_status"] == "same_frame"
    assert text.endswith("\n")

    empty = to_tsv([])
    assert empty.splitlines() == ["\t".join(FIELDS)]


def test_summary_counts_and_genes():
    feats = [
        # chr1: mismo marco
        _cds("a", 100, 300),
        _cds("b", 160, 400),
        # chr2: marco distinto
        _cds("c", 100, 300, contig="chr2"),
        _cds("d", 161, 330, contig="chr2"),
        # chr3: hebras opuestas
        _cds("e", 100, 300, contig="chr3"),
        _cds("f", 200, 400, strand="-", contig="chr3"),
    ]

    rows = analyze(feats, SEQUENCES)
    res = summary(rows)

    assert res["n_pairs"] == 3
    assert (res["n_same_strand"], res["n_opposite_strand"],
            res["n_unknown_strand"]) == (2, 1, 0)
    assert res["n_genes"] == 6
    assert res["genes"] == ["a", "b", "c", "d", "e", "f"]
    assert res["frame_same"] == 1
    assert res["frame_different"] == 1
    assert res["frame_opposite_strand"] == 1
    assert res["frame_unknown"] == 0
    assert (res["overlap_bp_total"], res["overlap_bp_max"]) == (141 + 140 + 101, 141)
    assert res["n_warnings"] == 2
    assert res["audit_only"] is True
    assert "ningún gen eliminado" in res["message"]


def test_summary_empty():
    res = summary([])
    assert res["n_pairs"] == 0
    assert res["n_genes"] == 0
    assert res["genes"] == []
    assert res["overlap_bp_total"] == 0
    assert res["overlap_bp_max"] == 0
    assert res["frame_unknown"] == 0
    assert res["audit_only"] is True


def test_ignores_non_cds_and_invalid_features():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 160, 500),
        {"contig": "chr1", "type": "tRNA", "start": 120, "end": 220,
         "strand": "+", "id": "t1"},
        _cds("bad_start", 0, 200),
        {"contig": "", "type": "CDS", "start": 100, "end": 300,
         "strand": "+", "id": "no_contig"},
        {"contig": "chr1", "type": "CDS", "start": "x", "end": 300,
         "strand": "+", "id": "bad_coord"},
        "not-a-mapping",
    ]

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    assert (rows[0]["gene_a"], rows[0]["gene_b"]) == ("a", "b")


def test_all_overlapping_pairs_reported_without_choosing():
    feats = [
        _cds("a", 100, 300),
        _cds("b", 150, 350),
        _cds("c", 200, 400),
    ]

    rows = analyze(feats, SEQUENCES)

    pairs = {(row["gene_a"], row["gene_b"]) for row in rows}
    assert pairs == {("a", "b"), ("a", "c"), ("b", "c")}
    res = summary(rows)
    assert res["n_genes"] == 3
    assert res["audit_only"] is True
    # El orden de salida es estable.
    assert rows == analyze(feats, SEQUENCES)


def test_unknown_strand_relation_is_flagged():
    feats = [
        _cds("a", 100, 300, strand="."),
        _cds("b", 160, 500, strand="+"),
    ]

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    assert rows[0]["strand_relation"] == "unknown"
    assert rows[0]["frame_shift"] == "unknown"
    assert rows[0]["frame_status"] == "unknown"
    assert "hebra desconocida" in rows[0]["frame_reason"]


def test_numeric_strings_are_normalised_and_input_not_mutated():
    feats = [
        {"contig": "chr1", "type": "CDS", "start": "100", "end": "300",
         "strand": "+", "id": "a"},
        {"contig": "chr1", "type": "CDS", "start": "160", "end": "500",
         "strand": "1", "id": "b"},
    ]
    snapshot = deepcopy(feats)

    rows = analyze(feats, SEQUENCES)

    assert len(rows) == 1
    assert rows[0]["overlap_bp"] == 141
    assert rows[0]["frame_shift"] == 0
    assert feats == snapshot
