"""Tests de la capa de regulación automática (operones/promotores por reglas).

Secuencias y features sintéticas: no se usan RegulonDB/DBTBS, herramientas
externas ni red. Se comprueba el contrato de columnas, la agrupación por
contig/hebra/hueco y el motivo -35/-10 en ambas hebras.
"""
from __future__ import annotations

import csv
from copy import deepcopy
from pathlib import Path

from webapp.regulation_layer import (
    DEFAULT_MAX_GAP,
    MINUS10,
    MINUS35,
    OPERON_FIELDS,
    PROMOTER_FIELDS,
    promoter_motif,
    revcomp,
    run,
)

# Ventana upstream canónica de 50 pb: -35 en 15..20, espaciador 17,
# -10 en 38..43 y 6 pb hasta el inicio del gen.
PROMOTER_TX = "A" * 15 + MINUS35 + "A" * 17 + MINUS10 + "A" * 6
assert len(PROMOTER_TX) == 50


def _features() -> list[dict]:
    return [
        {"contig": "chr", "type": "CDS", "start": 101, "end": 250,
         "strand": "+", "id": "cds_p"},
        {"contig": "chr", "type": "CDS", "start": 701, "end": 850,
         "strand": "-", "id": "cds_m"},
    ]


def _genome_with_promoters(length: int = 1000) -> str:
    """Contig sintético con un promotor en hebra + y otro en hebra -.

    Para la hebra - se guarda el reverso complementario: al extraer la ventana
    y reverse-complementarla debe recuperarse ``PROMOTER_TX``.
    """
    seq = list("A" * length)
    seq[50:100] = PROMOTER_TX                 # 1-based 51..100, gen + en 101
    seq[850:900] = revcomp(PROMOTER_TX)       # 1-based 851..900, gen - acaba en 850
    return "".join(seq)


def _read(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_grouping_same_contig_and_strand_by_gap(tmp_path: Path):
    feats = [
        {"contig": "chr", "type": "CDS", "start": 100, "end": 200,
         "strand": "+", "id": "g1"},
        {"contig": "chr", "type": "CDS", "start": 250, "end": 350,
         "strand": "+", "id": "g2"},   # hueco 49 -> mismo operón
        {"contig": "chr", "type": "CDS", "start": 1000, "end": 1100,
         "strand": "+", "id": "g3"},   # hueco 649 -> nuevo operón
        {"contig": "chr", "type": "CDS", "start": 1200, "end": 1300,
         "strand": "+", "id": "g4"},   # hueco 99 -> mismo operón
        {"contig": "chr", "type": "tRNA", "start": 400, "end": 470,
         "strand": "+", "id": "t1"},   # no es CDS: se ignora
    ]
    snapshot = deepcopy(feats)
    res = run(feats, tmp_path, "syn", max_gap=DEFAULT_MAX_GAP)

    assert res["ok"] is True
    assert res["n_operons"] == 2
    assert res["n_promoters"] == 2
    assert feats == snapshot, "run no debe mutar las features"

    ops = res["rows"]
    assert [r["n_genes"] for r in ops] == [2, 2]
    assert [r["genes"] for r in ops] == ["g1;g2", "g3;g4"]
    assert [r["internal_gaps"] for r in ops] == ["49", "99"]
    assert [r["max_internal_gap"] for r in ops] == [49, 99]
    assert ops[0]["start"] == 100 and ops[0]["end"] == 350
    assert ops[0]["confidence"] == "medium"
    assert ops[1]["confidence"] == "medium"
    for row in ops:
        assert row["evidence"] == "rule_based"
        assert row["status"] == "hypothesis_not_experimental"
        assert row["confidence"] in ("low", "medium")
        assert "no experimental" in row["limitation"].lower()
        assert row["promoter_id"].endswith(row["operon_id"].split(":op:")[-1])

    tsv = tmp_path / "syn.operons_auto.tsv"
    assert tsv.exists()
    with tsv.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert tuple(reader.fieldnames) == OPERON_FIELDS
        assert len(list(reader)) == 2


def test_gap_boundary_overlap_and_singleton(tmp_path: Path):
    feats = [
        {"contig": "chr", "type": "CDS", "start": 100, "end": 200,
         "strand": "-", "id": "a"},
        {"contig": "chr", "type": "CDS", "start": 501, "end": 600,
         "strand": "-", "id": "b"},    # hueco exactamente 300 -> une
        {"contig": "chr", "type": "CDS", "start": 902, "end": 1001,
         "strand": "-", "id": "c"},    # hueco 301 -> corta
        {"contig": "chr", "type": "CDS", "start": 2000, "end": 2100,
         "strand": "-", "id": "d"},
        {"contig": "chr", "type": "CDS", "start": 2050, "end": 2150,
         "strand": "-", "id": "e"},    # solapamiento -> mismo operón
    ]
    res = run(feats, tmp_path, "edge", max_gap=300)
    assert res["ok"] is True
    assert res["n_operons"] == 3
    assert [r["n_genes"] for r in res["rows"]] == [2, 1, 2]
    # En hebra - los genes se listan en orden de transcripción (5'->3').
    assert [r["genes"] for r in res["rows"]] == ["b;a", "c", "e;d"]
    assert res["rows"][1]["max_internal_gap"] == ""
    assert res["rows"][1]["confidence"] == "low"
    assert res["rows"][1]["confidence_reason"] == "gen_unico_sin_vecino_en_la_hebra"
    assert res["rows"][2]["max_internal_gap"] == -51

    # El umbral es configurable: con 100 pb, a/b se separan.
    strict = run(feats, tmp_path, "strict", max_gap=100)
    assert strict["n_operons"] == 4
    assert [r["n_genes"] for r in strict["rows"]] == [1, 1, 1, 2]


def test_opposite_strands_are_not_grouped(tmp_path: Path):
    feats = [
        {"contig": "chr", "type": "CDS", "start": 100, "end": 200,
         "strand": "+", "id": "plus"},
        {"contig": "chr", "type": "CDS", "start": 150, "end": 250,
         "strand": "-", "id": "minus"},
        {"contig": "other", "type": "CDS", "start": 100, "end": 200,
         "strand": "+", "id": "other_contig"},
    ]
    res = run(feats, tmp_path, "strands", max_gap=1000)
    assert res["ok"] is True
    assert res["n_operons"] == 3
    assert all(r["n_genes"] == 1 for r in res["rows"])
    assert {(r["contig"], r["strand"]) for r in res["rows"]} == {
        ("chr", "+"), ("chr", "-"), ("other", "+")}
    assert res["n_promoters"] == 3


def test_promoter_motif_both_strands(tmp_path: Path):
    res = run(_features(), tmp_path, "prom",
              sequences={"chr": _genome_with_promoters()})
    assert res["ok"] is True
    assert res["n_operons"] == 2
    assert res["n_promoters"] == 2

    proms = {p["strand"]: p for p in res["promoter_rows"]}
    for strand, expected_gene in (("+", "cds_p"), ("-", "cds_m")):
        p = proms[strand]
        assert p["gene_id"] == expected_gene
        assert p["sequence_status"] == "available"
        assert p["box_found"] == "yes"
        assert p["box35"] == MINUS35
        assert p["box10"] == MINUS10
        assert p["spacer"] == 17
        assert float(p["score"]) == 1.0
        assert p["confidence"] == "medium"
        assert p["evidence"] == "rule_based"
        assert p["status"] == "hypothesis_not_experimental"
        assert "confirmado" in p["limitation"]

    assert (proms["+"]["start"], proms["+"]["end"]) == (51, 100)
    assert (proms["-"]["start"], proms["-"]["end"]) == (851, 900)
    assert proms["+"]["length"] == 50

    tsv = tmp_path / "prom.promoters_auto.tsv"
    with tsv.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert tuple(reader.fieldnames) == PROMOTER_FIELDS
        assert len(list(reader)) == 2


def test_no_motif_is_low_confidence(tmp_path: Path):
    feats = [{"contig": "plain", "type": "CDS", "start": 101, "end": 250,
              "strand": "+", "id": "plain_gene"}]
    res = run(feats, tmp_path, "plain", sequences={"plain": "A" * 1000})
    p = res["promoter_rows"][0]
    assert p["box_found"] == "no"
    assert float(p["score"]) == 0.0
    assert p["confidence"] == "low"
    assert p["confidence_reason"] == "sin_motivo_-35/-10_plausible"


def test_no_sequence_writes_nominal_promoter(tmp_path: Path):
    res = run(_features(), tmp_path, "noseq")  # sin ``sequences``
    assert res["ok"] is True
    assert res["n_operons"] == 2
    assert res["n_promoters"] == 2
    for p in res["promoter_rows"]:
        assert p["sequence_status"] == "no_sequence_available"
        assert p["sequence"] == ""
        assert p["box_found"] == "no"
        assert float(p["score"]) == 0.0
        assert p["confidence"] == "low"
        assert p["confidence_reason"] == "sin_secuencia_para_evaluar_motivo"
    # Coordenadas nominales de la ventana de 50 pb en ambas hebras.
    plus = next(p for p in res["promoter_rows"] if p["strand"] == "+")
    minus = next(p for p in res["promoter_rows"] if p["strand"] == "-")
    assert (plus["start"], plus["end"]) == (51, 100)
    assert (minus["start"], minus["end"]) == (851, 900)


def test_empty_features_write_headers(tmp_path: Path):
    res = run([], tmp_path, "empty")
    assert res["ok"] is True
    assert res["n_operons"] == 0 and res["n_promoters"] == 0
    for name, fields in (("empty.operons_auto.tsv", OPERON_FIELDS),
                         ("empty.promoters_auto.tsv", PROMOTER_FIELDS)):
        path = tmp_path / name
        assert path.exists()
        with path.open(encoding="utf-8") as fh:
            assert tuple(csv.DictReader(fh, delimiter="\t").fieldnames) == fields


def test_promoter_motif_helper():
    motif = promoter_motif(PROMOTER_TX)
    assert motif["box35"] == MINUS35
    assert motif["box10"] == MINUS10
    assert motif["spacer"] == 17
    assert motif["score"] == 1.0
    assert motif["box_found"] == "yes"
    # Un motivo fuera de la posición promotora plausible no debe puntuar.
    shifted = "A" * 30 + MINUS35 + "A" * 17 + MINUS10 + "A" * 3
    assert promoter_motif(shifted)["box_found"] == "no"
    assert promoter_motif("")["box_found"] == "no"
