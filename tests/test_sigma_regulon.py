"""Pruebas de la capa sigma -> operón -> regulón.

Solo funciones puras de join/mapeo y un round-trip TSV en `tmp_path`.
No dependen de los 400 genomas ni de datos pesados.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "sigma" / "40_build_sigma_regulon.py"

_spec = importlib.util.spec_from_file_location("sigma_regulon_builder", MODULE_PATH)
assert _spec is not None and _spec.loader is not None
SR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SR)


# --------------------------------------------------------------------------
# Fixtures sintéticas
# --------------------------------------------------------------------------
@pytest.fixture
def sigma_row():
    return {
        "genome": "demo",
        "operon_id": "demo:1:+",
        "operon_label": "abc",
        "sigma": "sigma32",
        "sigma_class": "stress",
        "confidence_label": "sospecha_media",
        "confidence_score": "0.5",
        "sigma_status": "assigned",
        "score": "9.0",
        "first_gene_name": "abcA",
        "pathway_support": "yes",
        "pathway_hits": '["abcA:chaperone", "abcB:protease"]',
        "evidence": "PWM sigma32 z=2.2",
        "detected_by": "promoter_pwm:sigma32",
        "provenance": '{"window": 180}',
    }


@pytest.fixture
def table():
    operons = [{"operon_id": "demo:1:+", "label": "abc",
                "gene_names": "abcA,abcB"}]
    operon_genes = [
        {"operon_id": "demo:1:+", "position": "1", "gene_id": "g1",
         "gene_name": "abcA", "contig": "c1", "start": "10", "end": "20"},
        {"operon_id": "demo:1:+", "position": "2", "gene_id": "g2",
         "gene_name": "abcB", "contig": "c1", "start": "30", "end": "40"},
    ]
    genes = [
        {"gene_id": "g1", "gene_name": "abcA", "product": "chaperone A",
         "ec": "1.1.1.1"},
        {"gene_id": "g2", "gene_name": "abcB", "product": "protease B"},
    ]
    return {
        "operons": operons,
        "operon_genes": operon_genes,
        "genes": genes,
    }


def _build(genome, sigma_rows, table, edges=None, tf_available=True, rhea_index=None):
    operon_index = SR.build_operon_index(table["operons"])
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    tf_index, stats = SR.build_tf_index(edges or [], aliases)
    rows = SR.build_rows(
        genome, sigma_rows, operon_index, operon_gene_index, gene_index,
        tf_index, tf_available=tf_available, rhea_index=rhea_index or {},
    )
    return rows, stats


# --------------------------------------------------------------------------
# Esquema y campos derivados
# --------------------------------------------------------------------------
def test_output_schema_has_required_columns():
    required = {
        "genome", "operon_id", "operon_label", "sigma", "sigma_class",
        "confidence_label", "sigma_status", "target_genes", "target_functions",
        "tf_overlap", "tf_names", "pathway_context", "evidence",
        "detected_by", "provenance",
    }
    assert required <= set(SR.OUTPUT_COLS)
    assert len(SR.OUTPUT_COLS) == len(set(SR.OUTPUT_COLS))


def test_parse_pathway_hits():
    assert SR.parse_pathway_hits('["a:x", "b:y", "a:x"]') == ["a:x", "b:y"]
    assert SR.parse_pathway_hits("") == []
    assert SR.parse_pathway_hits("[]") == []
    assert SR.parse_pathway_hits("texto suelto") == ["texto suelto"]
    assert SR.parse_pathway_hits('"una cadena"') == ["una cadena"]
    assert SR.parse_pathway_hits('["a\\tb"]') == ["a b"]


def test_pathway_context_variants():
    assert SR.pathway_context(
        {"pathway_support": "yes", "pathway_hits": '["a:x"]'}
    ) == "yes: a:x"
    assert SR.pathway_context(
        {"pathway_support": "yes", "pathway_hits": "[]"}
    ) == "yes"
    assert SR.pathway_context(
        {"pathway_support": "no", "pathway_hits": '["a:x"]'}
    ) == "no"
    assert SR.pathway_context(
        {"pathway_support": "", "pathway_hits": "[]"}
    ) == "not_available"
    assert SR.pathway_context(
        {"pathway_support": "", "pathway_hits": '["a:x"]'}
    ) == "yes: a:x"


# --------------------------------------------------------------------------
# Índices
# --------------------------------------------------------------------------
def test_build_operon_gene_index_orders_and_dedupes():
    rows = [
        {"operon_id": "op", "position": "2", "gene_id": "g2", "gene_name": "b"},
        {"operon_id": "op", "position": "1", "gene_id": "g1", "gene_name": "a"},
        {"operon_id": "op", "position": "2", "gene_id": "g2", "gene_name": "b"},
        {"operon_id": "", "position": "1", "gene_id": "gx", "gene_name": "x"},
    ]
    index = SR.build_operon_gene_index(rows)
    assert [g["gene_name"] for g in index["op"]] == ["a", "b"]
    assert "" not in index


def test_build_gene_index_completes_missing_fields():
    rows = [
        {"gene_id": "g1", "gene_name": "a", "product": "", "ec": ""},
        {"gene_id": "g1", "gene_name": "", "product": "producto", "ec": "2.2.2.2"},
    ]
    index = SR.build_gene_index(rows)
    assert index["g1"] == {
        "gene_name": "a", "locus_tag": "", "product": "producto", "ec": "2.2.2.2",
    }


def test_build_operon_aliases_includes_lowercase(table):
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    assert aliases["abcA"] == {"demo:1:+"}
    assert aliases["abca"] == {"demo:1:+"}
    assert aliases["g1"] == {"demo:1:+"}


# --------------------------------------------------------------------------
# Funciones de destino
# --------------------------------------------------------------------------
def test_target_functions_with_ec_rhea_and_cap():
    genes = [
        {"gene_id": "g1", "gene_name": "a"},
        {"gene_id": "g2", "gene_name": "b"},
    ]
    gene_index = {
        "g1": {"gene_name": "a", "locus_tag": "", "product": "enzima A",
               "ec": "1.1.1.1"},
        "g2": {"gene_name": "b", "locus_tag": "", "product": "", "ec": ""},
    }
    rhea = {
        "by_gene": {"g1": ["RHEA:1", "RHEA:2", "RHEA:3"]},
        "by_coord": {("c1", "30", "40"): ["RHEA:9"]},
    }
    text = SR.build_target_functions(
        [
            genes[0],
            {"gene_id": "g2", "gene_name": "b", "contig": "c1",
             "start": "30", "end": "40"},
        ],
        gene_index,
        rhea,
    )
    assert "a:enzima A [EC 1.1.1.1] [RHEA:1,RHEA:2,+1]" in text
    assert "b [RHEA:9]" in text


def test_target_functions_truncates_explicitly():
    genes = [{"gene_id": f"g{i}", "gene_name": f"g{i}"} for i in range(12)]
    gene_index = {
        f"g{i}": {"gene_name": f"g{i}", "locus_tag": "", "product": f"p{i}",
                  "ec": ""}
        for i in range(12)
    }
    text = SR.build_target_functions(genes, gene_index, {})
    parts = text.split("; ")
    assert len(parts) == SR.MAX_FUNCTION_ENTRIES + 1
    assert parts[-1] == "(+4 más)"


def test_target_functions_not_available():
    assert SR.build_target_functions([], {}, {}) == "not_available"


# --------------------------------------------------------------------------
# TF -> operón
# --------------------------------------------------------------------------
def test_tf_index_uses_target_operon(table):
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    index, stats = SR.build_tf_index(
        [{"tf": "LexA", "target_operon": "demo:1:+", "target_genes": ""}],
        aliases,
    )
    assert [e["tf"] for e in index["demo:1:+"]] == ["LexA"]
    assert stats["n_tf_edges_mapped"] == 1


def test_tf_index_falls_back_to_target_genes(table):
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    index, stats = SR.build_tf_index(
        [{"tf": "CRP", "target_operon": "", "target_genes": "g2"}],
        aliases,
    )
    assert index["demo:1:+"][0]["tf"] == "CRP"
    assert stats["n_tf_edges_mapped"] == 1


def test_tf_index_ignores_ambiguous_and_unmapped(table):
    table["operon_genes"] = table["operon_genes"] + [
        {"operon_id": "demo:2:+", "position": "1", "gene_id": "g3",
         "gene_name": "abcA", "contig": "c1", "start": "50", "end": "60"},
    ]
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    index, stats = SR.build_tf_index(
        [
            {"tf": "Amb", "target_operon": "", "target_genes": "abcA"},
            {"tf": "NoMap", "target_operon": "", "target_genes": "desconocido"},
        ],
        aliases,
    )
    assert index == {}
    assert stats["n_tf_edges_ambiguous"] == 1
    assert stats["n_tf_edges_unmapped"] == 1


def test_tf_index_without_tf(table):
    operon_gene_index = SR.build_operon_gene_index(table["operon_genes"])
    gene_index = SR.build_gene_index(table["genes"])
    aliases = SR.build_operon_aliases(operon_gene_index, gene_index)
    index, stats = SR.build_tf_index(
        [{"tf": "", "target_operon": "demo:1:+"}], aliases
    )
    assert index == {}
    assert stats["n_tf_edges_without_tf"] == 1


# --------------------------------------------------------------------------
# Filas completas
# --------------------------------------------------------------------------
def test_build_rows_with_tf_and_pathway(sigma_row, table):
    rows, _ = _build(
        "demo", [sigma_row], table,
        edges=[{"tf": "LexA", "target_operon": "demo:1:+", "target_genes": ""}],
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["target_genes"] == "abcA,abcB"
    assert "abcA:chaperone A [EC 1.1.1.1]" in row["target_functions"]
    assert row["tf_overlap"] == 1
    assert row["tf_names"] == "LexA"
    assert row["pathway_context"] == "yes: abcA:chaperone; abcB:protease"
    assert row["sigma_class"] == "stress"
    assert row["confidence_label"] == "sospecha_media"
    assert row["sigma_status"] == "assigned"


def test_build_rows_without_regulon_marks_not_available(sigma_row, table):
    rows, _ = _build("demo", [sigma_row], table, tf_available=False)
    assert rows[0]["tf_overlap"] == "not_available"
    assert rows[0]["tf_names"] == "not_available"


def test_build_rows_falls_back_to_operon_gene_names(sigma_row):
    table = {
        "operons": [{"operon_id": "demo:1:+", "label": "abc",
                     "gene_names": "abcA,abcB"}],
        "operon_genes": [],
        "genes": [],
    }
    rows, _ = _build("demo", [sigma_row], table, tf_available=False)
    assert rows[0]["target_genes"] == "abcA,abcB"
    assert rows[0]["target_functions"] == "abcA; abcB"


def test_build_rows_uses_multiple_tfs_sorted(sigma_row, table):
    rows, _ = _build(
        "demo", [sigma_row], table,
        edges=[
            {"tf": "LexA", "target_operon": "demo:1:+"},
            {"tf": "CRP", "target_operon": "demo:1:+"},
            {"tf": "LexA", "target_operon": "demo:1:+"},
        ],
    )
    assert rows[0]["tf_overlap"] == 2
    assert rows[0]["tf_names"] == "CRP; LexA"


def test_build_rows_skips_empty_operon_id(sigma_row, table):
    bad = dict(sigma_row, operon_id="")
    rows, _ = _build("demo", [bad], table, tf_available=False)
    assert rows == []


# --------------------------------------------------------------------------
# Deduplicación y resumen
# --------------------------------------------------------------------------
def test_dedupe_sigma_rows_prefers_assigned():
    unresolved = {
        "operon_id": "op", "sigma_status": "unresolved",
        "confidence_score": "0.99", "score": "99",
    }
    assigned = {
        "operon_id": "op", "sigma_status": "assigned",
        "confidence_score": "0.1", "score": "1",
    }
    kept, dropped = SR._dedupe_sigma_rows([unresolved, assigned])
    assert dropped == 1
    assert kept == [assigned]


def test_dedupe_sigma_rows_keeps_order():
    first = {"operon_id": "op1", "sigma_status": "assigned"}
    second = {"operon_id": "op2", "sigma_status": "assigned"}
    kept, dropped = SR._dedupe_sigma_rows([first, second])
    assert dropped == 0
    assert kept == [first, second]


def test_summarize_counts(sigma_row, table):
    rows, stats = _build(
        "demo", [sigma_row], table,
        edges=[{"tf": "LexA", "target_operon": "demo:1:+"}],
    )
    without = {
        "sigma": "sigma70", "sigma_class": "housekeeping",
        "confidence_label": "sin_sospecha", "sigma_status": "ambiguous",
        "pathway_context": "no", "tf_overlap": 0, "tf_names": "",
    }
    unavailable = {
        "sigma": "sigma70", "sigma_class": "housekeeping",
        "confidence_label": "sin_sospecha", "sigma_status": "ambiguous",
        "pathway_context": "not_available", "tf_overlap": "not_available",
        "tf_names": "not_available",
    }
    rows = rows + [without, unavailable]
    summary = SR.summarize(
        "demo", rows, table["operons"], sigma_rows_input=3,
        duplicate_count=0, tf_available=True, tf_stats=stats, inputs={},
    )
    assert summary["n_edges"] == 3
    assert summary["n_operons_total"] == 1
    assert summary["n_operons_with_tf_overlap"] == 1
    assert summary["n_operons_sin_tf"] == 1
    assert summary["n_operons_tf_not_available"] == 1
    assert summary["sigma_counts"] == {"sigma32": 1, "sigma70": 2}
    assert summary["confidence_counts"] == {
        "sospecha_media": 1, "sin_sospecha": 2,
    }
    assert summary["pathway_context_counts"] == {
        "no": 1, "not_available": 1, "yes": 1,
    }
    assert summary["n_tfs"] == 1
    assert summary["tfs"] == ["LexA"]


# --------------------------------------------------------------------------
# E/S y autoprueba
# --------------------------------------------------------------------------
def test_write_tsv_is_raw_and_roundtrips(tmp_path):
    path = tmp_path / "tabla.tsv"
    rows = [{"a": "x\ty", "b": '{"k": "v"}', "c": "línea\nnueva"}]
    SR._write_tsv(path, rows, ["a", "b", "c"])
    text = path.read_text(encoding="utf-8")
    assert '{"k": "v"}' in text
    assert "\t".join(["a", "b", "c"]) in text
    read = SR._read_tsv(path)
    assert read == [{"a": "x y", "b": '{"k": "v"}', "c": "línea nueva"}]


def test_build_genome_without_regulon_end_to_end(tmp_path):
    data = tmp_path / "data"
    reg = data / "regulation" / "demo"
    sig = data / "sigma" / "demo"
    reg.mkdir(parents=True)
    sig.mkdir(parents=True)
    SR._write_tsv(sig / "sigma_promoters.tsv", [{
        "genome": "demo", "operon_id": "demo:1:+", "operon_label": "abc",
        "sigma": "sigma70", "sigma_class": "housekeeping",
        "confidence_label": "sospecha_media", "sigma_status": "assigned",
        "pathway_support": "no", "pathway_hits": "[]", "evidence": "ev",
        "detected_by": "pwm", "provenance": '{"x": 1}',
    }], [
        "genome", "operon_id", "operon_label", "sigma", "sigma_class",
        "confidence_label", "sigma_status", "pathway_support", "pathway_hits",
        "evidence", "detected_by", "provenance",
    ])
    SR._write_tsv(reg / "operons.tsv", [{
        "operon_id": "demo:1:+", "label": "abc", "gene_names": "a,b",
    }], ["operon_id", "label", "gene_names"])
    SR._write_tsv(reg / "operon_genes.tsv", [
        {"operon_id": "demo:1:+", "position": "1", "gene_id": "g1",
         "gene_name": "a", "contig": "c", "start": "1", "end": "10"},
        {"operon_id": "demo:1:+", "position": "2", "gene_id": "g2",
         "gene_name": "b", "contig": "c", "start": "20", "end": "30"},
    ], ["operon_id", "position", "gene_id", "gene_name", "contig", "start", "end"])
    SR._write_tsv(reg / "genes.tsv", [
        {"gene_id": "g1", "gene_name": "a", "product": "enzima A", "ec": ""},
        {"gene_id": "g2", "gene_name": "b", "product": "", "ec": ""},
    ], ["gene_id", "gene_name", "product", "ec"])

    out = tmp_path / "out"
    assert SR.discover_genomes(data) == ["demo"]
    summary = SR.build_genome(data, "demo", out)
    assert summary["tf_available"] is False
    assert summary["n_operons_tf_not_available"] == 1
    rows = SR._read_tsv(out / "demo" / "sigma_regulon.tsv")
    assert len(rows) == 1
    assert rows[0]["tf_overlap"] == "not_available"
    assert rows[0]["tf_names"] == "not_available"
    assert rows[0]["target_genes"] == "a,b"
    assert rows[0]["target_functions"] == "a:enzima A; b"


def test_selftest():
    assert "OK" in SR.selftest()
