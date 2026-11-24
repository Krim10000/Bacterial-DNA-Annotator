"""UX del navegador de genoma.

Cubre las piezas que se pueden probar sin levantar Streamlit:

- geometría pura de la vista circular (ángulos, ticks, ventana que cruza el
  origen y ventana que ocupa todo el contig);
- selección pura de ventana para genes, operones, promotores, horquillas,
  complejos, TSS/TTS y sospechas;
- tabla y contador de «elementos en pantalla», incluidos DataFrames vacíos.

Los SVG se validan con ``ElementTree``: si el dibujo no es XML válido, estos
tests fallan antes de que Streamlit lo incruste.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from browser import (
    ELEMENT_COLUMNS,
    WINDOW_SIZES,
    all_elements_on_screen,
    circular_geometry,
    circular_svg,
    contig_lengths,
    element_counts,
    element_span,
    elements_counter_text,
    parse_sel,
    selection_window,
    svg,
    viewport_for_span,
)


# ------------------------------------------------------------------ datos
def _tables() -> dict:
    """Genoma sintético con dos contigs y todas las pistas del navegador."""
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "gene_id": "g1", "locus_tag": "b0001", "gene_name": "thrA",
         "product": "aspartokinase", "ec": "2.7.2.4"},
        {"contig": "c1", "start": 500, "end": 800, "strand": "-1",
         "gene_id": "g2", "locus_tag": "b0002", "gene_name": "a<b>c",
         "product": "producto con <script>alert(1)</script>", "ec": ""},
        {"contig": "c2", "start": 1_000, "end": 1_200, "strand": "1",
         "gene_id": "g3", "locus_tag": "b0100", "gene_name": "otro",
         "product": "otro producto", "ec": ""},
    ])
    operons = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 800, "operon_id": "op1",
         "label": "thr", "n_genes": 2, "span_bp": 700, "confidence": "high",
         "strand": "1", "genes": "b0001,b0002"},
        {"contig": "c2", "start": 1_000, "end": 1_200, "operon_id": "op2",
         "label": "otro", "n_genes": 1, "span_bp": 200, "confidence": "low",
         "strand": "1", "genes": "b0100"},
    ])
    promoters = pd.DataFrame([
        {"operon_id": "op1", "label": "thr", "sigma": "70", "score": "0.9",
         "box35": "TTGACA", "box10": "TATAAT", "spacer": "17"},
    ])
    hairpins = pd.DataFrame([
        {"contig": "c1", "start": 820, "end": 860, "class": "A_fuerte",
         "dG": "-12.3", "stem": "8", "gc": "0.6", "polyU": "5"},
    ])
    complexes = pd.DataFrame([
        {"contig": "c1", "start": 200, "end": 600, "complex_id": "cx1",
         "name": "Complejo de prueba", "status": "complex_complete",
         "confidence": "high", "subunits": "b0001;b0002", "evidence": "vecindad"},
    ])
    return {"genes": genes, "operons": operons, "promoters": promoters,
            "hairpins": hairpins, "complexes": complexes}


def _tss() -> pd.DataFrame:
    return pd.DataFrame([
        {"pos": 150, "strand": "+", "source": "curado", "detail": "sigma70"},
        {"pos": 6_000, "strand": "+", "source": "alto rendimiento", "detail": "HT"},
    ])


def _tts() -> pd.DataFrame:
    return pd.DataFrame([
        {"pos": 950, "strand": "+", "source": "curado", "detail": "terminador"},
    ])


def _suspects() -> pd.DataFrame:
    return pd.DataFrame([
        {"contig": "c1", "start": 300, "end": 700, "n_genes": 2,
         "products": "enzima A; enzima B", "reactions": "R1;R2"},
    ])


# ------------------------------------------------------ geometría circular
def test_circular_geometry_fraction_angles_and_ticks():
    geo = circular_geometry(1_000, 250, 500)
    assert geo["contig_len"] == 1_000
    assert geo["fraction"] == pytest.approx(0.25)
    assert geo["start_angle"] == pytest.approx(90.0)
    assert geo["end_angle"] == pytest.approx(180.0)
    assert geo["sweep"] == pytest.approx(90.0)
    assert geo["full_circle"] is False
    assert geo["window_path"].startswith("M ")
    assert len(geo["ticks"]) == 12
    assert sum(1 for t in geo["ticks"] if t["major"]) == 4
    assert geo["ticks"][0]["pos"] == 0
    assert geo["ticks"][3]["pos"] == 250


def test_circular_geometry_window_wraps_around_origin():
    geo = circular_geometry(1_000, 900, 1_100)
    assert geo["span"] == 200
    assert geo["sweep"] == pytest.approx(72.0)
    assert geo["end_angle"] == pytest.approx(396.0)
    assert geo["wraps"] is True


def test_circular_geometry_full_contig_and_clamping():
    full = circular_geometry(1_000, 0, 1_000)
    assert full["full_circle"] is True
    assert full["fraction"] == 1.0
    assert full["sweep"] == 360.0
    assert full["window_path"] == ""

    bigger = circular_geometry(1_000, 0, 5_000)
    assert bigger["full_circle"] is True

    clamped = circular_geometry(1_000, -10, 10)
    assert clamped["start"] == 0
    assert clamped["span"] == 10
    assert clamped["fraction"] == pytest.approx(0.01)


def test_circular_svg_is_valid_xml_and_explains_the_portion():
    out = circular_svg("c<1>", 10_000, 2_000, 4_000)
    root = ET.fromstring(out)          # XML inválido -> falla aquí
    roles = {el.get("data-role") for el in root.iter() if el.get("data-role")}
    assert {"contig", "window", "window-start", "window-end",
            "ticks", "legend"} <= roles
    assert "porción del círculo" in out
    assert "c&lt;1&gt;" in out
    assert "2,000" in out and "4,000" in out


def test_circular_svg_full_circle_draws_a_ring():
    root = ET.fromstring(circular_svg("c1", 1_000, 0, 1_000))
    window = [el for el in root.iter() if el.get("data-role") == "window"]
    assert len(window) == 1
    assert window[0].tag.endswith("circle")


# ---------------------------------------------------- selección de ventana
def test_element_span_covers_every_clickable_kind():
    tables = _tables()

    gene = element_span("gene", "b0001", tables, contig="c1")
    assert gene is not None
    assert gene["id"] == "b0001"
    assert (gene["start"], gene["end"]) == (100, 400)

    operon = element_span("operon", "op1", tables)
    assert (operon["start"], operon["end"]) == (100, 800)

    promoter = element_span("promoter", "op1", tables)
    assert promoter is not None
    assert (promoter["start"], promoter["end"]) == (100, 800)

    hairpin = element_span("hairpin", "c1:820-860", tables)
    assert hairpin is not None
    assert (hairpin["start"], hairpin["end"]) == (820, 860)

    complex_ = element_span("complex", "cx1", tables)
    assert complex_ is not None
    assert (complex_["start"], complex_["end"]) == (200, 600)

    tss = element_span("tss", "150", tables, contig="c1")
    assert tss is not None
    assert (tss["start"], tss["end"]) == (150, 150)

    tts = element_span("tts", "950", tables, contig="c1")
    assert tts is not None
    assert (tts["start"], tts["end"]) == (950, 950)

    suspect = element_span("suspect", "300-700", tables,
                           suspects=_suspects(), contig="c1")
    assert suspect is not None
    assert (suspect["start"], suspect["end"]) == (300, 700)


def test_element_span_returns_none_for_unknown_or_bad_coordinates():
    tables = _tables()
    assert element_span("gene", "no-existe", tables) is None
    assert element_span("desconocido", "x", tables) is None
    assert element_span("hairpin", "sin-coordenadas", tables) is None
    empty = {k: pd.DataFrame() for k in tables}
    assert element_span("gene", "b0001", empty) is None


def test_viewport_for_span_centers_contains_and_clamps():
    centered = viewport_for_span(4_000, 4_200, 10_000, pad=500, min_size=1_000)
    assert centered["start"] <= 4_000
    assert centered["end"] >= 4_200
    assert centered["start"] >= 0
    assert centered["end"] <= 10_000
    assert centered["size_label"] in WINDOW_SIZES
    assert centered["contains"] is True

    at_start = viewport_for_span(0, 100, 10_000, pad=0, min_size=5_000)
    assert at_start["start"] == 0
    assert at_start["contains"] is True

    at_end = viewport_for_span(9_900, 10_000, 10_000, pad=0, min_size=5_000)
    assert at_end["end"] == 10_000
    assert at_end["contains"] is True


def test_selection_window_sets_contig_and_zoom():
    tables = _tables()
    inside = selection_window("gene", "b0001", tables, "c1", 0, 1_000,
                              contig_len=5_000)
    assert inside is not None
    assert inside["in_current_window"] is True
    assert inside["view_start"] <= inside["start"]
    assert inside["view_end"] >= inside["end"]
    assert inside["size_label"] in WINDOW_SIZES

    outside = selection_window("gene", "b0100", tables, "c1", 0, 1_000,
                               contig_len=5_000)
    assert outside is not None
    assert outside["contig"] == "c2"
    assert outside["in_current_window"] is False
    assert outside["view_start"] <= outside["start"]
    assert outside["view_end"] >= outside["end"]

    assert selection_window("gene", "no-existe", tables, "c1", 0, 1_000) is None


def test_contig_lengths_from_annotation_and_empty_tables():
    lengths = contig_lengths(_tables())
    assert lengths["c1"] >= 860
    assert lengths["c2"] >= 1_200
    assert contig_lengths({}) == {}
    assert contig_lengths({"genes": pd.DataFrame()}) == {}


# ------------------------------------------ tabla y contador en pantalla
def test_all_elements_on_screen_counts_all_tracks_in_window():
    table = all_elements_on_screen(
        _tables(), "c1", 0, 1_000,
        tss=_tss(), tts=_tts(), suspects=_suspects(),
    )
    assert not table.empty
    assert list(table.columns) == list(ELEMENT_COLUMNS)
    assert len(table) == 9
    counts = element_counts(table)
    assert counts["gene"] == 2
    assert counts["operon"] == 1
    assert counts["promoter"] == 1
    assert counts["hairpin"] == 1
    assert counts["complex"] == 1
    assert counts["tss"] == 1
    assert counts["tts"] == 1
    assert counts["suspect"] == 1
    # La TSS de 6.000 pb no pertenece a la ventana 0–1.000.
    assert set(table.loc[table["kind"] == "tss", "id"]) == {"150"}
    # La tabla está ordenada por coordenada y conserva ids del SVG.
    assert table["inicio"].is_monotonic_increasing
    assert "c1:820-860" in set(table["id"])


def test_all_elements_on_screen_empty_dataframes_do_not_break():
    empty = {k: pd.DataFrame() for k in
             ("genes", "operons", "promoters", "hairpins", "complexes")}
    table = all_elements_on_screen(empty, "c1", 0, 1_000)
    assert table.empty
    assert list(table.columns) == list(ELEMENT_COLUMNS)
    assert element_counts(table) == {}
    assert elements_counter_text(element_counts(table)) == "0 elementos en pantalla"

    also_empty = all_elements_on_screen({}, "c1", 0, 1_000,
                                        tss=pd.DataFrame(), tts=pd.DataFrame())
    assert also_empty.empty
    assert list(also_empty.columns) == list(ELEMENT_COLUMNS)


def test_all_elements_on_screen_skips_rows_without_coordinates():
    tables = {
        "genes": pd.DataFrame([
            {"contig": "c1", "start": "", "end": "", "locus_tag": "sin-coords"},
        ]),
        "operons": pd.DataFrame([{"operon_id": "op-sin-coords"}]),
    }
    table = all_elements_on_screen(tables, "c1", 0, 1_000)
    assert table.empty
    assert list(table.columns) == list(ELEMENT_COLUMNS)


def test_counter_text_singular_plural_and_breakdown():
    text = elements_counter_text({"gene": 2, "operon": 1})
    assert text.startswith("3 elementos en pantalla")
    assert "2 genes" in text
    assert "1 operón" in text

    one = elements_counter_text({"gene": 1})
    assert one.startswith("1 elemento en pantalla")
    assert "1 gen" in one

    assert element_counts(pd.DataFrame({"otra": [1, 2]})) == {}


def test_parse_sel_keeps_colons_inside_hairpin_id():
    assert parse_sel("gene:b0001") == {"kind": "gene", "id": "b0001"}
    assert parse_sel("hairpin:c1:820-860") == {
        "kind": "hairpin", "id": "c1:820-860"}
    assert parse_sel("tss:1234") == {"kind": "tss", "id": "1234"}
    assert parse_sel("") is None
    assert parse_sel("gene") is None
    assert parse_sel("raro:x") is None


# ------------------------------------------------------ SVG interactivo
def test_svg_marks_selected_and_all_elements():
    tables = _tables()
    out = svg(tables, "c1", 0, 1_000, focus="operon:op1",
              highlight_all=True, tss=_tss(), tts=_tts())
    ET.fromstring(out)
    assert 'class="hit all-hit sel-hit" data-kind="operon"' in out
    assert ".hit.sel-hit" in out
    # TSS/TTS también son clicables y comparten la identidad de la tabla.
    assert 'data-kind="tss"' in out
    assert 'data-kind="tts"' in out
    assert 'data-kind="complex"' in out
