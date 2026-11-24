"""Carriles de hebra y solapamientos de CDS del navegador.

Todo es puro: DataFrames sintéticos y SVG como texto/XML. No se levanta
Streamlit ni se llama a ``render``. Cubre:

- separación de genes en dos carriles (``+`` arriba, ``−`` abajo), con flecha
  en el sentido de lectura 5'→3' en ambos;
- detección de pares de CDS solapados en la misma ventana/contig, con hebra,
  marco (solo si viene en los datos) y locus compartido;
- tabla legible y resumen del expander «Solapamientos en pantalla».
"""
from __future__ import annotations

import inspect
import xml.etree.ElementTree as ET

import pandas as pd
import pytest

import browser
from browser import (
    GENE_LANE_COLUMNS,
    GENE_LANE_Y,
    LOCUS_SHARE_FRACTION,
    OVERLAP_COLUMNS,
    OVERLAP_VIEW_COLUMNS,
    OVERLAPS_TITLE,
    OVERLAP_WARNING,
    cds_overlaps,
    gene_lane_layout,
    gene_frame,
    normalize_frame,
    normalize_strand,
    overlap_summary,
    overlap_table_view,
    svg,
)


# ------------------------------------------------------------------ datos
def _strand_genes() -> pd.DataFrame:
    return pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "b0001", "gene_id": "g1", "gene_name": "thrA"},
        {"contig": "c1", "start": 500, "end": 800, "strand": "-1",
         "locus_tag": "b0002", "gene_id": "g2", "gene_name": "revB"},
        {"contig": "c1", "start": 900, "end": 1_000, "strand": "-",
         "locus_tag": "b0003", "gene_id": "g3", "gene_name": ""},
        {"contig": "c1", "start": 1_200, "end": 1_400, "strand": "+",
         "locus_tag": "b0004", "gene_id": "g4", "gene_name": "fwdD"},
        {"contig": "c1", "start": 1_500, "end": 1_600, "strand": "?",
         "locus_tag": "b0005", "gene_id": "g5", "gene_name": "sinHebra"},
        {"contig": "c2", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "b9001", "gene_id": "g9", "gene_name": "otroContig"},
    ])


def _overlap_genes() -> pd.DataFrame:
    """Cuatro pares solapados: mismo marco, distinto marco y marco unknown."""
    return pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "gA", "frame": "0"},
        {"contig": "c1", "start": 300, "end": 600, "strand": "1",
         "locus_tag": "gB", "frame": "0"},
        {"contig": "c1", "start": 500, "end": 800, "strand": "-1",
         "locus_tag": "gC", "frame": "1"},
        {"contig": "c1", "start": 700, "end": 1_000, "strand": "1",
         "locus_tag": "gD", "frame": "2"},
        {"contig": "c1", "start": 900, "end": 1_100, "strand": "1",
         "locus_tag": "gE"},
    ])


def _pair(table: pd.DataFrame, gen_a: str, gen_b: str) -> pd.Series:
    """Fila del par, sin depender del orden interno de ``gen_a``/``gen_b``."""
    mask = (((table["gen_a"] == gen_a) & (table["gen_b"] == gen_b))
            | ((table["gen_a"] == gen_b) & (table["gen_b"] == gen_a)))
    found = table[mask]
    assert len(found) == 1, (gen_a, gen_b, list(zip(table["gen_a"], table["gen_b"])))
    return found.iloc[0]


def _lane_group(root: ET.Element, lane: str) -> ET.Element:
    groups = [el for el in root.iter()
              if el.get("data-role") == "gene-lane" and el.get("data-lane") == lane]
    assert len(groups) == 1
    return groups[0]


def _first_gene(group: ET.Element) -> ET.Element:
    genes = [el for el in group.iter() if el.get("data-kind") == "gene"]
    assert genes
    return genes[0]


def _polygon_xs(gene: ET.Element) -> list[float]:
    polygons = [el for el in gene.iter() if el.tag.endswith("polygon")]
    assert polygons, "el gen debe llevar una flecha (polygon)"
    points = polygons[0].get("points", "").split()
    return [float(point.split(",")[0]) for point in points]


# ------------------------------------------------------- normalización
@pytest.mark.parametrize("value,expected", [
    ("1", "+"), ("+", "+"), ("+1", "+"), ("forward", "+"),
    ("-1", "−"), ("-", "−"), ("−1", "−"), ("reverse", "−"),
    ("", ""), (None, ""), ("raro", ""),
])
def test_normalize_strand_is_explicit(value, expected):
    assert normalize_strand(value) == expected


def test_normalize_frame_only_accepts_declared_frames():
    assert normalize_frame("0") == "0"
    assert normalize_frame(1) == "1"
    assert normalize_frame("2.0") == "2"
    assert normalize_frame("3") == "unknown"
    assert normalize_frame("-1") == "unknown"
    assert normalize_frame(None) == "unknown"
    assert normalize_frame("") == "unknown"
    # ``phase`` es el nombre GFF3 habitual: también cuenta como dato.
    assert gene_frame({"phase": "2"}) == "2"
    assert gene_frame({"frame": "no-vale", "phase": "1"}) == "1"
    assert gene_frame({"start": 10, "end": 20}) == "unknown"


# ------------------------------------------------------- layout de carriles
def test_gene_lane_layout_puts_plus_above_minus():
    layout = gene_lane_layout(_strand_genes(), "c1", 0, 2_000)
    assert list(layout.columns) == list(GENE_LANE_COLUMNS)
    assert GENE_LANE_Y["+"] < GENE_LANE_Y["−"]

    by_id = layout.set_index("id")
    assert by_id.loc["b0001", "lane"] == "+"
    assert by_id.loc["b0001", "direccion"] == "→"
    assert by_id.loc["b0001", "y"] == GENE_LANE_Y["+"]
    assert by_id.loc["b0002", "lane"] == "−"
    assert by_id.loc["b0002", "direccion"] == "←"
    assert by_id.loc["b0002", "y"] == GENE_LANE_Y["−"]
    assert by_id.loc["b0004", "lane"] == "+"
    assert by_id.loc["b0003", "hebra"] == "−"

    # primero todos los + y después todos los −: el dibujo es estable
    lanes = list(layout["lane"])
    assert lanes == sorted(lanes, key=lambda lane: 0 if lane == "+" else 1)


def test_gene_lane_layout_unknown_strand_is_honest_and_kept_in_plus_lane():
    layout = gene_lane_layout(_strand_genes(), "c1", 0, 2_000)
    row = layout[layout["id"] == "b0005"].iloc[0]
    assert row["lane"] == "+"
    assert row["hebra"] == "?"          # no se inventa la orientación
    assert row["direccion"] == "→"


def test_gene_lane_layout_filters_contig_and_window():
    layout = gene_lane_layout(_strand_genes(), "c1", 1_001, 1_300)
    assert set(layout["id"]) == {"b0004"}
    other = gene_lane_layout(_strand_genes(), "c2", 0, 2_000)
    assert set(other["id"]) == {"b9001"}
    assert gene_lane_layout(pd.DataFrame(), "c1", 0, 2_000).empty
    sin_coords = pd.DataFrame([{"contig": "c1", "locus_tag": "x"}])
    assert gene_lane_layout(sin_coords, "c1", 0, 2_000).empty


# ------------------------------------------------- click y flecha en el SVG
def test_svg_draws_two_lanes_with_five_to_three_arrows():
    out = svg({"genes": _strand_genes()}, "c1", 0, 2_000)
    root = ET.fromstring(out)          # XML inválido -> falla aquí

    plus = _lane_group(root, "+")
    minus = _lane_group(root, "−")
    assert plus.get("data-reading") == "left-to-right"
    assert minus.get("data-reading") == "right-to-left"
    assert plus.get("data-direction") == "→"
    assert minus.get("data-direction") == "←"

    plus_ids = {el.get("data-id") for el in plus.iter()
                if el.get("data-kind") == "gene"}
    minus_ids = {el.get("data-id") for el in minus.iter()
                 if el.get("data-kind") == "gene"}
    assert plus_ids == {"b0001", "b0004", "b0005"}
    assert minus_ids == {"b0002", "b0003"}

    plus_gene = _first_gene(plus)      # b0001
    minus_gene = _first_gene(minus)    # b0002
    plus_xs = _polygon_xs(plus_gene)
    minus_xs = _polygon_xs(minus_gene)
    # la punta va al final del rectángulo en + y al principio en −
    assert plus_xs[-1] > plus_xs[0]
    assert minus_xs[-1] < minus_xs[0]

    # el carril + está más arriba que el −
    def rect_y(gene: ET.Element) -> float:
        rects = [el for el in gene.iter() if el.tag.endswith("rect")]
        return float(rects[0].get("y"))

    assert rect_y(plus_gene) < rect_y(minus_gene)


def test_svg_lane_labels_and_legend_are_explicit():
    out = svg({"genes": _strand_genes()}, "c1", 0, 2_000)
    root = ET.fromstring(out)
    labels = [el.text for el in root.iter()
              if el.get("data-role") == "lane-label"]
    assert labels == ["GENES + · 5'→3' →", "GENES − · 5'→3' ←"]
    assert "5'→3'" in " ".join(labels)


def test_svg_keeps_click_id_and_focus_highlight_in_both_lanes():
    out = svg({"genes": _strand_genes()}, "c1", 0, 2_000,
              focus="gene:b0002")
    ET.fromstring(out)
    assert 'class="hit sel-hit" data-kind="gene" data-id="b0002"' in out
    plus = [el for el in ET.fromstring(out).iter()
            if el.get("data-role") == "gene-lane" and el.get("data-lane") == "+"]
    minus = [el for el in ET.fromstring(out).iter()
             if el.get("data-role") == "gene-lane" and el.get("data-lane") == "−"]
    assert plus and minus


# --------------------------------------------------------- solapamientos
def test_cds_overlaps_same_strand_and_same_frame():
    overlaps = cds_overlaps(_overlap_genes(), "c1", 0, 2_000)
    assert list(overlaps.columns) == list(OVERLAP_COLUMNS)
    row = _pair(overlaps, "gA", "gB")
    assert (row["inicio"], row["fin"]) == (300, 400)
    assert row["bp_solape"] == 101        # 300..400 inclusivo
    assert row["hebras"] == "misma (+/+)"
    assert bool(row["misma_hebra"]) is True
    assert (row["marco_a"], row["marco_b"]) == ("0", "0")
    assert bool(row["mismo_marco"]) is True
    assert row["relacion_marco"] == "mismo"
    assert bool(row["locus_compartido"]) is False
    assert row["cobertura_corta"] == pytest.approx(101 / 301, rel=1e-3)


def test_cds_overlaps_opposite_strands_and_different_frame():
    overlaps = cds_overlaps(_overlap_genes(), "c1", 0, 2_000)
    row = _pair(overlaps, "gB", "gC")
    assert row["hebras"] == "opuestas (+/−)"
    assert bool(row["misma_hebra"]) is False
    assert (row["marco_a"], row["marco_b"]) == ("0", "1")
    assert bool(row["mismo_marco"]) is False
    assert row["relacion_marco"] == "distinto"


def test_cds_overlaps_never_invents_frame():
    overlaps = cds_overlaps(_overlap_genes(), "c1", 0, 2_000)
    row = _pair(overlaps, "gD", "gE")
    assert bool(row["misma_hebra"]) is True
    assert row["marco_a"] == "2"
    assert row["marco_b"] == "unknown"
    assert row["mismo_marco"] is None
    assert row["relacion_marco"] == "unknown"

    # sin ninguna columna de marco, todos los pares quedan unknown
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 300, "strand": "1",
         "locus_tag": "x"},
        {"contig": "c1", "start": 200, "end": 400, "strand": "1",
         "locus_tag": "y"},
    ])
    no_frame = cds_overlaps(genes, "c1", 0, 1_000)
    assert set(no_frame["marco_a"]) == {"unknown"}
    assert set(no_frame["marco_b"]) == {"unknown"}
    assert all(value is None for value in no_frame["mismo_marco"])
    assert set(no_frame["relacion_marco"]) == {"unknown"}


def test_cds_overlaps_detects_shared_locus_by_containment_and_identifier():
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 1_000, "strand": "1",
         "locus_tag": "largo"},
        {"contig": "c1", "start": 200, "end": 900, "strand": "-1",
         "locus_tag": "corto"},
    ])
    nested = cds_overlaps(genes, "c1", 0, 2_000)
    assert len(nested) == 1
    assert bool(nested.iloc[0]["locus_compartido"]) is True
    assert nested.iloc[0]["cobertura_corta"] >= LOCUS_SHARE_FRACTION

    same_tag = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "mismo"},
        {"contig": "c1", "start": 300, "end": 600, "strand": "-1",
         "locus_tag": "mismo"},
    ])
    by_id = cds_overlaps(same_tag, "c1", 0, 2_000)
    assert bool(by_id.iloc[0]["locus_compartido"]) is True
    assert by_id.iloc[0]["cobertura_corta"] < LOCUS_SHARE_FRACTION


def test_cds_overlaps_requires_same_contig_and_window_intersection():
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "a"},
        {"contig": "c1", "start": 300, "end": 600, "strand": "1",
         "locus_tag": "b"},
        {"contig": "c2", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "c"},
    ])
    assert len(cds_overlaps(genes, "c1", 0, 2_000)) == 1
    assert cds_overlaps(genes, "c2", 0, 2_000).empty
    # ventana que solo toca uno de los dos genes: no hay par en pantalla
    assert cds_overlaps(genes, "c1", 450, 700).empty
    # ventana que toca a ambos: el par sí aparece
    assert len(cds_overlaps(genes, "c1", 350, 380)) == 1


def test_cds_overlaps_treats_coordinates_as_inclusive():
    touching = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 200, "strand": "1",
         "locus_tag": "a"},
        {"contig": "c1", "start": 200, "end": 300, "strand": "-1",
         "locus_tag": "b"},
    ])
    overlaps = cds_overlaps(touching, "c1", 0, 1_000)
    assert len(overlaps) == 1
    assert overlaps.iloc[0]["bp_solape"] == 1

    adjacent = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 200, "strand": "1",
         "locus_tag": "a"},
        {"contig": "c1", "start": 201, "end": 300, "strand": "-1",
         "locus_tag": "b"},
    ])
    assert cds_overlaps(adjacent, "c1", 0, 1_000).empty


def test_cds_overlaps_counts_cds_and_ignores_other_gene_types():
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "locus_tag": "cds1", "gene_type": "CDS"},
        {"contig": "c1", "start": 300, "end": 500, "strand": "1",
         "locus_tag": "trna", "gene_type": "tRNA"},
        {"contig": "c1", "start": 350, "end": 600, "strand": "-1",
         "locus_tag": "cds2", "gene_type": "CDS"},
        {"contig": "c1", "start": 400, "end": 700, "strand": "1",
         "locus_tag": "sin_tipo", "gene_type": ""},
    ])
    overlaps = cds_overlaps(genes, "c1", 0, 2_000)
    ids = set(overlaps["gen_a"]) | set(overlaps["gen_b"])
    assert "trna" not in ids
    assert {"cds1", "cds2", "sin_tipo"} <= ids
    # el tipo vacío se conserva como desconocido, no se descarta
    assert _pair(overlaps, "cds1", "sin_tipo")["bp_solape"] == 1


def test_cds_overlaps_empty_inputs_keep_stable_columns():
    for value in (None, pd.DataFrame(),
                  pd.DataFrame([{"contig": "c1", "locus_tag": "sin"}]),):
        table = cds_overlaps(value, "c1", 0, 1_000)
        assert table.empty
        assert list(table.columns) == list(OVERLAP_COLUMNS)


# ------------------------------------- tabla legible y resumen del expander
def test_overlap_table_view_is_human_readable_and_shows_unknown():
    overlaps = cds_overlaps(_overlap_genes(), "c1", 0, 2_000)
    view = overlap_table_view(overlaps)
    assert list(view.columns) == list(OVERLAP_VIEW_COLUMNS)
    assert len(view) == len(overlaps)

    known = view[view["Par"] == "gA ↔ gB"].iloc[0]
    assert "0 / 0" in known["Marco"]
    assert "mismo marco" in known["Marco"]
    assert known["Locus compartido"] == "no"

    unknown = view[view["Par"] == "gD ↔ gE"].iloc[0]
    assert "unknown" in unknown["Marco"]
    assert "2 / unknown" in unknown["Marco"]

    assert overlap_table_view(None).empty
    assert list(overlap_table_view(None).columns) == list(OVERLAP_VIEW_COLUMNS)


def test_overlap_summary_counts_strands_frames_and_locus():
    overlaps = cds_overlaps(_overlap_genes(), "c1", 0, 2_000)
    text = overlap_summary(overlaps)
    assert "4 solapamientos" in text
    assert "misma hebra: 2" in text
    assert "hebras opuestas: 2" in text
    assert "mismo marco: 1" in text
    assert "distinto marco: 2" in text
    assert "marco unknown: 1" in text
    assert "locus compartido: 0" in text
    assert overlap_summary(None) == "0 solapamientos entre CDS en pantalla."


def test_overlaps_expander_copy_and_warning():
    assert OVERLAPS_TITLE == "Solapamientos en pantalla"
    lower = OVERLAP_WARNING.lower()
    assert "pyrodigal" in lower
    assert "ambos" in lower
    assert "unknown" in lower
    # El expander existe en ``render``; se comprueba sin ejecutar Streamlit.
    source = inspect.getsource(browser.render)
    assert "OVERLAPS_TITLE" in source
    assert "st.expander" in source
    assert "overlap_table_view" in source
    assert "OVERLAP_WARNING" in source


def test_svg_groups_plus_and_minus_sections():
    import pandas as pd
    from webapp.browser import svg
    genes = pd.DataFrame([
        {"id": "g1", "locus_tag": "g1", "gene_name": "a", "product": "p",
         "contig": "c", "start": 100, "end": 300, "strand": "+",
         "gene_type": "CDS"},
        {"id": "g2", "locus_tag": "g2", "gene_name": "b", "product": "q",
         "contig": "c", "start": 200, "end": 400, "strand": "-",
         "gene_type": "CDS"},
    ])
    operons = pd.DataFrame([
        {"operon_id": "o1", "label": "o1", "contig": "c", "start": 100,
         "end": 300, "strand": "+"},
    ])
    out = svg({"genes": genes, "operons": operons}, "c", 0, 1000,
              tracks={"genes", "operones"})
    assert "HEBRA +" in out
    assert "HEBRA −" in out
    # el operón + se dibuja en la banda + (Y=132, rect en y=124),
    # no en la banda − (rect en y=372).
    assert 'y="124"' in out
    assert 'y="372"' not in out
