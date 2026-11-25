"""Pista sigma del navegador y vista Sigma (sin arrancar Streamlit).

Cubre, con datos sintéticos:

- carga defensiva de ``sigma_promoters.tsv`` (archivo ausente, vacío o raro);
- resolución de coordenadas por operón y fallback al ``operon_id``;
- forma estable de la tabla y funciones puras de color/confianza/JSON;
- SVG de la pista sigma (color por tipo, opacidad/trazo por confianza) y su
  integración en «elementos en pantalla», ``element_span`` y la ficha;
- métricas, filtros y tabla de ``webapp/sigma_view.py`` sin runtime de
  Streamlit.
"""
from __future__ import annotations

import inspect
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

import browser
from browser import (
    ELEMENT_COLUMNS,
    SIGMA_PROMOTER_COLUMNS,
    all_elements_on_screen,
    element_counts,
    element_span,
    elements_counter_text,
    format_sigma_json,
    load_sigma_promoters,
    sigma_color,
    sigma_confidence,
    sigma_promoters_path,
    svg,
)

SIGMA_COLUMNS = (
    "operon_id", "operon_label", "first_gene", "first_gene_name", "sigma",
    "sigma_class", "sigma_arch", "sigma_status", "confidence",
    "confidence_label", "evidence", "detected_by", "alternatives",
    "provenance",
)


# ------------------------------------------------------------------ datos
def _write_tsv(path: Path, rows: list[dict], columns: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=list(columns)).to_csv(
        path, sep="\t", index=False)


def _operons() -> list[dict]:
    return [
        {"operon_id": "g:ctg1:100:+", "label": "opA", "contig": "ctg1",
         "start": "100", "end": "400", "strand": "1"},
        {"operon_id": "g:ctg1:900:-", "label": "opB", "contig": "ctg1",
         "start": "900", "end": "1200", "strand": "-1"},
    ]


def _promoters() -> list[dict]:
    return [
        {"operon_id": "g:ctg1:100:+", "operon_label": "opA",
         "first_gene": "g1", "first_gene_name": "thrA", "sigma": "sigma70",
         "sigma_class": "housekeeping", "sigma_arch": "-35/-10",
         "sigma_status": "ambiguous", "confidence": "low",
         "confidence_label": "sin_sospecha", "evidence": "PWM sigma70",
         "detected_by": "promoter_pwm:sigma70;inventory:sigma70",
         "alternatives": '[{"sigma": "sigma38", "pvalue": 0.2}]',
         "provenance": '{"contig": "ctg1", "seed": 17}'},
        {"operon_id": "g:ctg1:900:-", "operon_label": "opB",
         "first_gene": "g2", "first_gene_name": "pyrD", "sigma": "sigma38",
         "sigma_class": "stationary", "sigma_arch": "-35/-10",
         "sigma_status": "assigned", "confidence": "medium",
         "confidence_label": "sospecha_media", "evidence": "PWM sigma38",
         "detected_by": "promoter_pwm:sigma38",
         "alternatives": "[]", "provenance": '{"contig": "ctg1"}'},
        {"operon_id": "g:ctg1:2000:-", "operon_label": "opC",
         "first_gene": "g3", "first_gene_name": "fliA", "sigma": "sigma28",
         "sigma_class": "flagellar", "sigma_arch": "-35/-10",
         "sigma_status": "assigned", "confidence": "high",
         "confidence_label": "sospecha_fuerte", "evidence": "PWM sigma28",
         "detected_by": "promoter_pwm:sigma28",
         "alternatives": "[]", "provenance": '{"contig": "ctg1"}'},
    ]


def _sigma_frame() -> pd.DataFrame:
    """Tabla ya resuelta, con la forma que devuelve ``load_sigma_promoters``."""
    return pd.DataFrame([
        {"operon_id": "s1", "operon_label": "opA", "first_gene_name": "thrA",
         "contig": "ctg1", "strand": "+", "sigma": "sigma70",
         "sigma_class": "housekeeping",
         "sigma_status": "assigned", "confidence": "high",
         "confidence_label": "sospecha_fuerte", "evidence": "PWM sigma70",
         "detected_by": "promoter_pwm", "alternatives": "[]",
         "provenance": '{"seed": 17}', "start": 100, "end": 400,
         "confidence_tier": "high", "sigma_color": sigma_color("sigma70")},
        {"operon_id": "s2", "operon_label": "opB", "first_gene_name": "pyrD",
         "contig": "ctg1", "strand": "−", "sigma": "sigma24",
         "sigma_class": "ecf",
         "sigma_status": "ambiguous", "confidence": "low",
         "confidence_label": "sospecha_baja", "evidence": "PWM sigma24",
         "detected_by": "promoter_pwm", "alternatives": "[]",
         "provenance": "{}", "start": 600, "end": 700,
         "confidence_tier": "low", "sigma_color": sigma_color("sigma24")},
        {"operon_id": "s3", "operon_label": "otro", "first_gene_name": "fliA",
         "contig": "ctg2", "strand": "+", "sigma": "sigma70",
         "sigma_class": "housekeeping",
         "sigma_status": "assigned", "confidence": "medium",
         "confidence_label": "sospecha_media", "evidence": "PWM sigma70",
         "detected_by": "promoter_pwm", "alternatives": "[]",
         "provenance": "{}", "start": 100, "end": 200,
         "confidence_tier": "medium", "sigma_color": sigma_color("sigma70")},
    ])


# ------------------------------------------------------- carga defensiva
def test_missing_sigma_file_returns_stable_empty_frame(tmp_path: Path):
    df = load_sigma_promoters(tmp_path, "sin_datos")
    assert df.empty
    assert set(SIGMA_PROMOTER_COLUMNS) <= set(df.columns)
    assert {"start", "end", "confidence_tier", "sigma_color"} <= set(df.columns)
    assert df.attrs["sigma_available"] is False
    assert "sigma_promoters.tsv" in df.attrs["sigma_warning"]
    assert sigma_promoters_path(tmp_path, "sin_datos").name == "sigma_promoters.tsv"


def test_load_sigma_resolves_operon_and_operon_id_fallback(tmp_path: Path):
    _write_tsv(tmp_path / "regulation" / "g" / "operons.tsv", _operons(),
               ("operon_id", "label", "contig", "start", "end", "strand"))
    _write_tsv(tmp_path / "sigma" / "g" / "sigma_promoters.tsv",
               _promoters(), SIGMA_COLUMNS)

    df = load_sigma_promoters(tmp_path, "g")
    assert len(df) == 3
    assert df.attrs["sigma_available"] is True
    assert df.attrs["sigma_warning"] == ""

    by_id = df.set_index("operon_id")
    # 1) operón encontrado: start/end/contig/hebra salen de operons.tsv
    row = by_id.loc["g:ctg1:900:-"]
    assert (row["start"], row["end"]) == (900, 1200)
    assert row["contig"] == "ctg1"
    assert row["strand"] == "−"
    assert row["confidence_tier"] == "medium"
    assert row["sigma_color"] == sigma_color("sigma38")
    # 2) operón ausente: se usa el start codificado en el operon_id
    fallback = by_id.loc["g:ctg1:2000:-"]
    assert (fallback["start"], fallback["end"]) == (2000, 2000)
    assert fallback["contig"] == "ctg1"
    assert fallback["strand"] == "−"
    assert fallback["confidence_tier"] == "high"
    # 3) los valores de confianza se normalizan sin inventar niveles
    assert by_id.loc["g:ctg1:100:+", "confidence_tier"] == "low"


def test_load_tables_attaches_sigma_pista(tmp_path: Path):
    _write_tsv(tmp_path / "regulation" / "g" / "genes.tsv",
               [{"contig": "ctg1", "start": "1", "end": "50",
                 "strand": "1", "locus_tag": "g1", "gene_name": "gen1"}],
               ("contig", "start", "end", "strand", "locus_tag", "gene_name"))
    _write_tsv(tmp_path / "regulation" / "g" / "operons.tsv", _operons(),
               ("operon_id", "label", "contig", "start", "end", "strand"))
    _write_tsv(tmp_path / "sigma" / "g" / "sigma_promoters.tsv",
               _promoters(), SIGMA_COLUMNS)
    tables = browser.load_tables(tmp_path, "g")
    assert "sigma" in tables
    assert len(tables["sigma"]) == 3
    assert tables["sigma"].attrs["sigma_available"] is True
    # La pista comparte identidad con el resto del navegador.
    assert "g:ctg1:100:+" in set(tables["sigma"]["operon_id"])


def test_sigma_loader_does_not_break_on_empty_or_odd_files(tmp_path: Path):
    path = tmp_path / "sigma" / "g" / "sigma_promoters.tsv"
    path.parent.mkdir(parents=True)

    path.write_text("", encoding="utf-8")
    empty = load_sigma_promoters(tmp_path, "g")
    assert empty.empty and empty.attrs["sigma_available"] is False

    path.write_text("genome\n", encoding="utf-8")
    header_only = load_sigma_promoters(tmp_path, "g")
    assert header_only.empty and header_only.attrs["sigma_available"] is False

    # TSV con columnas mínimas: carga, pero las filas sin coordenadas no se
    # dibujan ni entran en la tabla de pantalla.
    path.write_text("operon_id\tsigma\nsin:coords:1:+\tsigma70\n",
                    encoding="utf-8")
    odd = load_sigma_promoters(tmp_path, "g")
    assert len(odd) == 1
    assert odd.iloc[0]["start"] == 1
    assert odd.iloc[0]["confidence_tier"] == "unknown"
    out = svg({}, "ctg1", 0, 10, tracks={"sigma"}, sigma=odd)
    ET.fromstring(out)  # XML válido aunque la marca quede fuera de la ventana


def test_confidence_color_and_json_helpers_are_defensive():
    assert sigma_confidence({"confidence": "HIGH"}) == "high"
    assert sigma_confidence({"confidence_label": "sospecha_media"}) == "medium"
    assert sigma_confidence({"confidence": "raro"}) == "unknown"
    assert sigma_confidence({"confidence": ""}) == "unknown"

    assert sigma_color("sigma70") == "#177245"
    assert sigma_color("70") == "#177245"
    assert sigma_color("sigma_inventada") == browser.SIGMA_DEFAULT_COLOR

    pretty = format_sigma_json('[{"sigma": "sigma38"}]')
    assert '"sigma": "sigma38"' in pretty and "\n" in pretty
    assert format_sigma_json("no-json") == "no-json"
    assert format_sigma_json({"a": 1}) == '{\n  "a": 1\n}'
    assert format_sigma_json("x" * 50, max_chars=10).endswith("…")


# ------------------------------------------------------------- SVG / tabla
def _sigma_groups(out: str) -> dict[str, ET.Element]:
    root = ET.fromstring(out)
    groups = [el for el in root.iter()
              if el.get("data-kind") == "sigma"]
    return {el.get("data-id"): el for el in groups}


def test_svg_sigma_track_filters_window_and_styles_confidence():
    out = svg({}, "ctg1", 0, 1000, tracks={"sigma"}, sigma=_sigma_frame())
    ET.fromstring(out)
    groups = _sigma_groups(out)
    assert set(groups) == {"s1", "s2"}           # s3 es otro contig
    assert groups["s1"].get("data-sigma") == "sigma70"
    assert groups["s1"].get("data-confidence") == "high"
    assert groups["s2"].get("data-confidence") == "low"

    high_lines = [el for el in groups["s1"].iter() if el.tag.endswith("line")]
    low_lines = [el for el in groups["s2"].iter() if el.tag.endswith("line")]
    assert high_lines and low_lines
    assert float(high_lines[0].get("stroke-opacity")) > float(
        low_lines[0].get("stroke-opacity"))
    assert low_lines[0].get("stroke-dasharray") == "3 2"
    assert high_lines[0].get("stroke") == sigma_color("sigma70")
    assert high_lines[0].get("stroke") != low_lines[0].get("stroke")


def test_svg_sigma_track_respects_track_selector_and_handles_empty():
    data = _sigma_frame()
    without = svg({}, "ctg1", 0, 1000, tracks={"genes"}, sigma=data)
    assert "data-kind=\"sigma\"" not in without
    empty = svg({}, "ctg1", 0, 1000, tracks={"sigma"}, sigma=pd.DataFrame())
    ET.fromstring(empty)
    assert "data-kind=\"sigma\"" not in empty


def test_sigma_rows_enter_screen_elements_and_element_span():
    tables = {"sigma": _sigma_frame()}
    table = all_elements_on_screen(tables, "ctg1", 0, 1000,
                                   sigma=_sigma_frame())
    assert list(table.columns) == list(ELEMENT_COLUMNS)
    sigma_rows = table[table["kind"] == "sigma"]
    assert set(sigma_rows["id"]) == {"s1", "s2"}
    detail = sigma_rows[sigma_rows["id"] == "s1"].iloc[0]["detalle"]
    assert "σ70" in detail and "sospecha_fuerte" in detail and "high" in detail
    assert "2 anotaciones sigma" in elements_counter_text(element_counts(table))

    span = element_span("sigma", "s1", tables, contig="ctg1")
    assert span is not None
    assert (span["start"], span["end"]) == (100, 400)
    assert span["label"].startswith("σ70")
    assert element_span("sigma", "no-existe", tables) is None
    # Una tabla ausente no rompe la búsqueda del elemento.
    assert element_span("sigma", "s1", {}) is None


def test_screen_elements_without_sigma_keep_stable_shape():
    table = all_elements_on_screen({}, "ctg1", 0, 1000,
                                   sigma=pd.DataFrame())
    assert table.empty
    assert list(table.columns) == list(ELEMENT_COLUMNS)


# ------------------------------------------------------- integración GUI
def test_sigma_is_wired_into_render_and_pick_without_streamlit():
    assert browser.ELEMENT_KIND_LABELS["sigma"] == "Anotaciones sigma"
    assert "**Límite.**" in browser.SIGMA_HELP
    render_src = inspect.getsource(browser.render)
    for token in ("load_sigma_promoters", "br_t_sigma", "sigma=sigma_df",
                  "all_elements_on_screen", "sigma_legend_html"):
        assert token in render_src, token
    pick_src = inspect.getsource(browser._render_pick)
    assert "_sigma_detail" in pick_src
    detail_src = inspect.getsource(browser._sigma_detail)
    for token in ("evidence", "detected_by", "alternatives", "provenance"):
        assert token in detail_src, token


def test_sigma_view_summary_filters_and_table_view():
    import sigma_view

    df = _sigma_frame()
    stats = sigma_view.summary(df)
    assert stats["total"] == 3
    assert stats["assigned"] == 2
    assert stats["strong"] == 1
    assert stats["medium"] == 1
    assert stats["low"] == 1
    types = {item["sigma"]: item for item in stats["by_type"]}
    assert types["sigma70"]["total"] == 2
    assert types["sigma70"]["assigned"] == 2
    assert types["sigma24"]["low"] == 1

    only_70 = sigma_view.filter_promoters(df, sigma=["sigma70"])
    assert set(only_70["operon_id"]) == {"s1", "s3"}
    medium_70 = sigma_view.filter_promoters(
        df, sigma="sigma70", confidence_label=["sospecha_media"])
    assert list(medium_70["operon_id"]) == ["s3"]
    by_name = sigma_view.filter_promoters(df, query="pyrD")
    assert list(by_name["operon_id"]) == ["s2"]
    assert sigma_view.filter_promoters(df, query="no-existe").empty
    assert sigma_view.filter_promoters(None).empty

    view = sigma_view.table_view(df)
    assert list(view.columns) == [dst for _, dst in sigma_view.TABLE_COLUMNS]
    assert "σ70" in set(view["σ"])
    assert "sospecha fuerte" in set(view["Confianza"])
    assert view.loc[view["Operón"] == "opB", "Nivel"].iloc[0] == "baja"
    empty_view = sigma_view.table_view(None)
    assert empty_view.empty
    assert list(empty_view.columns) == [
        dst for _, dst in sigma_view.TABLE_COLUMNS]

    render_src = inspect.getsource(sigma_view.render)
    for token in ("load_sigma_promoters", "download_button", "_render_table"):
        assert token in render_src, token
    detail_src = inspect.getsource(sigma_view._detail)
    for token in ("evidence", "detected_by", "alternatives", "provenance"):
        assert token in detail_src, token
