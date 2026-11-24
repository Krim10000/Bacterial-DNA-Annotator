"""Ayuda in situ del navegador y limpieza de ``?genoma=demo``.

Estas pruebas no levantan Streamlit: importan ``browser`` (que solo define el
módulo) y comprueban textos, funciones puras y que ``render`` conecte ambas
cosas. Así se pueden ejecutar en CI sin servidor ni sesión.
"""
from __future__ import annotations

import inspect

import browser


# ------------------------------------------------------------ textos de ayuda
def test_help_expander_title_is_the_expected_one():
    assert browser.HELP_EXPANDER_TITLE == "❓ ¿Qué estoy viendo?"


def test_help_texts_cover_the_four_blocks():
    assert set(browser.HELP_TEXTS) >= {
        "elements", "hairpins", "circular", "operons",
    }
    for topic in ("elements", "hairpins", "circular", "operons"):
        text = browser.HELP_TEXTS[topic]
        assert isinstance(text, str)
        assert len(text) > 200, topic
        assert "**Límite.**" in text, topic


def test_elements_help_explains_scope_and_table_columns():
    text = browser.ELEMENTS_ON_SCREEN_HELP.lower()
    assert "ventana actual" in text
    assert "no el genoma completo" in text
    assert "tabla" in text
    assert "?sel=tipo:id" in browser.ELEMENTS_ON_SCREEN_HELP
    assert "confianza" in text


def test_hairpin_help_distinguishes_terminator_and_attenuator():
    text = browser.HAIRPIN_HELP.lower()
    assert "terminador" in text
    assert "atenuador" in text
    assert "final" in text
    assert "dentro" in text
    assert "poli-u" in text
    assert "δg" in text


def test_circular_help_contrasts_linear_and_circular_views():
    text = browser.CIRCULAR_HELP.lower()
    assert "ventana lineal" in text
    assert "contig completo" in text
    assert "círculo" in text or "circular" in text
    assert "porción" in text or "porcion" in text
    assert "naranja" in text
    assert "cromosoma" in text


def test_operon_help_explains_operon_and_promoter():
    text = browser.OPERON_HELP.lower()
    assert "operón" in text
    assert "promotor" in text
    assert "misma hebra" in text
    assert "−35" in browser.OPERON_HELP and "−10" in browser.OPERON_HELP
    assert "no es una prueba" in text


def test_all_elements_label_and_legend_are_explicit():
    assert browser.ALL_ELEMENTS_LABEL == "Todos los elementos en pantalla"
    legend = browser.ALL_ELEMENTS_LEGEND
    assert "Todos los elementos en pantalla" in legend
    assert "ventana actual" in legend
    assert "no añade" in legend.lower()


def test_help_expander_ignores_unknown_topic_without_streamlit_runtime():
    # El retorno temprano no toca Streamlit; el resto del helper lo haría.
    assert browser.help_expander("tema-inexistente") is None


def test_render_wires_each_help_block():
    source = inspect.getsource(browser.render)
    for topic in ("elements", "hairpins", "circular"):
        assert f'help_expander("{topic}")' in source, topic
    # Operones y promotores comparten la explicación en sus dos pestañas.
    assert source.count('help_expander("operons")') >= 2
    assert "ALL_ELEMENTS_LABEL" in source
    assert "ALL_ELEMENTS_LEGEND" in source


# ---------------------------------------------------- limpieza de ?genoma=demo
def test_clear_demo_query_removes_demo_and_preserves_the_rest():
    original = {"genoma": "demo", "pos": "1000", "w": "10 kb",
                "sel": "gene:b0001", "modo": "navegador"}
    cleaned = browser.clear_demo_query(original)

    assert cleaned == {"pos": "1000", "w": "10 kb",
                       "sel": "gene:b0001", "modo": "navegador"}
    assert "genoma" not in cleaned
    # Función pura: no muta la entrada y devuelve un mapping nuevo.
    assert original["genoma"] == "demo"
    assert cleaned is not original


def test_clear_demo_query_keeps_a_shareable_project_genome():
    query = {"genoma": "GCF_000020625.1", "pos": "500"}
    assert browser.clear_demo_query(query) == query
    assert browser.clear_demo_query({"genoma": "ecoli_k12"}) == {
        "genoma": "ecoli_k12"}


def test_clear_demo_query_handles_lists_and_other_values():
    assert browser.clear_demo_query({"genoma": ["demo"]}) == {}
    assert browser.clear_demo_query({"genoma": ["GCF_1"]}) == {
        "genoma": ["GCF_1"]}
    assert browser.clear_demo_query({"genoma": ""}) == {"genoma": ""}
    assert browser.clear_demo_query(None) == {}
    assert browser.clear_demo_query({}) == {}


def test_clear_demo_query_accepts_a_query_string():
    assert browser.clear_demo_query("?genoma=demo&pos=10") == "?pos=10"
    assert browser.clear_demo_query("genoma=demo") == ""
    assert browser.clear_demo_query(
        "?genoma=GCF_1&pos=10") == "?genoma=GCF_1&pos=10"
    assert browser.clear_demo_query("") == ""
    assert browser.clear_demo_query("?") == ""


def test_render_clears_the_demo_param_from_the_url():
    render_source = inspect.getsource(browser.render)
    assert "_clear_demo_query_from_url()" in render_source
    helper_source = inspect.getsource(browser._clear_demo_query_from_url)
    assert "clear_demo_query" in helper_source
    assert 'del st.query_params["genoma"]' in helper_source


def test_url_helper_removes_only_the_demo_param(monkeypatch):
    from types import SimpleNamespace

    query = {"genoma": "demo", "modo": "navegador", "pos": "10"}
    monkeypatch.setattr(
        browser, "st", SimpleNamespace(query_params=query))

    assert browser._clear_demo_query_from_url() is True
    assert query == {"modo": "navegador", "pos": "10"}


def test_url_helper_leaves_project_genomes_and_clean_urls_alone(monkeypatch):
    from types import SimpleNamespace

    project = {"genoma": "GCF_000020625.1"}
    monkeypatch.setattr(
        browser, "st", SimpleNamespace(query_params=project))
    assert browser._clear_demo_query_from_url() is False
    assert project == {"genoma": "GCF_000020625.1"}

    clean = {"modo": "navegador"}
    monkeypatch.setattr(browser, "st", SimpleNamespace(query_params=clean))
    assert browser._clear_demo_query_from_url() is False
    assert clean == {"modo": "navegador"}

