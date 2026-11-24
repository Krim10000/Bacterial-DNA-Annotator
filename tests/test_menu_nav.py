"""Pruebas puras del menú lateral y del alcance de rutas metabólicas.

No arrancan Streamlit: importan ``app`` en modo bare (sin runtime) y solo
ejercitan las funciones de agrupación y resolución de estado.
"""
from __future__ import annotations

import pytest

import app


# --------------------------------------------------------------------------
# Bloques del menú y claves de modo
# --------------------------------------------------------------------------
def test_menu_has_five_blocks_in_order():
    assert [g["title"] for g in app.MENU_GROUPS] == [
        "Anotar", "Explorar genoma", "Catálogo", "Reportes", "Créditos",
    ]
    assert len({g["id"] for g in app.MENU_GROUPS}) == len(app.MENU_GROUPS)


def test_menu_covers_every_mode_once():
    modes = app.flat_menu_modes()
    assert len(modes) == len(set(modes))
    assert set(modes) == set(app.MODE_LABELS)
    for mode in modes:
        assert app.menu_group_for(mode) is not None


def test_menu_group_for_each_mode():
    assert app.menu_group_for("annotate") == "anotar"
    assert app.menu_group_for("browser") == "explorar"
    assert app.menu_group_for("tree") == "explorar"
    assert app.menu_group_for("pathways") == "catalogo"
    assert app.menu_group_for("modules") == "catalogo"
    assert app.menu_group_for("reports") == "reportes"
    assert app.menu_group_for("hypotheses") == "reportes"
    assert app.menu_group_for("credits") == "creditos"
    assert app.menu_group_for("no-existe") is None


def test_catalog_block_says_it_does_not_change_species():
    catalog = next(g for g in app.MENU_GROUPS if g["id"] == "catalogo")
    assert "especie activa" in catalog["caption"].lower()
    assert "no cambian" in catalog["caption"].lower()


# --------------------------------------------------------------------------
# Enlaces ?modo= (compatibilidad) y estado mode_sel
# --------------------------------------------------------------------------
@pytest.mark.parametrize("query,expected", [
    ("navegador", "browser"),
    ("browser", "browser"),
    ("rutas", "pathways"),
    ("hipotesis", "hypotheses"),
    ("hipótesis", "hypotheses"),
    ("hypotheses", "hypotheses"),
    ("creditos", "credits"),
    ("créditos", "credits"),
    ("reportes", "reports"),
    ("descargas", "reports"),
    ("modulos", "modules"),
    ("módulos", "modules"),
    ("modules", "modules"),
    ("anotar", "annotate"),
    ("arbol", "tree"),
    ("árbol", "tree"),
    ("tree", "tree"),
])
def test_mode_from_query_keeps_legacy_urls(query, expected):
    assert app.mode_from_query(query) == expected


def test_mode_from_query_rejects_empty_and_unknown():
    assert app.mode_from_query(None) is None
    assert app.mode_from_query("") is None
    assert app.mode_from_query("no-existe") is None


def test_canonical_query_values_roundtrip():
    assert set(app.MODE_QUERY_VALUES) == set(app.MODE_LABELS)
    for mode_key, value in app.MODE_QUERY_VALUES.items():
        assert app.mode_from_query(value) == mode_key


def test_mode_key_from_label_roundtrip():
    labels = {key: f"etiqueta-{key}" for key in app.MODE_LABELS}
    for key in app.MODE_LABELS:
        assert app.mode_key_from_label(labels[key], labels) == key
    assert app.mode_key_from_label(None, labels) is None
    assert app.mode_key_from_label("", labels) is None
    assert app.mode_key_from_label("otra", labels) is None


# --------------------------------------------------------------------------
# Especie activa, demo y módulos de catálogo
# --------------------------------------------------------------------------
def test_resolve_active_genome_accepts_demo_and_falls_back():
    options = ["demo", "ecoli_k12", "GCF_000006985.1_ASM698v1"]
    assert app.resolve_active_genome("demo", options) == "demo"
    assert app.resolve_active_genome(options[2], options) == options[2]
    assert app.resolve_active_genome(None, options) == "demo"
    assert app.resolve_active_genome("no-anotado", options) == "demo"
    assert app.resolve_active_genome("demo", []) == ""


def test_active_species_label_is_explicit():
    assert app.active_species_label("demo") == "Especie activa: demo"
    assert app.active_species_label("GCF_1") == "Especie activa: GCF_1"
    assert app.active_species_label("") == "Especie activa: sin genoma"
    assert app.active_species_label(None) == "Especie activa: sin genoma"


def test_species_scope_note_demo_can_change_genome():
    note = app.species_scope_note("demo")
    assert "demo" in note
    assert "cambiarlo" in note.lower()
    assert "bloquea" in note.lower()


def test_catalog_species_note_is_explicit():
    note = app.CATALOG_SPECIES_NOTE.lower()
    assert "catálogo" in note
    assert "400" in note
    assert "no cambian la especie activa" in note


# --------------------------------------------------------------------------
# Rutas metabólicas: especie activa | catálogo de 400, sin mezclar
# --------------------------------------------------------------------------
def test_pathway_scope_order_and_labels():
    assert app.PATHWAY_SCOPE_ORDER == (
        app.PATHWAY_SCOPE_CATALOG, app.PATHWAY_SCOPE_SPECIES)
    options = app.pathway_scope_options("demo")
    assert list(options) == list(app.PATHWAY_SCOPE_ORDER)
    assert "400" in options[app.PATHWAY_SCOPE_CATALOG]
    assert "demo" in options[app.PATHWAY_SCOPE_SPECIES]


def test_normalize_pathway_scope_defaults_to_catalog():
    assert app.normalize_pathway_scope(None) == app.PATHWAY_SCOPE_CATALOG
    assert app.normalize_pathway_scope("basura") == app.PATHWAY_SCOPE_CATALOG
    assert app.normalize_pathway_scope(app.PATHWAY_SCOPE_SPECIES) == (
        app.PATHWAY_SCOPE_SPECIES)


def test_pathway_scope_caption_does_not_mix_scopes():
    species = app.pathway_scope_caption(app.PATHWAY_SCOPE_SPECIES, "demo", 310)
    catalog = app.pathway_scope_caption(app.PATHWAY_SCOPE_CATALOG, "demo", 310)

    assert "especie activa" in species.lower()
    assert "demo" in species
    assert "400" not in species

    assert "catálogo" in catalog.lower()
    assert "400" in catalog
    assert "no usa la especie activa" in catalog.lower()
    assert "demo" not in catalog


def render_catalog() -> None:
    pass


def render_kegg() -> None:
    pass


def render_curated() -> None:
    pass


def render_origin() -> None:
    pass


def render_suspects_tab() -> None:
    pass


def render_gaps() -> None:
    pass


def render_completeness() -> None:
    pass


def render_taxonomy() -> None:
    pass


@pytest.fixture()
def pathway_items():
    return [
        ("Catálogo", render_catalog),
        ("KEGG", render_kegg),
        ("Rutas curadas", render_curated),
        ("Origen", render_origin),
        ("Sospechas metabólicas", render_suspects_tab),
        ("Brechas", render_gaps),
        ("Completitud", render_completeness),
        ("Taxonomía", render_taxonomy),
    ]


def test_pathway_species_tabs_exist_in_app():
    assert app.PATHWAY_SPECIES_TABS == {
        "render_catalog", "render_kegg", "render_curated", "render_origin",
        "render_suspects_tab",
    }
    for name in app.PATHWAY_SPECIES_TABS:
        assert callable(getattr(app, name))


def test_pathway_scopes_partition_all_tabs(pathway_items):
    species = app.pathway_scope_items(pathway_items, app.PATHWAY_SCOPE_SPECIES)
    catalog = app.pathway_scope_items(pathway_items, app.PATHWAY_SCOPE_CATALOG)

    species_labels = [label for label, _ in species]
    catalog_labels = [label for label, _ in catalog]
    all_labels = [label for label, _ in pathway_items]

    assert species_labels == [
        "Catálogo", "KEGG", "Rutas curadas", "Origen",
        "Sospechas metabólicas",
    ]
    assert catalog_labels == ["Brechas", "Completitud", "Taxonomía"]
    assert set(species_labels).isdisjoint(catalog_labels)
    assert set(species_labels) | set(catalog_labels) == set(all_labels)


def test_pathway_unknown_scope_behaves_like_catalog(pathway_items):
    unknown = app.pathway_scope_items(pathway_items, "no-existe")
    catalog = app.pathway_scope_items(
        pathway_items, app.PATHWAY_SCOPE_CATALOG)
    assert [label for label, _ in unknown] == [label for label, _ in catalog]


def test_pathway_scope_groups_drop_empty_groups(pathway_items):
    groups = [
        ("Genoma", pathway_items[:4]),
        ("Brechas", pathway_items[5:]),
        ("Sospechas", [pathway_items[4]]),
    ]

    species = app.pathway_scope_groups(groups, app.PATHWAY_SCOPE_SPECIES)
    assert [title for title, _ in species] == ["Genoma", "Sospechas"]
    assert [label for label, _ in species[0][1]] == [
        "Catálogo", "KEGG", "Rutas curadas", "Origen",
    ]

    catalog = app.pathway_scope_groups(groups, app.PATHWAY_SCOPE_CATALOG)
    assert [title for title, _ in catalog] == ["Brechas"]
    assert [label for label, _ in catalog[0][1]] == [
        "Brechas", "Completitud", "Taxonomía",
    ]
