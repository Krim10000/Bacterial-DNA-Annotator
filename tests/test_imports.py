"""Pruebas de importación: la app debe poder cargar sus módulos sin Streamlit
en ejecución y sin dependencias opcionales instaladas."""
from __future__ import annotations

import importlib

import pytest

MODULES = [
    "i18n",
    "textos",
    "credits",
    "annotator",
    "requirements_ui",
    "browser",
    "closures",
    "operons",
    "pathways",
    "metabolism",
    "structural_filters",
    "structure_evidence",
]


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    mod = importlib.import_module(name)
    assert mod is not None


def test_credits_name_the_author_and_ai_help():
    import credits

    assert credits.AUTHOR["name"] == "Camilo Rojas Bravo"
    assert credits.AUTHOR["email"] == "Camilo.Rojas@usach.cl"
    assert "DeepSeek" in credits.AI_HELP["tool"]
    assert "2026" in credits.AI_HELP["period"]


def test_i18n_languages_include_spanish_english_chinese():
    import i18n

    assert {"es", "en", "zh"} <= set(i18n.LANGUAGES)
