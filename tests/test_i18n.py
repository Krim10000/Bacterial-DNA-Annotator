"""Comprueban que las traducciones estén completas y no rompan el formato.

Si añades una cadena nueva a la interfaz, ejecuta
``python tools/extract_ui_strings.py`` y traduce las claves que falten.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import i18n

ROOT = Path(__file__).resolve().parents[1]
LOCALES = ROOT / "webapp" / "locales"


def _json(name: str) -> dict:
    return json.loads((LOCALES / name).read_text(encoding="utf-8"))


@pytest.fixture(autouse=True)
def _reset_language():
    i18n.set_language("es")
    yield
    i18n.set_language("es")


def test_catalogs_cover_every_source_string():
    source = set(_json("_source.json"))
    for lang in ("en", "zh"):
        catalog = _json(f"{lang}.json")
        missing = source - set(catalog)
        assert not missing, f"{lang} sin traducir: {sorted(missing)[:10]}"
        assert not [k for k, v in catalog.items() if not str(v).strip()]


def test_translations_keep_placeholders_and_markup():
    source = _json("_source.json")
    for lang in ("en", "zh"):
        catalog = _json(f"{lang}.json")
        for key in source:
            value = catalog[key]
            # mismos placeholders {x}
            assert sorted(re.findall(r"\{\w+\}", key)) == sorted(
                re.findall(r"\{\w+\}", value)), f"{lang}: placeholders en {key!r}"
            # mismo número de backticks y de negritas
            assert key.count("`") == value.count("`"), f"{lang}: backticks en {key!r}"
            assert key.count("**") == value.count("**"), f"{lang}: negritas en {key!r}"


def test_chinese_is_actually_chinese():
    catalog = _json("zh.json")
    with_cjk = sum(1 for v in catalog.values() if re.search(r"[\u4e00-\u9fff]", v))
    assert with_cjk / len(catalog) > 0.35, "demasiadas entradas sin caracteres chinos"


def test_english_differs_from_spanish():
    catalog = _json("en.json")
    changed = sum(1 for k, v in catalog.items() if k != v)
    assert changed / len(catalog) > 0.9


def test_tr_and_textos_resolve_language():
    i18n.set_language("en")
    assert i18n.tr("Idioma") == "Language"
    assert i18n.tr("cadena que no existe") == "cadena que no existe"
    from textos import T

    assert T["mode_annotate"].startswith("🧬")
    i18n.set_language("es")
    assert i18n.tr("Idioma") == "Idioma"


def test_unknown_language_falls_back_to_spanish():
    i18n.set_language("xx")
    assert i18n.get_language() == "es"
