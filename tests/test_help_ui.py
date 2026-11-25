"""Pruebas de la ayuda contextual in situ (``webapp/help_ui.py``).

Cubren el contrato de temas, que cada idioma tenga contenido real (sin campos
vacíos), que el español explique los conceptos propios de cada tema y que
``render_help`` funcione dentro de un ``Streamlit AppTest`` real.

No se prueban datos del pipeline: todo el contenido es estático y vive en el
módulo de ayuda.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

import help_ui

WEBAPP = Path(__file__).resolve().parents[1] / "webapp"

EXPECTED_TOPICS = [
    "transporte",
    "membrana",
    "operones",
    "brechas",
    "complejos",
    "regulacion",
    "hipotesis",
    "taxonomia",
]

LANG_CODES = ["es", "en", "zh"]
REQUIRED_FIELDS = {
    "title",
    "summary",
    "details_heading",
    "details",
    "limits_heading",
    "limits",
}

#: Conceptos que la explicación en español debe nombrar para ser útil in situ.
TOPIC_KEYWORDS = {
    "transporte": [
        "transportador", "high", "medium", "low",
        "no es una función confirmada",
    ],
    "membrana": [
        "señal péptido", "tm", "lipobox", "tmhmm", "signalp",
    ],
    "operones": [
        "distancia", "cotranscripción", "tss", "tts",
    ],
    "brechas": [
        "brecha", "por qué no se cierra", "validated", "anomalous",
    ],
    "complejos": [
        "vecindad", "modelo", "transferencia", "amarillo no confirma",
    ],
    "regulacion": [
        "reglas", "rule_based", "experimental",
    ],
    "hipotesis": [
        "f1", "f2", "f3", "f4", "f5",
        "supported_hypothesis", "conflict",
    ],
    "taxonomia": [
        "precalculada", "ani", "no consulta",
    ],
}


def _joined(topic: str, lang: str = "es") -> str:
    """Todo el texto de un bloque, en minúsculas, para buscar conceptos."""
    block = help_ui.HELP[topic][lang]
    parts = [
        block["title"],
        block["summary"],
        block["details_heading"],
        block["limits_heading"],
        *block["details"],
        *block["limits"],
    ]
    return "\n".join(parts).lower()


def _apptest(script: str) -> AppTest:
    at = AppTest.from_string(script, default_timeout=60).run()
    return at


def _render_script(*calls: str) -> str:
    """Script de AppTest que añade ``webapp`` al path e importa el módulo."""
    prelude = (
        "import sys\n"
        f"sys.path.insert(0, {str(WEBAPP)!r})\n"
        "from help_ui import render_help\n"
    )
    return prelude + "\n".join(calls) + "\n"


# ---------------------------------------------------------------------------
# Contrato de temas y campos
# ---------------------------------------------------------------------------
def test_topics_are_exactly_the_expected_ones():
    assert list(help_ui.TOPICS) == EXPECTED_TOPICS
    assert set(help_ui.HELP) == set(EXPECTED_TOPICS)


def test_expander_label_is_the_requested_one():
    assert help_ui.EXPANDER_LABEL == "❓ ¿Qué significa esto?"


def test_render_help_api_signature():
    sig = inspect.signature(help_ui.render_help)
    assert list(sig.parameters) == ["topic", "key_prefix", "expanded"]
    assert sig.parameters["topic"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert (sig.parameters["key_prefix"].kind
            is inspect.Parameter.KEYWORD_ONLY)
    assert sig.parameters["expanded"].kind is inspect.Parameter.KEYWORD_ONLY
    assert sig.parameters["key_prefix"].default == ""
    assert sig.parameters["expanded"].default is False
    assert sig.return_annotation in (None, "None")


def test_languages_are_the_interface_languages():
    assert [code for code, _ in help_ui.LANGUAGES] == LANG_CODES
    assert [label for _, label in help_ui.LANGUAGES] == [
        "Español", "English", "中文",
    ]


@pytest.mark.parametrize("topic", EXPECTED_TOPICS)
def test_every_topic_has_all_languages_and_fields(topic):
    entry = help_ui.HELP[topic]
    assert set(entry) == set(LANG_CODES)
    for lang in LANG_CODES:
        block = entry[lang]
        assert set(block) == REQUIRED_FIELDS, f"{topic}/{lang}"


@pytest.mark.parametrize("topic", EXPECTED_TOPICS)
@pytest.mark.parametrize("lang", LANG_CODES)
def test_no_empty_fields(topic, lang):
    block = help_ui.HELP[topic][lang]
    for field in ("title", "summary", "details_heading", "limits_heading"):
        value = block[field]
        assert isinstance(value, str), f"{topic}/{lang}/{field}"
        assert value.strip(), f"{topic}/{lang}/{field} vacío"
        assert value.strip().lower() not in {"none", "null", "n/a"}

    for field in ("details", "limits"):
        value = block[field]
        assert isinstance(value, list), f"{topic}/{lang}/{field}"
        assert len(value) >= 3, f"{topic}/{lang}/{field} demasiado corto"
        for i, item in enumerate(value):
            assert isinstance(item, str), f"{topic}/{lang}/{field}[{i}]"
            assert item.strip(), f"{topic}/{lang}/{field}[{i}] vacío"


@pytest.mark.parametrize("topic", EXPECTED_TOPICS)
@pytest.mark.parametrize("lang", LANG_CODES)
def test_blocks_are_substantive(topic, lang):
    block = help_ui.HELP[topic][lang]
    # El chino concentra más contenido por carácter; el umbral es menor.
    min_summary = 60 if lang == "zh" else 120
    assert len(block["summary"]) >= min_summary
    assert len(block["details"]) >= 4
    assert len(block["limits"]) >= 3
    # No debe remitir a la documentación como única explicación.
    assert "docs/" not in _joined(topic, lang)


# ---------------------------------------------------------------------------
# Contenido específico de cada tema
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("topic", EXPECTED_TOPICS)
def test_spanish_explains_topic_specific_concepts(topic):
    text = _joined(topic, "es")
    missing = [kw for kw in TOPIC_KEYWORDS[topic] if kw not in text]
    assert not missing, f"{topic} no explica: {missing}"


def test_help_is_not_just_a_docs_pointer():
    for topic in EXPECTED_TOPICS:
        text = _joined(topic, "es")
        assert "ver documentación" not in text
        assert "consulte la documentación" not in text
        assert "see the docs" not in _joined(topic, "en")


def test_get_help_returns_requested_language():
    assert "transportador" in help_ui.get_help("transporte", "es")["summary"]
    assert "transporter" in help_ui.get_help("transporte", "en")["summary"]
    assert "转运蛋白" in help_ui.get_help("transporte", "zh")["summary"]


@pytest.mark.parametrize("topic,lang", [
    ("no_existe", "es"),
    ("transporte", "fr"),
])
def test_get_help_rejects_unknown_topic_or_language(topic, lang):
    with pytest.raises(KeyError):
        help_ui.get_help(topic, lang)


# ---------------------------------------------------------------------------
# Render con Streamlit AppTest
# ---------------------------------------------------------------------------
def test_render_help_does_not_raise_in_apptest():
    at = _apptest(_render_script(
        'render_help("transporte", key_prefix="at_", expanded=True)',
    ))
    assert not at.exception
    assert len(at.expander) == 1
    assert at.expander[0].label == help_ui.EXPANDER_LABEL
    assert [tab.label for tab in at.tabs] == ["Español", "English", "中文"]

    text = "\n".join(markdown.value for markdown in at.markdown)
    assert "transportador" in text
    assert "transporter" in text.lower()
    assert "转运蛋白" in text


def test_render_help_without_prefix_can_repeat_same_topic():
    at = _apptest(_render_script(
        'render_help("membrana")',
        'render_help("membrana", expanded=True)',
    ))
    assert not at.exception
    assert len(at.expander) == 2


def test_render_help_supports_distinct_key_prefixes():
    at = _apptest(_render_script(
        'render_help("operones", key_prefix="a_")',
        'render_help("brechas", key_prefix="b_")',
    ))
    assert not at.exception
    assert len(at.expander) == 2
    labels = [expander.label for expander in at.expander]
    assert labels == [help_ui.EXPANDER_LABEL, help_ui.EXPANDER_LABEL]


def test_render_unknown_topic_warns_without_raising():
    at = _apptest(_render_script('render_help("no_existe")'))
    assert not at.exception
    assert not at.expander
    assert at.warning
    assert "No hay ayuda" in at.warning[0].value
    assert "transporte" in at.warning[0].value


@pytest.mark.parametrize("topic", EXPECTED_TOPICS)
def test_render_every_topic_in_apptest(topic):
    at = _apptest(_render_script(
        f'render_help({topic!r}, key_prefix="all_")',
    ))
    assert not at.exception
    assert len(at.expander) == 1
    assert at.expander[0].label == help_ui.EXPANDER_LABEL
    text = "\n".join(markdown.value for markdown in at.markdown)
    assert help_ui.HELP[topic]["es"]["summary"] in text
