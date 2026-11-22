"""Internacionalización (i18n) de la interfaz.

Diseño
------
- El **texto en español es la clave** del catálogo. Así no hay que mantener
  identificadores paralelos y el código sigue siendo legible.
- Cada idioma vive en ``webapp/locales/<código>.json`` con el formato
  ``{"texto en español": "traducción"}``.
- ``tr()`` traduce una cadena exacta; si no hay traducción, devuelve el
  original. Nunca falla ni deja la interfaz en blanco.
- ``T`` (en ``textos.py``) resuelve sus claves en el idioma activo, de modo que
  las llamadas existentes ``T["clave"]`` siguen funcionando.
- ``install_streamlit_i18n()`` envuelve los métodos de Streamlit para traducir
  etiquetas y textos sin tocar las ~500 llamadas de la interfaz. Es un parche
  acotado: solo actúa sobre argumentos de texto, nunca sobre datos.

Para añadir un idioma: copia ``en.json``, traduce los valores y añade el código
a ``LANGUAGES``. El menú se construye solo.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

LOCALES_DIR = Path(__file__).resolve().parent / "locales"
DEFAULT_LANGUAGE = "es"

#: Código ISO 639-1 → nombre que se muestra en el menú (en su propio idioma).
LANGUAGES: dict[str, str] = {
    "es": "Español",
    "en": "English",
    "zh": "中文",
}

_STATE: dict[str, Any] = {"lang": None}
_CATALOGS: dict[str, dict[str, str]] = {}
_PATCHED = False


# ------------------------------------------------------------------ catálogo
def _catalog(lang: str) -> dict[str, str]:
    if lang == DEFAULT_LANGUAGE:
        return {}
    if lang not in _CATALOGS:
        path = LOCALES_DIR / f"{lang}.json"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            _CATALOGS[lang] = {str(k): str(v) for k, v in data.items() if v}
        except (OSError, ValueError):
            _CATALOGS[lang] = {}
    return _CATALOGS[lang]


def catalog_coverage(lang: str, source: list[str] | None = None) -> float:
    """Fracción de cadenas de la fuente que tienen traducción (para tests/CI)."""
    if lang == DEFAULT_LANGUAGE:
        return 1.0
    cat = _catalog(lang)
    if source is None:
        path = LOCALES_DIR / "_source.json"
        source = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    if not source:
        return 0.0
    return sum(1 for s in source if cat.get(s)) / len(source)


# ----------------------------------------------------------------- idioma
def _load_saved_language() -> str | None:
    path = Path(__file__).resolve().parents[1] / ".cache" / "ui_config.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        lang = data.get("lang")
        return lang if lang in LANGUAGES else None
    except (OSError, ValueError):
        return None


def _save_language(lang: str) -> None:
    path = Path(__file__).resolve().parents[1] / ".cache" / "ui_config.json"
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
        data["lang"] = lang
        path.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    except OSError:
        pass


def init_language() -> str:
    """Fija el idioma al arrancar: ``?lang=`` > sesión > archivo > español."""
    try:
        import streamlit as st

        qp = st.query_params
        lang = qp.get("lang") if hasattr(qp, "get") else None
        if lang in LANGUAGES:
            st.session_state["lang"] = lang
        elif st.session_state.get("lang") not in LANGUAGES:
            st.session_state["lang"] = _load_saved_language() or DEFAULT_LANGUAGE
        _STATE["lang"] = st.session_state["lang"]
    except Exception:  # sin Streamlit (tests, CLI)
        _STATE["lang"] = _STATE["lang"] or _load_saved_language() or DEFAULT_LANGUAGE
    return get_language()


def get_language() -> str:
    try:
        import streamlit as st

        lang = st.session_state.get("lang")
        if lang in LANGUAGES:
            return lang
    except Exception:
        pass
    if _STATE["lang"] in LANGUAGES:
        return _STATE["lang"]
    return DEFAULT_LANGUAGE


def set_language(lang: str) -> str:
    if lang not in LANGUAGES:
        return get_language()
    _STATE["lang"] = lang
    try:
        import streamlit as st

        st.session_state["lang"] = lang
    except Exception:
        pass
    _save_language(lang)
    return lang


# --------------------------------------------------------------- traducción
def tr(text: Any) -> Any:
    """Traduce una cadena exacta; si no es cadena o no hay entrada, la deja."""
    if not isinstance(text, str) or not text:
        return text
    lang = get_language()
    if lang == DEFAULT_LANGUAGE:
        return text
    return _catalog(lang).get(text, text)


def tr_tree(value: Any) -> Any:
    """Traduce cadenas dentro de listas, tuplas y dicts (etiquetas de menús)."""
    if isinstance(value, str):
        return tr(value)
    if isinstance(value, list):
        return [tr_tree(v) for v in value]
    if isinstance(value, tuple):
        return tuple(tr_tree(v) for v in value)
    if isinstance(value, dict):
        return {k: tr_tree(v) for k, v in value.items()}
    return value


# --------------------------------------------------------------- Streamlit
_LABEL_ARG0 = {
    "markdown", "caption", "subheader", "header", "title", "info", "warning",
    "error", "success", "write", "toast", "spinner", "button", "download_button",
    "link_button", "form_submit_button", "checkbox", "toggle", "radio",
    "selectbox", "multiselect", "text_input", "number_input", "text_area",
    "file_uploader", "expander", "popover", "slider", "metric", "status",
}
_LIST_ARG0 = {"tabs"}
_KW_LABEL = {"help", "label", "text", "placeholder"}


def _translate_args(args: tuple, kwargs: dict, *, list_arg0: bool) -> tuple:
    args = list(args)
    if args and isinstance(args[0], str):
        args[0] = tr(args[0])
    elif args and list_arg0 and isinstance(args[0], (list, tuple)):
        args[0] = tr_tree(args[0])
    for key in kwargs:
        if key in _KW_LABEL and isinstance(kwargs[key], str):
            kwargs[key] = tr(kwargs[key])
    return tuple(args)


def install_streamlit_i18n() -> None:
    """Envuelve los métodos de ``DeltaGenerator`` una sola vez.

    Se parchea la clase (no ``st``), porque ``st.sidebar`` es otra instancia de
    ``DeltaGenerator`` y de otro modo sus etiquetas no se traducirían.
    """
    global _PATCHED
    if _PATCHED:
        return
    try:
        from streamlit.delta_generator import DeltaGenerator
    except Exception:
        return

    for name in _LABEL_ARG0 | _LIST_ARG0:
        original = getattr(DeltaGenerator, name, None)
        if original is None:
            continue

        def make(orig, list_arg0: bool):
            def wrapper(self, *args, **kwargs):
                args = _translate_args(args, kwargs, list_arg0=list_arg0)
                return orig(self, *args, **kwargs)

            wrapper.__name__ = getattr(orig, "__name__", "wrapper")
            wrapper.__doc__ = getattr(orig, "__doc__", None)
            return wrapper

        setattr(DeltaGenerator, name, make(original, name in _LIST_ARG0))

    # ``st.header`` y compañía son métodos ya ligados guardados como atributos
    # del módulo: cambiar la clase no los afecta. Se parchean también aquí.
    try:
        import streamlit as st

        for name in _LABEL_ARG0 | _LIST_ARG0:
            original = getattr(st, name, None)
            if original is None or not callable(original):
                continue

            def make_module(orig, list_arg0: bool):
                def wrapper(*args, **kwargs):
                    args = _translate_args(args, kwargs, list_arg0=list_arg0)
                    return orig(*args, **kwargs)

                wrapper.__name__ = getattr(orig, "__name__", "wrapper")
                return wrapper

            setattr(st, name, make_module(original, name in _LIST_ARG0))
    except Exception:
        pass
    _PATCHED = True


# ------------------------------------------------------------ menú de idioma
def language_selector() -> str:
    """Dibuja el selector de idioma en la barra lateral y devuelve el código."""
    import streamlit as st

    codes = list(LANGUAGES)
    current = get_language()
    chosen = st.sidebar.selectbox(
        "🌐 " + tr("Idioma"),
        codes,
        index=codes.index(current),
        format_func=lambda c: LANGUAGES[c],
        key="lang_select",
        help=tr("Cambia el idioma de la interfaz. Los datos y los resultados "
                "no se traducen: se muestran tal como los produce el pipeline."),
    )
    if chosen != current:
        set_language(chosen)
        st.rerun()
    return chosen
