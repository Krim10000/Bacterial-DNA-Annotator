"""Extrae las cadenas visibles de la interfaz para el catálogo de traducción.

Uso:
    python tools/extract_ui_strings.py            # escribe webapp/locales/_source.json
    python tools/extract_ui_strings.py --report   # solo cuenta

El catálogo resultante se usa como clave: el texto en español es la clave y
cada idioma aporta su traducción en webapp/locales/<lang>.json.

Criterio: se extraen solo cadenas que llegan a la interfaz (argumentos de
funciones de Streamlit y diccionarios de etiquetas), no nombres de columnas ni
valores de datos. Las f-strings se omiten: sus partes se traducen en el código
con ``tr()`` cuando hace falta.
"""
from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UI_FILES = [
    "webapp/app.py",
    "webapp/credits.py",
    "webapp/browser.py",
    "webapp/esmfold_ui.py",
    "webapp/hypotheses_ui.py",
    "webapp/structural_filters.py",
    "webapp/structure_evidence.py",
    "webapp/requirements_ui.py",
    "webapp/i18n.py",
    "webapp/downloads_ui.py",
    "webapp/module_ui.py",
    "webapp/compare_ui.py",
    "webapp/structure_viewer.py",
    "webapp/gene_card.py",
    "webapp/taxonomy_tree_ui.py",
    "webapp/help_ui.py",
]
TEXTOS = "webapp/textos.py"

# Nombres cuyo valor es una etiqueta de interfaz (aunque se pase por variable).
LABEL_VARS = {"LAYOUTS", "LAYOUT_INTRO", "TOOL_LABEL", "TOOL_DESC",
              "RULE_LABEL", "PIECES", "GROUPS", "LICENSES_NOTE",
              "STEPS", "SYSTEM", "GLOSSARY", "STATUS_GLOSSARY",
              "AUTHOR", "AI_HELP", "HAIRPIN_MEANING", "CAPTION_COMPARE",
              "LIMITS", "HELP", "HELP_TEXTS"}

TEXT_METHODS = {
    "title", "header", "subheader", "caption", "markdown", "info", "warning",
    "error", "success", "button", "checkbox", "radio", "selectbox",
    "multiselect", "text_input", "number_input", "text_area", "download_button",
    "expander", "tabs", "metric", "spinner", "status", "progress", "write",
    "popover", "file_uploader", "toast", "link_button", "form_submit_button",
    "toggle", "slider", "tr",
}

SKIP_PREFIXES = ("<", "#", "http://", "https://", "python ", "streamlit ",
                 "pip ", "bash ", "git ", "--", "RHEA:", "GCF_", "data/",
                 "scripts/", "webapp/", "docs/", "SELECT ", "curl ", "sudo ",
                 "·", ")", ",", ":", "=", "+", "→", "/", "`", "*")
SKIP_CONTAINS = ("</", "/>", "background:", "font-family", "display:",
                 "margin", "padding", "localhost", "-->", "<!--", "\\n",
                 "printf ", "chmod ", "set -a")


def is_candidate(s: str) -> bool:
    s = s.strip()
    if len(s) < 3 or len(s) > 1400:
        return False
    if not re.search(r"[A-Za-zÁÉÍÓÚáéíóúÑñÜü]", s):
        return False
    if s.lower().endswith((".tsv", ".json", ".csv", ".py", ".md", ".sh",
                           ".png", ".txt")):
        return False
    if s.startswith(SKIP_PREFIXES):
        return False
    if any(k in s for k in SKIP_CONTAINS):
        return False
    if re.fullmatch(r"[\w./:-]+", s):
        # identificadores y valores de datos van en minúscula; las etiquetas
        # de una sola palabra ("Idioma", "Contigs", "Tamaño") empiezan mayúscula
        return s[:1].isupper()
    if s.count("`") % 2:
        return False
    return True


class _Collector(ast.NodeVisitor):
    """Constantes de cadena que no forman parte de una f-string."""

    def __init__(self) -> None:
        self.out: set[str] = set()

    def visit_JoinedStr(self, node: ast.JoinedStr) -> None:  # noqa: N802
        return

    def visit_Constant(self, node: ast.Constant) -> None:  # noqa: N802
        if isinstance(node.value, str) and is_candidate(node.value):
            self.out.add(node.value.strip())


def call_strings(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else (
                f.id if isinstance(f, ast.Name) else None)
            if name not in TEXT_METHODS:
                continue
            for arg in node.args:
                c = _Collector()
                c.visit(arg)
                out |= c.out
            for kw in node.keywords:
                if kw.arg in ("help", "label", "text", "placeholder"):
                    c = _Collector()
                    c.visit(kw.value)
                    out |= c.out
    return out


def label_dicts(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for t in targets:
                if isinstance(t, ast.Name) and t.id in LABEL_VARS:
                    c = _Collector()
                    c.visit(node.value)
                    out |= c.out
    return out


def credits_strings(path: Path) -> set[str]:
    """Solo los campos de CREDITS que se traducen (qué hace, uso y licencia).

    Autor, cita, URL y nombre de la herramienta se dejan tal cual: son datos
    bibliográficos y no deben traducirse.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if not any(isinstance(t, ast.Name) and t.id == "CREDITS" for t in targets):
            continue
        for item in ast.walk(node.value):
            if not isinstance(item, ast.Call):
                continue
            for kw in item.keywords:
                if kw.arg in ("what", "use", "license"):
                    c = _Collector()
                    c.visit(kw.value)
                    out |= c.out
    return out


def textos_values(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant):
            if isinstance(node.value.value, str) and len(node.value.value) > 2:
                out.add(node.value.value.strip())
    return out


def module_strings() -> set[str]:
    """Etiquetas, ayudas y textos de los módulos autodescubiertos."""
    out: set[str] = set()
    for path in sorted((ROOT / "modules").glob("*/module.py")):
        c = _Collector()
        c.visit(ast.parse(path.read_text(encoding="utf-8")))
        out |= c.out
    return out


def main() -> int:
    strings: set[str] = set()
    for rel in UI_FILES:
        p = ROOT / rel
        if p.exists():
            strings |= call_strings(p) | label_dicts(p)
    strings |= textos_values(ROOT / TEXTOS)
    strings |= credits_strings(ROOT / "webapp" / "credits.py")
    strings |= module_strings()
    ordered = sorted(strings, key=lambda s: (len(s), s))
    print(f"{len(ordered)} cadenas únicas")
    if "--report" not in sys.argv:
        dest = ROOT / "webapp" / "locales" / "_source.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(ordered, ensure_ascii=False, indent=1),
                        encoding="utf-8")
        print("escrito", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
