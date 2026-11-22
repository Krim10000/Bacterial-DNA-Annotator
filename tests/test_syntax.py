"""Todo el código Python del repositorio debe compilar.

Se ejecuta sobre ``webapp/`` y ``scripts/`` para detectar errores de sintaxis
que no cubren los tests funcionales (por ejemplo, un docstring mal puesto antes
de ``from __future__``).
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _files() -> list[Path]:
    out = []
    for base in ("webapp", "scripts", "tools", "tests"):
        for path in (ROOT / base).rglob("*.py"):
            if "__pycache__" in path.parts or "scripts_orginales" in path.parts:
                continue
            out.append(path)
    return out


def test_all_python_files_parse():
    errors = []
    for path in _files():
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            errors.append(f"{path.relative_to(ROOT)}:{exc.lineno}: {exc.msg}")
    assert not errors, "Errores de sintaxis:\n" + "\n".join(errors[:20])
