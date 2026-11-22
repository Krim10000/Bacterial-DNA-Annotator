"""Audita las referencias a rutas de `scripts/` en el repositorio.

Sirve antes y después de una reorganización: si un script, shell, doc o
workflow menciona `scripts/<algo>` y ese algo no existe, falla.

Uso:
    python tools/audit_script_paths.py
    python tools/audit_script_paths.py --quiet
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERN = re.compile(r"scripts/[A-Za-z0-9_./-]+")
SPLIT = re.compile(r"""["']scripts["']\s*/\s*["']([A-Za-z0-9_.-]+)["']""")
SKIP_DIRS = {".git", ".venv", ".venv_equilibrator", "venv", "dist", "data",
             "logs", "cache", ".cache", ".scratch", "__pycache__", "foldseek",
             "out", "out_fl",
             "checkpoints_descargados"}
SUFFIXES = {".py", ".sh", ".md", ".yml", ".yaml", ".toml", ".txt", ".cfg"}
NAMES = {"Dockerfile", "Makefile"}
# El plan de reorganización cita rutas *futuras*: no son referencias reales.
SKIP_FILES = {"docs/SCRIPTS_REORG_PLAN.md"}

# Referencias que no son rutas del repositorio (se crean en tiempo de ejecución
# o son fragmentos de un glob/plantilla).
ALLOW = {
    "scripts/data_prep/0",       # glob 0*_*.py en run_pipeline.sh
    "scripts/batter-tpe",           # BATTER descargado en data/terminators
    "scripts/crf.py",               # dentro del paquete BATTER descargado
    "scripts/constants.py",         # ídem
    "scripts/dataset.py",           # ídem
    "scripts/model.py",             # ídem
    "scripts/testdata",             # datos de prueba locales no versionados
}


def files() -> list[Path]:
    out = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        try:
            rel_parts = path.relative_to(ROOT).parts
        except ValueError:
            continue
        if any(part in SKIP_DIRS for part in rel_parts):
            continue
        if path.suffix not in SUFFIXES and path.name not in NAMES:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if str(path.relative_to(ROOT)) in SKIP_FILES:
            continue
        out.append(path)
    return out


def main() -> int:
    missing: dict[str, list[str]] = {}
    checked = 0
    for path in files():
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for comp in SPLIT.findall(text):
            token = f"scripts/{comp}"
            if token in ALLOW:
                continue
            checked += 1
            if not (ROOT / token).exists():
                rel = str(path.relative_to(ROOT))
                missing.setdefault(token, []).append(rel)
        for raw in PATTERN.findall(text):
            token = raw.rstrip(".,;:)\'\"`")
            if token in ALLOW:
                continue
            checked += 1
            if not (ROOT / token).exists():
                rel = str(path.relative_to(ROOT))
                missing.setdefault(token, []).append(rel)

    if not missing:
        print(f"OK: {checked} referencias a scripts/ existen.")
        return 0

    print(f"ERROR: {len(missing)} rutas referenciadas no existen "
          f"(de {checked} referencias):")
    for token in sorted(missing):
        where = ", ".join(sorted(set(missing[token]))[:4])
        print(f"  {token}  <-  {where}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
