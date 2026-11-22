"""Simulación de la reorganización de `scripts/` (no modifica archivos).

Lee `tools/reorg_manifest.csv`, comprueba colisiones de nombres y cuenta las
referencias que habría que reescribir. Es la herramienta de la fase 3 del plan
`docs/SCRIPTS_REORG_PLAN.md`.

Uso:
    python tools/reorg_dry_run.py
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tools" / "reorg_manifest.csv"
PATTERN = re.compile(r"scripts/[A-Za-z0-9_./-]+")
SKIP = {".git", ".venv", ".venv_equilibrator", "venv", "dist", "data", "logs",
        "cache", ".cache", ".scratch", "__pycache__", "foldseek",
        "checkpoints_descargados"}
SKIP_FILES = {"docs/SCRIPTS_REORG_PLAN.md"}


def read_manifest() -> list[dict]:
    with MANIFEST.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def text_files() -> list[Path]:
    out = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(p in SKIP for p in path.relative_to(ROOT).parts):
            continue
        if path.suffix not in {".py", ".sh", ".md", ".yml", ".yaml", ".toml",
                               ".txt", ".cfg"} and path.name not in {
                                   "Dockerfile", "Makefile"}:
            continue
        if str(path.relative_to(ROOT)) in SKIP_FILES:
            continue
        out.append(path)
    return out


def main() -> int:
    rows = read_manifest()
    print(f"Manifiesto: {len(rows)} orígenes")
    moves = [(r["source_dir"], r["target_dir"]) for r in rows
             if r["action"] != "keep"]
    print(f"Movimientos/ fusiones: {len(moves)}")

    # 1. Colisiones de nombres por destino
    by_target: dict[str, list[str]] = {}
    for src, dst in moves:
        for f in (ROOT / "scripts" / src).glob("*.py"):
            by_target.setdefault(dst, []).append(f.name)
    print("\nColisiones de nombre por destino:")
    any_collision = False
    for dst, names in sorted(by_target.items()):
        dup = [n for n, c in Counter(names).items() if c > 1]
        if dup:
            any_collision = True
            print(f"  {dst}: {sorted(dup)}")
    if not any_collision:
        print("  ninguna")

    # 2. Referencias que cambiarían
    prefixes = sorted(
        ((f"scripts/{src}/", f"scripts/{dst}/") for src, dst in moves),
        key=lambda p: len(p[0]), reverse=True,
    )
    counts: Counter[str] = Counter()
    files: Counter[str] = Counter()
    touched: set[str] = set()
    for path in text_files():
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for old, new in prefixes:
            n = text.count(old)
            if n:
                counts[f"{old} -> {new}"] += n
                files[f"{old} -> {new}"] += 1
                touched.add(str(path.relative_to(ROOT)))
    print("\nReescrituras previstas (ocurrencias / archivos):")
    for key in sorted(counts, key=lambda k: -counts[k]):
        print(f"  {counts[key]:4d} / {files[key]:3d}  {key}")

    total = sum(counts.values())
    print(f"\nTotal: {total} ocurrencias en {len(touched)} archivos")
    print("No se ha modificado ningún archivo (modo simulación).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
