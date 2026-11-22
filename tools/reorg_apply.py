"""Aplica la reorganización de `scripts/` descrita en el manifiesto.

Por defecto **no modifica nada**: muestra el plan. Con ``--confirm`` mueve las
carpetas, resuelve la colisión de `orchestrator.py` y reescribe las referencias.

Uso:
    python tools/reorg_apply.py            # simulación
    python tools/reorg_apply.py --confirm  # aplicar

Requiere que `python tools/audit_script_paths.py` esté en 0 antes de empezar.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MANIFEST = ROOT / "tools" / "reorg_manifest.csv"
REPORT = ROOT / ".scratch" / "reorg_applied.json"
SKIP_DIRS = {".git", ".venv", ".venv_equilibrator", "venv", "dist", "data",
             "logs", "cache", ".cache", ".scratch", "__pycache__", "foldseek",
             "checkpoints_descargados", "out", "out_fl"}
TEXT_SUFFIXES = {".py", ".sh", ".md", ".yml", ".yaml", ".toml", ".txt", ".cfg"}
TEXT_NAMES = {"Dockerfile", "Makefile"}
# El plan de reorganización cita rutas futuras: no se reescribe.
SKIP_FILES = {"docs/SCRIPTS_REORG_PLAN.md"}


def read_manifest() -> list[dict]:
    with MANIFEST.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def planned_rows(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r["action"] in {"move", "merge", "move+merge"}]


def text_files() -> list[Path]:
    out = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(p in SKIP_DIRS for p in path.parts):
            continue
        if path.suffix not in TEXT_SUFFIXES and path.name not in TEXT_NAMES:
            continue
        if str(path.relative_to(ROOT)) in SKIP_FILES:
            continue
        out.append(path)
    return out


def check_collisions(moves: list[dict]) -> list[str]:
    problems = []
    for row in moves:
        src = SCRIPTS / row["source_dir"]
        dst = SCRIPTS / row["target_dir"]
        if not src.exists():
            continue
        for item in src.iterdir():
            if item.name == "__pycache__":
                continue
            target = dst / item.name
            if (dst / item.name).exists():
                if item.name == "orchestrator.py" and row["source_dir"] == "04_defense":
                    continue
                if item.name != "orchestrator.py":
                    # README/run_*.sh: se renombran al mover, no son un error.
                    continue
                problems.append(f"{src.name}/{item.name} -> {dst.name}/{item.name}")
    return problems


def move_tree(row: dict, log: list[str]) -> None:
    src = SCRIPTS / row["source_dir"]
    dst = SCRIPTS / row["target_dir"]
    if not src.exists():
        log.append(f"(omitido, ya movido) {row['source_dir']}")
        return
    dst.mkdir(parents=True, exist_ok=True)
    for item in list(src.iterdir()):
        if item.name == "__pycache__":
            shutil.rmtree(item, ignore_errors=True)
            continue
        target = dst / item.name
        if target.exists() and item.name == "orchestrator.py" and row["source_dir"] == "04_defense":
            target = dst / "orchestrator_defense.py"
        if target.exists():
            # README/run_*.sh de módulos que se fusionan: se conservan ambos.
            target = dst / f"{target.stem}_{row['source_dir']}{target.suffix}"
            if target.exists():
                raise RuntimeError(f"colisión doble: {target}")
        shutil.move(str(item), str(target))
        log.append(f"{row['source_dir']}/{item.name} -> {row['target_dir']}/{target.name}")
    try:
        src.rmdir()
    except OSError:
        pass


def rewrite(moves: list[dict]) -> Counter:
    # Reemplazos ordenados de más largo a más corto para no partir prefijos.
    pairs = sorted(
        ((f"scripts/{r['source_dir']}", f"scripts/{r['target_dir']}") for r in moves),
        key=lambda p: len(p[0]), reverse=True,
    )
    counts: Counter = Counter()
    touched = set()
    for path in text_files():
        try:
            text = path.read_text(encoding="utf-8", errors="strict")
        except (UnicodeDecodeError, OSError):
            continue
        new = text
        for old, target in pairs:
            if old in new:
                n = new.count(old)
                new = new.replace(old, target)
                counts[f"{old} -> {target}"] += n
                touched.add(str(path.relative_to(ROOT)))
        if new != text:
            path.write_text(new, encoding="utf-8")
    print(f"Reescrituras: {sum(counts.values())} en {len(touched)} archivos")
    for key in sorted(counts, key=lambda k: -counts[k]):
        print(f"  {counts[key]:4d}  {key}")
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--confirm", action="store_true",
                    help="aplicar de verdad (sin esto solo simula)")
    args = ap.parse_args()

    rows = read_manifest()
    moves = planned_rows(rows)
    print(f"Movimientos/fusiones: {len(moves)}")
    problems = check_collisions(moves)
    if problems:
        print("COLISIONES:")
        for p in problems:
            print("  ", p)
        return 2
    print("Sin colisiones salvo la conocida de annotation/orchestrator.py.")

    if not args.confirm:
        for r in moves:
            print(f"  {r['source_dir']:18s} -> {r['target_dir']:16s} "
                  f"({r['py_files']} py, {r['action']})")
        print("\nSimulación. Usa --confirm para aplicar.")
        return 0

    log: list[str] = []
    for row in moves:
        move_tree(row, log)
    print(f"Movidos {len(log)} elementos.")
    counts = rewrite(moves)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps({
        "moves": moves, "moved": log,
        "rewrites": dict(counts), "total_rewrites": sum(counts.values()),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Informe: {REPORT}")
    print("Ejecuta ahora: python tools/audit_script_paths.py && "
          "python -m compileall -q scripts webapp && python -m pytest tests -q")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
