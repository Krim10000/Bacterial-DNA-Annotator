"""Utilidades para los módulos que resumen datos ya calculados.

No son un recalculo: leen tablas/JSON del proyecto y devuelven conteos
verificables. Cada salida declara el archivo de origen.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from modules.base import ModuleResult, RunContext


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def counts(rows: list[dict], column: str) -> dict[str, int]:
    values = [str(r.get(column, "")).strip() or "sin_dato" for r in rows]
    return dict(sorted(Counter(values).items(), key=lambda kv: (-kv[1], kv[0])))


def write_summary(ctx: RunContext, module: str, name: str,
                  data: dict) -> tuple[Path, ModuleResult]:
    out = ctx.out_root / module
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}_summary.json"
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return path, ModuleResult(module, True, outputs={"summary": path},
                              metrics=data)
