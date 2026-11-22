#!/usr/bin/env python3
"""doctor.py — verifica que el entorno pueda ejecutar el pipeline.

Pensado para un usuario nuevo: dice qué falta y para qué etapa se necesita,
sin lanzar errores crípticos a mitad de un proceso largo.

Uso
---
    python scripts/00_doctor.py
"""
from __future__ import annotations

import importlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TOOLS = [
    ("diamond", "anotación funcional y homología"),
    ("hmmscan", "dominios Pfam (cierres, sigma)"),
    ("hmmfetch", "extraer paneles de Pfam"),
    ("hmmpress", "indexar paneles de Pfam"),
    ("cmsearch", "rRNA (Infernal)"),
    ("tRNAscan-SE", "tRNA"),
    ("barrnap", "rRNA (alternativa rápida)"),
]

MODULES = [
    ("numpy", "núcleo", True),
    ("pandas", "tablas y GUI", True),
    ("scipy", "estadística", True),
    ("sklearn", "evaluación y modelos", True),
    ("networkx", "grafos metabólicos", True),
    ("sentencepiece", "tokenizador", False),
    ("streamlit", "GUI web", False),
    ("torch", "modelo de ADN (opcional)", False),
    ("RNA", "ViennaRNA: estructura de ARN", False),
]

DATA = [
    ("data/test_ecoli.gbff", "genoma de referencia K-12"),
    ("data/regulondb/gt/gt_operons.tsv", "verdad de referencia"),
    ("data/db/swissprot.dmnd", "anotación funcional"),
    ("data/annotations_ab_initio", "400 genomas anotados"),
]


def main() -> int:
    ok = True
    print("== herramientas externas ==")
    for tool, why in TOOLS:
        path = shutil.which(tool)
        mark = "OK   " if path else "FALTA"
        if not path:
            ok = False
        print(f"  {tool:14s} {mark}  ({why})")

    print("\n== módulos de Python ==")
    for mod, why, required in MODULES:
        try:
            importlib.import_module(mod)
            print(f"  {mod:14s} OK     ({why})")
        except Exception:
            print(f"  {mod:14s} FALTA  ({why})" + ("" if not required else "  [requerido]"))
            if required:
                ok = False

    print("\n== datos del proyecto ==")
    for rel, why in DATA:
        p = ROOT / rel
        print(f"  {rel:38s} {'OK' if p.exists() else 'FALTA'}  ({why})")

    print("\n== salidas que deberían existir tras correr el pipeline ==")
    outs = [
        ("data/operon_finder/ecoli_k12_operons_rules.tsv", "operones por reglas"),
        ("data/operon_finder/ecoli_k12_hairpins.tsv", "horquillas"),
        ("data/regulondb/gt/gt_operon_components.tsv", "componentes por operón"),
        ("data/metabolism/closure_final_top1000.tsv", "cierres de brechas"),
    ]
    for rel, why in outs:
        p = ROOT / rel
        print(f"  {rel:52s} {'OK' if p.exists() else 'pendiente'}  ({why})")

    print("\ninstalación rápida:")
    print("  python -m venv .venv && .venv/bin/pip install -r requirements.txt")
    print("  bash run_pipeline.sh components && bash run_pipeline.sh pwms")
    print("  bash run_pipeline.sh hairpins && bash run_pipeline.sh closure")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
