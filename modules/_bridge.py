"""Puente entre los módulos y el paquete `webapp`.

La app ejecuta `streamlit run webapp/app.py` con `webapp/` en `sys.path`; los
módulos viven fuera. Este helper reproduce esa condición sin duplicar lógica.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEBAPP = ROOT / "webapp"
for path in (str(WEBAPP), str(ROOT / "scripts" / "ghost")):
    if Path(path).is_dir() and path not in sys.path:
        sys.path.insert(0, path)
