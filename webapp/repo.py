"""Datos públicos del repositorio.

Cambia ``REPO_URL`` si el nombre final del repositorio es otro. El instalador y
la página de créditos lo usan para el enlace de «dale una estrella».
"""
from __future__ import annotations

import os

REPO_URL = os.environ.get(
    "BACTO_SCOPE_REPO",
    "https://github.com/Krim10000/Bacterial-DNA-Annotator",
)
STAR_URL = REPO_URL.rstrip("/") + "/stargazers"
