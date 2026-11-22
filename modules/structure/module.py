"""Estructura: Foldseek, P2Rank, AlphaFill y ESMFold."""
from __future__ import annotations

import shutil
from pathlib import Path

from modules.base import Module, ModuleStatus, Option


def _status(data_root: Path):
    tools = {t: shutil.which(t) is not None for t in ("foldseek", "prank")}
    if not any(tools.values()):
        return ModuleStatus(False, ("No se encontraron Foldseek ni P2Rank.",),
                            {"tools": tools})
    return ModuleStatus(True, detail={"tools": tools})


def _selftest() -> str:
    from structural_filters import run_filters
    from structure_evidence import evidence_for_rows
    assert callable(run_filters) and callable(evidence_for_rows)
    return "filtros estructurales y evidencia OK"


MODULE = Module(
    id="structure",
    name="Estructura y filtros",
    summary="Homología 3D, bolsillos, ligandos y sitio activo para priorizar "
            "candidatos.",
    order=40,
    category="estructura",
    options=(
        Option("foldseek", "F1 · Foldseek", "bool", True,
               help="Busca parecidos estructurales 3D."),
        Option("enzymm", "F2 · EnzyMM / M-CSA", "bool", True,
               help="Plantillas geométricas de residuos catalíticos."),
        Option("p2rank", "F3 · P2Rank", "bool", True,
               help="Predice bolsillos de unión."),
        Option("alphafill", "F4 · AlphaFill", "bool", True,
               help="Transfiere ligandos/cofactores desde PDB."),
        Option("active_site", "F5 · Sitio activo", "bool", True,
               help="Residuos catalíticos de UniProt/M-CSA."),
        Option("esmfold", "ESMFold (GPU o API)", "bool", False,
               help="Predice la estructura del ORF exacto; opcional y de pago."),
    ),
    docs="docs/FIVE_STRUCTURAL_FILTERS.md",
    selftest_fn=_selftest,
    status_fn=_status,
    estimate_fn=lambda cfg: "1–5 min por proteína (ESMFold: 1–2 min y crédito)",
)
