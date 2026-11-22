"""ORFs no anotados, desiertos y cierres por dominio."""
from __future__ import annotations

from modules.base import Module, ModuleStatus, Option


def _status(data_root):
    return ModuleStatus(True)


def _selftest() -> str:
    from ghost.orfs import merge_intervals
    from ghost.seqio import read_fasta
    assert callable(merge_intervals) and callable(read_fasta)
    return "motor de genes fantasma OK"


MODULE = Module(
    id="ghost",
    name="Genes fantasma",
    summary="Busca genes en los huecos que nadie anotó y regiones desiertas.",
    order=60,
    category="anotación",
    options=(
        Option("min_region_len", "Longitud mínima de región (nt)", "int", 150,
               min=30, max=5000, step=10,
               help="Huecos más cortos generan ruido."),
        Option("min_aa", "Longitud mínima de ORF (aa)", "int", 60,
               min=20, max=500, step=5,
               help="Por debajo de 60 aa es difícil distinguir un gen real."),
        Option("esm", "Puntuar con ESM-2", "bool", False,
               help="Modelo de lenguaje de proteínas; corre en CPU."),
    ),
    docs="scripts/ghost/README.md",
    status_fn=_status,
    selftest_fn=_selftest,
    estimate_fn=lambda cfg: "~5 s por genoma; ~6 h los 400",
)
