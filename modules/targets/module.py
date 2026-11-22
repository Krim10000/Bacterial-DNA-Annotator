"""Targets esenciales P1/P2/P3 y candidatos in silico."""
from __future__ import annotations

from pathlib import Path

from modules._summary import read_tsv, write_summary
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    if (data_root / "metabolism").exists():
        return ModuleStatus(True)
    return ModuleStatus(False, ("Faltan datos de metabolismo/targets.",))


def _run(ctx: RunContext, _inputs=None) -> ModuleResult:
    base = ctx.data_root / "metabolism" / "targets"
    if not base.exists():
        return ModuleResult("targets", False,
                            error="Falta data/metabolism/targets.")
    files = sorted(base.glob("*.tsv"))
    data = {
        "source": str(base),
        "files": {p.name: len(read_tsv(p)) for p in files},
        "n_files": len(files),
    }
    out, res = write_summary(ctx, "targets", "targets", data)
    res.outputs["source"] = base
    return res


MODULE = Module(
    id="targets",
    name="Targets esenciales",
    summary="Prioriza objetivos esenciales y busca candidatos en la bolsa de "
            "proteínas.",
    order=90,
    category="metabolismo",
    options=(
        Option("priority", "Prioridad", "select", "P1",
               choices=("P1", "P2", "P3"),
               help="P1/P2 son esenciales curados; P3 lleva excepción."),
        Option("top", "Candidatos por target", "int", 10, min=1, max=100,
               help="Cuántos candidatos mostrar por objetivo."),
        Option("include_exceptions", "Incluir excepciones documentadas", "bool", True,
               help="Una excepción no es una exención automática."),
    ),
    docs="docs/AGENT_CHAIN_BUILDER.md",
    run_label="Resumir targets",
    run_fn=_run,
    status_fn=_status,
    estimate_fn=lambda cfg: "depende del número de targets",
)
