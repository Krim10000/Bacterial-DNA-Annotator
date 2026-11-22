"""Membrana: hélices transmembrana, señal péptido, lipoproteínas y transporte."""
from __future__ import annotations

from pathlib import Path

from modules._summary import counts, read_tsv, write_summary
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    if (data_root / "membrane").exists():
        return ModuleStatus(True)
    return ModuleStatus(False, ("Faltan datos de membrana.",
                                "Genera con `python scripts/membrane/06_run_full.py`."))


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    path = ctx.data_root / "membrane" / f"{genome}.tsv"
    if not path.exists():
        return ModuleResult("membrane", False,
                            error=f"Sin datos de membrana para {genome}.")
    rows = read_tsv(path)
    data = {
        "genome": genome,
        "source": str(path.relative_to(ctx.data_root)),
        "rows": len(rows),
        "by_type": counts(rows, "type"),
        "by_confidence": counts(rows, "confidence"),
    }
    out, res = write_summary(ctx, "membrane", genome, data)
    res.outputs["source"] = path
    return res


MODULE = Module(
    id="membrane",
    name="Membrana y transporte",
    summary="Clasifica proteínas de membrana, señal péptido, lipoproteínas y "
            "transportadores.",
    order=80,
    category="genómica",
    options=(
        Option("tm", "Hélices transmembrana", "bool", True,
               help="Hidrofobicidad tipo Kyte-Doolittle; heurística, no TMHMM."),
        Option("signal", "Señal péptido", "bool", True,
               help="Predicción heurística; no sustituye SignalP."),
        Option("lipoprotein", "Lipoproteínas", "bool", True,
               help="Motivo lipobox y procesamiento."),
        Option("transporters", "Transportadores y sustrato", "bool", True,
               help="Solo asigna sustrato cuando es explícito."),
    ),
    docs="docs/AGENT_MEMBRANE_BUILDER.md",
    input_kind="genome",
    run_label="Resumir membrana",
    run_fn=_run,
    status_fn=_status,
    estimate_fn=lambda cfg: "~5 min por genoma",
)
