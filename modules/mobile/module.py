"""Elementos móviles, CRISPR, AMR y profagos."""
from __future__ import annotations

import json
from pathlib import Path

from modules._summary import counts, read_tsv, write_summary
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    if (data_root / "mobile_typing").exists():
        return ModuleStatus(True)
    return ModuleStatus(False, ("Faltan datos de móviles tipados.",
                                "Genera con `bash scripts/mobile/run_all.sh`."))


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    path = ctx.data_root / "mobile_typing" / f"{genome}.tsv"
    if not path.exists():
        return ModuleResult("mobile", False,
                            error=f"Sin datos de móviles para {genome}.")
    rows = read_tsv(path)
    data = {
        "genome": genome,
        "source": str(path.relative_to(ctx.data_root)),
        "rows": len(rows),
        "by_element": counts(rows, "element_type"),
        "by_subtype": counts(rows, "subtype"),
        "by_confidence": counts(rows, "confidence"),
    }
    out, res = write_summary(ctx, "mobile", genome, data)
    res.outputs["source"] = path
    return res


MODULE = Module(
    id="mobile",
    name="Móviles y defensa",
    summary="Tipado de transposones, CRISPR, AMR y profagos.",
    order=50,
    category="genómica",
    options=(
        Option("cas", "Sistemas CRISPR-Cas", "bool", True),
        Option("crispr", "Arrays y espaciadores", "bool", True),
        Option("amr", "Resistencia a antimicrobianos", "bool", True),
        Option("prophage", "Profagos", "bool", True),
    ),
    docs="docs/AGENT_MOBILE_BUILDER.md",
    input_kind="genome",
    run_label="Resumir móviles",
    run_fn=_run,
    status_fn=_status,
    estimate_fn=lambda cfg: "~25 min por barrido de 10 genomas",
)
