"""Identidad taxonómica: ANI, GTDB y especie más cercana."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from modules._summary import write_summary
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    if shutil.which("fastANI") is None:
        return ModuleStatus(False, ("Falta fastANI.",))
    if not (data_root / "taxonomy").exists():
        return ModuleStatus(False, ("Faltan datos de taxonomía del proyecto.",
                                    "Genera con `python scripts/taxonomy/07_assign_taxonomy_full.py`."))
    return ModuleStatus(True)


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    path = ctx.data_root / "taxonomy" / f"{genome}.json"
    if not path.exists():
        return ModuleResult("taxonomy", False,
                            error=f"Sin datos de taxonomía para {genome}.")
    data = json.loads(path.read_text(encoding="utf-8"))
    summary = {
        "genome": genome,
        "source": str(path.relative_to(ctx.data_root)),
        "status": data.get("status"),
        "species": data.get("species"),
        "closest_ref": data.get("closest_ref"),
        "ani": data.get("ani"),
        "confidence": data.get("confidence"),
        "gtdb_lineage": data.get("gtdb_lineage"),
    }
    out, res = write_summary(ctx, "taxonomy", genome, summary)
    res.outputs["source"] = path
    return res


MODULE = Module(
    id="taxonomy",
    name="Taxonomía y ANI",
    summary="Calcula el parecido global (ANI), la línea GTDB y la especie más "
            "cercana.",
    order=70,
    category="genómica",
    options=(
        Option("ani", "ANI con fastANI", "bool", True,
               help="Parecido global entre genomas; no es una probabilidad."),
        Option("gtdb", "Línea GTDB", "bool", True,
               help="Taxonomía GTDB cuando hay referencia disponible."),
        Option("ani_threshold", "Umbral de especie (ANI)", "float", 95.0,
               min=80.0, max=100.0, step=0.5,
               help="95 % es la convención habitual, no una ley universal."),
        Option("top", "Referencias a comparar", "int", 5, min=1, max=20,
               advanced=True),
    ),
    docs="docs/AGENT_TAXONOMY_BUILDER.md",
    input_kind="genome",
    run_label="Resumir taxonomía",
    run_fn=_run,
    status_fn=_status,
    estimate_fn=lambda cfg: "1–5 min por genoma",
)
