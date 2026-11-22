"""Validación contra RegulonDB, DBTBS y benchmarks internos."""
from __future__ import annotations

from pathlib import Path

from modules._summary import read_tsv, write_summary
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    reg = data_root / "regulondb" / "gt" / "gt_operons.tsv"
    if reg.exists():
        return ModuleStatus(True)
    return ModuleStatus(False, ("Falta la verdad de referencia de RegulonDB.",
                                "Genera con `python scripts/validation/regulondb_gt.py`."))


def _run(ctx: RunContext, _inputs=None) -> ModuleResult:
    gt = ctx.data_root / "regulondb" / "gt"
    files = {
        "operons": gt / "gt_operons.tsv",
        "promoters": gt / "gt_promoters.tsv",
        "components": gt / "gt_operon_components.tsv",
        "tf_pwms": gt / "tf_pwms.tsv",
    }
    if not files["operons"].exists():
        return ModuleResult("validation", False,
                            error="Falta data/regulondb/gt/gt_operons.tsv.")
    data = {"source": str(gt), "rows": {}}
    for key, path in files.items():
        data["rows"][key] = len(read_tsv(path)) if path.exists() else None
    out, res = write_summary(ctx, "validation", "regulondb", data)
    res.outputs["source"] = files["operons"]
    return res


MODULE = Module(
    id="validation",
    name="Validación y benchmarks",
    summary="Contrasta las predicciones con RegulonDB, DBTBS y benchmarks "
            "internos.",
    order=100,
    category="validación",
    options=(
        Option("regulondb", "RegulonDB (K-12)", "bool", True,
               help="Verdad curada de operones/promotores de E. coli."),
        Option("dbtbs", "DBTBS (B. subtilis)", "bool", False,
               help="Segunda especie curada; requiere descarga con registro."),
        Option("benchmarks", "Benchmarks internos", "bool", True,
               help="Comparaciones reproducibles del proyecto."),
        Option("strict", "Solo coincidencias estrictas", "bool", False,
               help="Excluye coincidencias parciales; reduce la cobertura."),
    ),
    docs="scripts/validation/README.md",
    run_label="Resumir validación",
    run_fn=_run,
    status_fn=_status,
    estimate_fn=lambda cfg: "2–10 min según el conjunto",
)
