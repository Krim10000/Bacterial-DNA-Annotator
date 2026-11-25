"""Sigma: factores sigma y promotores sigma-dependientes.

El módulo no importa Streamlit. La GUI genérica lo descubre y ejecuta sobre el
genoma seleccionado; las salidas incluyen la evidencia por anotación.
"""
from __future__ import annotations

import json
from pathlib import Path

from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path) -> ModuleStatus:
    pwm = data_root / "regulondb" / "gt" / "promoter_pwms.tsv"
    reg = data_root / "regulation"
    if not pwm.exists():
        return ModuleStatus(False, (
            "Falta data/regulondb/gt/promoter_pwms.tsv.",
            "Genera las PWMs con `bash run_pipeline.sh pwms`.",
        ))
    if not reg.exists() or not any((reg / d / "genes.tsv").exists()
                                   for d in (p.name for p in reg.iterdir() if p.is_dir())):
        return ModuleStatus(False, ("Faltan genomas en data/regulation/.",))
    return ModuleStatus(True, detail={"pwms": str(pwm), "genomes": str(reg)})


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    from scripts.sigma.sigma_annotator import annotate_genome

    if not genome:
        return ModuleResult("sigma", False, error="Selecciona un genoma.")
    try:
        res = annotate_genome(
            ctx.data_root, genome,
            ctx.out_root / "sigma" / genome,
            window=int(ctx.config.get("window", 180)),
            null_windows=int(ctx.config.get("null_windows", 250)),
            seed=int(ctx.config.get("seed", 17)),
            max_p=float(ctx.config.get("max_p", 0.05)),
            min_z=float(ctx.config.get("min_z", 2.0)),
            min_margin=float(ctx.config.get("min_margin", 1.0)),
            use_expression=bool(ctx.config.get("use_expression", True)),
        )
    except Exception as exc:  # noqa: BLE001
        return ModuleResult("sigma", False, error=str(exc))
    outputs = {
        "promoters": res["promoters"],
        "factors": res["factors"],
        "report": res["report"],
        "preview": res["promoters"],
        "summary": Path(res["out_dir"]) / "sigma_summary.json",
    }
    metrics = dict(res["summary"])
    warnings = [
        "Las anotaciones son hipótesis: `sospecha_fuerte` no es confirmación experimental.",
        "La expresión es opcional; si no hay tabla, se declara `not_available`.",
        "sigma54 se escanea con geometría -24/-12, no -35/-10.",
        "La agrupación por ortólogo es global: `python scripts/sigma/20_group_sigma.py --all`.",
        "Gold K-12: top-1 57,0 % (mayoría trivial 62,6 %); top-3 90,2 %. "
        "No es validación: usar como hipótesis priorizada.",
    ]
    return ModuleResult("sigma", True, outputs=outputs, metrics=metrics,
                        warnings=warnings)


def _selftest() -> str:
    from scripts.sigma.sigma_annotator import selftest
    return selftest()


def _estimate(cfg) -> str:
    return "10–30 s por genoma (CPU; sin GPU)"


MODULE = Module(
    id="sigma",
    name="Sigmas y promotores",
    summary="Anota factores sigma y promotores sigma-dependientes con evidencia, "
            "tipo y confianza; cruza función del operón y expresión opcional.",
    order=35,
    category="regulación",
    version="1.0.0",
    experimental=True,
    options=(
        Option("window", "Ventana aguas arriba (nt)", "int", 180,
               help="Región escaneada antes del primer gen del operón.",
               min=80, max=400, step=20, unit="nt"),
        Option("null_windows", "Ventanas nulas", "int", 250,
               help="Fondo para el p-valor empírico; más ventanas = más estable.",
               min=50, max=2000, step=50, advanced=True),
        Option("seed", "Semilla", "int", 17, min=0, max=10**9,
               advanced=True, help="Semilla del fondo nulo."),
        Option("max_p", "p-valor máximo", "float", 0.05,
               min=0.001, max=0.5, step=0.005, advanced=True),
        Option("min_z", "z mínimo", "float", 2.0,
               min=0.0, max=10.0, step=0.5, advanced=True),
        Option("min_margin", "Margen mínimo entre sigma", "float", 1.0,
               min=0.0, max=10.0, step=0.5, advanced=True),
        Option("use_expression", "Usar expresión si existe", "bool", True,
               help="Busca data/expression/<genoma>.tsv; si no existe, no falla."),
    ),
    docs="docs/SIGMA_MODULE.md",
    input_kind="genome",
    run_label="Anotar sigmas",
    run_fn=_run,
    selftest_fn=_selftest,
    status_fn=_status,
    estimate_fn=_estimate,
)
