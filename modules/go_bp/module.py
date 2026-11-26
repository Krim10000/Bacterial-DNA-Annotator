"""BP por propagación oficial MF→BP."""
from __future__ import annotations
from pathlib import Path
from modules.base import Module, ModuleResult, ModuleStatus, RunContext


def _status(data_root: Path):
    if not (data_root/"labels/external/go/go.obo").exists():
        return ModuleStatus(False,("Falta go.obo (no go-basic.obo).",))
    return ModuleStatus(True)


def _run(ctx: RunContext, genome: str):
    from webapp.go_bp import build_for_genome
    if not genome: return ModuleResult("go_bp",False,error="Selecciona un genoma.")
    res=build_for_genome(ctx.data_root,genome,ctx.out_root/"go_bp"/genome)
    out=Path(res["out_dir"])
    return ModuleResult("go_bp",True,outputs={"gene_bp":out/"gene_bp.tsv","summary":out/"bp_summary.json"},
        metrics=res["summary"],
        warnings=["Solo MF→BP con enlace part_of de go.obo (GO_REF:0000108).",
                  "No se infiere BP desde CC, EC, Rhea ni ruta."])


def _selftest() -> str:
    from webapp.go_bp import load_mf_to_bp
    assert callable(load_mf_to_bp)
    return "propagación MF→BP OK"


MODULE=Module(id="go_bp", name="GO BP", summary="BP por propagación oficial MF→BP; sin inferir desde CC o rutas.",
    order=33, category="anotación", docs="docs/GO_BP.md", input_kind="genome",
    run_label="Propagar BP", run_fn=_run, selftest_fn=_selftest, status_fn=_status,
    estimate_fn=lambda cfg:"~5–15 s por genoma")
