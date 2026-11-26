"""GO para rRNA/tRNA/ncRNA desde Infernal/GFF."""
from __future__ import annotations
from pathlib import Path
from modules.base import Module, ModuleResult, ModuleStatus, RunContext


def _status(data_root: Path):
    if not (data_root/"annotations").exists() and not (data_root/"regulation").exists():
        return ModuleStatus(False,("Faltan anotaciones o genes.",))
    return ModuleStatus(True)


def _run(ctx: RunContext, genome: str):
    from webapp.ncrna_go import build_for_genome
    if not genome:
        return ModuleResult("ncrna_go",False,error="Selecciona un genoma.")
    res=build_for_genome(ctx.data_root,genome,ctx.out_root/"ncrna_go"/genome)
    out=Path(res["out_dir"])
    return ModuleResult("ncrna_go",True,outputs={"ncrna_go":out/"ncrna_go.tsv","summary":out/"ncrna_summary.json"},metrics=res["summary"],
        warnings=["Anotación por tipo de feature/Rfam, no experimental.",
                  "ncRNA sin Rfam conocido queda fuera; no se inventa."])


def _selftest() -> str:
    from webapp.ncrna_go import build_for_genome
    assert callable(build_for_genome)
    return "motor ncRNA GO OK"


MODULE=Module(id="ncrna_go", name="GO ARN", summary="GO MF/CC/BP para rRNA, tRNA y ncRNA desde Infernal/GFF.",
    order=29, category="anotación", docs="docs/NCRNA_GO.md", input_kind="genome",
    run_label="Anotar GO ARN", run_fn=_run, selftest_fn=_selftest, status_fn=_status,
    estimate_fn=lambda cfg:"~2–10 s por genoma (GFF ya calculado)")
