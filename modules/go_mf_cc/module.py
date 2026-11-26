"""GO MF/CC por accesión UniProt."""
from __future__ import annotations
from pathlib import Path
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    if not (data_root/"regulation").exists():
        return ModuleStatus(False,("Faltan genes en data/regulation/.",))
    return ModuleStatus(True)


def _run(ctx: RunContext, genome: str):
    from webapp.go_mf_cc import build_for_genome
    if not genome:
        return ModuleResult("go_mf_cc",False,error="Selecciona un genoma.")
    res=build_for_genome(ctx.data_root,genome,ctx.out_root/"go_mf_cc"/genome)
    out=Path(res["out_dir"])
    return ModuleResult("go_mf_cc",True,outputs={
        "gene_go":out/"gene_go.tsv","summary":out/"go_summary.json"},
        metrics=res["summary"],
        warnings=["MF/CC por accesión UniProt; no es ortología ni prueba experimental.",
                  "Si no hay accesión o GO, el gen queda no_data."])


def _selftest() -> str:
    from webapp.go_mf_cc import load_index
    assert callable(load_index)
    return "índice GO MF/CC OK"


MODULE=Module(
    id="go_mf_cc", name="GO MF/CC", summary="Anota función molecular y componente celular por accesión UniProt.",
    order=31, category="anotación",
    docs="docs/GO_MF_CC.md", input_kind="genome", run_label="Anotar GO MF/CC",
    run_fn=_run, selftest_fn=_selftest, status_fn=_status,
    estimate_fn=lambda cfg:"~5 s por genoma con índice local",
)
