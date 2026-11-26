"""Etiquetas de proceso anidadas para genes y operones."""
from __future__ import annotations
from pathlib import Path
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path) -> ModuleStatus:
    reg=data_root/"regulation"
    if not reg.exists() or not any((d/"genes.tsv").exists() for d in reg.iterdir() if d.is_dir()):
        return ModuleStatus(False,("Faltan genes en data/regulation/.",))
    return ModuleStatus(True, detail={"external":str(data_root/"labels"/"external")})


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    from webapp.labels_engine import build
    if not genome:
        return ModuleResult("labels",False,error="Selecciona un genoma.")
    res=build(ctx.data_root,genome,ctx.out_root/"labels"/genome,
              use_lexicon=bool(ctx.config.get("use_lexicon",True)),
              use_go=bool(ctx.config.get("use_go",True)),
              use_kegg=bool(ctx.config.get("use_kegg",True)),
              allow_iea=bool(ctx.config.get("allow_iea",False)))
    out=Path(res["out_dir"])
    return ModuleResult("labels",True,outputs={
        "gene_labels":out/"gene_labels.tsv",
        "operon_labels":out/"operon_labels.tsv",
        "summary":out/"label_summary.json"},
        metrics=res["summary"],
        warnings=["Las etiquetas son evidencia contextual; no prueban función ni sigma.",
                  "GO/KEGG solo si están descargados; si faltan, se declara."])


def _selftest() -> str:
    from webapp.labels import labels_for_text
    assert "heat_shock" in labels_for_text("DnaK heat shock protein")
    return "léxico y motor de etiquetas OK"


MODULE=Module(
    id="labels",
    name="Etiquetas de proceso",
    summary="Cruce gen/operón con GO, KEGG, EC y léxico; insumo para sigma.",
    order=32,
    category="anotación",
    options=(
        Option("use_lexicon","Léxico de procesos","bool",True),
        Option("use_go","Usar GO/GOA","bool",True),
        Option("use_kegg","Usar KEGG gene→pathway","bool",True),
        Option("allow_iea","Incluir evidencia GO IEA","bool",False,
               help="Más cobertura, más ruido; se declara en el resumen."),
    ),
    docs="docs/LABELS_MODULE.md",
    input_kind="genome",
    run_label="Construir etiquetas",
    run_fn=_run, selftest_fn=_selftest, status_fn=_status,
    estimate_fn=lambda cfg:"~5–20 s por genoma (según catálogos)",
)
