"""Regulación transcripcional: operones, promotores, horquillas y fronteras."""
from __future__ import annotations

import json
from pathlib import Path

from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path):
    demo = Path(__file__).resolve().parents[2] / "webapp" / "demo" / "genes.tsv"
    if (data_root / "regulation" / "ecoli_k12" / "genes.tsv").exists() or demo.exists():
        return ModuleStatus(True)
    return ModuleStatus(False, ("Faltan datos de regulación; la demo sí está disponible.",
                                "Genera con `bash run_pipeline.sh components`."))


def _run(ctx: RunContext, genome: str) -> ModuleResult:
    from browser import load_tables, load_transcription

    tables = load_tables(ctx.data_root, genome)
    genes = tables.get("genes")
    if genes is None or genes.empty:
        return ModuleResult("regulation", False,
                            error=f"No hay genes para {genome}.")
    tss, tts = load_transcription(ctx.data_root, genome)
    summary = {
        "genome": genome,
        "genes": len(genes),
        "operons": len(tables.get("operons", [])),
        "promoters": len(tables.get("promoters", [])),
        "hairpins": len(tables.get("hairpins", [])),
        "tss": len(tss),
        "tts": len(tts),
    }
    out = ctx.out_root / "regulation"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{genome}_summary.json"
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return ModuleResult("regulation", True, outputs={"summary": path},
                        metrics=summary)


def _selftest() -> str:
    from browser import svg
    from operons import load_rule_operons
    assert callable(svg) and callable(load_rule_operons)
    return "navegador y lector de operones OK"


MODULE = Module(
    id="regulation",
    name="Regulación y operones",
    summary="Agrupa genes en operones, anota promotores y cruza TSS/TTS.",
    order=30,
    category="regulación",
    options=(
        Option("rules", "Operones por reglas deterministas", "bool", True,
               help="Fronteras por TSS, promotor y horquilla; no usa modelo."),
        Option("hairpins", "Horquillas de terminación", "bool", True,
               help="Pliegue del ARN y ΔG con ViennaRNA."),
        Option("suspects", "Sospechas por cadena metabólica", "bool", True,
               help="Genes contiguos que comparten metabolitos."),
        Option("tss_tts", "Fronteras TSS/TTS", "bool", True,
               help="Curación y alto rendimiento; el dato manda sobre la regla."),
    ),
    docs="docs/L2_INTERFACE.md",
    input_kind="genome",
    run_label="Resumir regulación",
    run_fn=_run,
    selftest_fn=_selftest,
    status_fn=_status,
    estimate_fn=lambda cfg: "~8 s por genoma (reglas) + ~1–10 min (horquillas)",
)
