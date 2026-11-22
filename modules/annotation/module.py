"""Anotación desde cero: CDS, tRNA, rRNA, ORFs no anotados y desiertos."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from modules._bridge import ROOT
from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root: Path) -> ModuleStatus:
    import importlib.util
    import shutil

    tools = {
        "pyrodigal": importlib.util.find_spec("pyrodigal") is not None,
        "trna": shutil.which("tRNAscan-SE") is not None,
        "rrna": shutil.which("barrnap") is not None or shutil.which("cmsearch") is not None,
        "diamond": shutil.which("diamond") is not None,
        "esm": importlib.util.find_spec("torch") is not None,
    }
    if not tools.get("pyrodigal"):
        return ModuleStatus(False, ("Falta pyrodigal (pip install pyrodigal).",),
                            {"tools": tools})
    reasons = []
    if not tools.get("trna"):
        reasons.append("tRNAscan-SE no está instalado: tRNA desactivado.")
    if not tools.get("rrna"):
        reasons.append("barrnap/cmsearch no está instalado: rRNA desactivado.")
    if not tools.get("diamond"):
        reasons.append("DIAMOND no está instalado: comparación desactivada.")
    return ModuleStatus(True, tuple(reasons), {"tools": tools})


def _estimate(cfg: dict[str, Any]) -> str:
    base = "1–3 min por genoma de ~5 Mb"
    if cfg.get("ghost"):
        base += " + 5 s de ORFs"
    if cfg.get("esm"):
        base += " + 0,3 s por ORF puntuado"
    if cfg.get("metabolic"):
        base += " + ~2 min de red metabólica"
    return base


def _safe(name: str) -> str:
    stem = Path(name).name
    for ext in (".gz", ".fna", ".fa", ".fasta", ".fas"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)[:80] or "genoma"


def _run(ctx: RunContext, inputs: Any) -> ModuleResult:
    from annotator import annotate

    if not inputs:
        return ModuleResult("annotation", False, error="No hay archivos FASTA.")
    opts = {
        "cds": ctx.config.get("cds", True),
        "trna": ctx.config.get("trna", True),
        "rrna": ctx.config.get("rrna", True),
        "ghost": ctx.config.get("ghost", True),
        "esm": ctx.config.get("esm", False),
        "esm_top": int(ctx.config.get("esm_top", 100)),
        "deserts": ctx.config.get("deserts", True),
        "min_region_len": int(ctx.config.get("min_region_len", 150)),
        "min_aa": int(ctx.config.get("min_aa", 60)),
    }
    outdir = ctx.out_root / "annotation"
    outdir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    metrics: dict[str, Any] = {}
    warnings: list[str] = []
    results = []
    for item in inputs:
        if isinstance(item, dict):
            name, data = item.get("name", "genoma.fna"), item.get("data", b"")
        else:
            name, data = item
        ctx.log(f"Anotando {name}…")
        res = annotate(str(name), bytes(data), opts)
        stem = _safe(str(name))
        (outdir / f"{stem}.gff3").write_text(res["gff"], encoding="utf-8")
        (outdir / f"{stem}.proteins.faa").write_text(
            res.get("proteins_fasta", ""), encoding="utf-8")
        (outdir / f"{stem}.summary.json").write_text(
            json.dumps({k: v for k, v in res.items()
                        if k not in {"gff", "proteins", "proteins_fasta",
                                     "candidates", "features", "deserts"}},
                       ensure_ascii=False, indent=2, default=str),
            encoding="utf-8")
        (outdir / f"{stem}.report.md").write_text(
            res.get("report_md", ""), encoding="utf-8")
        outputs[stem] = outdir / f"{stem}.gff3"
        metrics[stem] = {
            "bp": res["bp"], "gc": round(res["gc"], 4),
            "cds": res["stats"].get("cds", 0),
            "coverage": round(res.get("coverage", 0), 4),
            "candidates": len(res.get("candidates", [])),
            "deserts": len(res.get("deserts", [])),
        }
        results.append(res)
        if not res["stats"].get("cds"):
            warnings.append(f"{name}: no se encontraron CDS.")
    ctx.results["annotation"] = results
    if ctx.config.get("compare") and len(results) >= 2:
        from annotator import compare_results
        comp = compare_results(results)
        ctx.results["annotation_comparison"] = comp
        metrics["comparison"] = {k: comp.get(k) for k in
                                 ("clusters", "multi_genome") if k in comp}
        path = outdir / "comparison.json"
        path.write_text(json.dumps(comp, ensure_ascii=False, indent=2,
                                   default=str), encoding="utf-8")
        outputs["comparison"] = path
    return ModuleResult("annotation", True, outputs=outputs, metrics=metrics,
                        warnings=warnings)


def _selftest() -> str:
    from annotator import annotate

    sample = ROOT / "webapp" / "samples" / "demo.fna"
    res = annotate("demo.fna", sample.read_bytes(), {
        "cds": True, "trna": False, "rrna": False, "ghost": False,
        "esm": False, "deserts": False, "min_region_len": 150, "min_aa": 60,
    })
    if res["stats"].get("cds", 0) < 50:
        raise AssertionError("la demo produjo muy pocos CDS")
    return f"{res['stats']['cds']} CDS, cobertura {res['coverage']*100:.0f}%"


MODULE = Module(
    id="annotation",
    name="Anotación de FASTA",
    summary="Encuentra CDS, tRNA, rRNA, ORFs no anotados y desiertos en un "
            "genoma nuevo.",
    order=10,
    category="anotación",
    requires=("pyrodigal",),
    options=(
        Option("cds", "CDS con pyrodigal", "bool", True,
               help="Genes que producen proteína. Sin esto no hay análisis posterior."),
        Option("trna", "tRNA con tRNAscan-SE", "bool", True,
               help="ARN de transferencia; requiere el binario tRNAscan-SE."),
        Option("rrna", "rRNA con barrnap", "bool", True,
               help="ARN ribosómico; requiere barrnap o cmsearch."),
        Option("ghost", "ORFs en regiones no anotadas", "bool", True,
               help="Busca posibles genes en los huecos entre genes conocidos."),
        Option("esm", "Puntuar con ESM-2", "bool", False,
               help="Modelo de lenguaje de proteínas; corre en CPU pero tarda."),
        Option("esm_top", "ORFs a puntuar con ESM", "int", 100,
               min=0, max=2000, step=10, depends_on="esm",
               help="Puntuar todos tarda mucho; los mejores suelen estar al principio."),
        Option("deserts", "Desiertos de anotación", "bool", True,
               help="Regiones de ≥3 kb sin señales de genes."),
        Option("compare", "Comparar genomas subidos", "bool", True,
               help="Agrupa genes parecidos entre varios FASTA."),
        Option("min_region_len", "Longitud mínima de región (nt)", "int", 150,
               min=30, max=5000, step=10, advanced=True,
               help="Huecos más cortos generan ruido."),
        Option("min_aa", "Longitud mínima de ORF (aa)", "int", 60,
               min=20, max=500, step=5, advanced=True,
               help="Por debajo de 60 aa es difícil distinguir un gen real."),
    ),
    docs="README.md#anotar-fasta",
    run_fn=_run,
    selftest_fn=_selftest,
    estimate_fn=_estimate,
    status_fn=_status,
)
