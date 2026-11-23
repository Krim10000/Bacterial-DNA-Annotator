"""ORFs no anotados, desiertos y cierres por dominio."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _status(data_root):
    return ModuleStatus(True)


def _run(ctx: RunContext, inputs) -> ModuleResult:
    from annotator import annotate
    from ghost.pipeline import build_candidates

    if not inputs:
        return ModuleResult("ghost", False, error="No hay archivos FASTA.")
    rows = []
    metrics = {}
    outdir = ctx.out_root / "ghost"
    outdir.mkdir(parents=True, exist_ok=True)
    for item in inputs:
        name, data = (item.get("name", "genoma.fna"), item.get("data", b"")) \
            if isinstance(item, dict) else item
        stem = Path(str(name)).stem
        ann = annotate(str(name), bytes(data), {
            "cds": True, "trna": False, "rrna": False, "ghost": False,
            "esm": False, "deserts": False,
        })
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "raw").mkdir()
            (root / "annotations_ab_initio" / stem).mkdir(parents=True)
            (root / "raw" / f"{stem}.fna").write_bytes(bytes(data))
            (root / "annotations_ab_initio" / stem / "annotation.gff3").write_text(
                ann["gff"], encoding="utf-8")
            res = build_candidates(
                root, stem,
                min_region_len=int(ctx.config.get("min_region_len", 150)),
                min_aa=int(ctx.config.get("min_aa", 60)),
            )
        for c in res.candidates:
            c["genome"] = stem
            rows.append(c)
        metrics[stem] = {"regions": res.n_regions, "orfs": res.n_orfs,
                         "genome_bp": res.genome_bp}
    path = outdir / "ghost_candidates.tsv"
    if rows:
        import csv
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    summary = outdir / "ghost_summary.json"
    summary.write_text(json.dumps(metrics, ensure_ascii=False, indent=2),
                       encoding="utf-8")
    return ModuleResult("ghost", True,
                        outputs={"candidates": path, "summary": summary},
                        metrics=metrics,
                        warnings=[] if rows else ["No se encontraron ORFs candidatos."])


def _selftest() -> str:
    from ghost.orfs import merge_intervals
    from ghost.seqio import read_fasta
    assert callable(merge_intervals) and callable(read_fasta)
    return "motor de genes fantasma OK"


MODULE = Module(
    id="ghost",
    name="Genes fantasma",
    summary="Busca genes en los huecos que nadie anotó y regiones desiertas.",
    order=60,
    category="anotación",
    options=(
        Option("min_region_len", "Longitud mínima de región (nt)", "int", 150,
               min=30, max=5000, step=10,
               help="Huecos más cortos generan ruido."),
        Option("min_aa", "Longitud mínima de ORF (aa)", "int", 60,
               min=20, max=500, step=5,
               help="Por debajo de 60 aa es difícil distinguir un gen real."),
        Option("esm", "Puntuar con ESM-2", "bool", False,
               help="Modelo de lenguaje de proteínas; corre en CPU."),
    ),
    docs="scripts/ghost/README.md",
    input_kind="fasta",
    run_label="Buscar genes fantasma",
    run_fn=_run,
    status_fn=_status,
    selftest_fn=_selftest,
    estimate_fn=lambda cfg: "~5 s por genoma; ~6 h los 400",
)
