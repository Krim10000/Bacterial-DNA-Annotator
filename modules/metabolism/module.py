"""Red metabólica, reacciones, brechas y cierres."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from modules.base import Module, ModuleResult, ModuleStatus, Option, RunContext


def _has_data(data_root: Path) -> tuple[bool, list[str]]:
    need = [
        data_root / "metabolism" / "rhea2uniprot_sprot.tsv",
        data_root / "metabolism" / "rhea_reactions.pkl",
        data_root / "db" / "swissprot.dmnd",
    ]
    missing = [str(p.relative_to(data_root)) for p in need if not p.exists()]
    return (not missing, missing)


def _status(data_root: Path) -> ModuleStatus:
    ok, missing = _has_data(data_root)
    if not ok:
        return ModuleStatus(False, (
            "Faltan recursos: " + ", ".join(missing),
            "Ejecuta `python -m webapp.setup_data --data-root data`.",
        ))
    import shutil
    if shutil.which("diamond") is None:
        return ModuleStatus(False, ("Falta el binario DIAMOND.",))
    return ModuleStatus(True)


def _estimate(cfg: dict[str, Any]) -> str:
    return "~2 min por genoma (más si hay muchos ORFs)"


def _run(ctx: RunContext, inputs: Any) -> ModuleResult:
    from metabolism import build_network

    proteins = inputs
    if proteins is None:
        last = (ctx.results.get("annotation") or [None])[-1]
        proteins = (last or {}).get("proteins")
    if not proteins:
        return ModuleResult("metabolism", False,
                            error="No hay proteínas: ejecuta antes Anotación.")
    cfg = ctx.config
    data = ctx.data_root
    panel = data / "db" / "pfam" / "panel_focus.hmm" if cfg.get("use_panel") else None
    transport = data / "db" / "transport.dmnd" if cfg.get("use_transport") else None
    producer = data / "db" / "producer_family.dmnd" if cfg.get("use_producer") else None
    net = build_network(
        proteins, data, cfg.get("reference_gcf") or None,
        threads=int(cfg.get("threads", 4)),
        panel_hmm=panel if panel and panel.exists() else None,
        transport_db=transport if transport and transport.exists() else None,
        producer_db=producer if producer and producer.exists() else None,
    )
    out = ctx.out_root / "metabolism"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "network.json"
    path.write_text(json.dumps(net, ensure_ascii=False, indent=2, default=str),
                    encoding="utf-8")
    return ModuleResult(
        "metabolism", True, outputs={"network": path},
        metrics={"reactions": net.get("reactions", 0),
                 "gaps": net.get("n_gaps", 0)},
        warnings=[] if net.get("n_gaps") is not None else ["La red no reportó brechas."],
    )


def _selftest() -> str:
    from metabolism import build_network
    assert callable(build_network)
    return "import de build_network OK (los datos se comprueban en status)"


MODULE = Module(
    id="metabolism",
    name="Red metabólica y brechas",
    summary="Anota proteínas contra Swiss-Prot, las conecta con Rhea y busca "
            "los pasos que faltan.",
    order=20,
    category="metabolismo",
    requires=("diamond",),
    options=(
        Option("threads", "Hilos", "int", 4, min=1, max=64, step=1,
               help="Paralelismo de DIAMOND."),
        Option("reference_gcf", "GCF de referencia (opcional)", "text", "",
               help="Si el genoma ya está en el proyecto, reutiliza su red."),
        Option("use_panel", "Panel HMM de familias", "bool", True,
               help="Confirma familias con perfiles Pfam además de DIAMOND."),
        Option("use_transport", "Base de transportadores", "bool", True,
               help="Distingue brecha de transporte/importación."),
        Option("use_producer", "Base de productores", "bool", True,
               help="Busca el productor de la reacción en homólogos."),
    ),
    docs="docs/PIPELINE_LOGIC.md",
    run_fn=_run,
    selftest_fn=_selftest,
    estimate_fn=_estimate,
    status_fn=_status,
)
