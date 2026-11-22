"""GUI genérica para los módulos autodescubiertos.

No hay una pestaña escrita a mano por módulo: esta vista lee el contrato de
`modules/` y dibuja estado, opciones, tiempo estimado, autoprueba y ejecución.
Añadir un módulo no requiere tocar `app.py`.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules import discover, get  # noqa: E402
from modules.base import Module, Option, RunContext  # noqa: E402

from i18n import tr  # noqa: E402


@st.cache_data(show_spinner=False)
def _all_status(data_root: str) -> dict:
    """Calcula el estado una sola vez (importar torch/pyrodigal es lento)."""
    root = Path(data_root)
    return {m.id: m.status(root) for m in discover()}


def _status(module: Module, data_root: Path):
    try:
        return module.status(data_root)
    except Exception as exc:  # noqa: BLE001
        return type("S", (), {"available": False, "badge": lambda self: "⚠️",
                             "reasons": (str(exc),), "detail": {}})()


def _option_value(module: Module, opt: Option, values: dict[str, Any]) -> Any:
    key = f"mod_{module.id}_{opt.key}"
    label = tr(opt.label)
    help_ = tr(opt.help) if opt.help else None
    if opt.depends_on and not values.get(opt.depends_on):
        return opt.default
    if opt.kind == "bool":
        return st.checkbox(label, value=bool(opt.default), help=help_, key=key)
    if opt.kind == "int":
        return st.number_input(label, min_value=int(opt.min or 0),
                               max_value=int(opt.max or 10_000),
                               value=int(opt.default),
                               step=int(opt.step or 1), help=help_, key=key)
    if opt.kind == "float":
        return st.number_input(label, min_value=float(opt.min or 0.0),
                               max_value=float(opt.max or 1e9),
                               value=float(opt.default),
                               step=float(opt.step or 0.1), help=help_, key=key)
    if opt.kind == "select":
        idx = list(opt.choices).index(opt.default) if opt.default in opt.choices else 0
        return st.selectbox(label, opt.choices, index=idx, help=help_, key=key)
    if opt.kind == "multiselect":
        return st.multiselect(label, opt.choices, default=list(opt.default),
                              help=help_, key=key)
    return st.text_input(label, value=str(opt.default), help=help_, key=key)


def _render_options(module: Module) -> dict[str, Any]:
    values: dict[str, Any] = {}
    basic = [o for o in module.options if not o.advanced]
    advanced = [o for o in module.options if o.advanced]
    changed = False
    for opt in basic:
        val = _option_value(module, opt, values)
        values[opt.key] = val
        if val != opt.default:
            changed = True
    if advanced:
        with st.expander(tr("⚙️ Opciones avanzadas")):
            for opt in advanced:
                val = _option_value(module, opt, values)
                values[opt.key] = val
                if val != opt.default:
                    changed = True
    try:
        cfg = module.config(values)
    except ValueError as exc:
        st.error(str(exc))
        return values
    if changed:
        st.caption(tr("Configuración no predeterminada; se usará al ejecutar."))
    return cfg


def _store(module: Module, results: dict[str, Any]) -> None:
    """Guarda resultados para que los vea otro módulo (p. ej. metabolismo)."""
    st.session_state[f"module_results_{module.id}"] = results


def _render_result(res) -> None:
    if res.ok:
        st.success(tr("Módulo ejecutado en {s} s.").format(s=res.elapsed))
        if res.metrics:
            st.json(res.metrics, expanded=False)
        for name, path in res.outputs.items():
            p = Path(path)
            if p.exists() and p.stat().st_size < 5_000_000:
                st.download_button(f"⬇️ {name}", p.read_bytes(), file_name=p.name,
                                   key=f"dl_{res.module}_{name}")
        for w in res.warnings:
            st.warning(w)
    else:
        st.error(res.error or tr("El módulo falló sin detalle."))


def _run_inputs(module: Module, data_root: Path):
    if module.input_kind == "genome":
        genomes = []
        reg = data_root / "regulation"
        if reg.exists():
            genomes = sorted(d.name for d in reg.iterdir()
                             if d.is_dir() and (d / "genes.tsv").exists())
        if (ROOT / "webapp" / "demo" / "genes.tsv").exists():
            genomes = ["demo"] + [g for g in genomes if g != "demo"]
        if not genomes:
            st.info(tr("No hay genomas con datos para resumir."))
            return None
        return st.selectbox(tr("Genoma"), genomes,
                            key=f"mod_genome_{module.id}")
    if module.id == "annotation":
        uploaded = st.file_uploader(
            tr("FASTA para anotar"), type=["fa", "fna", "fasta", "fas", "gz"],
            accept_multiple_files=True, key="mod_annot_files",
        )
        use_demo = st.checkbox(tr("Usar el genoma de demostración"),
                               value=not uploaded, key="mod_annot_demo")
        inputs = []
        if use_demo:
            demo = ROOT / "webapp" / "samples" / "demo.fna"
            if demo.exists():
                inputs.append(("demo.fna", demo.read_bytes()))
        for f in uploaded or []:
            inputs.append((f.name, f.read()))
        return inputs
    if module.id == "metabolism":
        last = (st.session_state.get("module_results_annotation") or [None])[-1]
        if last:
            return last.get("proteins")
        st.info(tr("Ejecuta primero el módulo de Anotación en esta sesión."))
        return None
    return None


def render(data_root: Path, out_root: Path) -> None:
    st.subheader("🧩 " + tr("Módulos"))
    st.caption(tr(
        "Cada módulo declara sus opciones, requisitos, tiempo y autoprueba. "
        "La interfaz se construye sola: añadir una carpeta en `modules/` la "
        "hace aparecer aquí, sin editar la aplicación."))

    mods = discover()
    if not mods:
        st.warning(tr("No se encontró ningún módulo en `modules/`."))
        return

    statuses = _all_status(str(data_root))
    st.markdown("#### " + tr("Estado de todos los módulos"))
    rows = []
    for m in mods:
        st_m = statuses[m.id]
        rows.append({
            tr("Módulo"): tr(m.name),
            tr("Estado"): st_m.badge(),
            tr("Opciones"): len(m.options),
            tr("Hardware"): "GPU" if m.gpu else "CPU",
            tr("Tiempo estimado:"): tr(m.estimate(m.defaults())),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")

    labels = {m.id: f"{statuses[m.id].badge()} {tr(m.name)}" for m in mods}
    current = st.session_state.get("module_selected")
    if current not in labels:
        current = mods[0].id
    sel = st.selectbox(tr("Módulo"), list(labels), index=list(labels).index(current),
                       format_func=lambda i: labels[i], key="module_selected")

    module = get(sel)
    if module is None:
        st.error(tr("Módulo no encontrado."))
        return
    status = statuses.get(module.id) or _status(module, data_root)

    st.markdown(f"### {tr(module.name)}")
    st.markdown(tr(module.summary))
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(tr("Estado"), "✅ " + tr("disponible") if status.available
              else "⚠️ " + tr("incompleto"))
    c2.metric(tr("Hardware"), "GPU" if module.gpu else "CPU")
    c3.metric(tr("Opciones"), len(module.options))
    c4.metric(tr("Versión"), module.version)
    st.caption(f"⏱️ {tr('Tiempo estimado:')} {tr(module.estimate(module.defaults()))}")
    if status.reasons:
        for reason in status.reasons:
            st.warning(tr(reason))
    if module.experimental:
        st.warning(tr("Módulo experimental: no forma parte del núcleo estable."))

    cfg = _render_options(module)

    a, b = st.columns(2)
    with a:
        if st.button("🧪 " + tr("Probar módulo"), key=f"test_{module.id}"):
            with st.spinner(tr("Ejecutando autoprueba…")):
                try:
                    st.success(module.selftest())
                except Exception as exc:  # noqa: BLE001
                    st.error(str(exc))
    with b:
        st.caption(tr("La autoprueba no descarga datos ni modifica el proyecto."))

    st.divider()
    if module.run_fn is None:
        st.info(tr(
            "Este módulo ya declara sus opciones y su autoprueba, pero su "
            "ejecución desde esta vista todavía no está conectada. Puedes "
            "explorarlo en su pestaña correspondiente."))
    else:
        inputs = _run_inputs(module, data_root)
        disabled = not status.available or (module.input_kind is not None and inputs is None)
        if st.button("▶️ " + tr(module.run_label), type="primary",
                     disabled=disabled, key=f"run_{module.id}"):
            ctx = RunContext(
                data_root=data_root, out_root=out_root, config=cfg,
                log=lambda msg: st.write(msg),
                results={k: v for k, v in (
                    (m.id, st.session_state.get(f"module_results_{m.id}"))
                    for m in mods) if v},
            )
            with st.spinner(tr("Ejecutando {m}…").format(m=tr(module.name))):
                res = module.run(ctx, inputs)
            _store(module, ctx.results.get(module.id, []))
            _render_result(res)
            if module.id == "annotation" and res.ok:
                st.session_state["annotation_results"] = ctx.results.get("annotation", [])

    with st.expander(tr("Contrato del módulo")):
        st.markdown(tr(
            "Un módulo es una carpeta `modules/<id>/module.py` que exporta "
            "`MODULE` con opciones, estado, tiempo, autoprueba y, si aplica, "
            "una función de ejecución. Ver `docs/MODULES.md`."))
        st.code(f"python -m modules --describe {module.id}\n"
                f"python -m modules --selftest {module.id}", language="bash")
