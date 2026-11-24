"""Interfaz web (Streamlit) del anotador de ADN bacteriano.

    streamlit run webapp/app.py

Sube uno o varios FASTA, marca las anotaciones y descarga GFF3/TSV/FASTA.
El flujo es automático: el usuario decide qué analizar, no hay pasos
interactivos durante la ejecución.
"""
from __future__ import annotations

import io
import json
import math
import os
import struct
import subprocess
import sys
import time
import wave
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from i18n import (  # noqa: E402
    LANGUAGES,
    init_language,
    install_streamlit_i18n,
    language_selector,
    tr,
    tr_tree,
)
from metabolism import build_network  # noqa: E402
from requirements_ui import render_requirements  # noqa: E402
from active_sites import for_gap as active_sites_for_gap  # noqa: E402
from hypotheses_ui import render_hypotheses  # noqa: E402
from precomputed_ui import render as render_precomputed  # noqa: E402
import structure_viewer
from compare_ui import render_compare
from esmfold_ui import (  # noqa: E402
    render_esmfold_closure,
    render_structure_evidence,
    render_structural_annotation,
    render_five_filters,
)
from annotator import (  # noqa: E402
    annotate,
    build_report,
    compare_results,
    deserts_tsv,
    detect_tools,
)
from functional import annotate as annotate_functional  # noqa: E402
from regulation_layer import run as run_regulation  # noqa: E402
from complexes import infer as infer_complexes, enrich_rows as enrich_complex_rows  # noqa: E402
from complex_model import run as run_complex_model  # noqa: E402
from complex_transfer import load as load_complex_transfer  # noqa: E402
import oric as oric_engine  # noqa: E402
from closures import (  # noqa: E402
    LABEL as CLOSURE_LABEL,
    ORDER as CLOSURE_ORDER,
    available as closures_available,
    for_genome as closures_for_genome,
    load_closures,
    summary as closures_summary,
    table as closures_table,
)
from operons import (  # noqa: E402
    PIECES as OPERON_PIECES,
    hairpin_genomes,
    hairpin_summary,
    hairpin_table,
    load_hairpins,
    load_rule_boundaries,
    load_rule_operons,
    load_suspects,
    rule_genomes,
    rule_summary,
    RULE_LABEL,
    available as operons_available,
    coverage as operons_coverage,
    detail_for,
    load_components,
    load_detail,
    piece_table,
)
from textos import T  # noqa: E402
import credits  # noqa: E402
import downloads_ui  # noqa: E402
import module_ui  # noqa: E402
import browser  # noqa: E402
from pathways import (  # noqa: E402
    PATHWAYS,
    catalog,
    kegg_catalog,
    kegg_catalog_from_net,
    kegg_pathway_steps,
    local_context,
    pathway_view,
)

st.set_page_config(page_title="Anotador de ADN · DNA Annotator",
                   page_icon="🧬", layout="wide")

# La fuente de iconos (Material Symbols) puede no cargar por caché del navegador
# desincronizada o por bloqueo de fuentes. Cuando eso pasa, Streamlit dibuja el
# nombre de la ligadura ("arrow_forward") como texto encima de la etiqueta del
# desplegable. Este bloque lo evita: oculta ese icono y dibuja un triángulo
# propio, que no depende de ninguna fuente externa.
st.html(
    """
    <style>
      [data-testid="stExpanderIcon"] { display: none !important; }
      [data-testid="stExpander"] summary [data-testid="stIconMaterial"] {
        display: none !important;   /* ligadura rota si la fuente no carga */
      }
      [data-testid="stExpander"] summary::before {
        content: "\\25B8";            /* ▸ */
        margin-right: .45rem; opacity: .55; font-size: 1rem;
      }
      [data-testid="stExpander"] details[open] summary::before {
        content: "\\25BE";            /* ▾ */
      }
    </style>
    """
)

# Traduce etiquetas y textos de Streamlit sin tocar cada llamada.
install_streamlit_i18n()

try:
    from webapp import __version__ as APP_VERSION
except Exception:  # ejecutado como script suelto
    APP_VERSION = "1.0.0"


PROJECT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("DNA_DATA_ROOT", str(PROJECT / "data")))
OUT_ROOT = Path(os.environ.get("DNA_OUT_ROOT", str(PROJECT / "out_fl")))
UI_CFG = Path(__file__).resolve().parents[1] / ".cache" / "ui_config.json"


def _load_ui() -> dict:
    """Preferencias de la interfaz (vista de navegación, genoma)."""
    try:
        return json.loads(UI_CFG.read_text())
    except Exception:
        return {}


def _save_ui(cfg: dict) -> None:
    """Guarda preferencias sin perder las claves de otros módulos (p. ej. idioma)."""
    try:
        actual = _load_ui()
        actual.update(cfg)
        UI_CFG.parent.mkdir(parents=True, exist_ok=True)
        tmp = UI_CFG.with_suffix(".tmp")
        tmp.write_text(json.dumps(actual, ensure_ascii=False))
        tmp.replace(UI_CFG)
    except Exception:
        pass


_genomes_file = DATA_ROOT / "metabolism" / "free_living_genomes.json"
PROJECT_GENOMES = json.load(open(_genomes_file)) if _genomes_file.exists() else []


def render_catalog() -> None:
    st.caption(T["catalog_intro"])
    st.caption(T["catalog_read"])
    genome = _global_genome()
    cats = catalog(genome, OUT_ROOT)
    if not cats:
        st.info("Este genoma no tiene brechas refinadas.")
        return
    cat = st.selectbox("Categoría", list(cats), key="cat_cat")
    gaps = cats[cat]
    names = [g["gap"] for g in gaps]
    sel = st.selectbox("Brecha", names, key="cat_gap", help=T["help_brecha"])
    gap = gaps[names.index(sel)]
    st.caption(f"{len(gaps)} brechas en «{cat}» · mostramos un contexto local acotado.")
    cards = local_context(genome, gap, DATA_ROOT, OUT_ROOT)
    icon = {"presente": "🟢", "falta": "🔴", "sin datos": "⚪"}
    cols = st.columns(len(cards))
    for col, c in zip(cols, cards):
        with col:
            st.markdown(f"### {icon[c['status']]} {c['label']}")
            st.caption(c["role"])
            with st.expander("Detalle"):
                st.write(f"**Rhea:** {c['rhea'] or '—'}")
                st.write(f"**Reacción:** {c['equation'] or 'no disponible'}")
                st.write(f"**Evidencia:** {c['evidence']}")
    st.dataframe(pd.DataFrame(cards), width="stretch", hide_index=True)
    st.caption(T["catalog_scope"])
    st.caption(credits.line("rhea", "chebi"))


def render_kegg() -> None:
    st.caption(T["kegg_intro"])
    st.caption(T["kegg_read"])
    genome = _global_genome()
    cats = kegg_catalog(genome, DATA_ROOT, OUT_ROOT)
    if not cats:
        st.info("Catálogo KEGG no disponible. Descárgalo una vez con "
                "`python -m webapp.fetch_kegg --data-root data`.")
        return
    labels = [f"{c['name']} — faltan {c['missing']} / {c['total']}" for c in cats]
    idx = st.selectbox("Ruta KEGG", range(len(cats)),
                       format_func=lambda i: labels[i], key="kegg_path")
    cat = cats[idx]
    steps = kegg_pathway_steps(genome, cat["id"], DATA_ROOT, OUT_ROOT)
    per_page = 6
    pages = max(1, (len(steps) + per_page - 1) // per_page)
    page = st.number_input("Bloque de reacciones", 1, pages, 1, key="kegg_page") if pages > 1 else 1
    visible = steps[(page - 1) * per_page: page * per_page]
    st.caption(f"{cat['name']} · {cat['present']} presentes, {cat['missing']} faltantes "
               f"de {cat['total']} reacciones mapeadas.")
    st.caption(T["kegg_pages"])
    icon = {"presente": "🟢", "falta": "🔴", "sin datos": "⚪"}
    cols = st.columns(len(visible))
    for col, stp in zip(cols, visible):
        with col:
            st.markdown(f"### {icon[stp['status']]} {stp['rhea']}")
            st.caption(stp["status"])
            with st.expander("Detalle"):
                st.write(f"**Reacción:** {stp['equation'] or 'no disponible'}")
    st.dataframe(pd.DataFrame(visible), width="stretch", hide_index=True)
    st.caption(credits.line("kegg", "rhea"))


@st.cache_data(show_spinner=False)
def gaps_cached(path: str):
    return pd.read_csv(path, sep="\t")


@st.cache_data(show_spinner=False)
def csv_cached(path: str):
    return pd.read_csv(path, sep="\t")


def render_gaps() -> None:
    st.caption(T["gaps_intro"])
    path = DATA_ROOT / "metabolism" / "gap_interpretation.tsv"
    if not path.exists():
        st.info("Sin tabla de brechas interpretadas.")
        return
    df = csv_cached(str(path))
    genomes = sorted(df["genome"].unique())
    genome = st.selectbox("Genoma", genomes, key="gap_genome", help=T["help_genoma"])
    sub = df[df["genome"] == genome].copy()
    lpath = DATA_ROOT / "metabolism" / "essential" / "lifestyle_context.tsv"
    if lpath.exists():
        ldf = csv_cached(str(lpath))
        lrow = ldf[ldf["genome"] == genome]
        if len(lrow):
            lr = lrow.iloc[0]
            st.info(f"**Contexto del genoma:** {lr['message']} "
                    f"Tamaño {lr['size_mb']} Mb · expectativa `{lr['expectation']}` · "
                    f"confianza `{lr['confidence']}`.")
    # unir evidencia de cierre
    cl_path = DATA_ROOT / "metabolism" / "closure_annotations.tsv"
    if cl_path.exists():
        cl = csv_cached(str(cl_path))
        keep = [c for c in ("genome","gap","domain","query","seqid","start","end","closure_status","direction_status","support") if c in cl.columns]
        sub = sub.merge(cl[keep], on=["genome","gap"], how="left", suffixes=("","_cl"))
    statuses = sorted(df["status"].unique())
    sel = st.multiselect("Estado", statuses, default=statuses, key="gap_status", help=T["help_estado"])
    sub = sub[sub["status"].isin(sel)]
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Brechas", len(sub), help=T["gaps_m_total"])
    c2.metric("Soporte ≥10", int((sub["support_400"] >= 10).sum()),
              help=T["gaps_m_support"])
    c3.metric("Novel", int((sub["status"] == "novel").sum()), help=T["gaps_m_novel"])
    c4.metric("Con transporte",
              int((sub.get("transport_status", "") == "transport_reaction_only").sum()),
              help=T["gaps_m_transport"])
    c5.metric("Con cierre",
              int(sub.get("domain", pd.Series(dtype=str)).notna().sum()),
              help=T["gaps_m_closure"])
    st.caption(T["filter_hint"])
    origin_labels = {
        "enzymatic_missing": "enzima no anotada / brecha real",
        "enzymatic_present": "enzima presente, mapeo incompleto",
        "degradation": "degradación",
        "false_positive": "falso positivo",
        "mapped_homology": "mapeo por homología",
        "homology_present": "homólogo presente",
        "imported_or_entry": "importado / entrada",
        "lipid_remodeling": "remodelado lipídico",
        "spontaneous": "espontáneo",
        "promiscuous": "promiscuidad enzimática",
        "reannotated": "re-anotado",
    }
    sub["cause"] = sub["origin"].map(origin_labels).fillna(sub["origin"])
    st.markdown("#### Causas del hueco")
    st.dataframe(sub["cause"].value_counts().rename_axis("causa").reset_index(name="n"),
                 width="stretch", hide_index=True)
    with st.expander("Catálogo de causas y qué significan"):
        st.markdown("""
- **Anotación incompleta**: hay ORF/dominio, pero la enzima no está nombrada.
- **Dominio divergente**: homología débil o dominio parcial.
- **Transporte/importación**: el metabolito puede entrar; no se sintetiza.
- **Espontánea**: la reacción ocurre sin enzima.
- **Multifuncional/promiscua**: otra enzima cubre la reacción.
- **Ruta alternativa**: se llega por otro camino.
- **Pérdida génica/auxotrofía**: el gen no está; depende del medio/huésped.
- **Definición KEGG/Rhea**: reacción mal asignada o incompleta.
- **Sin referencia bacteriana**: Rhea solo tiene productor eucariota; no es ausencia.
- **Falso positivo**: el hueco no es real.
""")
    cause_cols = [c for c in ("gap","cause","review_status","transport_status","homolog_status","status_reason","review_reason") if c in sub.columns]
    st.dataframe(sub[cause_cols].head(100), width="stretch", hide_index=True)
    cols = ["gap","chebi","support_400","role","origin","status","status_reason",
            "transport_status","direction_status","review_status","hmm_evidence","domain","query"]
    st.dataframe(sub[[c for c in cols if c in sub.columns]],
                 width="stretch", hide_index=True)
    # panel de evidencia
    if len(sub):
        gap_sel = st.selectbox("Ver fuentes de evidencia de…", sub["gap"].tolist(), key="gap_detail")
        r = sub[sub["gap"] == gap_sel].iloc[0]
        st.markdown(f"### {gap_sel} · `{r['status']}`")
        st.caption(f"Origen: **{r.get('origin','')}** · Rol: **{r.get('role','')}** · Confianza: **{r.get('status_confidence','')}**")
        exf = DATA_ROOT / "metabolism" / "essential" / "gap_exception_flags.tsv"
        if exf.exists():
            exd = csv_cached(str(exf))
            hit = exd[(exd["genome"] == genome) & (exd["gap"] == gap_sel)]
            if len(hit):
                st.warning("Excepción documentada: este hueco se ha observado en casos de ejemplo. "
                           "No es una exención automática; revisa la implicación y la fuente.")
                st.dataframe(hit, width="stretch", hide_index=True)
        with st.expander("Sitios activos de referencias (bajo demanda)"):
            st.caption("Precalculado si es ruta esencial; si no, consulta UniProt on-demand y cachea.")
            if st.button("Buscar sitios activos", key=f"as_{genome}_{gap_sel}"):
                with st.spinner("Consultando UniProt..."):
                    res = active_sites_for_gap(DATA_ROOT, genome, gap_sel, fetch=True)
                st.caption(f"Fuente: {res['source']} · ruta: {res['map_id'] or 'sin ruta'}")
                if res["rows"]:
                    st.dataframe(pd.DataFrame(res["rows"]), width="stretch", hide_index=True)
                else:
                    st.info("Sin sitios activos anotados en UniProt para las referencias.")
        with st.expander("1 · Trazado y soporte (Rhea + 400 genomas)", expanded=True):
            st.caption(T["gaps_e1"])
            st.write(f"**Productores Rhea:** {r.get('producers','') or 'n/d'}")
            st.write(f"**Soporte:** {r.get('support_400','')}/400 genomas con productor")
            st.caption(T["gaps_support_help"])
            st.write(f"**Motivo:** {r.get('status_reason','') or 'n/d'}")
        with st.expander("2 · Transporte (Rhea + homólogo)"):
            st.caption(T["gaps_e2"])
            st.write(f"**Estado:** {r.get('transport_status','') or 'sin evidencia'}")
            st.write(f"**Reacción Rhea:** {r.get('transport_rhea','') or '—'}")
            st.write(f"**Transportador:** {r.get('transporter_acc','') or '—'} · pident {r.get('transport_pident','') or '—'} · qcov {r.get('transport_qcov','') or '—'} · e {r.get('transport_evalue','') or '—'}")
        with st.expander("3 · Dirección y termodinámica"):
            st.caption(T["gaps_e3"])
            st.write(f"**Estado:** {r.get('direction_status','') or 'sin dato'}")
            st.write(f"**Motivo:** {r.get('direction_reason','') or '—'}")
            st.write(f"**Termodinámica:** {r.get('thermo_status','') or 'pendiente'}")
            st.markdown("Fuentes: [Rhea: lados y direcciones](https://www.rhea-db.org/help/reaction%2Dside%2Ddirection) · "
                        "[eQuilibrator 3.0](https://doi.org/10.1093/nar/gkab1106) · [TECRDB](https://doi.org/10.1093/bioinformatics/bth314)")
        with st.expander("4 · Revisión manual"):
            st.caption(T["gaps_e4"])
            st.write(f"**Estado:** {r.get('review_status','') or 'sin revisar'}")
            st.write(f"**Motivo:** {r.get('review_reason','') or '—'}")
        with st.expander("5 · Cierre candidato"):
            st.caption(T["gaps_e5"])
            if str(r.get('domain','')) not in ('','nan'):
                st.write(f"**Dominio:** {r.get('domain')}")
                st.write(f"**ORF:** {r.get('query')}")
                st.write(f"**Locus:** {r.get('seqid')}:{r.get('start')}-{r.get('end')}")
                st.write(f"**Estado:** {r.get('closure_status') or r.get('closure_status_cl','candidate')}")
            else:
                st.write("Sin cierre por dominio.")
            render_esmfold_closure(r, None, key_prefix=f"gap_{genome}_{gap_sel}")
        with st.expander("6 · HMM"):
            st.caption(T["gaps_e6"])
            st.write(r.get('hmm_evidence','') or "Sin hit HMM.")
    rpath = DATA_ROOT / "metabolism" / "review" / "gap_review_current.tsv"
    if not rpath.exists():
        rpath = DATA_ROOT / "metabolism" / "review" / "gap_review.tsv"
    if rpath.exists():
        rdf = csv_cached(str(rpath))
        st.markdown("#### Clasificación por gap (capa de revisión)")
        cls = sorted(rdf["classification"].unique())
        csel = st.multiselect("Clase", cls, default=cls, key="gap_review_class", help=T["help_clase"])
        rview = rdf[rdf["classification"].isin(csel)]
        st.dataframe(rview.head(200), width="stretch", hide_index=True)
        st.download_button("⬇️ Revisión de gaps (TSV)", rview.to_csv(sep="\t", index=False),
                           file_name="gap_review.tsv", mime="text/plain")
        st.caption(T["gaps_review"])
    epath = DATA_ROOT / "metabolism" / "evidence" / "gap_evidence_all.tsv"
    if epath.exists():
        edf = csv_cached(str(epath))
        st.markdown("#### Evidencia por gap (lotes 01+02)")
        est = sorted(edf["status"].unique())
        esel = st.multiselect("Estado evidencia", est, default=est, key="gap_ev_status")
        eview = edf[(edf["genome"] == genome) & edf["status"].isin(esel)]
        st.dataframe(eview.head(200), width="stretch", hide_index=True)
        st.download_button("⬇️ Evidencia de gaps (TSV)", eview.to_csv(sep="\t", index=False),
                           file_name="gap_evidence.tsv", mime="text/plain")
        st.caption("Homología + dominio son ejes separados; `not_checked` no es ausencia. "
                   "3.956 gaps evaluados; 1.174 filas sin referencia bacteriana (`no_bacterial_reference`), no genomas eucariotas.")
    cpath = DATA_ROOT / "metabolism" / "evidence" / "mapping_candidates_all.tsv"
    if cpath.exists():
        cdf = csv_cached(str(cpath))
        st.markdown(f"#### Candidatos de mapeo ({len(cdf)})")
        st.dataframe(cdf.head(100), width="stretch", hide_index=True)
        st.download_button("⬇️ Candidatos de mapeo (TSV)", cdf.to_csv(sep="\t", index=False),
                           file_name="mapping_candidates.tsv", mime="text/plain")
        st.caption("Sugerencias de remapeo; no modifican la red ni la tabla base.")
    qpath = DATA_ROOT / "metabolism" / "REVIEW_QUEUE.tsv"
    if qpath.exists():
        q = csv_cached(str(qpath))
        st.markdown("#### Cola de revisión priorizada")
        st.caption(T["gaps_queue"])
        st.dataframe(q.head(20), width="stretch", hide_index=True)
        st.download_button("⬇️ Cola de revisión (TSV)", q.to_csv(sep="\t", index=False),
                           file_name="REVIEW_QUEUE.tsv", mime="text/plain")
    st.download_button("⬇️ Brechas interpretadas (TSV)",
                       sub.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_gaps_interpreted.tsv", mime="text/plain")
    st.caption(credits.line("rhea", "hmmer", "uniprot", "equilibrator"))


@st.cache_data(show_spinner=False)
def origins_cached(path: str):
    return pd.read_csv(path, sep="\t")


def render_origin() -> None:
    st.caption(T["origin_intro"])
    path = DATA_ROOT / "orca_predictions_fix" / "oriC_status.tsv"
    # v2 añade la evidencia de triaje (GC skew, cajas DnaA) sin cambiar estados
    v2 = DATA_ROOT / "orca_predictions_fix" / "oriC_status_v2.tsv"
    if v2.exists():
        path = v2
    if not path.exists():
        st.info("Sin tabla de estados de origen. Ejecuta el pipeline de orígenes.")
        return
    df = origins_cached(str(path))
    tri_path = DATA_ROOT / "orca_predictions_fix" / "oriC_triage_summary.tsv"
    if tri_path.exists():
        tri = origins_cached(str(tri_path))
        st.markdown("#### Triaje global de anomalías")
        st.dataframe(tri, width="stretch", hide_index=True)
    genome = _global_genome()
    row = df[df["genome"] == genome]
    if row.empty:
        st.info("Este genoma no tiene origen calculado.")
        return
    r = row.iloc[0]
    c1, c2, c3 = st.columns(3)
    c1.metric("Estado", str(r["status"]), help=T["origin_m_state"])
    c2.metric("Posición oriC",
              f"{int(float(r['position'])):,}" if str(r["position"]) else "n/d",
              help=T["origin_m_pos"])
    try:
        c3.metric("Confianza", f"{float(r['confidence']):.2f}", help=T["origin_m_conf"])
    except Exception:
        c3.metric("Confianza", "n/d", help=T["origin_m_conf"])
    st.write(f"**Validado por:** {r['validated_by'] or 'n/d'}")
    st.caption("«Validado por» lista los métodos que coinciden en esta posición.")
    st.write(f"**Motivo:** {r['reason'] or 'n/d'}")
    st.write(f"**Intentos:** {r['attempted'] or 'n/d'}")
    st.dataframe(row, width="stretch", hide_index=True)
    gff = path.parent / "oriC_all.gff3"
    if gff.exists():
        st.download_button("⬇️ Descargar GFF de orígenes", gff.read_text(),
                           file_name="oriC_all.gff3", mime="text/plain")
    if "triage" in row.columns and str(r.get("triage", "")):
        st.write(f"**Triaje de anomalía:** {r.get('triage')} — {r.get('triage_reason', '')}")
    st.caption(T["origin_states"])
    st.caption(T["origin_triage"])
    st.caption(credits.line("orca", "doric"))


def render_curated() -> None:
    st.caption(T["curated_intro"])
    st.caption(T["curated_read"])
    genome = _global_genome()
    pathway = st.selectbox("Ruta", list(PATHWAYS), index=0, key="cur_path", help=T["help_ruta"])
    steps = pathway_view(genome, pathway, DATA_ROOT, OUT_ROOT)
    per_page = 6
    pages = max(1, (len(steps) + per_page - 1) // per_page)
    page = st.number_input("Bloque de pasos", min_value=1, max_value=pages,
                           value=1, key="cur_page") if pages > 1 else 1
    visible = steps[(page - 1) * per_page: page * per_page]
    if pages > 1:
        st.caption(f"Mostrando pasos {(page-1)*per_page+1}–{min(page*per_page, len(steps))} "
                   f"de {len(steps)} (ruta larga).")
    cols = st.columns(len(visible))
    icon = {"presente": "🟢", "falta": "🔴", "sin datos": "⚪"}
    for col, s in zip(cols, visible):
        with col:
            st.markdown(f"### {icon[s['status']]} {s['enzyme']}")
            st.caption(f"EC {s['ec']} · {s['label']}")
            with st.expander("Detalle"):
                st.write(f"**Rhea:** {s['rhea'] or 'sin mapeo'}")
                st.write(f"**Reacción:** {s['equation'] or 'no disponible'}")
                st.write(f"**Locus ADN:** {s['locus'] or 'sin candidato localizado'}")
                st.write(f"**Evidencia:** {s.get('evidence') or 'sin evidencia'}")
    st.dataframe(pd.DataFrame(steps), width="stretch", hide_index=True)
    st.caption("Rutas curadas de referencia; el catálogo cubre todas las brechas del genoma.")
    st.caption(credits.line("rhea", "kegg"))


@st.cache_data(show_spinner=False)
def closures_cached(path: str):
    return load_closures(Path(path))


@st.cache_data(show_spinner=False)
def closures_ann_cached(path: str):
    return pd.read_csv(path, sep="\t")


def render_closures() -> None:
    st.caption(T["closures_intro"])
    path = DATA_ROOT / "metabolism" / "closure_annotations.tsv"
    if not path.exists():
        st.info("Sin tabla de cierres anotados.")
        return
    df = closures_ann_cached(str(path))
    if df.empty:
        st.info("La tabla de cierres está vacía.")
        return
    genomes = sorted(df["genome"].unique())
    genome = st.selectbox("Genoma", genomes, key="clo_genome", help=T["help_genoma"])
    sub = df[df["genome"] == genome].copy()
    states = sorted(df["closure_status"].dropna().unique())
    sel = st.multiselect("Estado", states, default=states, key="clo_state", help=T["help_estado"])
    sub = sub[sub["closure_status"].isin(sel)]
    c = st.columns(5)
    c[0].metric("Cierres", len(sub), help=T["closures_m_total"])
    _cl_help = {"candidate_thermo_forward": T["closures_m_forward"],
                "candidate_thermo_reversible": T["closures_m_reversible"],
                "candidate_thermo_reverse": T["closures_m_reverse"],
                "candidate_thermo_unknown": T["closures_m_unknown"]}
    for i, key in enumerate(["candidate_thermo_forward", "candidate_thermo_reversible",
                             "candidate_thermo_reverse", "candidate_thermo_unknown"], 1):
        c[i].metric(key.replace("candidate_", ""),
                    int((sub["closure_status"] == key).sum()), help=_cl_help[key])
    st.caption(T["filter_hint"])
    cols = [c for c in ("gap","domain","query","seqid","start","end","closure_status",
                         "direction_status","thermo_dg","review_status","source","support")
            if c in sub.columns]
    st.dataframe(sub[cols], width="stretch", hide_index=True)
    if len(sub):
        gap_sel = st.selectbox("Ver evidencia de…", sub["gap"].tolist(), key="clo_detail")
        r = sub[sub["gap"] == gap_sel].iloc[0]
        st.markdown(f"### {gap_sel} · `{r.get('closure_status','')}`")
        st.write(f"**Dominio:** {r.get('domain','')} · **ORF:** {r.get('query','')}")
        st.write(f"**Locus:** {r.get('seqid','')}:{r.get('start','')}-{r.get('end','')} · "
                 f"**Soporte 400:** {r.get('support','')}")
        st.write(f"**Fuente del cierre:** {r.get('source','')}")
        st.write(f"**Dirección:** {r.get('direction_status','')} · **ΔG′:** {r.get('thermo_dg','')} kJ/mol")
        st.write(f"**Motivo:** {r.get('direction_reason','')} · **Reacción:** {r.get('thermo_reaction','')}")
        st.write(f"**Revisión:** {r.get('review_status','') or 'sin revisar'}")
        gff = DATA_ROOT / "metabolism" / "closures_gff" / f"{genome}.gff3"
        st.code(str(gff))
        st.markdown("Fuentes: [Rhea](https://www.rhea-db.org/) · "
                    "[eQuilibrator 3.0](https://doi.org/10.1093/nar/gkab1106) · "
                    "[TECRDB](https://doi.org/10.1093/bioinformatics/bth314) · Pfam/HMM.")
    st.download_button("⬇️ Cierres (TSV)", sub.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_closures.tsv", mime="text/plain")
    st.caption(T["closures_limit"])
    st.caption(credits.line("hmmer", "rhea", "equilibrator"))


@st.cache_data(show_spinner=False)
def operons_cached(root: str):
    return load_components(Path(root))


@st.cache_data(show_spinner=False)
def operon_detail_cached(root: str):
    return load_detail(Path(root))


def render_operons() -> None:
    st.caption(T["operons_intro"])
    if not operons_available(DATA_ROOT):
        st.info(
            "Sin tabla de componentes. Genera con "
            "`python scripts/validation/operon_components.py`."
        )
        return
    df = operons_cached(str(DATA_ROOT))
    if df.empty:
        st.info("La tabla de componentes está vacía.")
        return

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Operones", len(df), help=T["operons_m_total"])
    c2.metric("Multigénicos", int((df["n_genes_i"] >= 2).sum()), help=T["operons_m_multi"])
    c3.metric("Con terminador", int((df.get("has_terminator") == "1").sum()),
              help=T["operons_m_term"])
    c4.metric("Completos (7/7)", int((df["complete"] == "yes").sum()),
              help=T["operons_m_complete"])
    st.caption(T["operons_pieces"])

    f1, f2, f3 = st.columns(3)
    multi = f1.checkbox("Solo multigénicos", value=True, key="op_multi")
    only_term = f2.checkbox("Solo con terminador", value=False, key="op_term")
    confs = sorted(x for x in df["confidence"].unique() if x)
    sel_conf = f3.multiselect("Confianza", confs, default=confs, key="op_conf")

    view = df.copy()
    if multi:
        view = view[view["n_genes_i"] >= 2]
    if only_term:
        view = view[view.get("has_terminator") == "1"]
    if sel_conf:
        view = view[view["confidence"].isin(sel_conf)]

    st.dataframe(piece_table(view), width="stretch", hide_index=True)
    st.download_button(
        "⬇️ Componentes (TSV)", view.to_csv(sep="\t", index=False),
        file_name="operon_components.tsv", mime="text/plain")

    cov = operons_coverage(view)
    if cov:
        st.caption("Cobertura en la selección actual: " +
                   " · ".join(f"{k} {v}/{len(view)}" for k, v in cov.items()))

    if view.empty:
        return
    st.subheader("Detalle de un operón")
    names = sorted(view["operon_name"])
    pick = st.selectbox("Operón", names, key="op_pick", help=T["help_operon"])
    row = view[view["operon_name"] == pick].iloc[0]
    st.markdown(
        f"**{pick}** · {row['n_genes']} genes · {row['n_tus']} unidades de "
        f"transcripción · evidencia `{row['evidence_codes'] or 'n/d'}` · "
        f"{row['n_pmids']} PMIDs"
    )
    detail = operon_detail_cached(str(DATA_ROOT))
    sub = detail_for(detail, row["operon_id"])
    if sub.empty:
        st.info("Este operón no tiene unidades de transcripción enlazadas.")
    else:
        st.dataframe(sub, width="stretch", hide_index=True)
    st.caption(T["operons_tu"])
    st.caption(T["operons_why"])
    st.caption(credits.line("regulondb"))
    render_rule_operons()
    render_hairpins()
    render_suspects()


@st.cache_data(show_spinner=False)
def rule_ops_cached(root: str, genome: str):
    return load_rule_operons(Path(root), genome)


@st.cache_data(show_spinner=False)
def rule_bnd_cached(root: str, genome: str):
    return load_rule_boundaries(Path(root), genome)


def render_rule_operons() -> None:
    st.subheader("⚖️ Operones por reglas (determinista, sin modelo)")
    st.caption(T["rules_intro"])
    genomes = rule_genomes(DATA_ROOT)
    if not genomes:
        st.info(
            "Sin operones por reglas. Genera con "
            "`python scripts/regulation/rule_operons.py`."
        )
        return
    genome = st.selectbox("Genoma (reglas)", genomes, key="rule_genome")
    if str(genome).endswith("(excluido)"):
        st.caption(
            "⚠️ **Genoma excluido**: posible parásito o simbionte intracelular. "
            "Está anotado para completar el catálogo, pero no entra en las "
            "métricas del conjunto libre porque su biología es distinta."
        )
    ops = rule_ops_cached(str(DATA_ROOT), genome)
    bnd = rule_bnd_cached(str(DATA_ROOT), genome)
    if ops.empty:
        st.info("Sin datos para ese genoma.")
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Operones", len(ops), help=T["rules_m_total"])
    c2.metric("Fronteras decididas", len(bnd), help=T["rules_m_borders"])
    c3.metric("Cortes", int((bnd["cut"] == "1").sum()) if "cut" in bnd else 0,
              help=T["rules_m_cuts"])

    summ = rule_summary(bnd)
    for k, v in sorted(summ.items(), key=lambda kv: -kv[1]):
        st.markdown(f"- {RULE_LABEL.get(k, k)}: **{v}**")
    for _k in ("rules_r1", "rules_r1b", "rules_r2", "rules_r3", "rules_r4", "rules_r5"):
        st.caption(T[_k])
    if "why" in bnd.columns:
        st.dataframe(bnd[["left", "right", "rule", "cut", "why"]].head(200),
                     width="stretch", hide_index=True)
        st.download_button("⬇️ Fronteras con motivo (TSV)",
                           bnd.to_csv(sep="\t", index=False),
                           file_name=f"{genome}_boundary_rules.tsv",
                           mime="text/plain")


@st.cache_data(show_spinner=False)
def hairpins_cached(root: str, genome: str):
    return load_hairpins(Path(root), genome)


def render_hairpins() -> None:
    st.subheader("🌀 Horquillas de terminación")
    st.caption(T["hairpins_intro"])
    genomes = hairpin_genomes(DATA_ROOT)
    if not genomes:
        st.info(
            "Sin horquillas anotadas. Genera con "
            "`bash run_pipeline.sh hairpins` (K-12) o "
            "`python scripts/regulation/annotate_hairpins.py --gff ... --fasta ...`."
        )
        return
    genome = st.selectbox("Genoma con horquillas", genomes, key="hp_genome")
    df = hairpins_cached(str(DATA_ROOT), genome)
    if df.empty:
        st.info("Sin datos para ese genoma.")
        return
    s = hairpin_summary(df)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Regiones con horquilla", len(df), help=T["hairpins_m_total"])
    c2.metric("Clase A (fuerte)", s.get("A_fuerte", 0), help=T["hairpins_m_a"])
    c3.metric("Clase B (media)", s.get("B_media", 0), help=T["hairpins_m_b"])
    c4.metric("Clase C (débil)", s.get("C_debil", 0), help=T["hairpins_m_c"])

    sel = st.multiselect("Clase", ["A_fuerte", "B_media", "C_debil"],
                         default=["A_fuerte"], key="hp_class",
                         help=T["hairpins_class_help"])
    view = df[df["class"].isin(sel)] if sel else df
    st.dataframe(hairpin_table(view), width="stretch", hide_index=True)
    st.download_button("⬇️ Horquillas (TSV)", view.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_hairpins.tsv", mime="text/plain")
    if "dG" in df.columns and view["dG_i"].notna().any():
        st.caption(f"ΔG mediana de la selección: "
                   f"{view['dG_i'].median():.1f} kcal/mol")
        st.caption(T["hairpins_dg"])
    st.caption(credits.line("vienna"))



def render_mobile() -> None:
    st.caption(T["mobile_intro"])
    path = DATA_ROOT / "mobile_defense" / "features.tsv"
    if not path.exists():
        st.info("Sin tabla de móviles/defensa.")
        return
    df = csv_cached(str(path))
    genomes = sorted(df["genome"].unique())
    genome = st.selectbox("Genoma", genomes, key="mob_genome", help=T["help_genoma"])
    sub = df[df["genome"] == genome].copy()
    types = sorted(df["type"].unique())
    sel = st.multiselect("Tipo", types, default=types, key="mob_types", help=T["help_tipo"])
    sub = sub[sub["type"].isin(sel)]
    c = st.columns(7)
    c[0].metric("Total", len(sub), help=T["mobile_types"])
    _mob = {"mobile_element": "Trozos de ADN que se copian y se pegan en otro sitio.",
            "prophage": "Virus bacteriano dormido dentro del genoma.",
            "defense": "Sistemas que protegen a la bacteria de virus y ADN extraño.",
            "crispr": "Sistema inmune que guarda memoria de virus pasados.",
            "amr": "Genes que dan resistencia a antibióticos.",
            "virulence": "Genes que ayudan a causar daño al hospedador."}
    for i, t in enumerate(["mobile_element","prophage","defense","crispr","amr","virulence"], 1):
        c[i].metric(t, int((sub["type"] == t).sum()), help=_mob[t])
    st.dataframe(sub.head(500), width="stretch", hide_index=True)
    st.download_button("⬇️ Descargar TSV", sub.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_mobile_defense.tsv", mime="text/plain")
    st.caption(T["mobile_limit"])



def render_family_explanations() -> None:
    st.caption(T["family_intro"])
    path = DATA_ROOT / "metabolism" / "pathway_explanations.tsv"
    if not path.exists():
        st.info("Sin tabla de explicaciones por familia.")
        return
    df = csv_cached(str(path))
    cols = [c for c in ("family","gaps","transport","t1","homolog_strong","homolog_weak",
                        "no_homolog","homology_present","mapped_homology","mapped_validated","high_support",
                        "no_homolog_conclusion","annotation_gap","low_support","conclusion","next_action")
            if c in df.columns]
    st.dataframe(df[cols], width="stretch", hide_index=True)
    fam = st.selectbox("Ver ejemplos de…", df["family"].tolist(), key="fam_sel")
    r = df[df["family"] == fam].iloc[0]
    st.markdown(f"### {fam}")
    st.write(f"**Brechas sin homólogo:** {r['gaps']} · **transporte:** {r['transport']} · "
             f"**homólogo presente:** {r.get('homology_present','')} · "
             f"**mapeo corregido:** {r.get('mapped_homology','')} · "
             f"**sin homólogo:** {r.get('no_homolog','')}")
    st.write(f"**Explicación:** {r['conclusion']}")
    st.write(f"**Acción:** {r.get('next_action','')}")
    if r.get('no_homolog_conclusion'):
        st.write(f"**Sin homólogo:** {r['no_homolog_conclusion']}")
    st.write(f"**Ejemplos:** {r['examples']}")
    st.download_button("⬇️ Explicaciones (TSV)", df.to_csv(sep="\t", index=False),
                       file_name="pathway_explanations.tsv", mime="text/plain")
    st.caption(T["family_limit"])
    st.caption(credits.line("rhea", "uniprot", "equilibrator", "kegg"))



def render_mapping() -> None:
    st.caption(T["mapping_intro"])
    path = DATA_ROOT / "metabolism" / "family_drilldown.tsv"
    if not path.exists():
        st.info("Sin tabla de drill-down de mapeo.")
        return
    df = csv_cached(str(path))
    fams = sorted(df["family"].unique())
    fam = st.selectbox("Familia", fams, key="map_fam")
    sub = df[df["family"] == fam].copy()
    states = sorted(sub["status"].unique())
    sel = st.multiselect("Estado", states, default=states, key="map_status", help=T["help_estado"])
    sub = sub[sub["status"].isin(sel)]
    c = st.columns(5)
    c[0].metric("Brechas", len(sub), help=T["mapping_m_total"])
    _map_help = {"homolog_strong": T["mapping_m_strong"],
                 "homolog_weak": T["mapping_m_weak"],
                 "homolog_low": T["mapping_m_low"],
                 "no_homolog": T["mapping_m_none"]}
    for i, k in enumerate(("homolog_strong","homolog_weak","homolog_low","no_homolog"), 1):
        c[i].metric(k, int((sub["status"] == k).sum()), help=_map_help[k])
    st.dataframe(sub[["genome","gap","support","transport","best_acc","pident","qcov","status"]],
                 width="stretch", hide_index=True)
    st.download_button("⬇️ Mapeo (TSV)", sub.to_csv(sep="\t", index=False),
                       file_name=f"{fam}_mapping.tsv", mime="text/plain")
    st.caption(credits.line("diamond"))
    corr = DATA_ROOT / "metabolism" / "mapping_corrections.tsv"
    if corr.exists():
        st.markdown("#### Correcciones de mapeo aplicadas")
        cdf = csv_cached(str(corr))
        val_path = DATA_ROOT / "metabolism" / "mapping_validation.tsv"
        if val_path.exists():
            vdf = csv_cached(str(val_path))
            cdf = cdf.merge(vdf[["genome","gap","closed"]], on=["genome","gap"], how="left")
            st.metric("Validados en red corregida", int((cdf["closed"] == 1).sum()),
                      help="Correcciones que, al aplicarse, hacen desaparecer la brecha de la red. Es la comprobación más directa.")
        ec_path = DATA_ROOT / "metabolism" / "mapping_ec.tsv"
        if ec_path.exists():
            edf = csv_cached(str(ec_path))
            cdf = cdf.merge(edf[["genome","gap","ec_status","producer_ec","current_ec"]],
                            on=["genome","gap"], how="left")
            st.metric("EC cambiado", int((cdf["ec_status"] == "ec_mismatch").sum()),
                      help="Casos donde el número EC propuesto para la enzima no coincidía con el esperado.")
        up_path = DATA_ROOT / "metabolism" / "mapping_ec_uniprot.tsv"
        if up_path.exists():
            udf = csv_cached(str(up_path))
            cdf = cdf.merge(udf[["genome","gap","uniprot_ec","status"]].rename(columns={"status":"uniprot_ec_status"}),
                            on=["genome","gap"], how="left")
            st.metric("EC UniProt coincide", int((cdf["uniprot_ec_status"] == "uniprot_match").sum()),
                      help="Casos donde la base UniProt confirma el número EC de la enzima propuesta.")
        rh_path = DATA_ROOT / "metabolism" / "mapping_rhea_uniprot.tsv"
        if rh_path.exists():
            rdf = csv_cached(str(rh_path))
            cdf = cdf.merge(rdf[["genome","gap","status"]].rename(columns={"status":"uniprot_rhea_status"}),
                            on=["genome","gap"], how="left")
            st.metric("Rhea UniProt respalda", int((cdf["uniprot_rhea_status"] == "rhea_supported").sum()),
                      help="Casos donde las proteínas de UniProt asociadas a la reacción de Rhea incluyen la enzima propuesta.")
        audit_path = DATA_ROOT / "metabolism" / "validation" / "mapping_audit.tsv"
        if audit_path.exists():
            adf = csv_cached(str(audit_path))
            cdf = cdf.merge(adf[["genome","gap","verdict","reason"]].rename(columns={"reason":"audit_reason"}),
                            on=["genome","gap"], how="left")
            st.metric("Auditoría support", int((cdf["verdict"] == "support").sum()),
                      help="Revisión automática que respalda la corrección.")
            st.metric("Auditoría reject", int((cdf["verdict"] == "reject").sum()),
                      help="Revisión automática que rechaza la corrección. Estas quedan en revisión.")
        ctx_path = DATA_ROOT / "metabolism" / "mapping_context.tsv"
        if ctx_path.exists():
            cdf2 = csv_cached(str(ctx_path))
            cdf = cdf.merge(cdf2[["genome","gap","status"]].rename(columns={"status":"mapping_context"}),
                            on=["genome","gap"], how="left")
            st.metric("Contexto misma ruta", int((cdf["mapping_context"] == "same_pathway").sum()),
                      help="Casos donde la enzima propuesta participa en la misma ruta que la brecha, lo que apoya la corrección.")
        show = [c for c in ("genome","gap","query","current_acc","producer_acc","producer_rhea",
                            "closed","ec_status","uniprot_ec","uniprot_ec_status",
                            "uniprot_rhea_status","mapping_context","verdict","audit_reason") if c in cdf.columns]
        st.dataframe(cdf[show].head(200), width="stretch", hide_index=True)
        st.download_button("⬇️ Correcciones (TSV)", cdf.to_csv(sep="\t", index=False),
                           file_name="mapping_corrections.tsv", mime="text/plain")
        st.caption(T["mapping_corrections"])



def render_completeness() -> None:
    st.caption(T["completeness_intro"])
    sp = DATA_ROOT / "metabolism" / "pathway_completeness_summary.tsv"
    fp = DATA_ROOT / "metabolism" / "pathway_completeness.tsv"
    if not sp.exists() or not fp.exists():
        st.info("Sin tabla de completitud. Ejecuta scripts/metabolism/40_pathway_completeness.py.")
        return
    summary = csv_cached(str(sp))
    full = csv_cached(str(fp))
    st.dataframe(summary.head(25), width="stretch", hide_index=True)
    pathway = st.selectbox("Ruta", summary["pathway"].tolist(), key="comp_path", help=T["help_ruta"])
    sub = full[full["pathway"] == pathway]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Genomas", len(sub), help=T["completeness_m_genomes"])
    c2.metric("Cobertura media", f"{sub['coverage'].mean():.3f}",
              help=T["completeness_m_cov"])
    c3.metric("Gaps", int(sub["gaps"].sum()), help=T["completeness_m_gaps"])
    c4.metric("Gaps esenciales", int(sub["essential_gaps"].sum()),
              help=T["completeness_m_ess"])
    st.dataframe(sub.sort_values("coverage"), width="stretch", hide_index=True)
    st.download_button("⬇️ Completitud (TSV)", sub.to_csv(sep="\t", index=False),
                       file_name=f"{pathway.split()[0]}_completeness.tsv", mime="text/plain")
    st.caption("Rutas Rhea-only no están incluidas: solo las que existen en KEGG.")
    st.caption(credits.line("kegg", "rhea"))



def render_mobile_typed() -> None:
    st.caption(T["mobile_typed_intro"])
    pilot = DATA_ROOT / "mobile_typing"
    if not pilot.exists():
        st.info("Sin capa tipada de móviles.")
        return
    base = pilot
    gff_base = pilot / "gff"
    sp_base = pilot / "crispr_spacers"
    full = pilot / "full"
    full_mode = False
    if full.exists() and list(full.glob("*.tsv")):
        src = st.radio("Fuente", ["Piloto (10)", "Full (400)"], horizontal=True, key="mobt_src")
        if src.startswith("Full"):
            base = full; gff_base = pilot / "gff" / "full"; sp_base = pilot / "crispr_spacers" / "full"
            full_mode = True
    files = sorted(p for p in base.glob("*.tsv") if p.name != "pilot_genomes.tsv")
    if not files:
        st.info("Sin TSV tipados.")
        return
    genome = st.selectbox("Genoma", [p.stem for p in files], key="mobt_genome", help=T["help_genoma"])
    path = base / f"{genome}.tsv"
    df = csv_cached(str(path))
    types = sorted(df["element_type"].unique())
    sel = st.multiselect("Tipo", types, default=types, key="mobt_types", help=T["help_tipo"])
    view = df[df["element_type"].isin(sel)]
    c = st.columns(5)
    c[0].metric("Filas", len(view))
    c[1].metric("Alta", int((view["confidence"] == "high").sum()))
    c[2].metric("Media", int((view["confidence"] == "medium").sum()))
    c[3].metric("Baja", int((view["confidence"] == "low").sum()))
    c[4].metric("Tipos", view["element_type"].nunique())
    st.dataframe(view, width="stretch", hide_index=True)
    gff = gff_base / f"{genome}.gff3"
    if gff.exists():
        st.download_button("⬇️ GFF móviles", gff.read_text(),
                           file_name=f"{genome}_mobile.gff3", mime="text/plain")
    st.download_button("⬇️ TSV móviles", view.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_mobile.tsv", mime="text/plain")
    sp = sp_base / f"{genome}.tsv"
    if sp.exists():
        sdf = csv_cached(str(sp))
        st.markdown(f"#### CRISPR spacers ({len(sdf)})")
        st.dataframe(sdf.head(50), width="stretch", hide_index=True)
    if full_mode:
        st.caption("Full 400: 44.832 filas; 3.170 `ncbi_only` en TSV no se unen por CDS; "
                   "90 genomas sin CRISPR evaluable; usar high/medium.")
    st.caption(T["mobile_typed_limit"])



def render_membrane_typed() -> None:
    st.caption(T["membrane_intro"])
    pilot = DATA_ROOT / "membrane"
    if not pilot.exists():
        st.info("Sin capa de membrana.")
        return
    base = pilot; gff_base = pilot / "gff"
    full_mode = False
    full = pilot / "full"
    if full.exists() and list(full.glob("GCF_*.tsv")):
        src = st.radio("Fuente", ["Piloto (10)", "Full (400)"], horizontal=True, key="mem_src")
        if src.startswith("Full"):
            base = full; gff_base = pilot / "gff" / "full"; full_mode = True
    files = sorted(p for p in base.glob("*.tsv") if p.name != "pilot_genomes.tsv")
    if not files:
        st.info("Sin TSV de membrana.")
        return
    genome = st.selectbox("Genoma", [p.stem for p in files], key="mem_genome", help=T["help_genoma"])
    df = csv_cached(str(base / f"{genome}.tsv"))
    types = sorted(df["type"].unique())
    sel = st.multiselect("Tipo", types, default=types, key="mem_types", help=T["help_tipo"])
    view = df[df["type"].isin(sel)]
    c = st.columns(5)
    c[0].metric("Filas", len(view))
    c[1].metric("Alta", int((view["confidence"] == "high").sum()))
    c[2].metric("Media", int((view["confidence"] == "medium").sum()))
    c[3].metric("Transportadores", int((view["type"] == "transporter").sum()))
    c[4].metric("Secreción", int((view["type"] == "secretion").sum()))
    st.dataframe(view, width="stretch", hide_index=True)
    gff = gff_base / f"{genome}.gff3"
    if gff.exists():
        st.download_button("⬇️ GFF membrana", gff.read_text(),
                           file_name=f"{genome}_membrane.gff3", mime="text/plain")
    st.download_button("⬇️ TSV membrana", view.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_membrane.tsv", mime="text/plain")
    if full_mode:
        st.caption("Full 400: 319.715 filas; 66,4 % low; 5.456 transportadores "
                   "con sustrato explícito; usar high/medium.")
    st.caption(T["membrane_limit"])
    st.caption(credits.line("uniprot"))



def render_taxonomy() -> None:
    st.caption(T["taxonomy_intro"])
    pilot = DATA_ROOT / "taxonomy"
    if not pilot.exists():
        st.info("Sin capa de taxonomía.")
        return
    base = pilot
    full = pilot / "full"
    full_mode = False
    if full.exists() and list(full.glob("GCF_*.json")):
        src = st.radio("Fuente", ["Piloto (10)", "Full (400)"], horizontal=True, key="tax_src")
        if src.startswith("Full"):
            base = full; full_mode = True
    records = []
    for p in sorted(base.glob("GCF_*.json")):
        d = json.load(open(p))
        records.append({
            "genome": d.get("genome", p.stem),
            "species_ncbi": d.get("species_ncbi", d.get("species", "")),
            "species_gtdb": d.get("species_gtdb", ""),
            "ani": d.get("ani"),
            "closest_ref": d.get("closest_ref"),
            "confidence": d.get("confidence"),
            "status": d.get("status"),
            "conflicts": bool((d.get("evidence") or {}).get("conflicts")),
            "_raw": d,
        })
    if not records:
        st.info("Sin JSON de taxonomía.")
        return
    st.dataframe(pd.DataFrame([{k: v for k, v in r.items() if k != "_raw"} for r in records]),
                 width="stretch", hide_index=True)
    genome = st.selectbox("Genoma", [r["genome"] for r in records], key="tax_genome", help=T["help_genoma"])
    r = next(x for x in records if x["genome"] == genome)
    d = r["_raw"]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Especie NCBI", r["species_ncbi"] or "n/d", help=T["taxonomy_m_ncbi"])
    c2.metric("Especie GTDB", r["species_gtdb"] or "n/d", help=T["taxonomy_m_gtdb"])
    c3.metric("ANI", "sin ANI" if r["ani"] is None else f"{float(r['ani']):.2f}",
              help=T["taxonomy_m_ani"])
    c4.metric("Confianza", r["confidence"] or "n/d", help=T["taxonomy_m_conf"])
    st.write(f"**Estado:** {r['status']} · **Referencia:** {r['closest_ref'] or 'n/d'} · "
             f"**Conflictos:** {'sí' if r['conflicts'] else 'no'}")
    with st.expander("Linaje NCBI", expanded=True):
        st.code(d.get("ncbi_lineage", ""))
    with st.expander("Linaje GTDB", expanded=True):
        st.code(d.get("gtdb_lineage", ""))
    if d.get("evidence"):
        with st.expander("Evidencia ANI/taxonomía"):
            st.json(d["evidence"])
    if full_mode:
        st.caption("Full 400: 117 sin ANI, 51 sin GTDB r232, 18 provisional. "
                   "`ani=null` = sin referencia, no ANI bajo.")
    st.caption(T["taxonomy_limit"])
    st.caption(credits.line("fastani") if "fastani" in credits.SHORT else "")



def render_chains() -> None:
    st.caption(T["chains_intro"])
    base = DATA_ROOT / "chains"
    files = sorted(base.glob("*.tsv")) if base.exists() else []
    files = [f for f in files if f.name not in ("families.tsv", "sources.tsv", "targets.tsv")]
    if not files:
        st.info("Sin cadenas.")
        return
    df = pd.concat([csv_cached(str(f)) for f in files], ignore_index=True)
    fams = sorted(df["family"].unique())
    fam = st.selectbox("Familia", fams, key="chain_fam")
    sub = df[df["family"] == fam].copy()
    c1, c2, c3 = st.columns(3)
    c1.metric("Cadenas", len(sub), help=T["chains_m_total"])
    c2.metric("Alternativas", int((sub["alternative_route"] == "yes").sum()),
              help=T["chains_m_alt"])
    c3.metric("Sin co-sustratos",
              int((sub["n_co_substrates"].astype(int) == 0).sum()),
              help=T["chains_m_noco"])
    only_alt = st.checkbox("Solo ruta alternativa", value=True, key="chain_alt",
                           help="Muestra solo caminos que evitan el paso que falta.")
    if only_alt:
        sub = sub[sub["alternative_route"] == "yes"]
    max_len = st.slider("Longitud máxima", 1, 6, 6, key="chain_len",
                        help=T["chains_len"])
    sub = sub[sub["length"].astype(int) <= max_len]
    cols = [c for c in ("source","target","length","reversible","alternative_route",
                         "n_co_substrates","path","path_metabolites") if c in sub.columns]
    st.dataframe(sub[cols].sort_values(["length","target"]), width="stretch", hide_index=True)
    st.download_button("⬇️ Cadenas (TSV)", sub.to_csv(sep="\t", index=False),
                       file_name=f"{fam}_chains.tsv", mime="text/plain")
    st.caption(T["chains_limit"])
    st.caption(credits.line("rhea"))



def render_transport_t1() -> None:
    st.caption(T["transport_intro"])
    base = DATA_ROOT / "metabolism" / "transport_t1"
    hits = base / "t1_hits.tsv"
    if not hits.exists():
        st.info("Sin capa T1.")
        return
    df = csv_cached(str(hits))
    c1, c2, c3 = st.columns(3)
    c1.metric("T1", len(df), help=T["transport_m_t1"])
    c2.metric("Genomas", df["genome"].nunique() if len(df) else 0,
              help=T["transport_m_genomes"])
    c3.metric("Alta confianza", int((df["confidence"] == "high").sum()) if len(df) else 0,
              help=T["transport_m_high"])
    st.dataframe(df, width="stretch", hide_index=True)
    st.download_button("⬇️ T1 (TSV)", df.to_csv(sep="\t", index=False),
                       file_name="t1_hits.tsv", mime="text/plain")
    summ = base / "summary.json"
    if summ.exists():
        with st.expander("Resumen y límites"):
            st.json(json.load(open(summ)))
    st.caption(T["transport_limit"])



def render_reannotation() -> None:
    st.caption(
        "Re-anotación propuesta por dominio (77 accionables). "
        "`reannotated` es candidato, no confirmado."
    )
    path = DATA_ROOT / "metabolism" / "reannotation" / "reannotation.tsv"
    if not path.exists():
        st.info("Sin re-anotaciones.")
        return
    df = csv_cached(str(path))
    c1, c2, c3 = st.columns(3)
    c1.metric("Re-anotar", int((df["decision"] == "reannotate").sum()))
    c2.metric("Validar dominio", int((df["decision"] == "validate_domain").sum()))
    c3.metric("Mantener brecha", int((df["decision"] == "keep_missing").sum()))
    dec = st.multiselect("Decisión", sorted(df["decision"].unique()),
                         default=sorted(df["decision"].unique()), key="reann_dec")
    view = df[df["decision"].isin(dec)]
    st.dataframe(view, width="stretch", hide_index=True)
    st.download_button("⬇️ Re-anotación (TSV)", view.to_csv(sep="\t", index=False),
                       file_name="reannotation.tsv", mime="text/plain")
    st.caption("HMM y ORF candidato; sin EC/estructura/experimento. "
               "46 keep_missing = sin evidencia para re-anotar, no ausencia confirmada.")


# genomas que la interfaz puede mostrar: los 400 del proyecto + K-12 curado
DEMO_GENOME = (Path(__file__).resolve().parent / "demo" / "genes.tsv").exists()
GENOME_CHOICES = (
    (["demo"] if DEMO_GENOME else [])
    + (["ecoli_k12"] if (DATA_ROOT / "regulation" / "ecoli_k12" / "genes.tsv").exists() else [])
    + list(PROJECT_GENOMES)
)


def _global_genome() -> str:
    """Genoma elegido en la barra lateral; lo comparten todas las pestañas."""
    g = st.session_state.get("genoma_global")
    if g in GENOME_CHOICES:
        return g
    return GENOME_CHOICES[0] if GENOME_CHOICES else ""


LAYOUTS = {
    "A": [
        ("🧪 Metabolismo y rutas", [
            ("📚 Catálogo", render_catalog),
            ("🗺️ KEGG", render_kegg),
            ("🧬 Brechas", render_gaps),
            ("🔬 Cierres", render_closures),
            ("🧭 Rutas curadas", render_curated),
            ("🚚 Transporte", render_transport_t1),
            ("🔗 Cadenas", render_chains),
            ("📊 Completitud", render_completeness),
        ]),
        ("🧩 Genoma y elementos", [
            ("🧩 Operones", render_operons),
            ("🦠 Móviles y defensa", render_mobile),
            ("🧬 Móviles por tipo", render_mobile_typed),
            ("🧫 Membrana", render_membrane_typed),
            ("🛠️ Reanotación", render_reannotation),
            ("🗺️ Mapeo de genes", render_mapping),
        ]),
        ("🧭 Linaje y origen", [
            ("🏷️ Taxonomía", render_taxonomy),
            ("📖 Explicaciones", render_family_explanations),
            ("🧭 Origen de replicación", render_origin),
        ]),
    ],
    "B": [
        ("🔎 Explorar los 400", [
            ("📚 Catálogo", render_catalog),
            ("🗺️ KEGG", render_kegg),
            ("📊 Completitud", render_completeness),
            ("🏷️ Taxonomía", render_taxonomy),
        ]),
        ("🧠 Entender un resultado", [
            ("🧬 Brechas", render_gaps),
            ("🔬 Cierres", render_closures),
            ("🧭 Rutas curadas", render_curated),
            ("📖 Explicaciones", render_family_explanations),
        ]),
        ("🧩 Genoma y regulación", [
            ("🧩 Operones", render_operons),
            ("🦠 Móviles y defensa", render_mobile),
            ("🧫 Membrana", render_membrane_typed),
            ("🗺️ Mapeo de genes", render_mapping),
            ("🛠️ Reanotación", render_reannotation),
        ]),
        ("🚚 Rutas y transporte", [
            ("🚚 Transporte", render_transport_t1),
            ("🔗 Cadenas", render_chains),
            ("🧭 Origen", render_origin),
        ]),
        ("📚 De dónde salen los datos", [
            ("🙏 Créditos y herramientas", credits.render),
        ]),
    ],
    "C": [
        ("🧬 Genoma", [
            ("📊 Completitud", render_completeness),
            ("🏷️ Taxonomía", render_taxonomy),
            ("🧭 Origen de replicación", render_origin),
        ]),
        ("🔤 Genes", [
            ("📚 Catálogo", render_catalog),
            ("🗺️ KEGG", render_kegg),
            ("🗺️ Mapeo de genes", render_mapping),
            ("🛠️ Reanotación", render_reannotation),
        ]),
        ("🧩 Operones", [
            ("🧩 Operones", render_operons),
            ("🧬 Móviles por tipo", render_mobile_typed),
        ]),
        ("🧪 Rutas", [
            ("🧭 Rutas curadas", render_curated),
            ("🚚 Transporte", render_transport_t1),
            ("🔗 Cadenas", render_chains),
        ]),
        ("🕳️ Brechas", [
            ("🧬 Brechas", render_gaps),
            ("🔬 Cierres", render_closures),
        ]),
        ("🎛️ Regulación", [
            ("📖 Explicaciones", render_family_explanations),
            ("🦠 Móviles y defensa", render_mobile),
            ("🧫 Membrana", render_membrane_typed),
        ]),
    ],
}

LAYOUT_INTRO = {
    "A": "Secciones temáticas: el mismo orden de siempre, con el genoma "
         "compartido entre pestañas.",
    "B": "Empieza por lo que quieres hacer, no por el nombre del análisis.",
    "C": "Recorre el genoma como lo haría un biólogo: genoma, genes, operones, "
         "rutas y brechas.",
}


@st.cache_data(show_spinner=False)
def suspects_cached(root: str, genome: str):
    return load_suspects(Path(root), genome)


def render_suspects() -> None:
    st.subheader("🔗 Sospecha por cadena metabólica")
    st.caption(
        "Genes **contiguos** cuyas enzimas participan en la **misma red de "
        "reacciones** (comparten metabolitos). Es una pista barata de operón "
        "metabólico: muchas rutas contiguas no son operones y muchos operones "
        "no tienen ruta asociada. Se marca como **sospecha**, nunca como hecho."
    )
    genome = _global_genome()
    df = suspects_cached(str(DATA_ROOT), genome) if genome else None
    if df is None or df.empty:
        st.info("Sin sospechas para este genoma. Genera con "
                "`bash run_pipeline.sh suspects`.")
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Grupos sospechosos", len(df))
    c2.metric("Genes implicados", int(df["n_genes"].sum()))
    c3.metric("Alta confianza", int((df.get("confidence") == "high").sum()))
    cols = [c for c in ("contig", "strand", "start", "end", "n_genes", "genes",
                        "gap_mean", "shared_metabolites", "confidence")
            if c in df.columns]
    st.dataframe(df[cols].head(200), width="stretch", hide_index=True)
    st.download_button("⬇️ Sospechas (TSV)", df.to_csv(sep="\t", index=False),
                       file_name=f"{genome}_operon_suspects.tsv", mime="text/plain")



def render_essential_pathways() -> None:
    st.caption(
        "Rutas curadas esenciales para vida libre en medio mínimo. "
        "Core = alta prioridad; conditional = depende de filo/estilo de vida. "
        "No implica que todos los gaps deban cerrarse."
    )
    base = DATA_ROOT / "metabolism" / "essential"
    comp = base / "essential_completeness.tsv"
    if not comp.exists():
        st.info("Sin panel de rutas esenciales.")
        return
    cal = base / "essential_calibration.tsv"
    cdf = csv_cached(str(cal)) if cal.exists() else csv_cached(str(comp))
    qdf = csv_cached(str(base / "essential_gap_queue.tsv")) if (base/"essential_gap_queue.tsv").exists() else None
    c1, c2, c3 = st.columns(3)
    c1.metric("Rutas", len(cdf))
    c2.metric("Core", int((cdf["tier"] == "core").sum()))
    c3.metric("Gaps en esenciales", 0 if qdf is None else len(qdf))
    st.dataframe(cdf, width="stretch", hide_index=True)
    rep = base / "COVERAGE_REPORT.md"
    if rep.exists():
        st.download_button("⬇️ Informe de cobertura esencial", rep.read_text(),
                           file_name="COVERAGE_REPORT.md", mime="text/markdown")
    pri = base / "essential_priority_targets.tsv"
    if pri.exists():
        st.download_button("⬇️ Objetivos prioritarios (TSV)", pri.read_text(),
                           file_name="essential_priority_targets.tsv", mime="text/plain")
    st.caption("Panel curado con literatura: 27 core / 8 conditional. "
               "Keio/LB no prueba esencialidad en medio mínimo; "
               "cobertura baja puede ser salvage/transporte.")
    asp = base / "ACTIVE_SITES_ESSENTIAL.tsv"
    if asp.exists():
        adf = csv_cached(str(asp))
        with st.expander(f"Sitios activos de referencias ({len(adf)} filas)"):
            at = st.selectbox("Target con sitios", sorted(adf["target"].unique()), key="as_target")
            av = adf[adf["target"] == at]
            st.dataframe(av, width="stretch", hide_index=True)
            st.download_button("⬇️ Sitios activos (TSV)", adf.to_csv(sep="\t", index=False),
                               file_name="ACTIVE_SITES_ESSENTIAL.tsv", mime="text/plain")
            st.caption("Residuos del homólogo, no del ORF; transferir solo por alineación.")
    clp = base / "ACTIVE_SITE_CLOSURE_READINESS.tsv"
    if clp.exists():
        with st.expander("Readiness de cierre por sitio activo"):
            cld = csv_cached(str(clp))
            st.dataframe(cld[cld["readiness"] == "structural_ready"], width="stretch", hide_index=True)
            st.download_button("⬇️ Readiness completo (TSV)", cld.to_csv(sep="\t", index=False),
                               file_name="ACTIVE_SITE_CLOSURE_READINESS.tsv", mime="text/plain")
            st.caption("No cierra gaps: indica qué falta para validar estructuralmente.")
    tp = base / "ACTIVE_SITE_TRANSFER.tsv"
    if tp.exists():
        with st.expander("Transferencia de sitios activos por secuencia"):
            tdf = csv_cached(str(tp))
            st.dataframe(tdf.head(200), width="stretch", hide_index=True)
            st.download_button("⬇️ Transferencia de sitios (TSV)", tdf.to_csv(sep="\t", index=False),
                               file_name="ACTIVE_SITE_TRANSFER.tsv", mime="text/plain")
            st.caption("Nivel secuencia; no estructura ni experimento. `conserved` = mismo residuo alineado.")
    exc = base / "essential_exceptions.tsv"
    if exc.exists():
        edf = csv_cached(str(exc))
        with st.expander(f"Excepciones documentadas ({len(edf)})"):
            st.dataframe(edf, width="stretch", hide_index=True)
    if qdf is not None and len(qdf):
        mid = st.selectbox("Ruta", cdf["map_id"].tolist(), key="ess_map")
        row = cdf[cdf["map_id"] == mid].iloc[0]
        st.markdown(f"### {row['pathway']} · `{mid}` · {row['tier']}")
        full = DATA_ROOT / "metabolism" / "pathway_completeness.tsv"
        if full.exists():
            pc = csv_cached(str(full))
            pv = pc[pc["pathway"].str.startswith(mid + " ")]
            st.dataframe(pv.sort_values("coverage"), width="stretch", hide_index=True)
            st.download_button("⬇️ Completitud por genoma", pv.to_csv(sep="\t", index=False),
                               file_name=f"{mid}_completeness.tsv", mime="text/plain")
        qv = qdf[qdf["map_id"] == mid].sort_values("support", ascending=False)
        st.markdown(f"#### Gaps en esta ruta ({len(qv)})")
        st.dataframe(qv.head(100), width="stretch", hide_index=True)
    st.caption("Límites: curaduría interna; no todas las rutas aplican a todos los filos; "
               "cobertura baja puede ser base Rhea incompleta, no ausencia biológica. "
               "El contexto de estilo de vida está en Brechas.")



def render_targets() -> None:
    st.caption(
        "Candidatos para los 17 targets esenciales P1/P2. "
        "Busca en toda la bolsa de proteínas, no solo en las anotadas. "
        "Score = dominio + homología + ORF no anotado; no es probabilidad."
    )
    base = DATA_ROOT / "metabolism" / "targets"
    p3 = DATA_ROOT / "metabolism" / "targets_p3"
    if p3.exists() and (p3/"target_candidates.tsv").exists():
        src = st.radio("Fuente", ["P1/P2", "P3 excepción"], horizontal=True, key="target_src")
        if src.startswith("P3"):
            base = p3
    path = base / "target_candidates_comparative.tsv"
    if not path.exists():
        path = base / "target_candidates.tsv"
    if not path.exists():
        st.info("Sin candidatos de targets.")
        return
    df = csv_cached(str(path))
    targets = sorted(df["target"].unique())
    t = st.selectbox("Target", targets, key="target_sel")
    sub = df[df["target"] == t].sort_values("score", ascending=False)
    st.dataframe(sub, width="stretch", hide_index=True)
    if "known_exception" in sub.columns and sub["known_exception"].astype(str).str.len().max() > 0:
        st.warning("P3: target con excepción documentada. Revisar la excepción antes de buscar; no es ausencia.")
    for doc, label in (("TOP3_REVIEW.md", "⬇️ Revisión top 3"), ("VALIDATION_6.md", "⬇️ Validación in silico"), ("STRUCTURAL_CANDIDATES.md", "⬇️ Candidatos estructura"), ("STRUCTURE_PIPELINE.md", "⬇️ Pipeline estructura"), ("ACTIVE_SITES.md", "⬇️ Sitios activos homólogos")):
        dp = base / doc
        if dp.exists():
            st.download_button(label, dp.read_text(), file_name=doc, mime="text/markdown")
    st.download_button("⬇️ Candidatos del target (TSV)", sub.to_csv(sep="\t", index=False),
                       file_name=f"{t.replace(' ','_')}_candidates.tsv", mime="text/plain")
    sf = base / "structural_candidates.faa"
    if sf.exists():
        st.download_button("⬇️ 10 ORFs estructurales (FASTA)", sf.read_text(),
                           file_name="structural_candidates.faa", mime="text/plain")
    fasta = base / "candidate_proteins.faa"
    if fasta.exists():
        st.download_button("⬇️ Proteínas candidatas (FASTA)", fasta.read_text(),
                           file_name="candidate_proteins.faa", mime="text/plain")
    st.caption("Límites: contexto = vecino anotado/same strand, no operón completo; "
               "co-ocurrencia 400 añadida; refs sin filtro eucariota; score no es probabilidad. "
               "AFDB solo cubre homólogos conocidos; el ORF exacto requiere ESMFold/ColabFold.")


def render_pathways() -> None:
    st.subheader("🧪 Rutas metabólicas")
    st.caption(
        "Análisis ya calculado de los 400 genomas del proyecto. Cada pestaña "
        "responde una pregunta distinta; empieza por «Brechas» si buscas qué le "
        "falta a una bacteria, y por «Operones» si te interesa cómo se organizan "
        "sus genes."
    )
    if not (DATA_ROOT / "reactions" / "tsv" / "rhea2ec.tsv").exists():
        st.warning("Faltan los recursos de Rhea. Se descargan una sola vez.")
        st.code("python -m webapp.setup_data --data-root data --only rhea")
        if st.button("⬇️ Descargar Rhea ahora"):
            with st.spinner("Descargando Rhea…"):
                subprocess.run([sys.executable, "-m", "webapp.setup_data",
                                "--data-root", str(DATA_ROOT), "--only", "rhea"], check=True)
            st.rerun()
        return
    if not PROJECT_GENOMES:
        st.info("Los modos de los 400 genomas necesitan datos locales. "
                "Puedes generarlos con los scripts de data/README.md; mientras "
                "tanto, la demo del navegador sí está disponible.")
        st.markdown("[🧬 Abrir la demo del navegador]"
                    "(?modo=navegador&genoma=demo)")
        return

    # ---- selector de vista: tres navegaciones sobre las mismas pantallas -----
    # enlace compartible: ?layout=B&genoma=GCF_...
    try:
        qp = st.query_params
        if qp.get("layout") in ("A", "B", "C"):
            st.session_state["layout"] = qp["layout"]
        if qp.get("genoma") and qp["genoma"] in PROJECT_GENOMES:
            st.session_state["genoma_global"] = qp["genoma"]
    except Exception:
        pass

    cfg = _load_ui()
    if "layout" not in st.session_state and cfg.get("layout") in ("A", "B", "C"):
        st.session_state["layout"] = cfg["layout"]
    if (PROJECT_GENOMES and "genoma_global" not in st.session_state
            and cfg.get("genoma") in PROJECT_GENOMES):
        st.session_state["genoma_global"] = cfg["genoma"]

    with st.sidebar.popover("⚙️ Configuración",
                            help="Vista de navegación y genoma activo. "
                                 "Se recuerdan en este equipo."):
        layout = st.radio(
            "Vista de navegación", ["A", "B", "C"], horizontal=True,
            key="layout",
            help="Tres formas de organizar los mismos datos: por secciones (A), "
                 "por tareas (B) o por objeto biológico (C).",
        )
        st.caption({
            "A": "**A · Mínima** — secciones temáticas, lo más parecido a hoy.",
            "B": "**B · Por tareas** — ¿qué quieres hacer?",
            "C": "**C · Por objeto** — Genoma → Genes → Operones → Rutas → Brechas.",
        }[layout])
        if PROJECT_GENOMES:
            cur = _global_genome()
            idx = PROJECT_GENOMES.index(cur) if cur in PROJECT_GENOMES else 0
            st.selectbox("Genoma activo", PROJECT_GENOMES, index=idx,
                         key="genoma_global", help=T["help_genoma"])

    if layout != cfg.get("layout") or _global_genome() != cfg.get("genoma"):
        _save_ui({"layout": layout, "genoma": _global_genome()})

    # Las etiquetas viven en español en LAYOUTS; se traducen al mostrarlas.
    groups = tr_tree(LAYOUTS[layout])
    names = [g[0] for g in groups]
    section = st.sidebar.radio(tr("Sección"), names, key=f"nav_{layout}")

    st.caption(tr(LAYOUT_INTRO[layout]))
    items = dict(groups)[section]
    tabs = st.tabs([label for label, _ in items])
    for tab, (_, fn) in zip(tabs, items):
        with tab:
            fn()


TOOL_DESC = {
    "pyrodigal": "encuentra los genes que producen proteína",
    "trna": "encuentra los ARN de transferencia",
    "rrna": "encuentra los ARN ribosómicos",
    "diamond": "busca parecidos en una base de proteínas ya conocidas",
    "esm": "mide si una secuencia parece una proteína real",
    "ghost": "busca genes en los huecos que nadie anotó",
}

TOOL_LABEL = {
    "pyrodigal": "pyrodigal (CDS)",
    "trna": "tRNAscan-SE (tRNA)",
    "rrna": "barrnap (rRNA)",
    "diamond": "DIAMOND (funcional)",
    "esm": "ESM-2 8M (sentido biológico)",
    "ghost": "Motor de ORFs no anotados",
}


@st.cache_data(show_spinner=False)
def tools_cached() -> dict:
    return detect_tools()


@st.cache_data(show_spinner=False)
def bell_wav() -> bytes:
    """Campana de dos tonos (880 y 1320 Hz) con caída exponencial."""
    rate = 22050
    frames = bytearray()
    for freq in (880.0, 1320.0):
        n = int(rate * 0.35)
        for i in range(n):
            t = i / rate
            env = min(1.0, t * 30.0) * math.exp(-4.0 * t)
            val = 0.4 * env * math.sin(2 * math.pi * freq * t)
            frames += struct.pack("<h", int(val * 32767))
    buf = io.BytesIO()
    with wave.open(buf, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(rate)
        w.writeframes(bytes(frames))
    return buf.getvalue()


def render_result(res: dict) -> None:
    st.subheader(f"🧬 {res['name']}")
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Contigs", res["contigs"], help=T["res_contigs"])
    c2.metric("Tamaño", f"{res['bp']:,} pb", help=T["res_size"])
    c3.metric("GC", f"{res['gc']*100:.1f}%", help=T["res_gc"])
    c4.metric("CDS", res["stats"].get("cds", 0), help=T["res_cds"])
    c5.metric("tRNA / rRNA",
              f"{res['stats'].get('trna', 0)} / {res['stats'].get('rrna', 0)}",
              help=T["res_rna"])
    c6.metric("Desiertos", len(res.get("deserts", [])), help=T["res_deserts"])

    coverage = res.get("coverage", 0.0)
    n_regions = res["stats"].get("regions", 0)
    n_deserts = len(res.get("deserts", []))
    st.progress(min(1.0, coverage),
                text=f"{coverage*100:.1f}% del genoma anotado")
    st.caption(T["res_coverage"])
    st.success(
        f"{coverage*100:.1f}% del genoma anotado. "
        f"Quedan {n_regions} posibles brechas y {n_deserts} desiertos ≥3 kb."
    )
    with st.expander("📄 Informe detallado (descargable)"):
        st.markdown(res.get("report_md", ""))
        st.download_button(
            "⬇️ Descargar informe (Markdown)", res.get("report_md", ""),
            file_name=f"{res['genome']}_report.md", mime="text/markdown",
        )
    st.caption(T["res_next"])
    st.caption(credits.line("pyrodigal", "diamond", "hmmer", "trnascan",
                            "barrnap", "esm"))

    tabs = st.tabs(["🔎 Genes fantasma", "🏜️ Desiertos", "🧪 Brechas",
                    "📥 Descargas", "🧾 Features", "🔬 Hipótesis estructural",
                    "🧬 Función y origen", "🧫 Membrana",
                    "📐 Regulación (reglas)", "🧊 Estructura 3D",
                    "🔀 Comparar genomas"])

    with tabs[0]:
        if res["candidates"]:
            df = pd.DataFrame(res["candidates"])
            cols = [c for c in ["contig", "start", "end", "strand", "aa_len",
                                "start_codon", "codon_llr", "rbs", "esm_llr",
                                "rank", "n_genomes", "label", "label_confidence"]
                    if c in df.columns]
            st.caption(
                f"{len(df):,} ORFs no anotados · "
                f"{sum(1 for r in res['candidates'] if r['rank'] >= 0.75)} de tier alto · "
                "cada fila lleva etiqueta, fuente y confianza"
            )
            st.caption(T["res_ghost"])
            st.dataframe(
                df[cols].sort_values("rank", ascending=False).head(100),
                width="stretch", hide_index=True,
            )
        else:
            st.info("No se generaron ORFs no anotados.")

    with tabs[1]:
        if res.get("deserts"):
            ddf = pd.DataFrame(res["deserts"])
            dcols = [c for c in ["contig", "start", "end", "length", "gc",
                                 "ghost_orfs", "conserved_orfs", "label",
                                 "label_confidence"] if c in ddf.columns]
            st.caption(T["res_deserts_tab"])
            st.dataframe(ddf[dcols].sort_values("length", ascending=False),
                         width="stretch", hide_index=True)
        else:
            st.info("No se encontraron desiertos ≥3 kb.")

    with tabs[2]:
        met = res.get("metabolic")
        if res.get("metabolic_error"):
            st.warning(res["metabolic_error"])
        elif not met:
            st.info("Activa «Construir red metabólica y brechas» para calcularlas.")
        else:
            m1, m2, m3 = st.columns(3)
            m1.metric("Reacciones", met["reactions"])
            m2.metric("Brechas", met["n_gaps"])
            m3.metric("Con soporte", sum(1 for g in met["gaps"] if g["support"] > 0))
            st.caption(f"Red: {met['source']} · {met['proteins']} proteínas.")
            if met["gaps"]:
                dfg = pd.DataFrame(met["gaps"])
                show = [c for c in ("gap","chebi","support","role","origin","status",
                                    "status_reason","direction_status","transport_status","homolog_status","review_class","closure_status","pathways") if c in dfg.columns]
                st.dataframe(dfg[show].head(100), width="stretch", hide_index=True)
                gap_sel = st.selectbox("Ver fuentes de evidencia de…", dfg["gap"].tolist(), key="up_gap_detail")
                r = dfg[dfg["gap"] == gap_sel].iloc[0]
                st.markdown(f"### {gap_sel} · `{r.get('status','')}`")
                with st.expander("1 · Trazado y soporte (Rhea + 400 genomas)", expanded=True):
                    st.write(f"**Soporte:** {r.get('support','')}/400")
                    st.write(f"**Motivo:** {r.get('status_reason','')}")
                    st.write(f"**Fuente:** {r.get('evidence_source','Rhea + 400 genomas')}")
                    st.write(f"**Familia KEGG:** {r.get('pathways','sin ruta')}")
                    st.write(f"**Homólogo de productor:** {r.get('homolog_status','no comprobado')} "
                             f"{r.get('homolog_acc','') or ''} {r.get('homolog_pident','') or ''}")
                with st.expander("2 · Transporte"):
                    st.write(f"**Estado:** {r.get('transport_status','no comprobado')}")
                    st.write(f"**Reacción Rhea:** {r.get('transport_rhea','') or '—'}")
                    st.write(f"**Transportador:** {r.get('transporter_acc','') or '—'} · "
                             f"pident {r.get('transport_pident','') or '—'} · "
                             f"qcov {r.get('transport_qcov','') or '—'}")
                with st.expander("3 · Dirección y termodinámica"):
                    st.write(f"**Estado:** {r.get('direction_status','direction_unspecified')}")
                    st.write(f"**Motivo:** {r.get('direction_reason','Rhea master UN')}")
                    st.write(f"**Termodinámica:** {r.get('thermo_status','pending')}")
                    st.markdown("Fuentes: [Rhea direcciones](https://www.rhea-db.org/help/reaction%2Dside%2Ddirection) · "
                                "[eQuilibrator](https://doi.org/10.1093/nar/gkab1106) · [TECRDB](https://doi.org/10.1093/bioinformatics/bth314)")
                with st.expander("4 · HMM de confirmación"):
                    if met.get("hmm_hits"):
                        st.dataframe(pd.DataFrame(met["hmm_hits"]).head(20), width="stretch", hide_index=True)
                    else:
                        st.write("Sin hits HMM para este genoma/brecha.")
                with st.expander("5 · Cierre por estructura (ESMFold)"):
                    render_esmfold_closure(r, res.get("proteins", []),
                                           key_prefix=f"up_{res['genome']}_{gap_sel}")
                    if res.get("esmfold_closure"):
                        st.markdown("##### Resultados de cierre por estructura")
                        st.dataframe(pd.DataFrame(res["esmfold_closure"]),
                                     width="stretch", hide_index=True)
                from family_logic import explain_families
                fr = explain_families(met["gaps"])
                if fr:
                    st.markdown("#### Familias de rutas")
                    st.dataframe(pd.DataFrame(fr[:12]), width="stretch", hide_index=True)
                sug = [g for g in met["gaps"] if g.get("homolog_rhea")]
                if sug:
                    st.markdown("#### Sugerencias de mapeo")
                    st.dataframe(pd.DataFrame(sug)[["gap","homolog_acc","homolog_rhea","homolog_ec","mapping_suggestion"]].head(20),
                                 width="stretch", hide_index=True)
                    st.caption("Sugerencias revisables: no modifican la red automáticamente.")
                st.download_button(
                    "⬇️ Brechas (TSV)", dfg.to_csv(sep="\t", index=False),
                    file_name=f"{res['genome']}_gaps.tsv", mime="text/plain",
                )
            st.caption(T["gaps_support_help"])
            kc = kegg_catalog_from_net(set(met.get("net", [])),
                                       set(met.get("gap_producers", [])),
                                       DATA_ROOT, OUT_ROOT)
            if kc:
                st.markdown("#### Rutas KEGG con señal en esta red")
                st.dataframe(pd.DataFrame(kc[:20]), width="stretch", hide_index=True)
            if met.get("hmm_hits"):
                st.markdown("#### HMM de confirmación")
                st.dataframe(pd.DataFrame(met["hmm_hits"]).head(50),
                             width="stretch", hide_index=True)
    with tabs[3]:
        st.caption(T["res_downloads"])
        d1, d2, d3, d4 = st.columns(4)
        d1.download_button(
            "⬇️ GFF3", res["gff"], file_name=f"{res['genome']}.gff3",
            mime="text/plain", width="stretch",
        )
        d2.download_button(
            "⬇️ Proteínas (FASTA)", res["proteins_fasta"],
            file_name=f"{res['genome']}_proteins.faa", mime="text/plain",
            width="stretch",
        )
        d3.download_button(
            "⬇️ Candidatos (TSV)", res["candidates_tsv"],
            file_name=f"{res['genome']}_ghost_candidates.tsv", mime="text/plain",
            width="stretch",
        )
        d4.download_button(
            "⬇️ Desiertos (TSV)", deserts_tsv(res.get("deserts", [])),
            file_name=f"{res['genome']}_deserts.tsv", mime="text/plain",
            width="stretch",
        )

    with tabs[4]:
        st.caption(T["res_features"])
        if res["features"]:
            st.dataframe(pd.DataFrame(res["features"]).head(500), width="stretch")
        else:
            st.info("No se anotaron features.")

    with tabs[5]:
        render_hypotheses(res)

    with tabs[6]:
        st.caption(
            "Caracterización automática de proteínas en tres caminos: "
            "A función nombrada, B uncharacterized/DUF, C sin hit Swiss-Prot. "
            "El origen se predice con ORCA + GC skew; DoriC solo aplica a "
            "genomas con accession conocido."
        )
        st.caption(
            "Etiquetas: A = función nombrada; B = homólogo no caracterizado/DUF; "
            "C = sin hit Swiss-Prot. `evidence_strength` = strong/moderate/weak/"
            "none. Una homología fuerte no es una función confirmada."
        )
        fun = res.get("functional")
        if res.get("functional_error"):
            st.warning(res["functional_error"])
        if fun:
            f1, f2, f3, f4 = st.columns(4)
            f1.metric("Proteínas", fun.get("n_proteins", 0))
            f2.metric("Con hit", fun.get("n_hits", 0))
            f3.metric("Cobertura", f"{100*fun.get('match_rate', 0):.1f}%")
            f4.metric("A / B / C", " / ".join(
                str(fun.get("paths", {}).get(k, 0)) for k in ("A", "B", "C")))
            if fun.get("warnings"):
                for w in fun["warnings"]:
                    st.warning(w)
            rows = fun.get("rows") or []
            if rows:
                with st.expander("Límites CPU/GPU de A/B/C y estructura"):
                    st.markdown(
                        "- **A/B/C:** clasificación por homología Swiss-Prot. "
                        "A no es función confirmada; C no es ausencia de función.\n"
                        "- **B estructural:** CPU busca modelos AFDB locales; "
                        "no pliega ni filtra automáticamente.\n"
                        "- **C dominios:** hmmscan contra un panel Pfam curado; "
                        "es una hipótesis de familia, no una función.\n"
                        "- **GPU:** AlphaFold-Multimer/Boltz/ESMFold local para "
                        "estructura; no cambia la etiqueta A/B/C.\n"
                        "- **Regla:** ninguna estructura cierra una brecha."
                    )
                st.caption(
                    "`subunit_hint=named_subunit` solo significa que el producto "
                    "del hit contiene 'subunit' o 'component'. No prueba que sea "
                    "una subunidad ni que pertenezca a un complejo. "
                    "Ver docs/SUBUNIT_NAMES.md."
                )
                hint_n = sum(1 for r in rows
                             if r.get("subunit_hint") == "named_subunit")
                b_afdb = sum(1 for r in rows
                             if r.get("structural_status") == "afdb_model_available")
                c_dom = sum(1 for r in rows if r.get("domain"))
                c_ctx = sum(1 for r in rows if r.get("context_hypothesis"))
                c_tmhmm = sum(1 for r in rows
                              if r.get("tmhmm_status") == "ok")
                c_remote = sum(1 for r in rows
                               if r.get("remote_status") == "hit")
                b_filtered = sum(1 for r in rows
                                 if str(r.get("filter_status", "")).startswith("evaluated"))
                st.caption(f"Con pista de nombre: **{hint_n}** proteínas. "
                           f"Es una pista débil, no una etiqueta.")
                st.caption(f"B con modelo AFDB local: **{b_afdb}**. "
                           f"B con filtros evaluados: **{b_filtered}**. "
                           f"C con dominio Pfam: **{c_dom}**. "
                           f"C con contexto: **{c_ctx}**; "
                           f"C con TMHMM: **{c_tmhmm}**. "
                           f"B/C con homología remota: **{c_remote}**. "
                           "Modelo/dominio/contexto no confirman función ni complejo.")
                dfp = pd.DataFrame(rows)
                f1, f2 = st.columns(2)
                paths = sorted(x for x in dfp.get("path", pd.Series()).dropna().unique())
                pick_paths = f1.multiselect("Filtrar camino", paths, default=paths,
                                            key="pp_paths")
                cstat = sorted(x for x in dfp.get("complex_status", pd.Series()).dropna().unique())
                pick_cstat = f2.multiselect("Filtrar complejo", cstat, default=cstat,
                                            key="pp_cstat")
                view = dfp
                if pick_paths:
                    view = view[view["path"].isin(pick_paths)]
                if pick_cstat and "complex_status" in view:
                    view = view[view["complex_status"].isin(pick_cstat)]
                st.dataframe(view.head(500), width="stretch", hide_index=True)
                st.download_button(
                    "⬇️ Caminos de proteínas (TSV)",
                    view.to_csv(sep="\t", index=False),
                    file_name=f"{res['genome']}_protein_paths.tsv",
                    mime="text/plain",
                )
            comp = res.get("complexes")
            if comp:
                st.markdown("#### 🧩 Complejos por vecindad, catálogo y ortología")
                st.caption(
                    "Solo genes del mismo contig y hebra dentro de la ventana. "
                    "La vecindad no demuestra interacción física. "
                    "Color: verde = completo; azul = transferido inferido; "
                    "amarillo = sospecha/incompleto; naranja = ambiguo; "
                    "gris = no evaluado."
                )
                crs = comp.get("complexes") or []
                suspected = sum(1 for c in crs
                                if c.get("status") == "complex_suspected")
                candidates = sum(1 for c in crs
                                 if c.get("status") == "complex_candidate")
                c1, c2, c3, c4, c5 = st.columns(5)
                c1.metric("Completos", comp.get("n_complete", 0))
                c2.metric("Incompletos", comp.get("n_incomplete", 0))
                c3.metric("Ambiguos", comp.get("n_ambiguous", 0))
                c4.metric("Sospecha fuerte", suspected)
                c5.metric("Candidatos", candidates)
                with st.expander("Cómo leer los complejos y el amarillo"):
                    st.markdown(
                        "- **Verde:** firma curada completa.\n"
                        "- **Azul:** complejo transferido desde K-12 por "
                        "ortología recíproca; nivel B, no experimental.\n"
                        "- **Amarillo:** se sospecha un complejo incompleto o "
                        "candidato; **no cierra la brecha**.\n"
                        "- **Naranja:** más de un complejo compatible.\n"
                        "- **Gris:** no evaluado.\n\n"
                        "Una brecha solo se asocia a un complejo si una EC/Rhea "
                        "del complejo produce exactamente su metabolito. "
                        "La estructura cuaternaria no se predice: se sospecha."
                    )
                if crs:
                    st.dataframe(pd.DataFrame(crs), width="stretch",
                                 hide_index=True)
                    st.download_button(
                        "⬇️ Complejos (TSV)",
                        pd.DataFrame(crs).to_csv(sep="\t", index=False),
                        file_name=f"{res['genome']}_complexes.tsv",
                        mime="text/plain",
                    )
                tax = res.get("taxonomy")
                if tax:
                    st.markdown("#### 🧭 Taxonomía local")
                    if tax.get("status") != "ok":
                        st.info(f"Taxonomía no disponible: {tax.get('reason','')}")
                    else:
                        trows = tax.get("rows") or []
                        tdf = pd.DataFrame(trows)
                        st.dataframe(tdf, width="stretch", hide_index=True)
                        st.caption(
                            "Normalización de un JSON local ya precalculado; "
                            "no se valida especie por ANI ni red. Un rango "
                            "ausente no es evidencia negativa."
                        )
                        st.download_button(
                            "⬇️ Taxonomía (TSV)",
                            tdf.to_csv(sep="\t", index=False),
                            file_name=f"{res['genome']}.taxonomy.tsv",
                            mime="text/plain",
                        )
        else:
            st.info("Anotación funcional no ejecutada: falta DIAMOND o Swiss-Prot.")
        o = res.get("oric")
        if res.get("oric_error"):
            st.warning(res["oric_error"])
        if o:
            st.markdown("#### 🧭 Origen de replicación")
            c1, c2, c3 = st.columns(3)
            c1.metric("Estado", str(o.get("status", "n/d")))
            c2.metric("Posición", f"{o.get('position'):,}" if o.get("position") else "n/d")
            c3.metric("Confianza", f"{float(o.get('confidence') or 0):.2f}")
            st.caption(f"{o.get('validated_by','')} · {o.get('reason','')}")
            st.download_button(
                "⬇️ oriC (GFF3)", oric_engine.gff_text(o),
                file_name=f"{res['genome']}.oric.gff3", mime="text/plain",
            )
        else:
            st.info("ORCA no ejecutado: falta el modelo o la dependencia orcapy.")

    with tabs[7]:
        mem = res.get("membrane")
        st.caption(
            "Membrana/transporte: hmmscan contra un panel Pfam curado, "
            "hidropatía tipo Kyte-Doolittle, señal péptido y lipobox. "
            "Son heurísticas; no sustituyen TMHMM/SignalP/LipoP."
        )
        if not mem:
            st.info("Membrana no ejecutada. Requiere hmmscan y el panel Pfam.")
        elif not mem.get("ok"):
            st.warning(f"No se pudo ejecutar membrana: {mem.get('reason', '')}")
        else:
            df = pd.DataFrame(mem.get("rows") or [])
            if df.empty:
                st.info("Sin filas de membrana para este genoma.")
            else:
                c1, c2, c3 = st.columns(3)
                c1.metric("Filas", len(df))
                c2.metric("Alta", int((df.get("confidence") == "high").sum())
                          if "confidence" in df else 0)
                c3.metric("Transportadores",
                          int((df.get("type") == "transporter").sum())
                          if "type" in df else 0)
                st.dataframe(df, width="stretch", hide_index=True)
                st.download_button(
                    "⬇️ Membrana (TSV)",
                    df.to_csv(sep="\t", index=False),
                    file_name=f"{res['genome']}_membrane.tsv", mime="text/plain",
                )
                if mem.get("gff"):
                    st.download_button(
                        "⬇️ Membrana (GFF3)", Path(mem["gff"]).read_text(),
                        file_name=f"{res['genome']}_membrane.gff3",
                        mime="text/plain",
                    )

    with tabs[8]:
        st.caption(
            "Operones y promotores **putativos** construidos por reglas sobre "
            "las features de este FASTA. `evidence=rule_based`: son hipótesis, "
            "no evidencia experimental. Un operón se agrupa por mismo "
            "contig/hebra con hueco ≤ 300 pb; el promotor es una ventana fija "
            "de 50 pb upstream y una coincidencia de consenso -35 TTGACA / "
            "-10 TATAAT (espaciador 15-19)."
        )
        reg = res.get("regulation")
        if res.get("regulation_error"):
            st.warning(f"Regulación automática no ejecutada: "
                       f"{res['regulation_error']}")
        elif not reg:
            st.info("Regulación automática no ejecutada para este genoma.")
        elif not reg.get("ok"):
            st.warning(f"No se pudo construir la regulación: "
                       f"{reg.get('reason', '')}")
        else:
            ops = reg.get("rows") or []
            proms = reg.get("promoter_rows") or []
            c1, c2, c3 = st.columns(3)
            c1.metric("Operones", len(ops))
            c2.metric("Promotores", len(proms))
            c3.metric("Con caja -10/-35",
                      sum(1 for p in proms if p.get("box_found") == "yes"))
            st.caption(
                "Cada fila declara `evidence=rule_based`, `confidence` "
                "(low/medium) y `limitation`. Ninguna etiqueta demuestra "
                "cotranscripción ni un TSS; no se usan RegulonDB/DBTBS ni red."
            )
            if ops:
                st.markdown("##### Operones")
                dfo = pd.DataFrame(ops)
                cols = [c for c in ("operon_id", "contig", "strand", "start",
                                    "end", "n_genes", "genes",
                                    "max_internal_gap", "confidence",
                                    "limitation") if c in dfo.columns]
                st.dataframe(dfo[cols], width="stretch", hide_index=True)
                st.download_button(
                    "⬇️ Operones (TSV)",
                    dfo.to_csv(sep="\t", index=False),
                    file_name=f"{res['genome']}.operons_auto.tsv",
                    mime="text/plain", key=f"dl_ops_{res['genome']}",
                )
            else:
                st.info("Sin operones putativos para este genoma.")
            if proms:
                st.markdown("##### Promotores putativos")
                dfp = pd.DataFrame(proms)
                cols = [c for c in ("promoter_id", "operon_id", "contig",
                                    "strand", "gene_id", "start", "end",
                                    "box35", "mm35", "box10", "mm10",
                                    "spacer", "score", "box_found",
                                    "confidence") if c in dfp.columns]
                st.dataframe(dfp[cols], width="stretch", hide_index=True)
                st.download_button(
                    "⬇️ Promotores (TSV)",
                    dfp.to_csv(sep="\t", index=False),
                    file_name=f"{res['genome']}.promoters_auto.tsv",
                    mime="text/plain", key=f"dl_proms_{res['genome']}",
                )
            else:
                st.info("Sin promotores putativos para este genoma.")

    with tabs[9]:
        st.caption(
            "Visor de coordenadas PDB. Puedes subir un PDB/mmCIF o indicar una "
            "ruta local. Una estructura experimental, AFDB o ESMFold no confirma "
            "función, sustrato, interfaz ni cierra brechas."
        )
        use_b = st.checkbox("Colorear por B-factor/RMSX", value=False,
                            key=f"pdb_b_{res['genome']}")
        up = st.file_uploader("Subir PDB/mmCIF", type=["pdb", "ent", "cif"],
                              key=f"pdb_up_{res['genome']}")
        if up is not None:
            structure_viewer.render(
                up.read().decode("utf-8", errors="replace"),
                key=f"pdb_upload_{res['genome']}",
                caption=f"Archivo subido: {up.name}", bfactor=use_b,
            )
        else:
            path_pdb = st.text_input(
                "Ruta local a PDB/mmCIF",
                key=f"pdb_path_{res['genome']}",
                placeholder="results/<genoma>_esmfold/modelo.pdb",
            )
            if path_pdb:
                structure_viewer.render(
                    path_pdb, key=f"pdb_local_{res['genome']}",
                    caption=f"Ruta: {path_pdb}", bfactor=use_b,
                )

    with tabs[10]:
        render_compare()


def main() -> None:
    init_language()
    language_selector()
    st.title(tr("🧬 Anotador de ADN bacteriano"))
    st.caption(T["app_subtitle"])
    render_requirements()
    with st.expander(tr("Anotar por lote (CLI)")):
        st.code("python -m webapp.pipeline genoma.fna --outdir results/ --threads 8\n"
                "python -m webapp.pipeline genoma.fna --outdir results/ \\\n"
                "  --esmfold \"neuros nap\" --esmfold-gaps 3 \\\n"
                "  --structure-evidence --structural-annotation --structural-filters \\\n"
                "  --structure-manifest manifest.tsv")
        st.caption("Misma ruta que la GUI, sin Streamlit: GFF, proteínas, brechas, "
                   "transporte, cierre por estructura, informe y summary.json. "
                   "Ver docs/NEUROSNAP_API.md.")
    with st.expander("¿Cómo funciona?"):
        st.markdown(
            "1. Subes uno o varios FASTA.\n"
            "2. Marcas los análisis y parámetros.\n"
            "3. Pulsas **Anotar**.\n"
            "4. Descargas GFF3, proteínas, candidatos y desiertos, todos etiquetados.\n\n"
            "El flujo es **automático**: no hay decisiones durante la ejecución. "
            "Cada etiqueta declara su fuente y su nivel de confianza."
        )

    # enlace directo: ?modo=rutas&layout=B&genoma=...
    try:
        _m = st.query_params.get("modo")
        if _m in ("navegador", "browser"):
            st.session_state["mode_sel"] = T["mode_browser"]
        elif _m == "rutas":
            st.session_state["mode_sel"] = T["mode_pathways"]
        elif _m in ("hipotesis", "hipótesis", "hypotheses"):
            st.session_state["mode_sel"] = T["mode_hypotheses"]
        elif _m in ("creditos", "créditos"):
            st.session_state["mode_sel"] = T["mode_credits"]
        elif _m in ("reportes", "descargas"):
            st.session_state["mode_sel"] = T["mode_reports"]
        elif _m in ("modulos", "módulos", "modules"):
            st.session_state["mode_sel"] = T["mode_modules"]
        elif _m == "anotar":
            st.session_state["mode_sel"] = T["mode_annotate"]
        if st.query_params.get("genoma"):
            st.session_state["genoma_global"] = st.query_params["genoma"]
    except Exception:
        pass
    mode = st.sidebar.radio(
        T["mode_label"], [T["mode_annotate"], T["mode_pathways"],
                           T["mode_browser"], T["mode_reports"],
                           T["mode_modules"], T["mode_hypotheses"], T["mode_credits"]],
        help=T["mode_help"], key="mode_sel",
    )
    if mode == T["mode_hypotheses"]:
        render_precomputed()
        return
    if mode == T["mode_credits"]:
        credits.render()
        return
    if mode == T["mode_reports"]:
        downloads_ui.render(PROJECT)
        return
    if mode == T["mode_modules"]:
        module_ui.render(DATA_ROOT, OUT_ROOT)
        return
    if mode == T["mode_pathways"]:
        render_pathways()
        return
    if mode == T["mode_browser"]:
        # se ofrece lo que tiene anotación de verdad, no la lista de "free living"
        reg = DATA_ROOT / "regulation"
        opciones = (sorted(d.name for d in reg.iterdir()
                           if d.is_dir() and (d / "genes.tsv").exists())
                    if reg.exists() else [])
        if DEMO_GENOME and "demo" not in opciones:
            opciones = ["demo"] + opciones
        if "ecoli_k12" in opciones:
            opciones = ["ecoli_k12"] + [o for o in opciones if o != "ecoli_k12"]
        if not opciones:
            st.error("No hay genomas del proyecto anotados todavía. "
                     "Ejecuta primero el pipeline de los 400.")
            return
        actual = st.session_state.get("genoma_global") or ""
        if actual not in opciones:
            actual = opciones[0] if opciones else ""
        st.sidebar.selectbox(
            "Genoma activo", opciones,
            index=opciones.index(actual) if actual in opciones else 0,
            key="genoma_global", help=T["help_genoma"],
        )
        elegido = st.session_state.get("genoma_global") or actual
        browser.render(DATA_ROOT, elegido)
        return

    if not (DATA_ROOT / "db" / "swissprot.dmnd").exists():
        st.info(
            "Primera vez: faltan Swiss-Prot y Rhea. Descárgalos una sola vez con "
            "`python -m webapp.setup_data --data-root data`."
        )

    with st.spinner("Detectando herramientas y preparando la interfaz…"):
        tools = tools_cached()
    with st.sidebar:
        st.header("1 · Herramientas disponibles")
        st.caption(T["sidebar_tools"])
        for key, ok in tools.items():
            label = tr(TOOL_LABEL.get(key, key))
            desc = tr(TOOL_DESC.get(key, ""))
            mark = "✅" if ok else "❌"
            # la descripción va aquí, no en un tooltip: en una lista de estado
            # el usuario necesita saber qué hace cada herramienta de un vistazo
            st.markdown(f"{mark} **{label}** — {desc}" if desc else f"{mark} {label}")
        st.caption("Créditos de cada herramienta, con sus autores y enlaces, "
                   "en **🙏 Créditos y herramientas**.")

        st.header("2 · Opciones de anotación")
        st.caption(T["sidebar_options"])
        cds = st.checkbox(
            "CDS (pyrodigal)", value=tools["pyrodigal"],
            disabled=not tools["pyrodigal"],
            help=T["opt_cds"],
        )
        trna = st.checkbox(
            "tRNA (tRNAscan-SE)", value=tools["trna"], disabled=not tools["trna"],
            help=T["opt_trna"],
        )
        rrna = st.checkbox(
            "rRNA (barrnap)", value=tools["rrna"], disabled=not tools["rrna"],
            help=T["opt_rrna"],
        )
        ghost = st.checkbox(
            "ORFs en regiones no anotadas", value=True,
            help=T["opt_ghost"],
        )
        esm = st.checkbox(
            "Puntuar con ESM-2 (sentido biológico)", value=tools["esm"],
            disabled=not tools["esm"],
            help=T["opt_esm"],
        )
        deserts = st.checkbox(
            "Desiertos de anotación", value=True,
            help=T["opt_deserts"],
        )
        compare = st.checkbox(
            "Comparar genomas subidos", value=True,
            help=T["opt_compare"],
        )
        bell = st.checkbox("🔔 Sonar campana al terminar", value=True,
                           help=T["bell"])
        if ghost and not cds:
            st.warning("Los ORFs no anotados necesitan CDS; se activará CDS.")

        st.header("3 · Metabolismo (opcional)")
        st.caption(T["sidebar_metabolism"])
        ref_gcf = st.text_input("GCF de referencia (opcional)", "",
                                help="Si el genoma ya está en el proyecto, reutiliza su red "
                                     "metabólica ya calculada en vez de recalcularla.")
        metabolic = st.checkbox(
            "Construir red metabólica y brechas", value=False,
            disabled=not (tools["diamond"] and (DATA_ROOT / "db" / "swissprot.dmnd").exists()),
            help=T["opt_metabolic"],
        )

        st.header("4 · Parámetros")
        st.caption(T["sidebar_params"])
        min_region_len = st.number_input("Longitud mínima de región (nt)", 30, 5000, 150, 10,
                                         help=T["param_min_region"])
        min_aa = st.number_input("Longitud mínima de ORF (aa)", 20, 500, 60, 5,
                                 help=T["param_min_aa"])
        esm_top = st.number_input("Puntuar con ESM los N mejores", 0, 2000, 100, 10,
                                  help=T["param_esm_top"])

    st.subheader("Genomas")
    uploaded = st.file_uploader(
        T["upload_label"],
        type=["fa", "fna", "fasta", "fas", "gz"],
        accept_multiple_files=True,
        help=T["upload_help"],
    )

    demo_paths = [p for p in (Path(__file__).parent / "samples" / "demo.fna",
                              Path(__file__).parent / "samples" / "demo1b.fna") if p.exists()]
    if demo_paths and st.button("▶️ Probar con 2 genomas de ejemplo", help=T["demo"]):
        st.session_state["demo"] = True
    if DEMO_GENOME and st.button(
            "🧬 Explorar la demo del navegador",
            help="Abre un genoma de demostración ya anotado; no descarga nada."):
        st.query_params["modo"] = "navegador"
        st.query_params["genoma"] = "demo"
        st.rerun()
    use_demo = bool(st.session_state.get("demo", False)) and bool(demo_paths)
    if use_demo:
        st.caption(f"{len(demo_paths)} genoma(s) de ejemplo cargados (misma especie). Pulsa Anotar.")

    if st.button(T["run_button"], type="primary",
                 disabled=not (uploaded or use_demo), help=T["run_help"]):
        inputs = []
        if use_demo:
            inputs += [(p.name, p.read_bytes()) for p in demo_paths]
        for f in uploaded or []:
            inputs.append((f.name, f.read()))

        opts = {
            "cds": cds or ghost,
            "trna": trna,
            "rrna": rrna,
            "ghost": ghost,
            "esm": esm,
            "esm_top": int(esm_top),
            "deserts": deserts,
            "min_region_len": int(min_region_len),
            "min_aa": int(min_aa),
        }
        results = []
        comparison = {}
        with st.status("Anotando genomas…", expanded=True) as status:
            for i, (name, data) in enumerate(inputs, 1):
                t0 = time.time()
                status.write(f"🧬 {name}: leyendo y aplicando herramientas ({i}/{len(inputs)})…")
                try:
                    res = annotate(name, data, opts)
                    res["elapsed"] = round(time.time() - t0, 1)
                    if res.get("features"):
                        status.write("📐 Operones y promotores por reglas…")
                        try:
                            res["regulation"] = run_regulation(
                                res["features"], OUT_ROOT, res["genome"],
                                sequences=res.get("sequences"))
                        except Exception as exc:  # noqa: BLE001
                            res["regulation_error"] = str(exc)
                    if res.get("proteins"):
                        status.write("🧾 Anotación funcional (DIAMOND/Swiss-Prot)…")
                        try:
                            res["functional"] = annotate_functional(
                                res["proteins"], DATA_ROOT, OUT_ROOT, res["genome"],
                                threads=4)
                            comp = infer_complexes(
                                res["features"],
                                res["functional"].get("hits") or {},
                                res["genome"])
                            model = run_complex_model(
                                res["features"],
                                res["functional"].get("hits") or {},
                                OUT_ROOT, res["genome"], DATA_ROOT)
                            if model.get("ok"):
                                comp["complexes"] += model["complexes"]
                                comp["by_protein"].update(model["by_protein"])
                                comp["n_candidate"] = model.get("n_candidate", 0)
                            transfer = load_complex_transfer(
                                DATA_ROOT, res["genome"], res["features"])
                            if transfer.get("ok"):
                                comp["complexes"] += transfer["complexes"]
                                comp["by_protein"].update(transfer["by_protein"])
                                comp["n_candidate"] = (comp.get("n_candidate", 0)
                                                       + transfer.get("n_candidate", 0))
                            enrich_complex_rows(
                                res["functional"].get("rows") or [], comp)
                            res["complexes"] = comp
                        except Exception as exc:  # noqa: BLE001
                            res["functional_error"] = str(exc)
                    if res.get("sequences"):
                        status.write("🧭 Predicción de origen (ORCA + GC skew)…")
                        try:
                            res["oric"] = oric_engine.run(res, DATA_ROOT)
                        except Exception as exc:  # noqa: BLE001
                            res["oric_error"] = str(exc)
                    try:
                        res["report_md"] = build_report(res)
                    except Exception:
                        pass
                    results.append(res)
                    status.write(
                        f"✅ {name}: {res['coverage']*100:.1f}% anotado · "
                        f"{len(res['candidates'])} ORFs · {len(res.get('deserts', []))} desiertos"
                    )
                    if metabolic and res.get("proteins"):
                        status.write("🧪 Construyendo red metabólica…")
                        try:
                            res["metabolic"] = build_network(
                                res["proteins"], DATA_ROOT,
                                ref_gcf.strip() or None, threads=4,
                                panel_hmm=DATA_ROOT / "db" / "pfam" / "panel_focus.hmm",
                                transport_db=DATA_ROOT / "db" / "transport.dmnd",
                                producer_db=DATA_ROOT / "db" / "producer_family.dmnd",
                            )
                            res["report_md"] = build_report(res)
                            status.write(
                                f"   red: {res['metabolic']['reactions']} reacciones · "
                                f"{res['metabolic']['n_gaps']} brechas"
                            )
                        except Exception as exc:  # noqa: BLE001
                            res["metabolic_error"] = str(exc)
                    elif metabolic:
                        res["metabolic_error"] = "Activa CDS para tener proteínas."
                except Exception as exc:  # noqa: BLE001
                    st.error(f"{name}: {exc}")
            if compare and len(results) >= 2:
                status.write("🔗 Comparando genomas subidos (95%)…")
                try:
                    comparison = compare_results(results)
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"No se pudo comparar: {exc}")
            status.update(label="Anotación completa", state="complete")
        st.session_state["results"] = results
        st.session_state["comparison"] = comparison
        st.session_state["bell_played"] = False

    results = st.session_state.get("results", [])
    comparison = st.session_state.get("comparison", {})
    if results:
        st.success(f"{len(results)} archivo(s) anotado(s).")
        if comparison:
            m1, m2, m3 = st.columns(3)
            m1.metric("Clusters entre genomas", comparison.get("clusters", 0))
            m2.metric("Clusters conservados", comparison.get("multi_genome", 0))
            m3.metric("Genomas comparados", len(results))
        if bell and not st.session_state.get("bell_played"):
            st.audio(bell_wav(), format="audio/wav", autoplay=True)
            st.session_state["bell_played"] = True
            st.caption("🔔 Análisis terminado.")
        for res in results:
            render_result(res)
            st.divider()

    st.caption(
        f"🧬 {tr('Anotador')} v{APP_VERSION} · "
        "👤 **Camilo Rojas Bravo** · "
        + tr("Ingeniero en Biotecnología (USACH)") + " · "
        "Camilo.Rojas@usach.cl  \n"
        "🤖 " + tr("Asistencia de desarrollo:") + " "
        + tr("**DeepSeek Harness (DSH)** con DeepSeek,") + " "
        + tr("noviembre de 2026") + ". "
        + tr("Los datos y las decisiones científicas son del autor.")
    )


if __name__ == "__main__":
    main()
