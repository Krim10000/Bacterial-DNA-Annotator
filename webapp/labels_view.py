"""Vista de etiquetas de proceso y GO MF/CC por genoma (solo lectura)."""
from __future__ import annotations
import csv
from collections import Counter
from pathlib import Path
import pandas as pd
import streamlit as st
from i18n import tr


def _find(data_root: Path, genome: str) -> Path | None:
    for root in (data_root/"labels"/genome, data_root.parent/"results"/"labels"/genome):
        p=root/"operon_labels.tsv"
        if p.exists(): return p
    return None


def _process(data_root: Path, genome: str) -> None:
    path=_find(data_root,genome)
    if path is None:
        st.info(tr("No hay etiquetas de proceso para este genoma. Ejecuta el módulo «Etiquetas de proceso» (orden 32)."))
        return
    with path.open(encoding="utf-8",errors="replace") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    labels=Counter(l for r in rows for l in (r.get("labels") or "").split(";") if l)
    c1,c2,c3=st.columns(3)
    c1.metric(tr("Operones"),len(rows))
    c2.metric(tr("Con etiqueta"),sum(1 for r in rows if r.get("labels")))
    c3.metric(tr("Etiquetas únicas"),len(labels))
    st.dataframe(pd.DataFrame(labels.most_common(),columns=[tr("Etiqueta"),tr("Operones")]),
                 width="stretch",hide_index=True)
    with st.expander(tr("Ver tabla por operón")):
        st.dataframe(pd.DataFrame(rows),width="stretch",hide_index=True)
    st.download_button(tr("⬇️ Descargar etiquetas (TSV)"),path.read_text(encoding="utf-8"),
                       file_name=f"{genome}_operon_labels.tsv",mime="text/tab-separated-values")
    st.caption(tr("Fuente: GO, KEGG, EC y léxico local. Si falta una fuente, se declara no_data; no se inventa."))


def _go(data_root: Path, genome: str) -> None:
    st.subheader("🧬 GO MF/CC por accesión")
    st.caption(tr("Función molecular y componente celular; fuente GAF/UniProt explícita. No es prueba experimental."))
    for root in (data_root/"labels"/"go_mf_cc"/genome, data_root.parent/"results"/"go_mf_cc"/genome):
        go_path=root/"gene_go.tsv"
        if not go_path.exists(): continue
        go=pd.read_csv(go_path,sep="\t",dtype=str).fillna("")
        if go.empty:
            st.info(tr("Sin términos GO MF/CC para este genoma.")); return
        c1,c2,c3=st.columns(3)
        c1.metric(tr("Genes anotados"),go["gene_id"].nunique())
        c2.metric(tr("Términos MF"),go.loc[go["aspect"]=="F","go_id"].nunique())
        c3.metric(tr("Términos CC"),go.loc[go["aspect"]=="C","go_id"].nunique())
        src=go.groupby("source").size().sort_values(ascending=False).head(10)
        st.caption(tr("Fuente explícita por GAF: ")+", ".join(f"{k} ({v})" for k,v in src.items()))
        st.dataframe(go.head(500),width="stretch",hide_index=True)
        st.download_button(tr("⬇️ Descargar GO MF/CC (TSV)"),go_path.read_text(encoding="utf-8"),
                           file_name=f"{genome}_gene_go.tsv",mime="text/tab-separated-values")
        return
    st.info(tr("Sin GO MF/CC. Ejecuta el módulo «GO MF/CC» (orden 31)."))


def render(data_root: Path, genome: str) -> None:
    st.subheader("🏷️ Etiquetas y GO")
    st.caption(tr("Evidencia contextual por gen/operón; fuente y límites explícitos."))
    _process(data_root, genome)
    st.divider()
    _go(data_root, genome)
