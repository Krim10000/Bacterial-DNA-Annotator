"""Vista de aristas TF–promotor curadas (RegulonDB/PRODORIC)."""
from __future__ import annotations
import csv
from pathlib import Path
import pandas as pd
import streamlit as st
from i18n import tr


def _load(data_root: Path) -> pd.DataFrame:
    p=data_root/"regulatory"/"evidence_edges.tsv"
    if not p.exists(): return pd.DataFrame()
    try: return pd.read_csv(p,sep="\t",dtype=str).fillna("")
    except Exception: return pd.DataFrame()


def render(data_root: Path, genome: str | None = None) -> None:
    st.subheader("🧩 Aristas TF–promotor")
    st.caption(tr("Solo aristas curadas/experimentales de RegulonDB y PRODORIC. No son inferencias."))
    df=_load(data_root)
    if df.empty:
        st.info(tr("No hay tabla de aristas. Ejecuta scripts/regulatory/02_build_edges.py."))
        return
    c1,c2,c3=st.columns(3)
    c1.metric(tr("Aristas"),len(df))
    c2.metric(tr("Genomas"),df["genome"].nunique())
    c3.metric(tr("Fuentes"),df["source_db"].nunique())
    genomes=sorted(df["genome"].unique())
    pick=st.multiselect(tr("Genoma"),genomes,default=[genome] if genome in genomes else genomes[:3])
    src=st.multiselect(tr("Fuente"),sorted(df["source_db"].unique()),default=sorted(df["source_db"].unique()))
    conf=st.multiselect(tr("Confianza"),sorted(df["confidence"].unique()),default=sorted(df["confidence"].unique()))
    sub=df[df["genome"].isin(pick) & df["source_db"].isin(src) & df["confidence"].isin(conf)]
    st.dataframe(sub,width="stretch",hide_index=True)
    st.download_button(tr("⬇️ Descargar aristas (TSV)"),sub.to_csv(sep="\t",index=False),
                       file_name="evidence_edges.tsv",mime="text/tab-separated-values")
    st.caption(tr("Fuente y confianza se muestran por fila; no se promedian."))
