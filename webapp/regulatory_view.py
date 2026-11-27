"""Vista de aristas TF–promotor y catálogo ENTRAF."""
from __future__ import annotations
import csv
from pathlib import Path
import pandas as pd
import streamlit as st
from i18n import tr


def _read(data_root: Path, rel: str) -> pd.DataFrame:
    p=data_root/rel
    if not p.exists(): return pd.DataFrame()
    try: return pd.read_csv(p,sep="\t",dtype=str).fillna("")
    except Exception: return pd.DataFrame()


def _edges(data_root: Path, genome: str | None) -> None:
    st.caption(tr("Aristas TF–promotor curadas de RegulonDB y PRODORIC. No son inferencias."))
    df=_read(data_root,"regulatory/evidence_edges.tsv")
    if df.empty:
        st.info(tr("No hay tabla de aristas. Ejecuta scripts/regulatory/02_build_edges.py.")); return
    c1,c2,c3=st.columns(3); c1.metric(tr("Aristas"),len(df)); c2.metric(tr("Genomas"),df["genome"].nunique()); c3.metric(tr("Fuentes"),df["source_db"].nunique())
    genomes=sorted(df["genome"].unique()); pick=st.multiselect(tr("Genoma"),genomes,default=[genome] if genome in genomes else genomes[:3])
    src=st.multiselect(tr("Fuente"),sorted(df["source_db"].unique()),default=sorted(df["source_db"].unique()))
    conf=st.multiselect(tr("Confianza"),sorted(df["confidence"].unique()),default=sorted(df["confidence"].unique()))
    sub=df[df["genome"].isin(pick) & df["source_db"].isin(src) & df["confidence"].isin(conf)]
    st.dataframe(sub,width="stretch",hide_index=True)
    st.download_button(tr("⬇️ Descargar aristas (TSV)"),sub.to_csv(sep="\t",index=False),file_name="evidence_edges.tsv",mime="text/tab-separated-values")


def _entraf(data_root: Path, genome: str | None) -> None:
    st.caption(tr("Catálogo ENTRAF: factores de transcripción curados, rol, familia, PDB y PubMed. Fuente: IIMAS-UNAM."))
    df=_read(data_root,"regulatory/entraf_entries.tsv")
    if df.empty:
        st.info(tr("No hay catálogo ENTRAF. Ejecuta scripts/regulatory/03_scrape_entraf.py.")); return
    c1,c2,c3=st.columns(3); c1.metric(tr("Entradas"),len(df)); c2.metric(tr("Mapeadas"),int((df["hit_type"]=="direct_species").sum())); c3.metric(tr("Con PubMed"),int((df["pubmed"]!="").sum()))
    genomes=sorted({g for s in df["genome_match"] for g in str(s).split(";") if g})
    pick=st.multiselect(tr("Genoma"),genomes,default=[genome] if genome in genomes else genomes[:3])
    roles=st.multiselect(tr("Rol"),sorted(df["role"].unique()),default=sorted(df["role"].unique()))
    sub=df[df["role"].isin(roles) & (df["genome_match"].apply(lambda s: any(g in str(s).split(";") for g in pick)) if pick else True)]
    st.dataframe(sub,width="stretch",hide_index=True)
    st.download_button(tr("⬇️ Descargar ENTRAF (TSV)"),sub.to_csv(sep="\t",index=False),file_name="entraf_entries.tsv",mime="text/tab-separated-values")
    st.caption("ENTRAF: An Encyclopedia of well-annotated DNA-binding TRanscription factors in Bacteria and Archaea. Computational Biology Lab, IIMAS-UNAM Sede Mérida.")


def render(data_root: Path, genome: str | None = None) -> None:
    st.subheader("🧩 Regulación")
    tab1,tab2=st.tabs([tr("Aristas TF–promotor"),tr("Catálogo ENTRAF")])
    with tab1: _edges(data_root,genome)
    with tab2: _entraf(data_root,genome)
