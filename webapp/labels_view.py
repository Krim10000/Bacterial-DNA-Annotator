"""Vista de etiquetas de proceso por genoma (solo lectura)."""
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


def render(data_root: Path, genome: str) -> None:
    st.subheader("🏷️ Etiquetas de proceso")
    st.caption(tr("Etiquetas GO/KEGG/EC/léxico por genoma. Son evidencia contextual, no prueba de función ni de sigma."))
    path=_find(data_root,genome)
    if path is None:
        st.info(tr("No hay etiquetas para este genoma. Ejecuta el módulo «Etiquetas de proceso» (orden 32)."))
        return
    rows=[]
    with path.open(encoding="utf-8",errors="replace") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    labels=Counter(l for r in rows for l in (r.get("labels") or "").split(";") if l)
    c1,c2,c3=st.columns(3)
    c1.metric(tr("Operones"),len(rows))
    c2.metric(tr("Con etiqueta"),sum(1 for r in rows if r.get("labels")))
    c3.metric(tr("Etiquetas únicas"),len(labels))
    st.markdown("**"+tr("Etiquetas más frecuentes")+"**")
    st.dataframe(pd.DataFrame(labels.most_common(),columns=[tr("Etiqueta"),tr("Operones")]),
                 width="stretch",hide_index=True)
    with st.expander(tr("Ver tabla por operón")):
        df=pd.DataFrame(rows)
        st.dataframe(df,width="stretch",hide_index=True)
    st.download_button(tr("⬇️ Descargar etiquetas (TSV)"),path.read_text(encoding="utf-8"),
                       file_name=f"{genome}_operon_labels.tsv",mime="text/tab-separated-values")
    st.caption(tr("Fuente: GO, KEGG, EC y léxico local. Si falta una fuente, se declara no_data; no se inventa."))
