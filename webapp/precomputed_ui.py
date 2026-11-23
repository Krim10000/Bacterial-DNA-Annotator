"""Vista GUI de hipótesis precalculadas."""
import streamlit as st

from hypotheses_ui import render_hypotheses
from precomputed_hypotheses import load


def render() -> None:
    st.subheader("🔬 Hipótesis estructurales precalculadas")
    st.caption(
        "1.192 accesiones con AFDB, 25.058 clusters y 176 hipótesis con tres fuentes. "
        "Son orientativas; no cierran brechas."
    )
    rows = load()
    if not rows:
        st.info("No está la tabla precalculada. Ejecuta UNCHAR_AFDB_PRESENT_filters.tsv.")
        return
    render_hypotheses({"genome": "UNCHAR_AFDB_PRESENT", "structural_hypotheses": rows})
