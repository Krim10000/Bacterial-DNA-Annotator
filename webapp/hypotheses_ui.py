"""GUI para revisar hipótesis estructurales sin confundirlas con anotación."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

import credits

DECISIONS = ["pending", "reviewed_orientative", "reviewed_supported",
             "conflict", "rejected", "needs_validation"]


def _decisions(genome: str) -> dict:
    return st.session_state.setdefault("review_decisions", {}).setdefault(genome, {})


def render_hypotheses(res: dict) -> None:
    rows = res.get("structural_hypotheses") or []
    st.caption(
        "Hipótesis computacionales a partir de Foldseek, EnzyMM/M-CSA, P2Rank y AlphaFill. "
        "No sustituyen la anotación base ni cierran brechas."
    )
    if not rows:
        st.info("Sin hipótesis estructurales. Ejecuta los 5 filtros o aporta un manifest de estructuras.")
        return
    df = pd.DataFrame(rows)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Hipótesis", len(df))
    c2.metric("Con ≥2 fuentes", int((df["n_sources"].astype(float) >= 2).sum()))
    c3.metric("Supported", int((df["confidence"] == "supported_hypothesis").sum()))
    c4.metric("Conflicto", int((df["confidence"] == "conflict").sum()))
    cols = [c for c in ("protein_id", "hypothesis", "ec", "cofactor", "pocket",
                        "template_mcsa", "foldseek_top", "n_sources", "confidence",
                        "review_decision") if c in df.columns]
    st.dataframe(df[cols], width="stretch", hide_index=True)
    st.download_button("⬇️ Hipótesis (TSV)", df.to_csv(sep="\t", index=False),
                       file_name=f"{res['genome']}_structural_hypotheses.tsv",
                       mime="text/plain")
    st.download_button("⬇️ Evidencia (JSON)",
                       json.dumps(rows, indent=2, ensure_ascii=False),
                       file_name=f"{res['genome']}_structural_hypotheses.json",
                       mime="application/json")

    key = st.selectbox("Revisar proteína", df["protein_id"].tolist(),
                       key=f"hyp_{res['genome']}")
    row = df[df["protein_id"] == key].iloc[0].to_dict()
    st.markdown(f"### {key} · `{row.get('hypothesis','')}`")
    st.caption(f"UniProt: `{row.get('uniprot','') or '—'}` · PDB: `{row.get('pdb','') or '—'}`")
    a, b, c = st.columns(3)
    a.metric("Confianza", row.get("confidence", ""))
    b.metric("Fuentes", row.get("n_sources", ""))
    c.metric("EC", row.get("ec", "") or "—")

    with st.expander("F1 · Foldseek", expanded=True):
        st.write(f"Top: **{row.get('foldseek_top','') or '—'}**")
        st.json(json.loads(row.get("evidence_json") or "{}").get("row", {}).get("foldseek_top", ""))
    with st.expander("F2 · EnzyMM / M-CSA"):
        st.write(f"Plantillas: `{row.get('template_mcsa','') or '—'}`")
        st.write(f"EC: `{row.get('ec','') or '—'}`")
    with st.expander("F3 · P2Rank"):
        st.write(f"Bolsillos: **{row.get('pocket','') or '—'}**")
    with st.expander("F4 · AlphaFill"):
        st.write(f"Ligandos/cofactores: `{row.get('cofactor','') or '—'}`")
    with st.expander("F5 · Sitio activo"):
        ev = json.loads(row.get("evidence_json") or "{}").get("row", {})
        st.write(f"Posiciones: `{ev.get('active_positions','') or '—'}`")
    with st.expander("Fuente y límites"):
        st.write(row.get("limitation", ""))
        st.caption(credits.line("foldseek", "enzymm", "p2rank", "alphafill", "mcsa", "afdb"))

    dec = _decisions(res["genome"]).get(key, row.get("review_decision", "pending"))
    new = st.selectbox("Decisión de revisión", DECISIONS,
                       index=DECISIONS.index(dec) if dec in DECISIONS else 0,
                       key=f"dec_{res['genome']}_{key}")
    if new != dec:
        _decisions(res["genome"])[key] = new
        st.success(f"Decisión guardada en la sesión: {new}")
    dec_rows = [{"genome": res["genome"], "protein_id": k, "decision": v}
                for k, v in _decisions(res["genome"]).items()]
    if dec_rows:
        st.download_button("⬇️ Decisiones de revisión (TSV)",
                           pd.DataFrame(dec_rows).to_csv(sep="\t", index=False),
                           file_name=f"{res['genome']}_review_decisions.tsv",
                           mime="text/plain")

    pdb = Path(str(row.get("pdb", "")))
    if pdb.exists():
        st.download_button("⬇️ PDB del candidato", pdb.read_bytes(),
                           file_name=pdb.name, mime="chemical/x-pdb")
