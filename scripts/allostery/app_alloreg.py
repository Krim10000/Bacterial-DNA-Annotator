#!/usr/bin/env python3
"""GUI del módulo: (1) genoma crudo → represores, (2) preparar estructura.

    streamlit run scripts/allostery/app_alloreg.py

Página 1 — Genoma: sube un FASTA, corre el pipeline (genes → dominios → homólogo →
efector → sitios) y muestra tabla, sitios, tiempos y RAM.

Página 2 — Preparar: sube un PDB, protona a pH 7 y minimiza con OpenMM **en local**
(sin GPU, sin Modal). Tarda minutos; se muestran tiempos medidos.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[2]
PIPE = ROOT / "scripts" / "allostery" / "18_genome_pipeline.py"
PREP = ROOT / "scripts" / "allostery" / "25_prepare_local.py"
CONTACT = ROOT / "scripts" / "allostery" / "26_contact_score.py"
MMPP = ROOT / "scripts" / "allostery" / "27_interaction_mm.py"
MMPY = ROOT / ".mamba_root" / "envs" / ".mamba_mm312" / "bin" / "python"
PY = str(ROOT / ".venv" / "bin" / "python")

st.set_page_config(page_title="AlloReg — genomas y estructuras", layout="wide")
st.title("AlloReg — represores en genomas y preparación de estructuras")
st.caption("Prototipo. Toda salida lleva `status`: `inherited_homolog`, "
           "`family_prior`, `regprecise_tf_name` o `unknown`. Nada es `confirmed`.")

pagina = st.sidebar.radio("Qué quieres hacer",
                          ["Genoma → represores", "Preparar estructura"])


def pagina_genoma():
    with st.sidebar:
        up = st.file_uploader("Genoma (FASTA)", type=["fa", "fna", "fasta", "txt"])
        max_tfs = st.slider("Máximo de candidatos", 5, 200, 40)
        solo_rep = st.checkbox("Buscar solo posibles represores", value=True)
        nivel = st.radio("Nivel de anotación", ["rapida", "estandar"], index=1,
                         help="rapida: sin sitios PWM (más rápido) · "
                              "estandar: incluye sitios en −250 pb")
        run = st.button("Ejecutar", type="primary", disabled=up is None)
        st.markdown("---")
        st.markdown("**Coste**: CPU local. ~22 s por genoma de 2 Mb; "
                    "~60-85 s para *E. coli* (4,6 Mb). Sin GPU.")

    if run and up is not None:
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / up.name
            fa.write_bytes(up.getbuffer())
            out = Path(tmp) / "run"
            cmd = [PY, str(PIPE), str(fa), "--out", str(out),
                   "--max-tfs", str(max_tfs), "--tier", nivel]
            if solo_rep:
                cmd.append("--find-repressors")
            with st.spinner("Corriendo pipeline…"):
                p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode != 0:
                st.error("El pipeline falló")
                st.code(p.stderr[-3000:])
                st.stop()
            summary = json.loads((out / "summary.json").read_text())
            df = pd.read_csv(out / "represores.tsv", sep="\t")
            st.session_state["genoma"] = {
                "summary": summary, "df": df.to_dict("records"),
                "report": (out / "informe.md").read_text(),
                "tsv": df.to_csv(sep="\t", index=False)}

    res = st.session_state.get("genoma")
    if not res:
        st.info("Sube un FASTA y pulsa **Ejecutar**. Para probar, un fragmento de "
                "100-200 kb tarda segundos.")
        return

    s = res["summary"]
    r = s.get("recursos", {})
    t = r.get("tiempos_s", {})
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Genes", s["genes"])
    c2.metric("Candidatos TF", s["tf_candidatos"])
    c3.metric("Represores", s["reportados"])
    c4.metric("Tiempo total", f"{t.get('total_s', '—')} s")
    c5.metric("RAM pico", f"{r.get('pico_ram_mb', '—')} MB")

    df = pd.DataFrame(res["df"])
    if not df.empty:
        st.subheader("Represores candidatos")
        st.dataframe(df.drop(columns=["sitios"]), use_container_width=True,
                     hide_index=True)
        pick = st.selectbox("Sitios de un TF", df["gen"].tolist())
        row = df[df["gen"] == pick].iloc[0]
        sites = json.loads(row["sitios"]) if isinstance(row["sitios"], str) else row["sitios"]
        st.table(pd.DataFrame(sites)) if sites else st.info("Sin sitios en este TF.")

    st.download_button("Descargar TSV", res["tsv"], file_name="represores.tsv")
    with st.expander("Informe completo"):
        st.markdown(res["report"])
    with st.expander("Recursos y pendientes"):
        st.json({"recursos": r, "profunda_bajo_demanda": s.get("profunda_bajo_demanda", [])})


def pagina_prep():
    with st.sidebar:
        pdb = st.file_uploader("Estructura (PDB)", type=["pdb", "ent"])
        lig = st.file_uploader("Ligando (SDF, opcional)", type=["sdf", "mol"],
                               help="Si lo subes, se calcula la puntuación de contacto "
                                    "en el bolsillo (no es energía).")
        ph = st.number_input("pH", 4.0, 10.0, 7.0, 0.5)
        iters = st.slider("Pasos de minimización", 0, 500, 100,
                          help="0 = solo protonar. 100 ≈ 3-6 min en este equipo; "
                               "500 ≈ 8-15 min (medido: 1,9 s/paso para 5.832 átomos).")
        run = st.button("Preparar", type="primary", disabled=pdb is None)
        st.markdown("---")
        st.markdown("**Coste**: $0. OpenMM en CPU local, sin GPU ni Modal.")
        st.markdown("**Límite**: la dinámica no es viable aquí "
                    "(medido 0,25 ns/día). Esto solo protona y minimiza.")

    if run and pdb is not None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / pdb.name
            src.write_bytes(pdb.getbuffer())
            out = Path(tmp) / (src.stem + "_preparado.pdb")
            with st.spinner("Protonando y minimizando… (puede tardar minutos)"):
                p = subprocess.run([PY, str(PREP), str(src), "--out", str(out),
                                    "--ph", str(ph), "--max-iter", str(iters)],
                                   capture_output=True, text=True)
            if p.returncode != 0:
                st.error("La preparación falló")
                st.code(p.stderr[-3000:])
                st.stop()
            res = json.loads(p.stdout)
            contacto = None
            if lig is not None:
                ligp = Path(tmp) / lig.name
                ligp.write_bytes(lig.getbuffer())
                q = subprocess.run([PY, str(CONTACT), str(out), str(ligp), "--decoy"],
                                   capture_output=True, text=True)
                if q.returncode == 0:
                    contacto = json.loads(q.stdout)
                if MMPY.exists():
                    import os
                    env = dict(os.environ)
                    env["PATH"] = str(MMPY.parent) + ":" + env.get("PATH", "")
                    m = subprocess.run([str(MMPY), str(MMPP), str(out), str(ligp),
                                        "--steps", "20"],
                                       capture_output=True, text=True, env=env)
                    if m.returncode == 0:
                        contacto["mm"] = json.loads(m.stdout)
            st.session_state["prep"] = {**res, "pdb": out.read_text(),
                                        "contacto": contacto}

    res = st.session_state.get("prep")
    if not res:
        st.info("Sube un PDB para protonar y minimizar. Funciona con receptores de "
                "proteína; los ligandos requieren parámetros GAFF/OpenFF aparte.")
        return

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Átomos (con H)", res["atomos"])
    c2.metric("Energía inicial", f"{res['energia_inicial_kj']:.0f}")
    c3.metric("Energía final", f"{res['energia_final_kj']:.0f}")
    c4.metric("Tiempo", f"{res['t_minimizacion_s']:.0f} s")
    st.download_button("Descargar PDB preparado", res["pdb"],
                       file_name="receptor_preparado.pdb")

    cont = res.get("contacto")
    if cont:
        st.subheader("Puntuación de contacto en el bolsillo")
        st.caption("Contactos geométricos, **no** energía. La columna *decoy* es la "
                   "misma molécula desplazada 25 Å: si la pose real no se entierra "
                   "más que el decoy, no hay bolsillo.")
        filas = [{"métrica": k, "pose": v, "decoy": cont["decoy_desplazado"].get(k)}
                 for k, v in cont["pose"].items()]
        st.dataframe(pd.DataFrame(filas), use_container_width=True, hide_index=True)
        st.success("La pose se separa del decoy" if cont.get("separa")
                   else "La pose NO se separa del decoy")
        mm = cont.get("mm")
        if mm:
            st.subheader("Energía de interacción MM (OpenFF)")
            c1, c2, c3 = st.columns(3)
            c1.metric("E interacción", f"{mm['e_interaccion_kj']:.0f} kJ/mol")
            c2.metric("LJ", f"{mm['e_lj_kj']:.0f}")
            c3.metric("Coulomb", f"{mm['e_coulomb_kj']:.0f}")
            st.caption("En **vacío**: el término Coulomb no está apantallado, así que "
                       "no sirve para ordenar ligandos. Es un paso intermedio hacia "
                       "MM/GBSA (con solvente implícito).")
        else:
            st.info("Energía MM no disponible: falta el entorno conda "
                    "(`.mamba_root/envs/.mamba_mm312`).")
    with st.expander("Detalle técnico"):
        st.json(res)


if pagina.startswith("Genoma"):
    pagina_genoma()
else:
    pagina_prep()
