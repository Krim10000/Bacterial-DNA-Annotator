"""UI compartida para el cierre de brechas con ESMFold."""
from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

import credits
import structure_viewer

from structure_evidence import evidence_for_rows
from structural_annotation import scan_pdbs, summarize as summarize_enzymm
from structural_filters import run_filters
from esmfold_tools import (
    LOCAL_LABEL,
    PROVIDER_LABEL,
    capability_markdown,
    capability_rows,
    detect_capability,
    has_key,
    local_command,
    sequence_for_gap,
    submit_to_neuros_nap,
)

NO = "No usar ESMFold"
LOCAL = f"{LOCAL_LABEL} (usa ESMFold)"
API = f"{PROVIDER_LABEL} (usa ESMFold)"


def render_esmfold_closure(gap: dict[str, Any], proteins: list[dict[str, Any]] | None = None,
                           key_prefix: str = "") -> None:
    """Muestra alternativas y, si se marca una, avisa que usa ESMFold."""
    gap = dict(gap)
    st.markdown("#### Cierre por estructura")
    st.caption(
        "Las alternativas de estructura de esta sección **usan ESMFold**. "
        "ESMFold predice una estructura 3D a partir de la secuencia. "
        "Una estructura de alta confianza **no confirma** función, sustrato ni cierre de brecha."
    )
    choice = st.radio(
        "Alternativa",
        [NO, LOCAL, API],
        key=f"{key_prefix}_esmfold_choice",
        help="Las dos alternativas de estructura usan ESMFold; se indica explícitamente.",
    )
    if choice == NO:
        st.caption("No se usará estructura para esta brecha.")
        return

    st.warning(
        f"Has marcado **{choice}**. Esta alternativa **usa ESMFold** para predecir la estructura. "
        "No es evidencia experimental ni cierra la brecha por sí sola."
    )

    cap = detect_capability()
    st.markdown("##### ¿Puede correr ESMFold en esta máquina?")
    st.caption("Se evalúa la máquina donde corre la GUI/pipeline, no el navegador.")
    st.dataframe(pd.DataFrame(capability_rows(cap)), width="stretch", hide_index=True)
    st.markdown(capability_markdown(cap))

    if choice == LOCAL:
        if not cap["can_run"]:
            st.error("Esta máquina no cumple los requisitos para ESMFold local.")
            st.info("Usa la alternativa **Neurosnap API** o una máquina con GPU ≥16 GB, RAM ≥32 GB y ≥20 GB libres.")
        else:
            st.code(local_command("candidatos.faa", "esmfold_local"), language="bash")
            st.caption("Local no envía secuencias a la nube; requiere pesos y OpenFold instalados.")
        return

    # API Neurosnap
    st.markdown("##### Neurosnap API")
    if has_key():
        st.success("`ESM_API_KEY` detectada en el entorno/.env.")
    else:
        st.warning("Falta `ESM_API_KEY`. Añádela a `.env`; no la pegues en el chat.")
        st.code(
            "printf 'ESM_API_KEY=%s\\n' 'TU_KEY' >> .env\n"
            "chmod 600 .env\n"
            "set -a; . ./.env; set +a",
            language="bash",
        )
    st.markdown(
        "- Coste: **~0,04 créditos** por ESMFold de ~300 aa; **$1/crédito**.\n"
        "- Consume crédito aunque la estructura salga con baja confianza.\n"
        "- Requiere autorización para enviar la secuencia a un servicio externo."
    )

    name, seq = sequence_for_gap(gap, proteins or [])
    if not seq:
        st.info("No hay secuencia de ORF asociada a esta brecha; usa el FASTA de candidatos.")
        st.code(
            "set -a; . ./.env; set +a\n"
            ".venv/bin/python scripts/structure/05_neurosnap_fold.py \\\n"
            "  --fasta candidatos.faa --out esmfold_neuros_nap --limit 1",
            language="bash",
        )
        return

    st.caption(f"ORF: `{name}` · {len(seq)} aa")
    if st.button("Enviar este ORF a Neurosnap (consume crédito)", key=f"{key_prefix}_ns_submit"):
        if not has_key():
            st.error("Falta ESM_API_KEY.")
            return
        with st.spinner("Enviando y esperando ESMFold..."):
            res = submit_to_neuros_nap(seq, name, outdir=st.session_state.get("_esmfold_out", "esmfold_neuros_nap"))
        if res.get("error"):
            st.error(res["error"])
        else:
            st.success(f"Job {res['job_id']} · {res['status']} · pLDDT {res.get('plddt')} · pTM {res.get('ptm')}")
            pdb_path = res.get("pdb", "")
            st.caption(f"PDB: {pdb_path or 'no descargado'}")
            st.session_state["_last_esmfold"] = res
            if pdb_path:
                structure_viewer.render(
                    pdb_path, key=f"{key_prefix}_esmfold_3d",
                    caption="ESMFold: estructura predicha; no confirmación.",
                )


def render_structure_evidence(gap: dict, genome: str, pdb: str, uniprot: str = "",
                              key_prefix: str = "") -> None:
    """Botón explícito para cruzar Foldseek + AlphaFill + sitios activos."""
    if not pdb:
        return
    st.markdown("##### Evidencia multi-fuente")
    st.caption("Foldseek PDB100 + AlphaFill + sitios activos locales. Consume tiempo de APIs; no cierra la brecha.")
    if st.button("Cruzar Foldseek + AlphaFill", key=f"{key_prefix}_evidence"):
        with st.spinner("Consultando Foldseek/AlphaFill..."):
            rows = evidence_for_rows([{
                "genome": genome, "gap": gap.get("gap", ""), "query": gap.get("query", ""),
                "pdb": pdb, "uniprot": uniprot,
            }])
        if not rows:
            st.warning("Sin resultados; revisa PDB/ruta.")
            return
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.caption("`multi_source_support` no es cierre; `closure=no` siempre.")
        st.session_state["_last_structure_evidence"] = rows
    structure_viewer.render(
        pdb, key=f"{key_prefix}_evidence_3d",
        caption=f"PDB: {pdb}",
    )


def render_structural_annotation(pdb: str, key_prefix: str = "") -> None:
    """EnzyMM/PyJess: plantillas catalíticas M-CSA; no es confirmación."""
    if not pdb:
        return
    st.markdown("##### Anotación automática EnzyMM/PyJess")
    st.caption(
        "Busca arreglos geométricos de residuos catalíticos de M-CSA. "
        "Un match es una hipótesis de plantilla, no una función confirmada."
    )
    if st.button("Escanear con EnzyMM/PyJess", key=f"{key_prefix}_enzymm"):
        with st.spinner("Escaneando plantillas M-CSA..."):
            res = scan_pdbs([pdb], "data/metabolism/targets/structures/enzymm_gui", jobs=4)
        rows = summarize_enzymm(res.get("hits") or [])
        if res.get("error"):
            st.error(res["error"])
        elif rows:
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            st.warning("Match de plantilla geométrica: no confirma reacción, sustrato ni cierre.")
        else:
            st.info("Sin plantilla M-CSA para este PDB. `no_template` no significa ausencia de función.")
        st.caption(f"Source: EnzyMM/PyJess + M-CSA · `{pdb}`")


def render_five_filters(pdb: str, genome: str = "", uniprot: str = "",
                        key_prefix: str = "") -> None:
    """Ejecuta los 5 filtros estructurales y muestra prioridad, no cierre."""
    if not pdb:
        return
    st.markdown("##### 5 filtros estructurales")
    st.caption(
        "F1 Foldseek · F2 EnzyMM/M-CSA · F3 P2Rank · F4 AlphaFill · F5 sitio activo. "
        "Reducen candidatos; no confirman función."
    )
    if st.button("Ejecutar 5 filtros", key=f"{key_prefix}_five"):
        with st.spinner("F1–F5: Foldseek, EnzyMM, P2Rank, AlphaFill, sitio activo..."):
            r = run_filters(pdb, genome=genome, uniprot=uniprot)
        st.dataframe(pd.DataFrame([r]), use_container_width=True, hide_index=True)
        st.markdown(f"**Prioridad:** `{r.get('priority')}` · **Pasan:** {r.get('n_pass')}/5")
        st.caption("Un filtro que no pasa no elimina por sí solo; `no_template` no es ausencia.")
        st.caption(credits.line("foldseek", "enzymm", "p2rank", "alphafill", "mcsa", "afdb"))
        st.session_state["_last_five_filters"] = r
        structure_viewer.render(
            pdb, key=f"{key_prefix}_five_3d", caption=f"PDB: {pdb}",
        )
