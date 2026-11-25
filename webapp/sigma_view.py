"""Vista **Sigma** del navegador: métricas, filtros, tabla y evidencia.

Se apoya en las funciones puras de ``browser.py`` (carga defensiva, color y
confianza) y en Streamlit solo para pintar la interfaz. No escribe datos ni
recalcula anotaciones: lee ``data/sigma/<genoma>/sigma_promoters.tsv``.

Uso desde ``app.py`` (una línea dentro de la vista de genoma)::

    import sigma_view
    sigma_view.render(DATA_ROOT, genome)

Los tests no arrancan Streamlit; ``summary``, ``filter_promoters`` y
``table_view`` son funciones puras y se pueden probar con DataFrames.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from browser import (
    format_sigma_json,
    load_sigma_promoters,
    sigma_color,
    sigma_confidence,
    sigma_short,
)

#: Etiquetas humanas de los niveles de confianza (``confidence``).
CONFIDENCE_ES = {
    "high": "fuerte",
    "medium": "media",
    "low": "baja",
    "unknown": "desconocida",
}

#: Etiquetas humanas de ``confidence_label``.
CONFIDENCE_LABEL_ES = {
    "sospecha_fuerte": "sospecha fuerte",
    "sospecha_media": "sospecha media",
    "sospecha_baja": "sospecha baja",
    "sin_sospecha": "sin sospecha",
}

#: Orden de las columnas de la tabla de pantalla.
TABLE_COLUMNS = (
    ("operon_label", "Operón"),
    ("first_gene_name", "Primer gen"),
    ("sigma", "σ"),
    ("sigma_class", "Clase"),
    ("sigma_arch", "Arquitectura"),
    ("spacer", "Espaciador"),
    ("confidence_label", "Confianza"),
    ("confidence", "Nivel"),
    ("sigma_status", "Estado σ"),
    ("status", "Estado"),
    ("score", "Score"),
    ("pvalue", "p-valor"),
    ("evidence", "Evidencia"),
    ("detected_by", "Detectado por"),
    ("start", "Inicio"),
    ("end", "Fin"),
)


def _empty_summary() -> dict:
    return {
        "total": 0, "assigned": 0, "ambiguous": 0, "unresolved": 0,
        "strong": 0, "medium": 0, "low": 0, "unknown": 0,
        "by_type": [], "confidence_labels": {},
    }


def _tiers(df: pd.DataFrame) -> pd.Series:
    if "confidence_tier" in df.columns:
        return df["confidence_tier"].astype(str)
    return df.apply(sigma_confidence, axis=1)


def summary(df: pd.DataFrame | None) -> dict:
    """Resumen de anotaciones sigma, sin Streamlit.

    - ``assigned`` cuenta ``sigma_status == "assigned"`` (asignación
      respaldada); ``ambiguous`` y ``unresolved`` se informan aparte.
    - ``strong``/``medium``/``low`` cuentan los niveles ``high``/``medium``/
      ``low`` de ``confidence`` (con fallback a ``confidence_label``).
    - ``by_type`` ordena por número de filas y resume cada tipo de σ.
    """
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return _empty_summary()
    status = (df["sigma_status"].astype(str).str.lower()
              if "sigma_status" in df.columns else pd.Series("", index=df.index))
    tiers = _tiers(df)
    out = {
        "total": int(len(df)),
        "assigned": int((status == "assigned").sum()),
        "ambiguous": int((status == "ambiguous").sum()),
        "unresolved": int((status == "unresolved").sum()),
        "strong": int((tiers == "high").sum()),
        "medium": int((tiers == "medium").sum()),
        "low": int((tiers == "low").sum()),
        "unknown": int((tiers == "unknown").sum()),
        "confidence_labels": {},
        "by_type": [],
    }
    if "confidence_label" in df.columns:
        counts = df["confidence_label"].astype(str).value_counts()
        out["confidence_labels"] = {str(k): int(v) for k, v in counts.items()}
    sigma_values = (df["sigma"].astype(str) if "sigma" in df.columns
                    else pd.Series("", index=df.index))
    types: list[dict] = []
    for name in sorted(set(sigma_values)):
        mask = sigma_values == name
        sub_tiers = tiers[mask]
        sub_status = status[mask]
        types.append({
            "sigma": name,
            "etiqueta": sigma_short(name),
            "total": int(mask.sum()),
            "assigned": int((sub_status == "assigned").sum()),
            "strong": int((sub_tiers == "high").sum()),
            "medium": int((sub_tiers == "medium").sum()),
            "low": int((sub_tiers == "low").sum()),
            "color": sigma_color(name),
        })
    types.sort(key=lambda item: (-item["total"], item["sigma"]))
    out["by_type"] = types
    return out


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = [value]
    return [str(v) for v in value if str(v).strip()]


def filter_promoters(df: pd.DataFrame | None, *, sigma=None,
                      confidence_label=None,
                      query: str = "") -> pd.DataFrame:
    """Filtra por tipo de σ, ``confidence_label`` y búsqueda de texto.

    ``sigma`` y ``confidence_label`` aceptan un valor o una lista; una lista
    vacía / ``None`` significa «todas». La búsqueda es por subcadena, sin
    distinguir mayúsculas, en operón, primer gen, etiqueta y tipo.
    """
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame() if df is None else df.copy()
    out = df
    sigmas = _as_list(sigma)
    if sigmas and "sigma" in out.columns:
        wanted = {s.strip().lower() for s in sigmas}
        out = out[out["sigma"].astype(str).str.lower().isin(wanted)]
    labels = _as_list(confidence_label)
    if labels and "confidence_label" in out.columns:
        wanted = {s.strip().lower() for s in labels}
        out = out[out["confidence_label"].astype(str).str.lower().isin(wanted)]
    q = str(query or "").strip()
    if q and not out.empty:
        cols = [c for c in ("operon_id", "operon_label", "first_gene",
                            "first_gene_name", "sigma", "sigma_class",
                            "evidence", "detected_by") if c in out.columns]
        if cols:
            mask = pd.Series(False, index=out.index)
            for col in cols:
                mask |= out[col].astype(str).str.contains(
                    q, case=False, na=False, regex=False)
            out = out[mask]
    return out.reset_index(drop=True)


def table_view(df: pd.DataFrame | None) -> pd.DataFrame:
    """Vista legible y estable de la tabla de promotores sigma."""
    columns = [dst for _, dst in TABLE_COLUMNS]
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame({col: pd.Series(dtype="object")
                             for col in columns})
    view = pd.DataFrame(index=df.index)
    for src, dst in TABLE_COLUMNS:
        view[dst] = df[src] if src in df.columns else ""
    if "σ" in view.columns:
        view["σ"] = view["σ"].map(sigma_short)
    if "Confianza" in view.columns and "confidence_label" in df.columns:
        view["Confianza"] = df["confidence_label"].map(
            lambda v: CONFIDENCE_LABEL_ES.get(str(v), str(v)))
    if "Nivel" in view.columns:
        view["Nivel"] = df.apply(sigma_confidence, axis=1).map(
            lambda v: CONFIDENCE_ES.get(v, v))
    return view.reset_index(drop=True)


def _safe_slug(value: str) -> str:
    text = "".join(ch if (ch.isalnum() or ch in "-_.") else "_"
                   for ch in str(value))
    return text.strip("._")[:60] or "genoma"


def _coord(value) -> str:
    text = str(value).strip()
    return "—" if text.lower() in ("", "nan", "none", "<na>") else text


def _json_records(df: pd.DataFrame) -> str:
    records = df.where(pd.notna(df), "").to_dict(orient="records")
    return json.dumps(records, ensure_ascii=False, indent=2, default=str)


def _detail(row: pd.Series) -> None:
    """Detalle de una fila con la evidencia que exige el contrato sigma."""
    sigma_name = str(row.get("sigma") or "")
    label = (str(row.get("operon_label") or "")
             or str(row.get("first_gene_name") or "")
             or str(row.get("first_gene") or row.get("operon_id") or ""))
    tier = sigma_confidence(row)
    st.markdown(f"### 🧪 {sigma_short(sigma_name)} · {label}")
    st.caption(
        "Anotación de promotor dependiente de σ. Es una **hipótesis "
        "priorizada**, no una confirmación experimental."
    )
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Confianza",
              f"{row.get('confidence_label') or '—'} ({tier})")
    c2.metric("Score", str(row.get("score") or "—"))
    c3.metric("p-valor", str(row.get("pvalue") or "—"))
    c4.metric("Estado σ", str(row.get("sigma_status") or "—"))
    st.markdown(
        f"**Tipo:** {row.get('sigma_class') or '—'} · "
        f"**Arquitectura:** {row.get('sigma_arch') or '—'} · "
        f"**Espaciador:** {row.get('spacer') or '—'} pb · "
        f"**Margen:** {row.get('score_margin') or '—'} · "
        f"**z:** {row.get('zscore') or '—'}"
    )
    st.markdown(
        f"**Cajas:** aguas arriba `{row.get('upstream_box') or '—'}` · "
        f"aguas abajo `{row.get('downstream_box') or '—'}` · "
        f"**dist. al inicio:** {row.get('dist_to_start') or '—'} pb"
    )
    st.markdown(
        f"**Posición dibujada:** {_coord(row.get('start'))}–"
        f"{_coord(row.get('end'))} pb (start del operón; no es un TSS)."
    )
    st.markdown("**Evidencia (`evidence`)**")
    st.code(str(row.get("evidence") or "—"), language="text")
    st.markdown("**Detectado por (`detected_by`)**")
    st.code(str(row.get("detected_by") or "—"), language="text")
    st.markdown(f"**Paso del pipeline:** `{row.get('pipeline_step') or '—'}`")
    with st.expander("Alternativas y procedencia (JSON)"):
        st.markdown("**Alternativas (`alternatives`)**")
        st.code(format_sigma_json(row.get("alternatives")) or "—",
                language="json")
        st.markdown("**Procedencia (`provenance`)**")
        st.code(format_sigma_json(row.get("provenance")) or "—",
                language="json")
        for key in ("pathway_hits", "expression_detail", "inventory_detail"):
            value = str(row.get(key) or "")
            if value:
                st.caption(f"{key}: {value}")


def _render_table(filtered: pd.DataFrame) -> None:
    view = table_view(filtered)
    key = "sigma_view_table"
    selected: int | None = None
    try:
        event = st.dataframe(
            view, width="stretch", hide_index=True, on_select="rerun",
            selection_mode="single-row", key=key,
        )
        rows = (event["selection"]["rows"] if isinstance(event, dict)
                else getattr(event, "selection", {}).get("rows", []))
        if rows:
            selected = int(rows[0])
    except Exception:  # noqa: BLE001 - Streamlit antiguo sin selection
        idx = st.selectbox(
            "Elegir fila", range(len(filtered)), key=key + "_sb",
            format_func=lambda i: " · ".join(
                str(view.iloc[i].get(c, ""))[:50]
                for c in ("σ", "Operón", "Primer gen")),
        )
        selected = int(idx)
    if selected is not None:
        _detail(filtered.iloc[selected])
    else:
        st.caption("Haz clic en una fila para ver su evidencia completa.")


def render(data_root: Path, genome: str) -> None:
    """Pinta la vista Sigma del genoma: métricas, filtros, tabla y detalle."""
    st.subheader("🧪 Sigma y promotores")
    st.caption(
        "Promotores dependientes de factor σ anotados por el módulo `sigma`. "
        "La tabla muestra el **tipo**, la **confianza** y la **evidencia** de "
        "cada operón; no se recalcula nada aquí."
    )
    df = load_sigma_promoters(Path(data_root), genome)
    if df.empty:
        st.warning(df.attrs.get("sigma_warning")
                   or f"No hay anotaciones sigma para {genome}.")
        st.caption(
            "Genera las salidas con el módulo **Sigmas y promotores** "
            "(o `python scripts/sigma/10_annotate_sigma.py --genome "
            f"{genome}`). Esta vista no se rompe si el archivo falta."
        )
        return

    stats = summary(df)
    st.caption(
        f"**{stats['total']:,}** operones anotados · "
        f"**{stats['assigned']:,} asignados** "
        f"(`sigma_status=assigned`) · "
        f"ambiguos {stats['ambiguous']:,} · sin resolver "
        f"{stats['unresolved']:,}. La confianza es del nivel `confidence`; "
        "`sospecha_fuerte` prioriza la hipótesis, **no** la confirma."
    )
    st.caption(
        "**Gold K-12 (2026-11-25):** top-1 57,0 % frente a una mayoría "
        "trivial de 62,6 %; top-3 90,2 %. La confianza ordena las hipótesis, "
        "pero **no es una probabilidad calibrada** ni una validación "
        "experimental. Ver `docs/SIGMA_BENCHMARK.md`."
    )
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Anotaciones", f"{stats['total']:,}")
    m2.metric("Asignados", f"{stats['assigned']:,}")
    m3.metric("Fuerte", f"{stats['strong']:,}")
    m4.metric("Media", f"{stats['medium']:,}")
    m5.metric("Baja", f"{stats['low']:,}")

    st.markdown("#### Por tipo de sigma")
    by_type = pd.DataFrame(stats["by_type"])
    if by_type.empty:
        st.caption("Sin tipos de σ en la tabla.")
    else:
        st.dataframe(
            by_type.rename(columns={
                "sigma": "σ", "etiqueta": "Etiqueta", "total": "Total",
                "assigned": "Asignados", "strong": "Fuerte",
                "medium": "Media", "low": "Baja", "color": "Color",
            })[["σ", "Etiqueta", "Total", "Asignados", "Fuerte", "Media",
                "Baja"]],
            width="stretch", hide_index=True,
        )
        st.caption(
            "«Asignados» = `sigma_status=assigned`. Los recuentos de "
            "confianza incluyen todas las filas del tipo, también las "
            "ambiguas y sin resolver."
        )

    st.markdown("### Filtros y tabla")
    sigma_options = sorted({str(v) for v in df.get("sigma", []) if str(v)})
    label_options = [k for k in CONFIDENCE_LABEL_ES
                     if "confidence_label" in df.columns
                     and (df["confidence_label"].astype(str) == k).any()]
    f1, f2, f3 = st.columns([1.4, 1.4, 2.2])
    sigma_sel = f1.multiselect(
        "Tipo de σ", sigma_options, format_func=sigma_short,
        help="Vacío = todas. El color de la pista sigue el tipo.")
    conf_sel = f2.multiselect(
        "Confianza (`confidence_label`)", label_options,
        format_func=lambda v: CONFIDENCE_LABEL_ES.get(v, v),
        help="Filtra por la etiqueta de sospecha que escribe el módulo.")
    query = f3.text_input(
        "Buscar por operón o gen", placeholder="p. ej. pyrD, thrLABC, rpoS")
    filtered = filter_promoters(
        df, sigma=sigma_sel, confidence_label=conf_sel, query=query)
    st.caption(
        f"**{len(filtered):,} de {len(df):,}** anotaciones coinciden con los "
        "filtros. La descarga incluye exactamente esta selección."
    )

    _slug = _safe_slug(genome)
    d1, d2 = st.columns(2)
    d1.download_button(
        "⬇️ Descargar TSV filtrado",
        data=filtered.to_csv(sep="\t", index=False).encode("utf-8"),
        file_name=f"{_slug}_sigma_promoters_filtrado.tsv",
        mime="text/tab-separated-values",
        key="sigma_dl_tsv",
    )
    d2.download_button(
        "⬇️ Descargar JSON filtrado",
        data=_json_records(filtered).encode("utf-8"),
        file_name=f"{_slug}_sigma_promoters_filtrado.json",
        mime="application/json",
        key="sigma_dl_json",
    )

    if filtered.empty:
        st.info("Ninguna anotación coincide con los filtros.")
        return
    _render_table(filtered)

    st.markdown("---")
    st.caption(
        "**Límites.** La posición dibujada usa el `start` del operón porque "
        "`sigma_promoters.tsv` no trae coordenadas propias del promotor; no "
        "es un TSS. La anotación depende del escaneo de PWMs y del fondo "
        "empírico del módulo, no de una medición experimental. "
        "`expression_support=not_available` significa que no había tabla de "
        "expresión, no que no se exprese."
    )
