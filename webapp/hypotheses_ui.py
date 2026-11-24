"""GUI para revisar hipótesis estructurales sin confundirlas con anotación.

La pestaña reúne en una sola vista: métricas de confianza explicadas, guía
in situ de F1–F5, tabla priorizada, selector de proteína, ficha con los cinco
filtros, visor 3D cuando hay PDB y descargas de tabla y estructura.

Ninguna etiqueta de confianza confirma función, sustrato ni cierre de brecha:
eso solo lo hace un experimento. Por eso las explicaciones viven aquí y no se
delegan a la documentación.
"""
from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

import credits
import structure_viewer

DECISIONS = ["pending", "reviewed_orientative", "reviewed_supported",
             "conflict", "rejected", "needs_validation"]

#: Columnas visibles de la tabla y del TSV descargable.
COLUMNS = ["protein_id", "hypothesis", "ec", "cofactor", "pocket",
           "template_mcsa", "foldseek_top", "n_sources", "confidence",
           "review_decision"]

_FLAG_FIELDS = ("f1_foldseek", "f2_enzymm", "f3_pocket", "f4_alphafill",
                "f5_active_site")

#: Orden de prioridad: primero lo apoyado, luego lo que exige revisión.
CONFIDENCE_ORDER = {
    "supported_hypothesis": 0,
    "conflict": 1,
    "high_priority": 2,
    "orientative": 3,
    "structural_only": 4,
    "uncharacterized": 5,
}

#: Explicación in situ de cada etiqueta de confianza (sin remitir a docs).
CONFIDENCE_INFO = {
    "supported_hypothesis": {
        "emoji": "🟢",
        "label": "Hipótesis apoyada por tres fuentes",
        "rule": "F2 EnzyMM/M-CSA + F3 P2Rank + F4 AlphaFill",
        "meaning": (
            "La plantilla enzimática (M-CSA), el bolsillo geométrico (P2Rank) "
            "y el ligando/cofactor transferido (AlphaFill) apuntan al mismo "
            "tipo de función."
        ),
        "caveat": (
            "Tres señales computacionales no son un ensayo: no confirma "
            "sustrato, reacción ni cierre de brecha."
        ),
    },
    "high_priority": {
        "emoji": "🟡",
        "label": "Alta prioridad para revisar",
        "rule": "F2 EnzyMM/M-CSA + F5 sitio activo",
        "meaning": (
            "Hay plantilla M-CSA y posiciones de sitio activo curadas que "
            "permiten revisar la hipótesis con residuos concretos."
        ),
        "caveat": (
            "El sitio activo es anotación curada transferida, no una medición "
            "hecha en esta proteína."
        ),
    },
    "orientative": {
        "emoji": "🔵",
        "label": "Orientativa",
        "rule": "Solo F2 EnzyMM/M-CSA",
        "meaning": (
            "La plantilla M-CSA sugiere una familia enzimática, pero faltan "
            "bolsillo, ligando o sitio activo independientes."
        ),
        "caveat": (
            "Una plantilla parecida puede compartir pliegue sin compartir la "
            "misma reacción."
        ),
    },
    "structural_only": {
        "emoji": "🟣",
        "label": "Solo evidencia estructural",
        "rule": "F1 Foldseek + (F3 P2Rank o F4 AlphaFill)",
        "meaning": (
            "Hay homología estructural y un bolsillo o un ligando compatible, "
            "pero no plantilla M-CSA con EC."
        ),
        "caveat": (
            "La geometría no identifica el sustrato: no permite nombrar la "
            "función."
        ),
    },
    "uncharacterized": {
        "emoji": "⚪",
        "label": "Sin señal estructural",
        "rule": "Sin fuentes F1–F5",
        "meaning": (
            "Los filtros aplicados no encontraron señal suficiente; la "
            "proteína queda sin hipótesis estructural."
        ),
        "caveat": (
            "La ausencia de señal no demuestra que la proteína no tenga "
            "función."
        ),
    },
    "conflict": {
        "emoji": "🟠",
        "label": "Conflicto por revisar",
        "rule": "La fila o la revisión declararon conflicto",
        "meaning": (
            "La evidencia estructural contradice la anotación base, o dos "
            "fuentes discrepan; hay que revisar antes de usar la hipótesis."
        ),
        "caveat": (
            "Un conflicto no se resuelve solo: no cambia la anotación base ni "
            "cierra la brecha."
        ),
    },
    "unknown": {
        "emoji": "⚫",
        "label": "Sin clasificar",
        "rule": "Etiqueta de confianza no reconocida",
        "meaning": (
            "La fila no trae una etiqueta de confianza válida; se muestra el "
            "valor original sin interpretarlo."
        ),
        "caveat": (
            "No se inventa una confianza: la hipótesis queda como no "
            "clasificada."
        ),
    },
}

#: Guía in situ de los cinco filtros estructurales.
FIELDS_GUIDE = [
    {
        "Filtro": "F1",
        "Fuente": "Foldseek",
        "Qué mide": "Homología estructural (pliegue) contra PDB100/AFDB.",
        "Qué no prueba": "Que la función sea la misma; pliegues parecidos "
                         "pueden tener reacciones distintas.",
    },
    {
        "Filtro": "F2",
        "Fuente": "EnzyMM / M-CSA",
        "Qué mide": "Plantilla enzimática M-CSA y EC asociada.",
        "Qué no prueba": "La EC es de la plantilla, no de esta proteína.",
    },
    {
        "Filtro": "F3",
        "Fuente": "P2Rank",
        "Qué mide": "Bolsillos y cavidades predichos sobre la superficie.",
        "Qué no prueba": "Qué ligando entra ni si hay catálisis.",
    },
    {
        "Filtro": "F4",
        "Fuente": "AlphaFill",
        "Qué mide": "Ligandos y cofactores transferidos desde homólogos.",
        "Qué no prueba": "Que ese cofactor se una de verdad a esta proteína.",
    },
    {
        "Filtro": "F5",
        "Fuente": "UniProt / estructura",
        "Qué mide": "Posiciones de sitio activo curadas, mapeadas al modelo.",
        "Qué no prueba": "Es anotación transferida, no una medición "
                         "experimental.",
    },
]

_UNPARSED = object()


def _text(value: Any) -> str:
    """Convierte un valor a texto corto sin lanzar excepciones."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    return str(value).strip()


def _json_text(value: Any) -> str:
    """Serializa un valor para mostrarlo como texto, sin fallar."""
    if isinstance(value, str):
        return value
    try:
        return json.dumps(value, ensure_ascii=False)
    except (TypeError, ValueError):
        return str(value)


def _as_int(value: Any) -> int:
    """Entero tolerante a ``"3"``, ``"3.0"``, ``None`` o basura."""
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


def is_on(value: Any) -> bool:
    """Interpreta booleanos de pipeline y de TSV (``True``/``"true"``/``1``)."""
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return value != 0
    return str(value).strip().lower() in ("true", "1", "yes", "si", "sí", "y")


def safe_foldseek_detail(row: Mapping[str, Any] | None) -> dict[str, Any]:
    """Resume el detalle F1 sin lanzar excepciones y sin confundir texto con JSON.

    ``evidence_json`` puede ser un ``dict``, JSON válido, texto plano
    (por ejemplo ``q9plu2``), una lista o basura. El resultado siempre trae
    ``detail`` como ``dict`` y ``text`` como ``str``: la GUI nunca pasa un
    string suelto a ``st.json``.

    Claves de salida:

    ``kind``
        ``structured``, ``text``, ``scalar``, ``top_only``, ``empty`` o
        ``invalid``.
    ``top``
        Homólogo Foldseek de la fila o del detalle.
    ``detail``
        Diccionario con la evidencia estructurada (vacío si no la hay).
    ``text``
        Texto plano que sí conviene mostrar cuando no hay estructura.
    ``message``
        Explicación breve para la interfaz.
    ``has_detail``
        ``True`` solo si ``detail`` tiene contenido.
    """
    out: dict[str, Any] = {
        "kind": "empty",
        "top": "",
        "detail": {},
        "text": "",
        "message": "",
        "has_detail": False,
    }
    if not isinstance(row, Mapping):
        out["kind"] = "invalid"
        out["message"] = (
            "Fila de hipótesis no válida: no hay detalle de Foldseek."
        )
        return out

    top = _text(row.get("foldseek_top"))
    raw = row.get("evidence_json")
    evidence: Mapping[str, Any] | None = None
    parsed: Any = _UNPARSED

    if isinstance(raw, Mapping):
        evidence = raw
    elif isinstance(raw, str):
        clean = raw.strip()
        if clean:
            try:
                parsed = json.loads(clean)
            except (TypeError, ValueError, RecursionError):
                parsed = _UNPARSED
            if parsed is _UNPARSED:
                out["kind"] = "text"
                out["text"] = clean
                out["message"] = (
                    "evidence_json es texto plano: se muestra tal cual y sin "
                    "detalle estructurado."
                )
            elif isinstance(parsed, Mapping):
                evidence = parsed
            elif isinstance(parsed, str):
                out["kind"] = "text"
                out["text"] = parsed
                out["message"] = (
                    "evidence_json contiene un texto, no un objeto: sin "
                    "detalle estructurado."
                )
            else:
                out["kind"] = "scalar"
                out["text"] = _json_text(parsed)
                out["message"] = (
                    "evidence_json no es un objeto: se muestra el valor tal "
                    "cual, sin detalle estructurado."
                )
    elif raw is not None:
        out["kind"] = "scalar"
        out["text"] = _json_text(raw)
        out["message"] = (
            "evidence_json no es un objeto: se muestra el valor tal cual, "
            "sin detalle estructurado."
        )

    if evidence is not None:
        if not top:
            top = _text(evidence.get("foldseek_top"))
        candidate = evidence.get("row")
        if isinstance(candidate, Mapping) and candidate:
            out["detail"] = dict(candidate)
        elif candidate is not None and not isinstance(candidate, Mapping):
            if not out["text"]:
                out["text"] = _text(candidate)
        if not out["detail"]:
            rest = {k: v for k, v in evidence.items() if k != "row"}
            if rest:
                out["detail"] = rest

    if (out["text"] and not out["detail"]
            and out["kind"] not in ("text", "scalar", "invalid")):
        out["kind"] = "text"
        out["message"] = (
            "La evidencia es texto plano: se muestra tal cual y sin detalle "
            "estructurado."
        )

    if out["detail"]:
        out["has_detail"] = True
        out["kind"] = "structured"
        out["message"] = (
            "Detalle estructurado de la evidencia; incluye la fila de filtros."
        )
        if not top:
            top = _text(out["detail"].get("foldseek_top"))
    elif out["kind"] not in ("text", "scalar", "invalid"):
        if top:
            out["kind"] = "top_only"
            out["message"] = (
                "Sin evidencia estructurada; solo se conoce el homólogo "
                "Foldseek."
            )
        else:
            out["kind"] = "empty"
            out["message"] = "Sin detalle de Foldseek para esta hipótesis."

    out["top"] = top
    return out


def confidence_info(value: Any) -> dict[str, Any]:
    """Devuelve la explicación de una etiqueta de confianza.

    Normaliza mayúsculas y espacios. Una etiqueta desconocida cae en
    ``unknown`` y conserva el valor original en ``raw``; nunca lanza.
    """
    raw = "" if value is None else str(value).strip()
    key = raw.lower()
    info = CONFIDENCE_INFO.get(key)
    if info is None:
        key = "unknown"
        info = CONFIDENCE_INFO["unknown"]
    out: dict[str, Any] = {"key": key, "raw": raw}
    out.update(info)
    return out


def classify_confidence(row: Mapping[str, Any] | None) -> str:
    """Clasifica la confianza de una fila con la misma regla del pipeline.

    F2 + F3 + F4 → ``supported_hypothesis``; F2 + F5 → ``high_priority``;
    F2 solo → ``orientative``; F1 + (F3 o F4) → ``structural_only``; si no,
    ``uncharacterized``. Un conflicto explícito de la fila o de la revisión
    tiene prioridad y devuelve ``conflict``.
    """
    if not isinstance(row, Mapping):
        return "uncharacterized"
    stored = _text(row.get("confidence")).lower()
    decision = _text(row.get("review_decision")).lower()
    if (stored == "conflict" or decision == "conflict"
            or is_on(row.get("conflict"))):
        return "conflict"
    f1, f2, f3, f4, f5 = (is_on(row.get(k)) for k in _FLAG_FIELDS)
    if f2 and f3 and f4:
        return "supported_hypothesis"
    if f2 and f5:
        return "high_priority"
    if f2:
        return "orientative"
    if f1 and (f3 or f4):
        return "structural_only"
    return "uncharacterized"


def confidence_label(row: Mapping[str, Any] | None) -> str:
    """Etiqueta de confianza de la fila, calculada si falta o no se reconoce."""
    if not isinstance(row, Mapping):
        return "uncharacterized"
    stored = _text(row.get("confidence")).lower()
    if stored in CONFIDENCE_INFO and stored != "unknown":
        return stored
    return classify_confidence(row)


def prioritize_rows(rows: Iterable[Any]) -> list[dict[str, Any]]:
    """Copia las filas-dict y las ordena por confianza y número de fuentes.

    No muta la entrada. Las filas que no son ``dict`` se descartan para que la
    interfaz no se rompa con datos parciales.
    """
    clean = [dict(row) for row in rows if isinstance(row, Mapping)]
    return sorted(
        clean,
        key=lambda row: (
            CONFIDENCE_ORDER.get(confidence_label(row), 99),
            -_as_int(row.get("n_sources")),
            _text(row.get("protein_id")),
            _text(row.get("hypothesis")),
        ),
    )


def _guide_rows() -> list[dict[str, str]]:
    return [dict(item) for item in FIELDS_GUIDE]


def _confidence_rows() -> list[dict[str, str]]:
    rows = []
    for key in ("supported_hypothesis", "high_priority", "conflict",
                "orientative", "structural_only", "uncharacterized"):
        info = CONFIDENCE_INFO[key]
        rows.append({
            "Confianza": f"{info['emoji']} {info['label']}",
            "Etiqueta": key,
            "Regla": info["rule"],
            "Significado": info["meaning"],
            "Límite": info["caveat"],
        })
    return rows


def _decisions(genome: str) -> dict:
    return st.session_state.setdefault("review_decisions", {}).setdefault(genome, {})


def _scalar(value: Any) -> Any:
    """Deja escalares como están y convierte objetos anidados a texto."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return _text(value)


def _frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = ""
        df[col] = [_scalar(value) for value in df[col]]
    return df


def _option_label(row: Mapping[str, Any]) -> str:
    pid = _text(row.get("protein_id")) or "—"
    hypothesis = _text(row.get("hypothesis")) or "sin hipótesis"
    label = confidence_info(confidence_label(row))["label"]
    return f"{pid} · {hypothesis} · {label}"


def _pdb_path(value: Any) -> Path | None:
    text = _text(value)
    if not text:
        return None
    path = Path(text)
    try:
        if path.is_file():
            return path
    except OSError:
        return None
    return None


def render_hypotheses(res: dict) -> None:
    """Dibuja la vista conjunta de hipótesis estructurales de un genoma."""
    genome = _text(res.get("genome")) or "genoma"
    raw_rows = res.get("structural_hypotheses") or []
    if isinstance(raw_rows, Mapping):
        raw_rows = [raw_rows]
    elif not isinstance(raw_rows, (list, tuple)):
        raw_rows = []

    st.caption(
        "Hipótesis computacionales a partir de Foldseek, EnzyMM/M-CSA, P2Rank "
        "y AlphaFill. No sustituyen la anotación base ni cierran brechas."
    )
    if not raw_rows:
        st.info(
            "Sin hipótesis estructurales. Ejecuta los 5 filtros o aporta un "
            "manifest de estructuras."
        )
        return

    ordered = prioritize_rows(raw_rows)
    if not ordered:
        st.warning(
            "Las hipótesis no tienen formato de fila (dict): no hay tabla que "
            "mostrar. Se conserva la descarga JSON de la evidencia original."
        )
        st.download_button(
            "⬇️ Evidencia (JSON)",
            json.dumps(raw_rows, indent=2, ensure_ascii=False, default=str),
            file_name=f"{genome}_structural_hypotheses.json",
            mime="application/json",
            key=f"hyp_json_raw_{genome}",
        )
        return

    df = _frame(ordered)
    labels = [confidence_label(row) for row in ordered]
    n_supported = sum(1 for label in labels if label == "supported_hypothesis")
    n_high = sum(1 for label in labels if label == "high_priority")
    n_conflict = sum(1 for label in labels if label == "conflict")
    n_two_sources = sum(1 for row in ordered
                        if _as_int(row.get("n_sources")) >= 2)

    st.markdown("##### Métricas de confianza")
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric(
        "Hipótesis", len(ordered),
        help="Filas revisables con al menos un filtro aplicado.",
    )
    c2.metric(
        "Con ≥2 fuentes", n_two_sources,
        help="Hipótesis donde coinciden dos o más de los filtros F1–F5. "
             "Coincidir no confirma función.",
    )
    c3.metric(
        "Apoyadas", n_supported,
        help="Etiqueta supported_hypothesis: F2 + F3 + F4. Sigue siendo una "
             "hipótesis, no un resultado experimental.",
    )
    c4.metric(
        "Alta prioridad", n_high,
        help="Etiqueta high_priority: F2 + F5; útil para revisar residuos del "
             "sitio activo.",
    )
    c5.metric(
        "Conflictos", n_conflict,
        help="Etiqueta conflict: la evidencia contradice la anotación o dos "
             "fuentes discrepan. Se revisa antes de usar la hipótesis.",
    )
    st.caption(
        "Estas cifras describen señales computacionales, no validación. "
        "`supported_hypothesis` no es una función confirmada y "
        "`high_priority` no cierra una brecha."
    )

    with st.expander("ℹ️ Qué mide F1–F5 y qué significa cada confianza"):
        st.markdown(
            "**Los cinco filtros son independientes.** Cada uno aporta una "
            "señal distinta; ninguno, por separado, confirma función."
        )
        st.dataframe(
            pd.DataFrame(_guide_rows()), width="stretch", hide_index=True,
        )
        st.markdown(
            "**Confianza:** combina los filtros que pasaron. Una etiqueta más "
            "fuerte significa más coincidencias, no una prueba."
        )
        st.dataframe(
            pd.DataFrame(_confidence_rows()), width="stretch", hide_index=True,
        )
        st.caption(
            "Conflicto es una señal de revisión (fila o decisión), no una "
            "categoría de apoyo."
        )

    st.markdown("##### Tabla priorizada")
    st.caption(
        "Orden: hipótesis apoyadas, conflictos, alta prioridad, orientativas, "
        "solo estructurales y sin señal; dentro de cada grupo, más fuentes "
        "primero."
    )
    show_cols = [col for col in COLUMNS if col in df.columns]
    st.dataframe(df[show_cols], width="stretch", hide_index=True)

    tsv = df[show_cols].to_csv(sep="\t", index=False)
    st.download_button(
        "⬇️ Hipótesis (TSV)", tsv,
        file_name=f"{genome}_structural_hypotheses.tsv",
        mime="text/plain",
        key=f"hyp_tsv_{genome}",
    )
    st.download_button(
        "⬇️ Evidencia (JSON)",
        json.dumps(raw_rows, indent=2, ensure_ascii=False, default=str),
        file_name=f"{genome}_structural_hypotheses.json",
        mime="application/json",
        key=f"hyp_json_{genome}",
    )

    st.divider()
    st.markdown("##### Ficha de la proteína seleccionada")
    selector_key = f"hyp_{genome}"
    options = list(range(len(ordered)))
    # Si el genoma se volvió a analizar con menos filas, la selección guardada
    # puede quedar fuera de rango; se reinicia para no romper el widget.
    if st.session_state.get(selector_key) not in options:
        st.session_state[selector_key] = options[0]
    selected = st.selectbox(
        "Revisar proteína",
        options=options,
        format_func=lambda index: _option_label(ordered[index]),
        key=selector_key,
    )
    row = ordered[int(selected)]
    pid = _text(row.get("protein_id")) or f"hipótesis {int(selected) + 1}"
    label = confidence_label(row)
    info = confidence_info(label)

    st.markdown(f"### {pid} · `{_text(row.get('hypothesis')) or 'sin hipótesis'}`")
    st.caption(
        f"UniProt: `{_text(row.get('uniprot')) or '—'}` · "
        f"PDB: `{_text(row.get('pdb')) or '—'}` · "
        f"Fuentes: {_as_int(row.get('n_sources'))}"
    )

    m1, m2, m3 = st.columns(3)
    m1.metric(
        "Confianza", f"{info['emoji']} {info['label']}",
        help=info["rule"],
    )
    m2.metric(
        "Fuentes", _as_int(row.get("n_sources")),
        help="Número de filtros F1–F5 que aportaron señal para esta fila.",
    )
    m3.metric(
        "EC", _text(row.get("ec")) or "—",
        help="EC de la plantilla M-CSA; no es una EC medida en esta proteína.",
    )
    st.markdown(f"{info['emoji']} **{info['label']}** · `{label}`")
    st.write(info["meaning"])
    st.caption(f"Regla: {info['rule']}")
    st.caption(f"Límite: {info['caveat']}")

    detail = safe_foldseek_detail(row)
    ev_row = detail["detail"]
    status = []
    for field_name, short in zip(_FLAG_FIELDS, ("F1", "F2", "F3", "F4", "F5")):
        status.append(f"{short} {'sí' if is_on(row.get(field_name)) else 'no'}")
    st.caption(
        "Filtros que pasaron: " + " · ".join(status)
        + " · el estado es señal computacional, no función confirmada."
    )

    with st.expander("F1 · Foldseek", expanded=True):
        st.caption(
            "Mide homología de pliegue contra PDB100/AFDB. No prueba que la "
            "función sea la misma."
        )
        st.write(f"Top: **{detail['top'] or '—'}**")
        if ev_row.get("foldseek_bits") not in (None, ""):
            st.write(f"Bits: `{ev_row['foldseek_bits']}`")
        if detail["detail"]:
            try:
                st.json(detail["detail"])
            except Exception as exc:  # pragma: no cover - datos no serializables
                st.warning(
                    "El detalle no se pudo mostrar como JSON "
                    f"({type(exc).__name__}); se muestra como texto."
                )
                st.text(_json_text(detail["detail"]))
        elif detail["text"]:
            st.text(detail["text"])
            st.caption(detail["message"])
        else:
            st.caption(detail["message"])

    with st.expander("F2 · EnzyMM / M-CSA"):
        st.caption(
            "Mide la plantilla enzimática M-CSA y su EC. La EC pertenece a la "
            "plantilla, no a esta proteína."
        )
        st.write(f"Plantillas: `{_text(row.get('template_mcsa')) or '—'}`")
        st.write(f"EC: `{_text(row.get('ec')) or '—'}`")

    with st.expander("F3 · P2Rank"):
        st.caption(
            "Mide bolsillos y cavidades predichos en la superficie. No dice "
            "qué ligando entra ni si hay catálisis."
        )
        st.write(f"Bolsillos: **{_text(row.get('pocket')) or '—'}**")

    with st.expander("F4 · AlphaFill"):
        st.caption(
            "Mide ligandos/cofactores transferidos desde estructuras "
            "homólogas. Un cofactor transferido no demuestra unión aquí."
        )
        st.write(
            f"Ligandos/cofactores: `{_text(row.get('cofactor')) or '—'}`"
        )

    with st.expander("F5 · Sitio activo"):
        st.caption(
            "Mide posiciones de sitio activo curadas (UniProt), mapeadas al "
            "modelo. Es anotación transferida, no una medición experimental."
        )
        positions = (ev_row.get("active_positions")
                     or row.get("active_positions") or "")
        st.write(f"Posiciones: `{_text(positions) or '—'}`")

    with st.expander("Fuente y límites"):
        st.write(_text(row.get("limitation"))
                 or "Hipótesis computacional; no confirma función.")
        st.caption(credits.line("foldseek", "enzymm", "p2rank", "alphafill",
                                "mcsa", "afdb"))

    dec = _decisions(genome).get(pid, _text(row.get("review_decision")) or "pending")
    new = st.selectbox(
        "Decisión de revisión", DECISIONS,
        index=DECISIONS.index(dec) if dec in DECISIONS else 0,
        key=f"dec_{genome}_{pid}",
        help="Se guarda solo en la sesión; no modifica la anotación base.",
    )
    if new != dec:
        _decisions(genome)[pid] = new
        st.success(f"Decisión guardada en la sesión: {new}")
    dec_rows = [{"genome": genome, "protein_id": k, "decision": v}
                for k, v in _decisions(genome).items()]
    if dec_rows:
        st.download_button(
            "⬇️ Decisiones de revisión (TSV)",
            pd.DataFrame(dec_rows).to_csv(sep="\t", index=False),
            file_name=f"{genome}_review_decisions.tsv",
            mime="text/plain",
            key=f"hyp_decisions_{genome}",
        )

    st.markdown("##### Estructura 3D del candidato")
    st.caption(
        "Una estructura predicha o experimental no confirma la hipótesis; "
        "solo permite inspeccionar la geometría."
    )
    pdb_path = _pdb_path(row.get("pdb"))
    if pdb_path is None:
        st.info(
            "No hay PDB asociado a esta hipótesis; no se inventa una "
            "estructura."
        )
        return
    try:
        structure_viewer.render(
            pdb_path, key=f"hyp_pdb_{genome}_{pid}",
            caption=f"PDB: {pdb_path.name}",
        )
    except Exception as exc:  # pragma: no cover - depende del entorno 3D
        st.warning(
            "No se pudo dibujar el visor 3D "
            f"({type(exc).__name__}: {exc}). La descarga del PDB sigue "
            "disponible."
        )
    try:
        st.download_button(
            "⬇️ PDB del candidato", pdb_path.read_bytes(),
            file_name=pdb_path.name, mime="chemical/x-pdb",
            key=f"hyp_pdb_dl_{genome}_{pid}",
        )
    except OSError as exc:
        st.warning(f"No se pudo leer el PDB para descargarlo: {exc}")
