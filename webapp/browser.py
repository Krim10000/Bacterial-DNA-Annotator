"""browser.py — navegador de genoma en Streamlit.

Dibuja una ventana del genoma con pistas apiladas (genes, operones, promotores,
horquillas) sobre un eje de coordenadas común, y permite buscar por coordenada,
gen, producto o reacción.

Principios de diseño
--------------------
- **Solo lectura**: no escribe datos del proyecto.
- **Sin dependencias nuevas**: el dibujo es SVG generado aquí; Streamlit solo lo
  incrusta. Nada de plotly ni componentes externos.
- **Un eje, varias pistas, una principal**: el usuario elige qué pista manda y
  las demás se atenúan, para no saturar.
- **Dos carriles de genes**: ``+`` arriba y ``−`` abajo; cada gen lleva una
  flecha en su sentido de lectura (5'→3') y la tabla de solapamientos avisa de
  los CDS que comparten coordenadas sin descartar ninguno.
- **Nada de decisiones automáticas**: aquí se muestra lo que otros pasos
  calcularon, con su clase y su confianza.
- **Lo que se ve es lo que se lista**: la tabla «elementos en pantalla» y su
  contador salen de la misma ventana que el dibujo; una fila se puede llevar
  al SVG, resaltarla y compartirla con ``?sel=tipo:id``.
- **Lineal + circular**: la vista circular complementaria recuerda que la
  ventana lineal es una porción del contig, no el contig completo.

Uso desde `app.py`
------------------
    import browser
    browser.render(data_root=DATA_ROOT, genome=selected_genome)
"""
from __future__ import annotations

from pathlib import Path
import csv
import math
import re

import pandas as pd
import streamlit as st
from streamlit.components.v1 import declare_component

from i18n import tr
import gene_card as _gene_card

# ---------------------------------------------------------------- constantes
GENE_FWD = "#4a7ebb"
GENE_REV = "#b07a3a"
HAIRPIN = {"A_fuerte": "#b30000", "B_media": "#c58a00", "C_debil": "#9aa0ad"}
HAIRPIN_MEANING = {
    "A_fuerte": "estructura clara, tallo largo y GC por encima del fondo del genoma",
    "B_media": "intermedia: no se distingue del azar",
    "C_debil": "débil, probablemente sin función de terminación",
}
PROMOTER = "#177245"
TSS_C = "#177245"    # medido (sólido) o alto rendimiento (atenuado)
TTS_C = "#8a5a22"
OPERON = "#8e6cc9"

WINDOW_SIZES = {"5 kb": 5_000, "10 kb": 10_000, "25 kb": 25_000,
                "50 kb": 50_000, "100 kb": 100_000}

#: Tamaños de ventana en pb, ordenados, para elegir el zoom de un elemento.
WINDOW_STEPS = tuple(sorted(WINDOW_SIZES.values()))

#: Etiquetas humanas de las pistas que se pueden resumir y resaltar.
ELEMENT_KIND_LABELS = {
    "gene": "Genes",
    "operon": "Operones",
    "promoter": "Promotores",
    "hairpin": "Horquillas",
    "tss": "TSS",
    "tts": "TTS",
    "complex": "Complejos",
    "suspect": "Sospechas",
}

ELEMENT_KIND_NAMES = {
    "gene": "Gen",
    "operon": "Operón",
    "promoter": "Promotor",
    "hairpin": "Horquilla",
    "tss": "TSS",
    "tts": "TTS",
    "complex": "Complejo",
    "suspect": "Sospecha",
}

ELEMENT_KIND_SINGULAR = {
    "gene": "gen",
    "operon": "operón",
    "promoter": "promotor",
    "hairpin": "horquilla",
    "tss": "TSS",
    "tts": "TTS",
    "complex": "complejo",
    "suspect": "sospecha",
}

#: Columnas de la tabla «elementos en pantalla». ``kind``/``id`` son la clave
#: que comparten la tabla, el SVG (``data-kind``/``data-id``) y ``_render_pick``.
ELEMENT_COLUMNS = ("tipo", "kind", "id", "etiqueta", "inicio", "fin", "detalle")

ELEMENTS_ON_SCREEN_NOTE = (
    "La tabla lista **solo los elementos dibujados en la ventana actual** "
    "(genes, operones, promotores, horquillas, TSS/TTS, complejos y sospechas). "
    "No es el genoma completo: al mover la ventana cambia el conjunto. La "
    "columna *tipo* indica la pista y el contador de arriba resume cuántos hay."
)

# --------------------------------------------------------------- ayuda in situ
#: Título común de los expanders de ayuda. Siempre es el mismo para que el
#: lector reconozca el patrón «aquí me explican lo que estoy viendo».
HELP_EXPANDER_TITLE = "❓ ¿Qué estoy viendo?"

ALL_ELEMENTS_LABEL = "Todos los elementos en pantalla"
ALL_ELEMENTS_LEGEND = (
    "**«Todos los elementos en pantalla»** reúne en una sola tabla todo lo "
    "dibujado en la ventana actual y lo resalta en el dibujo. No añade "
    "elementos de fuera de la ventana ni del resto del genoma."
)

#: Los genes se dibujan en dos carriles: ``+`` arriba y ``−`` abajo. Cada gen
#: lleva una flecha en su sentido de lectura (5'→3'), también en el carril ``−``.
GENE_LANE_ORDER = ("+", "−")
GENE_LANE_Y = {"+": 78.0, "−": 126.0}
GENE_LANE_DIRECTION = {"+": "→", "−": "←"}
GENE_LANE_READING = {"+": "left-to-right", "−": "right-to-left"}
GENE_LANE_LABELS = {
    "+": "GENES + · 5'→3' →",
    "−": "GENES − · 5'→3' ←",
}
GENE_LANE_COLUMNS = (
    "id", "etiqueta", "lane", "hebra", "start", "end", "y", "direccion",
    "color", "locus_tag", "gene_id", "gene_name",
)

#: Título y aviso del expander de solapamientos entre CDS.
OVERLAPS_TITLE = "Solapamientos en pantalla"
OVERLAP_WARNING = (
    "Aviso: **pyrodigal** (y el modo meta en general) puede conservar "
    "**ambos** CDS solapados. Este navegador no elimina ninguno ni decide "
    "cuál es correcto: solo muestra el solapamiento y sus datos. El marco "
    "(`frame`/`phase`) solo se muestra si viene en la anotación; si no, se "
    "marca `unknown` y no se infiere."
)

#: Un par de CDS comparte locus si coincide su identificador o si el solape
#: cubre al menos esta fracción del CDS más corto (llamadas anidadas/duplicadas).
LOCUS_SHARE_FRACTION = 0.8

#: Columnas estables de la tabla de solapamientos (pura, sin Streamlit).
OVERLAP_COLUMNS = (
    "gen_a", "gen_b", "etiqueta_a", "etiqueta_b",
    "inicio_a", "fin_a", "inicio_b", "fin_b",
    "inicio", "fin", "bp_solape",
    "hebra_a", "hebra_b", "hebras", "misma_hebra",
    "marco_a", "marco_b", "mismo_marco", "relacion_marco",
    "locus_compartido", "cobertura_corta",
)

#: Columnas de la vista legible que consume el expander de solapamientos.
OVERLAP_VIEW_COLUMNS = (
    "Par", "Gen A", "Gen B", "Solape", "bp", "Hebras", "Marco",
    "Locus compartido",
)

ELEMENTS_ON_SCREEN_HELP = (
    "**Qué lista.** Cada fila es un elemento **dibujado en la ventana actual**, "
    "no el genoma completo: genes, operones, promotores, horquillas, TSS/TTS, "
    "complejos y sospechas metabólicas. Al mover la ventana o cambiar de "
    "contig, la tabla cambia.\n\n"
    "**Cómo leerla.** *tipo* es la pista; *etiqueta* es el nombre visible; "
    "*inicio* y *fin* son coordenadas en pb; *detalle* resume la evidencia. El "
    "contador y la tabla salen de la misma consulta a la ventana, por eso "
    "coinciden con el dibujo.\n\n"
    "**Qué puedes hacer.** Haz clic en una fila para abrir su ficha y llevarla "
    "al dibujo; «Ver en el dibujo» fija `?sel=tipo:id`. Si algo no aparece, "
    "activa su pista arriba o mueve la ventana.\n\n"
    "**Límite.** Que un elemento esté en pantalla no significa que esté "
    "validado: es lo que el pipeline anotó para esa región, con su clase y su "
    "confianza."
)

HAIRPIN_HELP = (
    "**Qué es.** Una horquilla es un plegamiento del ARN (tallo y bucle) que "
    "se dibuja en las coordenadas del ADN donde está codificado. Un `ΔG` más "
    "negativo indica un tallo más estable; también se muestran la longitud del "
    "tallo, el GC y la cola poli-U.\n\n"
    "**Terminación.** Una horquilla al **final** de un operón, sobre todo si va "
    "seguida de poli-U, es candidata a terminador: corta la transcripción y "
    "cierra la unidad.\n\n"
    "**Atenuación.** Una horquilla **dentro** del operón puede ser un "
    "atenuador: termina el transcrito antes de tiempo según las condiciones, "
    "pero el operón sigue existiendo y puede regularse. Por eso la posición "
    "relativa al operón importa tanto como la estructura.\n\n"
    "**Clases A/B/C.** A = estructura clara, tallo largo y GC sobre el fondo "
    "del genoma; B = intermedia, no se distingue del azar; C = débil, "
    "probablemente sin función de terminación.\n\n"
    "**Límite.** Es una predicción termodinámica y de contexto, no una medida "
    "experimental: una horquilla aislada no prueba terminación ni atenuación."
)

CIRCULAR_HELP = (
    "**Lineal vs circular.** La vista de arriba es una **ventana lineal**: una "
    "porción del contig entre dos coordenadas. La vista circular dibuja el "
    "**contig completo** y marca en naranja la porción que ocupa esa ventana.\n\n"
    "**Cómo leerla.** El arco naranja es la ventana lineal, no un cromosoma "
    "aparte. Los ticks marcan posiciones; el centro resume qué fracción del "
    "contig se está viendo. Si la ventana cruza el origen, el arco da la "
    "vuelta.\n\n"
    "**Para qué sirve.** Para ubicar la ventana en el contexto del contig "
    "entero y notar cuánto queda fuera. La tabla, las fichas y la selección "
    "siguen usando la ventana lineal.\n\n"
    "**Límite.** La longitud del contig se estima desde la anotación (extremo "
    "del último gen), no desde el ensamblado completo; en genomas draft puede "
    "faltar material. Un contig no es necesariamente el cromosoma completo."
)

OPERON_HELP = (
    "**Operón.** Un operón es un grupo de genes contiguos, en la misma hebra, "
    "que se transcriben juntos. Aquí se propone por proximidad y hebra, con "
    "una confianza; no es una prueba experimental.\n\n"
    "**Promotor.** Es la región donde se une la ARN polimerasa (con su factor "
    "σ) para iniciar la transcripción. Se resumen las cajas **−35** y **−10**, "
    "el espaciador, la σ, el score y los desajustes.\n\n"
    "**Cómo leerlos juntos.** El promotor se asocia al operón y se dibuja como "
    "una marca vertical; al elegirlo se abre su ficha y el operón queda "
    "resaltado. Una confianza alta no sustituye a un TSS curado ni a un "
    "terminador experimental.\n\n"
    "**Límite.** En la demo los operones son una agrupación por proximidad y "
    "los genes, predicciones ab initio: sirven para probar el navegador, no "
    "para conclusiones biológicas."
)

HELP_TEXTS = {
    "elements": ELEMENTS_ON_SCREEN_HELP,
    "hairpins": HAIRPIN_HELP,
    "circular": CIRCULAR_HELP,
    "operons": OPERON_HELP,
}


def help_expander(topic: str) -> None:
    """Muestra la explicación in situ de ``topic`` bajo el título común.

    No depende de datos ni de estado: Streamlit solo pinta el expander. Si el
    tema no existe, no dibuja nada (útil para llamadas defensivas).
    """
    text = HELP_TEXTS.get(str(topic))
    if not text:
        return
    with st.expander(HELP_EXPANDER_TITLE):
        st.markdown(text)


def _is_demo_value(value) -> bool:
    """``True`` si el valor de un query param representa el genoma demo."""
    if isinstance(value, (list, tuple)):
        return any(str(v).strip() == "demo" for v in value)
    return str(value or "").strip() == "demo"


def clear_demo_query(query):
    """Devuelve ``query`` sin el parámetro ``genoma=demo``.

    Acepta un mapping (``dict`` o ``st.query_params``) y devuelve un ``dict``
    nuevo, o una cadena tipo ``?genoma=demo&pos=10`` y devuelve una cadena de
    la misma forma. No muta la entrada ni importa Streamlit: es una función
    pura y se puede probar sin runtime.

    El resto de parámetros (``pos``, ``w``, ``contig``, ``sel``, ``modo``…) se
    conserva tal cual. Solo se elimina ``genoma`` cuando su valor es ``demo``;
    un ``genoma=GCF_…`` compartible sigue intacto.
    """
    if query is None:
        return {}
    if isinstance(query, str):
        prefix = "?" if query.startswith("?") else ""
        raw = query[1:] if prefix else query
        pairs = [p for p in raw.split("&") if p]
        if not pairs:
            return ""
        kept = []
        for pair in pairs:
            key, _, value = pair.partition("=")
            if key.strip() == "genoma" and _is_demo_value(value):
                continue
            kept.append(pair)
        return prefix + "&".join(kept)
    if not hasattr(query, "items"):
        return {}
    cleaned: dict = {}
    for key, value in query.items():
        if str(key) == "genoma" and _is_demo_value(value):
            continue
        cleaned[key] = value
    return cleaned


def _clear_demo_query_from_url() -> bool:
    """Quita ``?genoma=demo`` de la URL tras consumirlo en esta ejecución.

    El enlace ``?modo=navegador&genoma=demo`` sigue abriendo la demo, pero al
    navegar a otro genoma el parámetro viejo no vuelve a forzar la selección en
    cada rerun. Devuelve ``True`` si eliminó el parámetro.
    """
    try:
        current = dict(st.query_params)
    except Exception:
        return False
    cleaned = clear_demo_query(current)
    if "genoma" not in current or "genoma" in cleaned:
        return False
    try:
        del st.query_params["genoma"]
        return True
    except (KeyError, TypeError):
        return False


# Componente propio: convierte el SVG en algo clicable sin añadir dependencias.
_COMPONENT_DIR = Path(__file__).resolve().parent / "components" / "genome_track"
_genome_track = declare_component("genome_track", path=str(_COMPONENT_DIR))


def genome_track(svg_markup: str, key: str = "genome_track", *,
                 pos: int | None = None, contig: str | None = None,
                 window: int | str | None = None) -> dict | None:
    """Dibuja el SVG y devuelve el elemento pulsado (o None).

    ``pos``/``contig``/``window`` permiten que el componente construya una URL
    de respaldo si Streamlit no procesa el valor del clic.
    """
    value = _genome_track(svg=svg_markup, show_hint=True,
                          pos=pos, contig=contig, window=window, key=key)
    return value if isinstance(value, dict) else None


# ------------------------------------------------------------------ lectura
def _p(data_root: Path, *parts: str) -> Path:
    return Path(data_root).joinpath(*parts)


DEMO_DIR = Path(__file__).resolve().parent / "demo"


def load_tables(data_root: Path, genome: str) -> dict:
    """Carga las pistas de un genoma. Las que falten se devuelven vacías.

    El genoma ``demo`` se lee de ``webapp/demo/`` si no hay datos locales, para
    que el navegador funcione en un clon recién descargado.
    """
    reg = _p(data_root, "regulation", genome)
    if genome == "demo" and not (reg / "genes.tsv").exists():
        reg = DEMO_DIR
    out: dict[str, pd.DataFrame] = {}
    for name in ("genes", "operons", "promoters"):
        path = reg / f"{name}.tsv"
        out[name] = (pd.read_csv(path, sep="\t", dtype=str).fillna("")
                     if path.exists() else pd.DataFrame())
    hp = _p(data_root, "operon_finder", f"{genome}_hairpins.tsv")
    if not hp.exists() and genome == "demo":
        hp = DEMO_DIR / "hairpins.tsv"
    out["hairpins"] = (pd.read_csv(hp, sep="\t", dtype=str).fillna("")
                       if hp.exists() else pd.DataFrame())
    cx = _p(data_root, "complexes", f"{genome}.complexes.tsv")
    if not cx.exists() and genome == "demo":
        cx = DEMO_DIR / "complexes.tsv"
    out["complexes"] = (pd.read_csv(cx, sep="\t", dtype=str).fillna("")
                        if cx.exists() else pd.DataFrame())
    if out["complexes"].empty:
        # Transferencia K-12: no trae coordenadas; se derivan de genes.tsv.
        tp = _p(data_root, "complexes", f"transfer_{genome}.tsv")
        if tp.exists() and not out["genes"].empty:
            tr = pd.read_csv(tp, sep="\t", dtype=str).fillna("")
            gmap = {str(r.get("gene_id") or r.get("locus_tag")): r
                    for _, r in out["genes"].iterrows()}
            rows = []
            for _, r in tr.iterrows():
                ids = [x for x in str(r.get("genes_genoma", "")).split(",")
                       if x in gmap]
                if len(ids) < 2:
                    continue
                coords = [(int(gmap[x]["start"]), int(gmap[x]["end"])) for x in ids]
                expected = int(r.get("expected_n") or 0)
                mapped = len(ids)
                status = ("complex_inferred_complete"
                          if expected and mapped >= expected
                          else "complex_inferred_incomplete")
                rows.append({
                    "complex_id": r.get("complex_id", ""),
                    "name": "K-12 transferido: " + str(r.get("complex_name", "")),
                    "contig": gmap[ids[0]].get("contig", ""),
                    "start": min(a for a, _ in coords),
                    "end": max(b for _, b in coords),
                    "subunits": ";".join(ids), "expected": r.get("expected_n", ""),
                    "missing": "", "status": status,
                    "evidence": (f"RBH K-12; frac={r.get('frac_mapped','')}; "
                                 f"phylo={r.get('phylo_jacc','')}; "
                                 f"masa={r.get('mass_kda','')} kDa"),
                    "confidence": "high" if float(r.get("frac_mapped") or 0) >= 0.8 else "medium",
                    "source": "transferencia por ortología K-12",
                    "flag": "blue" if status == "complex_inferred_complete" else "yellow",
                    "limitation": r.get("limitation", ""),
                })
            out["complexes"] = pd.DataFrame(rows)
    # Modelo K-12 en lote para los 400: se suma como sospecha amarilla.
    mp = _p(data_root, "complexes", f"model_{genome}.complexes.tsv")
    if mp.exists():
        try:
            model = pd.read_csv(mp, sep="\t", dtype=str).fillna("")
            if not model.empty:
                if out["complexes"].empty:
                    out["complexes"] = model
                else:
                    out["complexes"] = pd.concat([out["complexes"], model],
                                                 ignore_index=True)
                    out["complexes"] = out["complexes"].drop_duplicates(
                        subset=["complex_id"], keep="first")
        except Exception:
            pass
    for key in ("genes", "operons", "complexes"):
        df = out[key]
        for col in ("start", "end"):
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    df = out["hairpins"]
    for col in ("start", "end", "dG", "stem"):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return out


TRANS_FILES = {
    "ecoli_k12": {
        "tss": ("regulondb", "tss_k12_union.tsv"),
        "tts": ("regulondb", "TerminatorSet.tsv"),
    },
    "ecoli": {"tss": ("regulondb", "tss_k12_union.tsv"),
              "tts": ("regulondb", "TerminatorSet.tsv")},
}


def load_transcription(data_root: Path, genome: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """TSS y TTS del genoma, con su procedencia (medido vs inferido).

    La distinción importa: un TSS curado es una medición; los de alto
    rendimiento son mediciones de otro estudio; y una predicción sería una
    inferencia. El navegador lo dice en cada marca.
    """
    root = Path(data_root)
    tss_rows, tts_rows = [], []
    # K-12: unión de curados + alto rendimiento
    tss_path = root / "regulondb" / "tss_k12_union.tsv"
    if genome.startswith("ecoli") and tss_path.exists():
        for r in csv.DictReader(open(tss_path), delimiter="\t"):
            if not (r.get("pos") or "").isdigit():
                continue
            src = "alto rendimiento" if (r.get("sigma") or "").startswith("HT") else "curado"
            tss_rows.append({"pos": int(r["pos"]), "strand": r.get("strand", "+"),
                             "source": src, "detail": r.get("sigma", "")})
    # B. subtilis: BSGatlas
    bsg = root / "dbtbs" / "gt_tss.tsv"
    if genome in ("subtilis", "bacillus") and bsg.exists():
        for r in csv.DictReader(open(bsg), delimiter="\t"):
            if (r.get("pos") or "").isdigit():
                tss_rows.append({"pos": int(r["pos"]), "strand": r.get("strand", "+"),
                                 "source": "BSGatlas", "detail": r.get("sigma", "")})
    # TTS
    term = root / "regulondb" / "TerminatorSet.tsv"
    if genome.startswith("ecoli") and term.exists():
        for r in csv.DictReader(open(term), delimiter="\t"):
            v = (r.get("2)leftEndPos") or "").strip()
            if v.isdigit():
                st = (r.get("4)strand") or "").strip().lower()
                tts_rows.append({"pos": int(v), "strand": "-" if st.startswith("r") else "+",
                                 "source": "curado", "detail": (r.get("7)type") or "")})
    bterm = root / "dbtbs" / "gt_terminators.tsv"
    if genome in ("subtilis", "bacillus") and bterm.exists():
        for r in csv.DictReader(open(bterm), delimiter="\t"):
            if (r.get("start") or "").isdigit():
                tts_rows.append({"pos": int(r["start"]), "strand": r.get("strand", "+"),
                                 "source": "BSGatlas", "detail": ""})
    return pd.DataFrame(tss_rows), pd.DataFrame(tts_rows)


def genome_gc(data_root: Path, genome: str) -> float | None:
    """GC del genoma a partir de los genes (suficiente para un rótulo)."""
    path = _p(data_root, "regulation", genome, "genes.tsv")
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c in ("gc",))
    if "gc" in df:
        v = pd.to_numeric(df["gc"], errors="coerce").dropna()
        if len(v):
            return float(v.iloc[0])
    return None


# ----------------------------------------------------------------- búsqueda
def search(data_root: Path, genome: str, query: str, tables: dict) -> list[dict]:
    """Busca por coordenadas, gen, producto, EC o reacción.

    Devuelve entradas {etiqueta, contig, pos, tipo} listas para saltar.
    """
    q = (query or "").strip()
    if not q:
        return []
    hits: list[dict] = []
    genes = tables.get("genes", pd.DataFrame())
    contigs = genes["contig"].unique().tolist() if "contig" in genes else []

    # 1) coordenadas "1000000-1050000" o "1000000"
    parts = q.replace("–", "-").split("-")
    if parts and all(p.strip().isdigit() for p in parts if p.strip()):
        if len(parts) == 2:
            a, b = int(parts[0]), int(parts[1])
            return [{"etiqueta": f"Región {a:,}–{b:,}", "contig": None,
                     "pos": a, "fin": b, "tipo": "coordenada"}]
        return [{"etiqueta": f"Posición {int(parts[0]):,}", "contig": None,
                 "pos": int(parts[0]), "fin": None, "tipo": "coordenada"}]

    # 2) reacción Rhea
    if q.upper().startswith("RHEA:"):
        gm = _p(data_root, "metabolism", "genemaps", f"{genome}_gene_map.tsv")
        if gm.exists():
            df = pd.read_csv(gm, sep="\t", dtype=str).fillna("")
            for _, r in df[df["rhea_id"].str.upper() == q.upper()].head(20).iterrows():
                hits.append({"etiqueta": f"{r.get('gene_id','')} · {r.get('product','')[:60]}",
                             "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                             "fin": None, "tipo": "reacción"})
        return hits

    # 3) EC "2.7.1.39"
    if q.count(".") == 3 and q.replace(".", "").isdigit() and "ec" in genes:
        for _, r in genes[genes["ec"].astype(str).str.startswith(q)].head(20).iterrows():
            hits.append({"etiqueta": f"{r.get('locus_tag','')} · EC {r.get('ec','')}",
                         "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                         "fin": None, "tipo": "EC"})
        return hits

    # 4) nombre de gen, locus o producto
    if not genes.empty:
        mask = pd.Series(False, index=genes.index)
        for col in ("locus_tag", "gene_name", "gene_id", "product"):
            if col in genes:
                mask |= genes[col].astype(str).str.contains(q, case=False, na=False)
        for _, r in genes[mask].head(25).iterrows():
            etiqueta = r.get("gene_name") or r.get("locus_tag") or r.get("gene_id")
            hits.append({"etiqueta": f"{etiqueta} · {str(r.get('product',''))[:60]}",
                         "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                         "fin": None, "tipo": "gen"})
    return hits


# --------------------------------------------------------------------- svg
def _esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def _esc_attr(s) -> str:
    return _esc(s).replace('"', "&quot;")


def _split_focus(focus: str) -> tuple[str | None, str]:
    """Separa ``"gene:b0001"`` en ``("gene", "b0001")``.

    Si no hay prefijo de tipo se conserva el comportamiento histórico (un gen
    se puede seleccionar solo por su nombre). El identificador de una horquilla
    ya contiene ``:`` (``contig:inicio-fin``), así que solo se parte la primera
    aparición.
    """
    raw = str(focus or "")
    if ":" in raw:
        kind, _, ident = raw.partition(":")
        if kind in ELEMENT_KIND_LABELS:
            return kind, ident
    return None, raw


def svg(tables: dict, contig: str, start: int, end: int,
        tracks: set[str] | None = None, focus: str = "",
        suspects: pd.DataFrame | None = None,
        tss: pd.DataFrame | None = None, tts: pd.DataFrame | None = None,
        *, highlight_all: bool = False) -> str:
    """Dibuja la ventana como SVG ancho.

    Cada pista vive en su propia banda vertical, así que nunca se solapan: se
    pueden mostrar todas a la vez sin que se pisen. Las etiquetas de los genes
    se escalonan en dos filas cuando están muy juntas.

    ``focus`` acepta ``"tipo:id"`` (por ejemplo ``gene:b0001`` u
    ``operon:op1``) y añade la clase ``sel-hit`` al grupo correspondiente.
    ``highlight_all`` marca todos los elementos clicables con ``all-hit``.
    """
    tracks = tracks or {"genes", "operones", "promotores", "horquillas",
                        "transcripcion", "complejos"}
    span = max(1, end - start)
    width, height = 1180, 448
    pad_l, pad_r = 96, 12
    plot_w = width - pad_l - pad_r
    focus_kind, focus_id = _split_focus(focus)

    def x(bp: float) -> float:
        return pad_l + (bp - start) / span * plot_w

    def hit_class(kind: str, ident) -> str:
        cls = "hit"
        if highlight_all:
            cls += " all-hit"
        if focus_id and str(ident) == focus_id and focus_kind in (None, kind):
            cls += " sel-hit"
        return cls

    def hit_open(kind: str, ident) -> str:
        return (f'<g class="{hit_class(kind, ident)}" data-kind="{kind}" '
                f'data-id="{_esc_attr(ident)}">')


    genes = tables.get("genes", pd.DataFrame())
    operons = tables.get("operons", pd.DataFrame())
    promoters = tables.get("promoters", pd.DataFrame())
    hairpins = tables.get("hairpins", pd.DataFrame())
    complexes = tables.get("complexes", pd.DataFrame())
    if "contig" in genes and contig:
        genes = genes[genes["contig"] == contig]
        if "contig" in hairpins:
            hairpins = hairpins[hairpins["contig"] == contig]
    if not genes.empty:
        genes = genes[(genes["end"] >= start) & (genes["start"] <= end)]
    if not hairpins.empty:
        hairpins = hairpins[(hairpins["end"] >= start) & (hairpins["start"] <= end)]
    if not operons.empty and "contig" in operons:
        operons = operons[(operons["contig"] == contig)
                          & (operons["end"] >= start) & (operons["start"] <= end)]
    if not complexes.empty and "contig" in complexes:
        complexes = complexes[(complexes["contig"] == contig)
                              & (complexes["end"] >= start)
                              & (complexes["start"] <= end)]

    Y = {"ruler": 28,
         "genes_plus": GENE_LANE_Y["+"], "genes_minus": GENE_LANE_Y["−"],
         "operones": 188, "promotores": 232, "horquillas": 278,
         "transcripcion": 320, "complejos": 362, "sospecha": 404}
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" height="{height}" '
         f'viewBox="0 0 {width} {height}" font-family="sans-serif" '
         f'style="background:#fff;border-radius:8px">']
    if focus_id or highlight_all:
        o.append("<style>"
                 ".hit.sel-hit>rect,.hit.sel-hit>circle,.hit.sel-hit>line,"
                 ".hit.sel-hit>polygon{stroke:#ff4b4b;stroke-width:3}"
                 ".hit.all-hit{filter:drop-shadow(0 0 2px #f0a400)}"
                 "</style>")

    # regla
    o.append(f'<line x1="{pad_l}" y1="{Y["ruler"]}" x2="{width-pad_r}" y2="{Y["ruler"]}" stroke="#c9cfda"/>')
    for k in range(6):
        bp = start + span * k / 5
        xx = x(bp)
        o.append(f'<line x1="{xx:.1f}" y1="{Y["ruler"]-5}" x2="{xx:.1f}" y2="{Y["ruler"]+5}" stroke="#c9cfda"/>')
        anchor = "start" if k == 0 else ("end" if k == 5 else "middle")
        o.append(f'<text x="{xx:.1f}" y="{Y["ruler"]-10}" font-size="11" fill="#8b90a0" '
                 f'text-anchor="{anchor}">{bp:,.0f}</text>')

    # rótulos de banda (los de genes van dentro de cada carril)
    for key, label in (("operones", "OPERONES"),
                       ("promotores", "PROMOTORES"), ("horquillas", "HORQUILLAS"),
                       ("transcripcion", "TSS / TTS"),
                       ("complejos", "COMPLEJOS")):
        active = key in tracks
        o.append(f'<text x="8" y="{Y[key]-6}" font-size="9.5" fill="#808495" '
                 f'font-weight="700" opacity="{1 if active else 0.25}">{label}</text>')

    # genes en DOS carriles: + arriba, − abajo. La punta de cada gen apunta
    # siempre en el sentido de lectura 5'→3' (derecha en +, izquierda en −).
    if "genes" in tracks:
        layout = gene_lane_layout(genes)
        for lane in GENE_LANE_ORDER:
            yg = GENE_LANE_Y[lane]
            lane_rows = layout[layout["lane"] == lane]
            o.append(f'<g data-role="gene-lane" data-lane="{lane}" '
                     f'data-reading="{GENE_LANE_READING[lane]}" '
                     f'data-direction="{GENE_LANE_DIRECTION[lane]}">')
            # banda de fondo: separa visualmente los dos carriles
            o.append(f'<rect data-role="lane-band" x="{pad_l}" y="{yg-16:.1f}" '
                     f'width="{plot_w}" height="32" fill="#f7f9fc" rx="4"/>')
            o.append(f'<text data-role="lane-label" x="8" y="{yg-14:.1f}" '
                     f'font-size="9.5" fill="#808495" font-weight="700">'
                     f'{_esc(GENE_LANE_LABELS[lane])}</text>')
            labels: list[tuple[float, str]] = []
            for _, g in lane_rows.iterrows():
                gx0, gx1 = x(g["start"]), x(g["end"])
                w = max(3.0, gx1 - gx0)
                color = str(g["color"])
                name = str(g["etiqueta"])
                gid = str(g["id"])
                hebra = str(g["hebra"])
                rev = lane == "−"
                hebra_txt = hebra if hebra in ("+", "−") else "desconocida"
                tip_txt = (f"{name or gid} · {int(g['start']):,}–{int(g['end']):,} pb · "
                           f"hebra {hebra_txt} · 5'→3' {GENE_LANE_DIRECTION[lane]}")
                o.append(hit_open("gene", gid) + f'<title>{_esc(tip_txt)}</title>')
                # zona de clic algo mayor que el rectángulo visible
                o.append(f'<rect x="{gx0-2:.1f}" y="{yg-13:.1f}" '
                         f'width="{max(6.0, w+4):.1f}" height="26" fill="transparent"/>')
                o.append(f'<rect x="{gx0:.1f}" y="{yg-8:.1f}" width="{w:.1f}" '
                         f'height="16" fill="{color}"/>')
                if w > 10:
                    tip = min(8.0, w / 3)
                else:
                    # gen muy estrecho: una punta mínima igualmente marca 5'→3'
                    tip = 2.5
                if rev:
                    o.append(f'<polygon points="{gx0:.1f},{yg-8:.1f} {gx0:.1f},{yg+8:.1f} '
                             f'{gx0-tip:.1f},{yg:.1f}" fill="{color}"/>')
                else:
                    o.append(f'<polygon points="{gx1:.1f},{yg-8:.1f} {gx1:.1f},{yg+8:.1f} '
                             f'{gx1+tip:.1f},{yg:.1f}" fill="{color}"/>')
                gene_focus = (focus_kind in (None, "gene") and focus_id
                              and focus_id in (str(g.get("locus_tag") or ""),
                                               str(g.get("gene_id") or ""),
                                               str(g.get("gene_name") or "")))
                if gene_focus:
                    o.append(f'<rect x="{gx0-3:.1f}" y="{yg-12:.1f}" width="{w+6:.1f}" '
                             f'height="24" fill="none" stroke="#ff4b4b" stroke-width="2"/>')
                o.append("</g>")
                if w > 24 and name:
                    labels.append((gx0 + w / 2, name))
            # dos filas por carril para que las etiquetas no se pisen
            labels.sort()
            last_x, row = -999, 0
            for lx, name in labels:
                if lx - last_x < 52:
                    row = 1 - row
                else:
                    row = 0
                last_x = lx
                o.append(f'<text x="{lx:.1f}" y="{yg+20+row*10:.1f}" font-size="9.5" '
                         f'fill="#3b3f4a" text-anchor="middle">{_esc(name)[:16]}</text>')
            o.append("</g>")

    # operones
    if "operones" in tracks:
        for _, op in operons.iterrows():
            gx0, gx1 = x(op["start"]), x(op["end"])
            lab = op.get("label") or op.get("operon_id") or ""
            oid = op.get("operon_id") or lab or f"{int(op['start'])}-{int(op['end'])}"
            o.append(hit_open("operon", oid)
                     + f'<title>{_esc(lab or oid)} · operón</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["operones"]-8}" width="{max(4, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="none" stroke="{OPERON}" stroke-width="1.8"/>')
            if gx1 - gx0 > 70 and lab:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["operones"]+3}" font-size="9" '
                         f'fill="{OPERON}">{_esc(lab)[:30]}</text>')
            o.append("</g>")

    # promotores
    if "promotores" in tracks and not promoters.empty and not operons.empty:
        by_id = {str(r.get("operon_id")): r for _, r in operons.iterrows()}
        for _, p in promoters.iterrows():
            op = by_id.get(str(p.get("operon_id")))
            if op is None:
                continue
            gx = x(op["start"])
            sigma = p.get("sigma", "")
            pid = p.get("operon_id") or op.get("operon_id") or f"{int(op['start'])}"
            o.append(hit_open("promoter", pid)
                     + f'<title>Promotor {_esc(str(pid))} · σ{_esc(str(sigma))}</title>')
            o.append(f'<line x1="{gx:.1f}" y1="{Y["promotores"]-9}" x2="{gx:.1f}" '
                     f'y2="{Y["promotores"]+9}" stroke="{PROMOTER}" stroke-width="2"/>')
            o.append(f'<circle cx="{gx:.1f}" cy="{Y["promotores"]-11}" r="3" fill="{PROMOTER}"/>')
            o.append("</g>")
            # el sigma no se repite en el dibujo: va en la leyenda y en el detalle

    # horquillas
    if "horquillas" in tracks:
        for hi, (_, h) in enumerate(hairpins.iterrows()):
            gx = x(h["start"])
            color = HAIRPIN.get(str(h.get("class", "C_debil")), HAIRPIN["C_debil"])
            hid = f"{contig}:{int(h['start'])}-{int(h['end'])}"
            o.append(hit_open("hairpin", hid)
                     + f'<title>Horquilla {_esc(str(h.get("class", "")))} · '
                     f'{int(h["start"]):,}–{int(h["end"]):,} pb</title>')
            o.append(f'<line x1="{gx:.1f}" y1="{Y["horquillas"]-10}" x2="{gx:.1f}" '
                     f'y2="{Y["horquillas"]+10}" stroke="{color}" stroke-width="1.6"/>')
            o.append(f'<circle cx="{gx:.1f}" cy="{Y["horquillas"]-12}" r="5" fill="none" '
                     f'stroke="{color}" stroke-width="2"/>')
            o.append("</g>")

    # TSS (triángulo hacia arriba) y TTS (hacia abajo): medidos vs inferidos
    if "transcripcion" in tracks:
        yt = Y["transcripcion"]
        for _, t in (tss if tss is not None else pd.DataFrame()).iterrows():
            if not (start <= t["pos"] <= end):
                continue
            gx = x(t["pos"])
            full = t.get("source") in ("curado", "BSGatlas")
            op = 1.0 if full else 0.55
            o.append(hit_open("tss", int(t["pos"]))
                     + f'<title>TSS · {int(t["pos"]):,} pb · '
                     f'{_esc(str(t.get("source", "")))}</title>')
            o.append(f'<polygon points="{gx:.1f},{yt-16} {gx-4:.1f},{yt-9} '
                     f'{gx+4:.1f},{yt-9}" fill="{TSS_C}" opacity="{op}"/>')
            o.append(f'<line x1="{gx:.1f}" y1="{yt-9}" x2="{gx:.1f}" y2="{yt-3}" '
                     f'stroke="{TSS_C}" stroke-width="1" opacity="{op}"/>')
            o.append("</g>")
        for _, t in (tts if tts is not None else pd.DataFrame()).iterrows():
            if not (start <= t["pos"] <= end):
                continue
            gx = x(t["pos"])
            op = 1.0
            o.append(hit_open("tts", int(t["pos"]))
                     + f'<title>TTS · {int(t["pos"]):,} pb</title>')
            o.append(f'<polygon points="{gx:.1f},{yt+16} {gx-4:.1f},{yt+9} '
                     f'{gx+4:.1f},{yt+9}" fill="{TTS_C}" opacity="{op}"/>')
            o.append(f'<line x1="{gx:.1f}" y1="{yt+3}" x2="{gx:.1f}" y2="{yt+9}" '
                     f'stroke="{TTS_C}" stroke-width="1" opacity="{op}"/>')
            o.append("</g>")

    # complejos por vecindad
    if "complejos" in tracks:
        for _, c in complexes.iterrows():
            gx0, gx1 = x(c["start"]), x(c["end"])
            status = str(c.get("status", ""))
            color = {"complex_complete": "#2e7d32",
                     "complex_inferred_complete": "#1565c0",
                     "complex_incomplete": "#c58a00",
                     "complex_inferred_incomplete": "#c58a00",
                     "complex_ambiguous": "#b3261e"}.get(status, "#8b90a0")
            cid = str(c.get("complex_id") or f"{int(c['start'])}-{int(c['end'])}")
            tip = (f"{c.get('name','complejo')} · {status} · "
                   f"{int(c['start']):,}–{int(c['end']):,} pb")
            o.append(hit_open("complex", cid) + f'<title>{_esc(tip)}</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["complejos"]-8}" width="{max(6, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="{color}" opacity="0.22" '
                     f'stroke="{color}" stroke-width="1.6"/>')
            if gx1 - gx0 > 70:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["complejos"]+3}" font-size="8.5" '
                         f'fill="{color}">{_esc(str(c.get("name", "")))[:32]}</text>')
            o.append("</g>")

    # sospecha por cadena metabólica, como banda propia
    if suspects is not None and not suspects.empty:
        o.append(f'<text x="8" y="{Y["sospecha"]-6}" font-size="9.5" fill="#c58a00" '
                 f'font-weight="700">SOSPECHA</text>')
        for si, (_, s) in enumerate(suspects.iterrows()):
            gx0, gx1 = x(s["start"]), x(s["end"])
            sid = f"{int(s['start'])}-{int(s['end'])}"
            o.append(hit_open("suspect", sid)
                     + f'<title>Sospecha metabólica · {_esc(str(s.get("n_genes")))} genes</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["sospecha"]-8}" width="{max(6, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="#ffe9c7" stroke="#c58a00" stroke-width="1.6"/>')
            if gx1 - gx0 > 60:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["sospecha"]+3}" font-size="8.5" '
                         f'fill="#9a5b00">{_esc(str(s.get("n_genes")))} genes · cadena</text>')
            o.append("</g>")

    o.append("</svg>")
    return "".join(o)


# ------------------------------------------------------------------ página
def suspects_window(data_root: Path, genome: str, contig: str,
                    start: int, end: int) -> pd.DataFrame:
    """Sospechas por cadena metabólica que caen en la ventana."""
    path = Path(data_root) / "metabolism" / "operon_suspects.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df = df[df["genome"] == genome]
    if "contig" in df and contig:
        df = df[df["contig"] == contig]
    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df["end"] >= start) & (df["start"] <= end)]
    return df.sort_values("start").reset_index(drop=True)


# ------------------------------------------- geometría circular y selección
# Estas funciones son puras (no tocan Streamlit) para poder probarlas con
# DataFrames sintéticos. La UI solo las usa para fijar la ventana y el resalte.

def _as_int(value) -> int | None:
    """Convierte coordenadas/ids numéricos sin romperse con NaN o cadenas."""
    try:
        if value is None:
            return None
        if isinstance(value, float) and math.isnan(value):
            return None
        text = str(value).strip()
        if not text or text.lower() in ("nan", "none", "<na>"):
            return None
        return int(float(text))
    except (TypeError, ValueError, OverflowError):
        return None


def _first_match(df: pd.DataFrame | None, ident, cols: tuple[str, ...], *,
                 contig: str | None = None) -> dict | None:
    """Primera fila cuyo valor coincida en alguna de ``cols``.

    Tolera tablas vacías y columnas ausentes; si hay varias filas prefiere la
    del contig actual, para que un locus repetido no salte a otro contig.
    """
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return None
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if col in df.columns:
            mask |= df[col].astype(str) == str(ident)
    sub = df[mask]
    if sub.empty:
        return None
    if contig and "contig" in sub.columns:
        same = sub[sub["contig"].astype(str) == str(contig)]
        if not same.empty:
            sub = same
    return sub.iloc[0].to_dict()


def _coords(row: dict | None) -> tuple[int, int] | None:
    if not row:
        return None
    start, end = _as_int(row.get("start")), _as_int(row.get("end"))
    if start is None or end is None:
        return None
    return start, end


# ------------------------------------------- carriles de genes y solapamientos
# Todo lo de esta sección es puro: recibe DataFrames y devuelve DataFrames o
# cadenas. El SVG y los expanders solo consumen estos resultados.
def _text(value) -> str:
    """Texto limpio; convierte ``None``/``NaN``/``<NA>`` en cadena vacía."""
    if value is None:
        return ""
    try:
        if isinstance(value, float) and math.isnan(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "<na>"} else text


def _text_value(row, key: str) -> str:
    """``row[key]`` como texto, tolerando Series, dicts y columnas ausentes."""
    try:
        value = row.get(key)
    except (AttributeError, TypeError, ValueError):
        return ""
    return _text(value)


def normalize_strand(value) -> str:
    """Normaliza la hebra a ``"+"``, ``"−"`` o ``""`` (dato ausente).

    Acepta las convenciones de la anotación (``1``/``-1``) y variantes de texto.
    No adivina: un valor que no reconoce queda como desconocido.
    """
    text = _text(value).lower().replace("−", "-")
    if text in {"1", "+", "+1", "plus", "forward", "fwd", "directa"}:
        return "+"
    if text in {"-1", "-", "minus", "reverse", "rev", "reversa", "inversa"}:
        return "−"
    return ""


def normalize_frame(value) -> str:
    """Marco de lectura ``"0"``/``"1"``/``"2"`` o ``"unknown"``.

    Solo se acepta un marco explícito (columna ``frame``/``phase``). No se
    deduce de las coordenadas: si falta el dato, la respuesta es ``unknown``.
    """
    number = _as_int(value)
    if number in (0, 1, 2):
        return str(number)
    return "unknown"


def gene_frame(row) -> str:
    """Marco declarado en la fila; ``unknown`` si ninguna columna lo trae."""
    for key in ("frame", "phase", "reading_frame"):
        text = _text_value(row, key)
        if not text:
            continue
        frame = normalize_frame(text)
        if frame != "unknown":
            return frame
    return "unknown"


def _gene_ident(row) -> str:
    ident = (_text_value(row, "locus_tag") or _text_value(row, "gene_id")
             or _text_value(row, "gene_name"))
    if ident:
        return ident
    start, end = _as_int(row.get("start")), _as_int(row.get("end"))
    return f"{start}-{end}"


def _gene_label(row, ident: str) -> str:
    return (_text_value(row, "gene_name") or _text_value(row, "locus_tag")
            or _text_value(row, "gene_id") or ident)


def _window_genes(genes: pd.DataFrame, contig: str | None,
                  start: int | None, end: int | None) -> pd.DataFrame:
    """Filas de genes del contig que se solapan con la ventana."""
    sub = genes
    if contig is not None and "contig" in sub.columns:
        sub = sub[sub["contig"].astype(str) == str(contig)]
    if start is not None and end is not None:
        s_num = pd.to_numeric(sub["start"], errors="coerce")
        e_num = pd.to_numeric(sub["end"], errors="coerce")
        sub = sub[(e_num >= int(start)) & (s_num <= int(end))]
    return sub


def gene_lane_layout(genes, contig: str | None = None,
                     start: int | None = None,
                     end: int | None = None) -> pd.DataFrame:
    """Layout puro de los genes en dos carriles: ``+`` arriba, ``−`` abajo.

    Devuelve una fila por gen con ``lane``, ``y`` (menor es más arriba),
    ``direccion`` (``→``/``←``, sentido 5'→3') y los datos de identidad que
    comparten el SVG y la tabla. Un gen sin hebra reconocible se coloca en el
    carril ``+`` con ``hebra="?"`` y color neutro, sin inventar su orientación.
    """
    empty = pd.DataFrame({col: pd.Series(dtype="object")
                          for col in GENE_LANE_COLUMNS})
    if not isinstance(genes, pd.DataFrame) or genes.empty:
        return empty
    if "start" not in genes.columns or "end" not in genes.columns:
        return empty
    sub = _window_genes(genes, contig, start, end)
    rows: list[dict] = []
    for _, g in sub.iterrows():
        s_int, e_int = _as_int(g.get("start")), _as_int(g.get("end"))
        if s_int is None or e_int is None:
            continue
        hebra = normalize_strand(g.get("strand"))
        lane = hebra if hebra in GENE_LANE_ORDER else "+"
        ident = _gene_ident(g)
        color = (GENE_REV if lane == "−"
                 else (GENE_FWD if hebra == "+" else "#8b90a0"))
        rows.append({
            "id": ident,
            "etiqueta": _gene_label(g, ident),
            "lane": lane,
            "hebra": hebra or "?",
            "start": s_int,
            "end": e_int,
            "y": GENE_LANE_Y[lane],
            "direccion": GENE_LANE_DIRECTION[lane],
            "color": color,
            "locus_tag": _text_value(g, "locus_tag"),
            "gene_id": _text_value(g, "gene_id"),
            "gene_name": _text_value(g, "gene_name"),
        })
    if not rows:
        return empty
    table = pd.DataFrame(rows, columns=list(GENE_LANE_COLUMNS))
    order = {lane: i for i, lane in enumerate(GENE_LANE_ORDER)}
    table["_orden"] = table["lane"].map(order).fillna(len(order))
    table = table.sort_values(["_orden", "start", "end", "id"])
    return table.drop(columns="_orden").reset_index(drop=True)


def _locus_keys(row) -> set[str]:
    keys: set[str] = set()
    for col in ("locus_tag", "gene_id"):
        value = _text_value(row, col)
        if value:
            keys.add(f"{col}:{value}")
    return keys


def cds_overlaps(genes, contig: str | None = None,
                 start: int | None = None, end: int | None = None, *,
                 min_bp: int = 1,
                 locus_share_fraction: float = LOCUS_SHARE_FRACTION) -> pd.DataFrame:
    """Pares de CDS que comparten coordenadas en la misma ventana/contig.

    Reglas, sin inventar datos:

    - solo compara genes del mismo ``contig`` cuyas coordenadas se intersectan
      con la ventana pedida; el solape se mide sobre los CDS completos;
    - ``coordenadas`` son inclusivas (1-based), como en la anotación de
      pyrodigal: dos CDS que comparten un único pb sí se solapan;
    - ``mismo_marco`` solo es ``True``/``False`` si **ambos** traen marco;
      en caso contrario vale ``None`` y ``relacion_marco`` es ``unknown``;
    - ``locus_compartido`` es ``True`` si comparten ``locus_tag``/``gene_id`` o
      si el solape cubre ``locus_share_fraction`` del CDS más corto (llamadas
      anidadas o duplicadas sobre el mismo locus).
    """
    empty = pd.DataFrame({col: pd.Series(dtype="object")
                          for col in OVERLAP_COLUMNS})
    if not isinstance(genes, pd.DataFrame) or genes.empty:
        return empty
    if "start" not in genes.columns or "end" not in genes.columns:
        return empty
    sub = _window_genes(genes, contig, start, end)
    # Si la anotación distingue tipos, el solapamiento es entre CDS: tRNA/rRNA/
    # ncRNA no cuentan. Un tipo vacío se conserva (no se inventa que no sea CDS).
    if "gene_type" in sub.columns:
        types = sub["gene_type"].astype(str).str.strip().str.lower()
        non_cds = {"trna", "rrna", "ncrna", "tmrna", "rna", "pseudogene"}
        sub = sub[~types.isin(non_cds)]
    records: list[dict] = []
    for _, g in sub.iterrows():
        s_int, e_int = _as_int(g.get("start")), _as_int(g.get("end"))
        if s_int is None or e_int is None or e_int < s_int:
            continue
        ident = _gene_ident(g)
        records.append({
            "id": ident,
            "label": _gene_label(g, ident),
            "start": s_int,
            "end": e_int,
            "strand": normalize_strand(g.get("strand")),
            "frame": gene_frame(g),
            "locus": _locus_keys(g),
        })
    if len(records) < 2:
        return empty
    records.sort(key=lambda r: (r["start"], r["end"], r["id"]))
    rows: list[dict] = []
    for i, a in enumerate(records):
        for b in records[i + 1:]:
            # ordenados por inicio: en cuanto b empieza después del final de a,
            # ningún gen posterior puede solaparse con a.
            if b["start"] > a["end"]:
                break
            ov_start, ov_end = max(a["start"], b["start"]), min(a["end"], b["end"])
            if ov_end < ov_start:
                continue
            bp = ov_end - ov_start + 1
            if bp < int(min_bp):
                continue
            len_a = a["end"] - a["start"] + 1
            len_b = b["end"] - b["start"] + 1
            short = max(1, min(len_a, len_b))
            coverage = bp / short
            same_locus = bool(a["locus"] & b["locus"])
            locus_shared = bool(same_locus or coverage >= float(locus_share_fraction))

            hebra_a, hebra_b = a["strand"] or "?", b["strand"] or "?"
            if hebra_a != "?" and hebra_b != "?":
                same_strand = hebra_a == hebra_b
                hebras = (f"misma ({hebra_a}/{hebra_b})" if same_strand
                          else f"opuestas ({hebra_a}/{hebra_b})")
            else:
                same_strand = None
                hebras = "desconocida"

            frame_a, frame_b = a["frame"], b["frame"]
            if frame_a == "unknown" or frame_b == "unknown":
                same_frame, frame_relation = None, "unknown"
            else:
                same_frame = frame_a == frame_b
                frame_relation = "mismo" if same_frame else "distinto"

            rows.append({
                "gen_a": a["id"], "gen_b": b["id"],
                "etiqueta_a": a["label"], "etiqueta_b": b["label"],
                "inicio_a": a["start"], "fin_a": a["end"],
                "inicio_b": b["start"], "fin_b": b["end"],
                "inicio": ov_start, "fin": ov_end, "bp_solape": bp,
                "hebra_a": hebra_a, "hebra_b": hebra_b,
                "hebras": hebras, "misma_hebra": same_strand,
                "marco_a": frame_a, "marco_b": frame_b,
                "mismo_marco": same_frame, "relacion_marco": frame_relation,
                "locus_compartido": locus_shared,
                "cobertura_corta": round(coverage, 4),
            })
    if not rows:
        return empty
    table = pd.DataFrame(rows, columns=list(OVERLAP_COLUMNS))
    table = table.sort_values(
        ["inicio", "bp_solape", "gen_a", "gen_b"],
        ascending=[True, False, True, True],
    )
    return table.reset_index(drop=True)


def overlap_summary(overlaps: pd.DataFrame | None) -> str:
    """Frase-resumen de los solapamientos, con los contadores que importan."""
    if (overlaps is None or not isinstance(overlaps, pd.DataFrame)
            or overlaps.empty or "bp_solape" not in overlaps.columns):
        return "0 solapamientos entre CDS en pantalla."
    total = len(overlaps)

    def _true(value) -> bool:
        return value is not None and bool(value)

    def _count(column: str, predicate) -> int:
        if column not in overlaps.columns:
            return 0
        return sum(1 for value in overlaps[column] if predicate(value))

    same_strand = _count("misma_hebra", lambda v: v is not None and bool(v))
    opposite = _count("misma_hebra",
                      lambda v: v is not None and not bool(v))
    unknown_strand = _count("misma_hebra", lambda v: v is None)
    same_frame = _count("mismo_marco", lambda v: v is not None and bool(v))
    diff_frame = _count("mismo_marco",
                        lambda v: v is not None and not bool(v))
    unknown_frame = _count("mismo_marco", lambda v: v is None)
    shared_locus = _count("locus_compartido", _true)
    noun = "solapamiento" if total == 1 else "solapamientos"
    parts = [
        f"{total} {noun} entre CDS en pantalla",
        f"misma hebra: {same_strand}",
        f"hebras opuestas: {opposite}",
        f"mismo marco: {same_frame}",
        f"distinto marco: {diff_frame}",
        f"marco unknown: {unknown_frame}",
        f"locus compartido: {shared_locus}",
    ]
    if unknown_strand:
        parts.append(f"hebra desconocida: {unknown_strand}")
    return " · ".join(parts) + "."


def overlap_table_view(overlaps: pd.DataFrame | None) -> pd.DataFrame:
    """Vista legible (columnas en español) de la tabla de solapamientos.

    Es pura para poder probarla sin Streamlit. Muestra ``unknown`` cuando el
    marco no está en los datos y no rellena ese hueco con una suposición.
    """
    empty = pd.DataFrame({col: pd.Series(dtype="object")
                          for col in OVERLAP_VIEW_COLUMNS})
    if (overlaps is None or not isinstance(overlaps, pd.DataFrame)
            or overlaps.empty):
        return empty
    rows: list[dict] = []
    for _, r in overlaps.iterrows():
        frame_a = str(r.get("marco_a", "unknown"))
        frame_b = str(r.get("marco_b", "unknown"))
        relation = str(r.get("relacion_marco", "unknown"))
        if relation == "unknown":
            frame_txt = f"{frame_a} / {frame_b} · unknown"
        else:
            frame_txt = f"{frame_a} / {frame_b} · {relation} marco"
        rows.append({
            "Par": f"{r.get('gen_a', '')} ↔ {r.get('gen_b', '')}",
            "Gen A": (f"{r.get('etiqueta_a', '')} "
                      f"({int(r.get('inicio_a', 0)):,}–{int(r.get('fin_a', 0)):,})"),
            "Gen B": (f"{r.get('etiqueta_b', '')} "
                      f"({int(r.get('inicio_b', 0)):,}–{int(r.get('fin_b', 0)):,})"),
            "Solape": f"{int(r.get('inicio', 0)):,}–{int(r.get('fin', 0)):,}",
            "bp": int(r.get("bp_solape", 0)),
            "Hebras": str(r.get("hebras", "desconocida")),
            "Marco": frame_txt,
            "Locus compartido": "sí" if bool(r.get("locus_compartido")) else "no",
        })
    return pd.DataFrame(rows, columns=list(OVERLAP_VIEW_COLUMNS))


def element_span(kind: str, ident: str, tables: dict,
                 *, suspects: pd.DataFrame | None = None,
                 contig: str | None = None) -> dict | None:
    """Localiza un elemento de la tabla/SVG y devuelve su intervalo.

    Mantiene la misma identidad que el dibujo: genes por ``locus_tag`` →
    ``gene_id`` → ``gene_name``, horquillas por ``contig:inicio-fin``, etc.
    Los promotores heredan las coordenadas de su operón porque la tabla de
    promotores no trae coordenadas propias.
    """
    kind = str(kind or "").strip().lower()
    ident = str(ident or "").strip()
    tables = tables or {}
    genes = tables.get("genes", pd.DataFrame())
    operons = tables.get("operons", pd.DataFrame())

    if kind == "gene":
        row = _first_match(genes, ident,
                           ("locus_tag", "gene_id", "gene_name"), contig=contig)
        coords = _coords(row)
        if row is None or coords is None:
            return None
        start, end = coords
        gid = (row.get("locus_tag") or row.get("gene_id")
               or row.get("gene_name") or f"{start}-{end}")
        label = (row.get("gene_name") or row.get("locus_tag")
                 or row.get("gene_id") or gid)
        bits = [str(row.get("product") or "").strip()]
        if row.get("ec"):
            bits.append(f"EC {row['ec']}")
        return {"kind": kind, "id": str(gid), "label": str(label),
                "contig": str(row.get("contig") or contig or ""),
                "start": start, "end": end,
                "detail": " · ".join(b for b in bits if b)}

    if kind == "operon":
        row = _first_match(operons, ident, ("operon_id", "label"), contig=contig)
        coords = _coords(row)
        if row is None or coords is None:
            return None
        start, end = coords
        oid = row.get("operon_id") or row.get("label") or f"{start}-{end}"
        label = row.get("label") or row.get("operon_id") or oid
        detail = (f"{row.get('n_genes', '—')} genes · "
                  f"{row.get('span_bp', '—')} pb · {row.get('confidence', '—')}")
        return {"kind": kind, "id": str(oid), "label": str(label),
                "contig": str(row.get("contig") or contig or ""),
                "start": start, "end": end, "detail": detail}

    if kind == "promoter":
        promoters = tables.get("promoters", pd.DataFrame())
        p = _first_match(promoters, ident, ("operon_id", "label"))
        if p is None:
            return None
        oid = str(p.get("operon_id") or p.get("label") or "")
        op = (_first_match(operons, oid, ("operon_id",), contig=contig)
              if oid else None)
        coords = _coords(p)
        if coords is None:
            coords = _coords(op)
        if coords is None:
            return None
        start, end = coords
        pid = (p.get("operon_id") or (op.get("operon_id") if op is not None else None)
               or p.get("label") or f"{start}")
        p_contig = (p.get("contig") or (op.get("contig") if op is not None else None)
                    or contig or "")
        detail = (f"σ{p.get('sigma', '—')} · score {p.get('score', '—')} · "
                  f"−35 {p.get('box35', '—')} · −10 {p.get('box10', '—')}")
        return {"kind": kind, "id": str(pid),
                "label": f"Promotor · {oid or pid}",
                "contig": str(p_contig), "start": start, "end": end,
                "detail": detail}

    if kind == "hairpin":
        hairpins = tables.get("hairpins", pd.DataFrame())
        if (hairpins is None or hairpins.empty
                or "start" not in hairpins.columns or "end" not in hairpins.columns):
            return None
        hp_contig, coords_text = contig, ident
        if ":" in ident:
            hp_contig, _, coords_text = ident.partition(":")
        m = re.match(r"^(-?\d+)-(-?\d+)$", coords_text.strip())
        if not m:
            return None
        start, end = int(m.group(1)), int(m.group(2))
        sub = hairpins
        if "contig" in sub.columns and hp_contig:
            sub = sub[sub["contig"].astype(str) == str(hp_contig)]
        s_num = pd.to_numeric(sub["start"], errors="coerce")
        e_num = pd.to_numeric(sub["end"], errors="coerce")
        hit = sub[(s_num == start) & (e_num == end)]
        if hit.empty:
            return None
        row = hit.iloc[0].to_dict()
        hid = f"{hp_contig or contig}:{start}-{end}"
        detail = (f"clase {row.get('class', '—')} · ΔG {row.get('dG', '—')} · "
                  f"tallo {row.get('stem', '—')} pb")
        return {"kind": kind, "id": hid,
                "label": f"Horquilla {row.get('class', '')}".strip(),
                "contig": str(hp_contig or row.get("contig") or ""),
                "start": start, "end": end, "detail": detail}

    if kind == "complex":
        complexes = tables.get("complexes", pd.DataFrame())
        row = _first_match(complexes, ident, ("complex_id",), contig=contig)
        coords = _coords(row)
        if row is None or coords is None:
            return None
        start, end = coords
        cid = row.get("complex_id") or f"{start}-{end}"
        detail = (f"{row.get('status', '—')} · confianza "
                  f"{row.get('confidence', '—')}")
        return {"kind": kind, "id": str(cid),
                "label": str(row.get("name") or cid),
                "contig": str(row.get("contig") or contig or ""),
                "start": start, "end": end, "detail": detail}

    if kind == "suspect":
        row = _first_match(suspects, ident, ("id", "suspect_id"), contig=contig)
        if row is None and isinstance(suspects, pd.DataFrame) and not suspects.empty:
            m = re.match(r"^(-?\d+)-(-?\d+)$", ident)
            if m and "start" in suspects.columns and "end" in suspects.columns:
                start, end = int(m.group(1)), int(m.group(2))
                s_num = pd.to_numeric(suspects["start"], errors="coerce")
                e_num = pd.to_numeric(suspects["end"], errors="coerce")
                hit = suspects[(s_num == start) & (e_num == end)]
                row = hit.iloc[0].to_dict() if not hit.empty else None
        coords = _coords(row)
        if row is None or coords is None:
            return None
        start, end = coords
        sid = f"{start}-{end}"
        return {"kind": kind, "id": sid,
                "label": f"Sospecha · {row.get('n_genes', '—')} genes",
                "contig": str(row.get("contig") or contig or ""),
                "start": start, "end": end,
                "detail": str(row.get("products") or "")}

    if kind in ("tss", "tts"):
        pos = _as_int(ident)
        if pos is None:
            return None
        return {"kind": kind, "id": str(pos),
                "label": f"{kind.upper()} · {pos:,} pb",
                "contig": str(contig or ""), "start": pos, "end": pos,
                "detail": ""}

    return None


def contig_lengths(tables: dict) -> dict[str, int]:
    """Longitud aproximada de cada contig, a partir de la anotación.

    No se lee el FASTA (el navegador solo tiene tablas); se usa el mayor ``end``
    anotado. Es suficiente para la vista circular y se dice en la leyenda.
    """
    out: dict[str, int] = {}
    for key in ("genes", "operons", "complexes", "hairpins"):
        df = (tables or {}).get(key)
        if (df is None or not isinstance(df, pd.DataFrame) or df.empty
                or "contig" not in df.columns or "end" not in df.columns):
            continue
        sub = df[["contig", "end"]].copy()
        sub["end"] = pd.to_numeric(sub["end"], errors="coerce")
        sub = sub.dropna(subset=["contig", "end"])
        if sub.empty:
            continue
        grouped = sub.groupby(sub["contig"].astype(str))["end"].max()
        for name, value in grouped.items():
            out[str(name)] = max(out.get(str(name), 0), int(value))
    return out


def viewport_for_span(span_start: int, span_end: int, contig_len: int, *,
                      pad: int = 1_000, min_size: int = 1_000,
                      size_options: tuple[int, ...] = WINDOW_STEPS) -> dict:
    """Ventana de zoom centrada en un elemento y contenida en el contig.

    Elige el menor tamaño de ``size_options`` que deja ``pad`` pb a cada lado;
    si el elemento es enorme, usa el mayor tamaño disponible (el resalte sigue
    funcionando aunque no quepa entero).
    """
    length = max(1, int(contig_len))
    start = max(0, min(int(span_start), length))
    end = max(start, min(int(span_end), length))
    content = end - start
    needed = max(int(min_size), content + 2 * int(pad))
    options = sorted({int(o) for o in size_options if int(o) > 0}) or [needed]
    size = next((o for o in options if o >= needed), options[-1])
    center = (start + end) / 2.0
    lo = int(round(center - size / 2.0))
    lo = max(0, min(lo, max(0, length - size)))
    hi = lo + size
    label = next((lab for lab, bp in WINDOW_SIZES.items() if bp == size), None)
    if label is None:
        label = f"{size / 1000:g} kb" if size >= 1_000 else f"{size} pb"
    return {"start": lo, "end": hi, "size": size, "size_label": label,
            "contains": lo <= start and hi >= end}


def selection_window(kind: str, ident: str, tables: dict, contig: str,
                     start: int, end: int, *,
                     contig_len: int | None = None,
                     suspects: pd.DataFrame | None = None,
                     pad: int = 1_000, min_size: int = 1_000) -> dict | None:
    """Selección completa: elemento + ventana de zoom + si ya estaba visible.

    Es la función que usa el botón «Ver en el dibujo». Devuelve ``None`` si el
    elemento no está en las tablas cargadas (por ejemplo, una horquilla que ya
    no existe).
    """
    info = element_span(kind, ident, tables, suspects=suspects, contig=contig)
    if info is None:
        return None
    item_contig = str(info.get("contig") or contig or "")
    length = contig_len if contig_len else contig_lengths(tables).get(item_contig)
    if not length:
        length = max(int(info["end"]), int(end), int(start) + 1, 1)
    length = max(1, int(length), int(info["end"]))
    view = viewport_for_span(info["start"], info["end"], length,
                             pad=pad, min_size=min_size)
    in_current = (item_contig == str(contig)
                  and int(start) <= info["start"] and info["end"] <= int(end))
    return {**info, "contig_len": length,
            "view_start": view["start"], "view_end": view["end"],
            "size": view["size"], "size_label": view["size_label"],
            "in_current_window": in_current,
            "view_contains": view["contains"]}


# ------------------------------------- tabla y contador de elementos en pantalla
def _empty_elements() -> pd.DataFrame:
    return pd.DataFrame({col: pd.Series(dtype="object")
                         for col in ELEMENT_COLUMNS})


def _window_slice(df: pd.DataFrame | None, contig: str,
                  start: int, end: int) -> pd.DataFrame:
    """Filas cuyas coordenadas se solapan con la ventana; tolera tablas vacías."""
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    if "start" not in df.columns or "end" not in df.columns:
        return df.iloc[0:0].copy()
    sub = df
    if contig and "contig" in sub.columns:
        sub = sub[sub["contig"].astype(str) == str(contig)]
    s_num = pd.to_numeric(sub["start"], errors="coerce")
    e_num = pd.to_numeric(sub["end"], errors="coerce")
    return sub[(e_num >= int(start)) & (s_num <= int(end))]


def all_elements_on_screen(tables: dict, contig: str, start: int, end: int,
                           *, tss: pd.DataFrame | None = None,
                           tts: pd.DataFrame | None = None,
                           suspects: pd.DataFrame | None = None) -> pd.DataFrame:
    """Tabla única de todo lo que se dibuja en la ventana actual.

    Devuelve columnas estables incluso sin datos: ``tipo`` (etiqueta humana),
    ``kind``/``id`` (identidad que comparten SVG y ficha), ``etiqueta``,
    ``inicio``, ``fin`` y ``detalle``. No incluye nada fuera de la ventana.
    """
    tables = tables or {}
    rows: list[dict] = []

    def add(kind: str, ident, label, item_start, item_end, detail: str = "") -> None:
        s_int, e_int = _as_int(item_start), _as_int(item_end)
        if s_int is None or e_int is None:
            return
        rows.append({
            "tipo": ELEMENT_KIND_NAMES.get(kind, kind),
            "kind": kind,
            "id": str(ident),
            "etiqueta": str(label or ident),
            "inicio": s_int,
            "fin": e_int,
            "detalle": str(detail or ""),
        })

    genes = _window_slice(tables.get("genes"), contig, start, end)
    for _, r in genes.iterrows():
        s_int, e_int = _as_int(r.get("start")), _as_int(r.get("end"))
        gid = (r.get("locus_tag") or r.get("gene_id") or r.get("gene_name")
               or f"{s_int}-{e_int}")
        label = r.get("gene_name") or r.get("locus_tag") or r.get("gene_id") or gid
        detail = " · ".join(str(v) for v in
                            (r.get("product"), r.get("ec")) if v)
        add("gene", gid, label, s_int, e_int, detail)

    operons = _window_slice(tables.get("operons"), contig, start, end)
    for _, r in operons.iterrows():
        s_int, e_int = _as_int(r.get("start")), _as_int(r.get("end"))
        oid = r.get("operon_id") or r.get("label") or f"{s_int}-{e_int}"
        detail = (f"{r.get('n_genes', '—')} genes · "
                  f"{r.get('span_bp', '—')} pb · {r.get('confidence', '—')}")
        add("operon", oid, r.get("label") or oid, s_int, e_int, detail)

    promoters = tables.get("promoters")
    if isinstance(promoters, pd.DataFrame) and not promoters.empty:
        all_operons = tables.get("operons")
        op_by_id: dict[str, dict] = {}
        if isinstance(all_operons, pd.DataFrame) and not all_operons.empty \
                and "operon_id" in all_operons.columns:
            op_by_id = {str(op.get("operon_id")): op
                        for _, op in all_operons.iterrows()}
        for _, p in promoters.iterrows():
            oid = str(p.get("operon_id") or p.get("label") or "")
            op = op_by_id.get(oid)
            s_int, e_int = _as_int(p.get("start")), _as_int(p.get("end"))
            if s_int is None or e_int is None:
                s_int = _as_int(op.get("start")) if op is not None else None
                e_int = _as_int(op.get("end")) if op is not None else None
            if s_int is None or e_int is None:
                continue
            p_contig = (p.get("contig")
                        or (op.get("contig") if op is not None else None)
                        or contig or "")
            if contig and p_contig and str(p_contig) != str(contig):
                continue
            if e_int < int(start) or s_int > int(end):
                continue
            pid = (p.get("operon_id")
                   or (op.get("operon_id") if op is not None else None)
                   or p.get("label") or f"{s_int}")
            detail = (f"σ{p.get('sigma', '—')} · score {p.get('score', '—')} · "
                      f"−35 {p.get('box35', '—')} · −10 {p.get('box10', '—')}")
            add("promoter", pid, f"Promotor · {oid or pid}", s_int, e_int, detail)

    hairpins = _window_slice(tables.get("hairpins"), contig, start, end)
    for _, r in hairpins.iterrows():
        s_int, e_int = _as_int(r.get("start")), _as_int(r.get("end"))
        if s_int is None or e_int is None:
            continue
        hp_contig = str(r.get("contig") or contig or "")
        hid = f"{contig}:{s_int}-{e_int}"
        detail = (f"clase {r.get('class', '—')} · ΔG {r.get('dG', '—')} · "
                  f"tallo {r.get('stem', '—')} pb")
        add("hairpin", hid, f"Horquilla {r.get('class', '')}".strip(),
            s_int, e_int, detail)

    for kind, frame in (("tss", tss), ("tts", tts)):
        if (frame is None or not isinstance(frame, pd.DataFrame) or frame.empty
                or "pos" not in frame.columns):
            continue
        positions = pd.to_numeric(frame["pos"], errors="coerce")
        sub = frame[(positions >= int(start)) & (positions <= int(end))]
        for _, r in sub.iterrows():
            pos = _as_int(r.get("pos"))
            if pos is None:
                continue
            detail = " · ".join(str(v) for v in
                                (r.get("source"), r.get("detail")) if v)
            add(kind, str(pos), f"{kind.upper()} · {pos:,} pb",
                pos, pos, detail)

    complexes = _window_slice(tables.get("complexes"), contig, start, end)
    for _, r in complexes.iterrows():
        s_int, e_int = _as_int(r.get("start")), _as_int(r.get("end"))
        cid = r.get("complex_id") or f"{s_int}-{e_int}"
        detail = (f"{r.get('status', '—')} · confianza "
                  f"{r.get('confidence', '—')}")
        add("complex", cid, r.get("name") or cid, s_int, e_int, detail)

    suspect_rows = _window_slice(suspects, contig, start, end)
    for _, r in suspect_rows.iterrows():
        s_int, e_int = _as_int(r.get("start")), _as_int(r.get("end"))
        sid = f"{s_int}-{e_int}"
        add("suspect", sid, f"Sospecha · {r.get('n_genes', '—')} genes",
            s_int, e_int, str(r.get("products") or ""))

    if not rows:
        return _empty_elements()
    table = pd.DataFrame(rows, columns=list(ELEMENT_COLUMNS))
    return table.sort_values(["inicio", "fin", "tipo"]).reset_index(drop=True)


def element_counts(table: pd.DataFrame | None) -> dict[str, int]:
    """Conteo por ``kind`` de la tabla de elementos en pantalla."""
    if (table is None or not isinstance(table, pd.DataFrame) or table.empty
            or "kind" not in table.columns):
        return {}
    counts = table["kind"].astype(str).value_counts()
    return {str(k): int(v) for k, v in counts.items()}


def elements_counter_text(counts: dict[str, int] | None) -> str:
    """Frase del contador: total + desglose por pista, en español."""
    clean = {str(k): int(v) for k, v in (counts or {}).items() if int(v) > 0}
    total = sum(clean.values())
    if total == 0:
        return "0 elementos en pantalla"
    noun = "elemento" if total == 1 else "elementos"
    ordered = [k for k in ELEMENT_KIND_LABELS if k in clean]
    ordered += [k for k in clean if k not in ELEMENT_KIND_LABELS]
    parts = []
    for k in ordered:
        n = clean[k]
        label = (ELEMENT_KIND_SINGULAR.get(k, k) if n == 1
                 else ELEMENT_KIND_LABELS.get(k, k).lower())
        parts.append(f"{n} {label}")
    return f"{total} {noun} en pantalla: " + " · ".join(parts)


def parse_sel(value) -> dict | None:
    """Convierte ``?sel=gene:b0001`` en ``{"kind": "gene", "id": "b0001"}``."""
    raw = str(value or "").strip()
    if not raw or ":" not in raw:
        return None
    kind, _, ident = raw.partition(":")
    kind = kind.strip().lower()
    ident = ident.strip()
    if kind not in ELEMENT_KIND_LABELS or not ident:
        return None
    return {"kind": kind, "id": ident}


def _safe_slug(value: str) -> str:
    """Nombre de archivo seguro para el CSV de exportación."""
    text = "".join(ch if (ch.isalnum() or ch in "-_.") else "_"
                   for ch in str(value))
    return text.strip("._")[:60] or "genoma"



# ------------------------------------------------------- vista circular SVG
def _polar(cx: float, cy: float, radius: float,
           angle_deg: float) -> tuple[float, float]:
    """Punto en un círculo; 0° arriba y los ángulos crecen en sentido horario."""
    rad = math.radians(float(angle_deg))
    return (cx + radius * math.sin(rad), cy - radius * math.cos(rad))


def _arc_path(cx: float, cy: float, radius: float,
              angle0: float, angle1: float) -> str:
    """Camino SVG de un arco circular (sentido horario)."""
    sweep = float(angle1) - float(angle0)
    if abs(sweep) < 1e-9:
        return ""
    p0 = _polar(cx, cy, radius, angle0)
    p1 = _polar(cx, cy, radius, angle1)
    if abs(sweep) >= 359.999:
        mid = float(angle0) + sweep / 2.0
        pm = _polar(cx, cy, radius, mid)
        return (f"M {p0[0]:.2f} {p0[1]:.2f} "
                f"A {radius:.2f} {radius:.2f} 0 1 1 {pm[0]:.2f} {pm[1]:.2f} "
                f"A {radius:.2f} {radius:.2f} 0 1 1 {p1[0]:.2f} {p1[1]:.2f}")
    large = 1 if abs(sweep) > 180 else 0
    sweep_flag = 1 if sweep > 0 else 0
    return (f"M {p0[0]:.2f} {p0[1]:.2f} "
            f"A {radius:.2f} {radius:.2f} 0 {large} {sweep_flag} "
            f"{p1[0]:.2f} {p1[1]:.2f}")


def circular_geometry(contig_len: int, start: int, end: int, *,
                      cx: float = 210.0, cy: float = 220.0,
                      radius: float = 128.0, tick_count: int = 12) -> dict:
    """Geometría pura del contig circular y de la ventana lineal.

    - 0° está arriba y los ángulos crecen en sentido horario.
    - La ventana puede cruzar el origen (``start + span > contig_len``).
    - Si la ventana abarca todo el contig, ``full_circle`` es ``True`` y no hay
      un único arco SVG: se dibuja un anillo completo.
    """
    length = max(1, int(contig_len))
    s = max(0, int(start))
    e = max(s, int(end))
    raw_span = e - s
    full_circle = raw_span >= length
    span = length if full_circle else min(raw_span, length)
    start_mod = s % length
    start_angle = 360.0 * start_mod / length
    sweep = 360.0 if full_circle else 360.0 * span / length
    end_angle = start_angle + sweep
    window_path = ("" if full_circle or sweep <= 0
                   else _arc_path(cx, cy, radius, start_angle, end_angle))
    count = max(4, int(tick_count))
    major_ticks = {0, count // 4, count // 2, (3 * count) // 4}
    ticks: list[dict] = []
    for i in range(count):
        angle = 360.0 * i / count
        pos = length * i / count
        major = i in major_ticks
        x1, y1 = _polar(cx, cy, radius, angle)
        x2, y2 = _polar(cx, cy, radius + (10.0 if major else 5.0), angle)
        lx, ly = _polar(cx, cy, radius + (24.0 if major else 16.0), angle)
        ticks.append({"pos": int(round(pos)), "angle": angle, "major": major,
                      "x1": x1, "y1": y1, "x2": x2, "y2": y2,
                      "label_x": lx, "label_y": ly})
    start_point = _polar(cx, cy, radius, start_angle)
    end_point = _polar(cx, cy, radius, end_angle)
    return {
        "contig_len": length,
        "start": s, "end": e, "span": span,
        "fraction": 1.0 if full_circle else span / length,
        "cx": cx, "cy": cy, "radius": radius,
        "start_angle": start_angle, "end_angle": end_angle,
        "sweep": sweep, "full_circle": full_circle,
        "wraps": (not full_circle) and (start_mod + span) > length,
        "window_path": window_path,
        "start_point": start_point, "end_point": end_point,
        "ticks": ticks,
    }


def circular_svg(contig: str, contig_len: int, start: int, end: int, *,
                 width: int = 560, height: int = 500,
                 tick_count: int = 12) -> str:
    """SVG autónomo del contig circular con la ventana lineal resaltada.

    No usa librerías: es un dibujo geométrico con un anillo gris (contig), un
    arco naranja (ventana), ticks y una leyenda que explica que la ventana
    lineal es una porción del círculo.
    """
    width = max(360, int(width))
    height = max(380, int(height))
    cx = round(width * 0.38, 1)
    cy = round(height * 0.46, 1)
    radius = round(min(width * 0.24, height * 0.28), 1)
    geo = circular_geometry(contig_len, start, end, cx=cx, cy=cy,
                            radius=radius, tick_count=tick_count)
    length = geo["contig_len"]
    frac = geo["fraction"] * 100.0
    contig_safe = _esc_attr(contig)
    aria = (f"Contig {contig} de {length:,} pb en vista circular; la ventana "
            f"lineal {geo['start']:,}–{geo['end']:,} pb es la porción resaltada.")
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" height="auto" '
        f'viewBox="0 0 {width} {height}" font-family="sans-serif" role="img" '
        f'aria-label="{_esc_attr(aria)}" '
        f'style="max-width:{width}px;height:auto;display:block">',
        f'<title>{_esc_attr(aria)}</title>',
        f'<text x="12" y="22" font-size="14" font-weight="700" fill="#2f3440">'
        f'Contig {_esc(contig)}</text>',
        f'<text x="12" y="40" font-size="11" fill="#7a7f8c">'
        f'{length:,} pb · vista circular complementaria</text>',
        f'<circle data-role="contig" cx="{cx}" cy="{cy}" r="{radius}" '
        f'fill="#fbfcfe" stroke="#dfe3ea" stroke-width="18"/>',
    ]
    if geo["full_circle"]:
        out.append(f'<circle data-role="window" cx="{cx}" cy="{cy}" '
                   f'r="{radius}" fill="none" stroke="#e8710a" '
                   f'stroke-width="18" opacity="0.9"/>')
    elif geo["window_path"]:
        out.append(f'<path data-role="window" d="{geo["window_path"]}" '
                   f'fill="none" stroke="#e8710a" stroke-width="18" '
                   f'opacity="0.9" stroke-linecap="butt"/>')
        for role, point in (("start", geo["start_point"]),
                            ("end", geo["end_point"])):
            out.append(f'<circle data-role="window-{role}" cx="{point[0]:.2f}" '
                       f'cy="{point[1]:.2f}" r="4.5" fill="#b34700" '
                       f'stroke="#fff" stroke-width="1.5"/>')
    out.append('<g data-role="ticks">')
    for tick in geo["ticks"]:
        stroke = "#8f96a6" if tick["major"] else "#c3c9d4"
        out.append(f'<line x1="{tick["x1"]:.2f}" y1="{tick["y1"]:.2f}" '
                   f'x2="{tick["x2"]:.2f}" y2="{tick["y2"]:.2f}" '
                   f'stroke="{stroke}" '
                   f'stroke-width="{1.6 if tick["major"] else 0.9}"/>')
        if tick["major"]:
            out.append(f'<text x="{tick["label_x"]:.2f}" '
                       f'y="{tick["label_y"]:.2f}" text-anchor="middle" '
                       f'font-size="9.5" fill="#6f7482">'
                       f'{tick["pos"]:,}</text>')
    out.append("</g>")
    # centro: identidad del contig y fracción que ocupa la ventana
    out.append(f'<text x="{cx}" y="{cy - 16}" text-anchor="middle" '
               f'font-size="12" font-weight="700" fill="#2f3440">'
               f'{_esc(contig)}</text>')
    out.append(f'<text x="{cx}" y="{cy + 2}" text-anchor="middle" '
               f'font-size="11" fill="#5b6070">{length:,} pb</text>')
    out.append(f'<text x="{cx}" y="{cy + 20}" text-anchor="middle" '
               f'font-size="11" fill="#b34700">{frac:.1f}% del contig</text>')
    out.append(f'<text x="{cx}" y="{cy + 37}" text-anchor="middle" '
               f'font-size="10" fill="#7a7f8c">'
               f'ventana {geo["start"]:,}–{geo["end"]:,} pb</text>')
    # leyenda
    legend_y = height - 92
    out.append(f'<g data-role="legend" font-size="11" fill="#3b3f4a">')
    out.append(f'<line x1="28" y1="{legend_y}" x2="58" y2="{legend_y}" '
               f'stroke="#dfe3ea" stroke-width="8"/>')
    out.append(f'<text x="68" y="{legend_y + 4}">Contig completo · '
               f'{length:,} pb</text>')
    out.append(f'<line x1="28" y1="{legend_y + 24}" x2="58" '
               f'y2="{legend_y + 24}" stroke="#e8710a" stroke-width="8"/>')
    out.append(f'<text x="68" y="{legend_y + 28}">Ventana lineal actual · '
               f'{geo["start"]:,}–{geo["end"]:,} pb</text>')
    if geo["wraps"]:
        out.append(f'<text x="68" y="{legend_y + 48}" fill="#b34700">'
                   f'La ventana cruza el origen del contig circular.</text>')
    out.append("</g>")
    out.append(f'<text x="12" y="{height - 14}" font-size="11" fill="#5b6070">'
               f'La ventana lineal es una porción del círculo: el rectángulo de '
               f'arriba muestra solo ese tramo.</text>')
    out.append("</svg>")
    return "".join(out)


# ------------------------------------------------------- ficha al hacer clic
def _match_any(df: pd.DataFrame, value, cols: tuple[str, ...]) -> pd.Series:
    """Máscara booleana tolerante a columnas ausentes."""
    if df.empty:
        return pd.Series(dtype=bool)
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if col in df:
            mask |= df[col].astype(str) == str(value)
    return mask


def _row(df: pd.DataFrame, mask) -> dict | None:
    sub = df[mask] if not df.empty else df
    return sub.iloc[0].to_dict() if len(sub) else None


def load_genome_sequences(data_root: Path, genome: str) -> dict[str, str]:
    """Carga el FASTA local del genoma (raw o raw_excluded) al abrir una ficha."""
    import gzip
    for sub in ("raw", "raw_excluded"):
        path = Path(data_root) / sub / f"{genome}_genomic.fna.gz"
        if not path.exists():
            continue
        out: dict[str, str] = {}
        name, chunks = "", []
        try:
            with gzip.open(path, "rt", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if name:
                            out[name] = "".join(chunks)
                        name = line[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(line)
                if name:
                    out[name] = "".join(chunks)
        except OSError:
            return {}
        return out
    return {}


def _gene_detail(tables: dict, gid: str, data_root: Path | None = None,
                 genome: str = "") -> None:
    g = _row(tables["genes"], _match_any(
        tables["genes"], gid, ("locus_tag", "gene_id", "gene_name")))
    if not g:
        st.warning(tr("Ese gen ya no está en la ventana actual."))
        return
    data_root = Path(data_root) if data_root is not None else Path("data")
    seqs = load_genome_sequences(data_root, genome) if genome else {}
    _gene_card.render_gene_card(
        g, genome, data_root, sequences=seqs,
        key_prefix=f"browser_{genome}_{gid}",
    )


def _operon_detail(tables: dict, oid: str) -> None:
    o = _row(tables["operons"], _match_any(
        tables["operons"], oid, ("operon_id", "label")))
    if not o:
        st.warning(tr("Ese operón ya no está en la ventana actual."))
        return
    st.markdown(f"### 🧩 {o.get('label') or o.get('operon_id')}")
    c1, c2, c3 = st.columns(3)
    c1.metric(tr("Genes"), int(float(o.get("n_genes") or 0)))
    c2.metric(tr("Extensión"), f"{int(float(o.get('span_bp') or 0)):,} pb")
    c3.metric(tr("Hebra"), "−" if str(o.get("strand")) in ("-1", "-") else "+")
    st.markdown(f"**{tr('Genes')}:** {o.get('genes') or '—'}")
    st.markdown(f"**{tr('Nombres')}:** {o.get('gene_names') or '—'}")
    st.caption(tr("Hueco mínimo entre genes: {g} pb · evidencia: {e} · confianza: {c}").format(
        g=o.get("min_gap", "—"), e=o.get("evidence", "—"), c=o.get("confidence", "—")))


def _promoter_detail(tables: dict, oid: str) -> None:
    p = _row(tables["promoters"], _match_any(
        tables["promoters"], oid, ("operon_id", "label")))
    if not p:
        st.warning(tr("Ese promotor ya no está en la ventana actual."))
        return
    st.markdown(f"### 🎯 {tr('Promotor')} · {p.get('label') or p.get('operon_id')}")
    st.markdown(f"**σ:** {p.get('sigma') or '—'} · **score:** {p.get('score') or '—'}")
    st.markdown(f"**−35:** `{p.get('box35') or '—'}` ({tr('desajustes')} {p.get('mm35', '—')})")
    st.markdown(f"**−10:** `{p.get('box10') or '—'}` ({tr('desajustes')} {p.get('mm10', '—')})")
    st.caption(tr("Separación: {s} pb · distancia al inicio: {d} pb · sitios: {n}").format(
        s=p.get("spacer", "—"), d=p.get("dist_10_to_start", "—"), n=p.get("n_sites", "—")))


def _hairpin_detail(tables: dict, hid: str) -> None:
    hp = tables["hairpins"]
    mask = hp.apply(lambda r: f"{r.get('contig')}:{int(r['start'])}-{int(r['end'])}" == hid, axis=1)         if not hp.empty else pd.Series(dtype=bool)
    h = _row(hp, mask)
    if not h:
        st.warning(tr("Esa horquilla ya no está en la ventana actual."))
        return
    st.markdown(f"### 🌀 {tr('Horquilla')} {h.get('class', '')}")
    c1, c2, c3 = st.columns(3)
    c1.metric("ΔG", f"{h.get('dG', '—')} kcal/mol")
    c2.metric(tr("Tallo"), f"{h.get('stem', '—')} pb")
    c3.metric("GC", f"{h.get('gc', '—')}")
    st.markdown(f"**{tr('Bucle')}:** {h.get('loop', '—')} · **poli-U:** {h.get('polyU', '—')}")
    st.markdown(f"**{tr('Genes vecinos')}:** {h.get('left_gene') or '—'} ← → {h.get('right_gene') or '—'}")
    st.caption(HAIRPIN_MEANING.get(str(h.get("class")), ""))


def _tx_detail(tables: dict, kind: str, pos: str) -> None:
    label = "TSS" if kind == "tss" else "TTS"
    st.markdown(f"### 📡 {label} · {int(float(pos)):,} pb")
    st.caption(tr("Marca de {k}: {w}. Se muestra como triángulo en la pista de "
                  "transcripción.").format(
        k=label,
        w=tr("inicio de la transcripción (donde empieza el ARN)")
        if kind == "tss" else tr("fin de la transcripción (donde termina el ARN)")))


def _suspect_detail(data_root: Path, genome: str, sid: str) -> None:
    path = Path(data_root) / "metabolism" / "operon_suspects.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df = df[df["genome"] == genome]
    mask = df.apply(lambda r: f"{int(float(r['start']))}-{int(float(r['end']))}" == sid, axis=1)         if not df.empty else pd.Series(dtype=bool)
    r = _row(df, mask)
    if not r:
        st.warning(tr("Esa sospecha ya no está en la ventana actual."))
        return
    st.markdown(f"### 🔗 {tr('Sospecha metabólica')} · {r.get('n_genes')} genes")
    st.markdown(f"**{tr('Genes')}:** {r.get('genes') or '—'}")
    st.caption(tr("Genes contiguos que comparten metabolitos. Es una sospecha de "
                  "operón metabólico, no una prueba."))


def _complex_detail(tables: dict, ident: str) -> None:
    df = tables.get("complexes", pd.DataFrame())
    row = _row(df, _match_any(df, ident, ("complex_id",))) if not df.empty else None
    if not row:
        st.info(tr("Complejo no encontrado."))
        return
    st.markdown(f"### 🧩 {row.get('name', 'Complejo')}")
    st.caption(f"{row.get('status','')} · confianza {row.get('confidence','')}")
    st.write(f"**Subunidades:** {row.get('subunits','')}")
    st.write(f"**Esperadas:** {row.get('expected','')}")
    st.write(f"**Faltan:** {row.get('missing','') or '—'}")
    st.write(f"**Evidencia:** {row.get('evidence','')}")
    st.caption(row.get("limitation", ""))


def _query_value(name: str) -> str | None:
    """Lee un query param sin romperse fuera de una ejecución de Streamlit."""
    try:
        value = st.query_params.get(name)
    except Exception:
        return None
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    return None if value is None else str(value)


def _focus_payload(kind: str, ident: str, *, source: str = "tabla",
                   ts=None) -> dict:
    return {"kind": str(kind), "id": str(ident), "source": source, "ts": ts}


def _focus_key(focus: dict | None) -> str:
    if not isinstance(focus, dict):
        return ""
    kind, ident = focus.get("kind"), focus.get("id")
    if kind and ident is not None:
        return f"{kind}:{ident}"
    return ""


def _activate_element(data_root: Path, genome: str, tables: dict,
                      kind: str, ident: str, *,
                      suspects: pd.DataFrame | None = None) -> None:
    """Fija ``sel`` y, si hace falta, la ventana; luego resalta el elemento.

    Si el elemento ya está dentro de la ventana no se mueve el viewport: se
    resalta donde el usuario lo está viendo. Si está fuera, se deja preparada
    una vista de zoom que se aplica al inicio del siguiente rerun (antes de
    crear los widgets, para que Streamlit no conserve los valores viejos).
    """
    contig = str(st.session_state.get("br_contig") or "")
    try:
        cur_start = int(st.session_state.get("br_start") or 0)
    except (TypeError, ValueError):
        cur_start = 0
    cur_size = WINDOW_SIZES.get(str(st.session_state.get("br_size_sel")), 10_000)
    info = selection_window(kind, ident, tables, contig, cur_start,
                            cur_start + cur_size,
                            contig_len=contig_lengths(tables).get(contig),
                            suspects=suspects)
    if info is None:
        st.warning(tr("Ese elemento ya no está en las tablas cargadas."))
        return
    focus_id = str(info["id"])
    st.session_state["br_focus"] = _focus_payload(kind, focus_id, source="tabla")
    st.session_state["br_focus_key"] = f"{kind}:{focus_id}"
    st.session_state.pop("br_dismissed_key", None)
    st.session_state.pop("br_dismissed", None)
    if info["in_current_window"]:
        view_start = cur_start
        view_label = str(st.session_state.get("br_size_sel") or "")
    else:
        view_start = info["view_start"]
        view_label = info["size_label"]
        st.session_state["br_pending_view"] = {
            "contig": info["contig"], "start": view_start,
            "size_label": view_label,
        }
    try:
        st.query_params["sel"] = f"{kind}:{focus_id}"
        st.query_params["contig"] = str(info["contig"])
        st.query_params["pos"] = str(view_start)
        if view_label in WINDOW_SIZES:
            st.query_params["w"] = view_label
    except Exception:
        pass
    st.rerun()


def _clear_focus() -> None:
    """Quita el resalte y el ``sel`` de la URL, sin tocar la ventana."""
    focused = st.session_state.get("br_focus_key")
    if focused:
        st.session_state["br_auto_disabled"] = focused
    st.session_state.pop("br_focus", None)
    st.session_state.pop("br_focus_key", None)
    st.session_state.pop("br_dismissed_key", None)
    st.session_state.pop("br_url_sel_seen", None)
    st.session_state["br_dismissed"] = None
    try:
        if st.query_params.get("sel"):
            del st.query_params["sel"]
    except (KeyError, TypeError):
        pass
    st.rerun()


def _render_pick(data_root: Path, genome: str, tables: dict, pick: dict) -> None:
    kind, ident = pick.get("kind"), str(pick.get("id", ""))
    with st.container(border=True):
        head, close = st.columns([6, 1])
        with head:
            st.caption(tr("Elemento seleccionado (tabla, dibujo o enlace). "
                          "Esta ficha explica ese elemento."))
        with close:
            if st.button("✖ " + tr("Cerrar"), key="br_close_pick"):
                st.session_state["br_dismissed"] = pick.get("ts")
                st.session_state["br_dismissed_key"] = f"{kind}:{ident}"
                st.session_state.pop("br_focus", None)
                st.session_state.pop("br_focus_key", None)
                st.session_state.pop("br_url_sel_seen", None)
                if st.query_params.get("sel"):
                    try:
                        del st.query_params["sel"]
                    except (KeyError, TypeError):
                        pass
                st.rerun()
        if kind == "gene":
            _gene_detail(tables, ident, data_root, genome)
        elif kind == "operon":
            _operon_detail(tables, ident)
        elif kind == "promoter":
            _promoter_detail(tables, ident)
        elif kind == "hairpin":
            _hairpin_detail(tables, ident)
        elif kind in ("tss", "tts"):
            _tx_detail(tables, kind, ident)
        elif kind == "complex":
            _complex_detail(tables, ident)
        elif kind == "suspect":
            _suspect_detail(data_root, genome, ident)
        else:
            st.info(tr("Elemento sin ficha todavía."))


def render(data_root: Path, genome: str) -> None:
    """Dibuja la vista del navegador para un genoma."""
    # El enlace ?genoma=demo ya se consumió al elegir ``genome``; se retira de
    # la URL para que no vuelva a imponer la demo al cambiar de genoma.
    _clear_demo_query_from_url()
    st.subheader("🧬 Navegador de genoma")
    st.caption(
        "Ventana del genoma con sus pistas apiladas. El **eje es común**: cada "
        "posición en pb se alinea en todas las pistas, y cada pista tiene su "
        "banda para que no se pisen. Marca las que quieras ver."
    )
    st.caption("🧾 " + ALL_ELEMENTS_LEGEND)

    tables = load_tables(data_root, genome)
    if genome == "demo":
        st.info(tr("Demo incluida: 150 kb de *Chlorobaculum tepidum* anotados "
                   "desde cero. Los genes son predicciones ab initio y los "
                   "operones, una agrupación por proximidad; sirve para probar "
                   "el navegador, no para conclusiones biológicas."))
    tss_df, tts_df = load_transcription(data_root, genome)
    genes = tables["genes"]
    if genes.empty:
        st.warning(f"No hay anotación para **{genome}**. Genera primero la capa "
                   "regulatoria con `run_pipeline.sh components`.")
        return
    contigs = genes["contig"].dropna().unique().tolist()
    if not contigs:
        st.warning("El genoma no tiene contigs anotados.")
        return

    max_pos = max(0, int(genes["end"].max()))
    if st.session_state.get("br_genome") != genome:
        ini = 50_000 if genome == "ecoli_k12" else int(genes["start"].min())
        st.session_state.update({
            "br_genome": genome, "br_contig": contigs[0],
            "br_contig_sel": contigs[0], "br_start": ini, "br_pos": ini,
            "br_size": 10_000, "br_size_sel": "10 kb",
        })
        st.session_state.pop("br_focus", None)
        st.session_state.pop("br_focus_key", None)
        st.session_state.pop("br_auto_row", None)
        st.session_state.pop("br_auto_disabled", None)
        st.session_state["br_qp_token"] = None

    # Vista pendiente de un botón «Ver en el dibujo». Se aplica antes de crear
    # los widgets para que los controles no conserven el valor anterior.
    pending = st.session_state.pop("br_pending_view", None)
    if isinstance(pending, dict):
        p_contig = pending.get("contig")
        if p_contig in contigs:
            st.session_state["br_contig"] = p_contig
            st.session_state["br_contig_sel"] = p_contig
        p_start = pending.get("start")
        if p_start is not None:
            try:
                p_start = min(max(0, int(p_start)), max_pos)
            except (TypeError, ValueError):
                p_start = None
            if p_start is not None:
                st.session_state["br_start"] = p_start
                st.session_state["br_pos"] = p_start
        if pending.get("size_label") in WINDOW_SIZES:
            st.session_state["br_size_sel"] = pending["size_label"]

    # Enlace directo a una región: ?pos=123456&w=10+kb&contig=c1&sel=gene:b0048
    # El token evita reaplicar la URL en cada rerun y pisar lo que movió el
    # usuario; solo cambia cuando la URL cambia.
    qp_pos = _query_value("pos")
    qp_w = _query_value("w")
    qp_contig = _query_value("contig")
    qp_token = "|".join(x or "" for x in (qp_pos, qp_w, qp_contig))
    if qp_token != st.session_state.get("br_qp_token"):
        try:
            if qp_pos is not None:
                p_start = min(max(0, int(qp_pos)), max_pos)
                st.session_state["br_start"] = p_start
                st.session_state["br_pos"] = p_start
            if qp_w in WINDOW_SIZES:
                st.session_state["br_size_sel"] = qp_w
            if qp_contig in contigs:
                st.session_state["br_contig"] = qp_contig
                st.session_state["br_contig_sel"] = qp_contig
        except (TypeError, ValueError):
            pass
        st.session_state["br_qp_token"] = qp_token

    # ?sel=tipo:id deja el elemento resaltado al abrir un enlace compartido.
    # Se aplica solo cuando el valor de la URL cambia, para que un clic en el
    # dibujo no sea pisado en el siguiente rerun por el sel anterior.
    raw_sel = _query_value("sel")
    url_sel = parse_sel(raw_sel)
    if url_sel and str(raw_sel or "") != st.session_state.get("br_url_sel_seen"):
        st.session_state["br_url_sel_seen"] = str(raw_sel or "")
        info = element_span(url_sel["kind"], url_sel["id"], tables,
                            contig=st.session_state.get("br_contig"))
        f_kind = url_sel["kind"]
        f_id = str(info["id"]) if info else url_sel["id"]
        f_key = f"{f_kind}:{f_id}"
        if f_key != st.session_state.get("br_focus_key"):
            st.session_state["br_focus"] = _focus_payload(
                f_kind, f_id, source="url")
            st.session_state["br_focus_key"] = f_key
            st.session_state.pop("br_dismissed_key", None)
            if info and f_id != url_sel["id"]:
                try:
                    st.query_params["sel"] = f_key
                except Exception:
                    pass

    # --- búsqueda
    q = st.text_input(
        "🔎 Buscar por coordenada, gen, producto, EC o reacción (RHEA:…)",
        key="br_query", placeholder="p. ej. thrA, proteasa, 1000000-1050000, RHEA:12345",
    )
    if q:
        hits = search(data_root, genome, q, tables)
        if not hits:
            st.caption("Sin resultados en este genoma.")
        else:
            pick = st.selectbox(
                f"{len(hits)} resultado(s)", range(len(hits)),
                format_func=lambda i: f"{hits[i]['etiqueta']}  ·  {hits[i]['contig'] or ''} "
                                      f"{hits[i]['pos']:,}",
                key="br_hit",
            )
            if st.button("Saltar a este resultado", key="br_jump"):
                h = hits[pick]
                st.session_state["br_pending_view"] = {
                    "contig": h["contig"] or st.session_state.get("br_contig"),
                    "start": max(0, h["pos"] - 2_000),
                    "size_label": st.session_state.get("br_size_sel"),
                }
                st.rerun()

    # --- controles
    # Los valores viven en Session State (para que los botones puedan fijar la
    # vista); por eso el widget no recibe ``value``/``index`` a la vez.
    if st.session_state.get("br_contig") not in contigs:
        st.session_state["br_contig"] = contigs[0]
    if st.session_state.get("br_contig_sel") not in contigs:
        st.session_state["br_contig_sel"] = st.session_state["br_contig"]
    if st.session_state.get("br_size_sel") not in WINDOW_SIZES:
        st.session_state["br_size_sel"] = "10 kb"
    st.session_state.setdefault(
        "br_pos", min(int(st.session_state.get("br_start", 0)), max_pos))
    c1, c2, c3, c4, c5 = st.columns([2, 2, 1.6, 1, 1])
    contig = c1.selectbox("Contig", contigs, key="br_contig_sel")
    st.session_state["br_contig"] = contig
    size_label = c2.selectbox("Ventana", list(WINDOW_SIZES), key="br_size_sel")
    size = WINDOW_SIZES[size_label]
    start = c3.number_input("Posición (pb)", 0, max_pos, step=size, key="br_pos")
    st.session_state["br_start"] = int(start)
    if c4.button("◀", key="br_left"):
        st.session_state["br_pending_view"] = {
            "contig": contig, "start": max(0, int(start) - size),
            "size_label": size_label,
        }
        st.rerun()
    if c5.button("▶", key="br_right"):
        st.session_state["br_pending_view"] = {
            "contig": contig, "start": max(0, int(start) + size),
            "size_label": size_label,
        }
        st.rerun()

    # --- pistas visibles: casillas, para poder superponer sin que se pisen
    t1, t2, t3, t4, t5, t6, t7 = st.columns(7)
    tracks: set[str] = set()
    if t1.checkbox("Genes", True, key="br_t_genes"):
        tracks.add("genes")
    if t2.checkbox("Operones", True, key="br_t_op"):
        tracks.add("operones")
    if t3.checkbox("Promotores", True, key="br_t_prom"):
        tracks.add("promotores")
    if t4.checkbox("Horquillas", True, key="br_t_hp"):
        tracks.add("horquillas")
    ver_sospecha = t5.checkbox("Sospecha metabólica", True, key="br_t_susp")
    if t6.checkbox("TSS / TTS", True, key="br_t_tx"):
        tracks.add("transcripcion")
    if t7.checkbox("Complejos", True, key="br_t_complex"):
        tracks.add("complejos")

    start, end = int(start), int(start) + size
    wg = genes[(genes["contig"] == contig)
               & (genes["end"] >= start) & (genes["start"] <= end)].sort_values("start")
    overlaps = cds_overlaps(genes, contig, start, end)
    susp = (suspects_window(data_root, genome, contig, start, end)
            if ver_sospecha else pd.DataFrame())
    screen_elements = all_elements_on_screen(
        tables, contig, start, end, tss=tss_df, tts=tts_df, suspects=susp,
    )
    counts = element_counts(screen_elements)
    if "br_all_elements" not in st.session_state:
        st.session_state["br_all_elements"] = True
    highlight_all = bool(st.session_state.get("br_all_elements", True))

    st.caption(f"{genome} · {contig} · {start:,}–{end:,} pb · "
               f"{len(screen_elements)} elementos en pantalla "
               f"({len(wg)} genes en la ventana)")

    focus = st.session_state.get("br_focus")
    focus = focus if isinstance(focus, dict) else None
    focus_key = _focus_key(focus)
    focus_kind = str(focus.get("kind") or "") if focus else ""
    focus_id = str(focus.get("id") or "") if focus else ""

    clicked = genome_track(
        svg(tables, contig, start, end, tracks=tracks, suspects=susp,
            tss=tss_df, tts=tts_df, focus=focus_key,
            highlight_all=highlight_all),
        key="genome_track", pos=start, contig=contig, window=size_label,
    )
    st.caption(
        "🧬 **Carriles de genes:** **+** arriba (5'→3' →) y **−** abajo "
        "(5'→3' ←). Cada gen es una flecha en su sentido de lectura; el clic, "
        "el zoom y el resaltado funcionan igual en ambos carriles."
    )

    # Clic en el dibujo: se guarda la selección y se rerenderiza para que el
    # grupo reciba la clase de resalte. El ts evita bucles con el valor que el
    # componente conserva entre ejecuciones.
    if clicked:
        ts = clicked.get("ts")
        if (ts != st.session_state.get("br_dismissed")
                and ts != st.session_state.get("br_last_click_ts")):
            clicked_kind = str(clicked.get("kind") or "")
            clicked_id = str(clicked.get("id") or "")
            if clicked_kind and clicked_id:
                st.session_state["br_last_click_ts"] = ts
                new_key = f"{clicked_kind}:{clicked_id}"
                st.session_state["br_focus"] = _focus_payload(
                    clicked_kind, clicked_id, source="dibujo", ts=ts)
                st.session_state["br_focus_key"] = new_key
                st.session_state.pop("br_dismissed_key", None)
                try:
                    st.query_params["sel"] = new_key
                    st.query_params["contig"] = contig
                    st.query_params["pos"] = str(start)
                    if size_label in WINDOW_SIZES:
                        st.query_params["w"] = size_label
                    st.session_state["br_url_sel_seen"] = new_key
                except Exception:
                    pass
                st.rerun()

    if (focus and focus_key
            and st.session_state.get("br_dismissed_key") != focus_key):
        _render_pick(data_root, genome, tables, {
            "kind": focus_kind, "id": focus_id,
            "ts": focus.get("ts") if focus else None,
        })

    # --- distancias
    if len(wg) > 1:
        gaps = [wg.iloc[i]["start"] - wg.iloc[i - 1]["end"] - 1 for i in range(1, len(wg))]
        gaps = [g for g in gaps if g >= 0]
        if gaps:
            st.caption(f"Distancia entre genes contiguos: mediana "
                       f"{pd.Series(gaps).median():.0f} pb · mínimo {min(gaps)} pb · "
                       f"máximo {max(gaps)} pb")

    # --- solapamientos entre CDS de la ventana
    st.markdown(f"### 🔀 {OVERLAPS_TITLE}")
    with st.expander(f"🔀 {OVERLAPS_TITLE} ({len(overlaps)})",
                     expanded=not overlaps.empty):
        st.caption(
            "Pares de CDS del **mismo contig** cuyas coordenadas se cruzan en "
            "esta ventana. El solape se mide sobre los CDS completos, con "
            "coordenadas inclusivas (1-based). El marco solo se muestra si "
            "`frame`/`phase` viene en los datos; si no, aparece `unknown` y "
            "**no se infiere**. «Locus compartido» = mismo "
            "`locus_tag`/`gene_id` o solape ≥ "
            f"{LOCUS_SHARE_FRACTION:.0%} del CDS más corto. En hebras opuestas "
            "el marco se compara tal cual viene en los datos (cada CDS lo "
            "numera en su propia hebra): es un indicio, no una alineación de "
            "codones."
        )
        if overlaps.empty:
            st.info("No hay CDS solapados en esta ventana.")
        else:
            st.dataframe(
                overlap_table_view(overlaps),
                width="stretch", hide_index=True,
            )
            st.caption(overlap_summary(overlaps))
        st.warning(OVERLAP_WARNING)

    # --- vista circular complementaria
    c_len = contig_lengths(tables).get(str(contig))
    if c_len:
        st.markdown("#### 🧭 Vista circular del contig")
        st.caption(
            "El contig se dibuja como un **círculo completo**. La **ventana "
            "lineal** de arriba es la porción naranja resaltada: no es un "
            "cromosoma aparte, es el mismo contig visto entero. Los ticks "
            "marcan posiciones y el centro indica qué fracción ocupa la ventana."
        )
        help_expander("circular")
        st.markdown(circular_svg(contig, int(c_len), start, end),
                    unsafe_allow_html=True)
        st.caption(f"Longitud estimada desde la anotación: {int(c_len):,} pb. "
                   "La vista circular es complementaria; la tabla y las fichas "
                   "siguen usando la ventana lineal.")

    # --- leyenda
    st.markdown(
        "**Leyenda** · "
        f'<span style="color:{GENE_FWD}">▭</span> gen hebra + (5\'→3\' →) · '
        f'<span style="color:{GENE_REV}">▭</span> gen hebra − (5\'→3\' ←) · '
        '<span style="color:#8b90a0">▭</span> hebra desconocida · '
        '<span style="color:#b30000">●</span> horquilla clase A (fuerte) · '
        '<span style="color:#c58a00">●</span> clase B (intermedia) · '
        '<span style="color:#9aa0ad">●</span> clase C (débil) · '
        f'<span style="color:{PROMOTER}">│</span> promotor (la σ y las cajas se ven en el detalle) · '
        '<span style="color:#c58a00">▭</span> sospecha por cadena metabólica · '
        '<span style="color:#177245">▲</span> TSS medido · '
        '<span style="color:#177245;opacity:.55">▲</span> TSS alto rendimiento · '
        '<span style="color:#8a5a22">▼</span> TTS curado · '
        '<span style="color:#2e7d32">▭</span> complejo completo · '
        '<span style="color:#c58a00">▭</span> complejo incompleto',
        unsafe_allow_html=True,
    )
    if highlight_all:
        st.caption("🟡 Resaltado suave: todos los elementos de la tabla "
                   "«Todos los elementos en pantalla».")
    if focus_key:
        st.caption(f"🔴 Elemento seleccionado: `{focus_key}`. Pulsa "
                   "«Quitar resaltado» en la tabla o ✖ en la ficha para limpiarlo.")
        if st.button("✖ Quitar resaltado", key="br_clear_focus_top"):
            _clear_focus()

    def _elegir(df: pd.DataFrame, cols: list[str], key: str) -> int | None:
        """Tabla seleccionable; devuelve el índice de la fila elegida."""
        if df.empty:
            st.caption("Nada en esta ventana.")
            return None
        cols = [c for c in cols if c in df.columns] or list(df.columns)
        try:
            ev = st.dataframe(df[cols], width="stretch", hide_index=True,
                              on_select="rerun", selection_mode="single-row",
                              key=key)
            sel = ev["selection"]["rows"] if isinstance(ev, dict) else \
                  getattr(ev, "selection", {}).get("rows", [])
            if sel:
                return int(sel[0])
        except Exception:
            idx = st.selectbox("Elegir", range(len(df)), key=key + "_sb",
                               format_func=lambda i: " · ".join(
                                   str(df.iloc[i].get(c, ""))[:30] for c in cols[:3]))
            return int(idx)
        st.caption("Haz clic en una fila para ver el detalle y llevarla al dibujo.")
        return None

    def _ver(kind: str, ident: str, key: str) -> None:
        """Botón común de las pestañas: fija sel/viewport y resalta."""
        if st.button("🔎 Ver en el dibujo", key=key):
            _activate_element(data_root, genome, tables, kind, ident,
                              suspects=susp)


    # --- tabla resumen de elementos en pantalla
    st.markdown("### 🧾 Elementos en pantalla")
    st.caption(ELEMENTS_ON_SCREEN_NOTE)
    help_expander("elements")
    m1, m2 = st.columns([1, 2])
    with m1:
        st.metric("Elementos en pantalla", len(screen_elements))
    with m2:
        st.markdown(f"**{elements_counter_text(counts)}**")
        st.caption("Este contador y la tabla usan exactamente la misma fuente: "
                   "lo que se ve en la ventana lineal. No incluye el resto del genoma.")

    ver_todos = st.checkbox(
        ALL_ELEMENTS_LABEL,
        key="br_all_elements",
        help=("Incluye genes, operones, promotores, horquillas, TSS/TTS, "
              "complejos y sospechas de esta ventana; también los resalta en "
              "el dibujo superior."),
    )
    if ver_todos:
        if screen_elements.empty:
            st.caption("No hay elementos en esta ventana. Prueba a mover la "
                       "posición o a activar más pistas.")
        else:
            i = _elegir(
                screen_elements,
                ["tipo", "etiqueta", "inicio", "fin", "detalle"],
                "br_tab_all",
            )
            if i is not None:
                row = screen_elements.iloc[int(i)]
                row_key = f"{row['kind']}:{row['id']}"
                # Al elegir una fila se lleva al dibujo automáticamente. Se
                # recuerda la última fila para no repetir el rerun en cada
                # ejecución mientras la selección siga puesta.
                if (row_key != st.session_state.get("br_auto_row")
                        and row_key != st.session_state.get("br_auto_disabled")):
                    st.session_state["br_auto_row"] = row_key
                    _activate_element(data_root, genome, tables,
                                      str(row["kind"]), str(row["id"]),
                                      suspects=susp)
                st.markdown(
                    f"**{row.get('etiqueta', '')}** · {row.get('tipo', '')} · "
                    f"{int(row['inicio']):,}–{int(row['fin']):,} pb"
                )
                if str(row.get("detalle") or "").strip():
                    st.caption(str(row["detalle"]))
                b1, b2 = st.columns([1.3, 1])
                with b1:
                    if st.button("🔎 Ver en el dibujo", key="br_all_go"):
                        _activate_element(data_root, genome, tables,
                                          str(row["kind"]), str(row["id"]),
                                          suspects=susp)
                with b2:
                    if st.button("✖ Quitar resaltado", key="br_all_clear"):
                        _clear_focus()
            csv_bytes = screen_elements.to_csv(index=False).encode("utf-8")
            safe_contig = _safe_slug(contig)
            st.download_button(
                "⬇️ Exportar todos los elementos en pantalla (CSV)",
                data=csv_bytes,
                file_name=f"elementos_{genome}_{safe_contig}_{start}_{end}.csv",
                mime="text/csv",
                key="br_dl_all",
            )
            st.caption("La exportación incluye exactamente esta ventana, no el "
                       "genoma completo. Si eliges una fila, «Ver en el dibujo» "
                       "fija `sel` y, si hace falta, acerca la ventana.")
    else:
        st.caption("Activa «Todos los elementos en pantalla» para ver la tabla "
                   "combinada, exportarla y resaltar sus elementos. Las pestañas "
                   "de abajo siguen disponibles pista por pista.")

    # --- detalle por pestañas: aquí se "pincha" en cada elemento
    t_genes, t_ops, t_prom, t_hp, t_susp, t_tx, t_cx = st.tabs(
        ["🧬 Genes", "🧩 Operones", "🎯 Promotores", "🌀 Horquillas",
         "🔗 Sospechas", "📡 Transcripción", "🧩 Complejos"]
    )

    with t_genes:
        wg_view = wg.reset_index(drop=True)
        i = _elegir(wg_view,
                    ["locus_tag", "gene_name", "product", "start", "end", "strand", "ec"],
                    "br_tab_gen")
        if i is not None:
            r = wg_view.iloc[i]
            st.markdown(f"### 🧬 {r.get('gene_name') or r.get('locus_tag')}")
            a, b, c = st.columns(3)
            a.markdown(f"**Producto**  \n{r.get('product','') or '—'}")
            b.markdown(f"**Coordenadas**  \n{r.get('contig')}:{r.get('start')}–{r.get('end')} · "
                       f"hebra {'−' if str(r.get('strand')) in ('-1','-') else '+'}")
            c.markdown(f"**EC / familia**  \n{r.get('ec','') or '—'} · {r.get('fam','') or '—'}")
            ops = tables["operons"]
            if not ops.empty and "genes" in ops:
                hit = ops[ops["genes"].astype(str).str.contains(
                    str(r.get("locus_tag")), regex=False)]
                if not hit.empty:
                    o = hit.iloc[0]
                    st.caption(f"Pertenece al operón **{o.get('label') or o.get('operon_id')}** "
                               f"({o.get('n_genes')} genes, {o.get('span_bp','')} pb).")
            # promotor del operón si lo hay
            prom = tables["promoters"]
            if not prom.empty and not ops.empty and "operon_id" in prom:
                hit = ops[ops["genes"].astype(str).str.contains(
                    str(r.get("locus_tag")), regex=False)]
                if not hit.empty:
                    p2 = prom[prom["operon_id"] == hit.iloc[0].get("operon_id")]
                    if not p2.empty:
                        pr = p2.iloc[0]
                        st.markdown(f"**Promotor**: −35 `{pr.get('box35','')}` · "
                                    f"−10 `{pr.get('box10','')}` · σ **{pr.get('sigma','')}** "
                                    f"(score {pr.get('score','')})")
            gid = (r.get("locus_tag") or r.get("gene_id")
                   or r.get("gene_name"))
            if gid:
                _ver("gene", str(gid), "br_go_gen")

    with t_ops:
        help_expander("operons")
        wops = tables["operons"]
        if not wops.empty and {"contig", "start", "end"} <= set(wops.columns):
            wops = wops[(wops["contig"] == contig)
                        & (wops["end"] >= start) & (wops["start"] <= end)]
        wops_view = wops.reset_index(drop=True)
        i = _elegir(wops_view,
                    ["label", "n_genes", "start", "end", "strand", "span_bp", "confidence"],
                    "br_tab_op")
        if i is not None:
            r = wops_view.iloc[i]
            st.markdown(f"### 🧩 {r.get('label') or r.get('operon_id')}")
            a, b, c = st.columns(3)
            a.markdown(f"**Genes**  \n{r.get('n_genes')} · {r.get('span_bp','')} pb")
            b.markdown(f"**Coordenadas**  \n{r.get('start')}–{r.get('end')} · hebra "
                       f"{'−' if str(r.get('strand')) in ('-1','-') else '+'}")
            c.markdown(f"**Confianza**  \n{r.get('confidence','') or '—'}")
            st.caption(f"Genes: {r.get('gene_names') or r.get('genes','')}")
            st.caption(f"Evidencia: {r.get('evidence','') or '—'}")
            oid = r.get("operon_id") or r.get("label")
            if oid:
                _ver("operon", str(oid), "br_go_op")

    with t_prom:
        help_expander("operons")
        prom = tables["promoters"]
        if (not prom.empty and not tables["operons"].empty
                and {"contig", "start", "end", "operon_id"}
                <= set(tables["operons"].columns)
                and "operon_id" in prom.columns):
            ids = set(tables["operons"].loc[
                (tables["operons"]["contig"] == contig)
                & (tables["operons"]["end"] >= start)
                & (tables["operons"]["start"] <= end), "operon_id"].astype(str))
            prom = prom[prom["operon_id"].astype(str).isin(ids)]
        prom_view = prom.reset_index(drop=True)
        i = _elegir(prom_view,
                    ["operon_id", "sigma", "box35", "box10", "spacer", "score", "confidence"],
                    "br_tab_prom")
        if i is not None:
            r = prom_view.iloc[i]
            st.markdown(f"### 🎯 Promotor de {r.get('operon_id')}")
            st.markdown(f"**Caja −35** `{r.get('box35','')}` (desajustes {r.get('mm35','')})  \n"
                        f"**Caja −10** `{r.get('box10','')}` (desajustes {r.get('mm10','')})  \n"
                        f"**Espaciador** {r.get('spacer','')} pb · **σ** {r.get('sigma','')}  \n"
                        f"**Score** {r.get('score','')} · **Confianza** {r.get('confidence','')}")
            st.caption("Un promotor es la marca donde la ARN polimerasa empieza a "
                       "copiar. La caja −10 y la −35 son las dos secuencias que "
                       "reconoce el factor sigma.")
            pid = r.get("operon_id") or r.get("label")
            if pid:
                _ver("promoter", str(pid), "br_go_prom")

    with t_hp:
        help_expander("hairpins")
        hp = tables["hairpins"]
        if not hp.empty and {"contig", "start", "end"} <= set(hp.columns):
            hp = hp[(hp["contig"] == contig)
                    & (hp["end"] >= start) & (hp["start"] <= end)]
        hp_view = hp.reset_index(drop=True)
        i = _elegir(hp_view,
                    ["start", "end", "class", "dG", "stem", "gc", "polyU",
                     "left_gene", "right_gene"],
                    "br_tab_hp")
        if i is not None:
            r = hp_view.iloc[i]
            cls = str(r.get("class", ""))
            st.markdown(f"### 🌀 Horquilla en {r.get('start')}–{r.get('end')}")
            st.markdown(f"**Clase** {cls} — {HAIRPIN_MEANING.get(cls, '')}  \n"
                        f"**ΔG** {r.get('dG','')} kcal/mol · **tallo** {r.get('stem','')} nt · "
                        f"**GC** {r.get('gc','')} · **poli-U** {r.get('polyU','')}  \n"
                        f"**Entre** {r.get('left_gene','')} y {r.get('right_gene','')}")
            st.caption("La clase compara el tallo con el fondo de GC del genoma. "
                       "A: estructura clara; B: intermedia, no se distingue del azar; "
                       "C: débil.")
            hid = f"{contig}:{int(float(r['start']))}-{int(float(r['end']))}"
            _ver("hairpin", hid, "br_go_hp")

    with t_susp:
        if susp.empty:
            st.caption("Sin sospechas de cadena metabólica en esta ventana.")
        else:
            susp_view = susp.reset_index(drop=True)
            i = _elegir(susp_view,
                        ["n_genes", "products", "reactions", "gap_mean",
                         "shared_metabolites", "confidence"],
                        "br_tab_susp")
            if i is not None:
                r = susp_view.iloc[i]
                st.markdown(f"### 🔗 Sospecha: {r.get('n_genes')} genes contiguos")
                st.markdown(f"**Genes**: {r.get('genes','')}  \n"
                            f"**Enzimas**: {r.get('products','')}  \n"
                            f"**Reacciones**: {r.get('reactions','')}  \n"
                            f"**Metabolitos compartidos**: {r.get('shared_metabolites','')}")
                st.caption("Estas enzimas tocan los mismos metabolitos, así que "
                           "probablemente participan de la misma ruta. Al estar "
                           "contiguas en el ADN, es **sospecha** de operón "
                           "metabólico. No es prueba: hay rutas contiguas que no "
                           "son operones y operones sin ruta asociada.")
                sid = f"{int(float(r['start']))}-{int(float(r['end']))}"
                _ver("suspect", sid, "br_go_susp")

    with t_tx:
        wt = (tss_df[(tss_df.pos >= start) & (tss_df.pos <= end)]
              if not tss_df.empty else tss_df)
        wtt = (tts_df[(tts_df.pos >= start) & (tts_df.pos <= end)]
               if not tts_df.empty else tts_df)
        if wt.empty and wtt.empty:
            st.caption("Sin TSS ni TTS anotados en esta ventana.")
        else:
            st.caption(
                "Marcas de **transcripción**: dónde empieza (TSS) y dónde termina "
                "(TTS) el ARN. Las sólidas están curadas; las atenuadas vienen de "
                "alto rendimiento. Es el dato que decide las fronteras cuando existe."
            )
            if not wt.empty:
                st.markdown(f"**TSS en la ventana ({len(wt)})**")
                wt_view = wt.reset_index(drop=True)
                i = _elegir(wt_view, ["pos", "strand", "source", "detail"],
                            "br_tab_tss")
                if i is not None:
                    pos = int(float(wt_view.iloc[i]["pos"]))
                    st.caption(f"TSS en {pos:,} pb · {wt_view.iloc[i].get('source', '')}")
                    _ver("tss", str(pos), "br_go_tss")
            if not wtt.empty:
                st.markdown(f"**TTS en la ventana ({len(wtt)})**")
                wtt_view = wtt.reset_index(drop=True)
                i = _elegir(wtt_view, ["pos", "strand", "source", "detail"],
                            "br_tab_tts")
                if i is not None:
                    pos = int(float(wtt_view.iloc[i]["pos"]))
                    st.caption(f"TTS en {pos:,} pb · {wtt_view.iloc[i].get('source', '')}")
                    _ver("tts", str(pos), "br_go_tts")

    with t_cx:
        cxs = tables["complexes"]
        if not cxs.empty and {"contig", "start", "end"} <= set(cxs.columns):
            cxs = cxs[(cxs["contig"] == contig)
                      & (cxs["end"] >= start) & (cxs["start"] <= end)]
        cxs_view = cxs.reset_index(drop=True)
        i = _elegir(cxs_view,
                    ["complex_id", "name", "status", "start", "end",
                     "subunits", "confidence", "limitation"],
                    "br_tab_cx")
        if i is not None:
            r = cxs_view.iloc[i]
            st.markdown(f"### 🧩 {r.get('name') or r.get('complex_id')}")
            st.caption(f"{r.get('status','')} · confianza {r.get('confidence','')}")
            st.markdown(f"**Coordenadas** {r.get('start')}–{r.get('end')} pb  \n"
                        f"**Subunidades** {r.get('subunits','') or '—'}  \n"
                        f"**Evidencia** {r.get('evidence','') or '—'}")
            if r.get("limitation"):
                st.caption(r.get("limitation"))
            cid = r.get("complex_id")
            if cid:
                _ver("complex", str(cid), "br_go_cx")

    st.info(
        "**Terminador ≠ atenuador.** Una horquilla al final de un operón lo "
        "termina. Una horquilla *dentro* del operón puede ser un atenuador: "
        "corta el transcrito pero el operón continúa. Por eso una horquilla no "
        "basta para decidir dónde acaba un operón."
    )
