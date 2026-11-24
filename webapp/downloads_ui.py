"""Centro de descargas: reúne los resultados que un usuario puede reutilizar.

La idea es que no haya que buscar archivos por el disco: la app dice qué hay,
qué contiene, de dónde sale, cuánto pesa y con qué licencia se puede usar.
También arma un paquete ZIP pequeño con los resultados principales y la
información de citación.

Cada entrada del catálogo declara, además, el **nivel de evidencia** que
sostiene el archivo:

``interno``
    Auditoría o dato derivado del propio proyecto. No es una validación
    experimental independiente.
``not_evaluable``
    El informe describe un método que no se pudo evaluar con los datos o
    binarios disponibles. No es un cero ni un fracaso del método.
``externo``
    Resultado producido por una herramienta o un gold externo. Hoy no hay
    entradas de este tipo: el benchmark externo quedó ``not_evaluable`` porque
    Bakta, Prokka y DFAST no están instalados.

``module_reports(module, root)`` devuelve solo las entradas de un módulo para
que ``app.py`` pueda mostrarlas junto a sus resultados, sin duplicar la lista.
"""
from __future__ import annotations

import io
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

from i18n import tr

#: Niveles de evidencia admitidos. El valor se guarda tal cual en el catálogo.
EVIDENCE_LEVELS = ("interno", "not_evaluable", "externo")

#: Cómo se muestra el nivel en la interfaz (el valor interno no cambia).
EVIDENCE_LABELS = {
    "interno": "Interno",
    "not_evaluable": "No evaluable",
    "externo": "Externo",
}

#: Módulo canónico -> alias aceptados por :func:`module_reports`.
#: ``complex`` no es un módulo autodescubierto de ``modules/``, pero la GUI lo
#: usa para la sección de complejos; por eso se mantiene como clave propia.
MODULE_ALIASES = {
    "annotation": "annotation",
    "anotacion": "annotation",
    "anotación": "annotation",
    "mobile": "mobile",
    "movil": "mobile",
    "móvil": "mobile",
    "mobile_coverage": "mobile",
    "metabolism": "metabolism",
    "metabolismo": "metabolism",
    "gap": "metabolism",
    "gaps": "metabolism",
    "gap_closure": "metabolism",
    "complex": "complex",
    "complexes": "complex",
    "complejo": "complex",
    "complejos": "complex",
    "model": "complex",
    "clusters": "complex",
    "holdout": "complex",
    "regulation": "regulation",
    "regulacion": "regulation",
    "regulación": "regulation",
    "operon": "regulation",
    "operons": "regulation",
    "operones": "regulation",
    "operon_benchmark": "regulation",
    "validation": "validation",
    "validacion": "validation",
    "validación": "validation",
    "benchmark": "validation",
    "benchmarks": "validation",
    "layers": "validation",
    "membrane": "membrane",
    "membrana": "membrane",
    "structure": "structure",
    "estructura": "structure",
    "taxonomy": "taxonomy",
    "ghost": "ghost",
    "targets": "targets",
}

#: Catálogo de artefactos descargables.
#:
#: Cada entrada es un dict con:
#: - ``title``: título visible (clave de traducción);
#: - ``path``: ruta relativa a la raíz del proyecto;
#: - ``description``: qué contiene y qué no mide;
#: - ``license``: licencia y cita que hay que respetar;
#: - ``evidence``: ``interno``, ``not_evaluable`` o ``externo``;
#: - ``modules``: módulos donde se puede mostrar in situ (tupla, puede ir vacía).
CATALOG: list[dict[str, Any]] = [
    {
        "title": "Brechas metabólicas",
        "path": "data/metabolism/real_gaps.tsv",
        "description": "Una fila por brecha: reacción, enzima esperada, estado y soporte.",
        "license": "Derivado de Rhea/UniProt · cita Rhea (Bansal et al., 2022) y UniProt",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Interpretación de brechas",
        "path": "data/metabolism/gap_interpretation.tsv",
        "description": "Explicación por familia: ausencia real, transporte, mapeo o falso positivo.",
        "license": "Datos derivados del proyecto",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Cierres de brecha",
        "path": "data/metabolism/gap_closures.tsv",
        "description": "Cierres propuestos por dominio/homología, con veredicto.",
        "license": "Datos derivados del proyecto",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Cola de revisión",
        "path": "data/metabolism/REVIEW_QUEUE.tsv",
        "description": "Brechas priorizadas para revisión manual (alta/media/baja).",
        "license": "Datos derivados del proyecto",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Completitud por ruta",
        "path": "data/metabolism/pathway_completeness.tsv",
        "description": "Cobertura de cada ruta KEGG/Rhea en los 400 genomas.",
        "license": "Derivado de KEGG/Rhea · uso académico; no redistribuir KEGG",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Auditoría de mapeos",
        "path": "data/metabolism/mapping_audit.tsv",
        "description": "Cada mapeo con su veredicto (support/weak/reject) y contexto.",
        "license": "Datos derivados del proyecto",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Sospechas de operón",
        "path": "data/metabolism/operon_suspects.tsv",
        "description": "Genes contiguos que comparten metabolitos; es una sospecha, no una prueba.",
        "license": "Datos derivados del proyecto",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Operones por reglas (K-12)",
        "path": "data/operon_finder/ecoli_k12_operons_rules.tsv",
        "description": "Predicción determinista de operones para E. coli K-12.",
        "license": "Datos derivados del proyecto; verdad de referencia RegulonDB (no redistribuir)",
        "evidence": "interno",
        "modules": (),
    },
    {
        "title": "Paridad de anotación",
        "path": "reports/annotation_parity.md",
        "description": "Comparación de CDS/tRNA/rRNA contra una referencia; no mide función.",
        "license": "Auditoría interna · datos derivados",
        "evidence": "interno",
        "modules": ("annotation",),
    },
    {
        "title": "Paridad de anotación (JSON)",
        "path": "reports/annotation_parity.json",
        "description": ("Misma paridad en JSON: metadatos, conteos, coordenadas, "
                        "productos, ORFs pequeños, ncRNA y 200 diferencias."),
        "license": "Auditoría interna · datos derivados; cita el GFF3 de referencia",
        "evidence": "interno",
        "modules": ("annotation",),
    },
    {
        "title": "Cobertura móvil/AMR",
        "path": "reports/mobile_coverage.md",
        "description": "Cobertura real de móviles, defensa, CRISPR, AMR y virulencia.",
        "license": "Auditoría interna · datos derivados",
        "evidence": "interno",
        "modules": ("mobile",),
    },
    {
        "title": "Cobertura móvil/AMR (JSON)",
        "path": "reports/mobile_coverage.json",
        "description": ("Datos completos por genoma: inventario, cobertura, "
                        "brechas y advertencias de la auditoría móvil/AMR."),
        "license": ("Auditoría interna · datos derivados; fuentes móviles/AMR "
                    "con sus propias licencias (ver DATA_NOTICE.md)"),
        "evidence": "interno",
        "modules": ("mobile",),
    },
    {
        "title": "Cierre de brechas",
        "path": "reports/gap_closure_audit.md",
        "description": "Estados de brecha y límites de cierre; ninguna se cierra sin evidencia.",
        "license": "Auditoría interna · datos derivados",
        "evidence": "interno",
        "modules": ("metabolism",),
    },
    {
        "title": "Cierre de brechas (JSON)",
        "path": "reports/gap_closure_audit.json",
        "description": ("Instantánea completa de 8.655 brechas: estados, "
                        "evidencia, bloqueadores, asociaciones y detalle no evaluado."),
        "license": "Auditoría interna · datos derivados de Rhea/UniProt/ChEBI; cita Rhea y UniProt",
        "evidence": "interno",
        "modules": ("metabolism",),
    },
    {
        "title": "Benchmark por capas",
        "path": "reports/layer_benchmark.md",
        "description": "Métricas internas de A/B/C, membrana, complejos y oriC.",
        "license": "Auditoría interna · no es benchmark externo",
        "evidence": "interno",
        "modules": ("annotation", "membrane", "complex", "validation"),
    },
    {
        "title": "Benchmark por capas (JSON)",
        "path": "reports/layer_benchmark.json",
        "description": ("Datos por capa y sección del benchmark A/B/C, membrana, "
                        "complejos y oriC, con límites declarados."),
        "license": ("Auditoría interna · datos locales; cita este repositorio y "
                    "el gold de cada capa (ComplexPortal, DoriC)"),
        "evidence": "interno",
        "modules": ("annotation", "membrane", "complex", "validation"),
    },
    {
        "title": "Benchmark externo",
        "path": "reports/external_benchmark.md",
        "description": "Harness contra Bakta/Prokka/DFAST; not_evaluable si no están instalados.",
        "license": "Harness CPU · sin resultados externos si faltan binarios",
        "evidence": "not_evaluable",
        "modules": ("annotation", "validation"),
    },
    {
        "title": "Benchmark externo (JSON)",
        "path": "reports/external_benchmark.json",
        "description": ("Entorno, referencia y resultados del harness externo; "
                        "en esta corrida Bakta/Prokka/DFAST no están instalados "
                        "y todo queda not_evaluable."),
        "license": ("Harness interno · sin resultados externos; si se ejecuta, "
                    "cita Bakta, Prokka y DFAST con sus licencias"),
        "evidence": "not_evaluable",
        "modules": ("annotation", "validation"),
    },
    {
        "title": "Holdout de complejos",
        "path": "reports/complex_holdout.md",
        "description": "Holdout y cobertura de complejos; no hay gold independiente fuera de K-12.",
        "license": "Auditoría interna · hipótesis",
        "evidence": "not_evaluable",
        "modules": ("complex",),
    },
    {
        "title": "Holdout de complejos (JSON)",
        "path": "reports/complex_holdout.json",
        "description": ("Datos del holdout de complejos: procedencia del modelo, "
                        "catálogos K-12, evaluaciones y limitaciones."),
        "license": "Auditoría interna · ComplexPortal K-12; cita ComplexPortal",
        "evidence": "not_evaluable",
        "modules": ("complex",),
    },
    {
        "title": "Calibración de clusters",
        "path": "reports/model_clusters_calibration.md",
        "description": "Priorización P1/P2/P3 de los 281.090 clusters del modelo K-12.",
        "license": "Auditoría interna · no es probabilidad",
        "evidence": "interno",
        "modules": ("complex",),
    },
    {
        "title": "Calibración de clusters (JSON)",
        "path": "reports/model_clusters_calibration.json",
        "description": ("Agregados, evidencia y priorización P1/P2/P3 de los "
                        "281.090 clusters; la confianza no es una probabilidad."),
        "license": ("Auditoría interna · datos derivados K-12; cita este "
                    "repositorio; no es validación experimental"),
        "evidence": "interno",
        "modules": ("complex",),
    },
    {
        "title": "Benchmark de operones",
        "path": "reports/operon_benchmark.md",
        "description": "Operones por reglas contra RegulonDB K-12; no es holdout ciego.",
        "license": "Auditoría interna · gold K-12",
        "evidence": "interno",
        "modules": ("regulation", "validation"),
    },
    {
        "title": "Benchmark de operones (JSON)",
        "path": "reports/operon_benchmark.json",
        "description": ("Conteos, métricas por hebra, promotores y advertencias "
                        "del benchmark de operones contra RegulonDB K-12."),
        "license": "Auditoría interna · gold RegulonDB (uso académico, no redistribuir)",
        "evidence": "interno",
        "modules": ("regulation", "validation"),
    },
    {
        "title": "Evaluación del modelo",
        "path": "reports/eval_comprehensive.json",
        "description": "Métricas del modelo de ADN. Ojo: corrida corta en CPU, no es una evaluación completa.",
        "license": "Datos derivados del proyecto",
        "evidence": "not_evaluable",
        "modules": ("validation",),
    },
    {
        "title": "Estado de la evaluación",
        "path": "reports/README.md",
        "description": "Qué se evaluó, qué no y por qué.",
        "license": "Documentación del proyecto (MIT)",
        "evidence": "interno",
        "modules": ("validation",),
    },
]

BUNDLE = [
    "data/metabolism/real_gaps.tsv",
    "data/metabolism/gap_interpretation.tsv",
    "data/metabolism/REVIEW_QUEUE.tsv",
    "data/metabolism/mapping_audit.tsv",
    "data/metabolism/operon_suspects.tsv",
    "data/operon_finder/ecoli_k12_operons_rules.tsv",
    "reports/README.md",
]

README_BUNDLE = """# Paquete de resultados — Anotador de ADN bacteriano

Generado por la interfaz el {date}.

Contenido:
- real_gaps.tsv: brechas metabólicas detectadas.
- gap_interpretation.tsv: explicación por familia.
- REVIEW_QUEUE.tsv: cola de revisión priorizada.
- mapping_audit.tsv: auditoría de mapeos EC/Rhea.
- operon_suspects.tsv: sospechas de operón metabólico.
- ecoli_k12_operons_rules.tsv: operones predichos por reglas.
- reports/README.md: estado de la evaluación del modelo.

Cómo interpretar:
- Una brecha NO prueba que falte una enzima. Puede ser un gen no anotado,
  transporte, ruta alternativa o un falso positivo.
- Cada fila declara evidencia, confianza y estado. No hay `confirmed` sin
  validación experimental, EC o contexto.

Citación:
Rojas Bravo, C. (2026). Anotador de ADN bacteriano (v1.1.3). https://doi.org/10.5281/zenodo.24469732. Ingeniero en Biotecnología, titulado de la Universidad de Santiago de Chile (USACH); sin afiliación institucional actual.
Asistencia de desarrollo: DeepSeek Harness con DeepSeek, noviembre de 2026.

Cita también Rhea, UniProt, KEGG y RegulonDB según corresponda. Ver
THIRD_PARTY_LICENSES.md y DATA_NOTICE.md.

Licencia: el código es MIT; los datos derivados conservan las obligaciones de
sus fuentes. No redistribuyas RegulonDB ni KEGG.
"""


def _human(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{size} B"
        size /= 1024
    return f"{size:.1f} GB"


def _rows(path: Path, limit_mb: float = 8.0) -> int | None:
    """Filas de un TSV/CSV; ``None`` para informes o archivos demasiado grandes."""
    if path.suffix.lower() not in {".tsv", ".tab", ".csv"}:
        return None
    try:
        if path.stat().st_size > limit_mb * 1024 * 1024:
            return None
        with path.open("rb") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except OSError:
        return None


def _resolve(entry: dict[str, Any], root: Path) -> dict[str, Any]:
    """Copia una entrada y resuelve existencia, tamaño y filas contra ``root``."""
    path = root / entry["path"]
    try:
        exists = path.is_file()
        size = path.stat().st_size if exists else None
    except OSError:
        exists, size = False, None
    info = dict(entry)
    info.update({
        "abs_path": str(path),
        "exists": exists,
        "size_bytes": size,
        "size_human": _human(size) if size is not None else "—",
        "rows": _rows(path) if exists else None,
    })
    return info


def catalog_entries(root: Path | str, *, include_missing: bool = False) -> list[dict[str, Any]]:
    """Catálogo con datos de disco.

    Por defecto **filtra** las entradas cuya ruta no existe, para no ofrecer
    descargas rotas. Con ``include_missing=True`` las devuelve con
    ``exists=False``, ``size_bytes=None`` y ``size_human="—"``.
    """
    root = Path(root)
    entries = [_resolve(entry, root) for entry in CATALOG]
    if include_missing:
        return entries
    return [entry for entry in entries if entry["exists"]]


def module_reports(module: str, root: Path | str,
                   *, include_missing: bool = False) -> list[dict[str, Any]]:
    """Reportes del catálogo asociados a un módulo.

    Acepta los ids de ``modules/`` (``annotation``, ``mobile``, ``metabolism``,
    ``regulation``, ``membrane``, ``structure``, ``validation``, etc.) y alias
    de la GUI (``complex``, ``gap``, ``operon``, ``layers``). Un módulo
    desconocido, o uno conocido sin reportes, devuelve ``[]``. Por defecto solo
    devuelve archivos existentes bajo ``root``.

    Cada elemento agrega ``module`` (nombre canónico) y los datos de disco de
    :func:`catalog_entries` (``abs_path``, ``exists``, ``size_bytes``,
    ``size_human``, ``rows``).
    """
    key = MODULE_ALIASES.get(str(module).strip().lower())
    if key is None:
        return []
    reports = []
    for entry in catalog_entries(root, include_missing=include_missing):
        if key in entry.get("modules", ()):
            reports.append({**entry, "module": key})
    return reports


def render_module_reports(module: str, root: Path | str, *, key_prefix: str = "") -> None:
    """Muestra descargas de los reportes asociados a un módulo, si existen."""
    import streamlit as st
    entries = module_reports(module, root)
    if not entries:
        return
    with st.expander("📥 Reportes de este módulo"):
        for entry in entries:
            path = Path(entry["abs_path"])
            st.download_button(
                "⬇️ " + tr(entry["title"]),
                path.read_bytes(), file_name=path.name,
                mime=_mime(path),
                help=tr(entry.get("description", "")),
                key=f"{key_prefix}_dl_{entry['path']}",
            )


def _bundle(root: Path) -> bytes | None:
    buf = io.BytesIO()
    files = [rel for rel in BUNDLE if (root / rel).exists()]
    if not files:
        return None
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("LEEME.md", README_BUNDLE.format(
            date=datetime.now().strftime("%Y-%m-%d")))
        for rel in ("CITATION.cff", "LICENSE", "DATA_NOTICE.md",
                    "THIRD_PARTY_LICENSES.md"):
            if (root / rel).exists():
                zf.write(root / rel, rel)
        for rel in files:
            zf.write(root / rel, rel.split("/")[-1])
    return buf.getvalue()


def _mime(path: Path) -> str:
    return {
        ".tsv": "text/tab-separated-values",
        ".tab": "text/tab-separated-values",
        ".csv": "text/csv",
        ".json": "application/json",
        ".md": "text/markdown",
    }.get(path.suffix.lower(), "application/octet-stream")


def render(root: Path) -> None:
    import pandas as pd
    import streamlit as st

    st.subheader("📥 " + tr("Reportes y descargas"))
    st.caption(tr(
        "Aquí están los resultados que puedes reutilizar en tu investigación. "
        "Cada archivo dice qué contiene, de dónde sale y con qué licencia. "
        "Los datos de terceros con licencia académica no se incluyen: se "
        "descargan en tu equipo y no se redistribuyen."))

    entries = catalog_entries(root)
    rows = []
    for entry in entries:
        rows.append({
            tr("Archivo"): entry["path"],
            tr("Qué contiene"): tr(entry["description"]),
            tr("Filas"): entry["rows"] if entry["rows"] is not None else "—",
            tr("Tamaño"): entry["size_human"],
            tr("Evidencia"): tr(EVIDENCE_LABELS[entry["evidence"]]),
            tr("Licencia / cita"): tr(entry["license"]),
        })
    if not rows:
        st.info(tr("Todavía no hay resultados calculados en este equipo. "
                   "Anota un FASTA o genera los datos del proyecto."))
    else:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption(tr("Los archivos están en data/ y reports/; también puedes "
                      "descargarlos uno a uno abajo."))
        cols = st.columns(3)
        for i, entry in enumerate(entries):
            path = root / entry["path"]
            with cols[i % 3]:
                st.download_button(
                    "⬇️ " + tr(entry["title"]),
                    path.read_bytes(),
                    file_name=path.name,
                    mime=_mime(path),
                    help=tr(entry["description"]),
                    key=f"dl_{entry['path']}",
                )

    st.divider()
    st.markdown("### " + tr("Paquete para investigación"))
    st.caption(tr("ZIP con los resultados principales, la citación y las "
                  "licencias. Está pensado para adjuntar a un informe o "
                  "repositorio, no para sustituir las fuentes originales."))
    bundle = _bundle(root)
    if bundle:
        st.download_button(
            "⬇️ " + tr("Descargar paquete (ZIP)"), bundle,
            file_name="anotador_adn_resultados.zip",
            mime="application/zip", key="dl_bundle",
        )
        mb = len(bundle) / 1024 / 1024
        st.caption(f"{mb:.1f} MB" if mb >= 0.1 else f"{len(bundle)/1024:.0f} KB")
    else:
        st.info(tr("No hay archivos suficientes para armar el paquete."))

    st.divider()
    with st.expander("📚 " + tr("Cómo citar y qué licencia tiene cada cosa")):
        st.markdown(tr(
            "El **código** es MIT. Los **datos derivados** conservan las "
            "obligaciones de sus fuentes: Rhea y UniProt son CC BY; KEGG y "
            "RegulonDB son de uso académico y **no se pueden redistribuir**. "
            "Si publicas una tabla, cita la fuente original y este repositorio. "
            "El archivo CITATION.cff tiene la referencia lista."))
        for rel in ("CITATION.cff", "LICENSE", "DATA_NOTICE.md",
                    "THIRD_PARTY_LICENSES.md"):
            path = root / rel
            if path.exists():
                st.download_button("⬇️ " + rel, path.read_bytes(),
                                   file_name=rel, key=f"doc_{rel}")
