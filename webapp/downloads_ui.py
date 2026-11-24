"""Centro de descargas: reúne los resultados que un usuario puede reutilizar.

La idea es que no haya que buscar archivos por el disco: la app dice qué hay,
qué contiene, de dónde sale, cuánto pesa y con qué licencia se puede usar.
También arma un paquete ZIP pequeño con los resultados principales y la
información de citación.
"""
from __future__ import annotations

import io
import zipfile
from datetime import datetime
from pathlib import Path

from i18n import tr

#: (título, ruta relativa, descripción, licencia/cita)
CATALOG = [
    ("Brechas metabólicas", "data/metabolism/real_gaps.tsv",
     "Una fila por brecha: reacción, enzima esperada, estado y soporte.",
     "Derivado de Rhea/UniProt · cita Rhea (Bansal et al., 2022) y UniProt"),
    ("Interpretación de brechas", "data/metabolism/gap_interpretation.tsv",
     "Explicación por familia: ausencia real, transporte, mapeo o falso positivo.",
     "Datos derivados del proyecto"),
    ("Cierres de brecha", "data/metabolism/gap_closures.tsv",
     "Cierres propuestos por dominio/homología, con veredicto.",
     "Datos derivados del proyecto"),
    ("Cola de revisión", "data/metabolism/REVIEW_QUEUE.tsv",
     "Brechas priorizadas para revisión manual (alta/media/baja).",
     "Datos derivados del proyecto"),
    ("Completitud por ruta", "data/metabolism/pathway_completeness.tsv",
     "Cobertura de cada ruta KEGG/Rhea en los 400 genomas.",
     "Derivado de KEGG/Rhea · uso académico; no redistribuir KEGG"),
    ("Auditoría de mapeos", "data/metabolism/mapping_audit.tsv",
     "Cada mapeo con su veredicto (support/weak/reject) y contexto.",
     "Datos derivados del proyecto"),
    ("Sospechas de operón", "data/metabolism/operon_suspects.tsv",
     "Genes contiguos que comparten metabolitos; es una sospecha, no una prueba.",
     "Datos derivados del proyecto"),
    ("Operones por reglas (K-12)", "data/operon_finder/ecoli_k12_operons_rules.tsv",
     "Predicción determinista de operones para E. coli K-12.",
     "Datos derivados del proyecto; verdad de referencia RegulonDB (no redistribuir)"),
    ("Paridad de anotación", "reports/annotation_parity.md",
     "Comparación de CDS/tRNA/rRNA contra una referencia; no mide función.",
     "Auditoría interna · datos derivados"),
    ("Cobertura móvil/AMR", "reports/mobile_coverage.md",
     "Cobertura real de móviles, defensa, CRISPR, AMR y virulencia.",
     "Auditoría interna · datos derivados"),
    ("Cierre de brechas", "reports/gap_closure_audit.md",
     "Estados de brecha y límites de cierre; ninguna se cierra sin evidencia.",
     "Auditoría interna · datos derivados"),
    ("Benchmark por capas", "reports/layer_benchmark.md",
     "Métricas internas de A/B/C, membrana, complejos y oriC.",
     "Auditoría interna · no es benchmark externo"),
    ("Benchmark externo", "reports/external_benchmark.md",
     "Harness contra Bakta/Prokka/DFAST; not_evaluable si no están instalados.",
     "Harness CPU · sin resultados externos si faltan binarios"),
    ("Holdout de complejos", "reports/complex_holdout.md",
     "Holdout y cobertura de complejos; no hay gold independiente fuera de K-12.",
     "Auditoría interna · hipótesis"),
    ("Calibración de clusters", "reports/model_clusters_calibration.md",
     "Priorización P1/P2 de los 281.090 clusters del modelo K-12.",
     "Auditoría interna · no es probabilidad"),
    ("Benchmark de operones", "reports/operon_benchmark.md",
     "Operones por reglas contra RegulonDB K-12; no es holdout ciego.",
     "Auditoría interna · gold K-12"),
    ("Evaluación del modelo", "reports/eval_comprehensive.json",
     "Métricas del modelo de ADN. Ojo: corrida corta en CPU, no es una evaluación completa.",
     "Datos derivados del proyecto"),
    ("Estado de la evaluación", "reports/README.md",
     "Qué se evaluó, qué no y por qué.",
     "Documentación del proyecto (MIT)"),
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
    try:
        if path.stat().st_size > limit_mb * 1024 * 1024:
            return None
        with path.open("rb") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except OSError:
        return None


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


def render(root: Path) -> None:
    import pandas as pd
    import streamlit as st

    st.subheader("📥 " + tr("Reportes y descargas"))
    st.caption(tr(
        "Aquí están los resultados que puedes reutilizar en tu investigación. "
        "Cada archivo dice qué contiene, de dónde sale y con qué licencia. "
        "Los datos de terceros con licencia académica no se incluyen: se "
        "descargan en tu equipo y no se redistribuyen."))

    rows = []
    for title, rel, desc, lic in CATALOG:
        path = root / rel
        if not path.exists():
            continue
        n = _rows(path)
        rows.append({
            tr("Archivo"): rel,
            tr("Qué contiene"): tr(desc),
            tr("Filas"): n if n is not None else "—",
            tr("Tamaño"): _human(path.stat().st_size),
            tr("Licencia / cita"): tr(lic),
        })
    if not rows:
        st.info(tr("Todavía no hay resultados calculados en este equipo. "
                   "Anota un FASTA o genera los datos del proyecto."))
    else:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption(tr("Los archivos están en data/ y reports/; también puedes "
                      "descargarlos uno a uno abajo."))
        cols = st.columns(3)
        for i, (title, rel, desc, lic) in enumerate(CATALOG):
            path = root / rel
            if not path.exists():
                continue
            with cols[i % 3]:
                data = path.read_bytes()
                st.download_button(
                    "⬇️ " + tr(title),
                    data,
                    file_name=path.name,
                    mime="text/tab-separated-values" if path.suffix == ".tsv"
                    else "application/json" if path.suffix == ".json"
                    else "text/markdown",
                    help=tr(desc), key=f"dl_{rel}",
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
