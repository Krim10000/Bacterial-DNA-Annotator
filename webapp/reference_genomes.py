"""Panel de genomas de referencia para benchmark.

Los genomas de referencia son un panel **separado** de los 400 libres. El
sistema anota el FASTA y aquí se compara el conteo de features con la
anotación curada local. No se usa el gold para entrenar ni se declara
superioridad.
"""
from __future__ import annotations

import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Any

MANIFEST = "reference_genomes.json"


def load_manifest(data_root: Path | str) -> list[dict[str, Any]]:
    path = Path(data_root) / MANIFEST
    if not path.exists():
        # Copia empaquetada para instalaciones sin `data/` (release).
        path = Path(__file__).with_name(MANIFEST)
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return [r for r in data if isinstance(r, dict)] if isinstance(data, list) else []


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def count_gff(path: Path) -> dict[str, int]:
    """Cuenta CDS/tRNA/rRNA en GFF3; no interpreta pseudo."""
    counts: Counter[str] = Counter()
    if not path.exists():
        return {}
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                counts[f[2]] += 1
    return dict(counts)


def compare_genome(genome: str, data_root: Path | str,
                   out_root: Path | str) -> dict[str, Any]:
    """Compara el resumen del sistema con el GFF de referencia."""
    data_root, out_root = Path(data_root), Path(out_root)
    rec = next((r for r in load_manifest(data_root) if r.get("id") == genome), None)
    if rec is None:
        return {"genome": genome, "status": "not_available",
                "reason": "genoma no está en reference_genomes.json", "rows": []}
    summary = out_root / genome / f"{genome}.summary.json"
    if not summary.exists():
        summary = out_root / genome / "summary.json"
    if not summary.exists():
        return {"genome": genome, "status": "not_available",
                "reason": "sin resumen anotado por el sistema", "rows": [],
                "record": rec}
    try:
        system = json.loads(summary.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return {"genome": genome, "status": "invalid",
                "reason": f"resumen ilegible: {exc}", "rows": [], "record": rec}
    ref_path = (data_root / str(rec.get("annotation_path") or "")
                if rec.get("annotation_path") else Path(""))
    if not ref_path or not ref_path.exists():
        # fallbacks conocidos
        for cand in (data_root / "annotations" / f"{genome}.gff.gz",
                     data_root / "dbtbs" / "BSGatlas_v1.1.gff"):
            if cand.exists():
                ref_path = cand
                break
    ref_counts = count_gff(ref_path) if ref_path else {}
    sys_counts = (system.get("features") or {})
    rows = []
    for kind in ("CDS", "tRNA", "rRNA"):
        a, b = int(sys_counts.get(kind, 0) or 0), int(ref_counts.get(kind, 0) or 0)
        rows.append({"tipo": kind, "sistema": a, "referencia": b, "delta": a - b})
    return {"genome": genome, "status": "ok", "rows": rows,
            "system": system, "reference_counts": ref_counts,
            "reference_path": str(ref_path) if ref_path else "",
            "record": rec,
            "limitation": "Comparación de conteos; no evalúa función ni exactitud "
                          "de coordenadas. El gold hereda los criterios de la fuente."}


def render(data_root: Path | str, out_root: Path | str) -> None:
    import streamlit as st
    import pandas as pd

    st.subheader("🧬 Genomas de referencia")
    st.caption(
        "Panel separado de los 400 libres. Sirve para comparar la anotación "
        "del sistema contra una anotación curada; no se usa para entrenar."
    )
    rows = load_manifest(data_root)
    if not rows:
        st.warning("No hay `data/reference_genomes.json`.")
        return
    df = pd.DataFrame(rows)
    c1, c2 = st.columns([1, 2])
    level = c1.selectbox("Nivel", ["Todos", 1, 2], key="ref_level")
    if level != "Todos":
        df = df[df["level"] == level]
    st.dataframe(df[["level", "species", "strain", "resources", "status", "note"]],
                 width="stretch", hide_index=True)
    available = df[df["status"] == "available"]["id"].tolist()
    if not available:
        st.info("No hay genomas de referencia anotados todavía.")
        return
    selected = st.selectbox("Comparar genoma disponible", available,
                            key="ref_genome")
    comp = compare_genome(selected, data_root, out_root)
    if comp["status"] != "ok":
        st.info(f"`{selected}`: {comp['reason']}")
        return
    st.markdown(f"### {comp['record']['species']} {comp['record']['strain']}")
    st.caption(comp["limitation"])
    st.dataframe(pd.DataFrame(comp["rows"]), width="stretch", hide_index=True)
    st.caption(f"Referencia: `{comp['reference_path'] or 'no disponible'}`")
    with st.expander("Resumen del sistema"):
        st.json(comp["system"])
