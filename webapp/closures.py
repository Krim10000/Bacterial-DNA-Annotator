"""Cierres de brechas para la GUI.

Une el veredicto final (dominio × homología) con su contexto de soporte y lo
deja listo para la pestaña "Cierres".

Vocabulario que se muestra al usuario (ver `docs/PIPELINE_LOGIC.md`):

    validated   dos métodos independientes coinciden (dominio + homología)
    consensus   solo un método, pero coincide con el catálogo previo
    candidate   solo un método, sin cierre previo
    anomalous   el método nuevo y el previo discrepan de dominio
    unresolved  ninguna evidencia

Si los ficheros no existen (pipeline no ejecutado), se devuelve un DataFrame
vacío y la GUI muestra el aviso correspondiente; nunca rompe la página.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

FINAL = Path("metabolism") / "closure_final_top1000.tsv"
VALID = Path("metabolism") / "closure_validation_top1000.tsv"
ORDER = ["validated", "consensus", "candidate", "anomalous", "unresolved"]
LABEL = {
    "validated": "✅ validated — dos métodos independientes",
    "consensus": "🟢 consensus — un método, coincide con el previo",
    "candidate": "🟡 candidate — un método, sin cierre previo",
    "anomalous": "🔴 anomalous — los métodos discrepan",
    "unresolved": "⚪ unresolved — sin evidencia",
}
SHOW = ["gap", "domain", "domain_verdict", "homology_acc", "homology_pident",
        "homology_qcov", "same_protein", "verdict", "status"]


def available(data_root: Path) -> bool:
    return (Path(data_root) / FINAL).exists()


def load_closures(data_root: Path) -> pd.DataFrame:
    """Veredicto final + contexto (rol, soporte, cierre previo)."""
    path = Path(data_root) / FINAL
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    valid = Path(data_root) / VALID
    if valid.exists():
        v = pd.read_csv(valid, sep="\t", dtype=str).fillna("")
        keep = [c for c in ("genome", "gap", "role", "support_400", "prior")
                if c in v.columns]
        df = df.merge(v[keep].drop_duplicates(["genome", "gap"]),
                      on=["genome", "gap"], how="left")
    if "support_400" in df.columns:
        df["support_400"] = pd.to_numeric(df["support_400"], errors="coerce")
    return df


def for_genome(df: pd.DataFrame, genome: str) -> pd.DataFrame:
    if df.empty:
        return df
    return df[df["genome"] == genome].copy()


def summary(sub: pd.DataFrame) -> dict[str, int]:
    if sub.empty:
        return {k: 0 for k in ORDER}
    counts = sub["status"].value_counts().to_dict()
    return {k: int(counts.get(k, 0)) for k in ORDER}


def table(sub: pd.DataFrame) -> pd.DataFrame:
    """Columnas para mostrar, ordenadas por estado y soporte."""
    if sub.empty:
        return sub
    cols = [c for c in SHOW if c in sub.columns]
    key = {k: i for i, k in enumerate(ORDER)}
    tmp = sub.copy()
    tmp["_orden"] = tmp["status"].map(key).fillna(99) if "status" in tmp else 0
    tmp["_soporte"] = pd.to_numeric(tmp.get("support_400"), errors="coerce").fillna(-1)
    tmp = tmp.sort_values(["_orden", "_soporte"], ascending=[True, False])
    return tmp[cols].reset_index(drop=True)
