"""Componentes de operón para la GUI (RegulonDB, E. coli K-12).

Un operón no es solo una lista de genes: tiene promotor, TSS, cajas −10/−35,
factor sigma, reguladores, terminador y evidencia documental. Este módulo carga
la tabla por operón y el detalle por unidad de transcripción.

Fuentes (generadas por `scripts/validation/operon_components.py`):

    data/regulondb/gt/gt_operon_components.tsv   una fila por operón
    data/regulondb/gt/gt_operon_tu_detail.tsv    una fila por TU con sus piezas

Cobertura real del catálogo: 54% de los operones tienen las cinco piezas
regulatorias y solo 12% tienen terminador; por eso la GUI muestra la
completitud en vez de asumir que la pieza existe.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

COMPONENTS = Path("regulondb") / "gt" / "gt_operon_components.tsv"
DETAIL = Path("regulondb") / "gt" / "gt_operon_tu_detail.tsv"

PIECES = [
    ("has_promoter", "Promotor"),
    ("has_tss", "TSS"),
    ("has_minus10", "−10"),
    ("has_minus35", "−35"),
    ("has_sigma", "Sigma"),
    ("has_tf", "TF"),
    ("has_terminator", "Terminador"),
]
SHOW = ["operon_name", "n_genes", "n_tus", "confidence", "n_components",
        "complete", "n_pmids"]


def available(data_root: Path) -> bool:
    return (Path(data_root) / COMPONENTS).exists()


def load_components(data_root: Path) -> pd.DataFrame:
    path = Path(data_root) / COMPONENTS
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for key, _ in PIECES:
        if key in df.columns:
            df[f"pk_{key}"] = df[key].map({"1": "✓", "0": "—"})
    df["n_genes_i"] = pd.to_numeric(df["n_genes"], errors="coerce").fillna(0).astype(int)
    df["n_components_i"] = pd.to_numeric(df["n_components"], errors="coerce").fillna(0).astype(int)
    df["n_pmids_i"] = pd.to_numeric(df["n_pmids"], errors="coerce").fillna(0).astype(int)
    return df


def load_detail(data_root: Path) -> pd.DataFrame:
    path = Path(data_root) / DETAIL
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def detail_for(detail: pd.DataFrame, operon_id: str) -> pd.DataFrame:
    if detail.empty:
        return detail
    return detail[detail["operon_id"] == operon_id].copy()


def piece_table(df: pd.DataFrame) -> pd.DataFrame:
    """Tabla legible: una columna con ✓/— por pieza."""
    cols = [c for c in SHOW if c in df.columns]
    out = df[cols].copy()
    for key, label in PIECES:
        col = f"pk_{key}"
        if col in df.columns:
            out[label] = df[col].values
    out = out.rename(columns={"operon_name": "operón", "n_genes": "genes",
                              "n_tus": "TUs", "confidence": "confianza",
                              "n_components": "piezas", "complete": "completo",
                              "n_pmids": "PMIDs"})
    return out.reset_index(drop=True)


def coverage(df: pd.DataFrame) -> dict[str, int]:
    if df.empty:
        return {}
    return {label: int((df[key] == "1").sum()) for key, label in PIECES}


# ------------------------------------------------------------------ horquillas
HAIRPINS = Path("operon_finder")


def hairpin_genomes(data_root: Path) -> list[str]:
    """Genomas con anotación de horquillas disponible."""
    root = Path(data_root) / HAIRPINS
    if not root.exists():
        return []
    return sorted(p.name.replace("_hairpins.tsv", "")
                  for p in root.glob("*_hairpins.tsv"))


def load_hairpins(data_root: Path, genome: str) -> pd.DataFrame:
    path = Path(data_root) / HAIRPINS / f"{genome}_hairpins.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df["dG_i"] = pd.to_numeric(df.get("dG"), errors="coerce")
    df["stem_i"] = pd.to_numeric(df.get("stem"), errors="coerce")
    return df


HAIRPIN_COLS = ["contig", "start", "end", "class", "stem", "gc", "loop",
                "polyU", "dG", "gap", "left_gene", "right_gene"]


def hairpin_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    cols = [c for c in HAIRPIN_COLS if c in df.columns]
    order = {"A_fuerte": 0, "B_media": 1, "C_debil": 2}
    out = df.copy()
    out["_o"] = out["class"].map(order).fillna(3)
    out = out.sort_values(["_o", "dG_i"], ascending=[True, False])
    return out[cols].reset_index(drop=True)


def hairpin_summary(df: pd.DataFrame) -> dict[str, int]:
    if df.empty:
        return {}
    return df["class"].value_counts().to_dict()


# ------------------------------------------------- operones por reglas (sin ML)
#
# Tres carpetas, en este orden:
#   operon_finder/            K-12 y las especies de referencia
#   operon_finder/rules_all/  los 400 genomas libres del proyecto
#   operon_finder/excluded/   los 90 parásitos/simbiontes, fuera de las métricas
RULE_DIRS = [Path("operon_finder"), Path("operon_finder") / "rules_all",
             Path("operon_finder") / "excluded"]
MARCA_EXCLUIDO = " (excluido)"


def rule_genomes(data_root: Path) -> list[str]:
    """Todos los genomas con operones por reglas, excluidos al final."""
    out = []
    for d in RULE_DIRS:
        root = Path(data_root) / d
        if not root.exists():
            continue
        for p in root.glob("*_operons_rules.tsv"):
            g = p.name.replace("_operons_rules.tsv", "")
            if d.name == "excluded":
                g += MARCA_EXCLUIDO
            out.append(g)
    return sorted(out, key=lambda x: (x.endswith(MARCA_EXCLUIDO), x))


def _rule_path(data_root: Path, genome: str, suffix: str) -> Path | None:
    g = genome.replace(MARCA_EXCLUIDO, "")
    for d in RULE_DIRS:
        path = Path(data_root) / d / f"{g}{suffix}"
        if path.exists():
            return path
    return None


def load_rule_operons(data_root: Path, genome: str) -> pd.DataFrame:
    path = _rule_path(data_root, genome, "_operons_rules.tsv")
    if path is None:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def load_rule_boundaries(data_root: Path, genome: str) -> pd.DataFrame:
    path = _rule_path(data_root, genome, "_boundary_rules.tsv")
    if path is None:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def rule_summary(boundaries: pd.DataFrame) -> dict[str, int]:
    if boundaries.empty:
        return {}
    return boundaries["rule"].value_counts().to_dict()


RULE_LABEL = {
    "R1_atenuador": "R1 · péptido líder → unir",
    "R1b_atenuador_funcional": "R1b · horquilla + STRING → unir",
    "R2_terminador": "R2 · terminador fuerte → cortar",
    "R3_promotor": "R3 · promotor σ70 → cortar",
    "R4_distancia": "R4 · distancia → cortar",
    "R5_defecto": "R5 · sin señal → unir",
}


# ------------------------------------------- sospecha por cadena metabólica
SUSPECTS = Path("metabolism") / "operon_suspects.tsv"


def load_suspects(data_root: Path, genome: str) -> pd.DataFrame:
    """Grupos de genes contiguos cuyas reacciones comparten metabolitos."""
    path = Path(data_root) / SUSPECTS
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    sub = df[df["genome"] == genome].copy()
    for col in ("start", "end", "n_genes", "gap_mean"):
        if col in sub.columns:
            sub[col] = pd.to_numeric(sub[col], errors="coerce")
    return sub.sort_values(["contig", "start"]).reset_index(drop=True)
