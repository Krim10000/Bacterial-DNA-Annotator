#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Capa de contexto determinista y auditable para la asignación de sigma.

No entrena un modelo opaco. Cada evidencia contextual:

1. se extrae de una fuente local declarada (genoma, anotación, tabla curada);
2. se discretiza en categorías reproducibles;
3. se convierte en un *log-odds* (log-likelihood ratio) con suavizado
   Dirichlet, estimado **solo con las particiones de entrenamiento**;
4. se combina con el z del PWM mediante una suma lineal transparente;
5. se puede apagar (peso 0) para medir su contribución por ablación.

La ausencia de dato (`none`) es **neutral**: no aporta log-odds. No se usa la
falta de anotación como prior encubierto. La etiqueta no se fuerza: con
`scale=0` la capa reproduce exactamente el ranking del PWM, y el modo acotado
(`gap`) impide que el contexto desplace a un candidato que el motivo ya lidera
por un margen amplio.

Evidencias implementadas (ver `EVIDENCE_SPECS`):

  operon_position   ubicación del TSS respecto al gen/operón aguas abajo
  leader            distancia al gen aguas abajo (leader)
  operon_size       número de genes del operón que contiene/recibe el TSS
  pathway_ec        clase EC del gen aguas abajo (solo si `genes.tsv` trae EC)
  domain            familia/dominio del gen aguas abajo (solo si hay dato local)
  conservation      amplitud de conservación del gen aguas abajo (sintenía)
  expression        expresión medida del gen aguas abajo (solo si hay tabla)
  sigma54_geometry  geometría -24/-12 (o ventana proximal) del hit sigma54

Las fuentes que no existen para un genoma se declaran `available=False` y
cobertura 0; no se rellenan con datos de otra especie.
"""
from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

VERSION = "1.0.0"

SIGMA_ORDER = ("sigma70", "sigma54", "sigma32", "sigma38", "sigma28", "sigma24")


# --------------------------------------------------------------------------
# Registro de evidencias: fuente, peso declarado y ámbito
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class EvidenceSpec:
    """Metadatos de una evidencia contextual.

    `weight` es el peso declarado *a priori*. No se ajusta con el gold de
    evaluación; la ablación lo pone a 0 para medir su aporte. `source` debe
    apuntar a un fichero local; si no existe, la evidencia queda no disponible.
    """

    name: str
    source: str
    description: str
    weight: float
    applies_to: tuple[str, ...] = SIGMA_ORDER
    kind: str = "context"


EVIDENCE_SPECS: dict[str, EvidenceSpec] = {
    "operon_position": EvidenceSpec(
        name="operon_position",
        source="data/regulation/<genoma>/genes.tsv; "
               "subtilis: data/dbtbs/subtilis_tss_tts.tsv",
        description="TSS en borde de gen/operón, interno, antisentido o intergénico.",
        weight=1.0,
    ),
    "leader": EvidenceSpec(
        name="leader",
        source="data/regulation/<genoma>/genes.tsv; "
               "subtilis: data/dbtbs/subtilis_tss_tts.tsv:dist_gene",
        description="Distancia del TSS al inicio del gen aguas abajo, en bins.",
        weight=0.5,
    ),
    "operon_size": EvidenceSpec(
        name="operon_size",
        source="data/regulation/<genoma>/operons.tsv; "
               "subtilis: data/dbtbs/gt_operons.tsv",
        description="Número de genes del operón que contiene o recibe el TSS.",
        weight=0.5,
    ),
    "pathway_ec": EvidenceSpec(
        name="pathway_ec",
        source="data/regulation/<genoma>/genes.tsv:ec",
        description="Primer dígito de la clase EC del gen aguas abajo.",
        weight=1.0,
    ),
    "domain": EvidenceSpec(
        name="domain",
        source="data/sigma/structures/structure_features.tsv (no cubre K-12/BSu); "
               "data/regulation/<genoma>/genes.tsv:domain (no existe)",
        description="Dominio/estructura del gen aguas abajo. No se usa `fam`: es un "
                    "código casi único por gen (4212 valores/4262 genes en K-12) y "
                    "sobreajustaría.",
        weight=0.5,
    ),
    "conservation": EvidenceSpec(
        name="conservation",
        source="data/regulation/synteny_pairs.tsv + genes.tsv:acc; "
               "subtilis: data/dbtbs/bsu_vs_k12.tsv",
        description="Amplitud de conservación del gen aguas abajo (bins).",
        weight=1.0,
    ),
    "expression": EvidenceSpec(
        name="expression",
        source="data/expression/<genoma>.tsv (opcional)",
        description="Gen aguas abajo con expresión medida > 0.",
        weight=1.0,
    ),
    "sigma54_geometry": EvidenceSpec(
        name="sigma54_geometry",
        source="PWM sigma54 + ventana evaluada",
        description="Geometría del hit sigma54: canónica o desplazada.",
        weight=1.0,
        applies_to=("sigma54",),
        kind="geometry",
    ),
}

#: Evidencias que se usan por defecto en la capa (todas las disponibles).
DEFAULT_EVIDENCE_WEIGHTS: dict[str, float] = {
    name: spec.weight for name, spec in EVIDENCE_SPECS.items()
}


# --------------------------------------------------------------------------
# E/S
# --------------------------------------------------------------------------
def _read_tsv(path: Path) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _as_int(value: Any, default: int | None = None) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def _strand(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"+", "1", "1.0", "forward", "fwd", "plus"}:
        return "+"
    if text in {"-", "-1", "reverse", "rev", "minus"}:
        return "-"
    return ""


def _norm_strand(value: Any) -> int:
    return 1 if _strand(value) == "+" else -1


def _contig_of(row: Mapping[str, Any], default: str = "") -> str:
    for key in ("contig", "seqid", "chromosome", "replicon"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return default


def _bin_leader(value: int | None) -> str:
    if value is None or value < 0:
        return "none"
    if value <= 50:
        return "0-50"
    if value <= 150:
        return "50-150"
    if value <= 500:
        return "150-500"
    return ">500"


def _bin_operon_size(value: int | None) -> str:
    if value is None or value <= 0:
        return "none"
    if value == 1:
        return "1"
    if value == 2:
        return "2"
    if value <= 5:
        return "3-5"
    return ">5"


def _bin_conservation(value: int | None) -> str:
    if value is None or value <= 0:
        return "none"
    if value < 10:
        return "low"
    if value < 100:
        return "mid"
    return "high"


def _ec_class(value: str) -> str:
    text = str(value or "").strip()
    for part in re.split(r"[;,|/\s]+", text):
        part = part.strip()
        if part and part[0] in "1234567" and "." in part:
            return part[0]
    return "none"


def _default_contig(genes: Sequence[Mapping[str, Any]]) -> str:
    counts = Counter(_contig_of(g) for g in genes if _contig_of(g))
    return counts.most_common(1)[0][0] if counts else ""


# --------------------------------------------------------------------------
# Índice de genes y de operones
# --------------------------------------------------------------------------
class GeneIndex:
    """Índice determinista de genes para mapear un TSS al gen aguas abajo.

    No usa el gold: solo coordenadas, hebra y anotación del genoma.
    """

    def __init__(self, genes: Iterable[Mapping[str, Any]]):
        self.by: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
        for raw in genes:
            row = dict(raw)
            start, end = _as_int(row.get("start")), _as_int(row.get("end"))
            if start is None or end is None or start > end:
                continue
            strand = _norm_strand(row.get("strand"))
            row["_start"], row["_end"], row["_strand"] = start, end, strand
            self.by[(_contig_of(row), strand)].append(row)
        self._starts: dict[tuple[str, int], list[int]] = {}
        self._ends: dict[tuple[str, int], list[int]] = {}
        for key, rows in self.by.items():
            rows.sort(key=lambda r: (r["_start"], r["_end"]))
            self._starts[key] = [r["_start"] for r in rows]
            self._ends[key] = [r["_end"] for r in rows]

    def inside(self, contig: str, pos: int, strand: str) -> dict[str, Any] | None:
        key = (contig, _norm_strand(strand))
        rows = self.by.get(key, [])
        starts = self._starts.get(key, [])
        if not rows:
            return None
        i = int(np.searchsorted(starts, pos, side="right"))
        # Los genes pueden solaparse; se revisan unos pocos hacia atrás.
        for row in rows[max(0, i - 8):i + 1]:
            if row["_start"] <= pos <= row["_end"]:
                return row
        return None

    def downstream(self, contig: str, pos: int, strand: str) -> dict[str, Any] | None:
        key = (contig, _norm_strand(strand))
        rows = self.by.get(key, [])
        if not rows:
            return None
        if _norm_strand(strand) == 1:
            starts = self._starts[key]
            i = int(np.searchsorted(starts, pos, side="left"))
            return rows[i] if i < len(rows) else None
        ends = self._ends[key]
        i = int(np.searchsorted(ends, pos, side="right")) - 1
        return rows[i] if 0 <= i < len(rows) else None


class OperonIndex:
    """Índice de operones por coordenada (no por ID, para no depender de DBs)."""

    def __init__(self, rows: Iterable[Mapping[str, Any]], default_contig: str = ""):
        self.by: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
        for raw in rows:
            row = dict(raw)
            start, end = _as_int(row.get("start")), _as_int(row.get("end"))
            if start is None or end is None or start > end:
                continue
            strand = _norm_strand(row.get("strand"))
            size = _as_int(row.get("n_genes") or row.get("n_genes_tus"))
            row["_start"], row["_end"], row["_n"] = start, end, size
            self.by[(_contig_of(row, default_contig), strand)].append(row)
        for key in self.by:
            self.by[key].sort(key=lambda r: (r["_start"], r["_end"]))
        self._starts = {k: [r["_start"] for r in v] for k, v in self.by.items()}
        self._ends = {k: [r["_end"] for r in v] for k, v in self.by.items()}

    def containing_or_downstream(self, contig: str, pos: int,
                                 strand: str) -> dict[str, Any] | None:
        key = (contig, _norm_strand(strand))
        rows = self.by.get(key, [])
        if not rows:
            return None
        # Primero: operón que contiene el TSS.
        starts = self._starts[key]
        i = int(np.searchsorted(starts, pos, side="right"))
        for row in rows[max(0, i - 4):i + 1]:
            if row["_start"] <= pos <= row["_end"]:
                return row
        # Si no, el operón aguas abajo más cercano.
        if _norm_strand(strand) == 1:
            j = int(np.searchsorted(starts, pos, side="left"))
            return rows[j] if j < len(rows) else None
        ends = self._ends[key]
        j = int(np.searchsorted(ends, pos, side="right")) - 1
        return rows[j] if 0 <= j < len(rows) else None


# --------------------------------------------------------------------------
# Fuentes de conservación y expresión
# --------------------------------------------------------------------------
def _conservation_breadth(data_root: Path) -> dict[str, int]:
    """acc -> máximo n_genomes observado en `synteny_pairs.tsv`."""
    out: dict[str, int] = {}
    for row in _read_tsv(data_root / "regulation" / "synteny_pairs.tsv"):
        n = _as_int(row.get("n_genomes"), 0) or 0
        for acc in (row.get("acc_a") or "", row.get("acc_b") or ""):
            acc = acc.strip()
            if acc and n > out.get(acc, 0):
                out[acc] = n
    return out


def _bsu_to_acc(data_root: Path, genes_k12: Sequence[Mapping[str, Any]]) -> dict[str, str]:
    """BSU_* -> acc de K-12 vía bsu_vs_k12.tsv y el locus_tag b#### de K-12."""
    b_to_acc: dict[str, str] = {}
    for row in genes_k12:
        locus = str(row.get("locus_tag") or "").strip().lower()
        acc = str(row.get("acc") or "").strip()
        if locus and acc:
            b_to_acc[locus] = acc
    out: dict[str, str] = {}
    for row in _read_tsv(data_root / "dbtbs" / "bsu_vs_k12.tsv"):
        # Formato: BSU <TAB> b#### <TAB> identidad <TAB> evalue <TAB> longitud
        fields = [str(x).strip() for x in row.values() if str(x).strip()]
        if len(fields) >= 2:
            bsu, bnum = fields[0], fields[1].lower()
            if bsu.startswith("BSU_") and bnum in b_to_acc:
                out[bsu] = b_to_acc[bnum]
    return out


def _name_to_bsu(data_root: Path) -> dict[str, str]:
    """Nombre de gen de B. subtilis -> locus BSU_* (tabla curada DBTBS)."""
    out: dict[str, str] = {}
    for row in _read_tsv(data_root / "dbtbs" / "subtilis_tss_tts.tsv"):
        gene = str(row.get("gene") or "").strip()
        feature = str(row.get("feature") or "").strip()
        if gene and feature.startswith("BSU_"):
            out.setdefault(gene.lower(), feature)
    for row in _read_tsv(data_root / "dbtbs" / "gt_operons.tsv"):
        genes = [g.strip() for g in str(row.get("genes") or "").split(",")]
        names = [g.strip() for g in str(row.get("gene_names") or "").split(",")]
        for bsu, name in zip(genes, names):
            if bsu.startswith("BSU_") and name:
                out.setdefault(name.lower(), bsu)
    return out


def load_expression_map(data_root: Path, genome: str) -> dict[str, float]:
    """Carga expresión génica opcional; NO lee reglas de regulación.

    Formatos aceptados: TSV con `gene_id`/`locus_tag`/`gene` y
    `tpm`/`expression`/`value`. Si no existe, devuelve {}.
    """
    candidates = [
        data_root / "expression" / f"{genome}.tsv",
        data_root / "expression" / f"{genome}_expression.tsv",
        data_root / "regulation" / genome / "expression.tsv",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        return {}
    rows = _read_tsv(path)
    if not rows:
        return {}
    key = next((k for k in ("gene_id", "locus_tag", "gene", "id") if k in rows[0]), None)
    val = next((k for k in ("tpm", "expression", "value", "counts") if k in rows[0]), None)
    if not key or not val:
        return {}
    out: dict[str, float] = {}
    for row in rows:
        ident = str(row.get(key) or "").strip()
        if not ident:
            continue
        try:
            out[ident] = float(row.get(val) or 0.0)
        except (TypeError, ValueError):
            continue
    return out


# --------------------------------------------------------------------------
# Construcción de la tabla de contexto
# --------------------------------------------------------------------------
@dataclass
class ContextTable:
    genome: str
    sites: dict[tuple[str, int, str], dict[str, str]] = field(default_factory=dict)
    coverage: dict[str, float] = field(default_factory=dict)
    sources: dict[str, str] = field(default_factory=dict)
    available: dict[str, bool] = field(default_factory=dict)

    def row(self, site: tuple[str, int, str]) -> dict[str, str]:
        return self.sites.get(site, {})


def _curated_tss_rows(data_root: Path, genome: str) -> dict[tuple[str, int, str], dict[str, str]]:
    """Tabla curada TSS de B. subtilis (location/dist_gene), si existe."""
    if genome != "subtilis":
        return {}
    path = data_root / "dbtbs" / "subtilis_tss_tts.tsv"
    out: dict[tuple[str, int, str], dict[str, str]] = {}
    for row in _read_tsv(path):
        if (row.get("type") or "").strip().upper() != "TSS":
            continue
        pos = _as_int(row.get("pos"))
        strand = _strand(row.get("strand"))
        if pos is None or not strand:
            continue
        out[(_contig_of(row), pos, strand)] = row
    return out


def _operon_rows(data_root: Path, genome: str) -> tuple[list[dict[str, str]], str]:
    base = data_root / "regulation" / genome
    rows = _read_tsv(base / "operons.tsv")
    if rows:
        return rows, ""
    if genome == "subtilis":
        for candidate, default in (
            (data_root / "dbtbs" / "gt_operons.tsv", "basu168"),
            (data_root / "dbtbs" / "gt_operons_dbtbs.tsv", ""),
            (data_root / "operon_finder" / "subtilis_fixed" / "subtilis_operons_rules.tsv", ""),
        ):
            rows = _read_tsv(candidate)
            if rows:
                return rows, default
    return [], ""


def build_context(data_root: Path, genome: str,
                  sites: Iterable[tuple[str, int, str]],
                  use_curated: bool = False) -> ContextTable:
    """Extrae la tabla de contexto para los TSS indicados.

    `sites` son tuplas `(contig, pos, strand)` con strand `+`/`-`.
    Por defecto se mapea el TSS al gen aguas abajo con `genes.tsv`
    (método uniforme para cualquier genoma). En *B. subtilis* se puede
    pedir `use_curated=True` para usar `subtilis_tss_tts.tsv` (DBTBS), pero
    esa tabla comparte procedencia con el gold; se reporta como sensibilidad,
    no como fuente principal.

    Todas las fuentes son locales; si faltan, la categoría es `none` y la
    evidencia se declara no disponible.
    """
    data_root = Path(data_root)
    site_list = list(dict.fromkeys(sites))
    genes = _read_tsv(data_root / "regulation" / genome / "genes.tsv")
    default_contig = _default_contig(genes)
    gene_index = GeneIndex(genes)
    curated = _curated_tss_rows(data_root, genome) if use_curated else {}

    operon_rows, operon_contig = _operon_rows(data_root, genome)
    operon_index = OperonIndex(operon_rows, operon_contig or default_contig)

    cons_map = _conservation_breadth(data_root)
    genes_k12 = genes if genome == "ecoli_k12" else _read_tsv(
        data_root / "regulation" / "ecoli_k12" / "genes.tsv")
    bsu_to_acc = _bsu_to_acc(data_root, genes_k12) if genome == "subtilis" else {}
    name_to_bsu = _name_to_bsu(data_root) if genome == "subtilis" else {}
    expression = load_expression_map(data_root, genome)

    table = ContextTable(genome=genome)
    table.sources = {
        spec_name: spec.source for spec_name, spec in EVIDENCE_SPECS.items()
    }
    has_ec = bool(genes) and "ec" in (genes[0] if genes else {})
    # `genes.tsv:fam` es un código de familia casi único por gen (4212 valores
    # para 4262 genes en K-12): no es un dominio reutilizable y sobreajustaría.
    # Solo se acepta estructura/dominio si hay una tabla real que cubra el genoma.
    has_struct = (data_root / "sigma" / "structures" / "structure_features.tsv").exists()
    if has_struct:
        # La tabla de estructura cubre genomas GCF, no necesariamente estos dos.
        genomes_struct = {
            str(r.get("genome") or "").strip()
            for r in _read_tsv(data_root / "sigma" / "structures" / "structure_features.tsv")
        }
        has_struct = genome in genomes_struct
    table.available = {
        "operon_position": bool(curated) or bool(genes),
        "leader": bool(curated) or bool(genes),
        "operon_size": bool(operon_rows),
        "pathway_ec": has_ec,
        "domain": has_struct,
        "conservation": bool(cons_map) and (genome == "ecoli_k12" or bool(bsu_to_acc)),
        "expression": bool(expression),
        "sigma54_geometry": True,  # se completa en la evaluación con el hit PWM
    }

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for site in site_list:
        contig, pos, strand = site
        row: dict[str, str] = {}
        curated_row = curated.get(site)
        gene: dict[str, Any] | None = None

        # --- posición y leader ---
        if curated_row is not None:
            loc = (curated_row.get("location") or "").strip() or "none"
            row["operon_position"] = loc if loc in {
                "intergenic", "gene_edge", "inside_gene", "antisense"} else "none"
            leader = _as_int(curated_row.get("dist_gene"))
            row["leader"] = _bin_leader(leader)
            # Gen aguas abajo: el `feature` curado si es un locus BSU.
            feature = (curated_row.get("feature") or "").strip()
            if feature.startswith("BSU_"):
                gene = {"gene_name": (curated_row.get("gene") or "").strip(),
                        "gene_id": feature, "_start": 0, "_end": 0}
        else:
            inside = gene_index.inside(contig, pos, strand)
            anti = gene_index.inside(contig, pos, "-" if strand == "+" else "+")
            downstream = gene_index.downstream(contig, pos, strand)
            if inside is not None:
                row["operon_position"] = "inside_gene"
            elif anti is not None:
                row["operon_position"] = "antisense"
            elif downstream is None:
                row["operon_position"] = "no_downstream"
            else:
                leader = (downstream["_start"] - pos if strand == "+"
                          else pos - downstream["_end"])
                row["operon_position"] = "gene_edge" if leader <= 100 else "intergenic"
            if downstream is None:
                row["leader"] = "none"
            else:
                leader = (downstream["_start"] - pos if strand == "+"
                          else pos - downstream["_end"])
                row["leader"] = _bin_leader(leader)
                gene = downstream

        # Fallback: si la tabla curada no dio gen, usar el índice genómico.
        if gene is None and row.get("operon_position") not in {"inside_gene", "antisense"}:
            gene = gene_index.downstream(contig, pos, strand)

        # --- tamaño de operón ---
        operon = operon_index.containing_or_downstream(contig, pos, strand)
        row["operon_size"] = _bin_operon_size(
            _as_int(operon.get("n_genes") or operon.get("n_genes_tus")) if operon else None)

        # --- EC, dominio, conservación, expresión ---
        row["pathway_ec"] = "none"
        row["domain"] = "none"
        row["conservation"] = "none"
        row["expression"] = "none"
        if gene is not None:
            gene_name = str(gene.get("gene_name") or "").strip()
            gene_id = str(gene.get("gene_id") or "").strip()
            row["pathway_ec"] = _ec_class(str(gene.get("ec") or ""))
            # conservación
            acc = str(gene.get("acc") or "").strip()
            if not acc and genome == "subtilis":
                bsu = name_to_bsu.get(gene_name.lower()) or (
                    gene_id if gene_id.startswith("BSU_") else "")
                if bsu:
                    acc = bsu_to_acc.get(bsu, "")
            if acc:
                row["conservation"] = _bin_conservation(cons_map.get(acc))
            # expresión
            if expression:
                value = None
                for key in (gene_id, gene_name):
                    if key and key in expression:
                        value = expression[key]
                        break
                if value is not None:
                    row["expression"] = "on" if value > 0 else "off"

        row.setdefault("sigma54_geometry", "none")
        for name in EVIDENCE_SPECS:
            row.setdefault(name, "none")
            if row[name] != "none":
                counts[name][row[name]] += 1
        table.sites[site] = row

    for name in EVIDENCE_SPECS:
        table.coverage[name] = (len([s for s in site_list
                                     if table.sites[s].get(name, "none") != "none"])
                                / len(site_list)) if site_list else 0.0
        if not table.available.get(name, False):
            table.coverage[name] = 0.0
    table.sources["_mode"] = (
        "curated:subtilis_tss_tts.tsv" if use_curated and curated
        else "genomic:genes.tsv")
    return table


def active_weights(context: ContextTable, min_coverage: float = 0.05,
                   include_geometry: bool = True) -> dict[str, float]:
    """Pesos declarados que se activan para un genoma.

    Una evidencia sin fuente local o con cobertura inferior a `min_coverage`
    recibe peso 0.0: se conserva en la salida, pero no mueve el ranking.
    `sigma54_geometry` se completa en la evaluación; si se incluye, arranca
    con su peso declarado.
    """
    out: dict[str, float] = {}
    for name, spec in EVIDENCE_SPECS.items():
        if name == "sigma54_geometry" and not include_geometry:
            out[name] = 0.0
            continue
        coverage = float(context.coverage.get(name, 0.0))
        available = bool(context.available.get(name, False))
        out[name] = float(spec.weight) if available and coverage >= min_coverage else 0.0
    return out


# --------------------------------------------------------------------------
# Modelo de log-odds
# --------------------------------------------------------------------------
@dataclass
class LogOddsModel:
    """Tablas de log-likelihood ratio por evidencia y sigma.

    `tables[e][s][v] = log( P(v | sigma=s) / P(v) )` con suavizado Dirichlet.
    La categoría `none` no se modela: en scoring aporta 0.0.
    """

    tables: dict[str, dict[str, dict[str, float]]] = field(default_factory=dict)
    counts: dict[str, Any] = field(default_factory=dict)
    support: dict[str, dict[str, int]] = field(default_factory=dict)
    alpha: float = 1.0
    min_support: int = 5
    n_train: int = 0
    evidence_used: tuple[str, ...] = ()

    def llr(self, site_row: Mapping[str, str], sigma: str,
            weights: Mapping[str, float] | None = None) -> float:
        total = 0.0
        # `weights=None` usa los pesos declarados; un dict explícito es
        # EXCLUYENTE (las evidencias ausentes pesan 0). Así la ablación mide
        # una evidencia y no vuelve a sumar las demás por defecto.
        if weights is None:
            weights = {name: spec.weight for name, spec in EVIDENCE_SPECS.items()}
        for evidence, table in self.tables.items():
            weight = float(weights.get(evidence, 0.0))
            if weight == 0.0:
                continue
            spec = EVIDENCE_SPECS[evidence]
            if sigma not in spec.applies_to:
                continue
            value = str(site_row.get(evidence, "none"))
            if not value or value == "none":
                continue
            # Categorías con muy poco soporte en entrenamiento son neutrales:
            # evita que una familia/valor raro aprendido de 1-2 sitios mueva
            # el ranking (sobreajuste).
            if self.support.get(evidence, {}).get(value, 0) < self.min_support:
                continue
            total += weight * float(table.get(sigma, {}).get(value, 0.0))
        return total

    def to_rows(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for evidence, per_sigma in self.tables.items():
            spec = EVIDENCE_SPECS[evidence]
            for sigma, per_value in per_sigma.items():
                for value, llr in per_value.items():
                    rows.append({
                        "evidence": evidence, "sigma": sigma, "value": value,
                        "log_odds": round(float(llr), 6),
                        "n_train": int(self.support.get(evidence, {}).get(value, 0)),
                        "min_support": int(self.min_support),
                        "weight": spec.weight, "source": spec.source,
                    })
        return rows


def fit_logodds(context: ContextTable, labels: Mapping[tuple[str, int, str], Iterable[str]],
                evidence_names: Iterable[str] | None = None,
                alpha: float = 1.0, min_support: int = 5) -> LogOddsModel:
    """Ajusta las tablas de LLR usando SOLO los sitios de entrenamiento.

    `labels[site]` es el conjunto de sigma verdaderos. Un sitio multietiqueta
    cuenta para cada sigma verdadero en el numerador; en el denominador de
    fondo cuenta una vez.
    """
    names = tuple(evidence_names or [n for n, s in EVIDENCE_SPECS.items()
                                     if s.kind == "context"])
    model = LogOddsModel(alpha=float(alpha), min_support=int(min_support))
    model.n_train = len(labels)
    model.evidence_used = names
    for evidence in names:
        spec = EVIDENCE_SPECS.get(evidence)
        if spec is None:
            continue
        value_counts: Counter[str] = Counter()
        per_sigma: dict[str, Counter[str]] = defaultdict(Counter)
        n_sigma: Counter[str] = Counter()
        for site, sigmas in labels.items():
            row = context.row(site)
            value = str(row.get(evidence, "none"))
            if not value or value == "none":
                continue
            # En el fondo, un sitio cuenta una vez por valor observado.
            unique_values = {value}
            for v in unique_values:
                value_counts[v] += 1
            present = [s for s in sigmas if s in SIGMA_ORDER]
            for sigma in present:
                per_sigma[sigma][value] += 1
            for sigma in set(present):
                n_sigma[sigma] += 1
        vocab = sorted(value_counts)
        if not vocab:
            continue
        n_all = sum(value_counts.values())
        per_sigma_tables: dict[str, dict[str, float]] = {}
        for sigma in SIGMA_ORDER:
            n_s = n_sigma.get(sigma, 0)
            if n_s == 0:
                per_sigma_tables[sigma] = {}
                continue
            table: dict[str, float] = {}
            for value in vocab:
                p_s = ((per_sigma[sigma].get(value, 0) + alpha)
                       / (n_s + alpha * (len(vocab) + 1)))
                p_bg = ((value_counts.get(value, 0) + alpha)
                        / (n_all + alpha * (len(vocab) + 1)))
                table[value] = math.log(p_s / p_bg) if p_s > 0 and p_bg > 0 else 0.0
            per_sigma_tables[sigma] = table
        model.tables[evidence] = per_sigma_tables
        model.support[evidence] = dict(value_counts)
        model.counts[evidence] = {
            "vocab": vocab,
            "n_all": n_all,
            "n_sigma": dict(n_sigma),
            "min_support": int(min_support),
        }
    return model


def combine(motif_z: Mapping[str, float], model: LogOddsModel,
            site_row: Mapping[str, str], weights: Mapping[str, float] | None = None,
            scale: float = 0.5) -> dict[str, float]:
    """Score final = z del motivo + escala * suma ponderada de LLR.

    Con `scale=0` el resultado es exactamente el ranking del motivo.
    """
    out: dict[str, float] = {}
    for sigma, z in motif_z.items():
        out[sigma] = float(z) + float(scale) * model.llr(site_row, sigma, weights)
    return out


def select_top1(final: Mapping[str, float], motif_z: Mapping[str, float],
                gap: float | None = 1.0) -> str:
    """Top-1 con guardarraíl: el contexto solo puede promover candidatos
    dentro de `gap` unidades z del mejor motivo. `gap=None` desactiva el
    guardarraíl; `gap=0` reproduce el top-1 del motivo.
    """
    if not final:
        return ""
    if not motif_z:
        return max(final, key=final.get)
    zbest = max(motif_z.values())
    allowed = [s for s in final if gap is None or zbest - float(motif_z.get(s, -1e9)) <= gap]
    if not allowed:
        allowed = list(final)
    return max(allowed, key=lambda s: final[s])


def rank_topk(final: Mapping[str, float], k: int = 3) -> list[str]:
    return sorted(final, key=lambda s: (-float(final[s]), s))[:k]


# --------------------------------------------------------------------------
# Geometría sigma54
# --------------------------------------------------------------------------
def sigma54_geometry_category(hit: Mapping[str, Any] | None, genome: str) -> str:
    """Clasifica la geometría de un hit sigma54.

    - Dos cajas (K-12): canónica si el espaciador cae en 10-11 nt.
    - Una caja (B. subtilis): canónica si el extremo de la caja cae en la
      ventana proximal -32..-9 (región usada para construir la PWM).
    """
    if not hit:
        return "none"
    arch = str(hit.get("arch") or "")
    spacer = hit.get("spacer")
    if arch != "single-box" and spacer is not None:
        try:
            return "canonical" if 10 <= int(spacer) <= 11 else "shifted"
        except (TypeError, ValueError):
            return "none"
    dist = hit.get("dist_to_start")
    try:
        dist_i = int(dist)
    except (TypeError, ValueError):
        return "none"
    return "canonical" if 9 <= dist_i <= 32 else "shifted"


# --------------------------------------------------------------------------
# Calibración
# --------------------------------------------------------------------------
def softmax(scores: Mapping[str, float], temperature: float = 1.0,
            sigmas: Sequence[str] = SIGMA_ORDER) -> dict[str, float]:
    temp = float(temperature) if temperature and temperature > 0 else 1.0
    values = np.asarray([float(scores.get(s, -1e9)) for s in sigmas], dtype=float)
    values = values / temp
    values -= values.max()
    exp = np.exp(values)
    total = exp.sum()
    if not math.isfinite(total) or total <= 0:
        uniform = 1.0 / len(sigmas)
        return {s: uniform for s in sigmas}
    return {s: float(exp[i] / total) for i, s in enumerate(sigmas)}


def _site_nll(probs: Mapping[str, float], true: Iterable[str]) -> float:
    p_true = sum(float(probs.get(s, 0.0)) for s in true)
    return -math.log(max(p_true, 1e-12))


def _softmax_matrix(scores: np.ndarray, temperature: float) -> np.ndarray:
    temp = float(temperature) if temperature and temperature > 0 else 1.0
    values = np.asarray(scores, dtype=float) / temp
    values = values - values.max(axis=1, keepdims=True)
    exp = np.exp(values)
    total = exp.sum(axis=1, keepdims=True)
    bad = (~np.isfinite(total)) | (total <= 0)
    if bad.any():
        exp[bad.ravel()] = 1.0
        total[bad] = exp.shape[1]
    return exp / total


def _records_to_arrays(records: Sequence[Mapping[str, Any]],
                       score_key: str, true_key: str,
                       sigmas: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    n = len(records)
    scores = np.full((n, len(sigmas)), -1e9, dtype=float)
    truth = np.zeros((n, len(sigmas)), dtype=float)
    for i, record in enumerate(records):
        source = record[score_key]
        true = set(record[true_key])
        for j, sigma in enumerate(sigmas):
            scores[i, j] = float(source.get(sigma, -1e9))
            truth[i, j] = 1.0 if sigma in true else 0.0
    return scores, truth


def fit_temperature(sites: Sequence[Mapping[str, Any]],
                    score_key: str = "scores", true_key: str = "true",
                    grid: Sequence[float] | None = None,
                    sigmas: Sequence[str] = SIGMA_ORDER) -> float:
    """Ajusta una única temperatura por máxima verosimilitud (1 parámetro).

    `sites` debe traer `scores` (dict sigma->score) y `true` (conjunto).
    Es calibración escalar transparente, no un modelo opaco.
    """
    if not sites:
        return 1.0
    scores, truth = _records_to_arrays(sites, score_key, true_key, sigmas)
    candidates = list(grid) if grid is not None else [i / 20 for i in range(2, 121)]
    best_t, best_nll = 1.0, float("inf")
    for temp in candidates:
        probs = _softmax_matrix(scores, temp)
        p_true = (probs * truth).sum(axis=1)
        nll = float(-np.log(np.clip(p_true, 1e-12, None)).mean())
        if nll < best_nll:
            best_nll, best_t = nll, float(temp)
    return best_t


def calibration_metrics(records: Sequence[Mapping[str, Any]],
                        score_key: str = "probs", true_key: str = "true",
                        sigmas: Sequence[str] = SIGMA_ORDER,
                        n_bins: int = 10) -> dict[str, Any]:
    """ECE, Brier y NLL multietiqueta de un conjunto de probabilidades.

    Cada `record` trae `probs` (dict sigma->p) y `true` (conjunto).
    """
    if not records:
        return {"n": 0, "ece": None, "brier": None, "nll": None,
                "mean_confidence": None, "top1_accuracy": None, "bins": []}
    probs, truth = _records_to_arrays(records, score_key, true_key, sigmas)
    # Re-normalizar por si el registro trae probabilidades ligeramente sucias.
    probs = np.clip(probs, 0.0, None)
    total = probs.sum(axis=1, keepdims=True)
    total[total <= 0] = 1.0
    probs = probs / total
    best = probs.argmax(axis=1)
    rows = np.arange(len(probs))
    conf = probs[rows, best]
    correct = truth[rows, best]
    top1 = float(correct.mean())
    brier = float(((probs - truth) ** 2).mean())
    p_true = (probs * truth).sum(axis=1)
    nll = float(-np.log(np.clip(p_true, 1e-12, None)).mean())
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    bins: list[dict[str, Any]] = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        mask = (conf >= lo) & ((conf < hi) if i < n_bins - 1 else (conf <= hi))
        if not mask.any():
            bins.append({"lo": round(float(lo), 3), "hi": round(float(hi), 3),
                         "n": 0, "confidence": None, "accuracy": None})
            continue
        c = float(conf[mask].mean())
        a = float(correct[mask].mean())
        ece += float(mask.sum()) / len(conf) * abs(c - a)
        bins.append({"lo": round(float(lo), 3), "hi": round(float(hi), 3),
                     "n": int(mask.sum()), "confidence": round(c, 4),
                     "accuracy": round(a, 4)})
    return {
        "n": len(records),
        "ece": round(float(ece), 6),
        "brier": round(brier, 6),
        "nll": round(nll, 6),
        "mean_confidence": round(float(conf.mean()), 6),
        "top1_accuracy": round(top1, 6),
        "bins": bins,
    }


def bootstrap_delta(correct_baseline: Sequence[int], correct_context: Sequence[int],
                    n_boot: int = 1000, seed: int = 17) -> dict[str, Any]:
    """IC bootstrap pareado para la diferencia de acierto (contexto-base)."""
    a = np.asarray(correct_baseline, dtype=float)
    b = np.asarray(correct_context, dtype=float)
    n = len(a)
    if n == 0 or len(b) != n:
        return {"delta": None, "ci_low": None, "ci_high": None, "p_two_sided": None}
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        deltas[i] = b[idx].mean() - a[idx].mean()
    obs = float(b.mean() - a.mean())
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    # p bilateral aproximada: proporción de remuestreos con signo opuesto.
    if obs >= 0:
        p = float((deltas <= 0).mean())
    else:
        p = float((deltas >= 0).mean())
    return {"delta": round(obs, 6), "ci_low": round(float(lo), 6),
            "ci_high": round(float(hi), 6), "p_two_sided": round(p, 6),
            "n_boot": int(n_boot)}


# --------------------------------------------------------------------------
# Autoprueba
# --------------------------------------------------------------------------
def selftest() -> str:
    """Prueba sintética de índice, LLR, neutralidad, guardarraíl y calibración."""
    genes = [
        {"contig": "c", "start": 100, "end": 200, "strand": "+", "gene_name": "a"},
        {"contig": "c", "start": 300, "end": 400, "strand": "+", "gene_name": "b"},
        {"contig": "c", "start": 500, "end": 600, "strand": "-", "gene_name": "c"},
    ]
    idx = GeneIndex(genes)
    assert idx.downstream("c", 50, "+")["gene_name"] == "a"
    assert idx.downstream("c", 250, "+")["gene_name"] == "b"
    assert idx.downstream("c", 650, "-")["gene_name"] == "c"
    assert idx.inside("c", 150, "+")["gene_name"] == "a"
    assert idx.inside("c", 250, "+") is None
    ops = OperonIndex([
        {"contig": "c", "start": 90, "end": 200, "strand": "+", "n_genes": "2"},
        {"contig": "c", "start": 290, "end": 400, "strand": "+", "n_genes": "1"},
    ])
    assert ops.containing_or_downstream("c", 50, "+")["n_genes"] == "2"
    assert ops.containing_or_downstream("c", 150, "+")["n_genes"] == "2"

    table = ContextTable(genome="toy")
    sites = [("c", 50, "+"), ("c", 60, "+"), ("c", 70, "+"), ("c", 80, "+")]
    for i, site in enumerate(sites):
        table.sites[site] = {
            "operon_position": "gene_edge" if i < 2 else "inside_gene",
            "leader": "0-50", "pathway_ec": "1" if i < 2 else "none",
            "domain": "none", "conservation": "none", "expression": "none",
            "operon_size": "1", "sigma54_geometry": "none",
        }
    labels = {sites[0]: {"sigma70"}, sites[1]: {"sigma70"},
              sites[2]: {"sigma24"}, sites[3]: {"sigma24"}}
    model = fit_logodds(table, labels, ["operon_position", "pathway_ec"], min_support=1)
    assert model.llr(table.sites[sites[0]], "sigma70") > 0
    assert model.llr(table.sites[sites[2]], "sigma70") < 0
    # `none` es neutral.
    assert model.llr({"pathway_ec": "none"}, "sigma70") == 0.0
    # Categorías con soporte < min_support son neutrales (anti-sobreajuste).
    strict = fit_logodds(table, labels, ["operon_position"], min_support=5)
    assert strict.llr(table.sites[sites[0]], "sigma70") == 0.0
    # scale=0 reproduce el motivo.
    z = {"sigma70": 3.0, "sigma24": 2.0}
    assert combine(z, model, table.sites[sites[0]], scale=0.0) == z
    # guardarraíl: gap=0 no permite promover al segundo.
    final = {"sigma70": 3.0, "sigma24": 5.0}
    assert select_top1(final, z, gap=0.0) == "sigma70"
    assert select_top1(final, z, gap=None) == "sigma24"
    # geometría sigma54
    assert sigma54_geometry_category({"arch": "-24/-12", "spacer": 11}, "x") == "canonical"
    assert sigma54_geometry_category({"arch": "-24/-12", "spacer": 17}, "x") == "shifted"
    assert sigma54_geometry_category({"arch": "single-box", "dist_to_start": 12}, "x") == "canonical"
    assert sigma54_geometry_category({"arch": "single-box", "dist_to_start": 60}, "x") == "shifted"
    # calibración: predicciones perfectas y confiadas
    perfect = [{"probs": {"sigma70": 1.0, "sigma24": 0.0}, "true": {"sigma70"}},
               {"probs": {"sigma70": 0.0, "sigma24": 1.0}, "true": {"sigma24"}}]
    metrics = calibration_metrics(perfect, sigmas=("sigma70", "sigma24"))
    assert metrics["ece"] is not None and metrics["ece"] < 1e-6
    assert metrics["brier"] == 0.0
    return "contexto: índice, LLR, neutralidad, guardarraíl y calibración OK"


def main(argv: Sequence[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Autoprueba de la capa de contexto sigma.")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
