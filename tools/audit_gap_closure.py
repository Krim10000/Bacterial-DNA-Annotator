#!/usr/bin/env python3
"""Audita el cierre/interpretación de brechas metabólicas sin GPU ni red.

Regla central: una brecha **asociada** (a evidencia, a un complejo o a una
reacción de transporte) no es una brecha **cerrada**. Este script no modifica
redes, no escribe tablas de brechas y nunca emite ``confirmed``.

Uso:
    .venv/bin/python tools/audit_gap_closure.py
    .venv/bin/python tools/audit_gap_closure.py --data-root data \
        --out-md reports/gap_closure_audit.md \
        --out-json reports/gap_closure_audit.json

Salidas por defecto:
    reports/gap_closure_audit.md
    reports/gap_closure_audit.json
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

# Estados permitidos. ``confirmed`` no es un estado válido: haría falta un
# experimento (actividad enzimática, complementación, transporte medido, etc.).
ALLOWED_STATES = (
    "gap_evidence_supported",
    "gap_candidate",
    "complex_suspected",
    "complex_incomplete",
    "not_evaluated",
    "not_closable_cpu",
)
STATE_ORDER = ALLOWED_STATES

# Umbral usado por el proyecto para una transferencia de complejo "fuerte".
STRONG_TRANSFER_FRAC = 0.8

# Revisiones que aportan evidencia parcial de cierre. No incluye
# ``one_sided_not_transport`` (es una duda de transporte, no un cierre) ni
# ``keep_missing``/``no_domain_data`` (son ausencia de evidencia).
WEAK_REVIEW_STATUSES = {
    "evidence_high",
    "mapping_context_supported",
    "mapping_rhea_and_ec_supported",
    "mapping_rhea_supported_ec_unknown",
    "weak_rhea_supported",
    "high_priority_homology",
    "high_priority_low",
    "evidence_high_no_rhea",
    "mapping_context_no_producer",
}

# Columnas del contrato de datos. El lector tolera columnas faltantes; la lista
# se usa en los tests y documenta el esquema que se auditó.
GAP_FIELDS = [
    "genome", "gap", "chebi", "support_400", "organism", "phylum", "size_mb",
    "lifestyle", "multifunction", "entry_exit", "status",
    "status_validated_by", "status_reason", "status_attempted",
    "status_confidence", "status_evidence", "multifunction_ec", "role",
    "hmm_evidence", "origin", "audit_status", "audit_reason",
    "transport_status", "transport_reason", "direction_status",
    "direction_reason", "thermo_status", "review_status", "review_reason",
    "homolog_status", "homolog_acc", "homolog_pident", "homolog_qcov",
    "mapping_query", "mapping_rhea", "mapping_acc", "mapping_current_acc",
    "transport_rhea", "transporter_acc", "transport_pident",
    "transport_qcov", "transport_evalue", "mapping_validated",
    "mapping_ec_status", "uniprot_ec_status", "uniprot_ec",
    "uniprot_rhea_status", "uniprot_rhea", "homolog_rhea", "homolog_ec",
    "high_priority",
]

BLOCKER_TEXT = {
    "catabolic": (
        "catabolismo/degradación: no requiere un productor biosintético; "
        "no hay cierre de brecha que buscar"
    ),
    "false_positive_direction": (
        "falso positivo por dirección/moneda: la reacción consume el "
        "metabolito en sentido invertido; no es una brecha biosintética"
    ),
    "entry_point": (
        "punto de entrada/importado: el metabolito puede venir del medio; "
        "no es una brecha biosintética"
    ),
    "spontaneous": (
        "posible reacción espontánea/no enzimática; falta una base curada de "
        "reacciones no enzimáticas"
    ),
    "host_dependent": (
        "linaje endosimbionte/endoparásito: el hospedador puede proveer el "
        "metabolito; requiere experimento en el hospedador"
    ),
    "no_homolog": (
        "sin homólogo detectado: hace falta una base mayor (TrEMBL/UniRef), "
        "estructura o experimento"
    ),
    "novel_no_reference": (
        "química sin referencia en los 400: requiere caracterización "
        "experimental o una base externa"
    ),
    "transport_reaction_only": (
        "reacción de transporte en Rhea pero sin homólogo de transportador: "
        "la importación no está confirmada"
    ),
    "transporter_homolog": (
        "transportador detectado: la importación es una hipótesis; requiere "
        "ensayo de transporte"
    ),
    "direction_unresolved": (
        "dirección no resuelta (Rhea master UN): sin ΔG/experimento no se "
        "puede afirmar que la reacción produce el metabolito"
    ),
}

STATE_TEXT = {
    "gap_evidence_supported": (
        "hay evidencia directa (mapeo Rhea validado, homólogo fuerte o "
        "dominio confirmado) de un cierre candidato; NO es un cierre"
    ),
    "gap_candidate": (
        "hay evidencia parcial (dominio HMM, homólogo débil, revisión, "
        "multifunción same-EC o asociación débil a complejo); NO es un cierre"
    ),
    "complex_suspected": (
        "asociación a un complejo transferido fuerte (frac_mapped>=0.8); "
        "asociación, no cierre"
    ),
    "complex_incomplete": (
        "asociación a un complejo declarado incompleto (faltan subunidades o "
        "estequiometría); no puede cerrar la brecha"
    ),
    "not_evaluated": (
        "sin evidencia, sin complejo, sin transporte y sin bloqueador "
        "explícito: no hay hipótesis de cierre que evaluar"
    ),
    "not_closable_cpu": (
        "hay una capa evaluada que impide cerrar con CPU/offline: no es una "
        "brecha biosintética o requiere experimento/ΔG/base externa/GPU"
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _s(row: dict[str, Any], key: str) -> str:
    """Devuelve un campo como string limpio, tolerando None/columnas ausentes."""
    value = row.get(key, "")
    return str(value).strip() if value is not None else ""


def _f(value: Any) -> float | None:
    """Float tolerante; ``None`` si el campo está vacío o no se puede leer."""
    text = str(value).strip() if value is not None else ""
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _i(value: Any) -> int | None:
    number = _f(value)
    return int(number) if number is not None else None


def _split(value: Any) -> list[str]:
    """Separa listas del proyecto por ``;`` o ``,`` y descarta vacíos."""
    if value is None:
        return []
    return [x.strip() for x in re.split(r"[;,]", str(value)) if x.strip()]


def _rhea_ids(value: Any) -> set[int]:
    """Extrae IDs numéricos de ``RHEA:12345`` (acepta listas)."""
    out: set[int] = set()
    for item in _split(value):
        text = item.upper().replace("RHEA:", "")
        if text.isdigit():
            out.add(int(text))
    return out


def _pct(part: int, total: int) -> str:
    return f"{100.0 * part / total:.1f}%" if total else "0.0%"


def _display_path(path: Path) -> str:
    """Ruta relativa al repo cuando es posible; si no, la ruta original."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------
# Evidencia y clasificación
# ---------------------------------------------------------------------------
def evidence_basis(row: dict[str, Any]) -> tuple[str, list[str]]:
    """Clasifica la evidencia directa de cierre en ``strong``/``weak``/``none``.

    ``strong`` no significa ``confirmed``: significa que hay un candidato con
    mapeo validado, homólogo fuerte o dominio confirmado, todavía sin
    experimento.
    """
    strong: list[str] = []
    if _s(row, "mapping_validated") == "1":
        strong.append("mapping_rhea_validada")
    if _s(row, "homolog_status") == "homolog_strong" and _s(row, "homolog_acc"):
        strong.append("homologo_fuerte")
    if (
        _s(row, "review_status") in {"reannotated_hmm", "domain_confirmed"}
        and _s(row, "hmm_evidence")
    ):
        strong.append("dominio_confirmado")
    if strong:
        return "strong", strong

    weak: list[str] = []
    if _s(row, "hmm_evidence"):
        weak.append("dominio_hmm")
    if _s(row, "homolog_status") in {"homolog_weak", "homolog_low"} and _s(
        row, "homolog_acc"
    ):
        weak.append("homologo_debil")
    if _s(row, "status") == "candidate":
        weak.append("estado_candidate")
    if _s(row, "review_status") in WEAK_REVIEW_STATUSES:
        weak.append("revision_soporte")
    if _s(row, "origin") == "promiscuous":
        weak.append("multifuncion_misma_ec")
    if _s(row, "mapping_rhea"):
        weak.append("rhea_mapeada")
    if weak:
        return "weak", weak
    return "none", []


def gap_blockers(row: dict[str, Any]) -> list[str]:
    """Códigos que impiden cerrar la brecha con la evidencia CPU auditada.

    El orden es la prioridad del motivo primario que se reporta.
    """
    out: list[str] = []
    if (
        _s(row, "audit_status") == "degradation"
        or _s(row, "role") == "catabolic"
        or _s(row, "origin") == "degradation"
    ):
        out.append("catabolic")
    if (
        _s(row, "audit_status") == "false_positive"
        or _s(row, "origin") == "false_positive"
    ):
        out.append("false_positive_direction")
    if _s(row, "role") == "entry_point" or _s(row, "origin") == "imported_or_entry":
        out.append("entry_point")
    if _s(row, "origin") == "spontaneous":
        out.append("spontaneous")
    if _s(row, "status") == "lineage_specific":
        out.append("host_dependent")
    if _s(row, "homolog_status") == "no_homolog":
        out.append("no_homolog")
    if _s(row, "status") == "novel":
        out.append("novel_no_reference")
    if _s(row, "transport_status") == "transport_reaction_only":
        out.append("transport_reaction_only")
    elif _s(row, "transport_status") == "t1_bacterial":
        out.append("transporter_homolog")
    if _s(row, "direction_status") in {"direction_unspecified", "thermo_unknown"}:
        out.append("direction_unresolved")
    return out


def _reason_for_state(
    state: str,
    basis: list[str],
    blockers: list[str],
    associations: list[dict[str, Any]],
) -> str:
    if state == "gap_evidence_supported":
        return (
            "evidencia directa de cierre candidato ("
            + ", ".join(basis)
            + "); no confirmado: requiere experimento"
        )
    if state == "gap_candidate":
        if associations and not basis:
            weak = [a for a in associations if not a.get("incomplete")]
            names = ", ".join(sorted({str(a.get("complex_id", "")) for a in weak}))
            return (
                "asociación débil a complejo ("
                + (names or "sin id")
                + ", frac_mapped<0.8); es una hipótesis, no un cierre"
            )
        return (
            "evidencia parcial ("
            + (", ".join(basis) if basis else "asociación débil")
            + "); no confirmado: requiere validación"
        )
    if state == "complex_suspected":
        strong = [
            a for a in associations
            if a.get("strength") == "strong" and not a.get("incomplete")
        ]
        names = ", ".join(sorted({str(a.get("complex_id", "")) for a in strong}))
        return (
            "asociado a complejo transferido fuerte ("
            + (names or "sin id")
            + "); asociación, no cierre"
        )
    if state == "complex_incomplete":
        incomplete = [a for a in associations if a.get("incomplete")]
        names = ", ".join(
            sorted({str(a.get("complex_id", "")) for a in incomplete})
        )
        return (
            "asociado a complejo incompleto ("
            + (names or "sin id")
            + "); faltan subunidades o estequiometría, no cierra"
        )
    if state == "not_closable_cpu":
        return "; ".join(BLOCKER_TEXT.get(code, code) for code in blockers)
    return STATE_TEXT["not_evaluated"]


def classify_gap(
    row: dict[str, Any],
    associations: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Asigna uno de los seis estados a una brecha.

    Precedencia:
    1. evidencia directa fuerte -> ``gap_evidence_supported``;
    2. evidencia directa débil -> ``gap_candidate``;
    3. complejo fuerte no incompleto -> ``complex_suspected``;
    4. complejo incompleto -> ``complex_incomplete``;
    5. asociación débil a complejo -> ``gap_candidate``;
    6. bloqueador explícito -> ``not_closable_cpu``;
    7. nada -> ``not_evaluated``.

    Nunca devuelve ``confirmed``, aunque la fila de entrada lo diga.
    """
    associations = list(associations or [])
    tier, basis = evidence_basis(row)
    blockers = gap_blockers(row)

    if tier == "strong":
        state = "gap_evidence_supported"
        reason_code = "+".join(basis)
    elif tier == "weak":
        state = "gap_candidate"
        reason_code = "+".join(basis)
    else:
        strong_complete = [
            a for a in associations
            if a.get("strength") == "strong" and not a.get("incomplete")
        ]
        incomplete = [a for a in associations if a.get("incomplete")]
        weak_complete = [a for a in associations if not a.get("incomplete")]
        if strong_complete:
            state = "complex_suspected"
            reason_code = "complex_strong"
        elif incomplete:
            state = "complex_incomplete"
            reason_code = "complex_incomplete"
        elif weak_complete:
            state = "gap_candidate"
            reason_code = "complex_weak"
        elif blockers:
            state = "not_closable_cpu"
            reason_code = blockers[0]
        else:
            state = "not_evaluated"
            reason_code = "no_layers"

    if state not in ALLOWED_STATES:
        raise AssertionError(f"estado no permitido: {state!r}")
    return {
        "state": state,
        "reason_code": reason_code,
        "reason": _reason_for_state(state, basis, blockers, associations),
        "evidence_tier": tier,
        "evidence_basis": basis,
        "blockers": blockers,
        "associated": bool(associations),
        "closed": False,
    }


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------
def load_gaps(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Lee ``gap_interpretation.tsv``. No escribe ni normaliza la tabla."""
    if not path.exists():
        raise FileNotFoundError(f"no existe {path}")
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        rows = [dict(r) for r in csv.DictReader(fh, delimiter="\t")]
    stats = {
        "path": _display_path(path),
        "rows": len(rows),
        "rows_without_chebi": sum(1 for r in rows if not _s(r, "chebi")),
        "genomes": len({_s(r, "genome") for r in rows if _s(r, "genome")}),
        "duplicate_keys": len(rows) - len(
            {(_s(r, "genome"), _s(r, "gap"), _s(r, "chebi")) for r in rows}
        ),
    }
    return rows, stats


def load_bridge(data_root: Path) -> tuple[dict[str, set[str]], dict[str, set[int]], dict[str, Any]]:
    """Carga el puente EC -> Rhea -> productos (ChEBI).

    Solo usa archivos locales. Si faltan, devuelve diccionarios vacíos y deja
    constancia en ``stats``: no se inventan asociaciones.
    """
    reactions_path = data_root / "metabolism" / "rhea_reactions.pkl"
    rhea2ec_path = data_root / "reactions" / "tsv" / "rhea2ec.tsv"
    stats: dict[str, Any] = {
        "available": reactions_path.exists() and rhea2ec_path.exists(),
        "rhea_reactions_path": _display_path(reactions_path),
        "rhea2ec_path": _display_path(rhea2ec_path),
        "reactions": 0,
        "ecs_with_products": 0,
    }
    ec_products: dict[str, set[str]] = defaultdict(set)
    ec_rheas: dict[str, set[int]] = defaultdict(set)
    if not stats["available"]:
        return {}, {}, stats

    with reactions_path.open("rb") as fh:
        reactions = pickle.load(fh)
    stats["reactions"] = len(reactions)

    with rhea2ec_path.open(encoding="utf-8", errors="replace") as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4 or not cols[3] or not cols[0].isdigit():
                continue
            rhea_id = int(cols[0])
            ec = cols[3].strip()
            if not ec:
                continue
            products = reactions.get(rhea_id, {}).get("products", set())
            if products:
                ec_products[ec] |= set(products)
                ec_rheas[ec].add(rhea_id)
    stats["ecs_with_products"] = len(ec_products)
    return dict(ec_products), dict(ec_rheas), stats


def _complex_record(
    row: dict[str, Any],
    source: str,
    genome: str,
    ec_products: dict[str, set[str]],
    ec_rheas: dict[str, set[int]],
) -> dict[str, Any] | None:
    """Convierte una fila de complejo en un registro asociable (o None)."""
    ecs = _split(_s(row, "ecs"))
    if not ecs:
        return None
    products: set[str] = set()
    rheas: set[int] = set()
    for ec in ecs:
        products |= ec_products.get(ec, set())
        rheas |= ec_rheas.get(ec, set())
    if not products and not rheas:
        return None

    if source == "transfer":
        frac = _f(_s(row, "frac_mapped"))
        strength = "strong" if (frac is not None and frac >= STRONG_TRANSFER_FRAC) else "weak"
        expected = _i(_s(row, "expected_n"))
        observed = _i(_s(row, "observed_n"))
        incomplete = _s(row, "stoich_status") == "incomplete" or (
            expected is not None and observed is not None and observed < expected
        )
        complex_id = _s(row, "complex_id")
        complex_name = _s(row, "complex_name")
        evidence_class = _s(row, "evidence_class")
        limitation = _s(row, "limitation")
    else:  # model
        frac = None
        strength = "candidate"
        incomplete = bool(_s(row, "missing"))
        complex_id = _s(row, "complex_id")
        complex_name = _s(row, "name") or _s(row, "complex_name")
        evidence_class = _s(row, "evidence")
        limitation = _s(row, "limitation")

    return {
        "genome": genome,
        "complex_id": complex_id,
        "complex_name": complex_name,
        "source": source,
        "strength": strength,
        "frac_mapped": frac,
        "incomplete": incomplete,
        "stoich_status": _s(row, "stoich_status"),
        "expected_n": _i(_s(row, "expected_n")),
        "observed_n": _i(_s(row, "observed_n")),
        "ecs": ecs,
        "products": products,
        "rheas": rheas,
        "evidence_class": evidence_class,
        "limitation": limitation,
    }


def load_complexes(
    data_root: Path, ec_products: dict[str, set[str]], ec_rheas: dict[str, set[int]]
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Lee ``transfer_*.tsv`` y ``model_*.complexes.tsv``.

    Solo conserva complejos con EC/Rhea y productos: los que no tienen EC no
    pueden asociarse a una brecha por el puente actual. Igual se cuentan en las
    estadísticas para que el límite quede visible.
    """
    complexes_dir = data_root / "complexes"
    transfer_files = sorted(complexes_dir.glob("transfer_*.tsv"))
    model_files = sorted(complexes_dir.glob("model_*.complexes.tsv"))
    by_genome: dict[str, list[dict[str, Any]]] = defaultdict(list)

    transfer_stats: dict[str, Any] = {
        "files": 0, "rows": 0, "rows_with_ecs": 0, "rows_with_products": 0,
        "strong": 0, "weak": 0, "incomplete": 0, "stoich_complete": 0,
        "stoich_incomplete": 0, "stoich_unknown": 0,
    }
    for path in transfer_files:
        genome = path.name.removeprefix("transfer_").removesuffix(".tsv")
        transfer_stats["files"] += 1
        with path.open(encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                transfer_stats["rows"] += 1
                if _s(row, "ecs"):
                    transfer_stats["rows_with_ecs"] += 1
                stoich = _s(row, "stoich_status")
                expected = _i(_s(row, "expected_n"))
                observed = _i(_s(row, "observed_n"))
                incomplete_row = stoich == "incomplete" or (
                    expected is not None and observed is not None
                    and observed < expected
                )
                if stoich == "complete":
                    transfer_stats["stoich_complete"] += 1
                elif stoich == "incomplete":
                    transfer_stats["stoich_incomplete"] += 1
                else:
                    transfer_stats["stoich_unknown"] += 1
                if incomplete_row:
                    transfer_stats["incomplete"] += 1
                frac = _f(_s(row, "frac_mapped"))
                if frac is not None and frac >= STRONG_TRANSFER_FRAC:
                    transfer_stats["strong"] += 1
                else:
                    transfer_stats["weak"] += 1
                record = _complex_record(row, "transfer", genome, ec_products, ec_rheas)
                if record is None:
                    continue
                if record["products"]:
                    transfer_stats["rows_with_products"] += 1
                by_genome[genome].append(record)

    model_stats: dict[str, Any] = {
        "files": 0, "rows": 0, "rows_with_ecs": 0, "rows_with_products": 0,
        "candidate": 0, "with_missing": 0,
    }
    for path in model_files:
        genome = path.name.removeprefix("model_").removesuffix(".complexes.tsv")
        model_stats["files"] += 1
        with path.open(encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                model_stats["rows"] += 1
                if _s(row, "ecs"):
                    model_stats["rows_with_ecs"] += 1
                if _s(row, "missing"):
                    model_stats["with_missing"] += 1
                if _s(row, "status") == "complex_candidate":
                    model_stats["candidate"] += 1
                record = _complex_record(row, "model", genome, ec_products, ec_rheas)
                if record is None:
                    continue
                if record["products"]:
                    model_stats["rows_with_products"] += 1
                by_genome[genome].append(record)

    return dict(by_genome), {"transfer": transfer_stats, "model": model_stats}


def build_associations(
    gap_rows: list[dict[str, Any]],
    complexes_by_genome: dict[str, list[dict[str, Any]]],
) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Asocia brechas a complejos por producto ChEBI exacto.

    La asociación principal es ``EC del complejo -> Rhea -> producto`` y exige
    que el producto sea **exactamente** el ChEBI de la brecha, en el **mismo
    genoma**. El solapamiento de IDs Rhea se guarda como señal de apoyo, pero no
    crea la asociación por sí solo.
    """
    gaps_by_genome: dict[str, set[str]] = defaultdict(set)
    row_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for row in gap_rows:
        genome = _s(row, "genome")
        chebi = _s(row, "chebi")
        if not genome or not chebi:
            continue
        gaps_by_genome[genome].add(chebi)
        row_by_key[(genome, chebi)] = row

    associations: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for genome, complexes in complexes_by_genome.items():
        gaps = gaps_by_genome.get(genome)
        if not gaps:
            continue
        for complex_row in complexes:
            hit = gaps & complex_row.get("products", set())
            if not hit:
                continue
            complex_rheas = complex_row.get("rheas", set())
            for chebi in hit:
                gap_row = row_by_key[(genome, chebi)]
                gap_rheas = _rhea_ids(
                    _s(gap_row, "mapping_rhea") + ";" + _s(gap_row, "transport_rhea")
                )
                associations[(genome, chebi)].append({
                    "genome": genome,
                    "gap": _s(gap_row, "gap"),
                    "chebi": chebi,
                    "complex_id": complex_row.get("complex_id", ""),
                    "complex_name": complex_row.get("complex_name", ""),
                    "source": complex_row.get("source", ""),
                    "strength": complex_row.get("strength", ""),
                    "frac_mapped": complex_row.get("frac_mapped"),
                    "incomplete": bool(complex_row.get("incomplete")),
                    "stoich_status": complex_row.get("stoich_status", ""),
                    "ecs": list(complex_row.get("ecs", [])),
                    "rhea_overlap": bool(gap_rheas & complex_rheas),
                    "evidence_class": complex_row.get("evidence_class", ""),
                    "limitation": complex_row.get("limitation", ""),
                })
    for records in associations.values():
        records.sort(key=lambda r: (r["source"], r["complex_id"]))
    return dict(associations)


# ---------------------------------------------------------------------------
# Auditoría y reporte
# ---------------------------------------------------------------------------
def _association_summary(
    associations: dict[tuple[str, str], list[dict[str, Any]]],
) -> dict[str, Any]:
    strong = incomplete = both = 0
    rhea_overlap = 0
    for records in associations.values():
        has_strong_complete = any(
            r["strength"] == "strong" and not r["incomplete"] for r in records
        )
        has_incomplete = any(r["incomplete"] for r in records)
        if has_strong_complete:
            strong += 1
        if has_incomplete:
            incomplete += 1
        if has_strong_complete and has_incomplete:
            both += 1
        if any(r.get("rhea_overlap") for r in records):
            rhea_overlap += 1
    return {
        "gaps_with_any": len(associations),
        "gaps_with_strong": strong,
        "gaps_with_incomplete": incomplete,
        "gaps_with_strong_and_incomplete": both,
        "gaps_with_rhea_overlap": rhea_overlap,
        "association_count": sum(len(v) for v in associations.values()),
    }


def run_audit(
    gap_rows: list[dict[str, Any]],
    complexes_by_genome: dict[str, list[dict[str, Any]]],
    bridge_stats: dict[str, Any] | None = None,
    complex_stats: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Clasifica todas las brechas y devuelve el resultado estructurado."""
    associations = build_associations(gap_rows, complexes_by_genome)

    state_counts: Counter[str] = Counter()
    evidence_tier_counts: Counter[str] = Counter()
    basis_counts: Counter[str] = Counter()
    first_blocker_counts: Counter[str] = Counter()
    all_blocker_counts: Counter[str] = Counter()
    signals: Counter[str] = Counter()
    per_genome: dict[str, dict[str, Any]] = {}
    association_rows: list[dict[str, Any]] = []
    not_evaluated_rows: list[dict[str, Any]] = []
    confirmed_inputs_ignored = 0

    for row in gap_rows:
        genome = _s(row, "genome")
        chebi = _s(row, "chebi")
        records = associations.get((genome, chebi), []) if genome and chebi else []
        result = classify_gap(row, records)
        state = result["state"]
        state_counts[state] += 1
        evidence_tier_counts[result["evidence_tier"]] += 1
        for basis in result["evidence_basis"]:
            basis_counts[basis] += 1
        if state == "not_closable_cpu" and result["blockers"]:
            first_blocker_counts[result["blockers"][0]] += 1
            for blocker in result["blockers"]:
                all_blocker_counts[blocker] += 1
        if state == "not_evaluated":
            not_evaluated_rows.append({
                "genome": genome,
                "gap": _s(row, "gap"),
                "chebi": chebi,
                "status": _s(row, "status"),
                "role": _s(row, "role"),
                "origin": _s(row, "origin"),
                "audit_status": _s(row, "audit_status"),
                "transport_status": _s(row, "transport_status"),
                "direction_status": _s(row, "direction_status"),
                "review_status": _s(row, "review_status"),
            })
        if _s(row, "status") == "confirmed" or _s(row, "review_status") == "confirmed":
            confirmed_inputs_ignored += 1

        has_evidence = result["evidence_tier"] != "none"
        has_complex = bool(records)
        has_transport = _s(row, "transport_status") in {
            "transport_reaction_only", "t1_bacterial",
        }
        has_direction_resolved = _s(row, "direction_status") in {
            "thermo_forward", "thermo_reversible",
        }
        has_direction_unresolved = _s(row, "direction_status") in {
            "direction_unspecified", "thermo_unknown",
        }
        has_multifunction = _s(row, "origin") == "promiscuous"
        if has_evidence:
            signals["has_direct_evidence"] += 1
        if result["evidence_tier"] == "strong":
            signals["has_strong_evidence"] += 1
        if result["evidence_tier"] == "weak":
            signals["has_weak_evidence"] += 1
        if has_complex:
            signals["has_complex_association"] += 1
        if has_transport:
            signals["has_transport_signal"] += 1
        if _s(row, "transport_status") == "no_transport":
            signals["transport_evaluated_no_homolog"] += 1
        if has_direction_resolved:
            signals["has_direction_resolved"] += 1
        if has_direction_unresolved:
            signals["has_direction_unresolved"] += 1
        if has_multifunction:
            signals["has_multifunction_hypothesis"] += 1

        positive = (
            has_evidence or has_complex or has_transport
            or has_direction_resolved or has_multifunction
        )
        if not positive:
            signals["without_positive_signal"] += 1
        any_layer = any(
            _s(row, key)
            for key in (
                "audit_status", "transport_status", "direction_status",
                "review_status", "origin",
            )
        ) or has_evidence or has_complex
        if not any_layer:
            signals["without_any_layer"] += 1

        entry = per_genome.setdefault(
            genome,
            {
                "genome": genome,
                "gaps": 0,
                "states": Counter(),
                "with_evidence": 0,
                "with_complex": 0,
                "with_transport": 0,
                "with_direction": 0,
            },
        )
        entry["gaps"] += 1
        entry["states"][state] += 1
        entry["with_evidence"] += int(has_evidence)
        entry["with_complex"] += int(has_complex)
        entry["with_transport"] += int(has_transport)
        entry["with_direction"] += int(has_direction_resolved)

        for record in records:
            association_rows.append({
                **record,
                "state": state,
                "reason_code": result["reason_code"],
            })

    per_genome_rows = []
    for genome in sorted(per_genome):
        entry = per_genome[genome]
        per_genome_rows.append({
            **{k: v for k, v in entry.items() if k != "states"},
            "states": dict(entry["states"]),
        })

    total = len(gap_rows)
    closure = {
        "confirmed": 0,
        "closed": 0,
        "associated_not_closed": len(associations),
        "policy": (
            "ninguna brecha se declara cerrada ni confirmada; asociación != "
            "cierre; un cierre exige experimento (actividad, complementación, "
            "transporte medido) o validación equivalente"
        ),
        "confirmed_input_rows_ignored": confirmed_inputs_ignored,
    }
    return {
        "totals": {
            "gaps": total,
            "genomes": len(per_genome),
            "states_sum": sum(state_counts.values()),
        },
        "states": {state: state_counts.get(state, 0) for state in STATE_ORDER},
        "state_definitions": STATE_TEXT,
        "evidence": {
            "strong": evidence_tier_counts.get("strong", 0),
            "weak": evidence_tier_counts.get("weak", 0),
            "none": evidence_tier_counts.get("none", 0),
            "basis": dict(basis_counts),
        },
        "signals": {key: signals.get(key, 0) for key in (
            "has_direct_evidence", "has_strong_evidence", "has_weak_evidence",
            "has_complex_association", "has_transport_signal",
            "transport_evaluated_no_homolog", "has_direction_resolved",
            "has_direction_unresolved", "has_multifunction_hypothesis",
            "without_positive_signal", "without_any_layer",
        )},
        "blockers": {
            "first_reason": dict(first_blocker_counts),
            "all_reasons": dict(all_blocker_counts),
            "groups": {
                "not_a_biosynthetic_gap": sum(
                    first_blocker_counts.get(code, 0)
                    for code in ("catabolic", "false_positive_direction",
                                 "entry_point", "spontaneous")
                ),
                "biosynthetic_gap_not_closable": sum(
                    first_blocker_counts.get(code, 0)
                    for code in ("direction_unresolved", "transport_reaction_only",
                                 "transporter_homolog", "no_homolog",
                                 "novel_no_reference", "host_dependent")
                ),
            },
        },
        "complexes": {
            "bridge": bridge_stats or {},
            "transfer": (complex_stats or {}).get("transfer", {}),
            "model": (complex_stats or {}).get("model", {}),
            "associations": _association_summary(associations),
        },
        "closure": closure,
        "association_list": association_rows,
        "not_evaluated_detail": not_evaluated_rows,
        "per_genome": per_genome_rows,
        "pct": {state: _pct(state_counts.get(state, 0), total) for state in STATE_ORDER},
    }


def _write_markdown(result: dict[str, Any], path: Path) -> None:
    total = result["totals"]["gaps"]
    states = result["states"]
    signals = result["signals"]
    complexes = result["complexes"]
    transfer = complexes.get("transfer", {})
    model = complexes.get("model", {})
    assoc = complexes.get("associations", {})
    bridge = complexes.get("bridge", {})
    lines: list[str] = []
    add = lines.append

    add("# Auditoría de cierre/interpretación de brechas")
    add("")
    add(
        "Generado por `tools/audit_gap_closure.py`. Frente CPU, sin GPU, sin red "
        "y sin modificar redes ni tablas de entrada."
    )
    add("")
    add("## Regla")
    add("")
    add("- **Asociada no es cerrada.** Una brecha puede tener evidencia, un")
    add("  complejo o un transportador asociado y seguir abierta.")
    add("- **No se emite `confirmed`.** Un cierre real exige experimento")
    add("  (actividad enzimática, complementación, transporte medido) o una")
    add("  validación equivalente.")
    add(f"- Cierres declarados en esta auditoría: **{result['closure']['closed']}**.")
    add(
        f"- Brechas asociadas a algo pero no cerradas: "
        f"**{result['closure']['associated_not_closed']}**."
    )
    add("")
    add("## Instantánea")
    add("")
    add(f"- Brechas auditadas: **{total}** en **{result['totals']['genomes']}** genomas.")
    add(f"- Generado (UTC): `{result.get('generated_at', '')}`.")
    add(f"- Entrada: `data/metabolism/gap_interpretation.tsv`.")
    add(
        f"- Puente EC→Rhea→productos: "
        f"{'disponible' if bridge.get('available') else 'NO disponible'}; "
        f"{bridge.get('reactions', 0)} reacciones Rhea y "
        f"{bridge.get('ecs_with_products', 0)} EC con productos."
    )
    add("")
    add("## Estados")
    add("")
    add("| estado | brechas | % | criterio |")
    add("|---|---:|---:|---|")
    for state in STATE_ORDER:
        add(f"| `{state}` | {states[state]} | {result['pct'][state]} | "
            f"{STATE_TEXT[state]} |")
    add(f"| **total** | **{total}** | 100.0% | |")
    add("")
    add(
        "Precedencia: evidencia fuerte > evidencia débil > complejo fuerte > "
        "complejo incompleto > asociación débil > bloqueador > sin evaluar."
    )
    add("")
    add("### Evidencia directa")
    add("")
    ev = result["evidence"]
    strong_basis = {
        key: value for key, value in ev["basis"].items()
        if key in {"mapping_rhea_validada", "homologo_fuerte", "dominio_confirmado"}
    }
    weak_basis = {
        key: value for key, value in ev["basis"].items()
        if key not in {"mapping_rhea_validada", "homologo_fuerte", "dominio_confirmado"}
    }
    strong_text = ", ".join(f"{k}={v}" for k, v in strong_basis.items()) or "—"
    weak_text = ", ".join(f"{k}={v}" for k, v in weak_basis.items()) or "—"
    add("| nivel | brechas | base |")
    add("|---|---:|---|")
    add(f"| fuerte | {ev['strong']} | {strong_text} |")
    add(f"| débil | {ev['weak']} | {weak_text} |")
    add(f"| ninguna | {ev['none']} | — |")
    add("")
    add(
        "`fuerte` significa candidato con mapeo validado, homólogo fuerte o "
        "dominio confirmado; **no** significa función demostrada."
    )
    add("")
    ne_rows = result.get("not_evaluated_detail", [])
    if ne_rows:
        add(f"### Filas `not_evaluated` ({len(ne_rows)})")
        add("")
        add("| genoma | brecha | ChEBI | review | motivo |")
        add("|---|---|---|---|---|")
        for item in ne_rows:
            add(
                f"| `{item.get('genome', '')}` | {item.get('gap', '')} | "
                f"`{item.get('chebi', '')}` | `{item.get('review_status', '') or '—'}` "
                "| sin evidencia, complejo, transporte ni bloqueador explícito |"
            )
        add("")
    add("## Cobertura de señales")
    add("")
    add("| señal | brechas |")
    add("|---|---:|")
    labels = {
        "has_direct_evidence": "evidencia directa (fuerte o débil)",
        "has_strong_evidence": "evidencia fuerte",
        "has_weak_evidence": "evidencia débil",
        "has_complex_association": "asociación a complejo (transfer/model)",
        "has_transport_signal": "señal de transporte (reacción o transportador)",
        "transport_evaluated_no_homolog": "transporte evaluado sin homólogo",
        "has_direction_resolved": "dirección resuelta (forward/reversible)",
        "has_direction_unresolved": "dirección no resuelta (UN/unknown)",
        "has_multifunction_hypothesis": "hipótesis de multifunción (misma EC)",
        "without_positive_signal": "sin señal positiva de ningún tipo",
        "without_any_layer": "sin ninguna capa de interpretación",
    }
    for key, label in labels.items():
        add(f"| {label} | {signals.get(key, 0)} |")
    add("")
    add(
        "`sin señal positiva` no es lo mismo que `not_evaluated`: la fila puede "
        "tener un resultado negativo (p. ej. `no_transport`) o un bloqueador."
    )
    add("")
    add("## Complejos")
    add("")
    add("| insumo | archivos | filas | con EC | con productos | fuertes | incompletos |")
    add("|---|---:|---:|---:|---:|---:|---:|")
    add(
        f"| `transfer_*` | {transfer.get('files', 0)} | {transfer.get('rows', 0)} | "
        f"{transfer.get('rows_with_ecs', 0)} | {transfer.get('rows_with_products', 0)} | "
        f"{transfer.get('strong', 0)} | {transfer.get('incomplete', 0)} |"
    )
    add(
        f"| `model_*` | {model.get('files', 0)} | {model.get('rows', 0)} | "
        f"{model.get('rows_with_ecs', 0)} | {model.get('rows_with_products', 0)} | "
        f"— | {model.get('with_missing', 0)} (con `missing`) |"
    )
    add("")
    add("| asociación | brechas |")
    add("|---|---:|")
    add(f"| con cualquier complejo | {assoc.get('gaps_with_any', 0)} |")
    add(f"| con complejo fuerte (frac_mapped≥0.8) | {assoc.get('gaps_with_strong', 0)} |")
    add(f"| con complejo incompleto | {assoc.get('gaps_with_incomplete', 0)} |")
    add(f"| fuerte **y** además incompleto | {assoc.get('gaps_with_strong_and_incomplete', 0)} |")
    add(f"| con solapamiento Rhea | {assoc.get('gaps_with_rhea_overlap', 0)} |")
    add(f"| asociaciones totales (una brecha puede tener varias) | {assoc.get('association_count', 0)} |")
    add("")
    if assoc.get("gaps_with_any", 0):
        add("### Brechas asociadas a complejos")
        add("")
        add("| genoma | brecha | ChEBI | complejo | fuente | frac | estequiometría | estado |")
        add("|---|---|---|---|---|---:|---|---|")
        for record in result["association_list"]:
            frac = record.get("frac_mapped")
            frac_text = f"{frac:.2f}" if isinstance(frac, (int, float)) else "—"
            stoich = "incompleto" if record.get("incomplete") else (
                record.get("stoich_status") or "—"
            )
            add(
                f"| `{record.get('genome', '')}` | {record.get('gap', '')} | "
                f"`{record.get('chebi', '')}` | `{record.get('complex_id', '')}` "
                f"| {record.get('source', '')} | {frac_text} | {stoich} | "
                f"`{record.get('state', '')}` |"
            )
        add("")
    add(
        "**Hallazgo:** los `model_*` no aportan asociaciones en esta instantánea "
        f"porque tienen {model.get('rows_with_ecs', 0)} filas con EC. No se "
        "fuerza una asociación por vecindad: sin EC/Rhea no hay producto que "
        "comparar con el ChEBI de la brecha."
    )
    add("")
    add("## Bloqueadores de `not_closable_cpu`")
    add("")
    add("| motivo primario | brechas | significado |")
    add("|---|---:|---|")
    for code, count in sorted(
        result["blockers"]["first_reason"].items(), key=lambda kv: (-kv[1], kv[0])
    ):
        add(f"| `{code}` | {count} | {BLOCKER_TEXT.get(code, code)} |")
    add("")
    add("Conteo por todos los motivos (una brecha puede tener varios):")
    add("")
    add("| motivo | brechas |")
    add("|---|---:|")
    for code, count in sorted(
        result["blockers"]["all_reasons"].items(), key=lambda kv: (-kv[1], kv[0])
    ):
        add(f"| `{code}` | {count} |")
    add("")
    add(
        "`catabolic`, `false_positive_direction`, `entry_point` y `spontaneous` "
        "no son brechas biosintéticas: no hay cierre que buscar. El resto son "
        "brechas reales sin cierre CPU/offline posible."
    )
    groups = result["blockers"].get("groups", {})
    add("")
    add(
        f"Resumen: **{groups.get('not_a_biosynthetic_gap', 0)}** filas no son "
        f"brechas biosintéticas y **{groups.get('biosynthetic_gap_not_closable', 0)}** "
        "son brechas reales sin cierre CPU/offline."
    )
    add("")
    add("## Límites")
    add("")
    add("- El puente es **EC del complejo → Rhea → productos (ChEBI)**, con")
    add("  intersección exacta de ChEBI. No usa texto libre, estructura ni")
    add("  similitud de ligandos.")
    add("- La asociación puede incluir metabolitos moneda si una EC del complejo")
    add("  los produce; se reporta como asociación, nunca como cierre.")
    add("- `frac_mapped` es ortología K-12, no interacción física. `unknown` en")
    add("  estequiometría no es `complete`.")
    add("- `gap_evidence_supported` agrupa mapeo Rhea validado, homólogo fuerte y")
    add("  dominio confirmado; ninguno demuestra actividad en este genoma.")
    add("- No se calcula ΔG ni se ejecuta Equilibrator: `direction_unresolved`")
    add("  queda como bloqueador de esta auditoría.")
    add("- No hay GPU, AlphaFold-Multimer, co-folding ni dinámica; una predicción")
    add("  de estructura tampoco cerraría una brecha.")
    add("- `not_closable_cpu` no significa imposible: significa que con las capas")
    add("  CPU/offline auditadas no se puede emitir un cierre.")
    add("- La instantánea es de solo lectura; `gap_interpretation.tsv` y")
    add("  `data/complexes/` no se modifican.")
    add("")
    add("## Reproducir")
    add("")
    add("```bash")
    add(".venv/bin/python tools/audit_gap_closure.py \\")
    add("    --data-root data \\")
    add("    --out-md reports/gap_closure_audit.md \\")
    add("    --out-json reports/gap_closure_audit.json")
    add(".venv/bin/python -m pytest tests/test_audit_gap_closure.py -q")
    add("```")
    add("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def audit(
    data_root: Path | str = ROOT / "data",
    out_md: Path | str | None = ROOT / "reports" / "gap_closure_audit.md",
    out_json: Path | str | None = ROOT / "reports" / "gap_closure_audit.json",
) -> dict[str, Any]:
    """Ejecuta la auditoría completa. Solo lee las entradas."""
    data_root = Path(data_root)
    gap_rows, gap_stats = load_gaps(
        data_root / "metabolism" / "gap_interpretation.tsv"
    )
    ec_products, ec_rheas, bridge_stats = load_bridge(data_root)
    complexes_by_genome, complex_stats = load_complexes(
        data_root, ec_products, ec_rheas
    )
    result = run_audit(
        gap_rows, complexes_by_genome, bridge_stats=bridge_stats,
        complex_stats=complex_stats,
    )
    result["inputs"] = {
        "data_root": _display_path(data_root),
        "gap_interpretation": gap_stats,
        "bridge": bridge_stats,
    }
    result["generated_at"] = datetime.now(timezone.utc).isoformat()
    if out_md is not None:
        _write_markdown(result, Path(out_md))
    if out_json is not None:
        out_json = Path(out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(
            json.dumps(result, ensure_ascii=False, indent=2, default=str) + "\n",
            encoding="utf-8",
        )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audita el estado de cierre/interpretación de brechas sin GPU ni "
            "red. No modifica redes ni declara cierres."
        )
    )
    parser.add_argument(
        "--data-root", type=Path, default=ROOT / "data",
        help="directorio data/ (default: %(default)s)",
    )
    parser.add_argument(
        "--out-md", type=Path,
        default=ROOT / "reports" / "gap_closure_audit.md",
        help="reporte Markdown (default: %(default)s)",
    )
    parser.add_argument(
        "--out-json", type=Path,
        default=ROOT / "reports" / "gap_closure_audit.json",
        help="reporte JSON (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    result = audit(args.data_root, args.out_md, args.out_json)
    summary = {
        "gaps": result["totals"]["gaps"],
        "genomes": result["totals"]["genomes"],
        "states": result["states"],
        "complex_associations": result["complexes"]["associations"],
        "closure": result["closure"],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
