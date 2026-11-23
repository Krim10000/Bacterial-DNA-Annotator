"""Interpretación de brechas para cualquier genoma (proyecto o subido).

Asigna rol, origen, estado y fuente de evidencia. La dirección queda
`direction_unspecified` salvo evidencia termodinámica/fisiológica; no se
infiere de la reacción maestra de Rhea.
"""
from __future__ import annotations

import re
from pathlib import Path

ENTRY = re.compile(r"water|h2o|proton|h\(\+\)|carbon dioxide|co2|oxygen|o2|ammonia|"
                   r"ammonium|phosphate|sulfate|sulfite|nitrate|nitrite|nad|nadp|fad|fmn|"
                   r"\bcoa\b|atp|adp|amp|gtp|gdp|sodium|potassium|chloride|magnesium|"
                   r"calcium|iron|zinc|manganese|copper|molybdenum|nickel|cobalt", re.I)
EXOGENOUS = re.compile(r"stachyose|lactose|galactose|sucrose|maltose|trehalose|xylose|"
                       r"arabinose|mannose|fructose|lyxonate|gluconate|glucarate|curcumin|"
                       r"ochratoxin|verruculide|chrodrimanin|pheophorbide|octopamine|"
                       r"acetylcholine|coenzyme m|antibiotic|sulfogalactosyl", re.I)
LIPID = re.compile(r"acyl|glycero|phospho|eicos|octadec|hexadec|tetracosa|hexacos|"
                   r"octacos|triacont|phytanoyl|cholan|cucurbit|alkyl|epoxy|glycerol", re.I)
SPONT = re.compile(r"hydrogen peroxide|superoxide|radical|formaldehyde|glyoxal|"
                   r"methylglyoxal|carbamate|cyanate|isocyanate|maillard", re.I)


def interpret_gaps(gaps: list[dict], organism: str = "", ref_size_mb: float = 0) -> list[dict]:
    for g in gaps:
        name = g.get("gap", "")
        support = int(g.get("support", 0) or 0)
        if g.get("entry_exit") or ENTRY.search(name):
            role, origin = "entry_point", "imported_or_entry"
            status, why = "anomalous", "punto de entrada/salida no cuestionado"
        elif SPONT.search(name):
            role, origin = "unknown", "spontaneous"
            status, why = "anomalous", "química no enzimática conocida"
        elif LIPID.search(name):
            role, origin = "anabolic", "lipid_remodeling"
            status, why = ("validated" if support >= 10 else "consensus" if support else "novel",
                           "intermediario lipídico")
        elif EXOGENOUS.search(name):
            role, origin = "entry_point", "imported_or_entry"
            status, why = "anomalous", "compuesto exógeno/importado probable"
        else:
            role, origin = "anabolic", "enzymatic_missing"
            status = "validated" if support >= 10 else "consensus" if support else "novel"
            why = f"soporte {support}/400; enzima no mapeada"
        g.update({
            "role": role, "origin": origin, "status": status,
            "status_reason": why,
            "status_validated_by": "Rhea + 400 genomas",
            "status_attempted": "traza Rhea; soporte 400",
            "status_confidence": "high" if support >= 10 else "medium" if support else "low",
            "direction_status": "direction_unspecified",
            "direction_reason": "Rhea master UN: lados arbitrarios; requiere ΔG",
            "thermo_status": "pending",
            "transport_status": "not_checked",
            "evidence_source": "Rhea + 400 genomas",
            "review_status": "",
            # El pipeline actual evalua proteinas individuales: no ensambla
            # complejos ni exige estequiometria de subunidades.
            "complex_status": "not_evaluated",
            "complex_evidence": "",
            "quaternary_status": "unknown",
            "complex_note": (
                "No evaluado: la red se infiere de proteinas individuales; un "
                "complejo obligado puede aparecer como brecha sin que falte la "
                "capacidad enzimatica."
            ),
        })
    return gaps


def annotate_pathways(gaps: list[dict], data_root) -> list[dict]:
    """Añade familias KEGG a cada brecha a partir de sus productores Rhea."""
    path = Path(data_root) / "kegg" / "pathway_reactions.tsv"
    if not path.exists():
        return gaps
    pmap: dict[str, list[str]] = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0] != "pathway_id":
            pmap.setdefault(p[2], []).append(f"{p[0]} {p[1]}")
    for g in gaps:
        prod = g.get("producers", [])
        if isinstance(prod, str):
            prod = [x for x in prod.split(";") if x]
        paths = []
        for r in prod:
            paths += pmap.get(r, [])
        g["pathways"] = "; ".join(sorted(set(paths))) if paths else "sin ruta"
    return gaps


def classify_review(gaps: list[dict], data_root=None) -> list[dict]:
    """Clasificación de revisión para un FASTA nuevo (sin panel HMM completo)."""
    for g in gaps:
        tr = g.get("transport_status", "")
        hom = g.get("homolog_status", "")
        if tr in ("t1_bacterial", "t1_review", "transport_reaction_only"):
            cls = "transport"
        elif hom in ("homolog_strong", "homolog_weak"):
            cls = "annotation_gap"
        else:
            cls = "unknown"
        g["review_class"] = cls
    if data_root is not None:
        try:
            from exception_logic import annotate
            gaps = annotate(gaps, data_root)
        except Exception:
            pass
    return gaps
