"""Explicación por familia de rutas para cualquier genoma anotado desde cero."""
from __future__ import annotations

from collections import defaultdict


def explain_families(gaps: list[dict]) -> list[dict]:
    fam: dict[str, list[dict]] = defaultdict(list)
    for g in gaps:
        paths = str(g.get("pathways", "")).split("; ")
        for p in paths:
            if p:
                fam[p].append(g)
        if not paths or not any(paths):
            fam["sin ruta"].append(g)
    out = []
    for name, gs in sorted(fam.items(), key=lambda x: -len(x[1])):
        n = len(gs)
        strong = sum(1 for g in gs if g.get("homolog_status") == "homolog_strong")
        weak = sum(1 for g in gs if g.get("homolog_status") == "homolog_weak")
        no = sum(1 for g in gs if g.get("homolog_status") == "no_homolog")
        tr = sum(1 for g in gs if g.get("transport_status") in
                 ("t1_bacterial", "t1_review", "transport_reaction_only"))
        hom = strong + weak
        if "fatty acid" in name.lower():
            concl, action = "catabolismo/remodelación lipídica", "no priorizar como biosíntesis"
        elif name.startswith("sin ruta"):
            concl, action = "mapeo KEGG incompleto", "completar Rhea→KEGG"
        elif hom / n >= 0.5:
            concl, action = "enzima presente no mapeada", "corregir mapeo/dominio"
        elif tr / n >= 0.3:
            concl, action = "transporte/importación probable", "buscar transportador/T1"
        elif no / n >= 0.5:
            concl, action = "sin homólogo; ausencia/linaje probable", "revisar contexto"
        else:
            concl, action = "mixto; revisar casos", "revisión individual"
        out.append({"family": name, "gaps": n, "homolog": hom,
                    "no_homolog": no, "transport": tr,
                    "conclusion": concl, "action": action})
    return out
