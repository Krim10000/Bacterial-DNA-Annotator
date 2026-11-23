"""Complejos bacterianos por vecindad génica, sin GPU.

No predice interfaces ni estructura cuaternaria de novo. Agrupa CDS cercanos
(mismo contig y hebra), compara con firmas curadas de complejos y declara
completitud. La proximidad es una condición necesaria, no una prueba.

La tabla `complexes.tsv` y las columnas de complejo en brechas/proteínas son
una capa de interpretación: no cierran una brecha ni convierten una hipótesis
en `confirmed`.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any

# Firmas mínimas, basadas en nombres de genes y literatura clásica. No son un
# catálogo exhaustivo; se amplían con Complex Portal/EcoCyc/KEGG si la licencia
# lo permite. `essential` evita exigir subunidades accesorias.
SIGNATURES: list[dict[str, Any]] = [
    {
        "id": "ATP_SYNTHASE", "name": "ATP synthase (F1Fo)",
        "genes": ["atpA", "atpB", "atpC", "atpD", "atpE", "atpF", "atpG",
                  "atpH", "atpI"],
        "essential": ["atpA", "atpB", "atpC", "atpD"],
        "source": "curated literature; operon atp",
    },
    {
        "id": "NADH_DEHYDROGENASE", "name": "NADH dehydrogenase I (Nuo)",
        "genes": [f"nuo{c}" for c in "ABCDEFGHIJKLMN"],
        "essential": ["nuoB", "nuoC", "nuoD", "nuoE", "nuoF", "nuoG",
                      "nuoH", "nuoI", "nuoJ", "nuoK", "nuoL", "nuoM", "nuoN"],
        "source": "curated literature; operon nuo",
    },
    {
        "id": "SUCCINATE_DEHYDROGENASE", "name": "Succinate dehydrogenase",
        "genes": ["sdhA", "sdhB", "sdhC", "sdhD"],
        "essential": ["sdhA", "sdhB", "sdhC", "sdhD"],
        "source": "curated literature; operon sdh",
    },
    {
        "id": "CYTOCHROME_BC1", "name": "Cytochrome bc1 complex",
        "genes": ["petA", "petB", "petC"],
        "essential": ["petA", "petB", "petC"],
        "source": "curated literature; operon pet",
    },
    {
        "id": "CYTOCHROME_OXIDASE", "name": "Cytochrome c oxidase",
        "genes": ["cyoA", "cyoB", "cyoC", "cyoD"],
        "essential": ["cyoA", "cyoB", "cyoC"],
        "source": "curated literature; operon cyo",
    },
    {
        "id": "PYRUVATE_DEHYDROGENASE", "name": "Pyruvate dehydrogenase",
        "genes": ["aceE", "aceF", "lpdA"],
        "essential": ["aceE", "aceF", "lpdA"],
        "source": "curated literature; operon ace",
    },
    {
        "id": "OXOGLUTARATE_DEHYDROGENASE", "name": "2-oxoglutarate dehydrogenase",
        "genes": ["sucA", "sucB", "lpdA"],
        "essential": ["sucA", "sucB", "lpdA"],
        "source": "curated literature; operon suc",
    },
    {
        "id": "NITROGENASE", "name": "Nitrogenase",
        "genes": ["nifH", "nifD", "nifK"],
        "essential": ["nifH", "nifD", "nifK"],
        "source": "curated literature; operon nif",
    },
    {
        "id": "SEC_SYSTEM", "name": "Sec secretion system",
        "genes": ["secA", "secY", "secE", "secG"],
        "essential": ["secA", "secY", "secE"],
        "source": "curated literature; sec",
    },
    {
        "id": "TAT_SYSTEM", "name": "Tat secretion system",
        "genes": ["tatA", "tatB", "tatC"],
        "essential": ["tatA", "tatB", "tatC"],
        "source": "curated literature; tat",
    },
    {
        "id": "FTS_RING", "name": "Fts divisome ring",
        "genes": ["ftsA", "ftsZ", "ftsI"],
        "essential": ["ftsA", "ftsZ"],
        "source": "curated literature; fts",
    },
]


@dataclass
class Complex:
    complex_id: str
    name: str
    contig: str
    start: int
    end: int
    subunits: str
    expected: str
    missing: str
    status: str
    evidence: str
    confidence: str
    source: str
    limitation: str = (
        "Vecindad génica; no demuestra interacción física ni estequiometría. "
        "No usar como confirmación."
    )


def _genes_of(features: list[dict], hits: dict[str, dict]) -> list[dict]:
    out = []
    for f in features:
        if f.get("type") != "CDS":
            continue
        h = hits.get(f.get("id", ""), {})
        gene = str(h.get("gene", "") or "").lower()
        if gene:
            out.append({**f, "gene": gene, "hit": h})
    return out


def neighborhoods(features: list[dict], max_gap: int = 300) -> list[list[dict]]:
    """Agrupa CDS por contig/hebra; un RNA o un hueco grande corta el grupo."""
    by_key: dict[tuple[str, str], list[dict]] = {}
    for f in features:
        by_key.setdefault((str(f.get("contig", "")), str(f.get("strand", "+"))),
                          []).append(f)
    groups: list[list[dict]] = []
    for (_contig, _strand), items in by_key.items():
        items = sorted(items, key=lambda x: int(x.get("start", 0)))
        current: list[dict] = []
        last_end = None
        for f in items:
            if f.get("type") != "CDS":
                if current:
                    groups.append(current)
                current, last_end = [], None
                continue
            start, end = int(f.get("start", 0)), int(f.get("end", 0))
            if last_end is not None and start - last_end > max_gap:
                if current:
                    groups.append(current)
                current = []
            current.append(f)
            last_end = max(last_end or 0, end)
        if current:
            groups.append(current)
    return groups


def infer(features: list[dict], hits: dict[str, dict], genome: str = "",
          max_gap: int = 300) -> dict[str, Any]:
    """Devuelve complejos candidatos y asociación por proteína."""
    genes = _genes_of(features, hits)
    complexes: list[Complex] = []
    by_protein: dict[str, dict[str, str]] = {}
    gene_by_id = {g["id"]: g for g in genes}
    for group in neighborhoods(features, max_gap=max_gap):
        present = {gene_by_id[f["id"]]["gene"]: gene_by_id[f["id"]]
                   for f in group if f.get("id") in gene_by_id}
        for sig in SIGNATURES:
            expected = [g.lower() for g in sig["genes"]]
            found = [g for g in expected if g in present]
            if len(found) < 2:
                continue
            missing = [g for g in expected if g not in present]
            essential = [g.lower() for g in sig["essential"]]
            complete = all(g in present for g in essential)
            if complete and len(found) >= max(2, len(essential)):
                status, conf = "complex_complete", "high" if len(found) >= 3 else "medium"
            elif complete:
                status, conf = "complex_complete", "medium"
            elif len(found) >= 2:
                status, conf = "complex_incomplete", "low"
            else:
                continue
            coords = [(int(present[g]["start"]), int(present[g]["end"])) for g in found]
            start, end = min(x[0] for x in coords), max(x[1] for x in coords)
            cid = f"{sig['id']}:{group[0].get('contig','')}:{start}"
            subs = ";".join(f"{g}:{present[g].get('id','')}" for g in found)
            complexes.append(Complex(
                complex_id=cid, name=sig["name"],
                contig=str(group[0].get("contig", "")), start=start, end=end,
                subunits=subs, expected=";".join(expected),
                missing=";".join(missing), status=status,
                evidence=f"same contig/strand; gap<={max_gap} bp; genes={len(found)}/{len(expected)}",
                confidence=conf, source=sig["source"],
            ))
            for g in found:
                pid = str(present[g].get("id", ""))
                if pid:
                    by_protein[pid] = {
                        "complex_id": cid, "complex_name": sig["name"],
                        "complex_status": status, "complex_role": "subunit",
                        "complex_confidence": conf,
                    }
    # Una proteína puede estar en dos firmas (p. ej. lpdA); se conserva la primera
    # y se marca ambigüedad si hay más de una.
    counts: dict[str, set] = {}
    for c in complexes:
        for part in c.subunits.split(";"):
            if ":" in part:
                counts.setdefault(part.split(":", 1)[1], set()).add(c.complex_id)
    for pid, ids in counts.items():
        if len(ids) > 1:
            by_protein.setdefault(pid, {})["complex_status"] = "complex_ambiguous"
    return {
        "genome": genome,
        "complexes": [asdict(c) for c in complexes],
        "by_protein": by_protein,
        "n_complete": sum(c.status == "complex_complete" for c in complexes),
        "n_incomplete": sum(c.status == "complex_incomplete" for c in complexes),
        "n_candidate": 0,
        "n_ambiguous": sum(1 for x in by_protein.values()
                           if x.get("complex_status") == "complex_ambiguous"),
        "max_gap": max_gap,
    }


def write(result: dict[str, Any], outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{result.get('genome', 'genome')}.complexes.tsv"
    rows = result.get("complexes") or []
    cols = ["complex_id", "name", "contig", "start", "end", "subunits",
            "expected", "missing", "status", "evidence", "confidence",
            "source", "limitation"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    return path


def enrich_rows(rows: list[dict], result: dict[str, Any]) -> list[dict]:
    """Añade columnas de complejo a protein_paths.tsv."""
    by = result.get("by_protein") or {}
    for row in rows:
        c = by.get(row.get("protein_id", ""), {})
        row["complex_id"] = c.get("complex_id", "")
        row["complex_role"] = c.get("complex_role", "none")
        row["complex_status"] = c.get("complex_status", "not_evaluated")
        row["complex_confidence"] = c.get("complex_confidence", "")
    return rows


def enrich_gaps(gaps: list[dict], result: dict[str, Any],
                data_root: Path) -> list[dict]:
    """Añade complejo a brechas cuando una EC del complejo produce la reacción.

    Es una asociación por EC/Rhea; si no hay EC, deja ``not_evaluated``.
    """
    by = result.get("by_protein") or {}
    ec_by_complex: dict[str, set[str]] = {}
    # cargar hits desde el resultado de functional si está disponible
    rows = result.get("protein_rows") or []
    for row in rows:
        cid = row.get("complex_id")
        ec = row.get("ec")
        if cid and ec:
            ec_by_complex.setdefault(cid, set()).add(str(ec))
    r2e: dict[str, set[str]] = {}
    p = Path(data_root) / "reactions" / "tsv" / "rhea2ec.tsv" if data_root else None
    if p and p.exists():
        with p.open(encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[3]:
                    r2e.setdefault(f"RHEA:{c[0]}", set()).add(c[3])
    for g in gaps:
        ecs = set()
        for r in str(g.get("producers", "")).split(";"):
            ecs |= r2e.get(r.strip(), set())
        hit = next((cid for cid, es in ec_by_complex.items() if ecs & es), "")
        if hit:
            name = next((c.get("name", "") for c in result.get("complexes", [])
                         if c.get("complex_id") == hit), "")
            g["complex_id"] = hit
            g["complex_name"] = name
            g["complex_status"] = "complex_associated"
            g["complex_role"] = "subunit_reaction"
            g["complex_evidence"] = "EC/Rhea + vecindad"
            g["complex_confidence"] = "medium"
        else:
            g.setdefault("complex_id", "")
            g.setdefault("complex_name", "")
            g.setdefault("complex_status", "not_evaluated")
            g.setdefault("complex_role", "none")
            g.setdefault("complex_evidence", "")
            g.setdefault("complex_confidence", "")
    return gaps
