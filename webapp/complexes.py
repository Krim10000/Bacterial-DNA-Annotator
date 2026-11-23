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
        "ecs": ["7.1.2.2"], "source": "curated literature; operon atp",
    },
    {
        "id": "NADH_DEHYDROGENASE", "name": "NADH dehydrogenase I (Nuo)",
        "genes": [f"nuo{c}" for c in "ABCDEFGHIJKLMN"],
        "essential": ["nuoB", "nuoC", "nuoD", "nuoE", "nuoF", "nuoG",
                      "nuoH", "nuoI", "nuoJ", "nuoK", "nuoL", "nuoM", "nuoN"],
        "ecs": ["7.1.1.2"], "source": "curated literature; operon nuo",
    },
    {
        "id": "SUCCINATE_DEHYDROGENASE", "name": "Succinate dehydrogenase",
        "genes": ["sdhA", "sdhB", "sdhC", "sdhD"],
        "essential": ["sdhA", "sdhB", "sdhC", "sdhD"],
        "ecs": ["1.3.5.1"], "source": "curated literature; operon sdh",
    },
    {
        "id": "CYTOCHROME_BC1", "name": "Cytochrome bc1 complex",
        "genes": ["petA", "petB", "petC"],
        "essential": ["petA", "petB", "petC"],
        "ecs": ["7.1.1.8"], "source": "curated literature; operon pet",
    },
    {
        "id": "CYTOCHROME_OXIDASE", "name": "Cytochrome c oxidase",
        "genes": ["cyoA", "cyoB", "cyoC", "cyoD"],
        "essential": ["cyoA", "cyoB", "cyoC"],
        "ecs": ["7.1.1.9"], "source": "curated literature; operon cyo",
    },
    {
        "id": "PYRUVATE_DEHYDROGENASE", "name": "Pyruvate dehydrogenase",
        "genes": ["aceE", "aceF", "lpdA"],
        "essential": ["aceE", "aceF", "lpdA"],
        "ecs": ["1.2.4.1", "1.8.1.4"], "source": "curated literature; operon ace",
    },
    {
        "id": "OXOGLUTARATE_DEHYDROGENASE", "name": "2-oxoglutarate dehydrogenase",
        "genes": ["sucA", "sucB", "lpdA"],
        "essential": ["sucA", "sucB", "lpdA"],
        "ecs": ["1.2.4.2", "1.8.1.4"], "source": "curated literature; operon suc",
    },
    {
        "id": "NITROGENASE", "name": "Nitrogenase",
        "genes": ["nifH", "nifD", "nifK"],
        "essential": ["nifH", "nifD", "nifK"],
        "ecs": ["1.18.6.1"], "source": "curated literature; operon nif",
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
    {
        "id": "CYTOCHROME_BD", "name": "Cytochrome bd oxidase",
        "genes": ["cydA", "cydB"], "essential": ["cydA", "cydB"],
        "ecs": ["1.10.3.-"], "source": "curated literature; operon cyd",
    },
    {
        "id": "NITRATE_REDUCTASE", "name": "Nitrate reductase (Nar)",
        "genes": ["narG", "narH", "narI", "narJ"], "essential": ["narG", "narH"],
        "ecs": ["1.7.5.1"], "source": "curated literature; operon nar",
    },
    {
        "id": "NITRITE_REDUCTASE", "name": "Nitrite reductase (Nir)",
        "genes": ["nirB", "nirD"], "essential": ["nirB", "nirD"],
        "ecs": ["1.7.1.15"], "source": "curated literature; operon nir",
    },
    {
        "id": "HYDROGENASE", "name": "Hydrogenase",
        "genes": ["hyaA", "hyaB", "hybA", "hybB", "hydA"],
        "essential": ["hyaA", "hyaB"], "ecs": ["1.12.7.2"],
        "source": "curated literature; hya/hyb/hyd",
    },
    {
        "id": "FORMATE_DEHYDROGENASE", "name": "Formate dehydrogenase",
        "genes": ["fdhF", "fdhD", "fdhE"], "essential": ["fdhF"],
        "ecs": ["1.17.1.9", "1.2.1.2"], "source": "curated literature; fdh",
    },
    {
        "id": "PILI", "name": "Type IV pilus",
        "genes": ["pilA", "pilB", "pilC", "pilD"], "essential": ["pilA", "pilB"],
        "source": "curated literature; pil",
    },
    {
        "id": "FIMBRIAE", "name": "Type 1 fimbriae",
        "genes": ["fimA", "fimC", "fimD", "fimF", "fimG", "fimH"],
        "essential": ["fimA", "fimH"], "source": "curated literature; fim",
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
    ecs: str = ""
    flag: str = "grey"
    limitation: str = (
        "Vecindad génica; no demuestra interacción física ni estequiometría. "
        "No usar como confirmación."
    )


def flag(status: str) -> str:
    """Color de revisión: verde completo, amarillo sospecha, gris sin evaluar."""
    return {"complex_complete": "green",
            "complex_incomplete": "yellow",
            "complex_candidate": "yellow",
            "complex_ambiguous": "orange"}.get(status, "grey")


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
                ecs=";".join(sig.get("ecs", [])),
                flag=flag(status),
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
            "source", "ecs", "flag", "limitation"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    return path


def enrich_rows(rows: list[dict], result: dict[str, Any]) -> list[dict]:
    """Añade columnas de complejo a protein_paths.tsv."""
    by = result.get("by_protein") or {}
    missing = {c["complex_id"]: c.get("missing", "")
               for c in result.get("complexes", [])}
    for row in rows:
        c = by.get(row.get("protein_id", ""), {})
        cid = c.get("complex_id", "")
        row["complex_id"] = cid
        row["complex_name"] = c.get("complex_name", "")
        row["complex_role"] = c.get("complex_role", "none")
        row["complex_status"] = c.get("complex_status", "not_evaluated")
        row["complex_confidence"] = c.get("complex_confidence", "")
        row["complex_missing"] = missing.get(cid, "")
        row["complex_flag"] = flag(row.get("complex_status", ""))
    return rows


def enrich_gaps(gaps: list[dict], result: dict[str, Any],
                data_root: Path) -> list[dict]:
    """Asocia brechas a complejos por UniProt→Rhea o EC→Rhea.

    Si el complejo está incompleto, la brecha se marca como
    ``complex_incomplete``: la explicación es "puede faltar una subunidad",
    no "la enzima no existe". No cierra la brecha.
    """
    comps = {c["complex_id"]: c for c in result.get("complexes", [])}
    rows = result.get("protein_rows") or []
    by_pid = {r.get("protein_id"): r for r in rows}
    rhea_by_complex: dict[str, set[str]] = {}
    ec_by_complex: dict[str, set[str]] = {}
    for cid, c in comps.items():
        rhea, ec = set(), set()
        for part in str(c.get("subunits", "")).split(";"):
            if ":" not in part:
                continue
            pid = part.split(":", 1)[1]
            row = by_pid.get(pid, {})
            acc, e = str(row.get("uniprot", "")), str(row.get("ec", ""))
            if e:
                ec.add(e)
            if acc:
                # se completa después con la tabla Rhea→UniProt
                rhea.add(acc)
        rhea_by_complex[cid] = rhea
        ec_by_complex[cid] = ec | {x.strip() for x in str(c.get("ecs", "")).split(";") if x.strip()}
    # Rhea -> UniProt y Rhea -> EC
    acc_rhea: dict[str, set[str]] = {}
    p = Path(data_root) / "metabolism" / "rhea2uniprot_sprot.tsv" if data_root else None
    if p and p.exists():
        with p.open(encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[3]:
                    acc_rhea.setdefault(c[3], set()).add(f"RHEA:{c[0]}")
    r2e: dict[str, set[str]] = {}
    p = Path(data_root) / "reactions" / "tsv" / "rhea2ec.tsv" if data_root else None
    if p and p.exists():
        with p.open(encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[3]:
                    r2e.setdefault(f"RHEA:{c[0]}", set()).add(c[3])
    for cid in list(rhea_by_complex):
        rhea_by_complex[cid] = {r for acc in rhea_by_complex[cid]
                                for r in acc_rhea.get(acc, set())}
    for g in gaps:
        producers = {r.strip() for r in str(g.get("producers", "")).split(";") if r.strip()}
        ecs = {e for r in producers for e in r2e.get(r, set())}
        matches = [cid for cid in comps
                   if producers & rhea_by_complex.get(cid, set())
                   or ecs & ec_by_complex.get(cid, set())]
        if not matches:
            g.setdefault("complex_id", "")
            g.setdefault("complex_name", "")
            g.setdefault("complex_status", "not_evaluated")
            g.setdefault("complex_role", "none")
            g.setdefault("complex_evidence", "")
            g.setdefault("complex_confidence", "")
            g.setdefault("complex_missing", "")
            continue
        cid = matches[0]
        c = comps[cid]
        status = c.get("status", "complex_candidate")
        if len(matches) > 1:
            status = "complex_ambiguous"
        g["complex_id"] = cid
        g["complex_name"] = c.get("name", "")
        g["complex_status"] = status
        g["complex_role"] = "subunit_reaction"
        g["complex_evidence"] = "UniProt/Rhea + vecindad"
        g["complex_confidence"] = c.get("confidence", "low")
        g["complex_missing"] = c.get("missing", "")
        g["complex_flag"] = flag(status)
        if status != "complex_complete":
            g["complex_note"] = (
                "Se sospecha un complejo incompleto: no cierra la brecha. "
                "Faltan: " + (c.get("missing", "") or "subunidad no identificada"))
    return gaps
