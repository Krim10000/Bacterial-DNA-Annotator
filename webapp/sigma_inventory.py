"""Inventario de factores sigma por especie, sin panel fijo.

No traduce nombres nativos a una lista global de seis sigma. Cada especie
declara sus proteínas candidatas, su nombre nativo y una familia amplia
(``sigma70_like``, ``sigma54_like``, ``ecf_like``, ``unknown``). La familia es
una pista de dominio/anotación, no una equivalencia funcional exacta.
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable, Mapping

# Nombres que NO son factores sigma aunque contengan "sigma".
_ANTI = re.compile(r"anti[- ]?sigma|regulator of sigma|sigma[- ]?dependent|"
                   r"\brsd\b|\brse[ab]\b|\bflgm\b", re.I)
# Subunidades del core RNAP: no son sigma.
_CORE = re.compile(r"rna polymerase (subunit|omega|delta)|"
                   r"\brpo[abczy]\b", re.I)
# Candidatos por nombre nativo; cada especie puede tener otros.
_NAME = re.compile(r"^(sig[a-z0-9]+|rpo[dehns]|flia|feci)$", re.I)


@dataclass(frozen=True)
class SigmaFactor:
    species: str
    native_name: str
    protein_id: str = ""
    locus_tag: str = ""
    product: str = ""
    family: str = "unknown"
    evidence: str = "anotación"
    confidence: str = "medium"

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify_family(native_name: str, product: str = "") -> str:
    """Familia amplia, sin asignar un subtipo funcional inventado."""
    text = f"{native_name} {product}".lower()
    if "sigma54" in text or "sigma 54" in text or re.search(r"\b(rpon|sigl)\b", text):
        return "sigma54_like"
    if "ecf" in text or re.search(r"\b(rpoe|feci|sigv|sigw|sigx|sigy)\b", text):
        return "ecf_like"
    if "sigma" in text or re.match(r"^sig", native_name, re.I):
        return "sigma70_like"
    return "unknown"


def is_sigma_candidate(name: str, product: str) -> bool:
    """Filtro conservador: producto/nombre de sigma, sin anti-sigma ni core."""
    name = (name or "").strip()
    product = (product or "").strip()
    text = f"{name} {product}"
    if not name and not product:
        return False
    if _ANTI.search(text) or _CORE.search(text):
        return False
    if _NAME.match(name):
        return True
    return bool(re.search(r"rna polymerase sigma factor", product, re.I))


def inventory_from_rows(species: str,
                        rows: Iterable[Mapping[str, Any]]) -> list[SigmaFactor]:
    """Construye el inventario desde filas de genes; deduplica por nombre."""
    out: dict[str, SigmaFactor] = {}
    for row in rows:
        name = str(row.get("gene_name") or row.get("name") or "").strip()
        product = str(row.get("product") or "").strip()
        if not is_sigma_candidate(name, product):
            continue
        key = (name or str(row.get("locus_tag") or row.get("gene_id") or "")).lower()
        if not key or key in out:
            continue
        out[key] = SigmaFactor(
            species=species,
            native_name=name or key,
            protein_id=str(row.get("protein_id") or row.get("acc") or ""),
            locus_tag=str(row.get("locus_tag") or row.get("gene_id") or ""),
            product=product,
            family=classify_family(name, product),
            evidence="producto/nombre anotado",
            confidence="medium",
        )
    return sorted(out.values(), key=lambda f: f.native_name.lower())


def load_inventory(data_root: Path | str, genome: str) -> list[SigmaFactor]:
    """Lee `data/regulation/<genome>/genes.tsv`; no descarga nada."""
    path = Path(data_root) / "regulation" / genome / "genes.tsv"
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return inventory_from_rows(genome, rows)


def summary(factors: Iterable[SigmaFactor]) -> dict[str, Any]:
    items = list(factors)
    counts: dict[str, int] = {}
    for item in items:
        counts[item.family] = counts.get(item.family, 0) + 1
    return {
        "n_factors": len(items),
        "n_families": len(counts),
        "families": counts,
        "names": [item.native_name for item in items],
        "factors": [item.as_dict() for item in items],
    }


def to_tsv(factors: Iterable[SigmaFactor]) -> str:
    cols = ["species", "native_name", "protein_id", "locus_tag", "product",
            "family", "evidence", "confidence"]
    lines = ["\t".join(cols)]
    for item in factors:
        lines.append("\t".join(str(getattr(item, c)) for c in cols))
    return "\n".join(lines) + "\n"
