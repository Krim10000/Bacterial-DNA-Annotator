"""Esquema de etiquetado con procedencia explícita.

Toda fila de salida debe poder responder:
  - ¿qué etiqueta tiene? (label)
  - ¿de qué fuente sale? (source)
  - ¿con qué evidencia? (evidence)
  - ¿cómo se interpreta y con qué confianza? (interpretation, confidence)

No se permite una etiqueta sin fuente. Este módulo centraliza las fuentes,
las etiquetas válidas y la construcción de filas.
"""
from __future__ import annotations

from datetime import datetime, timezone

PIPELINE = "09_ghost 0.1.0"

SOURCES: dict[str, str] = {
    "orf": "ORF de novo, 6 marcos, motor 09_ghost 0.1.0",
    "pyrodigal": "pyrodigal 3.7.1 (predicción de CDS)",
    "infernal": "Infernal/cmsearch (rRNA)",
    "trnascan": "tRNAscan-SE 2.0.9 (tRNA)",
    "orca": "ORCA + DoriC (oriC)",
    "mobile": "mobile_elements.json (transposones)",
    "defense": "defense_systems.json (defensa/profagos)",
    "rna": "rna_annotations.json (ARN)",
    "diamond_swissprot": "DIAMOND 2.1.11 vs Swiss-Prot (data/db/swissprot.dmnd)",
    "pfam": "Pfam-A 2026 (hmmscan/hmmsearch, corte GA)",
    "rhea": "Rhea (rhea.rdf.gz) + rhea2uniprot_sprot.tsv",
    "esm": "ESM-2 8M (CPU, log-verosimilitud media)",
    "conservation": "Conservación entre genomas (diamond cluster 95%)",
    "phylum": "Conservación de reacción en el mismo filo (species_phylum.tsv)",
}

LABELS: dict[str, str] = {
    "annotated_known": "Gen anotado con función específica (homología/dominio).",
    "predicted_hypothetical": "Gen predicho sin función específica asignada.",
    "conserved_hypothetical": "Proteína 'hypothetical' presente en ≥2 genomas.",
    "candidate_novel": "ORF no anotado con evidencia de codificación.",
    "candidate_conserved": "ORF no anotado presente en ≥2 genomas.",
    "desert_region": "Región ≥3 kb sin ninguna señal de anotación.",
    "desert_candidate": "Desierto con ORFs conservados: candidato a gen/operón.",
    "gap_biosynthetic": "Brecha refinada: paso enzimático, no transporte, con demanda.",
    "gap_closure_candidate": "Candidato que podría cerrar una brecha (evidencia de dominio/reacción).",
    "gap_closure_confirmed": "Cierre validado tras verificar dirección y reacción.",
    "discarded_artifact": "Brecha descartada (ión, transporte, ROS o sin demanda).",
}

CONFIDENCE = ("high", "medium", "low")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def label_row(
    row: dict,
    *,
    label: str,
    source: str,
    evidence: str,
    confidence: str,
    interpretation: str,
) -> dict:
    """Añade los campos de procedencia a una fila."""
    if label not in LABELS:
        raise ValueError(f"etiqueta no registrada: {label}")
    if confidence not in CONFIDENCE:
        raise ValueError(f"confianza inválida: {confidence}")
    row.update({
        "label": label,
        "label_source": source,
        "label_evidence": evidence,
        "label_confidence": confidence,
        "label_interpretation": interpretation,
        "label_pipeline": PIPELINE,
        "label_date": now(),
    })
    return row
