"""Estados de decisión con evidencia, comunes a todos los pasos.

Complementa a `labels.py`: la etiqueta dice *qué es*; el estado dice *cómo se
decidió* y *qué falta*. Toda salida debería poder responder:

  - ¿qué estado tiene? (validated / consensus / lineage_specific / anomalous /
    unresolved / novel)
  - ¿quién lo valida? (validated_by)
  - ¿por qué se marca así? (reason)
  - ¿qué se intentó? (attempted)
"""
from __future__ import annotations

from dataclasses import dataclass, field

STATUSES: dict[str, str] = {
    "validated": "≥2 métodos independientes coinciden.",
    "consensus": "Los métodos coinciden, sin verdad externa.",
    "lineage_specific": "Desviación esperada del linaje (p. ej. sin dnaA).",
    "anomalous": "Los métodos discrepan; se declara el motivo.",
    "unresolved": "Ningún método dio candidato; se listan los intentos.",
    "novel": "Sin referencia conocida; candidato exploratorio.",
}
CONFIDENCE = ("high", "medium", "low")


@dataclass(frozen=True)
class Decision:
    status: str
    validated_by: str = ""
    reason: str = ""
    attempted: str = ""
    confidence: str | None = None
    evidence: str = ""

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError(f"estado no registrado: {self.status}")
        if self.confidence is not None and self.confidence not in CONFIDENCE:
            raise ValueError(f"confianza inválida: {self.confidence}")


def decide(status: str, *, validated_by: str = "", reason: str = "",
           attempted: str = "", confidence: str | None = None,
           evidence: str = "") -> Decision:
    return Decision(status=status, validated_by=validated_by, reason=reason,
                    attempted=attempted, confidence=confidence, evidence=evidence)


def apply_status(row: dict, decision: Decision) -> dict:
    """Añade los campos de estado a una fila de salida."""
    row.update({
        "status": decision.status,
        "status_validated_by": decision.validated_by,
        "status_reason": decision.reason,
        "status_attempted": decision.attempted,
        "status_confidence": decision.confidence or "",
        "status_evidence": decision.evidence,
    })
    return row
