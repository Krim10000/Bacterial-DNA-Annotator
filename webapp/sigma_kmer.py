"""Calibrador k-mer de sigma, específico de especie y serializable.

El modelo lo entrena ``scripts/sigma/63_train_kmer_calibrator.py`` sobre el
gold de una especie. No es universal: el JSON guarda genoma, ventana, n-gramo,
vocabulario, escala y coeficientes. Si falta, no se inventan probabilidades.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


@dataclass(frozen=True)
class KmerModel:
    genome: str
    labels: tuple[str, ...]
    vocab: dict[str, int]
    scale: tuple[float, ...]
    coef: tuple[tuple[float, ...], ...]
    intercept: tuple[float, ...]
    ngram: int = 4
    window_up: int = 40
    window_down: int = 5
    trained_on: str = ""
    schema: str = "sigma-kmer/v1"

    def counts(self, sequence: str) -> list[float]:
        """Cuenta n-gramos con el mismo tokenizador del entrenamiento."""
        seq = (sequence or "").lower()
        n = int(self.ngram)
        out = [0.0] * len(self.vocab)
        for i in range(max(0, len(seq) - n + 1)):
            idx = self.vocab.get(seq[i:i + n])
            if idx is not None:
                out[idx] += 1.0
        return out

    def probabilities(self, sequence: str) -> dict[str, float]:
        """Probabilidades one-vs-rest; no se normalizan entre sigma."""
        x = self.counts(sequence)
        if len(x) != len(self.scale):
            raise ValueError("vocabulario/escala incompatibles")
        x = [v / s if s else v for v, s in zip(x, self.scale)]
        out: dict[str, float] = {}
        for label, coef, b in zip(self.labels, self.coef, self.intercept):
            out[label] = round(_sigmoid(sum(a * c for a, c in zip(x, coef)) + b), 6)
        return out

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "genome": self.genome,
            "labels": list(self.labels),
            "vocab": dict(self.vocab),
            "scale": list(self.scale),
            "coef": [list(row) for row in self.coef],
            "intercept": list(self.intercept),
            "ngram": self.ngram,
            "window_up": self.window_up,
            "window_down": self.window_down,
            "trained_on": self.trained_on,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "KmerModel":
        if not isinstance(data, Mapping) or data.get("schema") != "sigma-kmer/v1":
            raise ValueError("esquema de calibrador k-mer no soportado")
        try:
            labels = tuple(str(x) for x in data["labels"])
            vocab = {str(k): int(v) for k, v in data["vocab"].items()}
            scale = tuple(float(x) for x in data["scale"])
            coef = tuple(tuple(float(x) for x in row) for row in data["coef"])
            intercept = tuple(float(x) for x in data["intercept"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"calibrador k-mer incompleto: {exc}") from exc
        if len(labels) != len(coef) or len(labels) != len(intercept):
            raise ValueError("labels/coef/intercept no coinciden")
        if any(len(row) != len(scale) for row in coef):
            raise ValueError("coef y escala no coinciden")
        return cls(
            genome=str(data.get("genome") or ""), labels=labels, vocab=vocab,
            scale=scale, coef=coef, intercept=intercept,
            ngram=int(data.get("ngram") or 4),
            window_up=int(data.get("window_up") or 40),
            window_down=int(data.get("window_down") or 5),
            trained_on=str(data.get("trained_on") or ""),
        )


def load_model(path: str | Path) -> KmerModel:
    p = Path(path)
    return KmerModel.from_dict(json.loads(p.read_text(encoding="utf-8")))


def save_model(model: KmerModel, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(model.as_dict(), ensure_ascii=False), encoding="utf-8")


__all__ = ["KmerModel", "load_model", "save_model"]
