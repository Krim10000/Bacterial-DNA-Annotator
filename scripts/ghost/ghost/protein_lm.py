"""Modelo de proteína ligero (ESM-2 8M) para filtrar péptidos no plegables.

Puntúa cada cadena con la log-verosimilitud media por residuo bajo ESM-2.
Las proteínas reales tienen secuencias mucho más probables que péptidos
aleatorios o composicionalmente sesgados.
"""
from __future__ import annotations

import math

_MODEL = None
_ALPHABET = None
_BC = None


def _load():
    global _MODEL, _ALPHABET, _BC
    if _MODEL is None:
        import esm  # type: ignore
        import torch

        _MODEL, _ALPHABET = esm.pretrained.esm2_t6_8M_UR50D()
        _MODEL.eval()
        _BC = _ALPHABET.get_batch_converter()
    return _MODEL, _ALPHABET, _BC


_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean(seq: str, max_len: int = 1022) -> str:
    """Normaliza para ESM: corta en la parada y sustituye símbolos raros."""
    seq = (seq or "").upper().split("*")[0]
    return "".join(c if c in _AA else "X" for c in seq)[:max_len]


def score_sequences(seqs: list[str], batch_size: int = 8) -> list[float]:
    """Log-verosimilitud media por residuo (mayor = más proteína-like)."""
    import torch

    model, alphabet, bc = _load()
    out_scores: list[float] = []
    for i in range(0, len(seqs), batch_size):
        chunk = [clean(x) for x in seqs[i:i + batch_size]]
        labels, strs, toks = bc([(str(j), s) for j, s in enumerate(chunk)])
        with torch.no_grad():
            logits = model(toks, repr_layers=[], return_contacts=False)["logits"]
            logprobs = torch.log_softmax(logits, dim=-1)
        for j, s in enumerate(strs):
            n = min(len(s), toks.shape[1] - 2)
            if n <= 0:
                out_scores.append(float("-inf"))
                continue
            ids = toks[j, 1:n + 1]          # 0 = BOS
            lp = logprobs[j, 1:n + 1, :]
            vals = lp.gather(1, ids.unsqueeze(-1)).squeeze(-1)
            out_scores.append(float(vals.mean().item()))
    return out_scores


def score_one(seq: str) -> float:
    return score_sequences([seq])[0]
