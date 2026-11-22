#!/usr/bin/env python3
"""Interfaz única para el modelo de ADN v2 (RoPE + ventana local).

Motivación: la arquitectura estaba copiada en 09_train_model_v2.py,
12_eval_probe_v2.py, 13_comprehensive_eval.py y 14_robust_eval.py. Eso hacía
imposible garantizar que todas las evaluaciones usaran el mismo modelo y el
mismo tokenizador.

Este módulo:
  - carga el checkpoint y el tokenizador una sola vez;
  - expone `encode()` (estados ocultos) y `forward()` (logits);
  - añade `key_padding_mask` para que el padding no contamine la atención
    (el entrenamiento original no lo tenía; para análisis es más limpio);
  - extrae embeddings por ventana en lotes, con pooling sobre tokens reales.

No modifica los pesos ni el checkpoint.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    import sentencepiece as spm
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Falta sentencepiece en el entorno (.venv).") from exc


# ---------------------------------------------------------------------------
# Arquitectura (debe coincidir exactamente con 09_train_model_v2.py)
# ---------------------------------------------------------------------------
class RoPE(nn.Module):
    def __init__(self, dim: int, max_seq_len: int = 32000):
        super().__init__()
        inv_freq = 1.0 / (10000 ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer("inv_freq", inv_freq)
        self.max_seq_len = max_seq_len

    def forward(self, x: torch.Tensor, seq_len: int):
        t = torch.arange(seq_len, device=x.device, dtype=self.inv_freq.dtype)
        freqs = torch.einsum("i,j->ij", t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        return emb.cos()[None, None, :, :], emb.sin()[None, None, :, :]


def apply_rotary_pos_emb(q, k, cos, sin):
    def rotate_half(x):
        x1, x2 = x[..., : x.shape[-1] // 2], x[..., x.shape[-1] // 2 :]
        return torch.cat((-x2, x1), dim=-1)

    q_embed = (q * cos) + (rotate_half(q) * sin)
    k_embed = (k * cos) + (rotate_half(k) * sin)
    return q_embed, k_embed


class RoPETransformerBlock(nn.Module):
    def __init__(self, d_model: int, n_heads: int, window_size: int = 4096):
        super().__init__()
        self.n_heads = n_heads
        self.d_head = d_model // n_heads
        self.window_size = window_size
        self.qkv = nn.Linear(d_model, 3 * d_model, bias=False)
        self.out = nn.Linear(d_model, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.mlp = nn.Sequential(
            nn.Linear(d_model, 4 * d_model),
            nn.GELU(),
            nn.Linear(4 * d_model, d_model),
        )

    def forward(self, x, cos, sin, key_padding_mask=None):
        B, T, C = x.shape
        x_norm = self.norm1(x)
        q, k, v = self.qkv(x_norm).chunk(3, dim=-1)
        q = q.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        k = k.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        v = v.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        q, k = apply_rotary_pos_emb(q, k, cos, sin)

        # Banda local (el entrenamiento original usaba esta misma máscara).
        mask = torch.ones(T, T, device=x.device, dtype=torch.bool)
        mask = mask.triu(diagonal=-self.window_size).tril(diagonal=self.window_size)

        # key_padding_mask: [B, T], True = posición de padding.
        # SDPA interpreta True como "puede atender".
        if key_padding_mask is not None and bool(key_padding_mask.any()):
            keep = ~key_padding_mask
            # [B, 1, T, T] para emitir con [B, H, T, T] de q/k.
            mask = mask.unsqueeze(0).unsqueeze(0) & keep.unsqueeze(1).unsqueeze(2)

        out = F.scaled_dot_product_attention(q, k, v, attn_mask=mask)
        out = out.transpose(1, 2).contiguous().view(B, T, C)
        x = x + self.out(out)
        x = x + self.mlp(self.norm2(x))
        return x


class DNAModelV2(nn.Module):
    def __init__(self, vocab_size: int, d_model: int = 384, n_layers: int = 6,
                 n_heads: int = 6, max_len: int = 32000):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, d_model)
        self.rope = RoPE(d_model // n_heads, max_len)
        self.blocks = nn.ModuleList(
            [RoPETransformerBlock(d_model, n_heads) for _ in range(n_layers)]
        )
        self.norm = nn.LayerNorm(d_model)
        self.lm_head = nn.Linear(d_model, vocab_size, bias=False)
        self.lm_head.weight = self.embed.weight
        self.gradient_checkpointing = False

        self.apply(self._init_weights)

    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            module.weight.data.normal_(mean=0.0, std=0.02)
            if module.bias is not None:
                module.bias.data.zero_()
        elif isinstance(module, nn.Embedding):
            module.weight.data.normal_(mean=0.0, std=0.02)

    def encode(self, input_ids: torch.Tensor, key_padding_mask=None):
        x = self.embed(input_ids)
        cos, sin = self.rope(x, x.shape[1])
        for block in self.blocks:
            x = block(x, cos, sin, key_padding_mask=key_padding_mask)
        return self.norm(x)

    def forward(self, input_ids: torch.Tensor, labels=None, key_padding_mask=None):
        x = self.encode(input_ids, key_padding_mask=key_padding_mask)
        logits = self.lm_head(x)
        if labels is not None:
            loss = F.cross_entropy(
                logits.view(-1, logits.size(-1)), labels.view(-1), ignore_index=-100
            )
            return loss, logits
        return logits


# ---------------------------------------------------------------------------
# Carga y tokenización
# ---------------------------------------------------------------------------
def load_v2(checkpoint, tokenizer_path, device: str = "cpu"):
    """Carga checkpoint + tokenizador y valida que las dimensiones coincidan."""
    checkpoint = Path(checkpoint)
    tokenizer_path = Path(tokenizer_path)
    if not checkpoint.exists():
        raise FileNotFoundError(f"No existe el checkpoint: {checkpoint}")
    if not tokenizer_path.exists():
        raise FileNotFoundError(f"No existe el tokenizador: {tokenizer_path}")

    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    state = torch.load(checkpoint, map_location="cpu")
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]

    model = DNAModelV2(vocab_size=sp.get_piece_size())
    model.load_state_dict(state, strict=True)
    model.to(torch.device(device)).eval()
    return model, sp


def tokenize(sp, seq: str, max_tokens: int | None = None) -> list[int]:
    """Tokeniza ADN. Devuelve [] para secuencias vacías."""
    if not seq:
        return []
    ids = sp.encode(seq, out_type=int)
    if max_tokens is not None:
        ids = ids[:max_tokens]
    return ids


def _pool(h: torch.Tensor, lengths: torch.Tensor, pool: str) -> torch.Tensor:
    rows = []
    for i, n in enumerate(lengths.tolist()):
        if n <= 0:
            rows.append(torch.zeros(h.shape[-1], device=h.device))
            continue
        v = h[i, :n]
        rows.append(v.mean(dim=0) if pool == "mean" else v.max(dim=0).values)
    return torch.stack(rows)


@torch.inference_mode()
def embed_sequences(model, sp, sequences: Sequence[str], batch_size: int = 32,
                    device: str = "cpu", max_tokens: int = 8192,
                    pool: str = "mean", show_progress: bool = False):
    """Extrae un vector [d_model] por secuencia (media de tokens reales).

    Devuelve (embeddings [N, d_model] float32, n_tokens [N] int32).
    Agrupa por longitud para minimizar padding y enmascara el padding en la
    atención, de modo que los vectores no dependan de los tokens de relleno.
    """
    if pool not in {"mean", "max"}:
        raise ValueError("pool debe ser 'mean' o 'max'")

    tokenized = [tokenize(sp, s, max_tokens=max_tokens) for s in sequences]
    lengths_all = np.array([len(t) for t in tokenized], dtype=np.int32)
    order = np.argsort(lengths_all, kind="stable")
    out = np.zeros((len(sequences), model.embed.embedding_dim), dtype=np.float32)
    dev = torch.device(device)

    n_batches = math.ceil(len(order) / batch_size)
    for b in range(n_batches):
        idx = order[b * batch_size : (b + 1) * batch_size]
        batch = [tokenized[i] for i in idx]
        max_len = max(len(t) for t in batch)
        if max_len == 0:
            continue

        input_ids = torch.zeros((len(batch), max_len), dtype=torch.long, device=dev)
        pad_mask = torch.ones((len(batch), max_len), dtype=torch.bool, device=dev)
        for j, ids in enumerate(batch):
            if ids:
                input_ids[j, : len(ids)] = torch.tensor(ids, dtype=torch.long, device=dev)
                pad_mask[j, : len(ids)] = False

        h = model.encode(input_ids, key_padding_mask=pad_mask)
        pooled = _pool(h, torch.tensor([len(t) for t in batch], device=dev), pool)
        out[idx] = pooled.float().cpu().numpy()

        if show_progress and (b + 1) % 10 == 0:
            print(f"  embeddings: lote {b + 1}/{n_batches} "
                  f"({min((b + 1) * batch_size, len(order))}/{len(order)} secuencias)")

    return out, lengths_all


@torch.inference_mode()
def sequence_perplexity(model, sp, seq: str, max_tokens: int = 8192,
                        device: str = "cpu") -> float:
    """Perplejidad de una secuencia sin padding (diagnóstico)."""
    ids = tokenize(sp, seq, max_tokens=max_tokens)
    if len(ids) < 2:
        return float("nan")
    dev = torch.device(device)
    x = torch.tensor([ids], dtype=torch.long, device=dev)
    logits = model(x)
    loss = F.cross_entropy(logits[0, :-1], x[0, 1:])
    return float(torch.exp(loss).item())


# ---------------------------------------------------------------------------
# CLI de autocomprobación
# ---------------------------------------------------------------------------
def _selftest():
    import argparse

    ap = argparse.ArgumentParser(description="Autocomprobación de dna_v2")
    ap.add_argument("--checkpoint", default="dna_model_v2_epoch9.pt")
    ap.add_argument("--tokenizer", default="data/tokenizer_v2/dna_sp_v2.model")
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    model, sp = load_v2(args.checkpoint, args.tokenizer, device=args.device)
    seqs = [
        "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA",
        "TTGACATATAATAGGAGGATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT",
        "NNNNNNNNNNACGTACGTNNNNNNNNNN",
    ]
    emb, lengths = embed_sequences(model, sp, seqs, batch_size=2, device=args.device)
    assert emb.shape == (len(seqs), 384), emb.shape
    assert np.isfinite(emb).all(), "embeddings no finitos"
    assert lengths.tolist() == [len(tokenize(sp, s)) for s in seqs]
    ppl = [sequence_perplexity(model, sp, s, device=args.device) for s in seqs]
    print(f"OK: embeddings {emb.shape}, dim={emb.shape[1]}, "
          f"tokens={lengths.tolist()}, ppl={[round(p, 2) for p in ppl]}")
    print(f"checkpoint={args.checkpoint} tokenizer_vocab={sp.get_piece_size()}")


if __name__ == "__main__":
    _selftest()
