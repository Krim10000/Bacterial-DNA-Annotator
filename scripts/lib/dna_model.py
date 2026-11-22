"""Definición del modelo de lenguaje de ADN y utilidades de entrenamiento."""
import torch
import torch.nn as nn
import torch.nn.functional as F
import math

class RoPE(nn.Module):
    def __init__(self, dim, max_seq_len=2048):
        super().__init__()
        inv_freq = 1.0 / (10000 ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer("inv_freq", inv_freq)
        self.max_seq_len = max_seq_len

    def forward(self, x, seq_len):
        t = torch.arange(seq_len, device=x.device, dtype=self.inv_freq.dtype)
        freqs = torch.einsum("i,j->ij", t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        return emb.cos()[None, None, :, :], emb.sin()[None, None, :, :]

def apply_rotary_pos_emb(q, k, cos, sin):
    def rotate_half(x):
        x1, x2 = x[..., :x.shape[-1]//2], x[..., x.shape[-1]//2:]
        return torch.cat((-x2, x1), dim=-1)
    q_embed = (q * cos) + (rotate_half(q) * sin)
    k_embed = (k * cos) + (rotate_half(k) * sin)
    return q_embed, k_embed

class RoPETransformerBlock(nn.Module):
    def __init__(self, d_model, n_heads, causal=False):
        super().__init__()
        self.n_heads = n_heads
        self.d_head = d_model // n_heads
        self.causal = causal
        self.qkv = nn.Linear(d_model, 3 * d_model, bias=False)
        self.out = nn.Linear(d_model, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.mlp = nn.Sequential(
            nn.Linear(d_model, 4 * d_model),
            nn.GELU(),
            nn.Linear(4 * d_model, d_model)
        )

    def forward(self, x, cos, sin):
        B, T, C = x.shape
        x_norm = self.norm1(x)
        q, k, v = self.qkv(x_norm).chunk(3, dim=-1)
        q = q.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        k = k.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        v = v.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        
        q, k = apply_rotary_pos_emb(q, k, cos, sin)
        
        out = F.scaled_dot_product_attention(q, k, v, is_causal=self.causal)
        out = out.transpose(1, 2).contiguous().view(B, T, C)
        x = x + self.out(out)
        x = x + self.mlp(self.norm2(x))
        return x

class DNAModel(nn.Module):
    def __init__(self, vocab_size, d_model=768, n_layers=12, n_heads=12, max_len=2048, causal=False):
        super().__init__()
        self.causal = causal
        self.embed = nn.Embedding(vocab_size, d_model)
        self.rope = RoPE(d_model // n_heads, max_len)
        self.blocks = nn.ModuleList([RoPETransformerBlock(d_model, n_heads, causal) for _ in range(n_layers)])
        self.norm = nn.LayerNorm(d_model)
        self.lm_head = nn.Linear(d_model, vocab_size, bias=False)
        self.lm_head.weight = self.embed.weight

    def forward(self, input_ids):
        x = self.embed(input_ids)
        cos, sin = self.rope(x, x.shape[1])
        for block in self.blocks:
            x = block(x, cos, sin)
        x = self.norm(x)
        return self.lm_head(x)
