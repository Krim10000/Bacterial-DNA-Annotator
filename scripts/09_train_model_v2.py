#!/usr/bin/env python3
"""Entrenamiento del Modelo v2 con Sliding Window y Gradient Checkpointing."""
import math, argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sentencepiece as spm
from pathlib import Path
from torch.utils.data import Dataset, DataLoader
from torch.utils.checkpoint import checkpoint
from tqdm import tqdm
import sys

sys.path.insert(0, str(Path(__file__).parent))
from lib.io import read_config

# ==============================================================================
# 1. DATASET (Manejo de secuencias largas de 32k y soporte IUPAC)
# ==============================================================================
class CircularDNADataset(Dataset):
    def __init__(self, npz_paths, tokenizer_path, max_len=32000, mask_prob=0.15):
        self.npz_paths = npz_paths
        self.sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
        self.max_len = max_len
        self.mask_prob = mask_prob

    def __len__(self):
        return len(self.npz_paths)

    def __getitem__(self, idx):
        data = np.load(self.npz_paths[idx])
        seq = data['nt'] 
        
        # Mapa de caracteres para decodificar de vuelta a string
        iupac_chars = "ACGTRYKMSWBDHVN"
        
        # Rotación circular
        if len(seq) > self.max_len:
            offset = np.random.randint(0, len(seq))
            seq = np.take(seq, np.arange(offset, offset + self.max_len), mode='wrap')
        else:
            seq = np.pad(seq, (0, self.max_len - len(seq)), constant_values=14) # 14 es N

        # Reverse Complement con soporte para IUPAC
        if np.random.rand() < 0.5:
            seq = seq[::-1] # Invertir
            # Reverso complementario biológico:
            # A(0)<->T(3), C(1)<->G(2), R(4)<->Y(5), K(6)<->M(7), S(8)<->S(8), W(9)<->W(9), B(10)<->V(13), D(11)<->H(12), N(14)<->N(14)
            rc_map = np.array([3, 2, 1, 0, 5, 4, 7, 6, 8, 9, 13, 12, 11, 10, 14], dtype=np.uint8)
            seq = rc_map[seq]

        # Convertir a texto usando el mapa IUPAC
        seq_str = "".join([iupac_chars[i] for i in seq])
        
        tokens = self.sp.encode(seq_str, out_type=int)[:self.max_len]
        
        if len(tokens) < self.max_len:
            tokens = tokens + [0] * (self.max_len - len(tokens))

        input_ids = torch.tensor(tokens, dtype=torch.long)
        labels = input_ids.clone()

        # Enmascaramiento MLM (15%)
        prob = torch.rand(input_ids.shape)
        mask = (prob < self.mask_prob) & (input_ids != 0)
        input_ids[mask] = 1
        labels[~mask] = -100 
            
        return {"input_ids": input_ids, "labels": labels}

# ==============================================================================
# 2. ARQUITECTURA (RoPE + Sliding Window + Checkpointing)
# ==============================================================================
class RoPE(nn.Module):
    def __init__(self, dim, max_seq_len=32000):
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
    def __init__(self, d_model, n_heads, window_size=4096):
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
        
        # Sliding Window Attention Mask (Máscara de banda)
        mask = torch.ones(T, T, device=x.device, dtype=torch.bool)
        mask = mask.triu(diagonal=-self.window_size).tril(diagonal=self.window_size)
        
        # FlashAttention con máscara
        out = F.scaled_dot_product_attention(q, k, v, attn_mask=mask)
        
        out = out.transpose(1, 2).contiguous().view(B, T, C)
        x = x + self.out(out)
        x = x + self.mlp(self.norm2(x))
        return x

class DNAModelV2(nn.Module):
    def __init__(self, vocab_size, d_model=384, n_layers=6, n_heads=6, max_len=32000):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, d_model)
        self.rope = RoPE(d_model // n_heads, max_len)
        self.blocks = nn.ModuleList([RoPETransformerBlock(d_model, n_heads) for _ in range(n_layers)])
        self.norm = nn.LayerNorm(d_model)
        self.lm_head = nn.Linear(d_model, vocab_size, bias=False)
        self.lm_head.weight = self.embed.weight
        
        # ACTIVAR GRADIENT CHECKPOINTING
        self.gradient_checkpointing = True
        self.apply(self._init_weights)

    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            module.weight.data.normal_(mean=0.0, std=0.02)
            if module.bias is not None: module.bias.data.zero_()
        elif isinstance(module, nn.Embedding):
            module.weight.data.normal_(mean=0.0, std=0.02)

    def forward(self, input_ids, labels=None):
        x = self.embed(input_ids)
        cos, sin = self.rope(x, x.shape[1])
        
        for block in self.blocks:
            if self.gradient_checkpointing and self.training:
                x = checkpoint(block, x, cos, sin, use_reentrant=False)
            else:
                x = block(x, cos, sin)
            
        x = self.norm(x)
        logits = self.lm_head(x)
        
        if labels is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), labels.view(-1), ignore_index=-100)
            return loss, logits
        return logits

# ==============================================================================
# 3. BUCLE DE ENTRENAMIENTO
# ==============================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--batch_size", type=int, default=2)
    parser.add_argument("--mode", type=str, default="bidirectional")
    parser.add_argument("--checkpoint_dir", type=str, default="checkpoints_v2")
    args = parser.parse_args()

    root = Path.cwd()
    tokenizer_path = root / "data" / "tokenizer_v2" / "dna_sp_v2.model"
    processed_dir = root / "data" / "processed_v2"
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(exist_ok=True)
    
    npz_files = list(processed_dir.glob("*.npz"))
    if not npz_files:
        print("❌ No hay archivos .npz v2.")
        return
    
    np.random.seed(42)
    np.random.shuffle(npz_files)
    split = int(len(npz_files) * 0.9)
    train_files, val_files = npz_files[:split], npz_files[split:]
    
    VOCAB_SIZE = spm.SentencePieceProcessor(model_file=str(tokenizer_path)).get_piece_size()
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Dispositivo: {device}")
    
    model = DNAModelV2(vocab_size=VOCAB_SIZE).to(device)
    total_params = sum(p.numel() for p in model.parameters())
    print(f"Parámetros del modelo v2: {total_params / 1e6:.1f}M")
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    
    train_dataset = CircularDNADataset(train_files, tokenizer_path)
    val_dataset = CircularDNADataset(val_files, tokenizer_path)
    
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=2)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, num_workers=2)
    
    for epoch in range(args.epochs):
        model.train()
        train_loss = 0
        pbar = tqdm(train_loader, desc=f"Época {epoch+1}/{args.epochs}")
        for batch in pbar:
            input_ids = batch["input_ids"].to(device)
            labels = batch["labels"].to(device)
            
            optimizer.zero_grad()
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                loss, _ = model(input_ids, labels=labels)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
            pbar.set_postfix({"loss": f"{loss.item():.4f}"})
            
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch in val_loader:
                input_ids = batch["input_ids"].to(device)
                labels = batch["labels"].to(device)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                    loss, _ = model(input_ids, labels=labels)
                val_loss += loss.item()
                
        avg_train = train_loss / len(train_loader)
        avg_val = val_loss / len(val_loader)
        print(f"Época {epoch+1} | Train Loss: {avg_train:.4f} | Val Loss: {avg_val:.4f}")
        
        torch.save(model.state_dict(), checkpoint_dir / f"dna_model_v2_epoch{epoch+1}.pt")

if __name__ == "__main__":
    main()
