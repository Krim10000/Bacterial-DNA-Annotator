#!/usr/bin/env python3
import torch, numpy as np, sentencepiece as spm, sys, random, re
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).parent))
from lib.io import read_config, read_fasta

# ==============================================================================
# ARQUITECTURA V2 (Embebida)
# ==============================================================================
class RoPE(nn.Module):
    def __init__(self, dim, max_seq_len=32000):
        super().__init__()
        inv_freq = 1.0 / (10000 ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer("inv_freq", inv_freq)
    def forward(self, x, seq_len):
        t = torch.arange(seq_len, device=x.device, dtype=self.inv_freq.dtype)
        freqs = torch.einsum("i,j->ij", t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        return emb.cos()[None, None, :, :], emb.sin()[None, None, :, :]

def apply_rotary_pos_emb(q, k, cos, sin):
    def rotate_half(x):
        x1, x2 = x[..., :x.shape[-1]//2], x[..., x.shape[-1]//2:]
        return torch.cat((-x2, x1), dim=-1)
    return (q * cos) + (rotate_half(q) * sin), (k * cos) + (rotate_half(k) * sin)

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
        self.mlp = nn.Sequential(nn.Linear(d_model, 4*d_model), nn.GELU(), nn.Linear(4*d_model, d_model))
    def forward(self, x, cos, sin):
        B, T, C = x.shape
        x_norm = self.norm1(x)
        q, k, v = self.qkv(x_norm).chunk(3, dim=-1)
        q = q.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        k = k.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        v = v.view(B, T, self.n_heads, self.d_head).transpose(1, 2)
        q, k = apply_rotary_pos_emb(q, k, cos, sin)
        mask = torch.ones(T, T, device=x.device, dtype=torch.bool).triu(diagonal=-self.window_size).tril(diagonal=self.window_size)
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
    def forward(self, input_ids):
        x = self.embed(input_ids)
        cos, sin = self.rope(x, x.shape[1])
        for block in self.blocks: x = block(x, cos, sin)
        x = self.norm(x)
        return self.lm_head(x)

# ==============================================================================
# UTILIDADES DE EVALUACIÓN
# ==============================================================================
def calculate_perplexity(model, sp, seq_str):
    ids = sp.encode(seq_str, out_type=int)[:2048]
    if len(ids) < 10: return 0.0
    input_ids = torch.tensor([ids])
    with torch.no_grad():
        logits = model(input_ids)
        shift_logits = logits[..., :-1, :].contiguous()
        shift_labels = input_ids[..., 1:].contiguous()
        loss_fct = torch.nn.CrossEntropyLoss()
        loss = loss_fct(shift_logits.view(-1, shift_logits.size(-1)), shift_labels.view(-1))
        return torch.exp(loss).item()

def get_windows(seq, motif, num_windows=100, window_size=2048):
    windows = []
    idxs = [m.start() for m in re.finditer(motif, seq)]
    random.shuffle(idxs)
    for i in idxs[:num_windows]:
        start = max(0, i - window_size // 2)
        w = seq[start:start + window_size]
        if len(w) == window_size: windows.append(w)
    return windows

def run_probe(model, sp, pos_windows, neg_windows):
    def get_embedding(seq_str):
        ids = sp.encode(seq_str, out_type=int)[:2048]
        real_len = len(ids)
        if real_len == 0: return np.zeros(384)
        ids_padded = ids + [0] * (2048 - real_len)
        with torch.no_grad():
            x = model.embed(torch.tensor([ids_padded]))
            cos, sin = model.rope(x, x.shape[1])
            for block in model.blocks: x = block(x, cos, sin)
            x = model.norm(x)
            return x[0, :real_len].mean(dim=0).numpy()

    X = np.array([get_embedding(w) for w in pos_windows + neg_windows])
    y = np.array([1]*len(pos_windows) + [0]*len(neg_windows))
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    probe = LogisticRegression(max_iter=1000)
    probe.fit(X_train, y_train)
    probs = probe.predict_proba(X_test)[:, 1]
    return roc_auc_score(y_test, probs)

# ==============================================================================
# MAIN
# ==============================================================================
def main():
    root = Path.cwd()
    tokenizer_path = root / "data/tokenizer_v2/dna_sp_v2.model"
    checkpoint_path = root / "dna_model_v2_epoch9.pt"
    
    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModelV2(vocab_size=sp.get_piece_size())
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
    model.eval()

    # Elegimos una bacteria al azar (distinta a la de la prueba anterior)
    fastas = list((root / "data/raw").glob("*.fna.gz"))
    test_fasta = random.choice(fastas[50:]) # Saltamos las primeras 50
    seq = next(read_fasta(test_fasta))[1]
    
    print(f"\n🧬 PRUEBAS DE ROBUSTEZ (Modelo v2 - Época 9)")
    print(f"Genoma de prueba: {test_fasta.name}")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # PRUEBA 1: GENERALIZACIÓN (Otra bacteria, otro filo)
    # -------------------------------------------------------------------------
    print("1. PRUEBA DE GENERALIZACIÓN (Cross-Species)")
    pos = get_windows(seq, "TATAAT", 100, 2048)
    neg = get_windows(seq, "GCGCGC", 100, 2048)
    auroc = run_probe(model, sp, pos, neg)
    print(f"   -> AUROC en bacteria desconocida: {auroc:.4f}")
    if auroc > 0.85:
        print("   ✅ ÉXITO: El modelo generaliza entre especies. No memorizó.")
    else:
        print("   ❌ ALERTA: El modelo podría estar sobreajustado a las bacterias de entrenamiento.")

    # -------------------------------------------------------------------------
    # PRUEBA 2: MUTAGÉNESIS IN SILICO (El modelo debe notar la mutación)
    # -------------------------------------------------------------------------
    print("\n2. PRUEBA DE MUTAGÉNESIS IN SILICO")
    # Secuencia real con promotor y RBS
    seq_real = "TTGACATATAATAGGAGGATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    seq_mutada = seq_real.replace("TATAAT", "GCGCGC")
    
    ppl_real = calculate_perplexity(model, sp, seq_real)
    ppl_mutada = calculate_perplexity(model, sp, seq_mutada)
    
    print(f"   -> Perplejidad Secuencia Real:   {ppl_real:.2f}")
    print(f"   -> Perplejidad Secuencia Mutada: {ppl_mutada:.2f}")
    
    if ppl_mutada > ppl_real * 1.1:
        print("   ✅ ÉXITO: El modelo se 'sorprende' al ver la mutación. Entiende el motivo.")
    else:
        print("   ❌ ALERTA: El modelo no distingue la mutación de la secuencia real.")

    # -------------------------------------------------------------------------
    # PRUEBA 3: GAP-FILLING (Genes vs Regiones Intergénicas)
    # -------------------------------------------------------------------------
    print("\n3. PRUEBA DE GAP-FILLING (Detección de Genes Ocultos)")
    # Tomamos un fragmento de un gen real y un fragmento aleatorio (basura)
    gen_real = seq[1000:3000]      # Asumimos que hay un gen aquí
    region_basura = seq[500000:502000] # Región probablemente intergénica
    
    ppl_gen = calculate_perplexity(model, sp, gen_real)
    ppl_basura = calculate_perplexity(model, sp, region_basura)
    
    print(f"   -> Perplejidad de un Gen Real:          {ppl_gen:.2f}")
    print(f"   -> Perplejidad de Región Intergénica:   {ppl_basura:.2f}")
    
    if ppl_basura > ppl_gen * 1.2:
        print("   ✅ ÉXITO: La IA detecta la diferencia. Es útil para encontrar genes ocultos.")
    else:
        print("   ⚠️ ADVERTENCIA: La perplejidad no discrimina bien. Necesitaremos más datos.")

    print("=" * 60)

if __name__ == "__main__":
    main()
