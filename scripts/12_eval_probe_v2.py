import torch, numpy as np, sentencepiece as spm, sys, random, re
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score

# Importar utilidades
sys.path.insert(0, str(Path(__file__).parent))
from lib.io import read_config, read_fasta

# ==============================================================================
# ARQUITECTURA V2 (Copiada para evitar problemas de importación)
# ==============================================================================
import torch.nn as nn
import torch.nn.functional as F

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
        
        mask = torch.ones(T, T, device=x.device, dtype=torch.bool)
        mask = mask.triu(diagonal=-self.window_size).tril(diagonal=self.window_size)
        
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
        self.gradient_checkpointing = False

# ==============================================================================
# EVALUACIÓN
# ==============================================================================
def get_windows(seq, motif, num_windows=200, window_size=2048):
    windows = []
    idxs = [m.start() for m in re.finditer(motif, seq)]
    random.shuffle(idxs)
    for i in idxs[:num_windows]:
        start = max(0, i - window_size // 2)
        w = seq[start:start + window_size]
        if len(w) == window_size: windows.append(w)
    return windows

def main():
    root = Path.cwd()
    
    # Rutas al v2
    tokenizer_path = root / "data" / "tokenizer_v2" / "dna_sp_v2.model"
    checkpoint_path = root / "dna_model_v2_epoch6.pt"
    
    if not tokenizer_path.exists():
        print(f"❌ Falta el tokenizador: {tokenizer_path}")
        return
    if not checkpoint_path.exists():
        print(f"❌ Falta el checkpoint: {checkpoint_path}")
        return

    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModelV2(vocab_size=sp.get_piece_size())
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
    model.eval()

    fasta = list((root / "data/raw").glob("*.fna.gz"))[0]
    seq = next(read_fasta(fasta))[1]

    print("Extrayendo ventanas (Promotores vs Aleatorias)...")
    pos_windows = get_windows(seq, "TATAAT", 200)
    neg_windows = get_windows(seq, "GCGCGC", 200) 

    def get_embedding(seq_str):
        ids = sp.encode(seq_str, out_type=int)[:2048]
        ids = ids + [0] * (2048 - len(ids))
        with torch.no_grad():
            x = model.embed(torch.tensor([ids]))
            cos, sin = model.rope(x, x.shape[1])
            for block in model.blocks:
                x = block(x, cos, sin)
            x = model.norm(x)
            return x[0].mean(dim=0).numpy()

    print("Generando embeddings (esto puede tardar unos minutos)...")
    X = np.array([get_embedding(w) for w in pos_windows + neg_windows])
    y = np.array([1]*len(pos_windows) + [0]*len(neg_windows))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print("Entrenando sonda lineal...")
    probe = LogisticRegression(max_iter=1000)
    probe.fit(X_train, y_train)

    preds = probe.predict(X_test)
    probs = probe.predict_proba(X_test)[:, 1]

    print("\n=== Resultados de la Sonda Lineal (Modelo v2) ===")
    print(f"Accuracy: {accuracy_score(y_test, preds):.4f}")
    print(f"AUROC:    {roc_auc_score(y_test, probs):.4f}")
    print("(Si AUROC > 0.80, el modelo capturó la gramática de los promotores).")

if __name__ == "__main__":
    main()
