"""Evalúa si los embeddings del modelo agrupan secuencias por función."""
import torch, numpy as np, sentencepiece as spm, sys, random
from pathlib import Path
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.io import read_config, read_fasta
from lib.dna_model import DNAModel

def get_windows(seq, motif, num_windows=100, window_size=512):
    windows = []
    idxs = [m.start() for m in __import__('re').finditer(motif, seq)]
    random.shuffle(idxs)
    for i in idxs[:num_windows]:
        start = max(0, i - window_size // 2)
        w = seq[start:start + window_size]
        if len(w) == window_size: windows.append(w)
    return windows

def main():
    root = Path.cwd()
    tokenizer_path = root / "data/tokenizer/dna_sp.model"
    checkpoint_path = root / "checkpoints_descargados/checkpoints_v2/dna_model_bidirectional_epoch25.pt"
    
    if not checkpoint_path.exists():
        print(f"❌ Falta el checkpoint: {checkpoint_path}")
        return

    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModel(vocab_size=sp.get_piece_size())
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
    model.eval()

    # Tomar un genoma al azar
    fasta = list((root / "data/raw").glob("*.fna.gz"))[0]
    seq = next(read_fasta(fasta))[1]

    print("Extrayendo ventanas de ADN...")
    pos_windows = get_windows(seq, "TATAAT", 100) # Promotores
    neg_windows = get_windows(seq, "GCGCGC", 100) # Regiones ricas en GC (distintas a promotores)

    def get_embedding(seq_str):
        ids = sp.encode(seq_str, out_type=int)[:2048]
        ids = ids + [0] * (2048 - len(ids))
        with torch.no_grad():
            out = model(torch.tensor([ids]))
            return out[0].mean(dim=0).numpy() # Mean pooling

    print("Generando embeddings (puede tardar unos minutos)...")
    emb_pos = [get_embedding(w) for w in pos_windows]
    emb_neg = [get_embedding(w) for w in neg_windows]

    X = np.array(emb_pos + emb_neg)
    y = np.array([1]*len(emb_pos) + [0]*len(emb_neg))

    print("Reduciendo dimensionalidad y graficando...")
    pca = PCA(n_components=2)
    X_2d = pca.fit_transform(X)

    plt.figure(figsize=(10, 8))
    plt.scatter(X_2d[y==1, 0], X_2d[y==1, 1], c='red', label='Promotores (TATAAT)', alpha=0.6)
    plt.scatter(X_2d[y==0, 0], X_2d[y==0, 1], c='blue', label='Regiones GC', alpha=0.6)
    plt.legend()
    plt.title("Espacio Latente del Modelo de ADN")
    plt.savefig(root / "reports/clustering_eval.png")
    print("✅ Gráfico guardado en reports/clustering_eval.png")
    print("Si ves dos grupos separados, el modelo aprendió biología.")

if __name__ == "__main__":
    main()
