import torch, numpy as np, sentencepiece as spm, sys, random
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
sys.path.insert(0, str(Path(__file__).parent))
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
    checkpoint_path = root / "checkpoints_descargados/checkpoints_v2/dna_model_bidirectional_epoch6.pt"
    
    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModel(vocab_size=sp.get_piece_size())
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
    model.eval()

    fasta = list((root / "data/raw").glob("*.fna.gz"))[0]
    seq = next(read_fasta(fasta))[1]

    print("Extrayendo ventanas...")
    pos_windows = get_windows(seq, "TATAAT", 200)
    neg_windows = get_windows(seq, "GCGCGC", 200)

    def get_embedding(seq_str):
        ids = sp.encode(seq_str, out_type=int)[:2048]
        ids = ids + [0] * (2048 - len(ids))
        with torch.no_grad():
            out = model(torch.tensor([ids]))
            return out[0].mean(dim=0).numpy()

    print("Generando embeddings...")
    X = np.array([get_embedding(w) for w in pos_windows + neg_windows])
    y = np.array([1]*len(pos_windows) + [0]*len(neg_windows))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print("Entrenando sonda lineal...")
    probe = LogisticRegression(max_iter=1000)
    probe.fit(X_train, y_train)

    preds = probe.predict(X_test)
    probs = probe.predict_proba(X_test)[:, 1]

    print("\n=== Resultados de la Sonda Lineal ===")
    print(f"Accuracy: {accuracy_score(y_test, preds):.4f}")
    print(f"AUROC:    {roc_auc_score(y_test, probs):.4f}")
    print("(Si AUROC > 0.85, el modelo capturó la gramática de los promotores).")

if __name__ == "__main__":
    main()
