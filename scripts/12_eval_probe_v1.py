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
    checkpoint_path = root / "checkpoints_descargados/checkpoints/dna_model_bidirectional_epoch6.pt"
    
    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModel(vocab_size=sp.get_piece_size(), d_model=256, n_layers=4, n_heads=4, max_len=2048)
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
    model.eval()

    fasta = list((root / "data/raw").glob("*.fna.gz"))[0]
    seq = next(read_fasta(fasta))[1]

    print("Extrayendo ventanas...")
    pos_windows = get_windows(seq, "TATAAT", 50)
    neg_windows = get_windows(seq, "GCGCGC", 50)

    def get_embedding(seq_str):
        ids = sp.encode(seq_str, out_type=int)[:512]
        ids = ids + [0] * (512 - len(ids))
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
    import json
    payload = {"checkpoint": str(checkpoint_path), "accuracy": float(accuracy_score(y_test, preds)),
               "auroc": float(roc_auc_score(y_test, probs)), "n_windows": len(y),
               "pos_motif": "TATAAT", "neg_motif": "GCGCGC"}
    (root / "reports").mkdir(exist_ok=True)
    (root / "reports/eval_probe_v1.json").write_text(json.dumps(payload, indent=2))
    print("escrito: reports/eval_probe_v1.json")

if __name__ == "__main__":
    main()
