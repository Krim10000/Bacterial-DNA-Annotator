#!/usr/bin/env python3
"""Entrena SentencePiece v2 con límite de memoria (Corregido OOM)."""
import sys, json, datetime
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from lib.io import read_config, setup_logging
import sentencepiece as spm

def main():
    cfg = read_config()
    log = setup_logging("06_train_tokenizer_v2", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    tok_dir = root / "data" / "tokenizer_v2"
    tok_dir.mkdir(parents=True, exist_ok=True)
    corpus = root / "data" / "tokenizer" / "corpus.txt"
    
    if not corpus.exists():
        print("❌ No existe corpus.txt.")
        return

    model_prefix = tok_dir / "dna_sp_v2"
    
    user_symbols = [
        "TATAAT", "TTGACA", "AGGAGG",
        "ATG", "TAA", "TAG", "TGA"
    ]

    with open(tok_dir / "user_symbols_v2.txt", "w") as f:
        f.write("\n".join(user_symbols) + "\n")

    print(f"Entrenando tokenizador v2 (vocab=4096, max_len=6)...")
    spm.SentencePieceTrainer.train(
        input=str(corpus),
        model_prefix=str(model_prefix),
        vocab_size=4096,
        model_type="bpe",
        user_defined_symbols=user_symbols,
        max_sentencepiece_length=6,
        max_sentence_length=2048,
        split_by_whitespace=False,
        split_by_unicode_script=False,
        split_by_number=False,
        split_digits=False,
        byte_fallback=False,
        character_coverage=1.0,
        # --- PARÁMETROS PARA EVITAR OOM ---
        input_sentence_size=200_000, 
        shuffle_input_sentence=True,
        # ----------------------------------
        num_threads=4,
    )

    meta = {
        "vocab_size": 4096,
        "max_token_length": 6,
        "user_symbols": user_symbols,
        "trained_at": datetime.datetime.utcnow().isoformat(),
    }
    with open(tok_dir / "tokenizer_meta_v2.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"✅ Tokenizador v2 guardado en {model_prefix}.model")

if __name__ == "__main__":
    main()
