#!/usr/bin/env python3
"""Entrena SentencePiece BPE con motivos biologicos protegidos."""
import sys, json, datetime
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.io import read_config, setup_logging
import sentencepiece as spm

def main():
    cfg = read_config()
    log = setup_logging("06_train_tokenizer", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    tok_dir = root / cfg["data"]["tokenizer_dir"]
    corpus = tok_dir / "corpus.txt"
    model_prefix = tok_dir / "dna_sp"
    model_file = Path(str(model_prefix) + ".model")

    if not corpus.exists():
        log.error(f"No existe corpus: {corpus}. Corre 05 primero.")
        sys.exit(1)
    if model_file.exists():
        log.info(f"Modelo ya existe: {model_file}. Borralo para reentrenar.")
        print(f"Ya existe {model_file}. Saltando.")
        return

    tcfg = cfg["tokenizer"]
    user_symbols = tcfg["user_symbols"]
    with open(tok_dir / "user_symbols.txt", "w") as f:
        f.write("\n".join(user_symbols) + "\n")

    log.info(f"SP: vocab={tcfg['vocab_size']}, type={tcfg['model_type']}, "
             f"{len(user_symbols)} motivos protegidos")
    print(f"Entrenando SentencePiece (puede tardar 15-60 min)...")
    print(f"  vocab={tcfg['vocab_size']}, {len(user_symbols)} motivos protegidos")

    spm.SentencePieceTrainer.train(
        input=str(corpus),
        model_prefix=str(model_prefix),
        vocab_size=tcfg["vocab_size"],
        model_type=tcfg["model_type"],
        user_defined_symbols=user_symbols,
        max_sentencepiece_length=tcfg["max_token_length"],
        max_sentence_length=2_048,
        split_by_whitespace=False,
        split_by_unicode_script=False,
        split_by_number=False,
        split_digits=False,
        byte_fallback=tcfg["byte_fallback"],
        character_coverage=1.0,
        num_threads=cfg["precompute"].get("sp_num_threads", cfg["precompute"]["n_workers"]),
    )

    meta = {
        "vocab_size": tcfg["vocab_size"],
        "model_type": tcfg["model_type"],
        "user_symbols": user_symbols,
        "trained_at": datetime.datetime.utcnow().isoformat(),
        "corpus_size_bytes": corpus.stat().st_size,
    }
    with open(tok_dir / "tokenizer_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    sp = spm.SentencePieceProcessor(model_file=str(model_file))
    test = "TTGACATATAATAGGAGGATGAAATAAGCGC"
    pieces = sp.encode(test, out_type=str)
    print(f"\nTest tokenizacion:")
    print(f"  in:  {test}")
    print(f"  out: {pieces}")
    atomic = [p for p in pieces if p in user_symbols]
    print(f"  motivos atomicos: {len(atomic)}/{len([m for m in user_symbols if m in test])}")

    log.info(f"Tokenizador guardado: {model_file}")

if __name__ == "__main__":
    main()
