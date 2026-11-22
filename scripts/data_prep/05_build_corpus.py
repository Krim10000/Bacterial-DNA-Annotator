#!/usr/bin/env python3
"""Concatena FASTA filtrados en corpus.txt para entrenar SentencePiece.
   Trocea a 1 Mbp por línea para evitar truncamiento en SP."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.io import read_config, read_fasta, setup_logging

CHUNK = 2_048

def main():
    cfg = read_config()
    log = setup_logging("05_build_corpus", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    raw = root / cfg["data"]["raw_dir"]
    tok = root / cfg["data"]["tokenizer_dir"]
    tok.mkdir(parents=True, exist_ok=True)
    out = tok / "corpus.txt"

    min_len = cfg["precompute"]["min_seq_length"]
    max_n = cfg["precompute"]["max_n_fraction"]
    fastas = sorted(raw.glob(cfg["data"]["fasta_pattern"]))
    if not fastas:
        log.error(f"No FASTA en {raw}")
        sys.exit(1)

    total_bp, n_seqs, n_chunks = 0, 0, 0
    with open(out, "w") as f:
        for fna in fastas:
            for name, seq in read_fasta(fna):
                if len(seq) < min_len:
                    continue
                if seq.count("N") / max(len(seq), 1) > max_n:
                    continue
                n_seqs += 1
                for i in range(0, len(seq), CHUNK):
                    chunk = seq[i:i+CHUNK]
                    if len(chunk) < min_len:
                        continue
                    f.write(chunk + "\n")
                    total_bp += len(chunk)
                    n_chunks += 1

    log.info(f"Corpus: {n_seqs} secuencias, {total_bp:,} bp, {n_chunks} chunks")
    log.info(f"Archivo: {out} ({out.stat().st_size/1e9:.2f} GB)")
    print(f"Corpus: {n_seqs:,} secuencias, {total_bp:,} bp → {out}")

if __name__ == "__main__":
    main()
