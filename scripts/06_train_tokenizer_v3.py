#!/usr/bin/env python3
"""Entrena el tokenizador v3 con las correcciones de diseño.

Cambios frente a v2:
  - default vocab 256, max_sentencepiece_length=4 (v2: 4096/6)
    (el barrido 1024/5, 512/5, 256/4 está en docs/TOKENIZER_V3_REPORT.md)
  - add_dummy_prefix=False (elimina el token ▁)
  - special tokens explícitos: <pad>, <unk>, <bos>, <eos>, <mask>
  - IUPAC ambiguo protegido como símbolo: R Y K M S W B D H V N U
  - motivos protegidos: TATAAT, TTGACA, AGGAGG, GGAGG, ATG, GTG, TTG,
    TAA, TAG, TGA
  - normalization_rule_name=identity (no NFKC sobre ADN)

Salida:
  data/tokenizer_v3/dna_sp_v3.model
  data/tokenizer_v3/dna_sp_v3.vocab
  data/tokenizer_v3/tokenizer_meta_v3.json

Uso:
  python scripts/06_train_tokenizer_v3.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import sentencepiece as spm

ROOT = Path(__file__).resolve().parents[1]

MOTIFS = [
    "TATAAT", "TTGACA", "AGGAGG", "GGAGG",
    "ATG", "GTG", "TTG", "TAA", "TAG", "TGA",
]
IUPAC = ["R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N", "U"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vocab-size", type=int, default=256)
    ap.add_argument("--max-piece", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "data" / "tokenizer_v3"))
    ap.add_argument("--corpus", default=str(ROOT / "data" / "tokenizer" / "corpus.txt"))
    args = ap.parse_args()

    corpus = Path(args.corpus)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "dna_sp_v3"

    if not corpus.exists():
        raise SystemExit(f"No existe el corpus: {corpus}")

    user_symbols = MOTIFS + IUPAC
    (out_dir / "user_symbols_v3.txt").write_text("\n".join(user_symbols) + "\n")

    print(f"Corpus: {corpus} ({corpus.stat().st_size / 1e6:.0f} MB)")
    print(f"Entrenando v3: vocab={args.vocab_size}, max_piece={args.max_piece}, "
          f"sin ▁, IUPAC completo ...")

    spm.SentencePieceTrainer.train(
        input=str(corpus),
        model_prefix=str(prefix),
        vocab_size=args.vocab_size,
        model_type="bpe",
        user_defined_symbols=user_symbols,
        control_symbols=["<mask>"],
        pad_piece="<pad>",
        unk_piece="<unk>",
        bos_piece="<bos>",
        eos_piece="<eos>",
        pad_id=0,
        unk_id=1,
        bos_id=2,
        eos_id=3,
        add_dummy_prefix=False,
        remove_extra_whitespaces=False,
        normalization_rule_name="identity",
        max_sentencepiece_length=args.max_piece,
        max_sentence_length=2048,
        split_by_whitespace=False,
        split_by_unicode_script=False,
        split_by_number=False,
        split_digits=False,
        byte_fallback=False,
        character_coverage=1.0,
        input_sentence_size=200_000,
        shuffle_input_sentence=True,
        num_threads=4,
    )

    sp = spm.SentencePieceProcessor(model_file=str(prefix) + ".model")
    meta = {
        "version": 3,
        "vocab_size": sp.get_piece_size(),
        "max_piece_length": args.max_piece,
        "user_symbols": user_symbols,
        "control_symbols": ["<mask>"],
        "special_tokens": {
            "pad": {"id": sp.pad_id(), "piece": sp.id_to_piece(sp.pad_id())},
            "unk": {"id": sp.unk_id(), "piece": sp.id_to_piece(sp.unk_id())},
            "bos": {"id": sp.bos_id(), "piece": sp.id_to_piece(sp.bos_id())},
            "eos": {"id": sp.eos_id(), "piece": sp.id_to_piece(sp.eos_id())},
        },
        "add_dummy_prefix": False,
        "normalization": "identity",
        "trained_at": datetime.now(timezone.utc).isoformat(),
        "corpus": str(corpus.relative_to(ROOT)),
        "corpus_bytes": corpus.stat().st_size,
    }
    # El id de <mask> no lo expone SentencePiece como campo; se busca.
    meta["special_tokens"]["mask"] = {
        "id": sp.piece_to_id("<mask>"),
        "piece": "<mask>",
    }
    (out_dir / "tokenizer_meta_v3.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False))

    # Comprobaciones mínimas.
    ok_motifs = all(len(sp.encode(m, out_type=int)) == 1 for m in MOTIFS)
    ok_iupac = all(sp.piece_to_id(c) != sp.unk_id() for c in IUPAC)
    no_prefix = all(sp.id_to_piece(i) != "▁" for i in sp.encode("ACGT", out_type=int))
    print("\nComprobaciones:")
    print(f"  motivos como token único: {ok_motifs}")
    print(f"  IUPAC sin <unk>:          {ok_iupac}")
    print(f"  sin ▁:                    {no_prefix}")
    print(f"  vocab:                    {sp.get_piece_size()}")
    print(f"  <mask> id:                {sp.piece_to_id('<mask>')}")
    print(f"  guardado en:              {prefix}.model")


if __name__ == "__main__":
    main()
