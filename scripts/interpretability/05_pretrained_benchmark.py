#!/usr/bin/env python3
"""Extrae embeddings de modelos pretrained para un windows.tsv dado.

Soporta:
  - prokbert  (neuralbioinfo/prokbert-mini-long, 384d, 4 kb)
  - ntv2      (InstaDeepAI/nucleotide-transformer-v2-50m-multi-species, 512d, 12 kb)
  - v2        (checkpoint local dna_model_v2_epoch9.pt, 384d)

Salida:
  <out>/windows.tsv
  <out>/embeddings.npy
  <out>/manifest.json

Uso:
  python scripts/interpretability/05_pretrained_benchmark.py \
      --model prokbert \
      --windows data/embeddings/_pretrained_subset/windows.tsv \
      --out data/embeddings/pretrained/prokbert_subset
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("HF_HOME", str(ROOT / ".cache" / "huggingface"))
os.environ.setdefault("TRANSFORMERS_CACHE", str(ROOT / ".cache" / "huggingface"))
sys.path.insert(0, str(ROOT / "scripts"))
from lib.io import read_fasta  # noqa: E402
from lib.dna_v2 import load_v2, embed_sequences  # noqa: E402

PROKBERT = "neuralbioinfo/prokbert-mini-long"
NTV2 = "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"


def mean_pool(h: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    m = mask.unsqueeze(-1).to(h.dtype)
    return (h * m).sum(dim=1) / m.sum(dim=1).clamp_min(1)


class ProkBERTAdapter:
    name = "prokbert"
    max_length = 4096

    def __init__(self):
        from transformers import AutoModel, AutoTokenizer
        self.tok = AutoTokenizer.from_pretrained(PROKBERT, trust_remote_code=True)
        self.model = AutoModel.from_pretrained(PROKBERT, trust_remote_code=True).eval()
        self.dim = self.model.config.hidden_size

    @torch.inference_mode()
    def embed(self, seqs: list[str], batch_size: int = 16) -> np.ndarray:
        out = []
        for i in range(0, len(seqs), batch_size):
            batch = self.tok(seqs[i:i + batch_size], return_tensors="pt",
                             padding=True, truncation=True,
                             max_length=self.max_length)
            h = self.model(**batch).last_hidden_state
            out.append(mean_pool(h, batch["attention_mask"]).float().cpu().numpy())
        return np.concatenate(out).astype(np.float32)


class NTV2Adapter:
    name = "ntv2"
    max_length = 2048

    def __init__(self):
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        self.tok = AutoTokenizer.from_pretrained(NTV2, trust_remote_code=True)
        self.model = AutoModelForMaskedLM.from_pretrained(
            NTV2, trust_remote_code=True).eval()
        self.dim = self.model.config.hidden_size

    @torch.inference_mode()
    def embed(self, seqs: list[str], batch_size: int = 32) -> np.ndarray:
        out = []
        for i in range(0, len(seqs), batch_size):
            batch = self.tok(seqs[i:i + batch_size], return_tensors="pt",
                             padding=True, truncation=True,
                             max_length=self.max_length)
            h = self.model(**batch, output_hidden_states=True).hidden_states[-1]
            out.append(mean_pool(h, batch["attention_mask"]).float().cpu().numpy())
        return np.concatenate(out).astype(np.float32)


class V2Adapter:
    name = "v2"
    max_length = 8192

    def __init__(self):
        self.model, self.sp = load_v2(
            ROOT / "dna_model_v2_epoch9.pt",
            ROOT / "data/tokenizer_v2/dna_sp_v2.model",
            device="cpu",
        )
        self.dim = self.model.embed.embedding_dim

    def embed(self, seqs: list[str], batch_size: int = 32) -> np.ndarray:
        emb, _ = embed_sequences(self.model, self.sp, seqs,
                                 batch_size=batch_size, device="cpu",
                                 max_tokens=self.max_length, pool="mean")
        return emb.astype(np.float32)


def get_adapter(name: str):
    if name == "prokbert":
        return ProkBERTAdapter()
    if name == "ntv2":
        return NTV2Adapter()
    if name == "v2":
        return V2Adapter()
    raise ValueError(name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["prokbert", "ntv2", "v2"], required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    meta = pd.read_csv(args.windows, sep="\t")
    if args.limit:
        meta = meta.head(args.limit).copy()
    print(f"Ventanas: {len(meta)} | modelo: {args.model}")

    # Si la tabla trae `sequence` (p. ej. negativos barajados), usarla.
    if "sequence" in meta.columns and meta["sequence"].notna().all():
        seqs = meta["sequence"].astype(str).tolist()
        genomes = sorted(meta["genome"].unique())
    else:
        genomes = sorted(meta["genome"].unique())
        seqs_by_genome = {}
        for g in genomes:
            for base in ("data/raw", "data/raw_excluded"):
                p = ROOT / base / f"{g}_genomic.fna.gz"
                if p.exists():
                    seqs_by_genome[g] = {n: s for n, s in read_fasta(p)}
                    break
            if g not in seqs_by_genome:
                raise SystemExit(f"Falta FASTA para {g}")
        seqs = [seqs_by_genome[r.genome][r.contig][r.start - 1:r.end]
                for r in meta.itertuples()]
    lengths = np.array([len(s) for s in seqs])

    t0 = time.time()
    adapter = get_adapter(args.model)
    load_s = time.time() - t0
    t1 = time.time()
    emb = adapter.embed(seqs, batch_size=args.batch_size)
    infer_s = time.time() - t1
    if emb.shape[0] != len(meta):
        raise RuntimeError(f"Desalineado: {emb.shape[0]} vs {len(meta)}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out / "windows.tsv", sep="\t", index=False)
    np.save(out / "embeddings.npy", emb)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "model": args.model,
        "model_repo": {"prokbert": PROKBERT, "ntv2": NTV2,
                       "v2": "dna_model_v2_epoch9.pt"}[args.model],
        "dim": int(emb.shape[1]),
        "n_windows": int(emb.shape[0]),
        "window_bp": int(meta["window_bp"].iloc[0]) if "window_bp" in meta else None,
        "genomes": genomes,
        "load_seconds": round(load_s, 2),
        "infer_seconds": round(infer_s, 2),
        "seconds_per_window": round(infer_s / max(1, len(meta)), 4),
        "batch_size": args.batch_size,
        "device": "cpu",
        "mean_nt": float(lengths.mean()),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"OK -> {out} | dim={emb.shape[1]} | "
          f"{manifest['seconds_per_window']:.3f}s/ventana | {infer_s:.1f}s total")


if __name__ == "__main__":
    main()
