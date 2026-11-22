#!/usr/bin/env python3
"""ESMFold local (fair-esm) para una o pocas secuencias.

Uso:
  python scripts/structure/04_local_esmfold.py --fasta candidatos.faa --out esmfold_local

Requiere GPU NVIDIA, PyTorch CUDA, fair-esm[esmfold] y OpenFold compatible.
No descarga ni envía secuencias a ningún servicio.
"""
import argparse
from pathlib import Path


def fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    import torch
    import esm

    if not torch.cuda.is_available():
        raise SystemExit("ESMFold local requiere CUDA; usa Neurosnap API o una GPU NVIDIA.")
    model = esm.pretrained.esmfold_v1().eval().cuda()
    model.set_chunk_size(128)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    seqs = fasta(a.fasta)
    if a.limit:
        seqs = seqs[:a.limit]
    for name, seq in seqs:
        seq = seq.replace("*", "").upper()
        pdb = model.infer_pdb(seq)
        safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)[:45] or "protein"
        (out / f"{safe}.pdb").write_text(pdb)
        print("saved", out / f"{safe}.pdb", len(seq))


if __name__ == "__main__":
    main()
