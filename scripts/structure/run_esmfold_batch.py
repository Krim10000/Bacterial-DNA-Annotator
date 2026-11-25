#!/usr/bin/env python3
"""ESMFold por lotes, reanudable, para Kaggle/Colab/GPU local gratuita."""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True, help="TSV con uniprot,length,sequence")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--chunk-size", type=int, default=64)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args()
    out = Path(args.outdir); pdbdir = out / "pdb"; pdbdir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(args.batch), delimiter="\t"))
    if args.limit: rows = rows[:args.limit]
    import torch
    import esm
    device = "cpu" if args.cpu or not torch.cuda.is_available() else "cuda"
    print("device", device, "batch", len(rows), flush=True)
    model = esm.pretrained.esmfold_v1().eval().to(device)
    if device == "cuda": model.set_chunk_size(args.chunk_size)
    metrics = out / "metrics.tsv"
    exists = metrics.exists()
    with open(metrics, "a", newline="") as fh:
        fields = ["uniprot", "status", "length", "plddt", "ptm", "seconds", "error"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        if not exists: w.writeheader()
        for i, r in enumerate(rows, 1):
            acc = r["uniprot"]; target = pdbdir / f"{acc}.pdb"
            if target.exists() and target.stat().st_size > 1000:
                print(i, acc, "skip", flush=True); continue
            seq = r["sequence"].strip()
            t0 = time.perf_counter()
            try:
                with torch.no_grad():
                    o = model.infer(seq)
                    pdb = model.output_to_pdb(o)[0]
                    plddt = float(o["plddt"].mean())
                    ptm = float(o.get("ptm", 0.0))
                target.write_text(pdb)
                w.writerow({"uniprot": acc, "status": "done", "length": len(seq),
                            "plddt": round(plddt, 2), "ptm": round(ptm, 4),
                            "seconds": round(time.perf_counter() - t0, 2), "error": ""})
            except Exception as exc:
                w.writerow({"uniprot": acc, "status": "error", "length": len(seq),
                            "seconds": round(time.perf_counter() - t0, 2),
                            "error": f"{type(exc).__name__}: {exc}"[:300]})
                print(i, acc, "ERROR", exc, flush=True)
            fh.flush()
            if i % 10 == 0: print(i, len(rows), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
