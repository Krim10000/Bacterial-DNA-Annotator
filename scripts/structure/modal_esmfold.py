"""ESMFold en Modal (GPU) — DEPRECADO.

> **Desde 2026-11-21 ESMFold se ejecuta vía neurosnap.ai**, no en Modal. Este
> script queda como referencia histórica: la app `dna-esmfold` está parada y no
> debe relanzarse. Ver `docs/EXTERNAL_SOURCES.md`.

Requiere `modal token new` y autorización de gasto.

POLÍTICA GPU: Modal limita a 10 GPU simultáneas por cuenta. Este app usa
max_containers=1. Anunciar "USANDO MODAL GPU" antes y cerrar con
`modal app stop <app-id>` + verificar `modal container list` después.

Uso:
  modal run scripts/structure/modal_esmfold.py \
    --fasta data/metabolism/targets/structural_candidates.faa \
    --out data/metabolism/targets/structures/esmfold --chunk 5
"""
import json
from pathlib import Path

import modal

app = modal.App("dna-esmfold")
cache = modal.Volume.from_name("dna-esmfold-cache", create_if_missing=True)

image = (
    modal.Image.from_registry("nvidia/cuda:11.8.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git", "wget", "libgl1", "build-essential", "clang")
    .pip_install("numpy==1.26.4", "Cython", "ninja", "setuptools<81", "wheel")
    .pip_install("torch==2.1.0", index_url="https://download.pytorch.org/whl/cu118")
    .pip_install("fair-esm[esmfold]==2.0.0", "modelcif")
    .add_local_file("scripts/structure/patch_fair_esm.py", "/tmp/patch_fair_esm.py", copy=True)
    .run_commands("python /tmp/patch_fair_esm.py")
    .env({"TORCH_CUDA_ARCH_LIST": "7.5"})
    .pip_install(
        "dllogger @ git+https://github.com/NVIDIA/dllogger.git",
        "openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307",
        extra_options="--no-build-isolation",
    )
    .run_commands("python -m pip uninstall -y deepspeed")
    .env({"TORCH_HOME": "/cache/torch", "HF_HOME": "/cache/hf", "CUDA_HOME": "/usr/local/cuda"})
)


@app.function(image=image, gpu="T4", timeout=1800, max_containers=1, volumes={"/cache": cache})
def fold_batch(seqs: dict[str, str]) -> dict[str, str]:
    import re
    from pathlib import Path
    root = Path("/usr/local/lib/python3.11/site-packages/esm/esmfold/v1")
    for fp in root.glob("*.py"):
        ts = fp.read_text()
        if "from dataclasses import dataclass" not in ts and "dataclass" not in ts:
            continue
        ts = ts.replace("from dataclasses import dataclass", "from dataclasses import dataclass, field")
        ts = re.sub(r"(\w+):\s*(\w+Config)\s*=\s*\2(?:\(\))?",
                    r"\1: \2 = field(default_factory=\2)", ts)
        ts = ts.replace("trunk: T.Any = FoldingTrunkConfig()",
                        "trunk: T.Any = field(default_factory=FoldingTrunkConfig)")
        fp.write_text(ts)
    for pyc in root.glob("__pycache__/*.pyc"):
        pyc.unlink(missing_ok=True)
    import torch
    import esm

    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    model.set_chunk_size(128)
    out = {}
    for name, seq in seqs.items():
        seq = seq.replace("*", "").strip().upper()
        if len(seq) < 30:
            out[name] = ""
            continue
        with torch.no_grad():
            out[name] = model.infer_pdb(seq)
    return out


def read_fasta(path: str) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip().replace("|", "__")
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


@app.local_entrypoint()
def main(fasta: str, out: str, chunk: int = 5):
    seqs = read_fasta(fasta)
    outdir = Path(out)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(seqs)
    for i in range(0, len(names), chunk):
        batch = {n: seqs[n] for n in names[i:i + chunk]}
        result = fold_batch.remote(batch)
        for name, pdb in result.items():
            if pdb:
                (outdir / f"{name}.pdb").write_text(pdb)
        print(f"folded {min(i + chunk, len(names))}/{len(names)}", flush=True)
    (outdir / "manifest.json").write_text(json.dumps({"n": len(names)}, indent=2))
