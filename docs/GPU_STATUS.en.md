🌐 [Español](GPU_STATUS.md) · **English**

# GPU status and pending processes

**Summary:** annotation, the browser, metabolic gaps, operons and the interface
work **without a GPU**. What remains GPU-dependent is optional or was closed out
in the cloud: full training of the DNA model, its evaluation, and local
ESMFold.

Reference machine for this version: `torch 2.14.1+cpu`,
`torch.cuda.is_available() == False`, no `nvidia-smi` and no `/dev/nvidia*`.
The GUI detects this and shows it in **⏱️ What do I need and how long does it take?**.

## 1 · What needs a GPU and what does not

| Stage | Status in v1.0 | GPU? | How it runs | Impact |
|---|---|---|---|---|
| Annotation (CDS, tRNA, rRNA, ghost genes) | ✅ works | No | GUI or `dna-annotate` | core |
| Gaps, closures, transport, operons | ✅ works | No | GUI / `run_pipeline.sh` | core |
| Genome browser | ✅ works | No | GUI | core |
| Structural filters (Foldseek, P2Rank, AlphaFill) | ✅ works | No | GUI / `scripts/structure/` | included |
| ESMFold via Neurosnap (API) | ☁️ optional | No (remote) | `--esmfold "neuros nap"` | ~0.04 credits/protein |
| Local ESMFold | ⛔ blocked | ≥16 GB VRAM | `scripts/structure/04_local_esmfold.py` | not included |
| ESMFold on Modal | ⛔ never run | A10G | `scripts/structure/modal_esmfold.py` (deprecated) | not included |
| Train DNA model v2 | ⚠️ partial (epoch 6 and 9 of 25) | A10G | `modal_train_v2.py` | not used by the annotator |
| Full model evaluation | ⛔ blocked | GPU | `scripts/13_comprehensive_eval.py --n 100` | not included |
| Boltz-2 at scale | ☁️ partial (24 predictions) | A10G | `scripts/allostery/04_boltz_modal.py` | evidence collected |
| MM/GBSA, ABFE, molecular dynamics | ⛔ outside v1.0 | GPU | `docs/ALOSTERIA_V1_CIERRE.md` | closed |

## 2 · Evidence of what remained blocked

- **DNA model.** `modal_train_v2.py:16` requires `gpu="A10G"`. Only the
  checkpoints from epochs 6 and 9 exist (`dna_model_v2_epoch6.pt`,
  `dna_model_v2_epoch9.pt`, 48.9 MB, ~13.8 M parameters). The 25-epoch plan
  did not finish: `checkpoints_descargados/checkpoints_v2/` is empty.
  `logs/09_train_model.log` is 0 bytes: it was never attempted locally.
- **Evaluation.** `reports/README.md:29-35` documents that the local run was
  cancelled after 45 minutes; `reports/eval_comprehensive.json` uses only 8 windows
  and `device: "cpu"`, so it **is not an evaluation**. `ROADMAP.md` claimed
  97.3 M parameters and 25 epochs; the actual weights are 11–14 M and 9–10 epochs.
  This version corrects the claim and does not present the model as evaluated.
- **Local ESMFold.** `scripts/structure/04_local_esmfold.py:38-40` aborts without
  CUDA; `webapp/esmfold_tools.py:78-88` reports "no NVIDIA GPU visible ·
  PyTorch without CUDA · missing openfold". The alternative is the Neurosnap API.
- **Boltz-2.** 24 predictions were made and the result was negative
  (`docs/NEGATIVE_RESULTS.md:11-24`); the campaign of ~900 co-folds was not launched
  (20 h of A10G, 25–40 USD).
- **Logs.** There is not a single occurrence of OOM/CUDA across the 913 logs: the heavy
  work was outsourced to Modal or Neurosnap, it did not fail locally.

## 3 · What this means for people using the app

- Without a GPU you get **complete annotation and prioritization**. What you do not get
  are free-energy predictions or the evaluation of the DNA model.
- The app **does not try** to load CUDA or download weights: if there is no GPU, the
  options that need one are shown as unavailable, with the reason.
- ESMFold via Neurosnap sends the sequence to an external service and consumes
  credit. The app warns about this before using it; it does not activate on its own.

## 4 · How to continue on a machine with a GPU

```bash
# 1. Verify that the GPU exists and that torch sees it
nvidia-smi
python -c "import torch; print(torch.cuda.is_available())"

# 2. Training (Modal; review docs/MODAL_GPU_POLICY.md first)
modal run modal_train_v2.py

# 3. Evaluation with the local GPU
python scripts/13_comprehensive_eval.py --n 100 --window 2048

# 4. Local ESMFold (requires weights and OpenFold)
python scripts/structure/04_local_esmfold.py --help
```

Documented indicative cost: Boltz-2 at scale 25–40 USD; ABFE 100–150 GPU-h
per system; MD 10 ns × 3 replicas 51–135 USD. None of this is necessary to
use the annotator.

## 4 · Kaggle as a first GPU filter

Kaggle is used to test ESMFold, co-folding and large models **before** spending
local GPU time or cloud budget. It is not a dependency of the annotator:

- sequences are uploaded manually and only with consent;
- results are imported as tables or structures;
- the app does not connect to Kaggle at startup;
- annotation, gaps and operons keep running on CPU.

The repressor line is detailed in
[`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md) and
[`ROADMAP.md`](ROADMAP.md).

## 5 · Decision for v1.0

The DNA model is released as **experimental and not evaluated**. The
checkpoints are not included in the repository: how to regenerate them is documented.
Annotation and the GUI are the stable product of this version.
