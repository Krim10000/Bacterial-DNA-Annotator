#!/usr/bin/env python3
"""F2d — especificidad del detector de bisagra: positivos vs negativos.

El hallazgo de LacI (modo DBD-dominado y antiparalelo) podría ser genérico: cualquier
proteína con un extremo N-terminal flexible daría "enriquecimiento alto" en la ventana
1-62. Para saberlo hace falta un **nulo interno por proteína**: comparar la ventana del
DBD contra 200 ventanas aleatorias del mismo tamaño.

Métrica:
  z_dbd(k) = (varianza_DBD(k) − media_ventanas_aleatorias(k)) / sd_ventanas_aleatorias(k)

  bisagra(k) = z_dbd > 2 y coseno antiparalelo entre cadenas < −0,5

Casos:
  POSITIVO  1LBG  LacI tetrámero apo (interruptor alostérico de manual)
  NEGATIVOS 1LMB λ cI · 5CRO λ Cro · 3CRO 434 Cro · 2OR1 434 cI
            1JHH LexA · 1F36 Fis · 1IHF IHF · 1B0N SinR
  (todos HTH, ninguno con efector pequeño alostérico)

    modal run scripts/allostery/08_hinge_specificity.py
"""

# NOTA MODAL: Modal admite 10 GPU simultáneas por workspace. Antes de correr esto,
# avisa que usas Modal y al terminar cierra con:
#     bash scripts/allostery/modal_guard.sh stop-all
#     bash scripts/allostery/modal_guard.sh status   # debe salir 0 activas
# Ver docs/MODAL.md.
from __future__ import annotations

import json
from pathlib import Path

import modal

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if len(_HERE.parents) > 2 else Path("/root")
OUT = ROOT / "data" / "allostery"

app = modal.App("alloreg-hinge")
vol = modal.Volume.from_name("alloreg-anm-results", create_if_missing=True)

image = modal.Image.debian_slim(python_version="3.11").pip_install("prody==2.6.1", "numpy")

CASES = [
    {"name": "LacI", "pdb": "1lbg", "role": "positivo"},
    {"name": "lambda_cI", "pdb": "1lmb", "role": "negativo"},
    {"name": "lambda_Cro", "pdb": "5cro", "role": "negativo"},
    {"name": "434_Cro", "pdb": "3cro", "role": "negativo"},
    {"name": "434_cI", "pdb": "2or1", "role": "negativo"},
    {"name": "LexA", "pdb": "1jhh", "role": "negativo"},
    {"name": "Fis", "pdb": "1f36", "role": "negativo"},
    {"name": "IHF", "pdb": "1ihf", "role": "negativo"},
    {"name": "SinR", "pdb": "1b0n", "role": "negativo"},
]


@app.function(image=image, gpu="A10G", timeout=1800, volumes={"/results": vol})
def analyze(case: dict) -> dict:
    import numpy as np
    import prody

    prot = prody.parsePDB(case["pdb"])
    ca = prot.select("name CA and protein")
    if ca is None or ca.numAtoms() < 40:
        return {**case, "error": "pocos CA"}
    if len(np.unique(ca.getChids())) > 1 and ca.numAtoms() > 1500:
        ca = ca.select("chain A")  # se conserva una cadena si es enorme

    anm = prody.ANM(case["name"])
    anm.buildHessian(ca)
    anm.calcModes(n_modes=5)
    vecs = np.asarray(anm.getEigvecs())
    eig = np.asarray(anm.getEigvals())[:5]

    resnums = ca.getResnums()
    chids = np.array([str(c) for c in ca.getChids()])
    n = len(resnums)
    L = min(62, n - 10)
    rng = np.random.default_rng(0)

    def window_mask(start):
        return (resnums >= start) & (resnums < start + L)

    real = np.zeros(n, dtype=bool)
    for ch in np.unique(chids):
        m = chids == ch
        real |= window_mask(resnums[m].min()) & m

    rand_fracs = {k: [] for k in range(vecs.shape[1])}
    for _ in range(200):
        start = int(rng.integers(int(resnums.min()), max(int(resnums.min()) + 1,
                                                        int(resnums.max()) - L + 2)))
        m = window_mask(start)
        if m.sum() < 10:
            continue
        for k in range(vecs.shape[1]):
            v = vecs[:, k].reshape(-1, 3)
            mag2 = (v ** 2).sum(1)
            rand_fracs[k].append(float(mag2[m].sum() / mag2.sum()))

    modes = []
    for k in range(vecs.shape[1]):
        v = vecs[:, k].reshape(-1, 3)
        mag2 = (v ** 2).sum(1)
        frac = float(mag2[real].sum() / mag2.sum())
        arr = np.array(rand_fracs[k])
        z = float((frac - arr.mean()) / arr.std()) if arr.std() > 0 else None
        # antiparalelismo entre cadenas (solo si hay >1)
        anti = None
        per_chain = {}
        for ch in np.unique(chids):
            m = real & (chids == ch)
            if m.sum():
                per_chain[ch] = v[m].mean(0)
        if len(per_chain) >= 2:
            ks = sorted(per_chain)
            a, b = per_chain[ks[0]], per_chain[ks[-1]]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            if na and nb:
                anti = float((a @ b) / (na * nb))
        modes.append({"modo": k + 1, "autovalor": float(eig[k]),
                      "frac_dbd": round(frac, 3), "z_dbd": None if z is None else round(z, 2),
                      "anti": None if anti is None else round(anti, 3)})

    hinge = any((m["z_dbd"] or -9) > 2 and (m["anti"] is not None and m["anti"] < -0.5)
                for m in modes)
    rec = {**case, "n_res": int(n), "ventana": L, "hinge": hinge, "modos": modes}
    Path("/results", f"hinge_{case['name']}.json").write_text(json.dumps(rec))
    vol.commit()
    return rec


@app.local_entrypoint()
def main(mode: str = "run"):
    if mode == "fetch":
        data = read_all.remote()
        (OUT / "hinge_specificity.json").write_text(json.dumps(data, indent=1, ensure_ascii=False))
        for r in data:
            best = max((m["z_dbd"] or -9) for m in r.get("modos", []))
            print(f"{r['name']:12s} {r['role']:9s} hinge={r.get('hinge')} z_max={best} "
                  f"anti={[m['anti'] for m in r.get('modos', [])][:3]}")
        return
    res = list(analyze.map(CASES))
    print(json.dumps(res, indent=1, ensure_ascii=False)[:2000])


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": vol}, timeout=300)
def read_all() -> list[dict]:
    import glob
    return [json.load(open(f)) for f in sorted(glob.glob("/results/hinge_*.json"))]
