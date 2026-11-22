#!/usr/bin/env python3
"""F2b — modos normales (ANM) sobre las estructuras AlphaFold, en Modal.

Idea: si el estado activo/inactivo de un represor es un movimiento de dominios
rígidos (no plegamiento), el modo más blando del ANM debe mover el dominio de
unión a ADN (DBD) respecto al cuerpo. Este script mide, para cada TF:

  - autovalores de los 5 modos más blandos
  - desplazamiento medio del DBD vs el resto (ratio)
  - fracción de la varianza del modo que cae en el DBD

Un represor con interruptor alostérico debería mostrar ratio > 1,5 y fracción
de DBD muy por encima de su fracción de residuos (~15%).

Se corre en Modal con GPU para no ocupar el CPU local (el ANM es CPU, va sobrado).

    modal run scripts/allostery/06_anm_modal.py            # lanza y recoge
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
AF = ROOT / "data" / "allostery" / "work" / "af"

app = modal.App("alloreg-anm")
results_vol = modal.Volume.from_name("alloreg-anm-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("prody==2.6.1", "numpy")
    .add_local_dir(str(AF), remote_path="/st")
)


@app.function(image=image, gpu="A10G", timeout=900,
              volumes={"/results": results_vol})
def anm_tf(job: dict) -> dict:
    import numpy as np
    import prody

    tf = job["tf"]
    pdb = Path("/st") / f"{tf}.pdb"
    if not pdb.exists():
        return {"tf": tf, "error": "sin estructura"}
    prot = prody.parsePDB(str(pdb))
    ca = prot.select("name CA and chain A")
    if ca is None or ca.numAtoms() < 30:
        return {"tf": tf, "error": f"pocos CA ({0 if ca is None else ca.numAtoms()})"}

    anm = prody.ANM(tf)
    anm.buildHessian(ca)
    anm.calcModes(n_modes=5)
    eig = np.asarray(anm.getEigvals())[:5]

    resnums = ca.getResnums()
    dbd = set(range(job["dbd_from"], job["dbd_to"] + 1)) if job.get("dbd_from") else set()
    in_dbd = np.array([r in dbd for r in resnums])

    modes = []
    eigvecs = np.asarray(anm.getEigvecs())
    for k in range(min(3, eigvecs.shape[1])):
        v = eigvecs[:, k].reshape(-1, 3)
        mag = np.linalg.norm(v, axis=1)
        frac_dbd = float((mag[in_dbd] ** 2).sum() / (mag ** 2).sum()) if in_dbd.any() else None
        ratio = (float(mag[in_dbd].mean() / mag[~in_dbd].mean())
                 if in_dbd.any() and (~in_dbd).any() else None)
        top = resnums[int(np.argmax(mag))]
        modes.append({"modo": k + 1, "autovalor": float(eig[k]),
                      "frac_dbd": frac_dbd, "ratio_dbd_core": ratio,
                      "residuo_max": int(top)})

    rec = {"tf": tf, "n_res": int(ca.numAtoms()),
           "dbd_from": job.get("dbd_from"), "dbd_to": job.get("dbd_to"),
           "frac_residuos_dbd": float(in_dbd.mean()) if len(in_dbd) else None,
           "eigvals": [float(x) for x in eig], "modos": modes}
    Path("/results", f"{tf}.json").write_text(json.dumps(rec))
    results_vol.commit()
    return rec


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": results_vol}, timeout=300)
def read_all() -> list[dict]:
    import glob
    return [json.load(open(f)) for f in sorted(glob.glob("/results/*.json"))]


@app.local_entrypoint()
def main(mode: str = "run"):
    if mode == "fetch":
        data = read_all.remote()
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "anm_results.json").write_text(json.dumps(data, indent=1, ensure_ascii=False))
        print(f"recogidos {len(data)}")
        for r in sorted(data, key=lambda x: -(x.get("modos", [{}])[0].get("ratio_dbd_core") or 0)):
            m = (r.get("modos") or [{}])[0]
            print(f"{r['tf']:8s} res={r.get('n_res')} "
                  f"ratio_m1={m.get('ratio_dbd_core')} frac_m1={m.get('frac_dbd')}")
        return

    ranges = json.load(open(OUT / "work" / "dbd_ranges.json"))
    tfs = sorted(p.stem for p in AF.glob("*.pdb"))
    jobs = [{"tf": tf, **ranges.get(tf, {})} for tf in tfs]
    print(f"lanzando {len(jobs)} TFs", flush=True)
    for r in anm_tf.map(jobs):
        if "error" in r:
            print("ERROR", r)


if __name__ == "__main__":
    pass
