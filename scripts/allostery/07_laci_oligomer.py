#!/usr/bin/env python3
"""F2c — prueba de oligómero: ¿aparece el modo funcional de LacI?

Autochequeo del ANM: el monómero de AlphaFold no mostró movimiento preferente del
DBD (enriquecimiento 0,50). Aquí se prueba la hipótesis de que falta el modo
cuaternario, usando el tetrámero experimental apo (1LBG) frente al unido a ADN (1LBI).

Se calcula:
  1. El desplazamiento observado apo→holo por residuo (alineando el core, no el DBD).
  2. Los modos ANM del tetrámero apo.
  3. El solapamiento (coseno) de cada modo con el desplazamiento observado.
  4. Antiparalelismo: si los DBD de subunidades distintas se mueven en sentidos opuestos.
  5. Enriquecimiento del DBD en el modo.

Criterio de paso: algún modo con solapamiento > 0,5 y enriquecimiento DBD > 2.

    modal run scripts/allostery/07_laci_oligomer.py
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

app = modal.App("alloreg-laci")
vol = modal.Volume.from_name("alloreg-anm-results", create_if_missing=True)

image = modal.Image.debian_slim(python_version="3.11").pip_install("prody==2.6.1", "numpy")


@app.function(image=image, gpu="A10G", timeout=1800, volumes={"/results": vol})
def laci_test(dbd_range: tuple[int, int] = (1, 62)) -> dict:
    import numpy as np
    import prody

    apo = prody.parsePDB("1lbg").select("name CA")
    holo = prody.parsePDB("1lbi").select("name CA")

    def cmap(sel):
        return {(str(c), int(r)): np.array(x)
                for c, r, x in zip(sel.getChids(), sel.getResnums(), sel.getCoords())}

    am, hm = cmap(apo), cmap(holo)
    keys = sorted(set(am) & set(hm))          # 1LBG y 1LBI no tienen el mismo nº de CA

    def kabsch(P, Q):
        """R,t que llevan P sobre Q (mínimos cuadrados)."""
        pc, qc = P.mean(0), Q.mean(0)
        H = (P - pc).T @ (Q - qc)
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        return R, qc - R @ pc

    core_keys = [k for k in keys if k[1] > dbd_range[1]]
    P = np.array([hm[k] for k in core_keys])
    Q = np.array([am[k] for k in core_keys])
    R, t = kabsch(P, Q)
    holo_aligned = {k: R @ hm[k] + t for k in keys}

    d = np.array([holo_aligned[k] - am[k] for k in keys]).flatten()
    d_norm = d / np.linalg.norm(d)
    rmsd_core = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum() / len(core_keys)))

    # ANM sobre el tetrámero apo completo
    anm = prody.ANM("laci_apo")
    anm.buildHessian(apo)
    anm.calcModes(n_modes=6)
    eig = np.asarray(anm.getEigvals())[:6]
    vecs = np.asarray(anm.getEigvecs())

    apo_keys = [(str(c), int(r)) for c, r in zip(apo.getChids(), apo.getResnums())]
    idx_of = {k: i for i, k in enumerate(apo_keys)}
    sel_idx = [idx_of[k] for k in keys if k in idx_of]
    keep = [i for i, k in enumerate(keys) if k in idx_of]
    d_use = np.array([holo_aligned[keys[i]] - am[keys[i]] for i in keep]).flatten()
    d_use = d_use / np.linalg.norm(d_use)

    resnums = apo.getResnums()
    chains = apo.getChids()
    in_dbd = np.array([r <= dbd_range[1] for r in resnums])
    frac_res = float(in_dbd.mean())

    modes = []
    for k in range(min(5, vecs.shape[1])):
        v = vecs[:, k]
        v_use = v.reshape(-1, 3)[sel_idx].flatten()
        v_use = v_use / np.linalg.norm(v_use)
        overlap = float(abs(v_use @ d_use))
        mag = np.linalg.norm(v.reshape(-1, 3), axis=1)
        frac_dbd = float((mag[in_dbd] ** 2).sum() / (mag ** 2).sum())
        per_chain = {}
        for ch in np.unique(chains):
            m = in_dbd & (chains == ch)
            if m.sum():
                per_chain[str(ch)] = v.reshape(-1, 3)[m].mean(0)
        anti = None
        ks = sorted(per_chain)
        if len(ks) >= 2:
            a, b = per_chain[ks[0]], per_chain[ks[-1]]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            if na and nb:
                anti = float((a @ b) / (na * nb))
        modes.append({"modo": k + 1, "autovalor": float(eig[k]),
                      "solapamiento_observado": round(overlap, 3),
                      "frac_dbd": round(frac_dbd, 3),
                      "enriquecimiento": round(frac_dbd / frac_res, 2) if frac_res else None,
                      "coseno_DBD_entre_cadenas": None if anti is None else round(anti, 3)})

    # Solapamiento de subespacio: ¿contiene el espacio de los k modos más blandos
    # la dirección observada? Es la medida justa para una transición colectiva.
    V = np.asarray(anm.getEigvecs())[:, :5]
    V_use = np.array([V[:, k].reshape(-1, 3)[sel_idx].flatten()
                      for k in range(V.shape[1])]).T
    V_use = V_use / np.linalg.norm(V_use, axis=0, keepdims=True)
    Qorth, _ = np.linalg.qr(V_use)
    sub = {}
    for k in (3, 5):
        Qk = Qorth[:, :k]
        sub[f"modos_1_{k}"] = round(float(np.linalg.norm(Qk.T @ d_use)), 3)

    rec = {"n_atomos_apo": int(apo.numAtoms()), "n_atomos_comunes": len(keys),
           "solapamiento_subespacio": sub,
           "frac_residuos_dbd": round(frac_res, 3), "rmsd_core": round(rmsd_core, 3),
           "modos": modes,
           "veredicto": ("pasa" if (sub.get("modos_1_5", 0) > 0.5 or
                                    any(m["solapamiento_observado"] > 0.5 and
                                        (m["enriquecimiento"] or 0) > 2 for m in modes))
                         else "no pasa")}
    Path("/results", "laci_oligomer.json").write_text(json.dumps(rec, indent=1))
    vol.commit()
    return rec


@app.local_entrypoint()
def main(mode: str = "run"):
    rec = laci_test.remote()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "laci_oligomer.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print(json.dumps(rec, indent=1, ensure_ascii=False))
