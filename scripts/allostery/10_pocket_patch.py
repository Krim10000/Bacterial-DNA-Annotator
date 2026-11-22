#!/usr/bin/env python3
"""F2f — ¿el modo bisagra deforma el bolsillo anotado?

Tras descartar P2Rank sobre conformeros extrapolados, se usa una medida robusta:
la **deformación interna** del parche del bolsillo (cambio en la matriz de distancias
entre sus residuos) a lo largo del modo bisagra, comparada contra 200 parches
aleatorios de igual tamaño y compactados igual (k vecinos de un residuo al azar).

Motivo: las distancias locales sobreviven a la extrapolación lineal; una superficie
entera no.

Casos: los TFs con sitio de unión anotado en UniProt fuera del DBD
(PurR, TreR, BirA, NemR, FrmR). Estructura: modelo AlphaFold (monómero — limitación).

    bash scripts/allostery/modal_guard.sh status
    modal run scripts/allostery/10_pocket_patch.py
    bash scripts/allostery/modal_guard.sh stop-all
"""
from __future__ import annotations

import json
from pathlib import Path

import modal

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if len(_HERE.parents) > 2 else Path("/root")
OUT = ROOT / "data" / "allostery"
AF = OUT / "work" / "af"

app = modal.App("alloreg-patch")
vol = modal.Volume.from_name("alloreg-anm-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("prody==2.6.1", "numpy")
    .add_local_dir(str(AF), remote_path="/st")
)


@app.function(image=image, gpu="A10G", timeout=1200, volumes={"/results": vol})
def patch(job: dict) -> dict:
    import numpy as np
    import prody

    tf = job["tf"]
    pdb = Path("/st") / f"{tf}.pdb"
    if not pdb.exists():
        return {"tf": tf, "error": "sin estructura"}
    prot = prody.parsePDB(str(pdb))
    ca = prot.select("name CA")
    anm = prody.ANM(tf)
    anm.buildHessian(ca)
    anm.calcModes(n_modes=5)
    eig = np.asarray(anm.getEigvals())[:5]
    vecs = np.asarray(anm.getEigvecs())
    resnums = ca.getResnums()
    coords = ca.getCoords()
    n = len(resnums)

    # modo con mayor fracción en el DBD (definido por el rango Pfam)
    lo, hi = job.get("dbd_from"), job.get("dbd_to")
    dbd = np.array([(r >= lo) & (r <= hi) for r in resnums]) if lo else np.zeros(n, bool)
    fracs = []
    for k in range(vecs.shape[1]):
        v = vecs[:, k].reshape(-1, 3)
        m2 = (v ** 2).sum(1)
        fracs.append(float(m2[dbd].sum() / m2.sum()) if dbd.any() else 0.0)
    k = int(np.argmax(fracs))
    v = vecs[:, k].reshape(-1, 3)
    v = v / np.linalg.norm(v, axis=1).max() * 1.5
    plus, minus = coords + v, coords - v

    # residuos del bolsillo anotado + vecinos (<6 Å) para incluir las paredes
    site = [s for s in job.get("sites", []) if any(abs(r - s) <= 0 for r in resnums)]
    if not site:
        return {"tf": tf, "error": "sitio anotado no mapea a la estructura"}
    idx_site = [i for i, r in enumerate(resnums) if r in site]
    # vecinos
    d = np.linalg.norm(coords[idx_site][:, None, :] - coords[None, :, :], axis=2)
    neigh = set(idx_site) | set(np.where(d.min(0) < 6.0)[0].tolist())
    idx_patch = sorted(neigh)
    kpatch = len(idx_patch)

    def deform(idx):
        a, b = plus[idx], minus[idx]
        da = np.linalg.norm(a[:, None] - a[None, :], axis=2)
        db = np.linalg.norm(b[:, None] - b[None, :], axis=2)
        iu = np.triu_indices(len(idx), 1)
        return float(np.abs(da - db)[iu].mean())

    real = deform(idx_patch)

    rng = np.random.default_rng(0)
    rand = []
    for _ in range(200):
        c = int(rng.integers(0, n))
        dvec = np.linalg.norm(coords - coords[c], axis=1)
        pick = np.argsort(dvec)[:kpatch]
        rand.append(deform(pick))
    rand = np.array(rand)
    z = float((real - rand.mean()) / rand.std()) if rand.std() > 0 else None

    rec = {"tf": tf, "modo": k + 1, "frac_dbd_modo": round(fracs[k], 3),
           "n_patch": kpatch, "sitios": site,
           "deform_real": round(real, 4),
           "deform_azar_media": round(float(rand.mean()), 4),
           "z": None if z is None else round(z, 2),
           "pasa": bool(z is not None and z > 2)}
    Path("/results", f"patch_{tf}.json").write_text(json.dumps(rec))
    vol.commit()
    return rec


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": vol}, timeout=300)
def read_all() -> list[dict]:
    import glob
    return [json.load(open(f)) for f in sorted(glob.glob("/results/patch_*.json"))]


@app.local_entrypoint()
def main(mode: str = "run"):
    if mode == "fetch":
        data = read_all.remote()
        (OUT / "pocket_patch.json").write_text(json.dumps(data, indent=1, ensure_ascii=False))
        for r in data:
            if "error" in r:
                print(f"{r['tf']:8s} ERROR {r['error']}")
            else:
                print(f"{r['tf']:8s} modo={r['modo']} n={r['n_patch']:3d} "
                      f"deform={r['deform_real']:.4f} azar={r['deform_azar_media']:.4f} "
                      f"z={r['z']} pasa={r['pasa']}")
        return

    sites = json.load(open(OUT / "work" / "uniprot_sites.json"))
    dbd = json.load(open(OUT / "work" / "dbd_ranges.json"))
    jobs = []
    for o in sites:
        if not o.get("features"):
            continue
        pts = [f["start"] for f in o["features"]
               if f["type"] in ("Binding site", "Site", "Metal binding") and f["start"]]
        if not pts:
            continue
        jobs.append({"tf": o["tf"], "sites": pts, **dbd.get(o["tf"], {})})
    print(f"lanzando {len(jobs)} TFs con sitio anotado", flush=True)
    for r in patch.map(jobs):
        print(r)


if __name__ == "__main__":
    pass
