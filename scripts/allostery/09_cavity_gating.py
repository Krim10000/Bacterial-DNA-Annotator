#!/usr/bin/env python3
"""F2e — gating de cavidad: ¿la bisagra abre/cierra un bolsillo?

El detector de bisagra tiene especificidad 5/8: dispara en λ cI, LexA y 434 Cro, que
tienen bisagra de linker pero **no** bolsillo de efector. La hipótesis es que lo que
distingue un interruptor alostérico es que la bisagra **gatee una cavidad**.

Diseño:
  1. ANM; elegir el modo con z_dbd > 2 y coseno antiparalelo < −0,5 (el modo bisagra).
  2. Generar dos conformeros x(±) = x0 ± s·v con desplazamiento máximo de 2,5 Å.
  3. P2Rank sobre cada conformero.
  4. Métrica: bolsillos presentes en + y ausentes en − (y viceversa), y diferencia de
     `sas_points` del bolsillo mayor.

Predicción: LacI (interruptor real) muestra un bolsillo que se abre/cierra; λ cI, LexA
y 434 Cro no.

USO DE MODAL: 1 contenedor A10G, 5 jobs. Cerrar al terminar con
    bash scripts/allostery/modal_guard.sh stop-all

    modal run scripts/allostery/09_cavity_gating.py
"""
from __future__ import annotations

import json
from pathlib import Path

import modal

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if len(_HERE.parents) > 2 else Path("/root")
OUT = ROOT / "data" / "allostery"

app = modal.App("alloreg-gating")
vol = modal.Volume.from_name("alloreg-anm-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("default-jre-headless", "curl", "tar")
    .run_commands(
        "curl -sL -o /tmp/p2rank.tgz "
        "https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz",
        "tar -xzf /tmp/p2rank.tgz -C /opt",
        "ls /opt/p2rank*",
    )
    .pip_install("prody==2.6.1", "numpy")
)

CASES = [
    {"name": "LacI", "pdb": "1lbg", "rol": "positivo"},
    {"name": "lambda_cI", "pdb": "1lmb", "rol": "falso_positivo"},
    {"name": "LexA", "pdb": "1jhh", "rol": "falso_positivo"},
    {"name": "434_Cro", "pdb": "3cro", "rol": "falso_positivo"},
    {"name": "Fis", "pdb": "1f36", "rol": "negativo"},
]


@app.function(image=image, gpu="A10G", timeout=1800, volumes={"/results": vol})
def cavity(case: dict) -> dict:
    import csv as _csv
    import subprocess

    import numpy as np
    import prody

    prot = prody.parsePDB(case["pdb"]).select("protein")
    ca = prot.select("name CA")
    anm = prody.ANM(case["name"])
    anm.buildHessian(ca)
    anm.calcModes(n_modes=5)
    vecs = np.asarray(anm.getEigvecs())
    resnums, chids = ca.getResnums(), np.array([str(c) for c in ca.getChids()])
    L = min(62, len(resnums) - 10)
    real = np.zeros(len(resnums), dtype=bool)
    for ch in np.unique(chids):
        m = chids == ch
        real |= (resnums >= resnums[m].min()) & (resnums < resnums[m].min() + L) & m

    # elegir el modo bisagra: DBD-dominado y antiparalelo entre cadenas
    k_best = 0
    for k in range(vecs.shape[1]):
        v = vecs[:, k].reshape(-1, 3)
        mag2 = (v ** 2).sum(1)
        frac = float(mag2[real].sum() / mag2.sum())
        per = {}
        for ch in np.unique(chids):
            m = real & (chids == ch)
            if m.sum():
                per[ch] = v[m].mean(0)
        anti = None
        if len(per) >= 2:
            ks = sorted(per)
            a, b = per[ks[0]], per[ks[-1]]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            if na and nb:
                anti = float((a @ b) / (na * nb))
        if anti is not None and anti < -0.5 and frac > 0.3:
            k_best = k
            break

    v = vecs[:, k_best].reshape(-1, 3)
    v = v / np.linalg.norm(v, axis=1).max() * 2.5   # 2,5 Å de desplazamiento máximo
    vmap = {(str(ch), int(rn)): v[i] for i, (ch, rn) in enumerate(zip(ca.getChids(), ca.getResnums()))}

    # diagnóstico: ¿P2Rank encuentra bolsillos en la estructura SIN perturbar?
    base = prot.getCoords().copy()
    subprocess.run(["/opt/p2rank_2.5.1/prank", "predict", "-f", str(prody.writePDB(
        "/tmp/orig.pdb", prot)), "-o", "/tmp/p2rank_orig"], capture_output=True, text=True)
    base_pockets = len(list(Path("/tmp/p2rank_orig").glob("*_predictions.csv")))
    allkeys = [(str(ch), int(rn)) for ch, rn in zip(prot.getChids(), prot.getResnums())]
    out = {}
    for sign, tag in ((+1, "plus"), (-1, "minus")):
        x = base + np.array([sign * vmap[k] for k in allkeys])
        moved = prot.copy()
        moved.setCoords(x)
        pdb = Path(f"/tmp/{case['name']}_{tag}.pdb")
        prody.writePDB(str(pdb), moved)
        odir = Path(f"/tmp/p2rank_{tag}_{case['name']}")
        pr = subprocess.run(["/opt/p2rank_2.5.1/prank", "predict", "-f", str(pdb),
                             "-o", str(odir)], capture_output=True, text=True)
        csvs = list(odir.glob("*_predictions.csv"))
        pockets = []
        header = []
        if csvs:
            with open(csvs[0]) as fh:
                rd = _csv.DictReader(fh)
                header = [h.strip() for h in (rd.fieldnames or [])]
                for row in rd:
                    row = {k.strip(): v for k, v in row.items()}
                    pockets.append({"rank": int(row["rank"]), "score": float(row["score"]),
                                    "prob": float(row["probability"]),
                                    "sas": int(float(row["sas_points"])),
                                    "x": float(row["center_x"]), "y": float(row["center_y"]),
                                    "z": float(row["center_z"])})
        out[tag] = pockets
        out[f"header_{tag}"] = header
        out[f"rc_{tag}"] = pr.returncode
        out[f"err_{tag}"] = (pr.stderr or pr.stdout)[-600:]

    def unmatched(a, b):
        res = []
        for p in a:
            d = [((p["x"] - q["x"]) ** 2 + (p["y"] - q["y"]) ** 2 + (p["z"] - q["z"]) ** 2) ** 0.5
                 for q in b]
            if not d or min(d) > 6.0:
                res.append(p)
        return res

    plus, minus = out.get("plus", []), out.get("minus", [])
    new_plus, new_minus = unmatched(plus, minus), unmatched(minus, plus)
    big_plus = max([p["sas"] for p in plus], default=0)
    big_minus = max([p["sas"] for p in minus], default=0)
    rec = {"name": case["name"], "rol": case["rol"], "modo": k_best + 1,
           "n_pockets_plus": len(plus), "n_pockets_minus": len(minus),
           "nuevos_plus": len(new_plus), "nuevos_minus": len(new_minus),
           "sas_max_plus": big_plus, "sas_max_minus": big_minus,
           "delta_sas": big_plus - big_minus,
           "nuevos_detalle": [{"sas": p["sas"], "prob": round(p["prob"], 2)} for p in new_plus],
           "gating": bool(new_plus or new_minus),
           "header": out.get("header_plus"),
           "raw_first": (plus[0] if plus else None),
           "p2rank_rc": out.get("rc_plus"), "p2rank_err": out.get("err_plus"),
           "baseline_p2rank_ok": base_pockets}
    Path("/results", f"gating_{case['name']}.json").write_text(json.dumps(rec))
    vol.commit()
    return rec


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": vol}, timeout=300)
def read_all() -> list[dict]:
    import glob
    return [json.load(open(f)) for f in sorted(glob.glob("/results/gating_*.json"))]


@app.local_entrypoint()
def main(mode: str = "run"):
    if mode == "fetch":
        data = read_all.remote()
        (OUT / "cavity_gating.json").write_text(json.dumps(data, indent=1, ensure_ascii=False))
        for r in data:
            print(f"{r['name']:12s} {r['rol']:16s} gating={str(r['gating']):5s} "
                  f"nuevos(+/-)={r['nuevos_plus']}/{r['nuevos_minus']} "
                  f"sas_max {r['sas_max_plus']} vs {r['sas_max_minus']} "
                  f"delta={r['delta_sas']}")
        return
    cases = CASES if mode == "run" else [c for c in CASES if c["name"] == mode]
    for r in cavity.map(cases):
        print(r["name"], r.get("gating"), r.get("header"))
