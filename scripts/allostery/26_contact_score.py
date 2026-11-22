#!/usr/bin/env python3
"""Puntuación de contacto receptor–ligando, local y sin campo de fuerzas.

Sustituye a la energía de intercción MM cuando no hay parámetros de ligando
(OpenFF/GAFF no están disponibles en este entorno). **No es una energía**: mide si el
ligando está efectivamente acomodado en el bolsillo.

Métricas por par (receptor preparado, ligando SDF):
  - contactos:               átomos pesados del ligando con algún átomo del receptor < 4,5 Å
  - fracción enterrada:      fracción de átomos del ligando con ≥3 contactos
  - distancia mínima:        control de choques (< 2,0 Å = mal)
  - puentes de H plausibles: pares N/O a 2,5-3,5 Å
  - complementariedad:       contactos / átomos del ligando

Se valida comparando la pose cristalográfica contra la misma molécula desplazada
(control negativo): la pose real debe enterrarse mucho más.

Uso:
    python scripts/allostery/26_contact_score.py receptor.pdb ligand.sdf
    python scripts/allostery/26_contact_score.py receptor.pdb ligand.sdf --decoy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

CONTACT = 4.5
HB_MIN, HB_MAX = 2.5, 3.5
CLASH = 2.0


def read_pdb_heavy(path: Path):
    xyz, el = [], []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        e = (line[76:78].strip() or line[12:14].strip()).upper()
        if e == "H":
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        el.append(e)
    return np.array(xyz), el


def read_sdf_heavy(path: Path):
    from rdkit import Chem

    mol = Chem.MolFromMolFile(str(path), removeHs=True)
    if mol is None:
        raise SystemExit("SDF ilegible")
    xyz = np.array(mol.GetConformer().GetPositions())
    el = [a.GetSymbol().upper() for a in mol.GetAtoms()]
    return xyz, el


def score(rec_xyz, rec_el, lig_xyz, lig_el) -> dict:
    d = np.linalg.norm(lig_xyz[:, None, :] - rec_xyz[None, :, :], axis=2)
    n_cont = int((d.min(1) < CONTACT).sum())
    n_multi = int(((d < CONTACT).sum(1) >= 3).sum())
    clash = int((d < CLASH).sum())
    hb = 0
    for i, e in enumerate(lig_el):
        if e not in ("N", "O"):
            continue
        for j, f in enumerate(rec_el):
            if f in ("N", "O") and HB_MIN <= d[i, j] <= HB_MAX:
                hb += 1
    n = len(lig_xyz)
    return {"atomos_ligando": n,
            "contactos": n_cont,
            "fraccion_contacto": round(n_cont / n, 3),
            "enterrados_3contactos": n_multi,
            "fraccion_enterrada": round(n_multi / n, 3),
            "dist_min_A": round(float(d.min()), 2),
            "choques_menos_2A": clash,
            "puentes_H_plausibles": hb,
            "complementariedad": round((n_cont + 0.5 * n_multi) / n, 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor")
    ap.add_argument("ligando")
    ap.add_argument("--decoy", action="store_true",
                    help="evalúa también una copia desplazada 25 Å (control negativo)")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    rec_xyz, rec_el = read_pdb_heavy(Path(args.receptor))
    lig_xyz, lig_el = read_sdf_heavy(Path(args.ligando))
    res = {"receptor": args.receptor, "ligando": args.ligando,
           "pose": score(rec_xyz, rec_el, lig_xyz, lig_el)}
    if args.decoy:
        res["decoy_desplazado"] = score(rec_xyz, rec_el, lig_xyz + 25.0, lig_el)
        a = res["pose"]["complementariedad"]
        b = res["decoy_desplazado"]["complementariedad"]
        res["separa"] = bool(a > b)
    if args.out:
        Path(args.out).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
