#!/usr/bin/env python3
"""Transfiere ligando/cofactor de un PDB holo a un PDB ESMFold.

Superpone CA por alineamiento Foldseek y aplica la transformación a HETATM.
No es docking: es transferencia estructural; mide contactos y clashes.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}


def parse(path, chain=None):
    atoms, het, seq, ca = [], [], [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            ch = line[21]
            if chain and ch != chain:
                continue
            resn = line[17:20].strip().upper()
            if line.startswith("ATOM  ") and resn in AA3:
                atoms.append((ch, int(line[22:26]), resn, line[12:16].strip(),
                              np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
                if line[12:16].strip() == "CA":
                    seq.append(AA3[resn]); ca.append(atoms[-1][4])
            elif line.startswith("HETATM"):
                het.append((ch, int(line[22:26]), resn, line[12:16].strip(),
                            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return atoms, het, "".join(seq), np.array(ca) if ca else np.zeros((0, 3))


def foldseek_map(q, t):
    with tempfile.TemporaryDirectory(prefix="holo_") as td:
        out = Path(td) / "a.m8"
        subprocess.run(["foldseek", "easy-search", str(q), str(t), str(out), td,
                        "--threads", "2", "--max-seqs", "1", "-v", "0",
                        "--format-output", "qstart,tstart,qend,tend,qaln,taln"],
                       capture_output=True, text=True, timeout=120)
        p = out.read_text().splitlines()[0].split("\t")
        qs, ts, qe, te, qa, ta = p[:6]
        pairs = []
        qi, ti = int(qs) - 1, int(ts) - 1
        for x, y in zip(qa.upper(), ta.upper()):
            if y != "-":
                if x != "-":
                    pairs.append((ti, qi, x, y))
                    qi += 1
                ti += 1
            elif x != "-":
                qi += 1
        return pairs


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--holo", required=True)
    ap.add_argument("--ligands", default="NAP,TLA")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--lys", type=int, default=171, help="Número del K catalítico en el candidato")
    a = ap.parse_args()
    ligs = set(a.ligands.split(","))
    catoms, _, cseq, cca = parse(a.candidate)
    _, het, tseq, tca = parse(a.holo, chain=a.chain)
    pairs = foldseek_map(a.candidate, a.holo)
    P, Q = [], []
    for ti, qi, _, _ in pairs:
        if ti < len(tca) and qi < len(cca):
            P.append(tca[ti]); Q.append(cca[qi])
    if len(P) < 10:
        raise SystemExit("alineamiento insuficiente")
    R, t = kabsch(np.array(P), np.array(Q))
    rmsd = float(np.sqrt(np.mean([np.sum((R @ p + t - q) ** 2) for p, q in zip(P, Q)])))
    moved = []
    for ch, ri, resn, atom, xyz in het:
        if resn in ligs:
            moved.append((ch, ri, resn, atom, R @ xyz + t))
    lines = [x for x in Path(a.candidate).read_text().splitlines() if x.startswith(("ATOM  ", "HETATM"))]
    for i, (ch, ri, resn, atom, xyz) in enumerate(moved, 1):
        lines.append(f"HETATM{i:5d} {atom:<3s} {resn:>3s} L{ri:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {atom[0]:>2s}")
    Path(a.out).write_text("\n".join(lines) + "\nEND\n")
    # contacts K171 (candidate numbering)
    k = [x for x in catoms if x[2] == "LYS" and x[1] == a.lys]
    kd = []
    for ch, ri, resn, atom, xyz in moved:
        for _, _, _, katom, kxyz in k:
            kd.append((round(float(np.linalg.norm(xyz - kxyz)), 2), resn, ri, atom, katom))
    kd.sort()
    # clashes ligand-protein (<2.0 A)
    clashes = 0
    for _, _, resn, _, xyz in moved:
        for _, _, _, _, pxyz in catoms:
            if np.linalg.norm(xyz - pxyz) < 2.0:
                clashes += 1
    print({"rmsd": round(rmsd, 2), "pairs": len(P), "ligand_atoms": len(moved),
           "lys": a.lys, "contacts": kd[:5], "clashes_lt2A": clashes, "out": a.out})


if __name__ == "__main__":
    main()
