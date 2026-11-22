#!/usr/bin/env python3
"""Energía de interacción receptor–ligando con OpenMM + OpenFF (entorno conda).

Requiere el entorno creado con micromamba (`.mamba_root/envs/.mamba_mm`), que sí tiene
`openff-toolkit` y `openmmforcefields` — no disponibles por pip en este equipo.

Qué calcula (en **vacío**, para evitar el artefacto de solvatación al desacoplar):
    E_int = E(complejo) − E(receptor con ligando desacoplado) − E(ligando desacoplado)
sobre el complejo minimizado. Es una energía de interacción MM, **no** una energía
libre: no incluye entropía ni solvatación.

Uso (con el intérprete del entorno conda):
    .mamba_root/envs/.mamba_mm312/bin/python scripts/allostery/27_interaction_mm.py \
        receptor_preparado.pdb ligand.sdf --steps 50
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from openff.toolkit import Molecule
from openmm import LangevinIntegrator, NonbondedForce, Platform, unit
from openmm.app import (CutoffNonPeriodic, Modeller, PDBFile, Simulation)
from openmmforcefields.generators import SystemGenerator


def build(receptor: Path, ligando: Path, ff_small: str):
    from rdkit import Chem

    rec = PDBFile(str(receptor))
    n_rec = rec.topology.getNumAtoms()
    modeller = Modeller(rec.topology, rec.positions)

    mol = Molecule.from_file(str(ligando))
    # el ligando se añade como residuo propio al modelo
    lig_top = mol.to_topology().to_openmm()
    positions = mol.conformers[0].to("nanometer").magnitude * unit.nanometer
    modeller.add(lig_top, positions)
    n_tot = modeller.topology.getNumAtoms()

    gen = SystemGenerator(
        forcefields=["amber14-all.xml"],
        small_molecule_forcefield=ff_small,
        molecules=[mol],
        nonperiodic_forcefield_kwargs={"nonbondedMethod": CutoffNonPeriodic},
    )
    system = gen.create_system(modeller.topology)
    return modeller, system, n_rec, n_tot, mol


def decouple(system, idx, zero_charge=True):
    import copy

    s = copy.deepcopy(system)
    for f in s.getForces():
        if isinstance(f, NonbondedForce):
            for i in idx:
                q, sigma, eps = f.getParticleParameters(i)
                f.setParticleParameters(i, 0.0 if zero_charge else q,
                                        sigma, 0.0)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor")
    ap.add_argument("ligando")
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--forcefield", default="openff-2.1.1")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.perf_counter()
    modeller, system, n_rec, n_tot, mol = build(
        Path(args.receptor), Path(args.ligando), args.forcefield)
    t_build = time.perf_counter() - t0

    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtosecond)
    plat = Platform.getPlatformByName("CPU")
    sim = Simulation(modeller.topology, system, integ, plat, {"Threads": "4"})
    sim.context.setPositions(modeller.positions)
    t1 = time.perf_counter()
    sim.minimizeEnergy(maxIterations=args.steps)
    t_min = time.perf_counter() - t1
    st = sim.context.getState(getEnergy=True, getPositions=True)
    pos = st.getPositions()
    e_tot = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    idx_lig = list(range(n_rec, n_tot))

    # Energía de interacción directa: suma de pares ligando–receptor con los
    # parámetros del propio NonbondedForce (evita artefactos de desacoplar).
    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    pos_nm = np.array([[p.x, p.y, p.z] for p in pos.value_in_unit(unit.nanometer)])
    COUL = 138.9354576  # kJ·nm/mol/e^2
    e_int = 0.0
    e_lj = 0.0
    e_coul = 0.0
    for i in idx_lig:
        qi, si, ei = nb.getParticleParameters(i)
        qi, si, ei = qi.value_in_unit(unit.elementary_charge), \
            si.value_in_unit(unit.nanometer), ei.value_in_unit(unit.kilojoule_per_mole)
        for k in range(n_rec):
            qk, sk, ek = nb.getParticleParameters(k)
            qk, sk, ek = qk.value_in_unit(unit.elementary_charge), \
                sk.value_in_unit(unit.nanometer), ek.value_in_unit(unit.kilojoule_per_mole)
            if ei == 0.0 and ek == 0.0 and qi == 0.0 and qk == 0.0:
                continue
            r = float(np.linalg.norm(pos_nm[i] - pos_nm[k]))
            if r < 1e-6 or r > 1.2:
                continue
            if ei and ek:
                s6 = (si * sk / r ** 2) ** 3
                e_lj += 4.0 * (ei * ek) ** 0.5 * (s6 * s6 - s6)
            if qi and qk:
                e_coul += COUL * qi * qk / r
    e_int = e_lj + e_coul

    res = {"receptor": args.receptor, "ligando": args.ligando,
           "forcefield_ligando": args.forcefield,
           "atomos_receptor": n_rec, "atomos_totales": n_tot,
           "smiles_ligando": mol.to_smiles(),
           "e_complejo_kj": round(e_tot, 1),
           "e_interaccion_kj": round(e_int, 1),
           "e_lj_kj": round(e_lj, 1),
           "e_coulomb_kj": round(e_coul, 1),
           "t_construccion_s": round(t_build, 1),
           "t_minimizacion_s": round(t_min, 1),
           "nota": "interacción MM en vacío; NO es energía libre ni ΔΔG"}
    if args.out:
        Path(args.out).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
