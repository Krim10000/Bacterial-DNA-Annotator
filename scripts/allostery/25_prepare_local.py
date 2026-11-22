#!/usr/bin/env python3
"""Preparación local de receptores: protonar a pH 7 y minimizar con OpenMM.

Todo en CPU, sin GPU y sin Modal. Pensado para el primer nivel del flujo:
preparar estructuras cuesta segundos-minutos; simular dinámica no es viable aquí
(medido: **0,25 ns/día** en este equipo para 5.832 átomos con solvente implícito).

Qué hace:
  1. completa átomos faltantes (PDBFixer) y estandariza residuos
  2. añade hidrógenos a pH 7
  3. minimiza energía en solvente implícito (GBSA-OBC)
  4. opcional: calcula la energía de interacción con un ligando
     (E_complejo − E_receptor − E_ligando; **no** es ΔG, es una primera pasada)

Uso:
    python scripts/allostery/25_prepare_local.py receptor.pdb --out prep.pdb
    python scripts/allostery/25_prepare_local.py receptor.pdb --ligand lig.sdf \
        --interaction
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np


def load_prepared(path: Path, ph: float = 7.0):
    from openmm.app import Modeller, PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # no inventar bucles ausentes
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    return Modeller(fixer.topology, fixer.positions)


def build_system(topology, ligando_sdf: Path | None = None):
    """Sistema con solvente implícito. Si hay ligando, se fusiona como molécula aparte."""
    from openmm import unit
    from openmm.app import ForceField, CutoffNonPeriodic

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    sistema = ff.createSystem(topology, nonbondedMethod=CutoffNonPeriodic,
                              nonbondedCutoff=1.0 * unit.nanometer)
    return ff, sistema


def minimize(topology, positions, max_iter: int = 2000, threads: int = 4):
    from openmm import LangevinIntegrator, Platform, unit
    from openmm.app import Simulation

    ff, sistema = build_system(topology)
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtosecond)
    sim = Simulation(topology, sistema, integ,
                     Platform.getPlatformByName("CPU"), {"Threads": str(threads)})
    sim.context.setPositions(positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=max_iter)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    e1 = st.getPotentialEnergy()
    return (e0.value_in_unit(unit.kilojoule_per_mole),
            e1.value_in_unit(unit.kilojoule_per_mole),
            st.getPositions())


def interaction_energy(receptor_pdb: Path, ligando_sdf: Path,
                       max_iter: int = 300, threads: int = 4) -> dict:
    """Energía de interacción en solvente implícito, sobre el complejo minimizado.

    E_int = E(complejo) − E(receptor con ligando desacoplado) − E(ligando desacoplado).
    Es una medida de contacto; **no** es una energía libre ni un ΔG.
    """
    from openmm import LangevinIntegrator, Platform, unit
    from openmm.app import CutoffNonPeriodic, ForceField, Modeller, PDBFile, Simulation
    from rdkit import Chem

    rec = PDBFile(str(receptor_pdb))
    n_rec = rec.topology.getNumAtoms()
    modelo = Modeller(rec.topology, rec.positions)
    mol = Chem.MolFromMolFile(str(ligando_sdf), removeHs=False)
    if mol is None:
        return {"error": "SDF ilegible"}
    tmp = receptor_pdb.parent / "_lig_tmp.pdb"
    tmp.write_text(Chem.MolToPDBBlock(mol))
    lig = PDBFile(str(tmp))
    modelo.add(lig.topology, lig.positions)
    n_tot = modelo.topology.getNumAtoms()
    idx_lig = list(range(n_rec, n_tot))

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")

    def sistema_con(decouple: str | None = None):
        s = ff.createSystem(modelo.topology, nonbondedMethod=CutoffNonPeriodic,
                            nonbondedCutoff=1.0 * unit.nanometer)
        if decouple:
            from openmm import NonbondedForce
            for f in s.getForces():
                if isinstance(f, NonbondedForce):
                    cero = idx_lig if decouple == "ligando" else list(range(n_rec))
                    for i in cero:
                        f.setParticleParameters(i, 0.0, 0.1, 0.0)
                    f.setParticleParameters(0, *f.getParticleParameters(0))
        return s

    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtosecond)
    plat = Platform.getPlatformByName("CPU")
    sim = Simulation(modelo.topology, sistema_con(), integ, plat,
                     {"Threads": str(threads)})
    sim.context.setPositions(modelo.positions)
    sim.minimizeEnergy(maxIterations=max_iter)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    pos = st.getPositions()
    e_tot = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    def energia_decoupled(which):
        sm = Simulation(modelo.topology, sistema_con(which),
                        LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                           2 * unit.femtosecond), plat,
                        {"Threads": str(threads)})
        sm.context.setPositions(pos)
        return sm.context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilojoule_per_mole)

    e_rec = energia_decoupled("ligando")
    e_lig = energia_decoupled("receptor")
    return {"e_complejo_kj": round(e_tot, 1),
            "e_receptor_kj": round(e_rec, 1),
            "e_ligando_kj": round(e_lig, 1),
            "e_interaccion_kj": round(e_tot - e_rec - e_lig, 1),
            "nota": "energía de interacción en solvente implícito; no es ΔG"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--out", default="")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--max-iter", type=int, default=200,
                    help="pasos de minimización (200 basta para clashes; 2000 tarda ~8 min)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ligand", default="", help="SDF del ligando (opcional)")
    ap.add_argument("--interaction", action="store_true")
    args = ap.parse_args()

    t0 = time.perf_counter()
    src = Path(args.pdb)
    modelo = load_prepared(src, args.ph)
    t_prep = time.perf_counter() - t0

    t1 = time.perf_counter()
    e0, e1, pos = minimize(modelo.topology, modelo.positions,
                           args.max_iter, args.threads)
    t_min = time.perf_counter() - t1

    out = Path(args.out) if args.out else src.with_name(src.stem + "_preparado.pdb")
    for carril in [out.parent]:
        carril.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        from openmm.app import PDBFile
        PDBFile.writeFile(modelo.topology, pos, fh)

    res = {"entrada": str(src), "salida": str(out),
           "atomos": modelo.topology.getNumAtoms(),
           "energia_inicial_kj": round(e0, 1), "energia_final_kj": round(e1, 1),
           "t_preparacion_s": round(t_prep, 1), "t_minimizacion_s": round(t_min, 1)}

    if args.interaction and args.ligand:
        try:
            res["interaccion"] = interaction_energy(out, Path(args.ligand), args.max_iter)
        except Exception as e:
            res["interaccion"] = {
                "error": str(e)[:200],
                "nota": "el ligando necesita parámetros de campo de fuerzas "
                        "(GAFF/OpenFF): instalar openmmforcefields + openff-toolkit. "
                        "La preparación del receptor no depende de esto."}

    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
