#!/usr/bin/env python3
"""Prepara receptores para el piloto de energía libre: protonar a pH 7 y minimizar.

Usa OpenMM + PDBFixer en Modal (contenedor CPU: la minimización de ~9k átomos es de
segundos y no necesita GPU; así no se ocupa el CPU local ni se gasta GPU).

Por cada sistema:
  1. completa residuos/átomos faltantes (PDBFixer)
  2. añade hidrógenos a pH 7
  3. minimiza energía en solvente implícito (GBSA-OBC)
  4. escribe `receptor_preparado.pdb`

Uso:
    modal run scripts/allostery/24_prepare_receptors.py
    modal run scripts/allostery/24_prepare_receptors.py --mode fetch
"""
from __future__ import annotations

import json
from pathlib import Path

import modal

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if len(_HERE.parents) > 2 else Path("/root")
PILOT = ROOT / "data" / "allostery" / "pilot"

app = modal.App("alloreg-prep")
vol = modal.Volume.from_name("alloreg-prep-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("openmm", "pdbfixer", "numpy")
    .add_local_dir(str(PILOT), remote_path="/pilot")
)


@app.function(image=image, timeout=1800, volumes={"/results": vol}, cpu=4, memory=8192)
def preparar(nombre: str, ph: float = 7.0) -> dict:
    from openmm import unit, LangevinIntegrator
    from openmm.app import (PDBFile, Modeller, ForceField, Simulation,
                            CutoffNonPeriodic)

    src = Path("/pilot") / nombre / "receptor.pdb"
    if not src.exists():
        return {"sistema": nombre, "error": "sin receptor.pdb"}

    info = {"sistema": nombre, "pasos": []}

    # 1) PDBFixer si está disponible; si no, Modeller directo
    try:
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(src))
        fixer.findMissingResidues()
        fixer.missingResidues = {}          # no inventar bucles ausentes
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        modelo = Modeller(fixer.topology, fixer.positions)
        info["pasos"].append("pdbfixer")
    except Exception as e:
        pdb = PDBFile(str(src))
        modelo = Modeller(pdb.topology, pdb.positions)
        modelo.addHydrogens(pH=ph)
        info["pasos"].append(f"modeller ({str(e)[:60]})")

    n_atomos = modelo.topology.getNumAtoms()

    # 2) minimización en solvente implícito
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    sistema = ff.createSystem(modelo.topology, nonbondedMethod=CutoffNonPeriodic,
                              nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=None)
    integrador = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    2 * unit.femtosecond)
    from openmm import Platform
    plataforma = Platform.getPlatformByName("CPU")
    sim = Simulation(modelo.topology, sistema, integrador, plataforma,
                     {"Threads": "4"})
    sim.context.setPositions(modelo.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=2000)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    e1 = st.getPotentialEnergy()

    salida = Path("/results") / nombre
    salida.mkdir(parents=True, exist_ok=True)
    out = salida / "receptor_preparado.pdb"
    with open(out, "w") as fh:
        PDBFile.writeFile(modelo.topology, st.getPositions(), fh)

    info.update({"n_atomos": n_atomos,
                 "energia_inicial_kj": round(e0.value_in_unit(unit.kilojoule_per_mole), 1),
                 "energia_final_kj": round(e1.value_in_unit(unit.kilojoule_per_mole), 1),
                 "delta_kj": round((e1 - e0).value_in_unit(unit.kilojoule_per_mole), 1)})
    vol.commit()
    return info


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": vol}, timeout=300)
def leer(nombre: str) -> str:
    return (Path("/results") / nombre / "receptor_preparado.pdb").read_text()


@app.local_entrypoint()
def main(mode: str = "run", sistemas: str = "LacI_IPTG,PurR_hipoxantina"):
    nombres = [s.strip() for s in sistemas.split(",") if s.strip()]
    if mode == "fetch":
        for n in nombres:
            try:
                txt = leer.remote(n)
                (PILOT / n / "receptor_preparado.pdb").write_text(txt)
                print(f"{n}: {len(txt)} bytes")
            except Exception as e:
                print(f"{n}: error {str(e)[:80]}")
        return
    for n in nombres:
        print(json.dumps(preparar.remote(n), ensure_ascii=False))
