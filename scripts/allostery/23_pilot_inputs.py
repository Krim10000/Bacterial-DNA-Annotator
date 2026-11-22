#!/usr/bin/env python3
"""Prepara los insumos del piloto de energía libre (LacI, PurR, TrpR).

Para cada sistema:
  - descarga la estructura holo de PDB (con el efector cristalizado)
  - separa receptor (proteína) y ligando (HETATM del efector)
  - escribe el ligando en SDF 3D con la química del componente CCD de PDB
    (orden de enlaces correcto) y las coordenadas de la pose cristalográfica
  - escribe metadatos (PDB, comp_id, SMILES de referencia, cadenas)

Salida: data/allostery/pilot/<sistema>/{receptor.pdb,ligand.sdf,metadata.json}
        data/allostery/pilot/README.md  (instrucciones OpenFE / NeuroSnap)

Uso: python scripts/allostery/23_pilot_inputs.py
"""
from __future__ import annotations

import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery" / "pilot"

SYSTEMS = [
    {"nombre": "LacI_IPTG", "pdb": "1lbh", "comp": "IPT",
     "proteina": "LacI (represor del operón lac)", "efector": "IPTG (inductor gratuito)",
     "modo": "inductor", "smiles_ref": "CC(C)S[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"},
    {"nombre": "PurR_hipoxantina", "pdb": "1pnr", "comp": "HPA",
     "proteina": "PurR (represor de purinas)", "efector": "hipoxantina",
     "modo": "correpresor", "smiles_ref": "O=c1[nH]cnc2nc[nH]c12"},
    {"nombre": "TrpR_Trp", "pdb": "1tro", "comp": "TRP",
     "proteina": "TrpR (represor del operón trp)", "efector": "L-triptófano",
     "modo": "correpresor",
     "smiles_ref": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"},
]


def fetch(url: str) -> bytes:
    req = Request(url, headers={"User-Agent": "alloreg"})
    with urlopen(req, timeout=60) as fh:
        return fh.read()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hechos = []
    for s in SYSTEMS:
        d = OUT / s["nombre"]
        d.mkdir(parents=True, exist_ok=True)
        pdb_txt = fetch(f"https://files.rcsb.org/download/{s['pdb'].upper()}.pdb").decode()
        (d / f"original_{s['pdb']}.pdb").write_text(pdb_txt)

        rec, het = [], []
        for line in pdb_txt.splitlines():
            if line.startswith(("ATOM", "TER")):
                rec.append(line)
            elif line.startswith("HETATM"):
                resn = line[17:20].strip()
                if resn == s["comp"]:
                    het.append(line)
                # aguas e iones se descartan
            elif line.startswith("END"):
                break
        # el cristal puede traer varias copias del ligando (una por subunidad):
        # se conserva la primera (cadena + número de residuo)
        grupos = {}
        for line in het:
            key = (line[21], line[22:27])
            grupos.setdefault(key, []).append(line)
        first = sorted(grupos)[0] if grupos else None
        het = grupos[first] if first else []
        (d / "receptor.pdb").write_text("\n".join(rec) + "\nEND\n")
        (d / "ligand_hetatm.pdb").write_text("\n".join(het) + "\nEND\n")
        copias = {f"{k[0]}:{k[1].strip()}": len(v) for k, v in grupos.items()}

        # química correcta desde el CCD, coordenadas desde el cristal
        ok_ccd = True
        modo_ccd = "asignacion_directa"
        try:
            ccd = fetch(f"https://files.rcsb.org/ligands/download/{s['comp']}_ideal.sdf").decode()
            ref = Chem.MolFromMolBlock(ccd.split("$$$$")[0], sanitize=True, removeHs=False)
            raw = Chem.MolFromPDBBlock("\n".join(het), sanitize=False, removeHs=False)
            raw = Chem.RemoveHs(raw, sanitize=False)
            # vía preferente: el CCD y el HETATM suelen compartir orden de átomos.
            ideal_heavy = Chem.RemoveHs(ref)
            cry = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het])
            sym_ccd = [a.GetSymbol() for a in ideal_heavy.GetAtoms()]
            sym_pdb = [(l[76:78].strip() or l[12:14].strip()) for l in het]
            if sym_ccd == sym_pdb:
                conf = Chem.Conformer(ideal_heavy.GetNumAtoms())
                for i, xyz in enumerate(cry):
                    conf.SetAtomPosition(i, [float(v) for v in xyz])
                ideal_heavy.RemoveAllConformers()
                ideal_heavy.AddConformer(conf)
                mol = Chem.AddHs(ideal_heavy, addCoords=True)
                modo_ccd = "coordenadas_cristal_sobre_ccd"
            elif len(sym_pdb) < len(sym_ccd):
                # faltan átomos en el cristal (p. ej. carboxilato desordenado):
                # no se inventa la pose; se entrega la molécula ideal del CCD
                raise ValueError(f"cristal incompleto ({len(sym_pdb)} de {len(sym_ccd)} átomos)")
            elif ref.GetNumHeavyAtoms() == raw.GetNumHeavyAtoms():
                mol = AllChem.AssignBondOrdersFromTemplate(ref, raw)
            else:
                # el cristal trae algún átomo desordenado/ausente: se alinea la molécula
                # ideal (química correcta) sobre los átomos presentes
                match = ref.GetSubstructMatch(raw)
                if not match:
                    raise ValueError("sin subestructura común")
                mol = Chem.Mol(ref)
                amap = [(i, m) for i, m in enumerate(match)]
                AllChem.AlignMol(mol, raw, atomMap=amap)
                modo_ccd = f"alineado_ideal_sobre_cristal ({ref.GetNumHeavyAtoms()} vs {raw.GetNumHeavyAtoms()} átomos)"
        except Exception as e:
            ok_ccd = False
            modo_ccd = f"pendiente_alinear ({str(e)[:40]})"
            mol = None
            # se guarda además la molécula ideal del CCD (química correcta, sin pose)
            try:
                ccd = fetch(f"https://files.rcsb.org/ligands/download/{s['comp']}_ideal.sdf").decode()
                ref2 = Chem.MolFromMolBlock(ccd.split("$$$$")[0], sanitize=True, removeHs=False)
                if ref2:
                    w2 = Chem.SDWriter(str(d / "ligand_ideal_ccd.sdf"))
                    w2.write(ref2)
                    w2.close()
                    modo_ccd += " + ligand_ideal_ccd.sdf"
            except Exception:
                pass
        n_heavy = mol.GetNumHeavyAtoms() if mol is not None else 0
        if mol is not None:
            mol.SetProp("_Name", f"{s['comp']}_{s['nombre']}")
            w = Chem.SDWriter(str(d / "ligand.sdf"))
            w.write(mol)
            w.close()

        (d / "metadata.json").write_text(json.dumps({
            **s, "pdb_file": f"original_{s['pdb']}.pdb",
            "natoms_receptor": len(rec), "natoms_ligando": len(het),
            "n_heavy_ligando": n_heavy, "ccd_usado": ok_ccd,
            "copias_ligando_en_cristal": len(copias), "copia_elegida": str(first),
            "modo_ccd": modo_ccd,
            "cadenas": sorted({l[21] for l in rec if l.startswith("ATOM")}),
            "pendiente": ["verificar protonación del ligando a pH 7",
                          "verificar pose en el bolsillo tras asignar órdenes de enlace",
                          "decidir si se mantiene el oligómero completo o solo el dímero"],
        }, indent=1, ensure_ascii=False))
        hechos.append(s["nombre"])
        print(f"  {s['nombre']}: {len(rec)} átomos receptor · {n_heavy} pesados ligando · ccd={ok_ccd}")

    (OUT / "README.md").write_text(f"""# Insumos del piloto de energía libre

Sistemas preparados: {", ".join(hechos)}

Cada carpeta contiene:

| archivo | contenido |
|---|---|
| `original_<pdb>.pdb` | estructura depositada completa |
| `receptor.pdb` | solo proteína (sin aguas, iones ni el ligando) |
| `ligand_hetatm.pdb` | el efector tal como está cristalizado |
| `ligand.sdf` | ligando 3D con la química del CCD y la pose del cristal |
| `metadata.json` | PDB, componente, cadenas, pendientes |

## Qué está resuelto

- Estructura holo real para los tres (1LBH/IPTG, 1PNR/hipoxantina, 1TRO/Trp).
- Pose cristalográfica del efector (no es docking: es la pose medida).
- Órdenes de enlace desde el componente químico de PDB.

## Qué falta antes de enviar

1. **Protonación a pH 7** del ligando y del receptor (PDBFixer/OpenMM). El SDF viene
   con la protonación del CCD, que puede no ser la fisiológica.
2. **Decidir la unidad funcional**: LacI es tetrámero y el bolsillo está en la interfaz.
   Para el piloto conviene quedarse con el dímero que contiene el ligando.
3. **Minimización** del receptor (los H de la estructura depositada faltan o están mal).
4. Elegir la variante de protocolo: ABFE (una pose) o RBFE (comparar dos ligandos).

## Comandos

**NeuroSnap** (OpenFE ABFE): subir `receptor.pdb` + `ligand.sdf`. El servicio pide
"a prepared protein structure and one positioned SDF ligand" — es exactamente lo que
hay aquí, tras protonar.

**Local (OpenFE CLI)**:
```bash
# requiere: conda env con openfe, openmm, openmmtools
openfe plan-rbfe-network -M ligand.sdf -p receptor.pdb -o network  # solo si hay >=2 ligandos
python -m openfe.protocols.openmm_afe   # ABFE vía API, no hay CLI estable
```

**Modal (OpenMM propio)**: instalar `openmm[cuda12]` y correr MM/GBSA primero;
ABFE con `openmmtools` + `perses` es un proyecto aparte.
""")
    print(json.dumps({"sistemas": hechos}, ensure_ascii=False))


if __name__ == "__main__":
    main()
