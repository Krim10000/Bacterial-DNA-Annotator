#!/usr/bin/env python3
"""F2 — co-folding con Boltz-2 en Modal (GPU).

Hipótesis operativa:
    Δ = afinidad(ligando | TF unido a ADN) − afinidad(ligando | TF libre)
    Δ < 0  ⇒ prefiere el estado unido  ⇒ correpresor
    Δ > 0  ⇒ prefiere el estado libre  ⇒ inductor

Cada (TF, ligando) genera dos predicciones: con operador (dos hebras) y sin él.
El ligando verdadero se compara contra un señuelo en el mismo TF, para que el
signo no dependa de la escala absoluta del modelo.

Flujo (pensado para sobrevivir a la desconexión del cliente):
    modal run --detach scripts/allostery/04_boltz_modal.py --mode pilot   # lanza
    modal run        scripts/allostery/04_boltz_modal.py --mode fetch     # recoge

Salida local: data/allostery/boltz_results.json
"""

# NOTA MODAL: Modal admite 10 GPU simultáneas por workspace. Antes de correr esto,
# avisa que usas Modal y al terminar cierra con:
#     bash scripts/allostery/modal_guard.sh stop-all
#     bash scripts/allostery/modal_guard.sh status   # debe salir 0 activas
# Ver docs/MODAL.md.
from __future__ import annotations

import json
import sys
from pathlib import Path

import modal

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if len(_HERE.parents) > 2 else Path("/root")
OUT = ROOT / "data" / "allostery"

app = modal.App("alloreg-boltz")
weights = modal.Volume.from_name("boltz-weights", create_if_missing=True)
results_vol = modal.Volume.from_name("alloreg-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch==2.6.0", index_url="https://download.pytorch.org/whl/cu124")
    .pip_install("boltz[cuda]")
    .env({"BOLTZ_CACHE": "/root/.boltz"})
)


def _slug(name: str) -> str:
    return name.replace("|", "__")


@app.function(image=image, gpu="A10G", timeout=1800,
              volumes={"/root/.boltz": weights, "/results": results_vol})
def run_boltz(job: dict) -> dict:
    """Corre `boltz predict` sobre un YAML y guarda afinidad + confianza en el volumen."""
    import glob
    import subprocess
    import tempfile

    work = Path(tempfile.mkdtemp())
    (work / "in.yaml").write_text(job["yaml"])
    cmd = ["boltz", "predict", str(work / "in.yaml"), "--out_dir", str(work / "out"),
           "--output_format", "pdb", "--diffusion_samples", "1"]
    if job.get("use_msa"):
        cmd.append("--use_msa_server")
    r = subprocess.run(cmd, capture_output=True, text=True)
    res = {"name": job["name"], "returncode": r.returncode,
           "stderr": r.stderr[-1500:], "meta": job["meta"]}
    for f in glob.glob(str(work / "out" / "**" / "affinity_*.json"), recursive=True):
        res["affinity"] = json.load(open(f))
    for f in glob.glob(str(work / "out" / "**" / "confidence_*.json"), recursive=True):
        c = json.load(open(f))
        res["confidence"] = {k: c.get(k) for k in
                             ("confidence_score", "ptm", "iptm", "ligand_iptm",
                              "protein_iptm", "complex_plddt")}
    Path("/results", _slug(job["name"]) + ".json").write_text(
        json.dumps(res, ensure_ascii=False))
    results_vol.commit()
    return res


@app.function(image=modal.Image.debian_slim(python_version="3.11"),
              volumes={"/results": results_vol}, timeout=300)
def read_all() -> list[dict]:
    import glob
    out = []
    for f in sorted(glob.glob("/results/*.json")):
        out.append(json.load(open(f)))
    return out


def yaml_protein_ligand(seq: str, smiles: str, dna: tuple[str, str] | None,
                        use_msa: bool = False) -> str:
    y = ["version: 1", "sequences:",
         "  - protein:", "      id: A", f"      sequence: {seq}"]
    if not use_msa:
        y.append("      msa: empty")
    y += ["  - ligand:", "      id: B", f"      smiles: '{smiles}'"]
    if dna:
        plus, minus = dna
        y += ["  - dna:", "      id: C", f"      sequence: {plus}",
              "  - dna:", "      id: D", f"      sequence: {minus}"]
    y += ["properties:", "  - affinity:", "      binder: B"]
    return "\n".join(y) + "\n"


def scramble(s: str, seed: str = "x") -> str:
    """Misma composición de bases, orden distinto: ADN no cognato de igual tamaño."""
    import hashlib
    h = int(hashlib.sha1((seed + s).encode()).hexdigest(), 16)
    bases = list(s)
    for i in range(len(bases) - 1, 0, -1):
        h = (1103515245 * h + 12345) % (2 ** 31)
        j = h % (i + 1)
        bases[i], bases[j] = bases[j], bases[i]
    out = "".join(bases)
    return out if out != s else out[::-1]


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def build_jobs(mode: str) -> list[dict]:
    import importlib.util

    sys.path.insert(0, str(ROOT / "scripts" / "allostery"))
    from _common import load_decisive
    spec = importlib.util.spec_from_file_location(
        "f2", ROOT / "scripts" / "allostery" / "03_pocket_specificity.py")
    f2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f2)

    ops = {  # operador consenso (hebra +) desde TF-RISet, mejor sitio S/C
        "TrpR": "TCGAAATTGTACTAGTTT",
        "CytR": "AAACGAAACATATTTTTGTGAGCAAT",
        "BetI": "TTGATTGGACGTTCAATATA",
        "PurR": "TCGCAAACGTTTGCTT",
        "GalR": "TTGTAATATTTTCATA",
        "MurR": "GGGATATATGATTCC",
    }
    decoy = "C[C@H]([NH3+])C(=O)[O-]"  # L-alanina
    cases = {"smoke": ["TrpR"],
             "pilot": ["TrpR", "CytR", "BetI"],
             "six": ["TrpR", "PurR", "GalR", "CytR", "BetI", "MurR"],
             "repaired": ["TrpR", "PurR", "GalR", "CytR", "BetI", "MurR"],
             "repaired_smoke": ["TrpR"]}[mode]
    use_msa = mode.startswith("repaired")

    jobs = []
    for d in load_decisive():
        if d["tf"] not in cases or not d["seq"]:
            continue
        chebi, smi = f2.resolve_ligand(d["tf"], d["effector"])
        if not smi:
            continue
        op = ops.get(d["tf"])
        if not op:
            continue
        # Diseño reparado: ambas condiciones son complejos de 4 cadenas del mismo
        # tamaño; solo cambia la identidad del ADN (cognato vs barajado). Así la
        # comparación no cruza composiciones distintas de complejo.
        states = ([("cognate", (op, revcomp(op))),
                   ("scrambled", (scramble(op, d["tf"]), revcomp(scramble(op, d["tf"]))))]
                  if use_msa else
                  [("dna", (op, revcomp(op))), ("free", None)])
        for kind, lig in (("true", smi), ("decoy", decoy)):
            for state, dna_arg in states:
                jobs.append({
                    "name": f"{d['tf']}|{kind}|{state}",
                    "yaml": yaml_protein_ligand(d["seq"], lig, dna_arg, use_msa),
                    "use_msa": use_msa,
                    "meta": {"tf": d["tf"], "mode": d["mode"], "kind": kind,
                             "state": state, "chebi": chebi},
                })
    return jobs


@app.function(image=modal.Image.debian_slim(python_version="3.11"), timeout=7200)
def orchestrate(jobs: list[dict]) -> int:
    """Mantiene vivo el run en modo detached: un solo contenedor orquesta el resto."""
    from concurrent.futures import ThreadPoolExecutor
    done = 0
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(run_boltz.remote, j) for j in jobs]
        for f in futs:
            try:
                f.result()
                done += 1
            except Exception as e:
                print("job falló:", e)
    return done


@app.local_entrypoint()
def main(mode: str = "smoke"):
    if mode == "fetch":
        data = read_all.remote()
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "boltz_results.json").write_text(json.dumps(data, indent=1, ensure_ascii=False))
        print(f"recogidos {len(data)} resultados -> data/allostery/boltz_results.json")
        for r in data:
            aff = r.get("affinity", {})
            print(f"{r['name']:22s} rc={r['returncode']} "
                  f"aff={aff.get('affinity_pred_value')} "
                  f"p={aff.get('affinity_probability_binary')}")
        return
    jobs = build_jobs(mode)
    print(f"lanzando {len(jobs)} jobs ({mode})", file=sys.stderr)
    orchestrate.spawn(jobs)
    print("lanzado en un solo orquestador; recoge con --mode fetch", file=sys.stderr)
