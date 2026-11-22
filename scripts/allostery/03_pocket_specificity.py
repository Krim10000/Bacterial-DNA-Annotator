#!/usr/bin/env python3
"""F2-reducido — ¿el bolsillo distingue al efector correcto?

Sin co-folding (no hay Boltz/AF3 local), se usa docking ciego con AutoDock Vina
sobre estructuras AlphaFold DB. El estadístico clave es **específico por ligando**:

    para cada (TF, efector verdadero): puntuar el MISMO ligando contra el bolsillo
    de su propia proteína y contra N proteínas sorteadas.

Si el bolsillo codifica su efector, el score propio debe ser mejor (más negativo)
que el de las proteínas sorteadas. Así se evita el sesgo de tamaño/hidrofobicidad
del ligando, que invalida comparar ligandos distintos entre sí.

Métricas:
  - rank del propio entre N+1 proteínas
  - fracción top-1
  - Δ = score(propio) − media(sorteado)

Salida: data/allostery/pocket_specificity.{json,md}
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import ROOT, WORK, load_decisive  # noqa: E402

OUT = ROOT / "data" / "allostery"
AF = WORK / "af"
DOCK = WORK / "dock"
BIN = ROOT / ".venv" / "bin"


def resolve_ligand(tf, eff_raw):
    """nombre de efector -> SMILES de ChEBI (con alias y fuzzy mínimo)."""
    import difflib
    import re
    import unicodedata
    from pathlib import Path as P

    rxn = P(ROOT) / "data" / "reactions" / "tsv"

    def norm(s):
        s = re.sub(r"<[^>]+>", "", s).replace("&alpha;", "alpha").replace("&beta;", "beta")
        s = unicodedata.normalize("NFKD", s)
        s = "".join(c for c in s if not unicodedata.combining(c))
        s = s.lower().replace("_", " ")
        s = re.sub(r"[^a-z0-9+ ]", " ", s)
        return re.sub(r"\s+", " ", s).strip()

    name2chebi = {}
    for line in open(rxn / "chebiId_name.tsv", errors="ignore"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            name2chebi.setdefault(norm(f[1]), f[0])
    chebi2smi = {}
    for line in open(rxn / "rhea-chebi-smiles.tsv", errors="ignore"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            chebi2smi[f[0]] = f[1]
    n = norm(re.sub(r"^[A-Za-z0-9]+\s*-\s*", "", re.sub(r"<[^>]+>", "", eff_raw)))
    n = re.sub(r"\s+(transcriptional\s+)?(activator|repressor|dual\s+regulator)$", "", n)
    c = name2chebi.get(n)
    if not c:
        cand = difflib.get_close_matches(n, name2chebi.keys(), n=1, cutoff=0.90)
        c = name2chebi[cand[0]] if cand else None
    return (c, chebi2smi.get(c)) if c else (None, None)


def prep_receptor(pdb: Path, name: str):
    out = DOCK / f"{name}.pdbqt"
    if out.exists() and out.stat().st_size > 1000:
        return out
    r = subprocess.run([str(BIN / "mk_prepare_receptor.py"), "--read_pdb", str(pdb),
                        "-o", str(DOCK / name), "-p", "--box_enveloping", str(pdb),
                        "--padding", "4"], capture_output=True, text=True)
    return out if out.exists() else None


def prep_ligand(smiles: str, tag: str):
    """SMILES -> SDF 3D -> PDBQT (cacheado por hash)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    h = hashlib.sha1(smiles.encode()).hexdigest()[:10]
    pdbqt = DOCK / f"lig_{h}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 200:
        return pdbqt
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=0xF00D) != 0:
        if AllChem.EmbedMolecule(m, randomSeed=0xF00D, useRandomCoords=True) != 0:
            return None
    AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    sdf = DOCK / f"lig_{h}.sdf"
    w = Chem.SDWriter(str(sdf)); w.write(m); w.close()
    subprocess.run([str(BIN / "mk_prepare_ligand.py"), "-i", str(sdf),
                    "-o", str(pdbqt)], capture_output=True)
    return pdbqt if pdbqt.exists() else None


def box_from_pdb(pdb: Path, size: float = 26.0):
    """Caja cúbica centrada en el centroide de la proteína (rápido y suficiente:
    un cubo de 26 Å cubre el bolsillo para monómeros de 200-400 aa)."""
    import numpy as _np
    xs = []
    for line in open(pdb):
        if line.startswith("ATOM"):
            xs.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    c = tuple(_np.mean(xs, axis=0))
    return c, (size, size, size)


def dock(job):
    """Job: (prot_pdbqt, prot_pdb, lig_pdbqt, exhaustiveness, box, cpu)."""
    prot_pdbqt, prot_pdb, lig_pdbqt, exh, box, cpu = job
    from vina import Vina
    v = Vina(sf_name="vina", cpu=cpu, seed=42, verbosity=0)
    v.set_receptor(str(prot_pdbqt))
    v.set_ligand_from_file(str(lig_pdbqt))
    c, s = box_from_pdb(Path(prot_pdb), box)
    v.compute_vina_maps(center=c, box_size=s)
    v.dock(exhaustiveness=exh, n_poses=1)
    return float(v.score()[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-ligand", type=int, default=4)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--box", type=float, default=26.0)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--serial", action="store_true")
    ap.add_argument("--job-timeout", type=int, default=150)
    args = ap.parse_args()

    DOCK.mkdir(parents=True, exist_ok=True)
    data = [d for d in load_decisive() if d["seq"]]
    random.seed(args.seed)

    # ligandos resolubles; se excluyen los muy flexibles (CoA, etc.):
    # el docking ciego con >12 enlaces rotables no converge y no es interpretable.
    from rdkit import Chem as _Chem, RDLogger as _RDLogger
    from rdkit.Chem import rdMolDescriptors as _rdMD
    _RDLogger.DisableLog("rdApp.*")
    cases, skipped_flex = [], []
    for d in data:
        chebi, smi = resolve_ligand(d["tf"], d["effector"])
        if not smi:
            continue
        mol = _Chem.MolFromSmiles(smi)
        nrot = _rdMD.CalcNumRotatableBonds(mol) if mol else 0
        if nrot > 12:
            skipped_flex.append({"tf": d["tf"], "n_rot": nrot})
            continue
        cases.append({**d, "chebi": chebi, "smiles": smi, "n_rot": nrot})
    if args.limit:
        cases = cases[:args.limit]

    # receptores
    protos = {d["tf"]: d for d in data}
    receptors = {}
    for tf in protos:
        rec = prep_receptor(AF / f"{tf}.pdb", tf)
        if rec:
            receptors[tf] = rec
    print(f"receptores {len(receptors)} · casos con ligando {len(cases)}", file=sys.stderr)

    cache = WORK / "dock_results.json"
    results = json.loads(cache.read_text()) if cache.exists() else {}

    import multiprocessing as _mp
    import os as _os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # OpenMP dentro de Vina ignora 'cpu=' y usa todos los núcleos; con varios
    # workers eso sobresuscribe la máquina (carga 39 en 8 núcleos) y no avanza.
    _os.environ["OMP_NUM_THREADS"] = "1"
    _os.environ["MKL_NUM_THREADS"] = "1"
    _ospawn = _mp.get_context("spawn")

    jobs = []
    for c in cases:
        lig = prep_ligand(c["smiles"], c["tf"])
        if not lig or c["tf"] not in receptors:
            continue
        others = [t for t in receptors if t != c["tf"]]
        picks = random.sample(others, min(args.per_ligand, len(others)))
        for prot in [c["tf"]] + picks:
            k = f"{c['tf']}|{prot}|{hashlib.sha1(c['smiles'].encode()).hexdigest()[:10]}"
            if k not in results:
                jobs.append((k, receptors[prot], str(AF / f"{prot}.pdb"), str(lig),
                             args.exhaustiveness, args.box, args.cpu))
    print(f"dockings pendientes: {len(jobs)}", file=sys.stderr)

    if args.serial:
        # Serial: el pool con vina/OpenMP se cuelga en esta máquina (BrokenProcessPool).
        import signal as _sig

        class _JobTimeout(Exception):
            pass

        def _alarm(_sig, _frm):
            raise _JobTimeout()

        _sig.signal(_sig.SIGALRM, _alarm)
        for i, j in enumerate(jobs, 1):
            k = j[0]
            _sig.setitimer(_sig.ITIMER_REAL, args.job_timeout)
            try:
                results[k] = dock(j[1:])
            except Exception:
                results[k] = None
            finally:
                _sig.setitimer(_sig.ITIMER_REAL, 0)
            if i % 1 == 0:
                cache.write_text(json.dumps(results))
                print(f"  {i}/{len(jobs)}", file=sys.stderr)
        cache.write_text(json.dumps(results))
        return _report(results, cases, args)

    with ProcessPoolExecutor(max_workers=args.workers, mp_context=_ospawn) as ex:
        futs = {ex.submit(dock, j[1:]): j[0] for j in jobs}
        done = 0
        for f in as_completed(futs):
            k = futs[f]
            try:
                results[k] = f.result()
            except Exception:
                results[k] = None
            done += 1
            if done % 5 == 0:
                cache.write_text(json.dumps(results))
                print(f"  {done}/{len(jobs)}", file=sys.stderr)
    cache.write_text(json.dumps(results))

    return _report(results, cases, args)


def _report(results, cases, args):
    import numpy as np
    # métricas
    rows = []
    for c in cases:
        h = hashlib.sha1(c["smiles"].encode()).hexdigest()[:10]
        own = results.get(f"{c['tf']}|{c['tf']}|{h}")
        others = [v for k, v in results.items()
                  if k.startswith(f"{c['tf']}|") and not k.startswith(f"{c['tf']}|{c['tf']}|")
                  and v is not None]
        if own is None or len(others) != args.per_ligand:
            continue
        better = sum(1 for o in others if own < o)
        rows.append({"tf": c["tf"], "mode": c["mode"], "chebi": c["chebi"],
                     "own": own, "others": others, "n_better": better,
                     "rank": 1 + sum(1 for o in others if o < own),
                     "rank_norm": (1 + sum(1 for o in others if o < own)) / (1 + len(others)),
                     "delta": own - float(np.mean(others))})
    if not rows:
        print(json.dumps({"error": "sin resultados"}, ensure_ascii=False))
        return
    ranks = [r["rank"] for r in rows]
    n_prot = max(r["rank"] for r in rows)
    summary = {
        "n_casos": len(rows),
        "n_proteínas_control": max(len(r["others"]) for r in rows),
        "top1": sum(1 for r in rows if r["rank"] == 1) / len(rows),
        "top2": sum(1 for r in rows if r["rank"] <= 2) / len(rows),
        "rank_medio": float(np.mean(ranks)),
        "rank_norm_medio": float(np.mean([r["rank_norm"] for r in rows])),
        "azar_top1": 1 / (1 + max(len(r["others"]) for r in rows)),
        "delta_medio": float(np.mean([r["delta"] for r in rows])),
        "delta_negativos": sum(1 for r in rows if r["delta"] < 0) / len(rows),
        "exhaustiveness": args.exhaustiveness,
    }
    (OUT / "pocket_specificity.json").write_text(json.dumps(
        {"summary": summary, "rows": rows}, indent=1, ensure_ascii=False))
    L = ["# F2-reducido — especificidad de bolsillo (docking ciego)", "",
         f"n={summary['n_casos']} · controles por ligando={summary['n_proteínas_control']} · "
         f"exhaustiveness={summary['exhaustiveness']}", "",
         "| métrica | valor | azar |", "|---|---:|---:|",
         f"| top-1 | {summary['top1']:.3f} | {summary['azar_top1']:.3f} |",
         f"| top-2 | {summary['top2']:.3f} | {2 * summary['azar_top1']:.3f} |",
         f"| rank medio | {summary['rank_medio']:.2f} | {(1 + summary['n_proteínas_control'] + 1) / 2:.2f} |",
         f"| Δ medio (propio − controles) | {summary['delta_medio']:.3f} kcal/mol | 0 |",
         f"| fracción con Δ<0 | {summary['delta_negativos']:.3f} | 0.5 |", ""]
    (OUT / "pocket_specificity.md").write_text("\n".join(L) + "\n")
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
