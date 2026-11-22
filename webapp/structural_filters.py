"""Cinco filtros estructurales para priorizar candidatos.

F1 Foldseek (PDB100/AFDB)
F2 EnzyMM/PyJess (plantillas M-CSA)
F3 P2Rank (bolsillo)
F4 AlphaFill (ligando/cofactor)
F5 sitio activo local (UniProt/estructura)

No cierra brechas. Reduce candidatos y deja procedencia.
"""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from structure_evidence import alphafill, foldseek_search, local_sites
from structural_annotation import scan_pdbs, summarize as summarize_enzymm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "metabolism" / "targets"
CACHE = DATA / "structures" / "filter_cache"


def _hash(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:20]


def _cache(key: str, value=None):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{key}.json"
    if value is None:
        return json.loads(f.read_text()) if f.exists() else None
    f.write_text(json.dumps(value, indent=2)[:2_000_000])
    return value


def p2rank(pdb: Path | str, timeout: int = 180) -> dict:
    pdb = Path(pdb)
    key = "p2rank_" + _hash(pdb)
    cached = _cache(key)
    if cached is not None:
        return cached
    exe = shutil.which("p2rank") or str(Path.home() / ".local/bin/p2rank")
    if not Path(exe).exists():
        return {"pockets": [], "error": "p2rank no instalado"}
    try:
        with tempfile.TemporaryDirectory(prefix="p2r_") as td:
            out = Path(td) / "out"
            r = subprocess.run([exe, "predict", "-threads", "2", "-f", str(pdb.resolve()), "-o", str(out)],
                               capture_output=True, text=True, timeout=timeout)
            csvs = list(out.glob("*_predictions.csv"))
            if not csvs:
                raise RuntimeError((r.stderr or r.stdout)[-300:])
            pockets = []
            with open(csvs[0]) as fh:
                reader = csv.reader(fh)
                header = next(reader)
                header = [h.strip() for h in header]
                for row in reader:
                    d = {header[i]: row[i].strip() for i in range(min(len(header), len(row)))}
                    pockets.append(d)
            res = {"pockets": pockets, "error": "", "csv": str(csvs[0])}
    except Exception as exc:
        res = {"pockets": [], "error": f"{type(exc).__name__}: {exc}"}
    return _cache(key, res)


def _active_positions(genome: str, uniprot: str = "", sites_file: Path | None = None) -> list[int]:
    if not genome or not uniprot:
        return []
    paths = [sites_file] if sites_file else [
        DATA / "ACTIVE_SITE_STRUCTURAL_MAP.tsv",
        DATA / "gap_glycerate" / "ACTIVE_SITE_STRUCTURAL_MAP.tsv",
    ]
    out = []
    for p in paths:
        if not p or not Path(p).exists():
            continue
        for line in Path(p).read_text().splitlines()[1:]:
            c = line.split("\t")
            if (len(c) > 10 and c[0] == genome and c[1] == uniprot
                    and c[2] == "Active site" and str(c[10]).lower() == "true"):
                try:
                    out.append(int(c[8]))
                except Exception:
                    pass
    return sorted(set(out))


def _pocket_hits(pockets: list[dict], positions: list[int]) -> list[str]:
    if not positions:
        return []
    hits = []
    for p in pockets[:3]:
        ids = p.get("residue_ids", "")
        for token in ids.replace("\t", " ").split():
            # A_171
            try:
                num = int(token.split("_")[-1])
            except Exception:
                continue
            if num in positions:
                hits.append(token)
    return sorted(set(hits))


def run_filters(pdb: Path | str, genome: str = "", uniprot: str = "",
                sites_file: Path | None = None, foldseek_hits: list[dict] | None = None) -> dict:
    pdb = Path(pdb)
    key = "filters_v5_" + _hash(pdb) + "_" + genome + "_" + uniprot + "_" + (str(sites_file) if sites_file else "auto") + "_" + ("local" if foldseek_hits is not None else "api")
    cached = _cache(key)
    if cached is not None:
        return cached
    fs = ({"source": "local_foldseek_swissprot", "hits": foldseek_hits or [], "error": ""}
          if foldseek_hits is not None else foldseek_search(pdb))
    enz_key = "enzymm_" + _hash(pdb)
    enz = _cache(enz_key)
    if enz is None:
        enz = summarize_enzymm(scan_pdbs([pdb], CACHE / f"enzymm_{_hash(pdb)}", jobs=2).get("hits") or [])
        _cache(enz_key, enz)
    pk = p2rank(pdb)
    af = alphafill(uniprot) if uniprot else {"hits": [], "error": "sin UniProt"}
    act = _active_positions(genome, uniprot, sites_file=sites_file)
    ph = _pocket_hits(pk.get("pockets", []), act)
    f = {
        "f1_foldseek": bool(fs.get("hits")),
        "f2_enzymm": bool(enz),
        "f3_pocket": bool(pk.get("pockets")),
        "f4_alphafill": bool(af.get("hits")),
        "f5_active_site": bool(act),
    }
    applicable = sum(1 for k in ("f1_foldseek", "f2_enzymm", "f4_alphafill", "f5_active_site") if True)
    passed = sum(int(f[k]) for k in ("f1_foldseek", "f2_enzymm", "f3_pocket", "f4_alphafill", "f5_active_site"))
    if f["f2_enzymm"] and f["f5_active_site"]:
        priority = "high"
    elif passed >= 2:
        priority = "review"
    else:
        priority = "low"
    out = {
        "pdb": str(pdb), "genome": genome, "uniprot": uniprot,
        "f1_foldseek": f["f1_foldseek"], "foldseek_top": (fs.get("hits") or [{}])[0].get("pdb_id", ""),
        "foldseek_bits": (fs.get("hits") or [{}])[0].get("bits", ""),
        "f2_enzymm": f["f2_enzymm"], "enzymm_hits": len(enz),
        "enzymm_templates": ";".join(x.get("template_pdb_id", "") for x in enz[:5]),
        "enzymm_ecs": ";".join(sorted({str(x.get("template_ec", "")) for x in enz if x.get("template_ec")})),
        "f3_pocket": f["f3_pocket"], "p2rank_pockets": len(pk.get("pockets", [])),
        "pocket_active_hits": ";".join(ph),
        "f4_alphafill": f["f4_alphafill"],
        "alphafill_ligands": ";".join(f"{x.get('ligand')}({x.get('pdb_id')})" for x in (af.get("hits") or [])[:5]),
        "f5_active_site": f["f5_active_site"], "active_positions": ";".join(map(str, act)),
        "n_pass": passed, "priority": priority,
        "source": "Foldseek + EnzyMM/PyJess + P2Rank + AlphaFill + UniProt/estructura",
        "limitation": "Filtros computacionales; no confirman función ni cierran brechas.",
        "error": "; ".join(x for x in (fs.get("error"), pk.get("error"), af.get("error")) if x),
    }
    return _cache(key, out)
