"""Evidencia estructural multi-fuente para una proteína.

Fuentes:
- Foldseek API (PDB100) o fallback local AFDB.
- AlphaFill (ligandos/cofactores transplantados).
- Sitios activos locales (ACTIVE_SITES.tsv / mapa estructural).

No declara cierre. Solo agrega evidencia y procedencia.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "metabolism" / "targets"
CACHE = DATA / "structures" / "evidence_cache"
AFDB = DATA / "structures" / "afdb"
FOLDSEEK_API = "https://search.foldseek.com/api"
ALPHAFILL_API = "https://alphafill.eu/v1/aff"


def _cache(key: str, value=None):
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{key}.json"
    if value is None:
        return json.loads(p.read_text()) if p.exists() else None
    p.write_text(json.dumps(value, indent=2)[:5_000_000])
    return value


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:20]


def foldseek_search(pdb: Path | str, database: str = "pdb100", max_hits: int = 10,
                    timeout: int = 180) -> dict:
    """Devuelve {source, hits, error}. Usa API y cae a AFDB local."""
    pdb = Path(pdb)
    key = "fs_" + _hash(pdb) + "_" + database
    cached = _cache(key)
    if cached:
        return cached
    try:
        with open(pdb, "rb") as fh:
            r = requests.post(
                f"{FOLDSEEK_API}/ticket",
                files={"q": fh},
                data={"database[]": database, "mode": "3diaa"},
                timeout=60,
            )
        r.raise_for_status()
        jid = r.json().get("id")
        status = "PENDING"
        deadline = time.time() + timeout
        while time.time() < deadline:
            s = requests.get(f"{FOLDSEEK_API}/ticket/{jid}", timeout=30).json()
            status = s.get("status", "")
            if status == "COMPLETE":
                break
            if status in ("ERROR", "FAILED"):
                raise RuntimeError(f"Foldseek {status}")
            time.sleep(3)
        if status != "COMPLETE":
            raise TimeoutError("Foldseek timeout")
        data = requests.get(f"{FOLDSEEK_API}/result/{jid}/0", timeout=60).json()
        alns = ((data.get("results") or [{}])[0].get("alignments") or [[]])[0] or []
        hits = []
        for a in alns[:max_hits]:
            target = a.get("target", "")
            m = re.match(r"([0-9A-Za-z]{4})", target)
            hits.append({
                "pdb_id": m.group(1).lower() if m else "",
                "target": target,
                "seq_id": a.get("seqId"),
                "bits": a.get("score"),
                "evalue": a.get("eval"),
                "aln_len": a.get("alnLength"),
                "q_start": a.get("qStartPos"), "q_end": a.get("qEndPos"),
                "t_start": a.get("dbStartPos"), "t_end": a.get("dbEndPos"),
            })
        out = {"source": "foldseek_api_pdb100", "hits": hits, "error": ""}
    except Exception as exc:
        out = _foldseek_local(pdb, max_hits)
        out["error"] = f"API: {type(exc).__name__}: {exc}"
    return _cache(key, out)


def _foldseek_local(pdb: Path, max_hits: int) -> dict:
    """Fallback: AFDB local (8 homólogos), no PDB100."""
    if not AFDB.exists():
        return {"source": "local_afdb", "hits": [], "error": "sin AFDB local"}
    try:
        with __import__("tempfile").TemporaryDirectory(prefix="fs_ev_") as td:
            out = Path(td) / "hits.m8"
            subprocess.run([
                "foldseek", "easy-search", str(pdb), str(AFDB), str(out), td,
                "--threads", "2", "--max-seqs", str(max_hits), "-v", "0",
                "--format-output", "target,fident,alnlen,evalue,bits",
            ], capture_output=True, text=True, timeout=120)
            hits = []
            if out.exists():
                for line in out.read_text().splitlines()[:max_hits]:
                    p = line.split("\t")
                    if len(p) >= 5:
                        hits.append({"pdb_id": Path(p[0]).stem.split("_")[0].lower(),
                                     "target": p[0], "seq_id": p[1], "bits": p[4],
                                     "evalue": p[3], "aln_len": p[2]})
            return {"source": "local_afdb", "hits": hits, "error": ""}
    except Exception as exc:
        return {"source": "local_afdb", "hits": [], "error": f"{type(exc).__name__}: {exc}"}


def alphafill(uniprot: str) -> dict:
    """Ligandos/cofactores transplantados por AlphaFill."""
    if not uniprot:
        return {"source": "alphafill", "hits": [], "error": "sin UniProt"}
    key = "af2_" + uniprot
    cached = _cache(key)
    if cached:
        return cached
    try:
        d = requests.get(f"{ALPHAFILL_API}/{uniprot}/json", timeout=12).json()
        hits = []
        for h in (d.get("hits") or []):
            for t in (h.get("transplants") or []):
                hits.append({
                    "pdb_id": h.get("pdb_id", ""),
                    "ligand": t.get("analogue_id", ""),
                    "rmsd": h.get("global_rmsd"),
                    "clash": (t.get("clash") or {}).get("clash_count"),
                })
        hits.sort(key=lambda x: (x["clash"] is None, x["clash"] or 0, x["rmsd"] or 99))
        out = {"source": "alphafill", "hits": hits[:20], "error": ""}
    except Exception as exc:
        out = {"source": "alphafill", "hits": [], "error": f"{type(exc).__name__}: {exc}"}
    return _cache(key, out)


def local_sites(genome: str, uniprot: str = "") -> dict:
    p = DATA / "ACTIVE_SITE_STRUCTURAL_MAP.tsv"
    if not p.exists():
        return {"source": "local_sites", "hits": [], "error": "sin mapa"}
    rows = []
    for line in p.read_text().splitlines()[1:]:
        c = line.split("\t")
        if len(c) >= 7 and c[0] == genome:
            rows.append({"acc": c[1], "feature": c[2], "start": c[4], "end": c[5],
                         "ref_pos": c[6], "qry_pos": c[8] if len(c) > 8 else "",
                         "conserved": c[10] if len(c) > 10 else ""})
    return {"source": "local_sites", "hits": rows, "error": ""}


def summarize(genome: str, pdb: str, uniprot: str = "") -> dict:
    fs = foldseek_search(pdb)
    af = alphafill(uniprot)
    sites = local_sites(genome, uniprot)
    top = (fs.get("hits") or [{}])[0]
    lig = [x for x in (af.get("hits") or []) if x.get("clash") is not None and x["clash"] <= 20][:5]
    active = [x for x in (sites.get("hits") or []) if x.get("feature") == "Active site"]
    sources = int(bool(fs.get("hits"))) + int(bool(af.get("hits"))) + int(bool(sites.get("hits")))
    return {
        "genome": genome, "pdb": pdb, "uniprot": uniprot,
        "foldseek_source": fs.get("source"), "foldseek_top": top.get("pdb_id", ""),
        "foldseek_bits": top.get("bits", ""), "foldseek_seq_id": top.get("seq_id", ""),
        "alphafill_ligands": ";".join(f"{x['ligand']}({x['pdb_id']},clash={x['clash']})" for x in lig),
        "active_site_rows": len(active),
        "active_site_conserved": sum(1 for x in active if str(x.get("conserved")).lower() == "true"),
        "n_sources": sources,
        "error": "; ".join(x for x in (fs.get("error"), af.get("error"), sites.get("error")) if x),
        "verdict": "multi_source_support" if sources >= 2 and active else "structure_only",
        "closure": "no",
        "limitation": "Evidencia computacional; no confirma actividad ni cierra la brecha.",
    }


def evidence_for_rows(rows: list[dict]) -> list[dict]:
    """Cruza Foldseek + AlphaFill + sitios activos para filas de cierre."""
    out = []
    for r in rows:
        pdb = Path(str(r.get("pdb") or ""))
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        if not pdb.exists():
            continue
        uniprot = str(r.get("uniprot") or r.get("best_acc") or r.get("homolog_acc") or "")
        s = summarize(str(r.get("genome") or ""), str(pdb), uniprot)
        s.update({"gap": r.get("gap", ""), "query": r.get("query", ""), "pdb": str(pdb)})
        out.append(s)
    return out
