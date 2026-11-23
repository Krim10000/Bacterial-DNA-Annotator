"""Dominios Pfam para el camino C (sin Swiss-Prot).

No asigna función ni complejo: aporta una hipótesis de familia/dominio con
e-value. Usa hmmscan contra un panel Pfam curado, no Pfam completo.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def available(data_root: Path) -> tuple[bool, str]:
    if shutil.which("hmmscan") is None:
        return False, "hmmscan no está en PATH"
    panel = Path(data_root) / "db" / "pfam" / "panel_focus.hmm"
    if not panel.exists():
        return False, f"falta {panel}"
    return True, ""


def _parse_domtbl(path: Path) -> dict[str, dict]:
    best: dict[str, dict] = {}
    if not path.exists():
        return best
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        pid = f[3]
        try:
            evalue = float(f[12])
            full_e = float(f[6])
            score = float(f[13])
        except ValueError:
            continue
        rec = {
            "domain": f[0], "accession": f[1], "evalue": evalue,
            "full_evalue": full_e, "score": score,
            "ali_from": f[17], "ali_to": f[18],
            "description": " ".join(f[22:]) if len(f) > 22 else "",
        }
        if pid not in best or evalue < best[pid]["evalue"]:
            best[pid] = rec
    return best


def run(rows: list[dict], proteins: list[dict], outdir: Path, genome: str,
        data_root: Path, cpu: int = 4) -> dict[str, Any]:
    ok, why = available(data_root)
    if not ok:
        return {"ok": False, "reason": why, "rows": [], "n": 0}
    c_ids = {r.get("protein_id") for r in rows if r.get("path") == "C"}
    seqs = {p.get("id"): p.get("aa", "") for p in proteins
            if p.get("id") in c_ids and p.get("aa")}
    if not seqs:
        return {"ok": True, "reason": "", "rows": [], "n": 0}
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="domains_") as td:
        tmp = Path(td)
        fa = tmp / "c.faa"
        with fa.open("w", encoding="utf-8") as fh:
            for pid, seq in seqs.items():
                fh.write(f">{pid}\n{seq}\n")
        domtbl = tmp / "domains.domtbl"
        panel = Path(data_root) / "db" / "pfam" / "panel_focus.hmm"
        cmd = ["hmmscan", "--domtblout", str(domtbl), "-o", "/dev/null",
               "--cpu", str(max(1, cpu)), str(panel), str(fa)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if proc.returncode != 0:
            return {"ok": False, "reason": (proc.stderr or proc.stdout)[-400:],
                    "rows": [], "n": 0}
        hits = _parse_domtbl(domtbl)
    out_rows = []
    for r in rows:
        pid = r.get("protein_id")
        if r.get("path") != "C" or pid not in hits:
            continue
        h = hits[pid]
        r["domain"] = h["domain"]
        r["domain_accession"] = h["accession"]
        r["domain_evalue"] = h["evalue"]
        r["domain_score"] = h["score"]
        r["domain_evidence"] = "hmmscan_panel_focus"
        r["domain_note"] = ("Hipótesis de dominio/familia; no confirma función "
                            "ni pertenencia a complejo.")
        out_rows.append({"genome": genome, "protein_id": pid, **h})
    path = outdir / f"{genome}.domains.tsv"
    cols = list(out_rows[0]) if out_rows else ["genome", "protein_id"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(out_rows)
    return {"ok": True, "reason": "", "rows": out_rows, "n": len(out_rows),
            "tsv": str(path)}
