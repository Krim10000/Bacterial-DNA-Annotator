#!/usr/bin/env python3
"""Mapea sitios activos UniProt sobre PDB ESMFold usando Foldseek.

No declara cierre: mapea residuos, comprueba conservación y mide dispersión.
"""
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "data" / "metabolism" / "targets"
RESULTS = TARGETS / "NEUROSNAP_4_RESULTS.tsv"
SITES = TARGETS / "ACTIVE_SITES.tsv"
AFDB = TARGETS / "structures" / "afdb"


def foldseek_alignment(query: Path, target: Path):
    with tempfile.TemporaryDirectory(prefix="fsmap_") as td:
        out = Path(td) / "aln.m8"
        cmd = [
            "foldseek", "easy-search", str(query), str(target), str(out), str(Path(td) / "tmp"),
            "--threads", "2", "--max-seqs", "1", "-v", "0",
            "--format-output", "qstart,tstart,qend,tend,qaln,taln",
        ]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode != 0 or not out.exists():
            return None
        lines = [x for x in out.read_text().splitlines() if x.strip()]
        if not lines:
            return None
        p = lines[0].split("\t")
        if len(p) < 6:
            return None
        return p


def ca_coords(path: Path) -> list[tuple[float, float, float]]:
    pts = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pts.append((0.0, 0.0, 0.0))
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--out", type=Path, default=TARGETS / "ACTIVE_SITE_STRUCTURAL_MAP")
    args = ap.parse_args()
    results = list(csv.DictReader(open(args.results), delimiter="\t"))
    sites = list(csv.DictReader(open(SITES), delimiter="\t"))
    rows, summary = [], []
    for r in results:
        acc = (r.get("esperado") or "").strip()
        pdb = Path((r.get("PDB") or "").strip())
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        plddt = float(r["pLDDT"]) if r.get("pLDDT") else 0.0
        if not acc or not pdb.exists():
            summary.append(dict(genoma=r["genoma"], acc=acc, plddt=plddt, features=0, mapped=0,
                                conserved=0, active_total=0, active_cons=0, max_ca="",
                                verdict="no_homolog" if not acc else "missing_pdb", closure="no"))
            continue
        ref = AFDB / f"{acc}.pdb"
        if not ref.exists():
            summary.append(dict(genoma=r["genoma"], acc=acc, plddt=plddt, features=0, mapped=0,
                                conserved=0, active_total=0, active_cons=0, max_ca="",
                                verdict="missing_ref", closure="no"))
            continue
        aln = foldseek_alignment(pdb, ref)
        if not aln:
            summary.append(dict(genoma=r["genoma"], acc=acc, plddt=plddt, features=0, mapped=0,
                                conserved=0, active_total=0, active_cons=0, max_ca="",
                                verdict="no_alignment", closure="no"))
            continue
        qstart, tstart, qend, tend, qaln, taln = aln[:6]
        qstart, tstart = int(qstart), int(tstart)
        coords = ca_coords(pdb)
        # map target (ref) 1-based -> query 0-based / aligned aa
        mapping: dict[int, tuple[int, str, str]] = {}
        qi = qstart - 1
        ti = tstart - 1
        for qa, ta in zip(qaln.upper(), taln.upper()):
            if ta != "-":
                if qa != "-":
                    mapping[ti + 1] = (qi, qa, ta)
                    qi += 1
                else:
                    mapping[ti + 1] = (-1, "-", ta)
                ti += 1
            elif qa != "-":
                qi += 1
        feats = []
        seen = set()
        for s in sites:
            if s.get("acc") != acc:
                continue
            try:
                a, b = int(s["start"]), int(s["end"])
            except Exception:
                continue
            key = (s.get("feature"), a, b)
            if key in seen:
                continue
            seen.add(key); feats.append((s, a, b))
        mapped = conserved = active_total = active_cons = 0
        pts = []
        for s, a, b in feats:
            for pos in range(a, b + 1):
                m = mapping.get(pos)
                qpos = qaa = refaa = ""
                if m:
                    qi2, qaa, refaa = m
                    qpos = qi2 + 1 if qi2 >= 0 else ""
                    if qi2 >= 0:
                        mapped += 1
                        if qaa == refaa:
                            conserved += 1
                            if s.get("feature") == "Active site":
                                active_cons += 1
                        if qi2 < len(coords):
                            pts.append(coords[qi2])
                if s.get("feature") == "Active site":
                    active_total += 1
                rows.append(dict(genoma=r["genoma"], acc=acc, feature=s.get("feature", ""),
                                 description=s.get("description", ""), start=a, end=b,
                                 ref_pos=pos, ref_aa=refaa, qry_pos=qpos, qry_aa=qaa,
                                 conserved=bool(qaa and qaa == refaa)))
        max_ca = ""
        if len(pts) >= 2:
            max_ca = round(max(math.dist(x, y) for i, x in enumerate(pts) for y in pts[i + 1:]), 1)
        if not feats:
            verdict = "no_active_site"
        elif mapped == 0:
            verdict = "unmapped"
        elif active_total and active_cons == active_total and plddt >= 80:
            verdict = "active_site_conserved"
        elif active_total and active_cons:
            verdict = "active_site_partial"
        elif conserved:
            verdict = "binding_partial"
        else:
            verdict = "mapped_not_conserved"
        summary.append(dict(genoma=r["genoma"], acc=acc, plddt=plddt, features=len(feats),
                            mapped=mapped, conserved=conserved, active_total=active_total,
                            active_cons=active_cons, max_ca=max_ca, verdict=verdict, closure="no"))
    cols = ["genoma", "acc", "feature", "description", "start", "end", "ref_pos", "ref_aa", "qry_pos", "qry_aa", "conserved"]
    out = Path(str(args.out) + ".tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    md = ["# Mapa estructural de sitios activos", "",
          "Foldseek: PDB ESMFold del ORF → AFDB del homólogo. No es cierre.", "",
          "| genoma | acc | pLDDT | features | mapped | conserved | active cons/total | max CA Å | veredicto | cierre |",
          "|---|---|---:|---:|---:|---:|---|---:|---|---|"]
    for s in summary:
        md.append(f"| {s['genoma']} | {s['acc']} | {s['plddt']} | {s['features']} | {s['mapped']} | {s['conserved']} | "
                  f"{s['active_cons']}/{s['active_total']} | {s['max_ca']} | {s['verdict']} | {s['closure']} |")
    md += ["", "## Criterio", "",
           "- `active_site_conserved`: todos los residuos `Active site` mapeados y conservados, pLDDT ≥80.",
           "- `active_site_partial`: algunos residuos catalíticos conservados.",
           "- `binding_partial`: solo sitios de unión, no catalíticos.",
           "- `unmapped` / `mapped_not_conserved`: no apoya función.",
           "- **`cierre = no`**: estructura y conservación no prueban actividad ni cierran la brecha.",
           "- `max CA Å` mide dispersión espacial; no es criterio de función."]
    Path(str(args.out) + ".md").write_text("\n".join(md) + "\n")
    for s in summary:
        print(s)


if __name__ == "__main__":
    main()
