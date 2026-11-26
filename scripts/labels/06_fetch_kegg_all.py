#!/usr/bin/env python3
"""Descarga gene→KO y gene→pathway de KEGG para el catálogo local."""
from __future__ import annotations
import argparse, csv, json, re, time, urllib.request
from collections import defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]


def norm(name: str) -> str:
    s = (name or "").lower()
    s = re.sub(r"\b(candidatus|subsp\.?|str\.?|strain)\b", " ", s)
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r"[^a-z0-9 ]+", " ", s)
    return " ".join(s.split()[:2])


def parse_list(path: Path):
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2: continue
        m = re.match(r"([A-Za-z0-9]+);\s*(.+)", parts[1].strip())
        if m:
            out.append((m.group(1), m.group(2).strip()))
    return out


def fetch(url: str, dest: Path):
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    req = urllib.request.Request(url, headers={"User-Agent": "dna-token-bac/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return "ok"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT/"data")
    ap.add_argument("--ext", type=Path, default=ROOT/"data"/"labels"/"external"/"kegg")
    ap.add_argument("--delay", type=float, default=0.25)
    a = ap.parse_args()
    kegg = parse_list(a.ext/"list_genome.tsv")
    by_species = defaultdict(list)
    for code, name in kegg:
        by_species[norm(name)].append((code, name))
    reg = a.data_root/"regulation"
    rows = []
    for d in sorted(reg.iterdir()):
        if not d.is_dir() or not (d/"genes.tsv").exists(): continue
        g = d.name
        t = a.data_root/"taxonomy"/"full"/f"{g}.json"
        species = ""
        if t.exists():
            try: species = json.loads(t.read_text()).get("species") or ""
            except Exception: species = ""
        key = norm(species)
        cands = by_species.get(key, [])
        # prefer K-12 / 168 when applicable, else first
        if key == "escherichia coli":
            cands = [x for x in cands if "K-12" in x[1] or x[0] == "eco"] or cands
        if key == "bacillus subtilis":
            cands = [x for x in cands if x[0] == "bsu"] or cands
        code, name = (cands[0] if cands else ("", ""))
        rows.append({"genome": g, "species": species, "kegg_code": code,
                     "kegg_name": name, "n_candidates": len(cands)})
    matched = [r for r in rows if r["kegg_code"]]
    print(f"genomas={len(rows)} matched={len(matched)} codes={len({r['kegg_code'] for r in matched})}")
    errors = []
    for code in sorted({r["kegg_code"] for r in matched}):
        for kind, endpoint in (("ko", "link/ko"), ("pathway", "link/pathway")):
            dest = a.ext/code/f"{kind}.tsv"
            try:
                fetch(f"https://rest.kegg.jp/{endpoint}/{code}", dest)
            except Exception as exc:
                errors.append({"code": code, "kind": kind, "error": str(exc)})
            time.sleep(a.delay)
    with (a.ext/"kegg_org_map.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    (a.ext/"kegg_fetch_summary.json").write_text(json.dumps({"genomes": len(rows), "matched": len(matched), "codes": sorted({r['kegg_code'] for r in matched}), "errors": errors}, indent=1), encoding="utf-8")
    print(json.dumps({"genomes": len(rows), "matched": len(matched), "errors": len(errors)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
