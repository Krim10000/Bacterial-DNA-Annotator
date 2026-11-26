#!/usr/bin/env python3
"""Perfil descriptivo de etiquetas por sigma en B. subtilis desde TSS gold."""
from __future__ import annotations
import argparse, csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def read(p: Path):
    if not p.exists(): return []
    with p.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def labs(path, col="labels"):
    return {r["operon_id"]: [x for x in (r.get(col) or "").replace(";", ",").split(",") if x] for r in read(path)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT/"data")
    ap.add_argument("--labels", type=Path, default=ROOT/"reports/labels_pilot")
    a = ap.parse_args()
    proc = defaultdict(set)
    for p in (a.labels/"subtilis_operon_labels.tsv", a.labels/"subtilis_go_operon_labels.tsv"):
        for oid, xs in labs(p).items(): proc[oid].update(xs)
    path = labs(a.labels/"subtilis_kegg_pathway_operons.tsv", "pathways")
    names = {}
    pn = a.data_root/"labels/external/kegg/bsu/pathway_names.tsv"
    for line in pn.read_text().splitlines():
        if "\t" in line:
            pid, name = line.split("\t",1); names[pid] = name.split(" - ")[0].strip()
    ops = read(a.data_root.parent/"results/subtilis/subtilis.operons_auto.tsv")
    gold = read(a.data_root/"sigma/gold/subtilis/tss_sigma.tsv")
    buckets = defaultdict(set)
    # map TSS to operon start on same contig/strand, within 1500 bp upstream
    by_key = defaultdict(list)
    for op in ops:
        by_key[(op.get("contig"), str(op.get("strand")))].append(op)
    for t in gold:
        try: pos = int(t["pos"])
        except Exception: continue
        key = (t.get("contig"), t.get("strand"))
        cand = []
        for op in by_key.get(key, []):
            try: s = int(op["start"])
            except Exception: continue
            if 0 <= s - pos <= 1500 or s <= pos <= int(op["end"]):
                cand.append((abs(s-pos), op["operon_id"]))
        if cand:
            buckets[t.get("sigma","")].add(sorted(cand)[0][1])
    out = {}
    for sig, oids in sorted(buckets.items()):
        out[sig] = {"n_tss_mapped": len(oids),
                    "process_labels": Counter(l for oid in oids for l in proc.get(oid,())).most_common(8),
                    "pathways": Counter(names.get(x,x) for oid in oids for x in path.get(oid,())).most_common(8)}
    (a.labels/"bsu_sigma_profiles.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    lines=["# Perfil de etiquetas por sigma (B. subtilis)","",
           "TSS gold mapeados a operones existentes; descriptivo, no validación independiente.",""]
    for sig,d in out.items():
        lines += [f"## {sig}", f"TSS mapeados: {d['n_tss_mapped']}", "", "| Etiqueta | n |", "|---|---:|"]
        lines += [f"| {k} | {v} |" for k,v in d["process_labels"]]
        lines += ["", "| Pathway KEGG | n |", "|---|---:|"]
        lines += [f"| {k} | {v} |" for k,v in d["pathways"]]
        lines.append("")
    (a.labels/"bsu_sigma_profiles.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(out, indent=1, ensure_ascii=False)[:4000])
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
