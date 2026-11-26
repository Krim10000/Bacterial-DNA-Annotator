#!/usr/bin/env python3
"""Perfil de etiquetas/pathways por sigma (K-12 gold; BSU descriptivo)."""
from __future__ import annotations
import argparse, csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def read(p: Path):
    if not p.exists(): return []
    with p.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def label_map(path: Path, col="labels"):
    return {r["operon_id"]: [x for x in (r.get(col) or "").replace(";", ",").split(",") if x]
            for r in read(path)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT/"data")
    ap.add_argument("--labels", type=Path, default=ROOT/"reports/labels_pilot")
    a = ap.parse_args()
    # K12 union: lexicon + GO + KEGG
    proc = defaultdict(set)
    for p in (a.labels/"ecoli_k12_operon_labels.tsv", a.labels/"ecoli_k12_go_operon_labels.tsv"):
        for oid, labs in label_map(p).items(): proc[oid].update(labs)
    path = label_map(a.labels/"ecoli_k12_kegg_pathway_operons.tsv", "pathways")
    names = {}
    pn = a.data_root/"labels"/"external"/"kegg"/"eco"/"pathway_names.tsv"
    if pn.exists():
        for line in pn.read_text(encoding="utf-8", errors="replace").splitlines():
            if "\t" in line:
                pid, name = line.split("\t", 1)
                names[pid] = name.split(" - ")[0].strip()
    ops = {r["label"]: r["operon_id"] for r in read(a.data_root/"regulation/ecoli_k12/operons.tsv")}
    net = read(a.data_root/"regulondb/NetworkSigmaTU.tsv")
    targets = defaultdict(set)
    for r in net:
        sig = (r.get("1)sigmaName") or "").strip()
        tu = re.split(r"[\[\s]", (r.get("2)regulatedTUName") or "").strip())[0]
        if sig and tu in ops: targets[sig].add(ops[tu])
    out = {}
    for sig, oids in sorted(targets.items()):
        pc = Counter(l for oid in oids for l in proc.get(oid, ()))
        kc = Counter(names.get(x, x) for oid in oids for x in path.get(oid, ()))
        out[sig] = {"n_targets": len(oids),
                    "process_labels": pc.most_common(8),
                    "pathways": kc.most_common(8)}
    (a.labels/"k12_sigma_profiles.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    lines = ["# Perfil de etiquetas por sigma (K-12)", "",
             "Fuente de dianas: `NetworkSigmaTU.tsv`. Etiquetas: léxico + GO + KEGG.", ""]
    for sig, d in out.items():
        lines += [f"## {sig}", f"Dianas mapeadas: {d['n_targets']}", "",
                  "| Etiqueta de proceso | n |", "|---|---:|"]
        lines += [f"| {k} | {v} |" for k, v in d["process_labels"]]
        lines += ["", "| Pathway KEGG | n |", "|---|---:|"]
        lines += [f"| {k} | {v} |" for k, v in d["pathways"]]
        lines.append("")
    (a.labels/"k12_sigma_profiles.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(out, indent=1, ensure_ascii=False)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
