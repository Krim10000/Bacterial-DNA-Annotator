#!/usr/bin/env python3
"""Agrega pathway labels (KEGG vía Rhea) a operones y valida K-12."""
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


def paths(s): return [x for x in (s or "").split(";") if x]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "labels_pilot")
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    cases = {
        "ecoli_k12": (a.data_root / "regulation" / "ecoli_k12" / "operons.tsv", "gene_ids", ","),
        "subtilis": (a.data_root.parent / "results" / "subtilis" / "subtilis.operons_auto.tsv", "genes", ";"),
    }
    summary = {}
    for genome, (op_path, col, sep) in cases.items():
        cw = {r["id"]: paths(r["pathways"]) for r in read(a.out / f"{genome}_catalog_crosswalk.tsv")}
        ops = read(op_path)
        rows = []
        for op in ops:
            ids = [x for x in (op.get(col) or "").split(sep) if x]
            ps = sorted({p for i in ids for p in cw.get(i, [])})
            rows.append({"operon_id": op.get("operon_id", ""), "label": op.get("label", ""),
                         "n_genes": len(ids), "n_pathways": len(ps),
                         "pathways": ";".join(ps)})
        summary[genome] = {"n_operons": len(rows),
                           "n_with_pathway": sum(1 for r in rows if r["pathways"]),
                           "n_pathways": len({p for r in rows for p in paths(r["pathways"])}),
                           "top": Counter(p for r in rows for p in paths(r["pathways"])).most_common(10)}
        with (a.out / f"{genome}_pathway_operons.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["operon_id"], delimiter="\t")
            w.writeheader(); w.writerows(rows)
    # K-12: enriquecimiento de pathways en operones diana de cada sigma.
    cw = {r["id"]: paths(r["pathways"]) for r in read(a.out / "ecoli_k12_catalog_crosswalk.tsv")}
    ops = {r["label"]: r["operon_id"] for r in read(a.data_root / "regulation" / "ecoli_k12" / "operons.tsv")}
    op_path = {r["operon_id"]: paths(r["pathways"]) for r in read(a.out / "ecoli_k12_pathway_operons.tsv")}
    net = read(a.data_root / "regulondb" / "NetworkSigmaTU.tsv")
    targets = defaultdict(set)
    for r in net:
        sig = (r.get("1)sigmaName") or "").strip()
        tu = re.split(r"[\[\s]", (r.get("2)regulatedTUName") or "").strip())[0]
        if sig and tu in ops: targets[sig].add(ops[tu])
    all_paths = Counter(p for ps in op_path.values() for p in ps)
    validation = {}
    for sig, tids in sorted(targets.items()):
        if not tids: continue
        target_counts = Counter(p for oid in tids for p in op_path.get(oid, []))
        bg = {p: all_paths[p] for p in all_paths}
        scored = []
        for p, n in target_counts.items():
            a1 = n; b1 = len(tids) - n
            c1 = bg.get(p, 0) - n; d1 = (len(op_path) - len(tids)) - c1
            if n < 2 or b1 < 0: continue
            try:
                from scipy.stats import fisher_exact
                odds, pv = fisher_exact([[a1, b1], [c1, d1]], alternative="greater")
            except Exception:
                odds, pv = None, None
            scored.append({"pathway": p, "n_target": n, "odds": round(float(odds), 3) if odds else None,
                           "pvalue": round(float(pv), 6) if pv else None})
        scored.sort(key=lambda x: (-(x["n_target"]), x["pvalue"] or 1))
        validation[sig] = {"n_targets_mapped": len(tids), "top": scored[:8]}
    (a.out / "k12_pathway_validation.json").write_text(json.dumps({"targets": validation}, indent=1), encoding="utf-8")
    (a.out / "crosswalk_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(json.dumps(summary, indent=1))
    print(json.dumps(validation, indent=1)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
