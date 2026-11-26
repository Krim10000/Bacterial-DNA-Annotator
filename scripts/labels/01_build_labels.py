#!/usr/bin/env python3
"""Prototipo de etiquetas de proceso para K-12 y B. subtilis."""
from __future__ import annotations
import argparse, csv, json, sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from webapp.labels import labels_for_text, is_regulator  # noqa: E402


def rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def k12(data_root: Path):
    genome = "ecoli_k12"
    genes = {r["gene_id"]: r for r in rows(data_root / "regulation" / genome / "genes.tsv")}
    ops = rows(data_root / "regulation" / genome / "operons.tsv")
    out = []
    for op in ops:
        ids = [x for x in (op.get("gene_ids") or "").split(",") if x]
        labels = set()
        for gid in ids:
            g = genes.get(gid, {})
            if is_regulator(g.get("gene_name", ""), g.get("product", "")):
                continue
            labels.update(labels_for_text(f"{g.get('gene_name','')} {g.get('product','')}"))
        out.append({"operon_id": op["operon_id"], "label": op.get("label", ""),
                    "labels": sorted(labels)})
    return out


def bsu(data_root: Path):
    genome = "subtilis"
    hits = {}
    fp = data_root.parent / "results" / genome / f"{genome}.functional_annotation.json"
    if fp.exists():
        hits = json.loads(fp.read_text()).get("hits", {})
    ops = rows(data_root.parent / "results" / genome / f"{genome}.operons_auto.tsv")
    out = []
    for op in ops:
        ids = [x for x in (op.get("genes") or "").split(";") if x]
        labels = set()
        for cid in ids:
            h = hits.get(cid, {})
            if is_regulator(h.get("gene", ""), h.get("product", "")):
                continue
            labels.update(labels_for_text(f"{h.get('gene','')} {h.get('product','')}"))
        out.append({"operon_id": op.get("operon_id", ""), "label": op.get("label", ""),
                    "labels": sorted(labels)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "labels_pilot")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for genome, fn in (("ecoli_k12", k12), ("subtilis", bsu)):
        ops = fn(a.data_root)
        counts = Counter(l for o in ops for l in o["labels"])
        summary[genome] = {"n_operons": len(ops),
                           "labels": {k: {"n": v,
                                          "fraction": round(v / len(ops), 4) if ops else 0}
                                      for k, v in sorted(counts.items())}}
        with (a.out / f"{genome}_operon_labels.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["operon_id", "label", "labels"], delimiter="\t")
            w.writeheader()
            for o in ops:
                w.writerow({**o, "labels": ",".join(o["labels"])})
        print(json.dumps(summary[genome], indent=1), flush=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
