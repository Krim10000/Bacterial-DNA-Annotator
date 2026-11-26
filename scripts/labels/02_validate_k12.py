#!/usr/bin/env python3
"""Valida etiquetas de proceso contra el sigmulón curado de RegulonDB K-12."""
from __future__ import annotations
import argparse, csv, json, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from webapp.labels import LABELS  # noqa: E402

EXPECTED = {"sigma32": "heat_shock", "sigma28": "flagellum",
            "sigma54": "nitrogen", "sigma38": "stationary",
            "sigma24": "envelope"}


def read(path):
    with path.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--labels", type=Path, default=ROOT / "reports" / "labels_pilot" / "ecoli_k12_operon_labels.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "labels_pilot")
    a = ap.parse_args()
    lab = {r["operon_id"]: set(x for x in (r["labels"] or "").split(",") if x)
           for r in read(a.labels)}
    op_rows = read(a.data_root / "regulation" / "ecoli_k12" / "operons.tsv")
    by_label = {r["label"]: r["operon_id"] for r in op_rows}
    net = read(a.data_root / "regulondb" / "NetworkSigmaTU.tsv")
    targets = {}
    for r in net:
        sigma = (r.get("1)sigmaName") or "").strip()
        tu = re.split(r"[\[\s]", (r.get("2)regulatedTUName") or "").strip())[0]
        if sigma and tu:
            targets.setdefault(sigma, set()).add(tu)
    out = {"n_operons": len(op_rows), "by_sigma": {}, "labels": sorted(LABELS)}
    for sigma, label in EXPECTED.items():
        tus = targets.get(sigma, set())
        mapped = {by_label[t] for t in tus if t in by_label}
        labeled = {oid for oid in mapped if label in lab.get(oid, set())}
        # Fondo: todos los operones con etiqueta.
        all_with = sum(1 for oid in by_label.values() if label in lab.get(oid, set()))
        non_target = set(by_label.values()) - mapped
        a1 = len(labeled)
        b1 = len(mapped) - a1
        c1 = sum(1 for oid in non_target if label in lab.get(oid, set()))
        d1 = len(non_target) - c1
        try:
            from scipy.stats import fisher_exact
            odds, pval = fisher_exact([[a1, b1], [c1, d1]], alternative="greater")
        except Exception:
            odds, pval = None, None
        out["by_sigma"][sigma] = {
            "expected_label": label,
            "n_targets": len(tus),
            "n_mapped_operons": len(mapped),
            "n_labeled": len(labeled),
            "recall": round(len(labeled) / len(mapped), 4) if mapped else None,
            "background_fraction": round(all_with / len(by_label), 4) if by_label else None,
            "odds_ratio": round(float(odds), 3) if odds is not None else None,
            "pvalue": round(float(pval), 6) if pval is not None else None,
        }
    # FDR Benjamini-Hochberg sobre los sigma evaluados.
    vals = [v for v in out["by_sigma"].values() if v.get("pvalue") is not None]
    order = sorted(range(len(vals)), key=lambda i: vals[i]["pvalue"])
    m = len(vals)
    for rank, idx in enumerate(order, 1):
        vals[idx]["fdr"] = round(min(1.0, vals[idx]["pvalue"] * m / rank), 6)
    (a.out / "k12_validation.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lines = ["# Validación de etiquetas K-12", "",
             f"Operones: {out['n_operons']}", "",
             "| Sigma | Etiqueta | Diana | Mapeados | Con etiqueta | Recall | Fondo | OR | FDR |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for s, v in out["by_sigma"].items():
        lines.append(f"| {s} | {v['expected_label']} | {v['n_targets']} | "
                     f"{v['n_mapped_operons']} | {v['n_labeled']} | "
                     f"{v['recall']} | {v['background_fraction']} | "
                     f"{v.get('odds_ratio')} | {v.get('fdr')} |")
    (a.out / "k12_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
