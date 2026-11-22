#!/usr/bin/env python3
"""boundary_classifier.py — ¿termina el transcrito en el borde 3' de este gen?

Encuadre
--------
No hace falta buscar terminadores en 4 Mb: hay un candidato por gen, en su borde
3'. La pregunta es binaria:

    ¿el transcrito TERMINA aquí, o sigue (read-through) al gen siguiente?

Ese "sí" es a la vez el final de la unidad de transcripción y la frontera del
operón. Clasificar por familia evita mezclar señales que no son comparables:

- **intrínseca**: horquilla clase A + poli-U detrás del gen
- **Rho**: ventana aguas arriba C-rica y pobre en G (Rho no reconoce secuencia,
  se engancha en regiones C-ricas sin estructura)
- **sin evidencia**: read-through

Uso
---
    # evaluar contra una curaduría de unidades de transcripción
    python scripts/regulation/boundary_classifier.py --gbff data/test_ecoli.gbff \
        --tu data/regulondb/gt/gt_operon_tu_detail.tsv --mode evaluate

    # correr sobre un genoma sin curaduría
    python scripts/regulation/boundary_classifier.py --gbff <genoma.gbff> \
        --mode scan --out <salida.tsv>
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))

C_FRAC_MIN = 0.27      # calibrado: rho-dependientes 0,274 vs independientes 0,238
G_FRAC_MAX = 0.26
TAIL = 150             # pb detrás del gen donde se busca la horquilla


def load_genes(gbff: Path):
    from Bio import SeqIO
    rec = next(SeqIO.parse(str(gbff), "genbank"))
    seq = str(rec.seq).upper()
    genes = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        st = "+" if f.location.strand == 1 else "-"
        tag = (f.qualifiers.get("locus_tag") or f.qualifiers.get("gene")
               or [f"{int(f.location.start)}"])[0]
        for p in f.location.parts:
            genes.append({"tag": tag, "start": int(p.start), "end": int(p.end),
                          "strand": st})
    genes.sort(key=lambda g: g["start"])
    return seq, genes


def hairpin_after(seq: str, end: int, strand: str, tail: int = TAIL):
    """Horquilla A con poli-U en la cola 3' del gen."""
    import annotate_hairpins as H
    if strand == "+":
        chunk = seq[end:end + tail]
    else:
        chunk = seq[max(0, end - tail):end][::-1].translate(str.maketrans("ACGT", "TGCA"))
    if len(chunk) < 20:
        return None
    hp = H.best_hairpin(chunk, use_dg=False)
    if hp and hp.get("cls") == "A_fuerte":
        return hp
    return None


def upstream_window(seq: str, end: int, strand: str, w: int = 100) -> str:
    if strand == "+":
        return seq[max(0, end - w):end]
    return seq[end:end + w][::-1].translate(str.maketrans("ACGT", "TGCA"))


def classify(seq: str, g: dict) -> dict:
    hp = hairpin_after(seq, g["end"], g["strand"])
    win = upstream_window(seq, g["end"], g["strand"])
    n = max(1, len(win))
    c = win.count("C") / n
    gg = win.count("G") / n
    if hp:
        call, why = "terminator_intrinsic", f"horquilla {hp['cls']} poli-U {hp['polyU']}"
    elif c >= C_FRAC_MIN and gg <= G_FRAC_MAX:
        call, why = "terminator_rho", f"C={c:.2f} G={gg:.2f}: región C-rica"
    else:
        call, why = "read_through", f"sin horquilla A; C={c:.2f} G={gg:.2f}"
    return {"tag": g["tag"], "strand": g["strand"], "end": g["end"],
            "c_frac": round(c, 3), "g_frac": round(gg, 3),
            "call": call, "why": why}


def tu_ends(path: Path) -> set[str]:
    """Último gen (3') de cada unidad de transcripción curada."""
    ends = set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
        if genes:
            ends.add(genes[-1])
    return ends


def evaluate(gbff: Path, tu_path: Path) -> None:
    seq, genes = load_genes(gbff)
    truth = tu_ends(tu_path)
    rows = [classify(seq, g) for g in genes if g["tag"] in truth or True]
    tp = fp = fn = tn = 0
    for r in rows:
        pred = r["call"].startswith("terminator")
        real = r["tag"] in truth
        if pred and real:
            tp += 1
        elif pred:
            fp += 1
        elif real:
            fn += 1
        else:
            tn += 1
    P = tp / (tp + fp) if tp + fp else 0
    R = tp / (tp + fn) if tp + fn else 0
    F = 2 * P * R / (P + R) if P + R else 0
    print(json.dumps({"n_genes": len(rows), "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                      "precision": round(P, 3), "recall": round(R, 3),
                      "f1": round(F, 3)}, indent=1))


def scan(gbff: Path, out: Path) -> None:
    seq, genes = load_genes(gbff)
    rows = [classify(seq, g) for g in genes]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    c = Counter(r["call"] for r in rows)
    print(f"{len(rows)} candidatos -> {out}")
    for k, v in c.most_common():
        print(f"  {k:22s} {v} ({100*v/len(rows):.0f}%)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["evaluate", "scan"], default="evaluate")
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--tu", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.mode == "evaluate":
        evaluate(Path(a.gbff), Path(a.tu))
    else:
        scan(Path(a.gbff), Path(a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
