#!/usr/bin/env python3
"""boundary_prior.py — clasificador de fronteras con prior calibrado.

El error del clasificador anterior
----------------------------------
Exigía evidencia local para llamar terminador, cuando en K-12 el 65% de los
genes ya son final de TU. La consecuencia: F1 0,473, **por debajo** de decir
"termina" siempre (0,786). Un modelo que no supera al predictor constante no
aporta nada.

Modelo
------
Se parte del prior y se actualiza con evidencia, en vez de exigirla:

    P(termina) = 0,65  →  × LR(horquilla)  × LR(rho)  × LR(hueco)  × LR(nombre)

- **horquilla A + poli-U**: evidencia fuerte a favor de terminar.
- **región C-rica (Rho)**: evidencia moderada.
- **hueco grande al gen siguiente**: a favor de terminar.
- **mismo prefijo de nombre con el gen siguiente**: a favor de read-through.

El modelo se ajusta sin rebalanceo de clases, para que el prior no se distorsione.
Se evalúa contra el baseline constante en el **mismo conjunto de prueba**.

Uso
---
    python scripts/regulation/boundary_prior.py --gbff data/test_ecoli.gbff \
        --tu data/regulondb/gt/gt_operon_tu_detail.tsv
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

TAIL = 150


def load(gbff: Path):
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
        name = (f.qualifiers.get("gene") or [tag])[0]
        for p in f.location.parts:
            genes.append({"tag": tag, "name": name, "start": int(p.start),
                          "end": int(p.end), "strand": st})
    genes.sort(key=lambda g: g["start"])
    return seq, genes


def tu_ends(path: Path) -> set[str]:
    out = set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        g = [x.strip() for x in (r.get("genes") or "").split(",") if x.strip()]
        if g:
            out.add(g[-1])
    return out


TSS: dict[str, list[int]] = {"+": [], "-": []}
EC: dict[str, str] = {}
SYN: dict[tuple[str, str], float] = {}


def load_ec(path: Path) -> None:
    if not path.exists():
        return
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r.get("locus_tag") and r.get("ec"):
            EC[r["locus_tag"]] = r["ec"]


def load_syn(path: Path) -> None:
    if not path.exists():
        return
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            SYN[(r["left"], r["right"])] = float(r["frac"])
        except (KeyError, ValueError):
            continue


def load_tss(path: Path) -> None:
    import bisect
    for r in csv.DictReader(open(path), delimiter="\t"):
        v = (r.get("pos") or "").strip()
        if not v.isdigit():
            continue
        st = "-" if (r.get("strand") or "+").strip().startswith("-") else "+"
        TSS[st].append(int(v))
    for v in TSS.values():
        v.sort()


def tss_between(a: int, b: int, strand: str) -> bool:
    import bisect
    arr = TSS.get(strand, [])
    i = bisect.bisect_left(arr, a)
    return i < len(arr) and arr[i] <= b


def feats(seq: str, genes: list[dict], i: int) -> dict:
    import annotate_hairpins as H
    g = genes[i]
    # siguiente gen en la misma hebra
    nxt = None
    for j in range(i + 1, len(genes)):
        if genes[j]["strand"] == g["strand"]:
            nxt = genes[j]
            break
    # cola 3'
    if g["strand"] == "+":
        tail = seq[g["end"]:g["end"] + TAIL]
        up = seq[max(0, g["end"] - 100):g["end"]]
    else:
        tail = seq[max(0, g["end"] - TAIL):g["end"]][::-1].translate(
            str.maketrans("ACGT", "TGCA"))
        up = seq[g["end"]:g["end"] + 100][::-1].translate(
            str.maketrans("ACGT", "TGCA"))
    hp = H.best_hairpin(tail, use_dg=False) if len(tail) >= 20 else None
    hairpin = bool(hp and hp.get("cls") == "A_fuerte" and hp.get("polyU", 0) >= 4)
    n = max(1, len(up))
    c = up.count("C") / n
    gg = up.count("G") / n
    rho = bool(c >= 0.27 and gg <= 0.26)
    if g["strand"] == "+":
        gap = (nxt["start"] - g["end"] - 1) if nxt else 10 ** 6
    else:
        gap = (g["end"] - nxt["end"] - 1) if nxt else 10 ** 6
    gap = max(-200, min(gap, 10 ** 6))
    same_name = False
    if nxt:
        a, b = g["name"].lower(), nxt["name"].lower()
        same_name = bool(a and b and (a[:4] == b[:4] or (len(a) > 3 and len(b) > 3
                                                          and a[:3] == b[:3])))
    # ¿hay un TSS documentado entre este gen y el siguiente? Si lo hay, el
    # transcrito actual termina obligatoriamente antes de él.
    tss_gap = 0
    if nxt:
        if g["strand"] == "+":
            tss_gap = int(tss_between(g["end"], nxt["start"], g["strand"]))
        else:
            tss_gap = int(tss_between(nxt["end"], g["end"], g["strand"]))
    ec_same = 0
    if nxt:
        a, b = EC.get(g["tag"], ""), EC.get(nxt["tag"], "")
        if a and b and ".".join(a.split(".")[:2]) == ".".join(b.split(".")[:2]):
            ec_same = 1
    syn = SYN.get((g["tag"], nxt["tag"]), 0.0) if nxt else 0.0
    # distancia al TSS siguiente (no solo presencia)
    import bisect
    arr = TSS.get(g["strand"], [])
    i2 = bisect.bisect_left(arr, g["end"])
    next_tss = (arr[i2] - g["end"]) if i2 < len(arr) else 10 ** 6
    return {"hairpin": int(hairpin), "rho": int(rho), "ec_same": ec_same,
            "synteny": syn, "next_tss": min(next_tss, 10 ** 6),
            "gap": gap, "same_name": int(same_name),
            "tss_between": tss_gap,
            "c_frac": round(c, 3), "g_frac": round(gg, 3)}


def main() -> int:
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split

    ap = argparse.ArgumentParser()
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--tu", required=True)
    ap.add_argument("--tss", default=str(ROOT / "data" / "regulondb" / "tss_k12.tsv"))
    ap.add_argument("--synteny", default=str(ROOT / "data" / "operon_finder" / "k12_synteny.tsv"))
    a = ap.parse_args()

    load_ec(ROOT / "data" / "regulation" / "ecoli_k12" / "genes.tsv")
    load_syn(Path(a.synteny))
    print(f"EC con anotacion: {len(EC)} | pares con sintenia: {len(SYN)}")
    if Path(a.tss).exists():
        load_tss(Path(a.tss))
        print(f"TSS cargados: {len(TSS['+'])+len(TSS['-'])}")
    seq, genes = load(Path(a.gbff))
    truth = tu_ends(Path(a.tu))
    X, y = [], []
    for i in range(len(genes) - 1):
        f = feats(seq, genes, i)
        X.append([f["hairpin"], f["rho"], f["gap"], f["same_name"],
                  f["tss_between"], f["ec_same"], f["synteny"],
                  f["next_tss"], f["c_frac"], f["g_frac"]])
        y.append(1 if genes[i]["tag"] in truth else 0)
    X = np.array(X, dtype=float)
    y = np.array(y)
    prior = y.mean()
    print(f"genes evaluados: {len(y)} | finales de TU: {int(y.sum())} "
          f"({100*prior:.0f}%)")

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3,
                                          random_state=0, stratify=y)
    clf = make_pipeline(StandardScaler(),
                        LogisticRegression(max_iter=2000)).fit(Xtr, ytr)
    pr = clf.predict_proba(Xte)[:, 1]
    pred = (pr >= 0.5).astype(int)
    tp = int(((pred == 1) & (yte == 1)).sum())
    fp = int(((pred == 1) & (yte == 0)).sum())
    fn = int(((pred == 1) & (yte == 1)).sum() == 0)  # placeholder
    fn = int(((pred == 0) & (yte == 1)).sum())
    P = tp / (tp + fp) if tp + fp else 0
    R = tp / (tp + fn) if tp + fn else 0
    F = 2 * P * R / (P + R) if P + R else 0
    base_P = yte.mean()
    base_F = 2 * base_P * 1.0 / (base_P + 1.0)
    coef = dict(zip(["hairpin", "rho", "gap", "same_name", "tss_between",
                     "ec_same", "synteny", "next_tss", "c_frac", "g_frac"],
                    [round(float(c), 3) for c in clf[-1].coef_[0]]))
    print(json.dumps({
        "modelo": {"P": round(P, 3), "R": round(R, 3), "F1": round(F, 3)},
        "baseline_constante": {"P": round(base_P, 3), "R": 1.0,
                               "F1": round(base_F, 3)},
        "gana_al_baseline": bool(F > base_F),
        "pesos": coef,
    }, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
