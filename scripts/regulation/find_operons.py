#!/usr/bin/env python3
"""find_operons.py — encontrar operones en un genoma NO anotado, desde la secuencia.

Por qué así
-----------
En un genoma sin anotar no hay componentes que leer: hay que **predecirlos**.
Solo *E. coli* K-12 tiene componentes curados, y se usa para dos cosas:
calibrar los detectores y entrenar el clasificador de fronteras. La aplicación
es al genoma 401.

Componentes de novo
-------------------
- `terminator_score`: horquilla GC + tracto poli-T (terminador intrínseco).
  Calibrado contra los 375 terminadores curados (`--validate`).
- `promoter_score`: PWM de σ70 −10/−35 sobre la región río arriba.
- `rbs_score`: motivo Shine-Dalgarno en la ventana de inicio.

Fronteras
---------
Regresión logística sobre [log distancia, solapamiento, terminador, promotor,
RBS, misma EC] entrenada en K-12; se corta donde P(frontera) ≥ umbral.

No lee `data/regulation/`: es un método independiente, comparable con L2.

Uso
---
    # validar los detectores de componentes contra RegulonDB
    python scripts/regulation/find_operons.py --validate

    # entrenar en K-12 y predecir
    python scripts/regulation/find_operons.py \
        --gbff data/test_ecoli.gbff --genome ecoli_k12 --train

    # aplicar a un genoma no anotado (GFF de pyrodigal + FASTA)
    python scripts/regulation/find_operons.py \
        --gff data/annotations_ab_initio/<GCF>/annotation.gff3 \
        --fasta data/raw/<GCF>_genomic.fna.gz --genome <GCF> --apply
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import pickle
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROM_PWM = ROOT / "data" / "regulondb" / "gt" / "promoter_pwms.tsv"
GT_OPERONS = ROOT / "data" / "regulondb" / "gt" / "gt_operons.tsv"
GT_TUS = ROOT / "data" / "regulondb" / "gt" / "gt_tus.tsv"
TERMINATORS = ROOT / "data" / "regulondb" / "TerminatorSet.tsv"
BASES = "ACGT"
RBS_PAT = re.compile(r"AGGAGG|GGAGG|AGGAG|AGGA")
COMP = str.maketrans("ACGT", "TGCA")
FEATURES = ["log_dist", "overlap", "terminator", "term_dist", "promoter",
            "promoter_adj", "prom_dist", "rbs", "same_ec", "string",
            "tts_near", "tss_near"]
TERMINI: dict[str, dict[str, list[int]]] = {"+": {"tts": [], "tss": []},
                                            "-": {"tts": [], "tss": []}}
# Los extremos experimentales NO mejoraron el F1 (0.807 vs 0.815 con STRING
# solo), así que quedan opcionales: cubren 61% (TSS) y 16% (TTS) de las
# fronteras y un TSS interno no implica frontera de operón.
_WITH_TERMINI = False
STRING_SCORES: dict[tuple[str, str], float] = {}
# prior empírico de la distancia promotor→inicio (pb), de gt_promoters σ70
PROM_BINS = [0, 20, 40, 60, 80, 120, 200, 300, 10 ** 9]
PROM_PRIOR = None


# ---------------------------------------------------------------- componentes
def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def hairpin_scan(seq: str, mismatches: int = 1, min_stem: int = 4,
                 loop_min: int = 2, loop_max: int = 10) -> tuple[float, int]:
    """Horquilla: tallo GC rico + bucle corto, tolerando 1 error.

    Exigir coincidencia exacta y además el tracto poli-T pegado bajaba el AUC a
    0.44 (peor que azar). Separadas, la horquilla sola da ~0.93 y el poli-T ~0.76.
    """
    seq = seq.upper()
    n = len(seq)
    best, best_i = 0.0, 0
    for i in range(n):
        for stem in range(min(12, n - i), min_stem - 1, -1):
            head = seq[i:i + stem]
            if any(b not in BASES for b in head):
                continue
            for loop in range(loop_min, loop_max + 1):
                j = i + stem + loop
                if j + stem > n:
                    continue
                tail = seq[j:j + stem]
                if any(b not in BASES for b in tail):
                    continue
                mm = sum(1 for x, y in zip(head, revcomp(tail)) if x != y)
                if mm > mismatches:
                    continue
                gc = sum(1 for b in head if b in "GC") / stem
                if stem + 3 * gc > best:
                    best, best_i = stem + 3 * gc, i
    return best, best_i


def poly_t_score(seq: str) -> float:
    run = best = 0
    for b in seq.upper() + "A":
        run = run + 1 if b == "T" else 0
        best = max(best, run)
    return float(best)


def terminator_scan(seq: str) -> tuple[float, int]:
    """(puntuación, posición de la horquilla)."""
    if not seq:
        return 0.0, 0
    hair, pos = hairpin_scan(seq)
    return round(hair + 1.5 * poly_t_score(seq), 2), pos


def load_promoter_pwm(sigma: str = "sigma70") -> list[dict]:
    if not PROM_PWM.exists():
        return []
    out = []
    for r in csv.DictReader(open(PROM_PWM), delimiter="\t"):
        if r["sigma"] == sigma:
            out.append({"box": r["box"], "width": int(r["width"]),
                        "pwm": json.loads(r["pwm"])})
    return out


def pwm_score(pwm: list[list[float]], seq: str) -> float:
    if len(seq) != len(pwm):
        return float("-inf")
    return sum(pwm[i][BASES.index(b)] if b in BASES else -6.0
               for i, b in enumerate(seq))


def promoter_scan(seq: str, boxes: list[dict],
                  spacer: tuple[int, int] = (14, 21)) -> tuple[float, int]:
    b35 = next((b for b in boxes if b["box"] == "-35"), None)
    b10 = next((b for b in boxes if b["box"] == "-10"), None)
    if not seq:
        return 0.0, 10 ** 6
    # matriz de una sola caja: hay organismos (B. subtilis) cuyos promotores
    # housekeeping no conservan caja -35; entonces se puntúa solo la que existe.
    if not b35 or not b10:
        box = b35 or b10
        w, pwm = box["width"], box["pwm"]
        up = seq.upper()
        best, best_dist = 0.0, 10 ** 6
        for i in range(0, max(0, len(up) - w) + 1):
            sc = pwm_score(pwm, up[i:i + w])
            if sc > best:
                best, best_dist = sc, len(up) - (i + w)
        return best, best_dist
    seq = seq.upper()
    best, best_dist = 0.0, 10 ** 6
    w35, w10 = b35["width"], b10["width"]
    for i in range(0, max(0, len(seq) - w35) + 1):
        s35 = pwm_score(b35["pwm"], seq[i:i + w35])
        if s35 < 0:
            continue
        for sp in range(spacer[0], spacer[1] + 1):
            j = i + w35 + sp
            if j + w10 > len(seq):
                break
            s10 = pwm_score(b10["pwm"], seq[j:j + w10])
            if s10 > 0 and s35 + s10 > best:
                best = s35 + s10
                best_dist = len(seq) - (j + w10)   # -10 -> inicio del gen
    return best, best_dist


def promoter_prob(score: float, dist: int) -> float:
    """Score de la PWM + log-prior empírico de la posición."""
    if PROM_PRIOR is None:
        return score
    import bisect
    idx = min(bisect.bisect_right(PROM_BINS, dist) - 1, len(PROM_PRIOR) - 1)
    return score + PROM_PRIOR[max(0, idx)]


def rbs_score(upstream: str) -> float:
    m = RBS_PAT.search(upstream.upper()[-20:])
    return float(len(m.group(0))) if m else 0.0


# --------------------------------------------------------------------- lectura
def read_gbff(path: Path) -> tuple[str, list[dict]]:
    """Secuencia + CDS (con coordenadas) desde un GenBank, sin Biopython."""
    seq_parts: list[str] = []
    genes: list[dict] = []
    cur: dict | None = None
    loc = ""
    in_seq = False
    with open(path, "rt") as fh:
        for line in fh:
            if line.startswith("ORIGIN"):
                in_seq = True
                continue
            if in_seq:
                if line.startswith("//"):
                    in_seq = False
                    continue
                seq_parts.append("".join(c for c in line if c.isalpha()))
                continue
            if line.startswith("     ") and len(line) > 5 and line[5] != " ":
                if cur and cur.get("start"):
                    genes.append(cur)
                key = line[5:21].strip()
                if key == "CDS":
                    cur = {"locus_tag": "", "gene": "", "ec": "", "strand": "+",
                           "start": 0, "end": 0, "seqid": ""}
                    loc = line[21:].strip()
                else:
                    cur, loc = None, ""
                continue
            if cur is None:
                continue
            body = line[21:].strip()
            if body.startswith("/"):
                q = body[1:]
                if "=" in q:
                    k, v = q.split("=", 1)
                    v = v.strip().strip('"')
                    if k == "locus_tag":
                        cur["locus_tag"] = v
                    elif k == "ID" and not cur["locus_tag"]:
                        # GenBank derivado de GFF: usa /ID en vez de /locus_tag
                        cur["locus_tag"] = v.strip('"')
                    elif k == "gene" and not cur["gene"]:
                        cur["gene"] = v
                    elif k == "EC_number":
                        cur["ec"] = v
            elif body and not body.startswith("/"):
                loc += body
            if ".." in loc and not cur["start"]:
                m = re.search(r"(\d+)\.\.(\d+)", loc)
                if m:
                    cur["start"], cur["end"] = int(m.group(1)), int(m.group(2))
                    cur["strand"] = "-" if "complement" in loc else "+"
    if cur and cur.get("start"):
        genes.append(cur)
    return "".join(seq_parts).upper(), genes


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = defaultdict(list)
    name = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def read_gff_genes(path: Path) -> list[dict]:
    genes = []
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        genes.append({"locus_tag": attrs.get("ID", ""), "gene": "", "ec": "",
                      "seqid": p[0], "start": int(p[3]), "end": int(p[4]),
                      "strand": p[6]})
    return genes


# ------------------------------------------------------------------- features
def ordered(genes: list[dict]) -> list[dict]:
    return sorted(genes, key=lambda g: (g.get("seqid", ""), g["start"], g["end"]))


def pairs(genes: list[dict]):
    """Pares consecutivos en la misma hebra (candidatos a frontera)."""
    for a, b in zip(genes, genes[1:]):
        if a["strand"] == b["strand"] and a.get("seqid") == b.get("seqid"):
            yield a, b


def features_of(a: dict, b: dict, seq: str, boxes: list[dict]) -> list[float]:
    gap = b["start"] - a["end"] - 1
    inter = seq[a["end"]:b["start"]] if gap > 0 else ""
    up = seq[max(0, b["start"] - 120):max(0, b["start"] - 1)]
    same_ec = 1.0 if (a.get("ec") and a.get("ec") == b.get("ec")) else 0.0
    term, tpos = terminator_scan(inter) if inter else (0.0, 0)
    pscore, pdist = promoter_scan(up, boxes)
    return [math.log1p(max(0, gap)), 1.0 if gap < 0 else 0.0,
            term, float(tpos),
            pscore, promoter_prob(pscore, pdist), float(pdist),
            rbs_score(seq[max(0, b["start"] - 20):max(0, b["start"] - 1)]),
            same_ec,
            STRING_SCORES.get((a["locus_tag"], b["locus_tag"]), 0.0),
            _near("tts", a["strand"], a["end"]) if _WITH_TERMINI else 0.0,
            _near("tss", b["strand"], b["start"]) if _WITH_TERMINI else 0.0]


def load_operon_of() -> dict[str, str]:
    out = {}
    if GT_OPERONS.exists():
        for r in csv.DictReader(open(GT_OPERONS), delimiter="\t"):
            for g in r["genes"].split(","):
                if g:
                    out[g] = r["operon_id"]
    return out


def segment(genes: list[dict], probs: list[float], threshold: float) -> list[list[dict]]:
    ops: list[list[dict]] = []
    cur: list[dict] = []
    it = iter(probs)
    for a, b in zip(genes, genes[1:]):
        same = a["strand"] == b["strand"] and a.get("seqid") == b.get("seqid")
        if not same:
            if cur:
                ops.append(cur)
            cur = [b]
            continue
        if not cur:
            cur = [a]
        p = next(it, 1.0)
        if p >= threshold:
            ops.append(cur)
            cur = [b]
        else:
            cur.append(b)
    if cur:
        ops.append(cur)
    return ops


def pair_f1(pred_ops: list[list[dict]], op_of: dict[str, str]) -> tuple[float, float, float]:
    """F1 de pares génicos contra el catálogo curado."""
    pairs_pred = set()
    for op in pred_ops:
        tags = [g["locus_tag"] for g in op if g["locus_tag"] in op_of]
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                pairs_pred.add(frozenset((tags[i], tags[j])))
    by_op: dict[str, set[str]] = defaultdict(set)
    for tag, oid in op_of.items():
        by_op[oid].add(tag)
    pairs_gt = set()
    for tags in by_op.values():
        t = sorted(tags)
        for i in range(len(t)):
            for j in range(i + 1, len(t)):
                pairs_gt.add(frozenset((t[i], t[j])))
    tp = len(pairs_pred & pairs_gt)
    p = tp / len(pairs_pred) if pairs_pred else 0.0
    r = tp / len(pairs_gt) if pairs_gt else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return round(p, 3), round(r, 3), round(f, 3)


def sweep_thresholds(genes: list[dict], probs: list[float],
                     op_of: dict[str, str]) -> float:
    print("barrido de umbral (pares):")
    best_f, best_t = -1.0, 0.5
    for t in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
        p, r, f = pair_f1(segment(genes, probs, t), op_of)
        mark = " <" if f > best_f else ""
        print(f"  t={t:.1f}  P={p:.3f} R={r:.3f} F1={f:.3f}{mark}")
        if f > best_f:
            best_f, best_t = f, t
    return best_t


def write_operons(ops: list[list[dict]], genome: str, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["operon_id", "genome", "contig", "strand", "start", "end",
                    "n_genes", "genes", "gene_names", "confidence", "status"])
        for i, op in enumerate(ops, 1):
            w.writerow([f"{genome}:op{i}", genome, op[0].get("seqid", ""),
                        op[0]["strand"], min(g["start"] for g in op),
                        max(g["end"] for g in op), len(op),
                        ",".join(g["locus_tag"] for g in op),
                        ",".join(g.get("gene", "") or "" for g in op),
                        "predicted", "candidate"])


# ------------------------------------------------------------------- validar
def load_termini(root: Path) -> None:
    """Extremos experimentales (TTS 3', TSS 5') para usar como evidencia."""
    import bisect
    for kind, fname, col in (("tts", "ht_tts.tsv", "rightEndPosition"),
                             ("tss", "ht_tss.tsv", "pos_1")):
        path = root / "regulondb" / "ht" / fname
        if not path.exists():
            print(f"sin {fname}: feature {kind} = 0")
            continue
        pos = {"+": set(), "-": set()}
        for r in csv.DictReader(open(path), delimiter="\t"):
            try:
                p = int(float(r.get(col) or r.get("leftEndPosition")))
            except (ValueError, TypeError):
                continue
            k = "+" if r.get("strand") in ("+", "1", "forward") else "-"
            pos[k].add(p)
        for k in pos:
            TERMINI[k][kind] = sorted(pos[k])
        print(f"{kind}: {len(pos['+'])} en +, {len(pos['-'])} en -")


def _near(kind: str, strand: str, target: int, tol: int = 150) -> float:
    """1.0 si hay un extremo experimental a <= tol pb del punto."""
    import bisect
    arr = TERMINI.get(strand, {}).get(kind, [])
    if not arr:
        return 0.0
    i = bisect.bisect_left(arr, target - tol)
    for x in arr[i:i + 4]:
        if abs(x - target) <= tol:
            return 1.0
    return 0.0


def load_string(path: Path) -> None:
    """scores de asociación funcional STRING por par (Taboada 2010: 94,6%)."""
    global STRING_SCORES
    if not path.exists():
        print("sin cache de STRING; feature = 0 para todos")
        return
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                v = float(r["score"])
            except (ValueError, KeyError):
                continue
            STRING_SCORES[(r["pair_a"], r["pair_b"])] = v
            STRING_SCORES[(r["pair_b"], r["pair_a"])] = v
    print(f"STRING: {len(STRING_SCORES)//2} pares con score")


def load_prom_prior() -> None:
    """log-prior empírico de la distancia -10 → inicio del gen (σ70)."""
    global PROM_PRIOR
    import math as _m
    src = ROOT / "data" / "regulondb" / "gt" / "gt_promoters.tsv"
    if not src.exists():
        PROM_PRIOR = None
        return
    counts = [1] * (len(PROM_BINS) - 1)          # Laplace
    import bisect
    for r in csv.DictReader(open(src), delimiter="\t"):
        if r.get("sigma_factor") != "sigma70":
            continue
        try:
            d = abs(int(float(r["dist_to_first_gene"]))) + 7  # TSS -> -10 aprox
        except (ValueError, KeyError):
            continue
        idx = min(bisect.bisect_right(PROM_BINS, d) - 1, len(counts) - 1)
        counts[max(0, idx)] += 1
    total = sum(counts)
    PROM_PRIOR = [_m.log(c / total) for c in counts]
    print("prior de posición promotor:", [round(x, 2) for x in PROM_PRIOR])


def validate_components() -> int:
    """Comprueba los detectores contra los datos curados de RegulonDB."""
    rnd = random.Random(42)
    rows = []
    with open(TERMINATORS, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = [c.split(")")[-1] for c in next(rd)]
        for r in rd:
            d = dict(zip(header, r))
            seq = d.get("sequence", "").strip().upper()
            if seq:
                rows.append(seq)
    if not rows:
        print("sin terminadores curados")
        return 1
    real = sorted(terminator_scan(s)[0] for s in rows)
    pool = "".join(rows)
    fake = sorted(terminator_scan("".join(rnd.choice(pool) for _ in range(len(rows[0]))))[0]
                  for _ in range(len(rows)))
    auc = sum(1 for r in real for f in fake if r > f) / (len(real) * len(fake))
    thr = sorted(real)[int(len(real) * 0.5)]
    sens = sum(1 for r in real if r >= 6) / len(real)
    spec = sum(1 for f in fake if f < 6) / len(fake)
    print(f"terminador: {len(rows)} curados")
    print(f"  score mediano real {sorted(real)[len(real)//2]:.1f} vs aleatorio {sorted(fake)[len(fake)//2]:.1f}")
    print(f"  AUC {auc:.3f} | corte >=6: sensibilidad {sens:.0%}, especificidad {spec:.0%}")
    return 0


# ---------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gbff", type=Path, help="GenBank con secuencia (K-12)")
    ap.add_argument("--gff", type=Path, help="GFF3 de CDS (genoma no anotado)")
    ap.add_argument("--fasta", type=Path, help="FASTA del genoma")
    ap.add_argument("--genome", default="genome")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "operon_finder")
    ap.add_argument("--train", action="store_true")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--with-termini", action="store_true",
                    help="añade features TTS/TSS experimentales (no mejora: 0.807 vs 0.815)")
    ap.add_argument("--sweep", action="store_true",
                    help="barre umbrales y elige el de mejor F1 de pares")
    a = ap.parse_args()

    if a.validate:
        return validate_components()

    load_prom_prior()
    load_string(a.out / "string_scores.tsv")
    global _WITH_TERMINI
    _WITH_TERMINI = a.with_termini
    if a.with_termini:
        load_termini(ROOT / "data")

    boxes = load_promoter_pwm()
    print(f"PWM σ70: {[b['box'] for b in boxes]}")
    if not boxes:
        print("faltan PWM: corre scripts/validation/build_tf_pwms.py --promoters",
              file=sys.stderr)
        return 1

    if a.gbff:
        seq, genes = read_gbff(a.gbff)
        genes = ordered([g for g in genes if g.get("locus_tag")])
        print(f"{a.genome}: {len(seq):,} pb, {len(genes)} CDS")
    else:
        if not (a.gff and a.fasta):
            print("para un genoma no anotado hacen falta --gff y --fasta",
                  file=sys.stderr)
            return 1
        seqs = read_fasta(a.fasta)
        genes = []
        for g in read_gff_genes(a.gff):
            g["_seqid"] = g["seqid"]
            genes.append(g)
        # se trabaja por contig: se concatenan con separadores para conservar
        # coordenadas locales por contig
        genes = ordered(genes)
        print(f"{a.genome}: {len(genes)} CDS en {len(seqs)} contigs")

    model_path = a.out / "boundary_model.pkl"
    if a.train:
        from sklearn.linear_model import LogisticRegression
        op_of = load_operon_of()
        X, y = [], []
        all_pairs = list(pairs(genes))
        for a1, b1 in all_pairs:
            # solo fronteras con etiqueta conocida: si un gen no está en el
            # catálogo curado, la etiqueta sería inventada
            if a1["locus_tag"] not in op_of or b1["locus_tag"] not in op_of:
                continue
            X.append(features_of(a1, b1, seq, boxes))
            y.append(1 if op_of[a1["locus_tag"]] != op_of[b1["locus_tag"]] else 0)
        print(f"fronteras: {len(y)} ({sum(y)} cortes, {len(y) - sum(y)} uniones)")
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        clf.fit(X, y)
        coef = dict(zip(FEATURES, [round(c, 3) for c in clf.coef_[0]]))
        print("coeficientes:", coef)
        a.out.mkdir(parents=True, exist_ok=True)
        pickle.dump(clf, open(model_path, "wb"))
        # IMPORTANTE: predecir sobre TODAS las fronteras, incluidas las que no
        # tienen etiqueta; si no, la lista queda desalineada con la segmentación
        X_all = [features_of(a1, b1, seq, boxes) for a1, b1 in all_pairs]
        probs = [clf.predict_proba([x])[0][1] for x in X_all]
        if a.sweep:
            best = sweep_thresholds(genes, probs, op_of)
            a.threshold = best
            print(f"umbral elegido por barrido: {best}")
        # probabilidad por frontera + etiqueta, para calibrar la cascada
        with open(a.out / f"{a.genome}_boundary_probs.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pair_a", "pair_b", "prob", "label", "log_dist",
                        "terminator", "promoter", "same_ec", "string"])
            for (x, z), pr, row in zip(all_pairs, probs, X_all):
                lab = ""
                if x["locus_tag"] in op_of and z["locus_tag"] in op_of:
                    lab = int(op_of[x["locus_tag"]] != op_of[z["locus_tag"]])
                w.writerow([x["locus_tag"], z["locus_tag"], f"{pr:.4f}", lab]
                           + [f"{v:.4f}" for v in row])
        ops = segment(genes, probs, a.threshold)
        out = a.out / f"{a.genome}_operons.tsv"
        write_operons(ops, a.genome, out)
        print(f"operones predichos: {len(ops)} -> {out.relative_to(ROOT)}")
        return 0

    if a.apply:
        if not model_path.exists():
            print(f"falta {model_path}; corre --train primero", file=sys.stderr)
            return 1
        clf = pickle.load(open(model_path, "rb"))
        # `seqs` solo se define en la ruta GFF+FASTA; con --gbff hay una sola
        # secuencia en `seq`. Se unifica aquí para que ambas rutas funcionen.
        seqs_map = locals().get("seqs") or {}
        if not seqs_map and seq:
            seqs_map = {g.get("seqid", "") or "": seq for g in genes}
        # multi-contig: cada contig con SU secuencia
        by_contig: dict[str, list[dict]] = defaultdict(list)
        for g in genes:
            by_contig[g.get("seqid", "")].append(g)
        ops: list[list[dict]] = []
        skipped = 0
        for seqid, grp in by_contig.items():
            s = seqs_map.get(seqid, "")
            if not s:
                skipped += len(grp)
                continue
            grp = ordered(grp)
            probs = [clf.predict_proba([features_of(x, z, s, boxes)])[0][1]
                     for x, z in pairs(grp)]
            ops.extend(segment(grp, probs, a.threshold))
        out = a.out / f"{a.genome}_operons.tsv"
        write_operons(ops, a.genome, out)
        print(f"operones: {len(ops)} -> {out.relative_to(ROOT)}"
              + (f" | CDS sin secuencia: {skipped}" if skipped else ""))
        return 0

    print("usa --validate, --train o --apply")
    return 1


if __name__ == "__main__":
    sys.exit(main())
