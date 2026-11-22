#!/usr/bin/env python3
"""Validación del motor de descubrimiento.

Control positivo
    Se ocultan K CDS anotados del GFF y se comprueba si el pipeline los
    recupera como ORFs en las regiones que quedan "vacías". Mide recall y la
    puntuación de los recuperados.

Control negativo
    Se barajan dinucleótidos de las regiones intergénicas (misma composición,
    sin estructura de codones) y se mide cuántos ORFs superan los umbrales.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
from pathlib import Path

import numpy as np

from .orfs import find_orfs_oriented, padded_regions
from .pipeline import _extract, _find_fasta
from .scoring import CodonModel, rank_score, rbs_score, start_weight
from .seqio import read_fasta, read_gff_features, revcomp


def _shuffle_dinuc(seq: str, rng: random.Random) -> str:
    """Baraja dinucleótidos preservando la composición local."""
    if len(seq) < 4:
        return seq
    pairs = [seq[i:i + 2] for i in range(0, len(seq) - 1, 2)]
    rng.shuffle(pairs)
    out = "".join(pairs)
    if len(seq) % 2:
        out += seq[-1]
    return out


def _identity(window: str, ref: str) -> float:
    return sum(1 for a, b in zip(window, ref) if a == b) / len(ref)


def peptide_match(orf_aa: str, ref_aa: str) -> float:
    """Mejor identidad de ``ref_aa`` contra ventanas de ``orf_aa``."""
    if not ref_aa:
        return 0.0
    if len(orf_aa) < len(ref_aa):
        return _identity(orf_aa.ljust(len(ref_aa), "-"), ref_aa)
    best = 0.0
    for off in range(0, len(orf_aa) - len(ref_aa) + 1):
        best = max(best, _identity(orf_aa[off:off + len(ref_aa)], ref_aa))
        if best == 1.0:
            break
    return best


def run(data_root: Path, genome: str, k: int, seed: int, out: Path) -> dict:
    rng = random.Random(seed)
    fasta = _find_fasta(data_root, genome)
    contigs = dict(read_fasta(fasta))
    _, feats = read_gff_features(
        data_root / "annotations_ab_initio" / genome / "annotation.gff3"
    )
    cds = [f for f in feats if f[4] == "CDS" and 300 <= f[2] - f[1] + 1 <= 1500]
    rng.shuffle(cds)
    held = cds[:k]
    held_keys = {(c, s, e, st) for c, s, e, st, _ in held}
    kept = [f for f in feats if (f[0], f[1], f[2], f[3]) not in held_keys]

    # Verdaderos péptidos (sin codón de parada final).
    truth = {}
    for c, s, e, st, _ in held:
        nt = _extract(contigs[c], s, e, st)
        aa = ""
        for i in range(0, len(nt) - 2, 3):
            from .seqio import translate_codon

            a, _ = translate_codon(nt[i:i + 3])
            if a == "*":
                break
            aa += a
        truth[(c, s, e, st)] = aa

    # Modelo con los CDS restantes.
    model = CodonModel()
    for c, s, e, st, t in kept:
        if t == "CDS":
            model.add_cds(_extract(contigs[c], s, e, st))
    regs_by_c: dict[str, list[tuple[int, int]]] = {}
    for c, seq in contigs.items():
        regs = padded_regions(
            len(seq), [(s, e) for cc, s, e, _, _ in kept if cc == c], min_len=150, pad=30
        )
        regs_by_c[c] = regs
        for s, e in regs:
            model.add_background(seq[s - 1:e])
    model.fit()

    ann_pref = {}
    for c, seq in contigs.items():
        arr = np.zeros(len(seq) + 2, dtype=np.int32)
        for cc, s, e, _, _ in kept:
            if cc == c:
                arr[s:e + 1] = 1
        ann_pref[c] = np.cumsum(arr)

    orfs = []
    for c, seq in contigs.items():
        for s, e in regs_by_c[c]:
            sub = seq[s - 1:e]
            for strand in ("+", "-"):
                oriented = sub if strand == "+" else revcomp(sub)
                for hit in find_orfs_oriented(oriented, min_aa=60):
                    i0, i1 = hit["start_idx"], hit["end_idx"]
                    gs, ge = (s + i0, s + i1 - 1) if strand == "+" else (e - i1 + 1, e - i0)
                    pref = ann_pref[c]
                    if int(pref[ge] - pref[gs - 1]) / (ge - gs + 1) >= 0.95:
                        continue
                    sc = model.score(oriented[i0:i1])
                    rb, _ = rbs_score(oriented, i0)
                    rec = {
                        "contig": c, "start": gs, "end": ge, "strand": strand,
                        "aa": hit["aa"], "aa_len": len(hit["aa"]),
                        "codon_llr": sc["codon_llr"], "frac_bad": sc["frac_bad"],
                        "rbs": rb, "start_w": start_weight(hit["start_codon"]),
                    }
                    rec["rank"] = rank_score(rec)
                    orfs.append(rec)

    # Recall
    recovered = 0
    details = []
    for (c, s, e, st), ref in truth.items():
        best = 0.0
        best_orf = None
        for o in orfs:
            if o["contig"] != c or o["strand"] != st:
                continue
            ov = max(0, min(o["end"], e) - max(o["start"], s) + 1)
            if ov / (e - s + 1) < 0.8:
                continue
            ident = peptide_match(o["aa"], ref)
            if ident > best:
                best, best_orf = ident, o
        ok = best >= 0.9
        recovered += int(ok)
        details.append({"contig": c, "start": s, "end": e, "identity": round(best, 3),
                        "recovered": ok, "rank": round(best_orf["rank"], 3) if best_orf else None})

    # Control negativo: barajar las regiones
    neg_orfs = []
    for c, seq in contigs.items():
        for s, e in regs_by_c[c]:
            sh = _shuffle_dinuc(seq[s - 1:e], rng)
            for strand in ("+", "-"):
                oriented = sh if strand == "+" else revcomp(sh)
                for hit in find_orfs_oriented(oriented, min_aa=60):
                    i0, i1 = hit["start_idx"], hit["end_idx"]
                    sc = model.score(oriented[i0:i1])
                    rb, _ = rbs_score(oriented, i0)
                    rec = {"aa_len": len(hit["aa"]), "codon_llr": sc["codon_llr"],
                           "frac_bad": sc["frac_bad"], "rbs": rb,
                           "start_w": start_weight(hit["start_codon"])}
                    rec["rank"] = rank_score(rec)
                    neg_orfs.append(rec)

    pos_llr = [o["codon_llr"] for o in orfs]
    neg_llr = [o["codon_llr"] for o in neg_orfs]
    neg_kb = sum(e - s + 1 for regs in regs_by_c.values() for s, e in regs) / 1000.0

    doc = {
        "genome": genome,
        "held_out": len(held),
        "recovered": recovered,
        "recall": round(recovered / len(held), 4) if held else None,
        "mean_identity": round(statistics.mean(d["identity"] for d in details), 4) if details else None,
        "orfs_in_regions": len(orfs),
        "llr_real_cds_hidden": {
            "mean": round(statistics.mean(pos_llr), 4) if pos_llr else None,
            "max": round(max(pos_llr), 4) if pos_llr else None,
        },
        "negative_control": {
            "orfs": len(neg_orfs),
            "per_kb": round(len(neg_orfs) / neg_kb, 3) if neg_kb else None,
            "mean_llr": round(statistics.mean(neg_llr), 4) if neg_llr else None,
            "frac_rank_ge_075": round(sum(1 for o in neg_orfs if o["rank"] >= 0.75) / max(1, len(neg_orfs)), 4),
            "frac_rank_ge_055": round(sum(1 for o in neg_orfs if o["rank"] >= 0.55) / max(1, len(neg_orfs)), 4),
        },
        "details": details,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{genome}_validation.json").write_text(json.dumps(doc, indent=1))
    return doc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--k", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    d = run(a.data_root, a.genome, a.k, a.seed, a.out)
    print(json.dumps({k: v for k, v in d.items() if k != "details"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
