"""Puntuaciones ligeras (CPU, numpy opcional) para priorizar ORFs.

Capas
-----
1. Potencial codificante: modelo de Markov de orden 1 sobre codones entrenado
   con los CDS anotados del propio genoma y comparado contra el fondo
   intergénico (log-odds). Detecta estructura de codones típica de proteína.
2. Inicio de traducción: motivo Shine-Dalgarno (RBS) y tipo de codón de inicio.
3. Composición proteica: GRAVY, fracción hidrofóbica, desorden estimado.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

from .seqio import CODON_TABLE, STOP_CODONS, gc

# Kyte-Doolittle
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
DISORDER = set("PESQKR")
HYDROPHOBIC = set("AILMFWV")


class CodonModel:
    """Modelo log-odds de codones (orden 1 con backoff a orden 0)."""

    def __init__(self, alpha: float = 20.0, pseudocount: float = 1.0):
        self.alpha = alpha
        self.pc = pseudocount
        self.cds0: Counter = Counter()
        self.bg0: Counter = Counter()
        self.cds1: dict[str, Counter] = defaultdict(Counter)
        self.bg1: dict[str, Counter] = defaultdict(Counter)
        self.llr0: dict[str, float] = {}
        self.llr1: dict[str, dict[str, float]] = {}
        self.fitted = False

    # -- entrenamiento -----------------------------------------------------
    def add_cds(self, nt: str) -> None:
        nt = nt.upper()
        prev = None
        for i in range(0, len(nt) - 2, 3):
            codon = nt[i:i + 3]
            if set(codon) - set("ACGT"):
                prev = None
                continue
            aa = CODON_TABLE[codon]
            if aa == "*":
                prev = None
                continue
            self.cds0[codon] += 1
            if prev is not None:
                self.cds1[prev][codon] += 1
            prev = codon

    def add_background(self, nt: str) -> None:
        nt = nt.upper()
        prev = None
        for i in range(0, len(nt) - 2, 3):
            codon = nt[i:i + 3]
            if set(codon) - set("ACGT"):
                prev = None
                continue
            self.bg0[codon] += 1
            if prev is not None:
                self.bg1[prev][codon] += 1
            prev = codon

    def fit(self) -> "CodonModel":
        codons = list(CODON_TABLE)
        total_c = sum(self.cds0.values())
        total_b = sum(self.bg0.values())
        for c in codons:
            p_c = (self.cds0[c] + self.pc) / (total_c + self.pc * 64)
            p_b = (self.bg0[c] + self.pc) / (total_b + self.pc * 64)
            self.llr0[c] = math.log2(p_c / p_b)
        for prev in set(self.cds1) | set(self.bg1):
            tc = sum(self.cds1[prev].values())
            tb = sum(self.bg1[prev].values())
            row: dict[str, float] = {}
            for c in codons:
                p0_c = (self.cds0[c] + self.pc) / (total_c + self.pc * 64)
                p0_b = (self.bg0[c] + self.pc) / (total_b + self.pc * 64)
                p_c = (self.cds1[prev][c] + self.alpha * p0_c) / (tc + self.alpha)
                p_b = (self.bg1[prev][c] + self.alpha * p0_b) / (tb + self.alpha)
                row[c] = math.log2(p_c / p_b)
            self.llr1[prev] = row
        self.fitted = True
        return self

    # -- puntuación --------------------------------------------------------
    def score(self, nt: str) -> dict[str, float]:
        """Devuelve LLR medio, fracción de codones desfavorables y GC3."""
        nt = nt.upper()
        if not self.fitted:
            raise RuntimeError("CodonModel no entrenado")
        vals: list[float] = []
        gc3 = 0
        n = 0
        prev = None
        for i in range(0, len(nt) - 2, 3):
            codon = nt[i:i + 3]
            if set(codon) - set("ACGT"):
                prev = None
                continue
            aa = CODON_TABLE[codon]
            if aa == "*":
                prev = None
                continue
            v = self.llr1.get(prev, {}).get(codon, self.llr0.get(codon, 0.0))
            vals.append(v)
            if codon[2] in "GC":
                gc3 += 1
            n += 1
            prev = codon
        if not vals:
            return {"codon_llr": 0.0, "frac_bad": 1.0, "gc3": 0.0, "n_codons": 0}
        mean = sum(vals) / len(vals)
        bad = sum(1 for v in vals if v < -1.0) / len(vals)
        return {
            "codon_llr": mean,
            "frac_bad": bad,
            "gc3": gc3 / n if n else 0.0,
            "n_codons": n,
        }


# --------------------------------------------------------------------------
# RBS
# --------------------------------------------------------------------------

_RBS = {
    "AGGAGG": 1.00,
    "AAGGAG": 0.95,
    "GGAGG": 0.85,
    "AGGAG": 0.85,
    "AGGA": 0.70,
    "GGAG": 0.60,
    "GAGG": 0.55,
    "AGGT": 0.45,
}


def rbs_score(oriented: str, start_idx: int) -> tuple[float, str]:
    """Puntúa el Shine-Dalgarno en los 20 nt río arriba del inicio."""
    lo = max(0, start_idx - 20)
    up = oriented[lo:start_idx].upper()
    best = 0.0
    best_motif = ""
    for motif, w in _RBS.items():
        pos = up.rfind(motif)
        if pos < 0:
            continue
        spacer = len(up) - (pos + len(motif))
        if not (3 <= spacer <= 12):
            continue
        # Óptimo ~7 nt de separación.
        factor = max(0.0, 1.0 - abs(spacer - 7) / 12.0)
        val = w * factor
        if val > best:
            best, best_motif = val, f"{motif}@{spacer}nt"
    return best, best_motif


# --------------------------------------------------------------------------
# Composición
# --------------------------------------------------------------------------

def protein_features(aa: str) -> dict[str, float]:
    clean = [a for a in aa if a in KD]
    n = len(clean)
    if not n:
        return {"gravy": 0.0, "hydro": 0.0, "disorder": 0.0, "aromatic": 0.0}
    return {
        "gravy": sum(KD[a] for a in clean) / n,
        "hydro": sum(1 for a in clean if a in HYDROPHOBIC) / n,
        "disorder": sum(1 for a in clean if a in DISORDER) / n,
        "aromatic": sum(1 for a in clean if a in "FWY") / n,
    }


def start_weight(codon: str) -> float:
    return {"ATG": 1.0, "GTG": 0.75, "TTG": 0.55}.get(codon, 0.35)


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------

def rank_score(c: dict) -> float:
    """Puntaje 0-1 combinando señales ligeras (sin estructura)."""
    llr = c.get("codon_llr", 0.0)
    # LLR típico: CDS real > 0.5; intergénico ~ -0.3. Escala suave.
    s_coding = 1.0 / (1.0 + math.exp(-(llr - 0.15) * 3.0))
    s_len = min(1.0, math.log10(max(c.get("aa_len", 1), 10) / 60.0))
    s_rbs = c.get("rbs", 0.0)
    s_start = c.get("start_w", 0.5)
    s_pure = 1.0 - min(1.0, c.get("frac_bad", 0.5))
    return (
        0.42 * s_coding
        + 0.12 * s_len
        + 0.16 * s_rbs
        + 0.12 * s_start
        + 0.18 * s_pure
    )
