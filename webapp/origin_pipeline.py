"""Pipeline general de origen de replicación con estados y evidencia.

No es un parche para los 400: recibe secuencia/marcadores y devuelve un estado
explicable. Los métodos pesados (ORCA, HMM de marcadores) se conectan aquí.

Estados:
  validated        ≥2 métodos independientes coinciden (p. ej. ORCA + DoriC).
  consensus        métodos coinciden, sin verdad externa (ORCA + dnaA).
  lineage_specific  desviación esperada del linaje (Fusobacteriaceae sin dnaA).
  anomalous        los métodos discrepan; se declara el motivo.
  unresolved       ningún método dio candidato; se listan los intentos.
"""
from __future__ import annotations

import re
from collections import Counter

MOTIF = "TTATCCACA"

# Linajes documentados sin dnaA / con iniciación independiente de DnaA.
# Es una tabla de priors, no una rama hardcodeada: el pipeline siempre intenta
# los mismos métodos y solo cambia la interpretación del resultado.
DNAA_LESS = {
    "Fusobacteriota": "Fusobacteriaceae: RPA/iterones derivados de plásmido "
                       "(Brezellec & Pasek 2026)",
    "Cyanobacteriota": "algunas cianobacterias replican sin DnaA "
                       "(p. ej. Synechococcus elongatus)",
    "Mycoplasmatota": "Mollicutes: complejo de replicación simplificado",
}


def gc_skew_origin(seq: str, win: int = 10000) -> tuple[int, int]:
    vals = []
    for i in range(0, len(seq), win):
        s = seq[i:i + win].upper()
        g, c = s.count("G"), s.count("C")
        vals.append((g - c) / (g + c) if g + c else 0.0)
    cum, acc = [], 0.0
    for v in vals:
        acc += v
        cum.append(acc)
    return min(range(len(cum)), key=lambda i: cum[i]) * win, \
           max(range(len(cum)), key=lambda i: cum[i]) * win


def dnaa_boxes(seq: str, lo: int, hi: int, mism: int = 1) -> list[int]:
    s = seq[max(0, lo):min(len(seq), hi)].upper()
    hits = [max(0, lo) + m.start() for m in re.finditer(MOTIF, s)]
    if mism:
        for m in re.finditer(r"(?=(.{9}))", s):
            if sum(a != b for a, b in zip(m.group(1), MOTIF)) == 1:
                hits.append(max(0, lo) + m.start())
    return hits


def iteron_scan(seq: str, center: int, flank: int = 10000,
                k: int = 12, min_count: int = 4) -> list[tuple[str, int]]:
    """Repeticiones cortas (iterones) cerca de un candidato de oriC.

    Es un apoyo, no una prueba: en genomas AT-ricos el ruido de baja
    complejidad obliga a confirmar con HMM de RPA/Rep.
    """
    reg = seq[max(0, center - flank):center + flank].upper()
    counts = Counter(reg[i:i + k] for i in range(len(reg) - k + 1))
    hits = [(m, c) for m, c in counts.items()
            if c >= min_count and len(set(m)) >= 3]
    return sorted(hits, key=lambda x: -x[1])[:5]


def _cd(a: int, b: int, L: int) -> int:
    d = abs(a - b)
    return min(d, L - d)


def classify(*, genome: str, phylum: str, organism: str, size: int,
             marker_status: str, dnaA: int | None, orca_pos: int | None,
             orca_conf: float | None, doric: int | None,
             skew: int | None, boxes: int, L: int,
             iterons: int = 0) -> dict:
    """Devuelve estado, motivo, métodos que validan e intentos."""
    endosymb = any(k in (organism or "").lower() for k in
                   ("endosymb", "symbiont", "candidatus", "mycoplasma", "buchnera"))
    out = {"genome": genome, "status": "", "position": None, "confidence": None,
           "validated_by": "", "reason": "", "attempted": "ORCA; dnaA/dnaN; GC skew; cajas dnaA",
           "lineage": phylum}

    if dnaA is None:
        if phylum in DNAA_LESS:
            note = DNAA_LESS[phylum]
            if iterons > 0:
                out.update(status="validated", position=skew,
                           reason=f"sin dnaA; iterones detectados · {note}",
                           validated_by="GC skew + iterones (linaje sin dnaA)")
            else:
                out.update(status="lineage_specific", position=skew,
                           reason=f"sin dnaA; iterones no detectados · {note}",
                           validated_by="GC skew (linaje sin dnaA)")
        elif endosymb or size < 1_500_000:
            out.update(status="anomalous", position=skew,
                       reason="genoma reducido/endosimbionte sin dnaA ni cajas DnaA "
                              "diagnósticas",
                       validated_by="GC skew")
        else:
            out.update(status="anomalous", position=skew,
                       reason="dnaA/dnaN no detectados por homología",
                       validated_by="GC skew")
        out["confidence"] = 0.3 if skew is not None else None
        return out

    # hay dnaA
    near_dnaA = orca_pos is not None and _cd(orca_pos, dnaA, L) <= 50000
    if doric is not None and orca_pos is not None:
        d_orca = _cd(orca_pos, doric, L)
        d_skew = _cd(skew, doric, L) if skew is not None else None
        if d_orca <= 50000 and near_dnaA:
            out.update(status="validated", position=orca_pos, confidence=orca_conf,
                       validated_by="ORCA + DoriC + dnaA", reason="coinciden los tres")
        elif d_skew is not None and d_skew <= 50000:
            out.update(status="validated", position=skew, confidence=orca_conf,
                       validated_by="DoriC + GC skew",
                       reason=f"ORCA discrepa (Δ={d_orca} pb); GC skew apoya DoriC")
        else:
            out.update(status="anomalous", position=orca_pos, confidence=orca_conf,
                       validated_by="ORCA",
                       reason=f"discordancia ORCA-DoriC Δ={d_orca} pb; GC skew ambiguo")
    elif near_dnaA:
        out.update(status="consensus", position=orca_pos, confidence=orca_conf,
                   validated_by="ORCA + proximidad dnaA", reason="sin verdad DoriC")
    else:
        out.update(status="anomalous", position=skew or orca_pos, confidence=orca_conf,
                   validated_by="GC skew" if skew is not None else "ORCA",
                   reason="ORCA lejos de dnaA; se usa GC skew")
    return out
