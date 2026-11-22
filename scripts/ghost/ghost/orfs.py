"""Regiones no anotadas y enumeración de ORFs con ambigüedad.

Concepto
--------
1. Las features anotadas (CDS/tRNA/rRNA/...) se fusionan y se calcula su
   complemento: las regiones intergénicas.
2. En cada región se enumeran ORFs en los 6 marcos (3 por hebra). Un ORF es
   el tramo sin codón de parada que va desde el primer inicio en fase hasta
   la parada siguiente.
3. Si un codón es ambiguo (IUPAC/N), se resuelve a un conjunto de
   aminoácidos. Si el número de combinaciones es pequeño se emiten varias
   cadenas; si es grande se emite una cadena con X en las posiciones
   ambiguas (y se marca la bandera).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Iterable

from .seqio import (
    ALT_START_CODONS,
    START_CODONS,
    revcomp,
    translate_codon,
)


# --------------------------------------------------------------------------
# Intervalos
# --------------------------------------------------------------------------

def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Fusiona intervalos 1-based inclusivos solapados o adyacentes."""
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if not out or s > out[-1][1] + 1:
            out.append([s, e])
        else:
            out[-1][1] = max(out[-1][1], e)
    return [(s, e) for s, e in out]


def intergenic_regions(
    contig_len: int,
    intervals: Iterable[tuple[int, int]],
    min_len: int = 150,
) -> list[tuple[int, int]]:
    """Huecos 1-based entre features; incluye extremos del contig."""
    merged = merge_intervals(intervals)
    regs: list[tuple[int, int]] = []
    cur = 1
    for s, e in merged:
        if s - cur >= min_len:
            regs.append((cur, s - 1))
        cur = max(cur, e + 1)
    if contig_len - cur + 1 >= min_len:
        regs.append((cur, contig_len))
    return regs


# --------------------------------------------------------------------------
# ORFs
# --------------------------------------------------------------------------

@dataclass
class ORF:
    region_id: str
    contig: str
    strand: str
    frame: int
    start: int          # 1-based, genómico, inclusivo
    end: int            # 1-based, genómico, inclusivo
    aa: str             # cadena consenso (X en ambigüedad)
    variants: list[str] = field(default_factory=list)
    ambiguous_positions: list[tuple[int, list[str]]] = field(default_factory=list)
    start_codon: str = ""
    stop_codon: str = ""
    n_ambiguous: int = 0
    ambig_stop: bool = False
    aas: list[str] = field(default_factory=list)

    @property
    def aa_len(self) -> int:
        return len(self.aa)

    @property
    def nt_len(self) -> int:
        return self.end - self.start + 1


def _build_variants(
    aas: list[str],
    poss: list[frozenset[str]],
    max_variants: int,
) -> tuple[str, list[str], list[tuple[int, list[str]]], int]:
    """Construye cadena consenso y variantes a partir de conjuntos posibles."""
    consensus = "".join(aas)
    amb = [(i, sorted(p)) for i, p in enumerate(poss) if len(p) > 1]
    n_amb = sum(1 for p in poss if len(p) > 1)
    if not amb:
        return consensus, [], [], 0

    # Producto cartesiano acotado.
    total = 1
    for _, options in amb:
        total *= len(options)
        if total > max_variants:
            break

    if total <= max_variants:
        variants = []
        for combo in product(*[options for _, options in amb]):
            seq = list(consensus)
            for (idx, _), aa in zip(amb, combo):
                seq[idx] = aa
            variants.append("".join(seq))
        return consensus, variants, amb, n_amb
    return consensus, [], amb, n_amb


def _orfs_in_frame(
    seq: str,
    frame: int,
    *,
    min_aa: int,
    starts: frozenset[str],
    max_variants: int,
) -> list[dict]:
    """ORFs de un único marco de lectura (frame = 0, 1 o 2)."""
    codons: list[tuple[int, str, str, frozenset[str]]] = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        aa, poss = translate_codon(codon)
        codons.append((i, codon, aa, poss))

    results: list[dict] = []
    # Segmentos delimitados por paradas. Una parada ambigua ("*" entre las
    # posibilidades, p. ej. TGR = TGA/* o TGG/W) también corta, pero se marca.
    seg_start = 0
    sentinel = (None, "", "*", frozenset({"*"}))
    for k, (idx, codon, aa, poss) in enumerate(codons + [sentinel]):
        is_stop = "*" in poss
        if not is_stop and k < len(codons):
            continue
        segment = codons[seg_start:k]
        seg_start = k + 1
        if not segment:
            continue
        # primer inicio en fase del segmento
        start_k = None
        for m, (i2, c2, a2, p2) in enumerate(segment):
            if c2 in starts:
                start_k = m
                break
        if start_k is None:
            continue
        sub = segment[start_k:]
        if len(sub) < min_aa:
            continue
        stop_codon = codon.upper() if is_stop else ""
        stop_poss = poss if is_stop else frozenset()
        aas = [a for _, _, a, _ in sub]
        # Una parada ambigua no aporta aminoácido a las variantes.
        poss_sub = [frozenset(p - {"*"}) for _, _, _, p in sub]
        consensus, variants, amb, n_amb = _build_variants(aas, poss_sub, max_variants)
        first = sub[0]
        last = sub[-1]
        results.append(
            dict(
                start_idx=first[0],
                end_idx=last[0] + 3,          # exclusivo
                frame=frame,
                aa=consensus,
                variants=variants,
                ambiguous_positions=amb,
                start_codon=first[1],
                stop_codon=stop_codon,
                ambig_stop=is_stop and len(stop_poss) > 1,
                partial=(codon == ""),
                n_ambiguous=n_amb,
            )
        )
    return results


def find_orfs_oriented(
    seq: str,
    *,
    min_aa: int = 60,
    starts: frozenset[str] = START_CODONS,
    max_variants: int = 8,
) -> list[dict]:
    """Encuentra ORFs en los tres marcos de una secuencia orientada 5'->3'."""
    out: list[dict] = []
    for frame in range(3):
        out.extend(
            _orfs_in_frame(
                seq, frame, min_aa=min_aa, starts=starts, max_variants=max_variants
            )
        )
    return out


def enumerate_regions(
    contig: str,
    seq: str,
    regions: list[tuple[int, int]],
    *,
    region_prefix: str,
    min_aa: int = 60,
    include_revcomp: bool = True,
    starts: frozenset[str] = START_CODONS,
    max_variants: int = 8,
) -> list[ORF]:
    """Enumera ORFs en todas las regiones de un contig (ambas hebras)."""
    out: list[ORF] = []
    for ridx, (rs, re_) in enumerate(regions):
        region_id = f"{region_prefix}_r{ridx}"
        sub = seq[rs - 1:re_]
        strands = ("+", "-") if include_revcomp else ("+",)
        for strand in strands:
            oriented = sub if strand == "+" else revcomp(sub)
            for hit in find_orfs_oriented(
                oriented, min_aa=min_aa, starts=starts, max_variants=max_variants
            ):
                i0, i1 = hit["start_idx"], hit["end_idx"]
                if strand == "+":
                    gstart = rs + i0
                    gend = rs + i1 - 1
                else:
                    gstart = re_ - i1 + 1
                    gend = re_ - i0
                out.append(
                    ORF(
                        region_id=region_id,
                        contig=contig,
                        strand=strand,
                        frame=hit.get("frame", i0 % 3),
                        start=gstart,
                        end=gend,
                        aa=hit["aa"],
                        variants=hit["variants"],
                        ambiguous_positions=hit["ambiguous_positions"],
                        start_codon=hit["start_codon"],
                        stop_codon=hit["stop_codon"],
                        n_ambiguous=hit["n_ambiguous"],
                        ambig_stop=hit.get("ambig_stop", False),
                        aas=list(hit["aa"]),
                    )
                )
    return out


def padded_regions(
    contig_len: int,
    intervals: Iterable[tuple[int, int]],
    *,
    min_len: int = 150,
    pad: int = 30,
) -> list[tuple[int, int]]:
    """Regiones intergénicas extendidas ``pad`` nt hacia los genes vecinos.

    Los genes bacterianos suelen solaparse unos pocos nucleótidos; sin pad, el
    codón de inicio o de parada de un gen no anotado puede quedar fuera de la
    región y el ORF no se genera. El pad se recorta a los límites del contig.
    """
    merged = merge_intervals(intervals)
    regs: list[tuple[int, int]] = []
    cur = 1
    for s, e in merged:
        gap_len = s - cur
        if gap_len > 0 or pad > 0:
            s2 = max(1, cur - pad)
            e2 = min(contig_len, s - 1 + pad)
            if e2 - s2 + 1 >= min_len and e2 >= s2:
                regs.append((s2, e2))
        cur = max(cur, e + 1)
    s2 = max(1, cur - pad)
    e2 = contig_len
    if e2 - s2 + 1 >= min_len:
        regs.append((s2, e2))
    return regs
