"""E/S de secuencias y traducción con notación ambigua.

Sin dependencias externas (solo biblioteca estándar). Soporta FASTA/FASTQ
básico, gzip, GFF3 y traducción de codones IUPAC.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator


# --------------------------------------------------------------------------
# Lectura
# --------------------------------------------------------------------------

def open_text(path: str | Path):
    """Abre texto plano o gzip de forma transparente."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Itera (nombre, secuencia) de un FASTA (posiblemente .gz)."""
    name = None
    buf: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def fasta_index(path: str | Path) -> dict[str, int]:
    """Devuelve {contig: longitud} sin cargar las secuencias completas."""
    out: dict[str, int] = {}
    name = None
    n = 0
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    out[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        out[name] = n
    return out


# --------------------------------------------------------------------------
# Complemento reverso y tablas
# --------------------------------------------------------------------------

_COMP = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


IUPAC: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# Tabla estándar (código genético 11) en orden TCAG.
_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODON_TABLE: dict[str, str] = {}
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            CODON_TABLE[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]

STOP_CODONS = frozenset(c for c, a in CODON_TABLE.items() if a == "*")
START_CODONS = frozenset({"ATG", "GTG", "TTG"})
# Inicios alternativos documentados en bacterias (se pueden activar por CLI).
ALT_START_CODONS = frozenset({"ATT", "CTG", "ATC", "ATA"})


def translate_codon(codon: str) -> tuple[str, frozenset[str]]:
    """Traduce un codón, resolviendo ambigüedad IUPAC.

    Devuelve (aa, posibles). Si el codón es inequívoco, ``posibles`` tiene un
    solo elemento. Si es ambiguo, ``aa`` es "X" salvo que todas las
    resoluciones den el mismo aminoácido (p. ej. GCN -> A).
    """
    codon = codon.upper()
    if len(codon) != 3:
        return "X", frozenset({"X"})
    if set(codon) <= {"A", "C", "G", "T"}:
        aa = CODON_TABLE[codon]
        return aa, frozenset({aa})

    options: list[str] = []
    for b in codon:
        options.append(IUPAC.get(b, "ACGT"))
    # Producto cartesiano acotado: peor caso NNN = 64.
    aas: set[str] = set()
    for a in options[0]:
        for b in options[1]:
            for c in options[2]:
                aas.add(CODON_TABLE[a + b + c])
    if len(aas) == 1:
        aa = next(iter(aas))
        return aa, frozenset({aa})
    return "X", frozenset(aas)


def translate(seq: str, to_stop: bool = False) -> str:
    """Traducción simple de una secuencia en marco 0."""
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa, _ = translate_codon(seq[i:i + 3])
        if aa == "*":
            if to_stop:
                break
            out.append("*")
        else:
            out.append(aa)
    return "".join(out)


def gc(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


# --------------------------------------------------------------------------
# GFF3
# --------------------------------------------------------------------------

def read_gff_features(
    path: str | Path,
    types: frozenset[str] = frozenset({"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"}),
) -> tuple[dict[str, int], list[tuple[str, int, int, str, str]]]:
    """Lee un GFF3.

    Devuelve (longitudes de contig, features) donde cada feature es
    (contig, inicio 1-based, fin 1-based, strand, tipo).
    Las longitudes salen de las directivas ##sequence-region.
    """
    lengths: dict[str, int] = {}
    feats: list[tuple[str, int, int, str, str]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        lengths[parts[1]] = int(parts[3])
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != 9:
                continue
            if cols[2] not in types:
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            feats.append((cols[0], start, end, cols[6], cols[2]))
    return lengths, feats
