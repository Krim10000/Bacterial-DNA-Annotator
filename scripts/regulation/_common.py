#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilidades comunes de la capa regulatoria (scripts/regulation).

Contrato:
  * determinista: orden estable, sin aleatoriedad, floats redondeados;
  * sin red;
  * solo escribe bajo `data/regulation/`.

Vocabulario de estados: se importa el módulo existente
`scripts/ghost/ghost/status.py` (no se duplica).
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import re
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # dna_poc/
DATA = ROOT / "data"
REG = DATA / "regulation"
AB = DATA / "annotations_ab_initio"
NCBI_GFF = DATA / "annotations"
RAW = DATA / "raw"
TMP = REG / ".tmp"


def _load_status():
    p = ROOT / "scripts" / "ghost" / "ghost" / "status.py"
    spec = importlib.util.spec_from_file_location("ghost_status", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ghost_status"] = mod          # requerido por dataclasses
    spec.loader.exec_module(mod)
    return mod


_st = _load_status()
STATUSES = _st.STATUSES
decide = _st.decide

# --------------------------------------------------------------------------- #
# secuencia / motivos
# --------------------------------------------------------------------------- #
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
_ALLOWED = {k: set(re.sub(r"[\[\]]", "", v)) for k, v in IUPAC.items()}
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def rx(motif: str) -> str:
    return "".join(IUPAC.get(c, c) for c in motif.upper())


def _seed(motif: str) -> tuple[str, int]:
    """Fragmento contiguo no degenerado más largo (ancla de búsqueda)."""
    best, cur, best_at, at = "", "", 0, 0
    for i, c in enumerate(motif.upper()):
        if c in "ACGT":
            cur += c
            if len(cur) > len(best):
                best, best_at = cur, i - len(cur) + 1
        else:
            cur = ""
    return best, best_at


def scan(motif: str, window: str) -> list[tuple[int, int, str]]:
    """Sitios del motivo en `window` como [(offset, mismatches, seq)].

    Búsqueda por semilla exacta + verificación IUPAC. No es un PWM: no
    puntúa afinidad, solo coincidencia. Se declara así en la evidencia.
    """
    motif = motif.upper()
    window = window.upper()
    seed, seed_at = _seed(motif)
    if not seed:
        return []
    pat = re.compile(rx(motif))
    out, pos = [], 0
    while True:
        j = window.find(seed, pos)
        if j < 0:
            break
        pos = j + 1
        s = j - seed_at
        if s < 0 or s + len(motif) > len(window):
            continue
        sub = window[s:s + len(motif)]
        if not pat.fullmatch(sub):
            continue
        mm = sum(1 for a, m in zip(sub, motif) if a not in _ALLOWED[m])
        out.append((s, mm, sub))
    return out


def best_site(motif: str, window: str):
    hits = scan(motif, window)
    if not hits:
        return None
    return min(hits, key=lambda h: (h[1], h[0]))


# --------------------------------------------------------------------------- #
# genes / genomas
# --------------------------------------------------------------------------- #
@dataclass
class Gene:
    contig: str
    start: int          # 1-based, inclusivo
    end: int            # 1-based, inclusivo
    strand: int         # +1 / -1
    gene_id: str
    gene_type: str = "CDS"
    locus_tag: str = ""
    gene_name: str = ""
    product: str = ""
    ec: str = ""
    acc: str = ""
    fam: str = ""        # familia UniProt (stem del entry name)
    partial: str = "00"
    rbs_motif: str = ""
    rbs_source: str = ""
    start_type: str = ""
    conf: str = ""
    source: str = ""

    @property
    def is_rna(self) -> bool:
        return self.gene_type != "CDS"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Genome:
    genome_id: str
    kind: str                       # ab_initio | gbff
    genes: list[Gene] = field(default_factory=list)
    contigs: dict[str, int] = field(default_factory=dict)
    seqs: dict[str, str] = field(default_factory=dict)
    source: dict = field(default_factory=dict)

    def cds(self) -> list[Gene]:
        return [g for g in self.genes if not g.is_rna]

    def rna(self) -> list[Gene]:
        return [g for g in self.genes if g.is_rna]

    def has_seq(self) -> bool:
        return bool(self.seqs)

    def sequence(self, contig: str) -> str:
        return self.seqs.get(contig, "")


def _tx_key(g: Gene):
    """Clave de orden de transcripción estable dentro de (contig, hebra)."""
    return -g.end if g.strand == -1 else g.start


def _overlaps(a: Gene, b: Gene) -> bool:
    return a.contig == b.contig and a.start <= b.end and b.start <= a.end


def transcription_units(genome: Genome) -> list[list[Gene]]:
    """Genes en orden de transcripción, cortando en cambio de contig/hebra y en ARN.

    Un gen de ARN corta la unidad SOLO si está intercalado (misma hebra y sin
    solapar ningún CDS). Los ARN anidados dentro de un CDS o solapantes son
    elementos internos del transcrito (p. ej. xylZ dentro de xylA) y no cortan:
    tratarlos como frontera fragmenta operones reales.
    """
    units: list[list[Gene]] = []
    contigs = sorted(genome.contigs) or sorted({g.contig for g in genome.genes})
    for contig in contigs:
        for strand in (1, -1):
            ordered = sorted((x for x in genome.genes
                              if x.contig == contig and x.strand == strand),
                             key=lambda x: (-x.end if strand == -1 else x.start))
            cds = [x for x in ordered if not x.is_rna]
            cur: list[Gene] = []
            for g in ordered:
                if g.is_rna:
                    # ncRNA/sRNA dentro de un operón es un miembro procesado del
                    # mismo transcrito (RegulonDB: xylZ en xylAB, fimR2 en
                    # fimAICDFGH, rbsZ en rbsDACBZKR). Solo tRNA/rRNA cortan.
                    if g.gene_type == "ncRNA":
                        continue
                    if cur:
                        units.append(cur)
                    cur = []
                else:
                    cur.append(g)
            if cur:
                units.append(cur)
    return units


def gap_between(up: Gene, down: Gene) -> int:
    """Distancia del 3' de `up` al 5' de `down` (negativo = solapamiento)."""
    if up.strand == 1:
        return down.start - up.end
    return up.start - down.end


def upstream_seq(seq: str, start: int, end: int, strand: int, n: int) -> str:
    """`n` nt río arriba del 5' del gen (orientación de transcripción)."""
    if not seq:
        return ""
    if strand == 1:
        a = max(0, start - 1 - n)
        return seq[a:start - 1]
    a = min(len(seq), end)
    return revcomp(seq[a:a + n])


def downstream_seq(seq: str, start: int, end: int, strand: int, n: int) -> str:
    if not seq:
        return ""
    if strand == 1:
        return seq[end:end + n]
    a = max(0, start - 1 - n)
    return revcomp(seq[a:start - 1])


RBS_PAT = re.compile(r"GGAGG|AGGAG|AGGA|GGAG")   # más largo primero


def detect_rbs(seq: str, gene: Gene) -> str:
    """RBS tipo Shine-Dalgarno en la ventana -20..-4 del inicio.

    Reimplementación simple (no es pyrodigal): devuelve el motivo encontrado
    o "" si no hay. Se marca `rbs_source=computed`.
    """
    if not seq:
        return ""
    w = upstream_seq(seq, gene.start, gene.end, gene.strand, 20)
    if not w:
        return ""
    best = ""
    for m in RBS_PAT.finditer(w):
        if len(m.group(0)) > len(best):
            best = m.group(0)
    return best


def rbs_strength(motif: str) -> int:
    """0 sin RBS, 1 débil (4-mer), 2 fuerte (5-6-mer Shine-Dalgarno)."""
    m = (motif or "").upper()
    if not m or m == "NONE":
        return 0
    joined = re.sub(r"[^ACGT]", "", m)
    if any(x in joined for x in ("GGAGG", "AGGAG", "AGGAGG")):
        return 2
    if "AGGA" in joined or "GGAG" in joined:
        return 1
    return 0


def has_strong_rbs(motif: str) -> bool:
    return rbs_strength(motif) >= 2


# --------------------------------------------------------------------------- #
# parsers
# --------------------------------------------------------------------------- #
def _attrs(s: str) -> dict[str, str]:
    out = {}
    for kv in s.strip().strip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_gff3(path: Path) -> tuple[list[Gene], dict[str, int]]:
    genes, contigs = [], {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    contigs[p[1]] = int(p[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = _attrs(f[8])
            rbs = a.get("rbs_motif", "")
            genes.append(Gene(
                contig=f[0], start=int(f[3]), end=int(f[4]),
                strand=1 if f[6] == "+" else -1,
                gene_id=a.get("ID", ""),
                partial=a.get("partial", ""),
                rbs_motif=rbs, rbs_source="pyrodigal" if rbs else "",
                start_type=a.get("start_type", ""), conf=a.get("conf", ""),
                source="pyrodigal"))
    return genes, contigs


def parse_functional(path: Path) -> dict[str, dict]:
    """functional.tsv (outfmt 6, 13 col): qseqid -> acc/product/gene."""
    out = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            st = p[12]
            parts = p[1].split("|")
            acc = parts[1] if p[1].startswith("sp|") and len(parts) > 1 else p[1]
            entry = parts[2] if len(parts) > 2 else ""
            fam = entry.rsplit("_", 1)[0] if "_" in entry else entry
            prod = st.split(" OS=")[0].strip()
            tok = prod.split(" ", 1)
            if tok[0].startswith("sp|"):
                prod = tok[1] if len(tok) > 1 else ""
            gn = re.search(r"\bGN=(\S+)", st)
            ecs = re.findall(r"\bEC[=:\s]+([0-9]+\.[0-9.]+n?|1\.[0-9.]+)", st)
            out[p[0]] = {"acc": acc, "fam": fam, "product": prod,
                         "gene_name": gn.group(1) if gn else "",
                         "ec": ecs[0] if ecs else "",
                         "bitscore": p[11], "evalue": p[10]}
    return out


_RNA_TYPES = {"tRNA": "tRNA", "rRNA": "rRNA", "ncRNA": "ncRNA",
              "RNase_P_RNA": "ncRNA", "SRP_RNA": "ncRNA", "tmRNA": "ncRNA"}


def parse_ncbi_rna(genome_id: str) -> list[Gene]:
    """ARN del GFF de NCBI (el GFF ab initio no los tiene). Trae hebra."""
    p = NCBI_GFF / f"{genome_id}.gff.gz"
    if not p.exists():
        return []
    out, seen = [], set()
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in _RNA_TYPES:
                continue
            a = _attrs(f[8])
            key = (f[0], f[3], f[4], f[6], _RNA_TYPES[f[2]])
            if key in seen:
                continue
            seen.add(key)
            prod = a.get("product", "")
            out.append(Gene(
                contig=f[0], start=int(f[3]), end=int(f[4]),
                strand=1 if f[6] == "+" else -1,
                gene_id=a.get("ID", key[0] + ":%s-%s" % (f[3], f[4])),
                gene_type=_RNA_TYPES[f[2]],
                locus_tag=a.get("locus_tag", ""),
                gene_name=prod or a.get("locus_tag", ""),
                product=prod, source="ncbi_gff"))
    return out


def load_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    if not path.exists():
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def load_genome(genome_id: str, with_seq: bool = False) -> Genome | None:
    """Carga un genoma del set ab initio (CDS pyrodigal + ARN NCBI + funcional)."""
    d = AB / genome_id
    if not d.is_dir():
        return None
    g = Genome(genome_id, "ab_initio",
               source={"gff": str(d / "annotation.gff3"),
                       "functional": str(d / "functional.tsv")})
    cds, contigs = parse_gff3(d / "annotation.gff3")
    func = parse_functional(d / "functional.tsv")
    for c in cds:
        hit = func.get(c.gene_id, {})
        c.acc = hit.get("acc", "")
        c.fam = hit.get("fam", "")
        c.product = hit.get("product", "")
        c.gene_name = hit.get("gene_name", "")
        c.ec = hit.get("ec", "")
    rna = parse_ncbi_rna(genome_id)
    for r in rna:
        contigs.setdefault(r.contig, r.end)
        contigs[r.contig] = max(contigs.get(r.contig, 0), r.end)
    for c in cds:
        contigs.setdefault(c.contig, c.end)
        contigs[c.contig] = max(contigs.get(c.contig, 0), c.end)
    g.genes = cds + rna
    g.contigs = contigs
    if with_seq:
        g.seqs = load_fasta(RAW / f"{genome_id}_genomic.fna.gz")
        # recorta a la longitud real del contig (el GFF puede exceder por ARN)
        g.contigs = {k: len(v) for k, v in g.seqs.items()} or contigs
    # RBS calculado cuando pyrodigal no lo dio
    for c in g.genes:
        if not c.is_rna and not c.rbs_motif and g.seqs:
            c.rbs_motif = detect_rbs(g.seqs.get(c.contig, ""), c)
            c.rbs_source = "computed" if c.rbs_motif else ""
    return g


def load_gbff(path: Path, genome_id: str, diamond: dict | None = None) -> Genome:
    """Carga un GenBank completo (p. ej. E. coli K-12) a la misma estructura."""
    from Bio import SeqIO
    rec = SeqIO.read(str(path), "genbank")
    contig = rec.id
    seq = str(rec.seq).upper()
    g = Genome(genome_id, "gbff",
               source={"gbff": str(path), "contigs": 1})
    contigs = {contig: len(seq)}
    genes = []
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"):
            continue
        q = f.qualifiers
        start = min(int(p.start) for p in f.location.parts) + 1
        end = max(int(p.end) for p in f.location.parts)
        strand = -1 if f.location.strand == -1 else 1
        lt = (q.get("locus_tag") or [""])[0]
        name = (q.get("gene") or [""])[0] or lt
        gid = (q.get("protein_id") or [""])[0] or lt or f"{contig}:{start}"
        partial = "00"
        if any(getattr(p, "start", 0) == 0 or getattr(p, "end", 0) == len(seq)
               for p in f.location.parts):
            partial = "01"
        if "pseudo" in q or "pseudogene" in q:
            partial = partial if partial != "00" else "10"
        gt = "CDS" if f.type == "CDS" else (
            "ncRNA" if f.type in ("ncRNA", "tmRNA") else f.type)
        genes.append(Gene(
            contig=contig, start=start, end=end, strand=strand, gene_id=gid,
            gene_type=gt, locus_tag=lt, gene_name=name,
            product=(q.get("product") or [""])[0],
            ec=(q.get("EC_number") or [""])[0],
            partial=partial, source="genbank"))
    g.genes, g.contigs, g.seqs = genes, contigs, {contig: seq}
    if diamond:
        for c in g.genes:
            hit = diamond.get(c.gene_id) or diamond.get(c.locus_tag) or {}
            c.acc = hit.get("acc", "")
            c.fam = hit.get("fam", "")
            if not c.product:
                c.product = hit.get("product", "")
            if not c.gene_name:
                c.gene_name = hit.get("gene_name", "")
            if not c.ec:
                c.ec = hit.get("ec", "")
    for c in g.genes:
        if not c.is_rna and not c.rbs_motif:
            c.rbs_motif = detect_rbs(seq, c)
            c.rbs_source = "computed" if c.rbs_motif else ""
    return g


# --------------------------------------------------------------------------- #
# evidencia funcional
# --------------------------------------------------------------------------- #
_STOP = {"protein", "subunit", "chain", "putative", "probable", "uncharacterized",
         "hypothetical", "family", "domain", "precursor", "the", "and", "for",
         "with", "type", "system", "component", "membrane", "transport",
         "binding", "permease", "atpase", "synthase", "synthetase", "protein.",
         "dehydrogenase", "reductase", "kinase", "transferase", "isomerase",
         "hydrolase", "lyase", "ligase", "oxidase", "epimerase", "mutase",
         "phosphatase", "decarboxylase", "aminotransferase", "processing",
         "assembly", "precursor.", "unknown", "duf"}


def tokens(product: str) -> set[str]:
    out = set()
    for t in re.findall(r"[A-Za-z]{4,}", (product or "").lower()):
        if t not in _STOP and not t.endswith("ase") or t in ("arabinose", "maltose",
                                                            "lactose", "xylose",
                                                            "ribose", "galactose",
                                                            "tryptophan", "histidine"):
            if t not in _STOP:
                out.add(t)
    return out


def gene_prefix(a: Gene, b: Gene) -> str:
    n1, n2 = (a.gene_name or ""), (b.gene_name or "")
    if not n1 or not n2:
        return ""
    k = 0
    for x, y in zip(n1, n2):
        if x != y:
            break
        k += 1
    pre = n1[:k]
    if len(pre) < 3:
        return ""
    if pre.lower() in _STOP:
        return ""
    if len(pre) == 3 and pre.lower().startswith("y"):
        return ""
    return pre


def _specificity(count: int, n: int) -> float:
    """Cuán informativo es un token compartido según lo raro que sea."""
    if count <= 3:
        return 1.0
    if count <= 10:
        return 0.8
    if count <= 30:
        return 0.5
    if count <= 100:
        return 0.2
    return 0.0


def coherent(a: Gene, b: Gene, stats: dict | None = None):
    """Coherencia funcional débil (texto de anotación) + especificidad.

    Devuelve (bool, motivo, especificidad 0-1). Un prefijo génico raro
    (`ara`, `mal`) es evidencia mucho más fuerte que una palabra común.
    """
    stats = stats or {}
    n = stats.get("n", 0)
    if a.ec and b.ec:
        x, y = a.ec.split("."), b.ec.split(".")
        # 3 niveles: 2 niveles mete cualquier transportador (7.5) con otro
        if len(x) >= 3 and len(y) >= 3 and x[:3] == y[:3]:
            key = ".".join(x[:3])
            return True, f"ec:{key}", _specificity(stats.get("ec", {}).get(key, 0), n)
    pre = gene_prefix(a, b)
    if pre:
        return True, f"prefix:{pre}", _specificity(stats.get("prefix", {}).get(pre, 0), n)
    shared = tokens(a.product) & tokens(b.product)
    if shared:
        t = sorted(shared)[0]
        return True, "kw:" + t, _specificity(stats.get("token", {}).get(t, 0), n)
    return False, "", 0.0


def cooccurrence_stats(genes) -> dict:
    """Frecuencia de prefijos génicos, tokens de producto y clase EC."""
    from collections import Counter
    pre: Counter = Counter()
    tok: Counter = Counter()
    ec: Counter = Counter()
    for g in genes:
        nm = (g.gene_name or "").lower()
        if len(nm) >= 3:
            pre[nm[:3]] += 1
        for t in tokens(g.product):
            tok[t] += 1
        if g.ec:
            ec[".".join(g.ec.split(".")[:3])] += 1
    return {"prefix": pre, "token": tok, "ec": ec, "n": len(genes)}


# --------------------------------------------------------------------------- #
# salida / determinismo
# --------------------------------------------------------------------------- #
def _clean(v) -> str:
    s = "" if v is None else str(v)
    return s.replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_clean(r.get(c, "")) for c in cols) + "\n")


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]


def dump_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=True,
                               ensure_ascii=False, default=str) + "\n")


def sha256_files(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for p in sorted(paths):
        if not p.exists():
            continue
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()


def fnum(x: float, nd: int = 4) -> str:
    return f"{x:.{nd}f}"


def status_row(row: dict, decision) -> dict:
    """Inyecta status/evidence/attempted/confidence (vocabulario de 09_ghost)."""
    row["status"] = decision.status
    row["evidence"] = decision.evidence
    row["attempted"] = decision.attempted
    row["confidence"] = decision.confidence or ""
    row["validated_by"] = decision.validated_by
    row["reason"] = decision.reason
    return row


def load_genes_tsv(genome_id: str, out_dir: Path | None = None) -> Genome | None:
    """Reconstruye el genoma desde `genes.tsv` (interfaz estable L1->L2)."""
    d = (out_dir or (REG / genome_id))
    gp = d / "genes.tsv"
    mp = d / "genes_meta.json"
    if not mp.exists():
        mp = d / "meta.json"          # compatibilidad con corridas antiguas
    if not gp.exists():
        return None
    meta = json.loads(mp.read_text()) if mp.exists() else {}
    g = Genome(genome_id, meta.get("kind", ""),
               contigs=dict(meta.get("contigs", {})), source={"meta": str(mp)})
    for r in read_tsv(gp):
        g.genes.append(Gene(
            contig=r["contig"], start=int(r["start"]), end=int(r["end"]),
            strand=int(r["strand"]), gene_id=r["gene_id"],
            gene_type=r.get("gene_type") or "CDS",
            locus_tag=r.get("locus_tag", ""), gene_name=r.get("gene_name", ""),
            product=r.get("product", ""), ec=r.get("ec", ""),
            acc=r.get("acc", ""), fam=r.get("fam", ""),
            partial=r.get("partial", ""),
            rbs_motif=r.get("rbs_motif", ""), rbs_source=r.get("rbs_source", ""),
            start_type=r.get("start_type", ""), conf=r.get("conf", ""),
            source=r.get("source", "")))
    return g


def load_input(genome_id: str) -> Genome | None:
    """genes.tsv si existe; si no, parseo directo del bundle ab initio."""
    return load_genes_tsv(genome_id) or load_genome(genome_id)


def load_sequence(genome_id: str) -> dict[str, str]:
    """Secuencia del genoma según meta.json (gbff o raw fna.gz)."""
    mp = REG / genome_id / "genes_meta.json"
    if not mp.exists():
        mp = REG / genome_id / "meta.json"
    if not mp.exists():
        return {}
    meta = json.loads(mp.read_text())
    src = meta.get("seq_source") or {}
    kind, path = src.get("kind"), src.get("path")
    if not kind or not path or not Path(path).exists():
        return {}
    if kind == "raw":
        return load_fasta(Path(path))
    if kind == "gbff":
        from Bio import SeqIO
        out = {}
        for rec in SeqIO.parse(str(path), "genbank"):
            out[rec.id] = str(rec.seq).upper()
        return out
    return {}
