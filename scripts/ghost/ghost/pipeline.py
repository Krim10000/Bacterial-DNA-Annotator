"""Pipeline de descubrimiento de genes no anotados.

Flujo por genoma:
    FASTA + GFF3
      -> regiones intergénicas
      -> ORFs en 6 marcos (con ambigüedad)
      -> scoring ligero (codones, RBS, composición)
      -> [opcional] DIAMOND vs Swiss-Prot / proteomas
      -> TSV + JSON priorizados
"""
from __future__ import annotations

import json
import math
import re
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .orfs import ORF, enumerate_regions, padded_regions
from .scoring import (
    CodonModel,
    protein_features,
    rank_score,
    rbs_score,
    start_weight,
)
from .seqio import (
    ALT_START_CODONS,
    START_CODONS,
    gc,
    read_fasta,
    read_gff_features,
    revcomp,
)


@dataclass
class GenomeResult:
    genome: str
    n_contigs: int = 0
    genome_bp: int = 0
    n_features: int = 0
    n_regions: int = 0
    region_bp: int = 0
    n_orfs: int = 0
    candidates: list[dict] = field(default_factory=list)
    model: dict = field(default_factory=dict)
    contigs: set[str] = field(default_factory=set)

    def stats(self) -> dict:
        return {
            "genome": self.genome,
            "n_contigs": self.n_contigs,
            "genome_bp": self.genome_bp,
            "n_features": self.n_features,
            "n_regions": self.n_regions,
            "region_bp": self.region_bp,
            "region_frac": round(self.region_bp / self.genome_bp, 4) if self.genome_bp else 0,
            "n_orfs": self.n_orfs,
            "n_candidates": len(self.candidates),
            "model": self.model,
        }


def _find_fasta(data_root: Path, genome: str) -> Path | None:
    for sub in ("raw", "raw_excluded"):
        raw = data_root / sub
        if not raw.is_dir():
            continue
        for pat in (f"{genome}_genomic.fna.gz", f"{genome}*.fna.gz", f"{genome}*.fna"):
            hits = sorted(raw.glob(pat))
            if hits:
                return hits[0]
    return None


def _extract(seq: str, start: int, end: int, strand: str) -> str:
    sub = seq[start - 1:end]
    return sub if strand == "+" else revcomp(sub)


def _nearest_genes(
    feats: list[tuple[str, int, int, str, str]],
    contig: str,
) -> tuple[list[tuple[int, int]], list[str]]:
    ordered = sorted(
        (s, e, t) for c, s, e, _, t in feats if c == contig
    )
    return [(s, e) for s, e, _ in ordered], [t for _, _, t in ordered]


def build_candidates(
    data_root: str | Path,
    genome: str,
    *,
    min_region_len: int = 150,
    min_aa: int = 60,
    max_variants: int = 8,
    alt_starts: bool = False,
    max_regions: int | None = None,
    pad: int = 30,
    exclude_partial: bool = False,
) -> GenomeResult:
    data_root = Path(data_root)
    fasta = _find_fasta(data_root, genome)
    if fasta is None:
        raise FileNotFoundError(f"FASTA no encontrado para {genome} en {data_root/'raw'}")
    gff = data_root / "annotations_ab_initio" / genome / "annotation.gff3"
    if not gff.exists():
        raise FileNotFoundError(f"GFF no encontrado: {gff}")

    contigs = dict(read_fasta(fasta))
    lengths, feats = read_gff_features(gff)
    starts = START_CODONS | ALT_START_CODONS if alt_starts else START_CODONS

    res = GenomeResult(genome=genome)
    res.n_contigs = len(contigs)
    res.genome_bp = sum(len(s) for s in contigs.values())
    res.n_features = len(feats)
    res.contigs = set(contigs)

    # --- regiones intergénicas ------------------------------------------
    regions_by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, seq in contigs.items():
        iv = [(s, e) for c, s, e, _, _ in feats if c == contig]
        regs = padded_regions(len(seq), iv, min_len=min_region_len, pad=pad)
        if max_regions is not None:
            regs = regs[:max_regions]
        regions_by_contig[contig] = regs
        res.n_regions += len(regs)
        res.region_bp += sum(e - s + 1 for s, e in regs)

    # --- modelo de codones ----------------------------------------------
    model = CodonModel()
    for contig, s, e, strand, typ in feats:
        if typ != "CDS":
            continue
        model.add_cds(_extract(contigs[contig], s, e, strand))
    for contig, regs in regions_by_contig.items():
        for s, e in regs:
            model.add_background(contigs[contig][s - 1:e])
    model.fit()
    res.model = {
        "cds_codons": int(sum(model.cds0.values())),
        "bg_codons": int(sum(model.bg0.values())),
        "bg_fraction": round(
            sum(model.bg0.values()) / max(1, sum(model.cds0.values()) + sum(model.bg0.values())), 4
        ),
    }

    # --- máscara de lo ya anotado ---------------------------------------
    # Con pad, las regiones invaden los genes vecinos. Un ORF que cae casi
    # entero dentro de una feature anotada es un fragmento de un gen conocido.
    ann_pref: dict[str, np.ndarray] = {}
    for contig, seq in contigs.items():
        arr = np.zeros(len(seq) + 2, dtype=np.int32)
        for c, s, e, _, _ in feats:
            if c == contig:
                arr[s:e + 1] = 1
        ann_pref[contig] = np.cumsum(arr)

    # --- ORFs ------------------------------------------------------------
    orfs: list[ORF] = []
    for contig, seq in contigs.items():
        regs = regions_by_contig[contig]
        # guardamos la secuencia orientada por región para RBS
        for ridx, (rs, re_) in enumerate(regs):
            region_id = f"{contig}_r{ridx}"
            sub = seq[rs - 1:re_]
            for strand in ("+", "-"):
                oriented = sub if strand == "+" else revcomp(sub)
                for hit in _find_orfs(oriented, min_aa, starts, max_variants):
                    if hit.get("partial") and exclude_partial:
                        continue
                    i0, i1 = hit["start_idx"], hit["end_idx"]
                    if strand == "+":
                        gstart, gend = rs + i0, rs + i1 - 1
                    else:
                        gstart, gend = re_ - i1 + 1, re_ - i0
                    pref = ann_pref[contig]
                    covered = int(pref[gend] - pref[gstart - 1])
                    if covered / (gend - gstart + 1) >= 0.95:
                        continue
                    o = ORF(
                        region_id=region_id,
                        contig=contig,
                        strand=strand,
                        frame=i0 % 3,
                        start=gstart,
                        end=gend,
                        aa=hit["aa"],
                        variants=hit["variants"],
                        ambiguous_positions=hit["ambiguous_positions"],
                        start_codon=hit["start_codon"],
                        stop_codon=hit["stop_codon"],
                        n_ambiguous=hit["n_ambiguous"],
                        ambig_stop=hit["ambig_stop"],
                        aas=list(hit["aa"]),
                    )
                    orfs.append(o)
                    # RBS y composición
                    r_score, r_motif = rbs_score(oriented, i0)
                    sc = model.score(oriented[i0:i1])
                    feats_p = protein_features(o.aa)
                    rec = {
                        "genome": genome,
                        "contig": contig,
                        "start": gstart,
                        "end": gend,
                        "strand": strand,
                        "frame": o.frame,
                        "nt_len": o.nt_len,
                        "aa_len": o.aa_len,
                        "start_codon": o.start_codon,
                        "stop_codon": o.stop_codon,
                        "partial": hit.get("partial", False),
                        "n_ambiguous": o.n_ambiguous,
                        "ambig_stop": o.ambig_stop,
                        "variants": len(o.variants),
                        "ann_frac": round(covered / (gend - gstart + 1), 3),
                        "aa": o.aa,
                        "rbs": round(r_score, 4),
                        "rbs_motif": r_motif,
                        "start_w": start_weight(o.start_codon),
                        **{k: round(v, 4) if isinstance(v, float) else v for k, v in sc.items()},
                        **{k: round(v, 4) for k, v in feats_p.items()},
                        "rank": 0.0,
                    }
                    # Un ORF que invade mucho un gen conocido es sospechoso.
                    rec["rank"] = round(rank_score(rec) * (1.0 - 0.6 * rec["ann_frac"]), 4)
                    if rec["partial"]:
                        # ORF que corre hasta el borde: conserva recall pero baja prioridad.
                        rec["rank"] = round(rec["rank"] * 0.7, 4)
                    res.candidates.append(rec)

    res.n_orfs = len(res.candidates)
    return res


def _find_orfs(oriented: str, min_aa: int, starts, max_variants: int):
    from .orfs import find_orfs_oriented

    return find_orfs_oriented(
        oriented, min_aa=min_aa, starts=starts, max_variants=max_variants
    )


# --------------------------------------------------------------------------
# ESM-2 (filtro ligero de proteína)
# --------------------------------------------------------------------------

def add_esm_scores(candidates: list[dict], top_n: int | None = None, batch_size: int = 16) -> int:
    """Puntúa con ESM-2 8M y ajusta el ranking (requiere torch + fair-esm)."""
    from .protein_lm import score_sequences

    sel = sorted(candidates, key=lambda c: c["rank"], reverse=True)
    if top_n:
        sel = sel[:top_n]
    if not sel:
        return 0
    scores = score_sequences([c["aa"] for c in sel], batch_size=batch_size)
    for c, s in zip(sel, scores):
        c["esm_llr"] = round(s, 4)
        # Calibrado con CDS reales (media -0.42) vs barajados (-0.53).
        s_norm = 1.0 / (1.0 + math.exp(-(s + 0.48) * 12.0))
        c["esm_score"] = round(s_norm, 3)
        c["rank"] = round(min(1.0, 0.85 * c["rank"] + 0.15 * s_norm), 4)
    return len(sel)


# --------------------------------------------------------------------------
# DIAMOND
# --------------------------------------------------------------------------

FMT = "qseqid sseqid pident length qlen slen evalue bitscore stitle"


def _acc(name: str) -> str:
    """GCF_000006985.1_ASM698v1 -> GCF_000006985.1"""
    m = re.match(r"(GC[AF]_\d+\.\d+)", name or "")
    return m.group(1) if m else (name or "")


def contig_to_genome(data_root: Path) -> dict[str, str]:
    """Mapa contig -> genoma leyendo la primera directiva ##sequence-region."""
    out: dict[str, str] = {}
    for gff in (data_root / "annotations_ab_initio").glob("*/annotation.gff3"):
        genome = gff.parent.name
        with open(gff) as fh:
            for line in fh:
                if line.startswith("##sequence-region"):
                    out[line.split()[1]] = genome
                    break
    return out


def load_genome_phylum(data_root: Path) -> dict[str, str]:
    path = data_root / "species_phylum.tsv"
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ip = header.index("phylum")
        except ValueError:
            return out
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > ip:
                out[_acc(p[0])] = p[ip]
    return out


def write_faa(candidates: list[dict], path: Path) -> int:
    seen = set()
    n = 0
    with open(path, "w") as fh:
        for c in candidates:
            key = (c["contig"], c["start"], c["end"], c["strand"])
            if key in seen:
                continue
            seen.add(key)
            pid = f"{c['genome']}|{c['contig']}|{c['start']}|{c['end']}|{c['strand']}"
            fh.write(f">{pid}\n")
            aa = c["aa"]
            for i in range(0, len(aa), 60):
                fh.write(aa[i:i + 60] + "\n")
            n += 1
    return n


def run_diamond(query: Path, db: Path, out: Path, threads: int = 4, evalue: float = 1e-5) -> None:
    cmd = [
        "diamond", "blastp",
        "--query", str(query),
        "--db", str(db),
        "--out", str(out),
        "--outfmt", "6", *FMT.split(),
        "--evalue", str(evalue),
        "--max-target-seqs", "1",
        "--threads", str(threads),
        "--quiet",
    ]
    subprocess.run(cmd, check=True)


def parse_diamond(path: Path) -> dict[str, dict]:
    hits: dict[str, dict] = {}
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            q, s, pid, ln, ql, sl, ev, bs = p[:8]
            stitle = "\t".join(p[8:])
            if q in hits:
                continue
            hits[q] = {
                "subject": s,
                "pident": float(pid),
                "length": int(ln),
                "qlen": int(ql),
                "slen": int(sl),
                "evalue": float(ev),
                "bitscore": float(bs),
                "title": stitle,
                "qcov": round(int(ln) / max(1, int(ql)), 3),
            }
    return hits


def attach_hits(
    candidates: list[dict],
    hits: dict[str, dict],
    self_contigs: set[str] | None = None,
    contig_genome: dict[str, str] | None = None,
    phylum: dict[str, str] | None = None,
    query_genome: str | None = None,
) -> None:
    for c in candidates:
        pid = f"{c['genome']}|{c['contig']}|{c['start']}|{c['end']}|{c['strand']}"
        h = hits.get(pid)
        if not h:
            continue
        c["hit"] = h["subject"]
        c["hit_pident"] = h["pident"]
        c["hit_qcov"] = h["qcov"]
        c["hit_title"] = h["title"][:180]
        c["hit_evalue"] = h["evalue"]
        c["hit_alnlen"] = h["length"]
        c["hit_slen"] = h["slen"]
        c["len_ratio"] = round(h["slen"] / max(1, h["qlen"]), 2)
        self_hit = bool(self_contigs) and h["subject"].split("__")[0] in self_contigs
        c["self_hit"] = self_hit
        donor = (contig_genome or {}).get(h["subject"].split("__")[0], "")
        donor_phylum = (phylum or {}).get(_acc(donor), "")
        query_phylum = (phylum or {}).get(_acc(query_genome or ""), "")
        c["hit_genome"] = donor
        c["hit_phylum"] = donor_phylum
        c["same_phylum"] = bool(donor_phylum and query_phylum and donor_phylum == query_phylum)
        if self_hit:
            # Parálogo/pseudogén del propio genoma: no es conservación.
            c["conservation"] = 0.0
            c["rank"] = round(c["rank"] * 0.8, 4)
            continue
        if h["evalue"] > 1e-10 or h["qcov"] < 0.6:
            ev = 0.0
        elif h["pident"] >= 90 and h["qcov"] >= 0.8:
            ev = 1.0
        elif h["pident"] >= 70 and h["qcov"] >= 0.6:
            ev = 0.6
        else:
            ev = 0.3
        if not c["same_phylum"]:
            # Un hit a otro filo para un ORF corto suele ser composicional.
            ev *= 0.25
        # Si el homólogo es mucho más largo, el ORF es un fragmento/dominio.
        if c["len_ratio"] >= 1.5:
            ev *= 0.5
        if c["len_ratio"] >= 2.5:
            ev = 0.0
        c["conservation"] = ev
        c["rank"] = round(min(1.0, c["rank"] + 0.15 * ev), 4)


# --------------------------------------------------------------------------
# Salida
# --------------------------------------------------------------------------

def write_outputs(res: GenomeResult, outdir: str | Path) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / f"{res.genome}_ghost_candidates.tsv"
    js = outdir / f"{res.genome}_ghost_candidates.json"
    summ = outdir / f"{res.genome}_ghost_summary.json"

    cols = [
        "genome", "contig", "start", "end", "strand", "frame", "nt_len", "aa_len",
        "start_codon", "stop_codon", "partial", "n_ambiguous", "ambig_stop", "variants",
        "ann_frac", "codon_llr",
        "frac_bad", "gc3", "rbs", "rbs_motif", "start_w", "gravy", "hydro",
        "disorder", "aromatic", "rank", "hit", "hit_pident", "hit_qcov", "esm_llr", "esm_score", "hit_alnlen", "hit_slen", "len_ratio", "hit_evalue", "hit_title", "self_hit", "hit_genome", "hit_phylum", "same_phylum", "conservation",
    ]
    cands = sorted(res.candidates, key=lambda c: c["rank"], reverse=True)
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in cands:
            row = []
            for k in cols:
                v = c.get(k, "")
                if isinstance(v, float):
                    v = f"{v:.4f}"
                row.append(str(v).replace("\t", " "))
            fh.write("\t".join(row) + "\n")

    with open(js, "w") as fh:
        json.dump({"stats": res.stats(), "candidates": cands}, fh, indent=1)

    by_rank = Counter()
    for c in cands:
        if c["rank"] >= 0.75:
            by_rank["high"] += 1
        elif c["rank"] >= 0.55:
            by_rank["medium"] += 1
        else:
            by_rank["low"] += 1
    doc = res.stats()
    doc["tiers"] = dict(by_rank)
    doc["top"] = [
        {k: c.get(k) for k in ("contig", "start", "end", "strand", "aa_len", "codon_llr", "rank")}
        for c in cands[:20]
    ]
    with open(summ, "w") as fh:
        json.dump(doc, fh, indent=1)
    return {"tsv": tsv, "json": js, "summary": summ}
