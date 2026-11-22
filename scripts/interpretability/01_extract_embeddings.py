#!/usr/bin/env python3
"""Extrae embeddings de 384 dimensiones por ventana genómica y los une a
anotaciones de ADN.

Salida (por defecto en data/embeddings/v2_epoch9_pilot/):
  windows.tsv        una fila por ventana + metadatos + etiquetas
  embeddings.npy     matriz float32 [N, 384] (mismo orden que windows.tsv)
  embeddings.tsv     windows.tsv + dim_000..dim_383 (tabla legible)
  manifest.json      checkpoint/tokenizer (sha256), genomas, parámetros, conteos

Diseño:
  - La ventana es el unidad de análisis (por defecto 1.024 pb).
  - Se muestrean clases: CDS, tRNA, rRNA, ncRNA, intergénico, upstream
    (200 pb río arriba del inicio del CDS, proxy de promotor), oriC,
    CRISPR/profago/transposón y ventanas aleatorias.
  - Las etiquetas `is_*` describen el contenido en el centro de la ventana;
    `sample_class` describe cómo se eligió. No son lo mismo.
  - El embedding es la media de los tokens reales (el padding no participa
    en la atención ni en el pooling).

Uso:
  python scripts/interpretability/01_extract_embeddings.py \
      --n-genomes 6 --max-per-class 120 --device cpu

  # Unos genomas concretos
  python scripts/interpretability/01_extract_embeddings.py \
      --genomes GCF_000006985.1_ASM698v1,GCF_000006685.1_ASM668v1
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import random
import re
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "lib"))

from lib.dna_v2 import embed_sequences, load_v2, tokenize  # noqa: E402
from lib.io import read_fasta  # noqa: E402

FEATURE_TYPES = {
    "CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA",
    "riboswitch", "pseudogene",
}
PRIORITY = ["tRNA", "rRNA", "ncRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA",
            "riboswitch", "pseudogene", "CDS"]


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def parse_attrs(field: str) -> dict[str, str]:
    out = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        else:
            bits = part.split(" ", 1)
            if len(bits) != 2:
                continue
            k, v = bits
        out[k.strip()] = v.strip()
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def parse_gff(path: Path) -> list[dict]:
    """Devuelve features de ADN con coordenadas 0-based half-open."""
    feats: list[dict] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in FEATURE_TYPES:
                continue
            try:
                start = int(p[3]) - 1
                end = int(p[4])
            except ValueError:
                continue
            if end <= start:
                continue
            attrs = parse_attrs(p[8])
            feats.append({
                "seqid": p[0],
                "type": p[2],
                "start": start,
                "end": end,
                "strand": 1 if p[6] == "+" else (-1 if p[6] == "-" else 0),
                "gene": attrs.get("gene") or attrs.get("Name", ""),
                "locus_tag": attrs.get("locus_tag", ""),
                "product": attrs.get("product", ""),
                "ec": attrs.get("EC_number", ""),
            })
    return feats


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(int(s), int(e)) for s, e in merged]


def overlap_frac(start: int, end: int, intervals: list[tuple[int, int]]) -> float:
    if not intervals:
        return 0.0
    total = 0
    for s, e in intervals:
        if e <= start:
            continue
        if s >= end:
            break
        total += max(0, min(e, end) - max(s, start))
    return total / max(1, end - start)


def interval_at(center: int, feats: list[dict]) -> dict | None:
    hits = [f for f in feats if f["start"] <= center < f["end"]]
    if not hits:
        return None
    for typ in PRIORITY:
        for f in hits:
            if f["type"] == typ:
                return f
    return hits[0]


def load_selected(data_root: Path) -> dict[str, dict]:
    path = data_root / "selected.tsv"
    out = {}
    if not path.exists():
        return out
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx.get("phylum", 4):
                continue
            acc = p[idx["assembly_accession"]]
            out[acc] = {
                "phylum": p[idx["phylum"]],
                "organism": p[idx.get("organism_name", 2)],
            }
    return out


def discover_genomes(data_root: Path, include_excluded: bool,
                     gff_source: str = "auto") -> list[dict]:
    raw = {p.name.replace("_genomic.fna.gz", ""): p
           for p in (data_root / "raw").glob("*_genomic.fna.gz")}
    if include_excluded:
        raw.update({p.name.replace("_genomic.fna.gz", ""): p
                    for p in (data_root / "raw_excluded").glob("*_genomic.fna.gz")})
    selected = load_selected(data_root)
    out = []
    for acc, fasta in sorted(raw.items()):
        gff_ncbi = data_root / "annotations" / f"{acc}.gff.gz"
        gff_ab = data_root / "annotations_ab_initio" / acc / "annotation.gff3"
        if gff_source == "ncbi":
            gff = gff_ncbi
        elif gff_source == "ab_initio":
            gff = gff_ab
        else:
            gff = gff_ncbi if gff_ncbi.exists() else gff_ab
        if not gff.exists():
            continue
        # selected.tsv usa el accession base (GCF_000006985.1); los FASTA
        # añaden el sufijo _ASM... o similar.
        m = re.match(r"(GC[AF]_\d+\.\d+)", acc)
        base = m.group(1) if m else acc
        meta = selected.get(base, selected.get(acc, {}))
        split = "train" if (data_root / "raw" / fasta.name).exists() else "excluded"
        out.append({
            "genome": acc,
            "fasta": fasta,
            "gff": gff,
            "phylum": meta.get("phylum", "unknown"),
            "organism": meta.get("organism", ""),
            "split": split,
        })
    return out


def choose_genomes(genomes: list[dict], n: int, per_phylum: int,
                   include_excluded: bool, seed: int) -> list[dict]:
    rng = random.Random(seed)
    pool = [g for g in genomes if include_excluded or g["split"] == "train"]
    by_phylum: dict[str, list[dict]] = defaultdict(list)
    for g in pool:
        by_phylum[g["phylum"]].append(g)
    for v in by_phylum.values():
        rng.shuffle(v)
    chosen = []
    # Ronda 1: hasta per_phylum por filo.
    for phylum in sorted(by_phylum):
        chosen.extend(by_phylum[phylum][:per_phylum])
    # Ronda 2: completar con el resto.
    rest = [g for g in pool if g not in chosen]
    rng.shuffle(rest)
    chosen.extend(rest[: max(0, n - len(chosen))])
    # Preferir diversidad de filo y repartir train/excluded si aplica.
    seen, final = set(), []
    for g in chosen:
        if g["genome"] in seen:
            continue
        seen.add(g["genome"])
        final.append(g)
        if len(final) >= n:
            break
    return final


# ---------------------------------------------------------------------------
# Muestreo de ventanas
# ---------------------------------------------------------------------------
class WindowCollector:
    def __init__(self, window: int, max_n: float = 0.10):
        self.window = window
        self.max_n = max_n
        self.rows: dict[tuple, dict] = {}

    def add(self, genome: dict, contig: str, seq: str, start: int, end: int,
            sample_class: str, center_feature: dict | None = None,
            extra: dict | None = None) -> bool:
        if start < 0 or end > len(seq) or end - start != self.window:
            return False
        w = seq[start:end]
        if w.count("N") / self.window > self.max_n:
            return False
        key = (genome["genome"], contig, start, end)
        row = self.rows.get(key)
        if row is None:
            row = {
                "window_id": f"{genome['genome']}:{contig}:{start + 1}-{end}",
                "genome": genome["genome"],
                "organism": genome["organism"],
                "phylum": genome["phylum"],
                "split": genome["split"],
                "contig": contig,
                "start": start + 1,
                "end": end,
                "window_bp": self.window,
                "sample_class": sample_class,
                "strand": 0,
                "gene": "",
                "locus_tag": "",
                "product": "",
                "ec": "",
                "is_upstream": 0,
                "is_oric": 0,
                "is_crispr": 0,
                "is_prophage": 0,
                "is_transposon": 0,
            }
            self.rows[key] = row
        # Une etiquetas si la ventana ya existía.
        if extra:
            for k, v in extra.items():
                if v:
                    row[k] = 1
        if sample_class not in row["sample_class"].split("|"):
            row["sample_class"] += "|" + sample_class
        if center_feature:
            row["strand"] = center_feature["strand"]
            for k in ("gene", "locus_tag", "product", "ec"):
                if center_feature.get(k):
                    row[k] = center_feature[k]
        return True


def free_runs(covered: np.ndarray, width: int) -> list[tuple[int, int]]:
    """Runs de bases no cubiertas con longitud >= width."""
    # Transiciones False->True y True->False sobre un vector booleano.
    padded = np.concatenate(([True], covered, [True]))
    diff = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(diff == -1)      # True -> False
    ends = np.flatnonzero(diff == 1)         # False -> True
    runs = []
    for s, e in zip(starts, ends):
        if e - s >= width:
            runs.append((int(s), int(e)))
    return runs


def annotate_window(row: dict, seq: str, feats: list[dict],
                    intervals_by_type: dict[str, list[tuple[int, int]]],
                    dist_map: np.ndarray | None) -> None:
    start0, end0 = row["start"] - 1, row["end"]
    center = (start0 + end0) // 2
    feat = interval_at(center, feats)
    row["feature_type"] = feat["type"] if feat else "intergenic"
    row["is_cds"] = int(row["feature_type"] == "CDS")
    row["is_trna"] = int(row["feature_type"] == "tRNA")
    row["is_rrna"] = int(row["feature_type"] == "rRNA")
    row["is_ncrna"] = int(row["feature_type"] in {"ncRNA", "tmRNA", "SRP_RNA",
                                                  "RNase_P_RNA", "riboswitch"})
    row["is_intergenic"] = int(row["feature_type"] == "intergenic")
    if feat:
        row["strand"] = feat["strand"]
        for k in ("gene", "locus_tag", "product", "ec"):
            if feat.get(k) and not row.get(k):
                row[k] = feat[k]

    w = seq[start0:end0]
    gc = (w.count("G") + w.count("C")) / len(w)
    row["gc"] = round(gc, 6)
    g = w.count("G"); c = w.count("C"); a = w.count("A"); t = w.count("T")
    row["gc_skew"] = round((g - c) / (g + c), 6) if g + c else 0.0
    row["at_skew"] = round((a - t) / (a + t), 6) if a + t else 0.0
    counts = np.array([w.count(x) for x in "ACGT"], dtype=float)
    p = counts / counts.sum() if counts.sum() else counts
    nz = p[p > 0]
    row["entropy"] = round(float(-(nz * np.log2(nz)).sum()), 6)
    row["n_frac"] = round(w.count("N") / len(w), 6)

    for typ in ("CDS", "tRNA", "rRNA", "ncRNA"):
        iv = intervals_by_type.get(typ, [])
        row[f"overlap_{typ.lower()}"] = round(overlap_frac(start0, end0, iv), 6)

    if dist_map is not None and 0 <= center < len(dist_map):
        row["dist_to_feature"] = int(dist_map[center])
    else:
        row["dist_to_feature"] = -1


def load_mobile(data_root: Path) -> dict[str, dict]:
    """CRISPR/profago/transposón por genoma, si existen."""
    out: dict[str, dict] = {}
    for name, keys in [
        ("defense_systems.json", ("crispr", "prophage")),
        ("mobile_elements.json", ("transposon",)),
    ]:
        path = data_root / name
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text())
        except Exception:
            continue
        for genome, payload in data.items():
            slot = out.setdefault(genome, {"crispr": [], "prophage": [], "transposon": []})
            if isinstance(payload, list):
                slot["transposon"] = payload
            elif isinstance(payload, dict):
                for k in keys:
                    slot[k] = payload.get(k, [])
    return out


def load_oric(data_root: Path) -> dict[str, int]:
    path = data_root / "orca_predictions_fix" / "oriC_status.tsv"
    out = {}
    if not path.exists():
        return out
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                pos = int(float(p[idx["position"]]))
            except Exception:
                continue
            if pos > 0:
                out[p[idx["genome"]]] = pos
    return out


def build_contig_context(contigs: dict[str, str], feats: list[dict]) -> dict[str, dict]:
    """Precalcula, por contig, secuencia, features, intervalos y distancia."""
    by_contig: dict[str, list[dict]] = defaultdict(list)
    for f in feats:
        by_contig[f["seqid"]].append(f)

    ctx: dict[str, dict] = {}
    for contig, seq in contigs.items():
        cfeats = by_contig.get(contig, [])
        if not cfeats:
            continue
        intervals: dict[str, list[tuple[int, int]]] = {}
        for typ in ("CDS", "tRNA", "rRNA"):
            intervals[typ] = merge_intervals(
                [(f["start"], f["end"]) for f in cfeats if f["type"] == typ]
            )
        intervals["ncRNA"] = merge_intervals([
            (f["start"], f["end"]) for f in cfeats
            if f["type"] in {"ncRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA", "riboswitch"}
        ])
        all_iv = merge_intervals([(f["start"], f["end"]) for f in cfeats])
        covered = np.zeros(len(seq), dtype=bool)
        for s, e in all_iv:
            covered[s:e] = True
        try:
            from scipy.ndimage import distance_transform_edt
            dist = distance_transform_edt(~covered).astype(np.int32)
        except Exception:
            dist = None
        ctx[contig] = {"seq": seq, "feats": cfeats,
                       "intervals": intervals, "dist": dist}
    return ctx


def sample_genome(genome: dict, contigs: dict[str, str], collector: WindowCollector,
                  rng: random.Random, max_per_class: int,
                  oric: dict[str, int], mobile: dict[str, dict]) -> None:
    feats = parse_gff(genome["gff"])
    if not feats:
        print(f"  [skip] sin features: {genome['genome']}")
        return
    ctx = build_contig_context(contigs, feats)
    if not ctx:
        print(f"  [skip] sin contigs anotados: {genome['genome']}")
        return

    half = collector.window // 2

    def contig_at(pos: int):
        for contig, c in ctx.items():
            if 0 <= pos < len(c["seq"]):
                return contig, c
        return None, None

    # 1) Features anotadas + 2) upstream, por contig.
    for contig, c in ctx.items():
        seq = c["seq"]
        by_type: dict[str, list[dict]] = defaultdict(list)
        for f in c["feats"]:
            by_type[f["type"]].append(f)

        def add_centered(f: dict, sample_class: str, extra=None):
            mid = (f["start"] + f["end"]) // 2
            collector.add(genome, contig, seq, mid - half, mid + half,
                          sample_class, center_feature=f, extra=extra)

        for typ in ("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"):
            candidates = by_type.get(typ, [])
            if not candidates:
                continue
            rng.shuffle(candidates)
            for f in candidates[:max_per_class]:
                add_centered(f, typ.lower())

        cds = list(by_type.get("CDS", []))
        rng.shuffle(cds)
        for f in cds[:max_per_class]:
            tss = f["start"] if f["strand"] >= 0 else f["end"]
            center = tss - 200 if f["strand"] >= 0 else tss + 200
            collector.add(genome, contig, seq, center - half, center + half,
                          "upstream", center_feature=f, extra={"is_upstream": 1})

        # 3) Intergénico: centro sin ninguna feature. La ventana puede
        # solapar genes vecinos; `overlap_*` lo cuantifica.
        covered = np.zeros(len(seq), dtype=bool)
        for s, e in merge_intervals([(f["start"], f["end"]) for f in c["feats"]]):
            covered[s:e] = True
        added, attempts = 0, 0
        while added < max_per_class and attempts < max_per_class * 80:
            start = rng.randrange(0, max(1, len(seq) - collector.window))
            center = start + half
            if 0 <= center < len(seq) and not covered[center]:
                if collector.add(genome, contig, seq, start, start + collector.window,
                                 "intergenic"):
                    added += 1
            attempts += 1

        # 4) Ventanas aleatorias (pueden caer en cualquier parte).
        for _ in range(min(max_per_class, 50)):
            start = rng.randrange(0, max(1, len(seq) - collector.window))
            collector.add(genome, contig, seq, start, start + collector.window,
                          "random")

    # 5) oriC (una posición; se asigna al contig que la contiene).
    pos = oric.get(genome["genome"])
    if pos:
        contig, c = contig_at(pos)
        if c is not None:
            collector.add(genome, contig, c["seq"], pos - half, pos + half,
                          "oric", extra={"is_oric": 1})

    # 6) CRISPR / profago / transposón.
    mob = mobile.get(genome["genome"], {})
    for key, cls, flag in [
        ("crispr", "crispr", "is_crispr"),
        ("prophage", "prophage", "is_prophage"),
        ("transposon", "transposon", "is_transposon"),
    ]:
        items = list(mob.get(key, []) or [])
        rng.shuffle(items)
        for item in items[:max_per_class]:
            try:
                s = int(item["start"]) - 1
                e = int(item["end"])
            except Exception:
                continue
            contig, c = contig_at((s + e) // 2)
            if c is None:
                continue
            center = (s + e) // 2
            collector.add(genome, contig, c["seq"], center - half, center + half,
                          cls, extra={flag: 1})

    # Anotar todas las ventanas de este genoma usando el contexto precalculado.
    for row in collector.rows.values():
        if row["genome"] != genome["genome"] or row.get("feature_type"):
            continue
        c = ctx.get(row["contig"])
        if c is not None:
            annotate_window(row, c["seq"], c["feats"], c["intervals"], c["dist"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(ROOT / "data"))
    ap.add_argument("--out", default=str(ROOT / "data" / "embeddings" / "v2_epoch9_pilot"))
    ap.add_argument("--checkpoint", default=str(ROOT / "dna_model_v2_epoch9.pt"))
    ap.add_argument("--tokenizer", default=str(ROOT / "data" / "tokenizer_v2" / "dna_sp_v2.model"))
    ap.add_argument("--window", type=int, default=1024)
    ap.add_argument("--max-per-class", type=int, default=120)
    ap.add_argument("--n-genomes", type=int, default=6, help="0 = todos los disponibles")
    ap.add_argument("--per-phylum", type=int, default=1)
    ap.add_argument("--genomes", default="", help="Lista separada por comas; vacío = automático")
    ap.add_argument("--include-excluded", action="store_true",
                    help="Incluye parásitos/simbiontes de data/raw_excluded")
    ap.add_argument("--gff-source", choices=["ncbi", "ab_initio", "auto"],
                    default="auto", help="Anotación de entrada (default: NCBI, fallback ab initio)")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--max-tokens", type=int, default=8192)
    ap.add_argument("--max-n", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-tsv", action="store_true", help="No escribir embeddings.tsv (solo .npy)")
    args = ap.parse_args()

    t0 = time.time()
    data_root = Path(args.data_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Descubriendo genomas en {data_root} ...")
    genomes = discover_genomes(data_root, args.include_excluded, args.gff_source)
    print(f"  genomas con FASTA+GFF: {len(genomes)}")

    if args.genomes.strip():
        wanted = [x.strip() for x in args.genomes.split(",") if x.strip()]
        index = {g["genome"]: g for g in genomes}
        missing = [g for g in wanted if g not in index]
        if missing:
            raise SystemExit(f"No encontrados (FASTA+GFF): {missing}")
        selected = [index[g] for g in wanted]
    else:
        n = len(genomes) if args.n_genomes <= 0 else args.n_genomes
        selected = choose_genomes(genomes, n, args.per_phylum,
                                  args.include_excluded, args.seed)
    if not selected:
        raise SystemExit("No hay genomas seleccionables.")

    print("Genomas seleccionados:")
    for g in selected:
        print(f"  - {g['genome']}  [{g['phylum']}]  split={g['split']}")

    print(f"Cargando modelo: {Path(args.checkpoint).name} ({args.device}) ...")
    model, sp = load_v2(args.checkpoint, args.tokenizer, device=args.device)

    print("Cargando secuencias de los genomas seleccionados ...")
    genome_seqs: dict[str, dict[str, str]] = {}
    for g in selected:
        genome_seqs[g["genome"]] = {name: seq for name, seq in read_fasta(g["fasta"])}
        n_bp = sum(len(s) for s in genome_seqs[g["genome"]].values())
        print(f"  {g['genome']}: {len(genome_seqs[g['genome']])} contig(s), {n_bp:,} pb")

    collector = WindowCollector(window=args.window, max_n=args.max_n)
    oric = load_oric(data_root)
    mobile = load_mobile(data_root)

    for i, g in enumerate(selected, 1):
        t = time.time()
        rng = random.Random(args.seed + i)
        before = len(collector.rows)
        sample_genome(g, genome_seqs[g["genome"]], collector, rng,
                      args.max_per_class, oric, mobile)
        print(f"[{i}/{len(selected)}] {g['genome']}: "
              f"+{len(collector.rows) - before} ventanas ({time.time() - t:.1f}s)")

    rows = sorted(collector.rows.values(), key=lambda r: (r["genome"], r["contig"], r["start"]))
    if not rows:
        raise SystemExit("No se generó ninguna ventana.")
    print(f"Total ventanas: {len(rows)}")

    seqs = []
    for r in rows:
        try:
            seq = genome_seqs[r["genome"]][r["contig"]]
        except KeyError as exc:
            raise RuntimeError(f"Contig no encontrado: {r['contig']} en {r['genome']}") from exc
        seqs.append(seq[r["start"] - 1:r["end"]])

    print(f"Extrayendo embeddings (batch={args.batch_size}, device={args.device}) ...")
    emb, lengths = embed_sequences(
        model, sp, seqs,
        batch_size=args.batch_size, device=args.device,
        max_tokens=args.max_tokens, pool="mean", show_progress=True,
    )
    for r, n in zip(rows, lengths.tolist()):
        r["n_tokens"] = int(n)

    meta = pd.DataFrame(rows)
    # Orden de columnas legible.
    front = ["window_id", "genome", "organism", "phylum", "split", "contig",
             "start", "end", "window_bp", "sample_class", "feature_type",
             "strand", "gene", "locus_tag", "product", "ec",
             "is_cds", "is_trna", "is_rrna", "is_ncrna", "is_intergenic",
             "is_upstream", "is_oric", "is_crispr", "is_prophage",
             "is_transposon", "overlap_cds", "overlap_trna", "overlap_rrna",
             "overlap_ncrna", "dist_to_feature", "gc", "gc_skew", "at_skew",
             "entropy", "n_frac", "n_tokens"]
    meta = meta[[c for c in front if c in meta.columns] +
                [c for c in meta.columns if c not in front]]
    meta.to_csv(out_dir / "windows.tsv", sep="\t", index=False)

    np.save(out_dir / "embeddings.npy", emb.astype(np.float32))

    if not args.no_tsv:
        dims = pd.DataFrame(
            emb.astype(np.float32),
            columns=[f"dim_{i:03d}" for i in range(emb.shape[1])],
        )
        pd.concat([meta.reset_index(drop=True), dims], axis=1).to_csv(
            out_dir / "embeddings.tsv", sep="\t", index=False,
            float_format="%.6g",
        )

    counts = meta["sample_class"].str.split("|").explode().value_counts().to_dict()
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "checkpoint_sha256": sha256(Path(args.checkpoint)),
        "tokenizer": str(Path(args.tokenizer).resolve()),
        "tokenizer_sha256": sha256(Path(args.tokenizer)),
        "window_bp": args.window,
        "max_per_class": args.max_per_class,
        "max_n": args.max_n,
        "seed": args.seed,
        "device": args.device,
        "n_windows": int(len(meta)),
        "n_genomes": int(meta["genome"].nunique()),
        "genomes": sorted(meta["genome"].unique().tolist()),
        "split_counts": meta["split"].value_counts().to_dict(),
        "sample_class_counts": counts,
        "feature_type_counts": meta["feature_type"].value_counts().to_dict(),
        "label_counts": {
            c: int(meta[c].sum()) for c in
            ["is_cds", "is_trna", "is_rrna", "is_ncrna", "is_intergenic",
             "is_upstream", "is_oric", "is_crispr", "is_prophage", "is_transposon"]
            if c in meta.columns
        },
        "seconds": round(time.time() - t0, 1),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))

    print(f"\nOK -> {out_dir}")
    print(f"  windows.tsv      {len(meta)} filas")
    print(f"  embeddings.npy   {emb.shape}")
    if not args.no_tsv:
        print(f"  embeddings.tsv   {emb.shape[0]} x {emb.shape[1] + len(meta.columns)}")
    print(f"  clases: {counts}")
    print(f"  tiempo: {manifest['seconds']}s")


if __name__ == "__main__":
    main()
