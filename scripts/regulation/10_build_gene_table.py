#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10_build_gene_table.py — tabla génica completa por genoma (cimiento de L2/L3).

Fuentes:
  * `data/annotations_ab_initio/<G>/annotation.gff3`  CDS de pyrodigal
  * `data/annotations_ab_initio/<G>/functional.tsv`   DIAMOND vs Swiss-Prot
  * `data/annotations/<G>.gff.gz`                     ARN de NCBI (con hebra)
  * `data/raw/<G>_genomic.fna.gz`                     secuencia (opcional)
  * `--gbff`                                          GenBank completo (E. coli)

Salida: `data/regulation/<genoma>/genes.tsv` + `meta.json`.

Uso:
  python scripts/regulation/10_build_gene_table.py --all
  python scripts/regulation/10_build_gene_table.py --genome GCF_000006685.1_ASM668v1
  python scripts/regulation/10_build_gene_table.py --gbff data/test_ecoli.gbff
  python scripts/regulation/10_build_gene_table.py --selftest
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

COLS = ["genome", "contig", "start", "end", "strand", "gene_id", "locus_tag",
        "gene_name", "product", "ec", "acc", "fam", "partial", "rbs_motif",
        "rbs_source", "start_type", "conf", "rna", "gene_type", "source",
        "status", "evidence", "attempted", "confidence"]


def _row(genome_id: str, g: C.Gene) -> dict:
    if g.is_rna:
        d = C.decide("consensus", validated_by="ncbi_gff",
                     reason="ARN anotado por NCBI; una sola fuente",
                     attempted="ncbi_gff", confidence="high",
                     evidence=f"global:{g.contig}:{g.start}-{g.end}:{g.strand:+d}")
    elif g.partial not in ("", "00"):
        d = C.decide("anomalous", validated_by="pyrodigal",
                     reason=f"CDS parcial (partial={g.partial})",
                     attempted="pyrodigal", confidence="low",
                     evidence=f"rbs:{g.rbs_motif or 'none'}")
    elif g.acc:
        d = C.decide("validated", validated_by="pyrodigal;diamond_swissprot",
                     reason="coordenadas y función con métodos distintos",
                     attempted="pyrodigal;diamond", confidence="high",
                     evidence=f"sp:{g.acc}")
    else:
        d = C.decide("consensus", validated_by="pyrodigal",
                     reason="CDS predicho sin homología específica",
                     attempted="pyrodigal;diamond", confidence="medium",
                     evidence=f"rbs:{g.rbs_motif or 'none'}")
    r = {"genome": genome_id, "contig": g.contig, "start": g.start,
         "end": g.end, "strand": g.strand, "gene_id": g.gene_id,
         "locus_tag": g.locus_tag, "gene_name": g.gene_name,
         "product": g.product, "ec": g.ec, "acc": g.acc, "fam": g.fam,
         "partial": g.partial,
         "rbs_motif": g.rbs_motif, "rbs_source": g.rbs_source,
         "start_type": g.start_type, "conf": g.conf,
         "rna": "true" if g.is_rna else "false", "gene_type": g.gene_type,
         "source": g.source}
    return C.status_row(r, d)


def _diamond(genome_id: str, faa: Path, cache: Path) -> dict:
    """Anota con DIAMOND vs Swiss-Prot (BD local) y cachea el resultado."""
    if cache.exists():
        return C.parse_functional(cache)
    db = C.DATA / "db" / "swissprot.dmnd"
    if not db.exists() or not faa.exists():
        return {}
    cache.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["diamond", "blastp", "--query", str(faa), "--db", str(db),
                    "--out", str(cache), "--outfmt", "6", "qseqid", "sseqid",
                    "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                    "sstart", "send", "evalue", "bitscore", "stitle",
                    "--evalue", "1e-5", "--max-target-seqs", "1",
                    "--threads", "4", "--quiet"], check=True)
    return C.parse_functional(cache)


def _faa_from_gbff(gbff: Path, out: Path) -> Path:
    from Bio import SeqIO
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for rec in SeqIO.parse(str(gbff), "genbank"):
            for f in rec.features:
                if f.type == "CDS" and "translation" in f.qualifiers:
                    q = f.qualifiers
                    gid = (q.get("protein_id") or [""])[0] or \
                          (q.get("locus_tag") or ["x"])[0]
                    fh.write(f">{gid}\n{q['translation'][0]}\n")
    return out


def build_one(genome_id: str, out_dir: Path, with_seq: bool,
              gbff: Path | None = None, use_diamond: bool = True) -> dict:
    t0 = time.time()
    if gbff:
        diamond = {}
        if use_diamond:
            faa = _faa_from_gbff(gbff, out_dir / "proteins.faa")
            cache = out_dir / "diamond_swissprot.tsv"
            diamond = _diamond(genome_id, faa, cache)
        g = C.load_gbff(gbff, genome_id, diamond)
    else:
        g = C.load_genome(genome_id, with_seq=with_seq)
        if g is None:
            raise SystemExit(f"genoma no encontrado: {genome_id}")
    rows = [_row(genome_id, x) for x in g.genes]
    rows.sort(key=lambda r: (r["contig"], int(r["start"]), int(r["end"])))
    C.write_tsv(out_dir / "genes.tsv", rows, COLS)
    n_cds = sum(1 for x in g.genes if not x.is_rna)
    n_rna = len(g.genes) - n_cds
    meta = {
        "genome": genome_id, "kind": g.kind, "contigs": g.contigs,
        "n_genes": len(g.genes), "n_cds": n_cds, "n_rna": n_rna,
        "n_rna_with_strand": sum(1 for x in g.genes if x.is_rna and x.strand),
        "n_cds_with_acc": sum(1 for x in g.genes if not x.is_rna and x.acc),
        "n_cds_with_ec": sum(1 for x in g.genes if not x.is_rna and x.ec),
        "n_cds_with_name": sum(1 for x in g.genes if not x.is_rna and x.gene_name),
        "has_sequence": g.has_seq(),
        "seq_source": ({"kind": "gbff", "path": str(Path(gbff).resolve())} if gbff else
                       ({"kind": "raw", "path": str((C.RAW / f"{genome_id}_genomic.fna.gz").resolve())}
                        if g.has_seq() else {"kind": "none", "path": ""})),
        "golden": C.AB / genome_id / "functional.tsv",
        "source": g.source,
        "elapsed_s": round(time.time() - t0, 2),
        "sha256_genes": C.sha256_files([out_dir / "genes.tsv"]),
        "status_counts": _counts(rows),
    }
    C.dump_json(out_dir / "genes_meta.json", meta)
    return meta


def _counts(rows: list[dict]) -> dict:
    out: dict[str, int] = {}
    for r in rows:
        out[r["status"]] = out.get(r["status"], 0) + 1
    return dict(sorted(out.items()))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--gbff", type=Path)
    ap.add_argument("--name", default="ecoli_k12")
    ap.add_argument("--with-seq", action="store_true")
    ap.add_argument("--no-diamond", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        return selftest()
    if a.gbff:
        m = build_one(a.name, C.REG / a.name, True, gbff=a.gbff,
                      use_diamond=not a.no_diamond)
        print(f"{m['genome']}: {m['n_cds']} CDS, {m['n_rna']} ARN, "
              f"acc={m['n_cds_with_acc']}, ec={m['n_cds_with_ec']}")
        return 0

    if a.all:
        dirs = sorted(p for p in C.AB.glob("GCF_*") if p.is_dir())
        if a.limit:
            dirs = dirs[:a.limit]
        t0 = time.time()
        ok = 0
        for i, d in enumerate(dirs, 1):
            try:
                m = build_one(d.name, C.REG / d.name, with_seq=a.with_seq)
                ok += 1
                if i % 50 == 0:
                    print(f"  [{i}/{len(dirs)}] {d.name}: {m['n_cds']} CDS",
                          flush=True)
            except Exception as e:                        # noqa: BLE001
                print(f"  FALLO {d.name}: {e}", file=sys.stderr)
        print(f"genes.tsv: {ok}/{len(dirs)} genomas en {time.time()-t0:.1f}s")
        return 0 if ok else 1

    if not a.genome:
        ap.error("indica --genome, --all o --gbff")
    m = build_one(a.genome, C.REG / a.genome, a.with_seq)
    print(f"{m['genome']}: {m['n_cds']} CDS, {m['n_rna']} ARN "
          f"(hebra {m['n_rna_with_strand']}), acc={m['n_cds_with_acc']}")
    return 0


def selftest() -> int:
    import tempfile
    ok = 0
    # reverse complement
    assert C.revcomp("ACGTN") == "NACGT", C.revcomp("ACGTN")
    # scan con degeneración
    w = "GGGGTTGACATATAATGGGG"
    hit = C.best_site("TTGACA", w)
    assert hit and hit[0] == 4 and hit[1] == 0, hit
    assert C.best_site("TTGTCA", w) is None          # sin semilla no hay sitio
    hit2 = C.best_site("TTGWCA", w)
    assert hit2 and hit2[0] == 4, hit2
    # gap en ambas hebras
    a = C.Gene("c", 100, 200, 1, "a")
    b = C.Gene("c", 230, 300, 1, "b")
    assert C.gap_between(a, b) == 30
    am = C.Gene("c", 100, 200, -1, "a")
    bm = C.Gene("c", 40, 85, -1, "b")
    assert C.gap_between(am, bm) == 15
    # orden de transcripción
    g = C.Genome("x", "gbff")
    g.contigs = {"c": 400}
    g.genes = [C.Gene("c", 1, 100, 1, "1"), C.Gene("c", 120, 200, 1, "2"),
               C.Gene("c", 300, 350, 1, "r", gene_type="tRNA"),
               C.Gene("c", 360, 400, 1, "3")]
    units = C.transcription_units(g)
    assert len(units) == 2 and [x.gene_id for x in units[0]] == ["1", "2"], units
    # RBS
    rbs = C.detect_rbs("A" * 50 + "GGAGG" + "A" * 10 + "ATG" + "A" * 30,
                       C.Gene("c", 66, 95, 1, "g"))
    assert rbs in ("GGAGG", "AGGAG"), rbs
    assert not C.has_strong_rbs("None") and C.has_strong_rbs("AGGA/GGAG")
    # prefijo génico
    g1, g2 = C.Gene("c", 1, 9, 1, "1", gene_name="lacZ"), \
             C.Gene("c", 1, 9, 1, "2", gene_name="lacY")
    assert C.coherent(g1, g2)[0] and C.gene_prefix(g1, g2) == "lac"
    assert not C.coherent(C.Gene("c", 1, 9, 1, "3", gene_name="ybaB"),
                          C.Gene("c", 1, 9, 1, "4", gene_name="ybaC"))[0]
    # estados
    assert "validated" in C.STATUSES and "novel" in C.STATUSES
    with tempfile.TemporaryDirectory(dir=str(C.REG)) as td:
        p = Path(td) / "t.tsv"
        C.write_tsv(p, [{"a": 1, "b": "x\ty"}], ["a", "b"])
        assert C.read_tsv(p)[0]["b"] == "x y"
        ok = 1
    print(f"selftest 10_build_gene_table: OK ({ok} bloque)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
