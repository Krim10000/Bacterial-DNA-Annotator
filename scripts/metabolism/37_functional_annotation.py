#!/usr/bin/env python3
"""
37_functional_annotation.py — Anotación funcional por homología (DIAMOND vs Swiss-Prot).

Recibe el directorio de un genoma ya anotado ab initio (salida de 32):
    data/annotations_ab_initio/<prefix>/
        proteins.faa       (CDS predichos por Pyrodigal)
        annotation.gbk     (GenBank estructural, solo coordenadas)

Y produce en el mismo directorio:
    functional_annotation.json  mapa CDS -> {gene, product, uniprot, ec, ...} + cobertura
    annotation_func.gbk         copia de annotation.gbk con qualifiers gene=/product=/
                                db_xref=/note= añadidos a cada CDS con hit.
                                NO modifica annotation.gbk.

Motivación (hallazgo del test de rotación): ORCA necesita dnaA/dnaN con nombre
para anclar el oriC; el GenBank ab initio solo trae coordenadas y 35 cae en
needs_functional_annotation. El GenBank enriquecido de este script es el que
35 debe consumir con --gbk-name annotation_func.gbk.

Requisitos (una sola vez):
    sudo apt install diamond-aligner      # o binario de GitHub
    mkdir -p data/db
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -P data/db/
    gunzip data/db/uniprot_sprot.fasta.gz
    diamond makedb --in data/db/uniprot_sprot.fasta -d data/db/swissprot

Uso:
    python scripts/metabolism/37_functional_annotation.py data/annotations_ab_initio/GCF_000006685.1_ASM668v1
    python scripts/metabolism/37_functional_annotation.py <dir> --threads 4
    python scripts/metabolism/37_functional_annotation.py --selftest    # prueba del parser
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("func_annot")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = PROJECT_ROOT / "data" / "db" / "swissprot"

# --- Parser de encabezados UniProt (fasta header = DE line compactada) ----

RE_OS_OX = re.compile(r" OS=(.*?) OX=")
RE_OS_END = re.compile(r" OS=(.*)$")
RE_GN = re.compile(r" GN=([^\s]+)")
RE_EC_PAREN = re.compile(r"\(EC[ :]*([\d\.]+\.-?[\d\-\.]*)\)")
RE_EC_EQ = re.compile(r"\bEC[=:]\s*([\d\.]+\.-?[\d\-\.]*)")
RE_ACC = re.compile(r"\b((?:G[CF]A?_\d+)|(?:[A-Z]{2}_\d+))\.\d+\b")


@dataclass
class Hit:
    cds_id: str
    uniprot: str = ""
    entry_name: str = ""
    gene: str = ""
    product: str = ""
    organism: str = ""
    ec: str = ""
    pident: float = 0.0
    aln_len: int = 0
    evalue: float = 1.0
    bitscore: float = 0.0


def parse_sprot_hit(qseqid: str, sseqid: str, stitle: str,
                    pident: float, aln_len: int, evalue: float,
                    bitscore: float) -> Hit:
    """Convierte una fila m8 de DIAMOND vs uniprot_sprot en un Hit estructurado."""
    hit = Hit(cds_id=qseqid, pident=pident, aln_len=aln_len,
              evalue=evalue, bitscore=bitscore)

    parts = sseqid.split("|")
    if len(parts) >= 3:
        hit.uniprot = parts[1]
        hit.entry_name = parts[2]
    else:
        hit.entry_name = sseqid

    m = RE_OS_OX.search(stitle) or RE_OS_END.search(stitle)
    if m:
        hit.organism = m.group(1).strip()

    m = RE_GN.search(stitle)
    if m:
        # GN puede traer sinonimos: "dnaA, rhaA" -> primero
        hit.gene = m.group(1).split(",")[0].strip()
    if not hit.gene and "_" in hit.entry_name:
        # DNAA_ECOLI -> dnaa (ORCA compara en minusculas)
        hit.gene = hit.entry_name.split("_")[0].lower()

    cut = stitle.find(" OS=")
    hit.product = stitle[:cut].strip() if cut > 0 else stitle.strip()

    m = RE_EC_PAREN.search(stitle) or RE_EC_EQ.search(stitle)
    if m:
        hit.ec = m.group(1)
    return hit


def selftest() -> int:
    samples = [
        ("NC_002620.2__NC_002620.2_1", "sp|P0304|DNAA_ECOLI",
         "Chromosomal replication initiator protein DnaA "
         "OS=Escherichia coli (strain K12) OX=83333 GN=dnaA PE=1 SV=4",
         {"gene": "dnaA", "uniprot": "P0304", "ec": ""}),
        ("x__2", "sp|P28990|DNAN_BACSU",
         "DNA polymerase sliding clamp subunit DnaN "
         "OS=Bacillus subtilis (strain 168) OX=224308 GN=dnaN PE=1 SV=1",
         {"gene": "dnaN", "uniprot": "P28990", "ec": ""}),
        ("x__3", "sp|X12345|TUS_ECOLI",
         "Transcription termination protein Tus "
         "OS=Escherichia coli (strain K12) OX=83333 PE=1 SV=2",
         {"gene": "tus", "uniprot": "X12345", "ec": ""}),
        ("x__4", "sp|X99999|AHHX_PSEAE",
         "Anthranilate synthase component II (EC 4.1.3.27) "
         "OS=Pseudomonas aeruginosa PAO1 OX=208964 GN=trpG PE=1 SV=1",
         {"gene": "trpG", "uniprot": "X99999", "ec": "4.1.3.27"}),
    ]
    ok = True
    for q, sid, title, want in samples:
        h = parse_sprot_hit(q, sid, title, 99.0, 400, 0.0, 800.0)
        got = {"gene": h.gene, "uniprot": h.uniprot, "ec": h.ec}
        status = "OK " if got == want else "FAIL"
        if got != want:
            ok = False
        log.info("%s gene=%s uniprot=%s ec=%s product=%r org=%r",
                 status, h.gene, h.uniprot, h.ec, h.product[:40], h.organism[:30])
    return 0 if ok else 1


# --- DIAMOND ---------------------------------------------------------------

def run_diamond(faa: Path, db: Path, threads: int, evalue: str,
                max_target_seqs: int, out_tsv: Path) -> None:
    if shutil.which("diamond") is None:
        raise RuntimeError(
            "diamond no esta en PATH. Instalalo: sudo apt install diamond-aligner "
            "(o descarga el binario de GitHub)"
        )
    db_path = db if str(db).endswith(".dmnd") else Path(str(db) + ".dmnd")
    if not db_path.exists():
        raise RuntimeError(
            f"No existe la DB DIAMOND: {db_path}. Creala con: "
            "diamond makedb --in data/db/uniprot_sprot.fasta -d data/db/swissprot"
        )
    cmd = [
        "diamond", "blastp",
        "-d", str(db_path), "-q", str(faa), "-o", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "evalue", "bitscore", "stitle",
        "--evalue", evalue,
        "--max-target-seqs", str(max_target_seqs),
        "--threads", str(threads),
        "--quiet",
    ]
    log.info("  diamond blastp: %s vs %s", faa.name, db_path.name)
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if res.returncode != 0:
        raise RuntimeError("diamond fallo: " + (res.stderr or res.stdout or "")[-500:])


def parse_m8(tsv: Path) -> dict[str, Hit]:
    hits: dict[str, Hit] = {}
    with open(tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            q, s = cols[0], cols[1]
            try:
                hit = parse_sprot_hit(
                    q, s, cols[6],
                    pident=float(cols[2]),
                    aln_len=int(cols[3]),
                    evalue=float(cols[4]),
                    bitscore=float(cols[5]),
                )
            except ValueError:
                continue
            # max-target-seqs 1 => 1 fila por query, pero por si acaso:
            if q not in hits or hit.bitscore > hits[q].bitscore:
                hits[q] = hit
    return hits


def count_faa(faa: Path) -> int:
    n = 0
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


# --- GenBank ---------------------------------------------------------------

def enrich_genbank(gbk_in: Path, hits: dict[str, Hit], gbk_out: Path) -> tuple[int, int]:
    """Copia annotation.gbk anadiendo gene=/product=/db_xref=/note= a los CDS con hit."""
    try:
        from Bio import SeqIO
    except ImportError as e:
        raise RuntimeError("Biopython no instalado (pip install biopython)") from e

    records = list(SeqIO.parse(str(gbk_in), "genbank"))
    n_cds = 0
    n_enriched = 0
    for rec in records:
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            n_cds += 1
            cds_id = (feat.qualifiers.get("ID") or [None])[0]
            hit = hits.get(cds_id) if cds_id else None
            if hit is None:
                continue
            if hit.gene:
                feat.qualifiers["gene"] = [hit.gene]
            if hit.product:
                feat.qualifiers["product"] = [hit.product]
            else:
                feat.qualifiers["product"] = ["hypothetical protein"]
            xrefs = list(feat.qualifiers.get("db_xref", []))
            if hit.uniprot:
                xrefs.append(f"UniProtKB:{hit.uniprot}")
            if hit.ec:
                xrefs.append(f"EC:{hit.ec}")
            feat.qualifiers["db_xref"] = xrefs
            notes = list(feat.qualifiers.get("note", []))
            notes.append(
                f"DIAMOND vs Swiss-Prot: {hit.uniprot} "
                f"pident={hit.pident:.1f}% evalue={hit.evalue:.1e}"
            )
            feat.qualifiers["note"] = notes
            n_enriched += 1

    SeqIO.write(records, str(gbk_out), "genbank")
    return n_cds, n_enriched


# --- CLI -------------------------------------------------------------------

@dataclass
class Report:
    genome_dir: str
    timestamp: str
    db: str
    n_cds: int
    n_hits: int
    match_rate: float
    with_gene: int
    with_product: int
    with_ec: int
    mean_pident: float
    gbk_cds: int
    gbk_enriched: int
    func_gbk: str | None
    warnings: list = None  # type: ignore[assignment]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Anotacion funcional DIAMOND vs Swiss-Prot para un genoma",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("genome_dir", type=Path, nargs="?",
                   help="Directorio del genoma (con proteins.faa y annotation.gbk)")
    p.add_argument("--db", type=Path, default=DEFAULT_DB,
                   help="Ruta de la DB DIAMOND (sin .dmnd o con el)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--evalue", default="1e-5")
    p.add_argument("--max-target-seqs", type=int, default=1)
    p.add_argument("--skip-gbk", action="store_true",
                   help="Solo DIAMOND + JSON; no genera annotation_func.gbk")
    p.add_argument("--force", action="store_true",
                   help="Re-corre aunque exista functional_annotation.json")
    p.add_argument("--selftest", action="store_true",
                   help="Prueba el parser de encabezados UniProt y sale")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    if args.genome_dir is None:
        log.error("Falta <genome_dir> (o usa --selftest)")
        return 2

    gdir = args.genome_dir.resolve()
    faa = gdir / "proteins.faa"
    gbk = gdir / "annotation.gbk"
    out_json = gdir / "functional_annotation.json"
    out_gbk = gdir / "annotation_func.gbk"

    if not gdir.is_dir():
        log.error("No existe: %s", gdir)
        return 2
    if not faa.exists():
        log.error("Falta %s (salida de 32)", faa)
        return 2
    if out_json.exists() and not args.force:
        log.info("SKIP (ya existe %s; usa --force)", out_json.name)
        return 0
    if not args.skip_gbk and not gbk.exists():
        log.warning("Falta %s: se omite la generacion de annotation_func.gbk", gbk)
        args.skip_gbk = True

    warnings: list[str] = []
    tmp = tempfile.mkdtemp(prefix="func_annot_")
    try:
        tsv = Path(tmp) / "diamond.tsv"
        run_diamond(faa, args.db, args.threads, args.evalue,
                    args.max_target_seqs, tsv)
        hits = parse_m8(tsv)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    n_cds = count_faa(faa)
    n_hits = len(hits)
    rate = n_hits / n_cds if n_cds else 0.0
    with_gene = sum(1 for h in hits.values() if h.gene)
    with_prod = sum(1 for h in hits.values() if h.product)
    with_ec = sum(1 for h in hits.values() if h.ec)
    mean_pid = (sum(h.pident for h in hits.values()) / n_hits) if n_hits else 0.0

    log.info("  CDS=%d | hits=%d (%.1f%%) | gene=%.1f%% product=%.1f%% EC=%.1f%% | pident medio=%.1f%%",
             n_cds, n_hits, 100 * rate, 100 * with_gene / n_cds if n_cds else 0,
             100 * with_prod / n_cds if n_cds else 0,
             100 * with_ec / n_cds if n_cds else 0, mean_pid)

    if rate < 0.5:
        warnings.append(
            f"match rate bajo ({100*rate:.1f}%): los IDs de proteins.faa no "
            "coinciden con los qualifiers ID= del GenBank; revisar el prefijo "
            "contig__ que anade 32"
        )
        log.warning("  %s", warnings[-1])

    gbk_cds = gbk_enriched = 0
    if not args.skip_gbk:
        gbk_cds, gbk_enriched = enrich_genbank(gbk, hits, out_gbk)
        log.info("  GenBank enriquecido: %d/%d CDS -> %s",
                 gbk_enriched, gbk_cds, out_gbk.name)

    rep = Report(
        genome_dir=str(gdir),
        timestamp=datetime.now().isoformat(timespec="seconds"),
        db=str(args.db),
        n_cds=n_cds,
        n_hits=n_hits,
        match_rate=round(rate, 4),
        with_gene=with_gene,
        with_product=with_prod,
        with_ec=with_ec,
        mean_pident=round(mean_pid, 2),
        gbk_cds=gbk_cds,
        gbk_enriched=gbk_enriched,
        func_gbk=str(out_gbk) if not args.skip_gbk else None,
        warnings=warnings or [],
    )
    out_json.write_text(json.dumps(
        {"report": asdict(rep), "hits": {k: asdict(v) for k, v in hits.items()}},
        indent=2))
    log.info("  Escrito: %s", out_json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
