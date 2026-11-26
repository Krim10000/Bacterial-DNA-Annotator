#!/usr/bin/env python3
"""
ab_initio_annotator.py — Pipeline "cero a héroe" para genomas bacterianos.

Recibe un FASTA de genoma (comprimido o no) SIN anotación previa y produce:
  - <out>/annotation.gff3 : GFF3 consolidado (CDS + rRNA + tRNA)
  - <out>/proteins.faa    : proteínas predichas (solo CDS)
  - <out>/annotation.gbk  : GenBank multi-registro
  - <out>/report.json     : métricas + versiones de herramientas

Herramientas:
  - Pyrodigal   : CDS ab initio (wrapper Python de Prodigal)
  - barrnap     : rRNA (5S, 16S, 23S)          [opcional]
  - tRNAscan-SE : tRNA                          [opcional]
  - Biopython   : escritura GenBank             [opcional]

Uso básico:
    python scripts/annotation/00_annotate_single.py genome.fna.gz -o out_dir/
    python scripts/annotation/00_annotate_single.py genome.fna -o out/ --prefix GCF_xxx --mode meta
    python scripts/annotation/00_annotate_single.py --check-deps
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("ab_initio")


# ---------------------------------------------------------------------------
# 1. Verificación de dependencias
# ---------------------------------------------------------------------------

def _tool_version(cmd: str, flag: str) -> str:
    try:
        out = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=10)
        combined = (out.stdout or "") + (out.stderr or "")
        for line in combined.splitlines():
            s = line.strip()
            if s and not s.lower().startswith(("unknown", "usage", "error", "option")):
                return s
        return "(version not detected)"
    except Exception:
        return "unknown"


def check_dependencies(
    skip_rrna: bool = False, skip_trna: bool = False
) -> dict[str, str]:
    versions: dict[str, str] = {}

    try:
        import pyrodigal
        versions["pyrodigal"] = getattr(pyrodigal, "__version__", "unknown")
    except ImportError as e:
        raise RuntimeError(
            "pyrodigal no instalado. Ejecuta: pip install pyrodigal"
        ) from e

    if not skip_rrna:
        has_barrnap = shutil.which("barrnap") is not None
        has_cmsearch = shutil.which("cmsearch") is not None
        if not has_barrnap and not has_cmsearch:
            raise RuntimeError(
                "Se necesita barrnap o cmsearch (Infernal). Instálalos o usa --skip-rrna"
            )
        if has_barrnap:
            versions["barrnap"] = _tool_version("barrnap", "--version")
        if has_cmsearch:
            versions["cmsearch"] = "infernal (cmsearch)"

    if not skip_trna:
        if shutil.which("tRNAscan-SE") is None:
            raise RuntimeError(
                "tRNAscan-SE no encontrado. Instálalo o usa --skip-trna"
            )
        versions["tRNAscan-SE"] = _tool_version("tRNAscan-SE", "--version")

    try:
        import Bio
        versions["biopython"] = getattr(Bio, "__version__", "unknown")
    except ImportError:
        versions["biopython"] = "MISSING (GenBank se omitirá)"

    return versions


# ---------------------------------------------------------------------------
# 2. Lectura de FASTA
# ---------------------------------------------------------------------------

def validate_input_type(input_path: Path) -> None:
    """Reject protein FASTA; accept IUPAC nucleotide genomes."""
    name = input_path.name.lower()
    if name.endswith((".faa", ".faa.gz", ".pep", ".pep.gz", ".faa.bz2")):
        raise ValueError(
            f"Input '{input_path.name}' parece proteinas (.faa); "
            f"se necesita .fna.gz (nucleotidos)."
        )

    DNA_CHARS = set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv-*.")
    PROTEIN_ONLY = set("EFILPQZefilpqz")
    opener = gzip.open if name.endswith(".gz") else open

    buf: list[str] = []
    with opener(input_path, "rt") as fh:
        in_seq = False
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if in_seq:
                    break
                in_seq = True
                continue
            if in_seq:
                buf.append(s)
                if sum(len(x) for x in buf) >= 1000:
                    break

    first_seq = "".join(buf)
    if not first_seq:
        raise ValueError(f"Input '{input_path.name}' no tiene secuencia")

    chars = set(first_seq[:2000])
    protein_hits = chars & PROTEIN_ONLY
    dna_hits = chars & DNA_CHARS

    if len(protein_hits) >= 3 or (not dna_hits and protein_hits):
        raise ValueError(
            f"Input '{input_path.name}' parece proteinas "
            f"(caracteres exclusivos: {sorted(protein_hits)[:8]}). "
            f"Se necesita ADN (.fna.gz)."
        )


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    contigs: dict[str, str] = {}
    name: Optional[str] = None
    chunks: list[str] = []

    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            contigs[name] = "".join(chunks).upper()

    if not contigs:
        raise ValueError(f"No se encontraron secuencias en {path}")
    return contigs


def infer_prefix(input_path: Path) -> str:
    name = input_path.name
    for suffix in (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_path.stem


def ensure_uncompressed_fasta(input_path: Path, work_dir: Path) -> Path:
    """Devuelve un FASTA plano para herramientas externas que no leen gzip."""
    if not str(input_path).endswith(".gz"):
        return input_path
    out = work_dir / "genome_uncompressed.fna"
    with gzip.open(input_path, "rt") as fin, open(out, "w") as fout:
        shutil.copyfileobj(fin, fout)
    return out


# ---------------------------------------------------------------------------
# 3. Predicción de CDS (Pyrodigal)
# ---------------------------------------------------------------------------

def _prefix_attrs(attrs: str, contig: str) -> str:
    """Añade <contig>__ a ID= y Parent= para unicidad multi-contig."""
    return re.sub(
        r"\b(ID|Parent)=([^;\s]+)",
        lambda m: f"{m.group(1)}={contig}__{m.group(2)}",
        attrs,
    )


def _sanitize_for_pyrodigal(seq: str) -> bytes:
    """Convert IUPAC ambiguity codes (and anything else) to N; keep ACGT."""
    seq_up = seq.upper()
    if seq_up.isascii() is False:
        seq_up = "".join(c if ord(c) < 128 else "N" for c in seq_up)
    if all(c in "ACGT" for c in seq_up):
        return seq_up.encode("ascii")
    return "".join(c if c in "ACGT" else "N" for c in seq_up).encode("ascii")


def predict_cds(contigs: dict[str, str], mode: str = "single") -> dict:
    import pyrodigal

    valid: dict[str, str] = {}
    for name, seq in contigs.items():
        if len(seq) < 200:
            log.warning("  %s muy corto (%d bp), omitido", name, len(seq))
        else:
            valid[name] = seq
    if not valid:
        raise ValueError("No hay contigs >= 200 bp para anotar")

    meta = mode == "meta"
    try:
        finder = pyrodigal.GeneFinder(meta=meta)
    except AttributeError:
        finder = pyrodigal.OrfFinder(meta=meta)

    if not meta:
        training_seq = "N".join(valid.values())
        if len(training_seq) < 20_000:
            log.warning("  Genoma total <20 kb; cambiando a meta mode")
            try:
                finder = pyrodigal.GeneFinder(meta=True)
            except AttributeError:
                finder = pyrodigal.OrfFinder(meta=True)
            meta = True
        else:
            log.info("  Entrenando Pyrodigal en %s bp", f"{len(training_seq):,}")
            finder.train(_sanitize_for_pyrodigal(training_seq))

    gff_lines: list[str] = []
    proteins: list[tuple[str, str]] = []

    for contig_name, seq in valid.items():
        t0 = time.time()
        genes = finder.find_genes(_sanitize_for_pyrodigal(seq))

        with tempfile.NamedTemporaryFile(mode="r", suffix=".gff", delete=False) as tf:
            gff_tmp = Path(tf.name)
        with tempfile.NamedTemporaryFile(mode="r", suffix=".faa", delete=False) as tf:
            faa_tmp = Path(tf.name)

        try:
            with open(gff_tmp, "w") as dst:
                genes.write_gff(dst, sequence_id=contig_name)
            with open(faa_tmp, "w") as dst:
                genes.write_translations(dst, sequence_id=contig_name)

            n_cds = 0
            with open(gff_tmp) as fh:
                for line in fh:
                    if line.startswith("##") or not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9:
                        continue
                    cols[0] = contig_name
                    cols[8] = _prefix_attrs(cols[8], contig_name)
                    gff_lines.append("\t".join(cols))
                    if cols[2] == "CDS":
                        n_cds += 1

            cur_id = None
            cur_seq = []
            with open(faa_tmp) as fh:
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if cur_id is not None:
                            proteins.append((cur_id, "".join(cur_seq)))
                        raw = line[1:].split()[0]
                        cur_id = f"{contig_name}__{raw}"
                        cur_seq = []
                    elif line:
                        cur_seq.append(line)
                if cur_id is not None:
                    proteins.append((cur_id, "".join(cur_seq)))

            log.info("  %s: %s bp, %d CDS (%.1fs)", contig_name, f"{len(seq):,}", n_cds, time.time() - t0)
        finally:
            gff_tmp.unlink(missing_ok=True)
            faa_tmp.unlink(missing_ok=True)

    return {"gff_lines": gff_lines, "proteins": proteins}


def predict_rrna(fasta_path: Path, kingdom: str = "bac") -> list[str]:
    cmd = [
        "barrnap", "--kingdom", kingdom,
        "--threads", "2", "--quiet",
        # El E-value por defecto (1e-6) descarta copias de 5S que sí pasan el
        # corte GA de Rfam; con 10 el propio modelo decide qué es un hit.
        "--evalue", "10", str(fasta_path),
    ]
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120,
            check=True, stdin=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as e:
        log.warning("barrnap falló: %s", (e.stderr or "")[:200])
        return []
    except subprocess.TimeoutExpired:
        log.warning("barrnap timeout >30 min; se omite rRNA")
        return []

    lines = [
        line for line in out.stdout.splitlines()
        if line and not line.startswith("#") and len(line.split("\t")) >= 9
    ]
    log.info("  barrnap: %d features rRNA", len(lines))
    return lines


# ---------------------------------------------------------------------------
# 5. tRNA (tRNAscan-SE)
# ---------------------------------------------------------------------------

def predict_trna(fasta_path: Path) -> list[str]:
    workdir = Path(tempfile.mkdtemp(prefix="trna_"))
    trna_out = workdir / "trna.out"

    cmd = [
        "tRNAscan-SE", "-B", "--thread", "2", "-q",
        "-o", str(trna_out), str(fasta_path),
    ]
    try:
        subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800, check=True
        )
    except subprocess.CalledProcessError as e:
        log.warning("tRNAscan-SE falló: %s", (e.stderr or "")[:200])
        shutil.rmtree(workdir, ignore_errors=True)
        return []
    except subprocess.TimeoutExpired:
        log.warning("tRNAscan-SE timeout >30 min; se omite tRNA")
        shutil.rmtree(workdir, ignore_errors=True)
        return []

    gff_lines: list[str] = []
    try:
        with open(trna_out) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if (not line or line.startswith("---")
                        or line.startswith("Sequence")
                        or line.startswith("--------")):
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    seq_name = parts[0]
                    begin = int(parts[2])
                    end = int(parts[3])
                    aa = parts[4] or "Undet"
                except (ValueError, IndexError):
                    continue
                strand = "+" if begin < end else "-"
                start, stop = min(begin, end), max(begin, end)
                attrs = (
                    f"ID={seq_name}__tRNA_{start}_{stop};"
                    f"product=tRNA-{aa}"
                )
                gff_lines.append("\t".join([
                    seq_name, "tRNAscan-SE", "tRNA",
                    str(start), str(stop), ".", strand, ".", attrs,
                ]))
        log.info("  tRNAscan-SE: %d features tRNA", len(gff_lines))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    return gff_lines


# ---------------------------------------------------------------------------
def predict_rrna_infernal(fasta_path: Path, cm_path: Path) -> list[str]:
    """Detect rRNA with Infernal cmsearch against bacterial Rfam models."""
    workdir = Path(tempfile.mkdtemp(prefix="infernal_"))
    tbl = workdir / "out.tbl"

    cmd = [
        "cmsearch", "--cpu", "2", "--cut_ga", "--noali",
        "--tblout", str(tbl), str(cm_path), str(fasta_path),
    ]
    try:
        subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600,
            check=True, stdin=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as e:
        log.warning("cmsearch falló: %s", (e.stderr or "")[:200])
        shutil.rmtree(workdir, ignore_errors=True)
        return []
    except subprocess.TimeoutExpired:
        log.warning("cmsearch timeout >10 min; se omite rRNA")
        shutil.rmtree(workdir, ignore_errors=True)
        return []

    MODEL_TO_PRODUCT = {
        "SSU_rRNA_bacteria": "16S_rRNA",
        "LSU_rRNA_bacteria": "23S_rRNA",
        "5S_rRNA": "5S_rRNA",
    }

    gff_lines = []
    try:
        if tbl.exists():
            with open(tbl) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    c = line.split()
                    if len(c) < 10:
                        continue
                    contig = c[0]
                    model = c[2]
                    rf_acc = c[3]
                    try:
                        start = int(c[7])
                        end = int(c[8])
                    except (ValueError, IndexError):
                        continue
                    strand = c[9] if c[9] in ("+", "-") else "+"
                    if start > end:
                        start, end = end, start
                    product = MODEL_TO_PRODUCT.get(model, model)
                    attrs = (
                        f"ID={contig}__rRNA_{product}_{start}_{end};"
                        f"product={product};Dbxref=Rfam:{rf_acc}"
                    )
                    gff_lines.append("\t".join([
                        contig, "cmsearch", "rRNA",
                        str(start), str(end), ".", strand, ".", attrs,
                    ]))
        log.info("  cmsearch: %d features rRNA", len(gff_lines))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    return gff_lines


def _run_rrna(args, fasta_path: Path) -> list[str]:
    """Dispatch rRNA detection: Infernal preferred, barrnap fallback."""
    if args.rna_tool == "barrnap":
        log.info("─── [3/6] rRNA (barrnap) ───")
        return predict_rrna(fasta_path, kingdom=args.kingdom)

    cm_path = args.rfam_cm
    if cm_path is None:
        # parents[2] = raíz del proyecto (scripts/annotation/<archivo>.py).
        project_root = Path(__file__).resolve().parents[2]
        cm_path = project_root / "data" / "rfam_models" / "bac_rRNA.cm"
    cm_path = Path(cm_path)
    # cmsearch contra un CM no necesita los índices .i1*; solo cmscan los usa.
    can_infernal = (
        cm_path.exists() and shutil.which("cmsearch") is not None
    )

    if args.rna_tool == "infernal" and not can_infernal:
        raise RuntimeError(
            f"Infernal requerido pero falta cmsearch o modelos en {cm_path}"
        )

    if can_infernal:
        log.info("─── [3/6] rRNA (Infernal: %s) ───", cm_path.name)
        return predict_rrna_infernal(fasta_path, cm_path)

    log.info("─── [3/6] rRNA (barrnap, fallback) ───")
    return predict_rrna(fasta_path, kingdom=args.kingdom)


# 6. Escritura consolidada
# ---------------------------------------------------------------------------

def write_gff3(
    path: Path,
    contig_lengths: dict[str, int],
    features: list[str],
) -> None:
    header = ["##gff-version 3"]
    for name, length in contig_lengths.items():
        header.append(f"##sequence-region {name} 1 {length}")

    def key(line: str):
        c = line.split("\t")
        return (c[0], int(c[3]), int(c[4]))

    features_sorted = sorted(features, key=key)
    with open(path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write("\n".join(features_sorted) + "\n")


def write_faa(path: Path, proteins: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for pid, seq in proteins:
            fh.write(f">{pid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def write_genbank(
    path: Path, contigs: dict[str, str], features: list[str]
) -> bool:
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        log.warning("Biopython no disponible; se omite GenBank")
        return False

    per_contig: dict[str, list[str]] = {c: [] for c in contigs}
    for line in features:
        c = line.split("\t")[0]
        per_contig.setdefault(c, []).append(line)

    records = []
    for name, seq in contigs.items():
        rec = SeqRecord(
            Seq(seq), id=name, name=name[:16],
            description="ab_initio_annotator",
        )
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["organism"] = "bacterium"
        rec.annotations["data_file_division"] = "BCT"
        # Feature source (requerida por ORCA y por convencion GenBank)
        source_feature = SeqFeature(
            FeatureLocation(0, len(seq)),
            type="source",
            qualifiers={
                "organism": ["bacterium"],
                "mol_type": ["genomic DNA"],
            },
        )
        rec.features.append(source_feature)
        for line in per_contig.get(name, []):
            cols = line.split("\t")
            try:
                sf = SeqFeature(
                    FeatureLocation(
                        int(cols[3]) - 1,
                        int(cols[4]),
                        strand=-1 if cols[6] == "-" else 1,
                    ),
                    type=cols[2],
                    qualifiers=_parse_attrs(cols[8]),
                )
                rec.features.append(sf)
            except (ValueError, IndexError):
                continue
        records.append(rec)

    try:
        SeqIO.write(records, str(path), "genbank")
        log.info("  GenBank: %d registro(s)", len(records))
        return True
    except Exception as e:
        log.warning("Falló GenBank: %s", e)
        return False


def _parse_attrs(attrs: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for pair in attrs.split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            out.setdefault(k, []).append(v)
    return out


# ---------------------------------------------------------------------------
# 7. Reporte
# ---------------------------------------------------------------------------

@dataclass
class Report:
    input_path: str
    output_dir: str
    timestamp: str
    mode: str
    contigs: int
    genome_bp: int
    cds_features: int
    proteins: int
    rrna_features: int
    trna_features: int
    tool_versions: dict = field(default_factory=dict)
    files: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    elapsed_seconds: float = 0.0


# ---------------------------------------------------------------------------
# 8. CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Anotación ab initio de un genoma bacteriano crudo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", type=Path, nargs="?",
                   help="FASTA (.fna / .fasta / .fa, opcionalmente .gz)")
    p.add_argument("-o", "--output-dir", type=Path, default=None)
    p.add_argument("--prefix", default=None)
    p.add_argument("--mode", choices=["single", "meta"], default="single",
                   help="single=aislado, meta=metagenoma")
    p.add_argument("--kingdom", default="bac",
                   choices=["bac", "arc", "euk", "mito"])
    p.add_argument("--skip-rrna", action="store_true")
    p.add_argument("--rna-tool", choices=["auto", "infernal", "barrnap"],
                   default="auto",
                   help="Motor rRNA: auto=infernal->barrnap fallback")
    p.add_argument("--rfam-cm", type=Path, default=None,
                   help="Modelos Rfam (default: data/rfam_models/bac_rRNA.cm)")
    p.add_argument("--skip-trna", action="store_true")
    p.add_argument("--skip-genbank", action="store_true")
    p.add_argument("--check-deps", action="store_true",
                   help="Solo verifica dependencias y sale")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)

    if args.check_deps:
        try:
            v = check_dependencies(args.skip_rrna, args.skip_trna)
        except RuntimeError as e:
            log.error(str(e))
            return 2
        log.info("Dependencias OK:")
        for k, val in v.items():
            log.info("  %-15s %s", k, val)
        return 0

    if args.input is None:
        log.error("Falta el argumento 'input' (o usa --check-deps)")
        return 2

    input_path = args.input.resolve()
    if not input_path.is_file():
        log.error("No existe: %s", input_path)
        return 2

    try:
        validate_input_type(input_path)
    except ValueError as e:
        log.error(str(e))
        return 2

    prefix = args.prefix or infer_prefix(input_path)
    if args.output_dir:
        output_dir = args.output_dir
    else:
        project_root = Path(__file__).resolve().parent.parent
        output_dir = project_root / "data" / "annotations_ab_initio" / prefix
    output_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    warnings: list[str] = []

    log.info("═══ ab_initio_annotator ═══")
    log.info("Input  : %s", input_path)
    log.info("Output : %s", output_dir)
    log.info("Mode   : %s", args.mode)

    try:
        versions = check_dependencies(args.skip_rrna, args.skip_trna)
    except RuntimeError as e:
        log.error(str(e))
        return 2
    log.info("Tools  : %s", versions)

    log.info("─── [1/6] Genoma ───")
    contigs = read_fasta(input_path)
    total_bp = sum(len(s) for s in contigs.values())
    log.info("  %d contig(s), %s bp", len(contigs), f"{total_bp:,}")
    short = sum(1 for s in contigs.values() if len(s) < 200)
    if short:
        warnings.append(f"{short} contig(s) <200 bp")

    work_dir = output_dir / ".work"
    work_dir.mkdir(exist_ok=True)
    fasta_for_ext = ensure_uncompressed_fasta(input_path, work_dir)

    try:
        log.info("─── [2/6] CDS (Pyrodigal) ───")
        cds = predict_cds(contigs, mode=args.mode)
        cds_features = sum(
            1 for ln in cds["gff_lines"] if ln.split("\t")[2] == "CDS"
        )
        log.info("  Total CDS: %d, proteínas: %d",
                 cds_features, len(cds["proteins"]))

        rrna_lines: list[str] = []
        if not args.skip_rrna:
            rrna_lines = _run_rrna(args, fasta_for_ext)
        else:
            log.info("─── [3/6] rRNA: omitido ───")

        trna_lines: list[str] = []
        if not args.skip_trna:
            log.info("─── [4/6] tRNA (tRNAscan-SE) ───")
            trna_lines = predict_trna(fasta_for_ext)
        else:
            log.info("─── [4/6] tRNA: omitido ───")

        log.info("─── [5/6] Consolidando GFF3 ───")
        all_features = cds["gff_lines"] + rrna_lines + trna_lines
        gff_path = output_dir / "annotation.gff3"
        write_gff3(
            gff_path,
            {c: len(s) for c, s in contigs.items()},
            all_features,
        )
        log.info("  %s (%d features)", gff_path.name, len(all_features))

        faa_path = output_dir / "proteins.faa"
        write_faa(faa_path, cds["proteins"])
        log.info("  %s (%d proteínas)", faa_path.name, len(cds["proteins"]))

        gbk_ok = False
        gbk_path = output_dir / "annotation.gbk"
        if not args.skip_genbank:
            log.info("─── [6/6] GenBank ───")
            gbk_ok = write_genbank(gbk_path, contigs, all_features)

        report = Report(
            input_path=str(input_path),
            output_dir=str(output_dir),
            timestamp=datetime.now().isoformat(timespec="seconds"),
            mode=args.mode,
            contigs=len(contigs),
            genome_bp=total_bp,
            cds_features=cds_features,
            proteins=len(cds["proteins"]),
            rrna_features=len(rrna_lines),
            trna_features=len(trna_lines),
            tool_versions=versions,
            files={
                "gff3": str(gff_path),
                "faa": str(faa_path),
                "gbk": str(gbk_path) if gbk_ok else None,
            },
            warnings=warnings,
            elapsed_seconds=round(time.time() - t_start, 2),
        )
        report_path = output_dir / "report.json"
        with open(report_path, "w") as fh:
            json.dump(asdict(report), fh, indent=2)

        log.info("═══ Listo en %.1fs ═══", report.elapsed_seconds)
        log.info("  GFF3    : %s", gff_path)
        log.info("  FAA     : %s", faa_path)
        if gbk_ok:
            log.info("  GenBank : %s", gbk_path)
        log.info("  Reporte : %s", report_path)
        return 0
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
