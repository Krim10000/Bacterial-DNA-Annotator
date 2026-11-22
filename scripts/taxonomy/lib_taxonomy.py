"""Utilidades compartidas del modulo de taxonomia (ANI/GTDB).

Reglas del modulo:
- Solo se escribe en `data/taxonomy/` y en el directorio de trabajo cacheado
  (`$TAXONOMY_WORK`, por defecto `~/.cache/dna_poc/taxonomy`).
- Ningun valor taxonomico se inventa: todo sale de
  (a) `data/assembly_summary.txt` + `data/taxdump/` (NCBI),
  (b) el fichero de taxonomia GTDB descargado (release pineada),
  (c) ANI medido con fastANI.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = DATA / "taxonomy"
ASSEMBLY_SUMMARY = DATA / "assembly_summary.txt"
TAXDUMP = DATA / "taxdump"
PILOT_TSV = DATA / "mobile_typing" / "pilot_genomes.tsv"
RAW_DIRS = [DATA / "raw", DATA / "raw_excluded"]

# Directorio de trabajo pesado (FASTA descomprimidos, referencias descargadas).
WORK = Path(os.environ.get("TAXONOMY_WORK", Path.home() / ".cache" / "dna_poc" / "taxonomy"))
REF_DIR = WORK / "refs"
Q_DIR = WORK / "queries"
ANI_DIR = WORK / "ani"
GTDB_DIR = WORK / "gtdb"

GTDB_RELEASE = os.environ.get("GTDB_RELEASE", "latest")
GTDB_BASE = f"https://data.gtdb.ecogenomic.org/releases/{GTDB_RELEASE}"
FASTANI_MIN_FRACTION = "0.05"  # bajamos de 0.20 para no perder pares lejanos
ANI_SPECIES_THRESHOLD = 95.0    # umbral clasico de especie (ANI >= 95 %)


def ensure_dirs() -> None:
    for d in (OUT, REF_DIR, Q_DIR, ANI_DIR, GTDB_DIR, WORK / "bin"):
        d.mkdir(parents=True, exist_ok=True)


def fastani_bin() -> str:
    """Devuelve la ruta a fastANI: $FASTANI, bin local o PATH."""
    env = os.environ.get("FASTANI")
    if env:
        return env
    local = WORK / "bin" / "fastANI"
    if local.exists():
        return str(local)
    found = shutil.which("fastANI")
    if found:
        return found
    raise SystemExit(
        "fastANI no encontrado. Ejecuta scripts/taxonomy/00_setup_fastani.sh "
        "o define $FASTANI."
    )


def genome_key(path: Path) -> str:
    """`GCF_xxx_ASMxxx_genomic.fna[.gz]` -> `GCF_xxx_ASMxxx`."""
    name = path.name
    for suffix in (".fna.gz", ".fna", ".fa.gz", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_genomic"):
        name = name[: -len("_genomic")]
    return name


def pilot_accession(genome: str) -> str:
    """`GCF_045784305.1_DRML_24019` -> `GCF_045784305.1`."""
    parts = genome.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return genome


def read_pilot_tsv() -> list[dict]:
    with open(PILOT_TSV, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_assembly_summary(accessions: set[str] | None = None) -> dict[str, dict]:
    """Filas de assembly_summary.txt indexadas por accession (solo 'latest')."""
    wanted = accessions
    cols = None
    rows: dict[str, dict] = {}
    with open(ASSEMBLY_SUMMARY, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        cols = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(cols):
                continue
            if wanted is not None and p[0] not in wanted:
                continue
            row = dict(zip(cols, p))
            if row["version_status"] != "latest":
                continue
            rows[row["assembly_accession"]] = row
    return rows


def taxdump_lineage() -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    """(parent, rank, scientific_name) del taxdump NCBI."""
    parent: dict[int, int] = {}
    rank: dict[int, str] = {}
    with open(TAXDUMP / "nodes.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.split("\t|\t")
            if len(p) < 3:
                continue
            tid = int(p[0])
            parent[tid] = int(p[1])
            rank[tid] = p[2]
    name: dict[int, str] = {}
    with open(TAXDUMP / "names.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.split("\t|\t")
            if len(p) < 4 or p[3].rstrip("\t|\n") != "scientific name":
                continue
            name[int(p[0])] = p[1]
    return parent, rank, name


def lineage_for(taxid: int, parent, rank, name) -> list[tuple[str, str]]:
    """Lista [(rank, name)] de raiz a taxid (sin 'cellular root')."""
    out: list[tuple[str, str]] = []
    cur = taxid
    seen = 0
    while cur and cur != 1 and seen < 64:
        out.append((rank.get(cur, "no rank"), name.get(cur, f"taxid:{cur}")))
        cur = parent.get(cur, 1)
        seen += 1
    out.reverse()
    # recorta el prefijo "cellular root" que no aporta
    return [(r, n) for r, n in out if n != "cellular root"]


def ncbi_ranks(lineage: list[tuple[str, str]]) -> dict[str, str]:
    """Primer valor por rango principal (superkingdom/domain, phylum, class...)."""
    want = {
        "domain": ["domain", "superkingdom"],
        "phylum": ["phylum"],
        "class": ["class"],
        "order": ["order"],
        "family": ["family"],
        "genus": ["genus"],
        "species": ["species"],
    }
    flat: dict[str, str] = {}
    for rank, name in lineage:
        for key, accepted in want.items():
            if rank in accepted and key not in flat:
                flat[key] = name
    return flat


def gtdb_taxonomy_path() -> Path | None:
    """Fichero de taxonomia GTDB ya descargado (bac120 + ar53 combinados)."""
    combined = GTDB_DIR / "gtdb_taxonomy.tsv"
    return combined if combined.exists() else None


def load_gtdb_taxonomy() -> dict[str, str]:
    """accession (RS_GCF_/GB_GCA_) -> lineage GTDB."""
    path = gtdb_taxonomy_path()
    if path is None:
        return {}
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out[p[0]] = p[1]
    return out


def gtdb_lineage_for(accession: str, gtdb: dict[str, str]) -> str | None:
    """Busca `accession` con prefijo RS_/GB_ en la taxonomia GTDB."""
    if not gtdb:
        return None
    for prefix in ("RS_", "GB_"):
        hit = gtdb.get(prefix + accession)
        if hit:
            return hit
    return None


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False, sort_keys=False)
        fh.write("\n")


def download(url: str, dest: Path, timeout: int = 600) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = ["curl", "-sSL", "--fail", "--retry", "3", "--max-time", str(timeout),
           "-o", str(tmp), url]
    rc = subprocess.run(cmd).returncode
    if rc != 0 or not tmp.exists() or tmp.stat().st_size == 0:
        tmp.unlink(missing_ok=True)
        return False
    tmp.replace(dest)
    return True


def gunzip_to(src: Path, dst: Path) -> None:
    if dst.exists() and dst.stat().st_size > 0:
        return
    tmp = dst.with_suffix(dst.suffix + ".part")
    with gzip.open(src, "rb") as inp, open(tmp, "wb") as out:
        shutil.copyfileobj(inp, out, length=1 << 20)
    tmp.replace(dst)


def fasta_stats(path: Path) -> dict:
    """n contigs + bp, sin cargar todo en memoria."""
    n = 0
    bp = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
            else:
                bp += len(line.strip())
    return {"contigs": n, "bp": bp}


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)
