#!/usr/bin/env python3
"""fetch_gold_k12.py — gold de complejos de *E. coli* K-12 y caché de STRING.

Qué descarga
------------
- ComplexPortal, tabla de la especie 83333 (K-12): complejos curados con
  subunidades y estequiometría.
- STRING v12.0, `protein.links.full` de 511145 (K-12 MG1655): seis canales por
  par, incluido `experiments` (evidencia experimental) y `database` (curada).

Cómo mapea
----------
ComplexPortal identifica subunidades por acceso UniProt; el pipeline trabaja con
locus_tag (b-number). El GenBank de K-12 trae `db_xref="UniProtKB/Swiss-Prot:..."`,
así que el mapeo sale del propio genoma, sin tabla externa.

Salidas (bajo `data/complexes/`)
--------------------------------
- `k12_genes.tsv`     genes del GenBank con UniProt (para features y mapeo)
- `gold_k12.tsv`      complejos con >=2 genes distintos y mapeo completo
- `raw/`              descargas originales
- `manifest.json`     URLs, fecha y conteos (provenance)

Uso
---
    python scripts/complexes/fetch_gold_k12.py
    python scripts/complexes/fetch_gold_k12.py --force
    python scripts/complexes/fetch_gold_k12.py --gbff data/test_ecoli.gbff
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import complexes_lib as L  # noqa: E402

CP_URL = ("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/"
          "complextab/83333.tsv")
STRING_URL = ("https://stringdb-downloads.org/download/protein.links.full.v12.0/"
              "511145.protein.links.full.v12.0.txt.gz")
SPECIES = "511145"
RAW = L.CX / "raw"


def download(url: str, dest: Path, force: bool = False, tries: int = 3) -> Path:
    if dest.exists() and not force:
        print(f"  ya existe: {dest.relative_to(L.ROOT)}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, tries + 1):
        try:
            print(f"  bajando {url}")
            with urllib.request.urlopen(url, timeout=180) as fh, open(tmp, "wb") as out:
                while True:
                    chunk = fh.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(dest)
            print(f"  ok: {dest.stat().st_size / 1e6:.1f} MB")
            return dest
        except Exception as exc:                       # red, 5xx, corte
            print(f"  intento {attempt}/{tries} falló: {exc}")
            if attempt < tries:
                time.sleep(3 * attempt)
    raise SystemExit(f"no se pudo descargar {url}")


def build_gold(gbff: Path, cp_tsv: Path, out: Path) -> tuple[dict, dict]:
    genes = L.read_genes_gbff(gbff)
    if not genes:
        raise SystemExit(f"sin CDS en {gbff}")
    L.write_genes_tsv(genes, L.CX / "k12_genes.tsv")
    acc2loc = {g.uniprot: g.locus_tag for g in genes.values() if g.uniprot}
    loc2acc = {g.locus_tag: g.uniprot for g in genes.values() if g.uniprot}
    print(f"  genes: {len(genes)}  con UniProt: {len(acc2loc)}")

    complexes, stats = L.read_complexportal(cp_tsv, acc2loc, drop_incomplete=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["complex_id", "complex_name", "subunits", "n_genes",
                    "locus_tags", "gene_names", "evidence_code",
                    "evidence_class", "properties", "source"])
        for c in sorted(complexes, key=lambda x: (-len(x["locus_tags"]), x["complex_id"])):
            names = [genes[t].gene_name or t for t in c["locus_tags"]]
            w.writerow([c["complex_id"], c["complex_name"], c["subunits"],
                        len(c["locus_tags"]), ",".join(c["locus_tags"]),
                        ",".join(names), c["evidence_code"],
                        c["evidence_class"], c["properties"], c["source"]])
    print(f"  ComplexPortal: {stats['total']} complejos, {stats['kept']} usables "
          f"(sin mapear: {stats['unmapped']}, homómeros: {stats['homomers']})")
    return stats, {"n_genes": len(genes), "n_with_uniprot": len(acc2loc)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gbff", type=Path, default=L.DATA / "test_ecoli.gbff")
    ap.add_argument("--out-dir", type=Path, default=L.CX)
    ap.add_argument("--force", action="store_true", help="re-descarga aunque exista")
    ap.add_argument("--no-download", action="store_true",
                    help="usa solo lo que ya esté en data/complexes/raw")
    a = ap.parse_args(argv)

    print("1) descargas")
    cp = RAW / "complexportal_83333.tsv"
    st = RAW / "string_511145_full.txt.gz"
    if a.no_download:
        for p in (cp, st):
            if not p.exists():
                raise SystemExit(f"falta {p} y --no-download está activo")
    else:
        download(CP_URL, cp, a.force)
        download(STRING_URL, st, a.force)

    print("2) genes y gold")
    if not a.gbff.exists():
        raise SystemExit(
            f"falta {a.gbff}. Descárgalo de NCBI (NC_000913.3) o pasa --gbff.")
    stats, gene_stats = build_gold(a.gbff, cp, a.out_dir / "gold_k12.tsv")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "species_string": SPECIES,
        "sources": {"complexportal": CP_URL, "string": STRING_URL,
                    "genbank": str(a.gbff)},
        "complexportal": stats,
        "genes": gene_stats,
    }
    (a.out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"3) listo: {a.out_dir / 'gold_k12.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
