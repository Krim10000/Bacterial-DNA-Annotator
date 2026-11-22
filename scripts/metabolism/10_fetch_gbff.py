#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10_fetch_gbff.py — descarga los GenBank PGAP (*_genomic.gbff.gz) de NCBI.

POR QUÉ: los GFF3 de RefSeq NO traen números EC (ni en el product ni en
atributos ec=/Dbxref), así que el gene-map del 09 sale vacío con ellos. Los
EC sí están en el GenBank (/EC_number por CDS) y el 09 usa el .gbff
automáticamente cuando está (prioridad gbk > bakta > gff), con
/translation incluida (secuencias sin FAA).

Fuente: data/assembly_summary.txt (RefSeq) — columna ftp_path ->
    https://ftp.ncbi.nlm.nih.gov/genomes/all/<...>/{acc}_{name}_genomic.gbff.gz
Si no tienes assembly_summary.txt, la URL se construye igualmente desde el
accesión (GCF_000007985.2_ASM798v2 -> genomes/all/GCF/000/007/985/...).
Cada descarga se verifica (gzip que empieza por LOCUS); un .gbff corrupto
o un 404/HTML se descarta y se re-descarga en la siguiente pasada.

Uso:
  # probar con un genoma:
  python scripts/metabolism/10_fetch_gbff.py --genome GCF_000006685.1_ASM668v1

  # el lote: solo los genomas que tienen red del 04 (~1-2 MB gz por genoma)
  python scripts/metabolism/10_fetch_gbff.py --net-dir data/metabolism

  # ver qué bajaría sin bajar nada:
  python scripts/metabolism/10_fetch_gbff.py --net-dir data/metabolism --dry-run

Los ficheros caen en --out-dir (data/annotations por defecto), junto a los
.gff/.faa que ya tienes. Es idempotente: salta los que ya existen salvo
--force. Después, borra los gene maps vacíos o simplemente reinicia el
portal: el 07 regenera los gene maps vacíos al vuelo con el 09.
"""
import argparse
import gzip
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
m04 = importlib.import_module("04_network_trace")
log = m04.log
PROJECT_ROOT = m04.PROJECT_ROOT

UA = {"User-Agent": "dna_token_bac-genemap/1.0 (metabolism pipeline)"}


def parse_assembly_summary(path):
    """{nombre_assembly: ftp_path} de assembly_summary.txt.

    Claves: 'GCF_000006685.1_ASM668v1' (nombre = último tramo del ftp_path,
    como llaman tus .gff/.faa/genemaps) y también el accesoión corto
    'GCF_000006685.1'. La columna ftp_path se localiza por la cabecera
    (línea '# ...') para no depender del número exacto de columnas."""
    out_acc, out_name = {}, {}
    header_idx = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                fields = [x.strip() for x in line.lstrip("# \t").split("\t")]
                if "ftp_path" in fields and header_idx is None:
                    header_idx = fields.index("ftp_path")
                continue
            cols = line.rstrip("\n").split("\t")
            idx = header_idx if header_idx is not None else 19
            if len(cols) <= idx or not cols[idx]:
                continue
            ftp = cols[idx].strip()
            name = ftp.rstrip("/").rsplit("/", 1)[-1]
            out_name[name] = ftp
            out_acc[cols[0]] = ftp
    return out_name, out_acc


def gbff_url(ftp_path):
    """URL https del *_genomic.gbff.gz de un ftp_path de assembly_summary."""
    name = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    base = ftp_path.rstrip("/")
    if base.startswith("ftp://"):
        base = "https://" + base[len("ftp://"):]
    return f"{base}/{name}_genomic.gbff.gz"


_ID_RE = re.compile(r"^(GC[AF])_(\d{9})\.\d+_.+$")


def base_url_from_id(genome_id):
    """https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/685/{id}
    construido desde el accesoión — para cuando no hay assembly_summary.txt.
    None si el id no tiene forma GC[AF]_9dígitos.ver_nombre."""
    m = _ID_RE.match(genome_id)
    if not m:
        return None
    d = m.group(2)
    return (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{m.group(1)}"
            f"/{d[0:3]}/{d[3:6]}/{d[6:9]}/{genome_id}")


def verify_gbff_gz(path):
    """True si el fichero es un GenBank descomprimible que empieza por
    LOCUS (descarta un HTML/404 guardado por error)."""
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            chunk = f.read(65536)
    except Exception:  # noqa: BLE001 — cualquier error = fichero inválido
        return False
    return any(line.startswith("LOCUS") for line in chunk.splitlines())


def net_genomes(net_dir):
    return sorted(p.name[:-len("_network.json")]
                  for p in Path(net_dir).glob("*_network.json"))


def download(url, dest, retries=3, timeout=60):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r, \
                    open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 16)
                    if not chunk:
                        break
                    f.write(chunk)
            if not verify_gbff_gz(tmp):
                log.error(f"  {url}: descargado pero NO es un GenBank "
                          "(¿404/HTML?) — no se guarda")
                tmp.unlink()
                return False
            tmp.rename(dest)
            return True
        except Exception as e:  # noqa: BLE001 — reintentar y contar
            last = e
            if attempt < retries:
                time.sleep(2 * attempt)
    log.error(f"  {url}: {last}")
    if tmp.exists():
        tmp.unlink()
    return False


def main():
    ap = argparse.ArgumentParser(
        description="Descarga {genome}_genomic.gbff.gz de NCBI (fuente de "
                    "los EC para el gene-map del 09).")
    ap.add_argument("--genome", action="append", default=[],
                    help="genoma (nombre o accesoión); repetible")
    ap.add_argument("--net-dir", type=Path, default=None,
                    help="descarga solo los genomas con {genome}_network.json")
    ap.add_argument("--summary", type=Path, default=None,
                    help="assembly_summary.txt (por defecto "
                         "data/assembly_summary.txt)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="destino (por defecto data/annotations, donde el "
                         "09/07 autodescubre)")
    ap.add_argument("--limit", type=int, default=None,
                    help="descarga como mucho N genomas (prueba a mano)")
    ap.add_argument("--force", action="store_true",
                    help="re-descarga aunque ya exista el .gbff")
    ap.add_argument("--dry-run", action="store_true",
                    help="lista las URLs sin descargar")
    args = ap.parse_args()

    summary = args.summary or (PROJECT_ROOT / "data" / "assembly_summary.txt")
    out_dir = args.out_dir or (PROJECT_ROOT / "data" / "annotations")
    by_name, by_acc = {}, {}
    if summary.exists():
        by_name, by_acc = parse_assembly_summary(summary)
        log.info(f"assembly_summary: {len(by_name)} assemblies con ftp_path")
    else:
        log.warning(f"no hay {summary}: construyo las URLs desde el "
                    "accesión (GCF_xxx.yyy_NOMBRE -> ftp.ncbi.nlm.nih.gov/"
                    "genomes/all/...). Para el mapeo exacto usa --summary")

    wanted = []
    if args.net_dir:
        wanted += net_genomes(args.net_dir)
    wanted += args.genome
    if not wanted:
        wanted = sorted(by_name)          # todo el summary

    jobs, missing = [], []
    for g in dict.fromkeys(wanted):       # únicos, conservando orden
        ftp = by_name.get(g) or by_acc.get(g)
        base = ftp.rstrip("/") if ftp else base_url_from_id(g)
        if not base:
            missing.append(g)
            continue
        if not ftp:
            log.info(f"  {g}: sin summary — URL construida desde el accesión")
        jobs.append((g, gbff_url(base), out_dir / f"{g}_genomic.gbff.gz"))
    if missing:
        log.warning(f"{len(missing)} genomas sin URL (ni summary ni accesión "
                    f"reconocible, p. ej. {missing[:3]}) — se omiten")

    if args.limit:
        jobs = jobs[:args.limit]
    # ya descargado = existe, no vacío Y es un GenBank de verdad (los
    # .part corruptos o HTML de un intento anterior se re-descargan)
    todo = [(g, u, d) for g, u, d in jobs
            if args.force or not (d.exists() and d.stat().st_size > 0
                                  and verify_gbff_gz(d))]
    skipped = len(jobs) - len(todo)
    log.info(f"genomas: {len(jobs)} ({skipped} ya descargados, "
             f"{len(todo)} por bajar) -> {out_dir}")

    if args.dry_run:
        for g, u, _ in todo:
            print(f"  {g}: {u}")
        print(f"DRY-RUN: {len(todo)} ficheros, nada descargado")
        return

    if not todo:
        log.info("nada que descargar")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    for i, (g, u, d) in enumerate(todo, 1):
        log.info(f"[{i}/{len(todo)}] {g}")
        if download(u, d):
            ok += 1
        else:
            fail += 1
    log.info(f"listo: {ok} descargados, {fail} fallos, {skipped} ya estaban")
    if fail:
        log.warning("re-lanza el mismo comando para reintentar los fallos "
                    "(es idempotente)")
    if ok:
        log.info("siguiente paso: reinicia el portal (el 07 regenera los "
                 "gene maps vacíos con el 09 usando los .gbff)")


if __name__ == "__main__":
    main()
