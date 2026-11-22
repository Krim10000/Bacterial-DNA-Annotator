#!/usr/bin/env python3
"""01 — Construye el panel de referencia y descomprime los FASTA.

Panel (todo acotado y documentado):

1. `local_panel`: los genomas del proyecto ya presentes en `data/raw/` y
   `data/raw_excluded/` (400 FASTA). Coste: cero descargas.
2. `conspecific`: hasta `MAX_CONSPECIFIC` ensamblados RefSeq (Complete/Chromosome)
   de la **misma especie** que cada piloto, distintos del piloto.
3. `congeneric_type`: ensamblados RefSeq de tipo ("assembly from type material" /
   "neotype"), Complete Genome, del **mismo genero** que cada piloto
   (hasta `MAX_PER_GENUS` por genero, submuestreados de forma determinista).

No se descarga ninguna base de datos grande: solo genomas individuales que
aparecen en `data/assembly_summary.txt` y solo para los generos de los pilotos.

Salidas:
  data/taxonomy/refpanel.tsv     manifiesto (una fila por genoma de referencia)
  data/taxonomy/pilot_checks.tsv contigs/bp medidos vs. los declarados
  <WORK>/refs/*.fna              FASTA listos para fastANI
"""

from __future__ import annotations

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import csv

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

MAX_CONSPECIFIC = int(os.environ.get("MAX_CONSPECIFIC", "4"))
MAX_PER_GENUS = int(os.environ.get("MAX_PER_GENUS", "200"))
WORKERS = int(os.environ.get("TAXONOMY_WORKERS", "8"))


def pick_evenly(rows: list[dict], limit: int) -> list[dict]:
    """Submuestreo determinista y repartido en todo el genero."""
    rows = sorted(rows, key=lambda r: r["assembly_accession"])
    if len(rows) <= limit:
        return rows
    if limit <= 1:
        return rows[:1]
    step = (len(rows) - 1) / (limit - 1)
    idx = sorted({int(round(i * step)) for i in range(limit)})
    return [rows[i] for i in idx]


def collect_candidates(pilots: list[dict]) -> tuple[dict[str, dict], list[dict]]:
    """Una sola pasada por assembly_summary.txt recogiendo lo relevante."""
    pilot_acc = {L.pilot_accession(p["genome"]) for p in pilots}
    genus_by_species_taxid: dict[str, str] = {}

    # 1a pasada: obtener species_taxid y genero de cada piloto
    rows_by_acc: dict[str, dict] = {}
    with open(L.ASSEMBLY_SUMMARY, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        cols = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(cols):
                continue
            if p[0] not in pilot_acc:
                continue
            row = dict(zip(cols, p))
            if row["version_status"] != "latest":
                continue
            rows_by_acc[row["assembly_accession"]] = row
            genus_by_species_taxid[row["species_taxid"]] = row["organism_name"].split()[0]

    wanted_genera = set(genus_by_species_taxid.values())
    wanted_species = set(genus_by_species_taxid)
    local_keys = [L.genome_key(f) for d in L.RAW_DIRS if d.exists() for f in d.glob("*.fna.gz")]
    local_acc = {L.pilot_accession(k) for k in local_keys}
    local_meta: dict[str, dict] = {}

    # 2a pasada: candidatos conspecificos y congenericos de tipo
    conspecific: dict[str, list[dict]] = {}
    congeneric: dict[str, list[dict]] = {}
    with open(L.ASSEMBLY_SUMMARY, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        cols = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(cols):
                continue
            row = dict(zip(cols, p))
            if row["version_status"] != "latest":
                continue
            if row["assembly_accession"] in local_acc and row["assembly_accession"] not in pilot_acc:
                local_meta[row["assembly_accession"]] = row
            if not row["assembly_accession"].startswith("GCF_"):
                continue
            if row["assembly_level"] not in ("Complete Genome", "Chromosome"):
                continue
            genus = row["organism_name"].split()[0]
            if row["species_taxid"] in wanted_species and row["assembly_accession"] not in pilot_acc:
                conspecific.setdefault(row["species_taxid"], []).append(row)
            if (
                genus in wanted_genera
                and row["relation_to_type_material"] != "na"
                and row["assembly_accession"] not in pilot_acc
            ):
                congeneric.setdefault(genus, []).append(row)

    return genus_by_species_taxid, (conspecific, congeneric, rows_by_acc, local_meta)


def main() -> int:
    L.ensure_dirs()
    pilots = L.read_pilot_tsv()
    rows = []
    log = L.log

    # --- 1. panel local ------------------------------------------------------
    local = []
    for d in L.RAW_DIRS:
        if not d.exists():
            continue
        for f in sorted(d.glob("*.fna.gz")):
            local.append(f)
    log(f"[panel] FASTA locales: {len(local)}")
    jobs = []
    for f in local:
        key = L.genome_key(f)
        dst = L.REF_DIR / f"{key}.fna"
        jobs.append((f, dst, key))
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        list(ex.map(lambda j: L.gunzip_to(j[0], j[1]), jobs))

    for f, dst, key in jobs:
        rows.append(
            {
                "genome": key,
                "role": "local_panel",
                "source": "data/" + f.parent.name,
                "organism_name": "",
                "taxid": "",
                "species_taxid": "",
                "assembly_level": "",
                "relation_to_type_material": "",
                "fasta_path": str(dst),
                "ftp_url": "",
            }
        )
    log(f"[panel] panel local listo: {len(jobs)}")

    # --- 2. metadatos de los pilotos y candidatos ---------------------------
    species_of_genus, (consp, congeneric, pilot_rows, local_meta) = collect_candidates(pilots)
    pilot_acc = set(pilot_rows)
    # enriquece el panel local con la metadata NCBI (para saber que especie es cada hit)
    for r in rows:
        m = local_meta.get(L.pilot_accession(r["genome"]))
        if m:
            r["organism_name"] = m["organism_name"]
            r["taxid"] = m["taxid"]
            r["species_taxid"] = m["species_taxid"]
            r["assembly_level"] = m["assembly_level"]
            r["relation_to_type_material"] = m["relation_to_type_material"]
    missing_pilot_meta = [L.pilot_accession(p["genome"]) for p in pilots
                          if L.pilot_accession(p["genome"]) not in pilot_rows]
    if missing_pilot_meta:
        log(f"[warn] pilotos sin fila en assembly_summary: {missing_pilot_meta}")

    ext: dict[str, dict] = {}
    for pilot in pilots:
        acc = L.pilot_accession(pilot["genome"])
        prow = pilot_rows.get(acc)
        if not prow:
            continue
        # conspecificos
        for r in pick_evenly(consp.get(prow["species_taxid"], []), MAX_CONSPECIFIC):
            ext.setdefault(r["assembly_accession"], {**r, "_role": "conspecific"})
        # congenericos de tipo
        genus = prow["organism_name"].split()[0]
        for r in pick_evenly(congeneric.get(genus, []), MAX_PER_GENUS):
            ext.setdefault(r["assembly_accession"], {**r, "_role": "congeneric_type"})

    log(f"[panel] referencias externas a descargar: {len(ext)}")
    to_dl = []
    for acc, r in sorted(ext.items()):
        url = r["ftp_path"].rstrip("/") + "/" + r["ftp_path"].rstrip("/").split("/")[-1] + "_genomic.fna.gz"
        gz = L.REF_DIR.parent / "downloads" / f"{acc}_genomic.fna.gz"
        dst = L.REF_DIR / f"{acc}.fna"
        to_dl.append((acc, r, url, gz, dst))

    def fetch(job):
        acc, r, url, gz, dst = job
        ok = L.download(url, gz)
        if not ok:
            return acc, False, url
        L.gunzip_to(gz, dst)
        return acc, True, url

    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(fetch, j) for j in to_dl]
        for fut in as_completed(futs):
            acc, ok, url = fut.result()
            done += 1
            if not ok:
                log(f"[warn] fallo descarga {acc}: {url}")
            elif done % 25 == 0:
                log(f"[panel] descargados {done}/{len(to_dl)}")

    for acc, r, url, gz, dst in to_dl:
        if not dst.exists():
            continue
        rows.append(
            {
                "genome": acc,
                "role": r["_role"],
                "source": "ncbi_refseq_ftp",
                "organism_name": r["organism_name"],
                "taxid": r["taxid"],
                "species_taxid": r["species_taxid"],
                "assembly_level": r["assembly_level"],
                "relation_to_type_material": r["relation_to_type_material"],
                "fasta_path": str(dst),
                "ftp_url": url,
            }
        )

    # --- 3. manifiesto -------------------------------------------------------
    fields = ["genome", "role", "source", "organism_name", "taxid", "species_taxid",
              "assembly_level", "relation_to_type_material", "fasta_path", "ftp_url"]
    with open(L.OUT / "refpanel.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    # --- 4. comprobacion de los FASTAS de los pilotos -----------------------
    checks = []
    for pilot in pilots:
        key = pilot["genome"]
        f = L.REF_DIR / f"{key}.fna"
        rec = {
            "genome": key,
            "declared_contigs": pilot.get("contigs", ""),
            "declared_bp": pilot.get("genome_bp", ""),
            "fasta_found": f.exists(),
        }
        if f.exists():
            st = L.fasta_stats(f)
            rec["measured_contigs"] = st["contigs"]
            rec["measured_bp"] = st["bp"]
            rec["bp_delta"] = st["bp"] - int(pilot["genome_bp"])
            rec["contigs_match"] = str(st["contigs"]) == str(pilot.get("contigs"))
            rec["bp_within_1pct"] = abs(rec["bp_delta"]) <= 0.01 * int(pilot["genome_bp"])
        checks.append(rec)
    with open(L.OUT / "pilot_checks.tsv", "w", newline="", encoding="utf-8") as fh:
        fields2 = ["genome", "declared_contigs", "measured_contigs", "declared_bp",
                   "measured_bp", "bp_delta", "contigs_match", "bp_within_1pct", "fasta_found"]
        w = csv.DictWriter(fh, fieldnames=fields2, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for c in checks:
            w.writerow(c)

    from collections import Counter

    counts = Counter(r["role"] for r in rows)
    log(f"[panel] panel final: {dict(counts)} -> {L.OUT / 'refpanel.tsv'}")
    missing = [c["genome"] for c in checks if not c.get("contigs_match") or not c.get("bp_within_1pct")]
    if missing:
        log(f"[warn] pilotos con FASTA incoherente: {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
