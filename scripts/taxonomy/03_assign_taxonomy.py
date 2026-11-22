#!/usr/bin/env python3
"""03 — Asigna especie y linaje y escribe un JSON por genoma.

Fuentes (ninguna se inventa):

* **NCBI**: `data/assembly_summary.txt` (metadata del ensamblado) + `data/taxdump/`
  (linaje real por taxid, incluidos los rangos kingdom/phylum que NCBI ya usa
  con la nomenclatura de GTDB: Bacillati, Pseudomonadati, ...).
* **GTDB**: fichero de taxonomia GTDB de la release pineada
  (`bac120_taxonomy.tsv.gz` + `ar53_taxonomy.tsv.gz`, ~10 MB en total), buscado
  por accession (`RS_GCF_...`). No se usa GTDB-Tk (no esta instalado) ni su
  base de datos completa (~100 GB).
* **ANI**: `data/taxonomy/ani_hits.tsv` generado por fastANI en el paso 02.

El ANI nunca se estima a partir de la taxonomia: si no hubo hit, `ani` es null
y el linaje se queda con la evidencia puramente meta-datos.

Salidas:
  data/taxonomy/<genome>.json   (10)
  data/taxonomy/summary.json
"""

from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

SPECIES_ANI = L.ANI_SPECIES_THRESHOLD
VERSION = "1.0"


# --------------------------------------------------------------------------- GTDB
def ensure_gtdb() -> tuple[str, dict[str, str]]:
    """Descarga (una vez) y carga la taxonomia GTDB. Devuelve (release, mapa)."""
    combined = L.GTDB_DIR / "gtdb_taxonomy.tsv"
    release_file = L.GTDB_DIR / "release.txt"
    if combined.exists() and release_file.exists():
        release = release_file.read_text(encoding="utf-8").strip()
    else:
        L.GTDB_DIR.mkdir(parents=True, exist_ok=True)
        release = "unknown"
        for domain, name in (("bac120", "bac120_taxonomy.tsv.gz"), ("ar53", "ar53_taxonomy.tsv.gz")):
            gz = L.GTDB_DIR / name
            if not gz.exists():
                url = f"{L.GTDB_BASE}/{name}"
                L.log(f"[gtdb] descargando {url}")
                if not L.download(url, gz, timeout=900):
                    L.log(f"[gtdb] no se pudo descargar {url}")
            if gz.exists():
                with gzip.open(gz, "rt", encoding="utf-8", errors="replace") as inp, \
                        open(combined.with_suffix(".part"), "a", encoding="utf-8") as out:
                    shutil.copyfileobj(inp, out, length=1 << 20)
        if combined.with_suffix(".part").exists():
            combined.with_suffix(".part").replace(combined)
        release = detect_gtdb_release()
        release_file.write_text(release + "\n", encoding="utf-8")

    gtdb: dict[str, str] = {}
    if combined.exists():
        with open(combined, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    gtdb[p[0]] = p[1]
    return release, gtdb


def detect_gtdb_release() -> str:
    """Lee la release real de la pagina de GTDB (sin inventarla)."""
    try:
        import subprocess

        html = subprocess.run(
            ["curl", "-sSL", "--max-time", "60", "https://data.gtdb.ecogenomic.org/releases/"],
            capture_output=True, text=True, timeout=90,
        ).stdout
        hits = sorted(set(re.findall(r"release(\d{3})/", html)))
        if hits:
            return "r" + hits[-1]
    except Exception:
        pass
    return "unknown"


# ------------------------------------------------------------------------- helpers
def gtdb_species(lineage: str) -> str | None:
    if not lineage:
        return None
    for part in lineage.split(";"):
        if part.startswith("s__") and len(part) > 3:
            return part[3:]
    return None


def gtdb_genus(lineage: str) -> str | None:
    if not lineage:
        return None
    for part in lineage.split(";"):
        if part.startswith("g__") and len(part) > 3:
            return part[3:]
    return None


def strip_gtdb_suffix(name: str | None) -> str | None:
    """`Treponema_D` -> `Treponema`, `Escherichia_A` -> `Escherichia`."""
    if not name:
        return None
    return re.sub(r"_[A-Z]+$", "", name)


def is_placeholder_species(name: str | None) -> bool:
    if not name:
        return True
    return bool(re.search(r"\bsp\d{5,}\b|\bsp\.", name))


def split_binomial(name: str | None) -> tuple[str | None, str | None]:
    if not name:
        return None, None
    parts = name.split()
    if len(parts) < 2:
        return parts[0] if parts else None, None
    return parts[0], parts[1]


def fmt_ani(value: float) -> float:
    """Redondea a 2 decimales sin convertir 99.9999 en 100.0."""
    v = round(value, 2)
    if v >= 100.0 and value < 100.0:
        v = 99.99
    return v


def epithet_key(epithet: str | None) -> str | None:
    """Normaliza la terminacion latina del epiteto para poder comparar."""
    if not epithet:
        return None
    e = epithet.lower()
    for suffix in ("ensis", "ense", "ensis", "um", "us", "us", "a", "e", "i", "o", "u"):
        if e.endswith(suffix) and len(e) > len(suffix) + 2:
            return e[: -len(suffix)]
    return e


def ncbi_lineage_string(lineage: list[tuple[str, str]]) -> str:
    """Formato tipo GTDB con los rangos NCBI: d__;p__;c__;o__;f__;g__;s__."""
    mapping = [
        ("domain", ("domain", "superkingdom"), "d__"),
        ("phylum", ("phylum",), "p__"),
        ("class", ("class",), "c__"),
        ("order", ("order",), "o__"),
        ("family", ("family",), "f__"),
        ("genus", ("genus",), "g__"),
        ("species", ("species",), "s__"),
    ]
    parts: list[str] = []
    for _label, ranks, prefix in mapping:
        val = None
        for rank, name in lineage:
            if rank in ranks:
                val = name
                break
        parts.append(prefix + (val if val else ""))
    return ";".join(parts)


# --------------------------------------------------------------------------- main
def main() -> int:
    L.ensure_dirs()
    pilots = L.read_pilot_tsv()
    pilot_keys = [p["genome"] for p in pilots]
    pilot_acc = {p["genome"]: L.pilot_accession(p["genome"]) for p in pilots}

    # panel y ANI
    panel: dict[str, dict] = {}
    with open(L.OUT / "refpanel.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panel[r["genome"]] = r

    hits: dict[str, list[dict]] = {}
    hits_path = L.OUT / "ani_hits.tsv"
    if hits_path.exists():
        with open(hits_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["ani_f"] = float(r["ani"])
                hits.setdefault(r["genome"], []).append(r)
        for v in hits.values():
            v.sort(key=lambda r: -r["ani_f"])

    # metadatos NCBI
    meta = L.load_assembly_summary(set(pilot_acc.values()))
    L.log(f"[ncbi] filas de assembly_summary para los pilotos: {len(meta)}")

    L.log("[ncbi] leyendo taxdump (puede tardar ~2 min)")
    parent, rank, names = L.taxdump_lineage()

    release, gtdb = ensure_gtdb()
    L.log(f"[gtdb] release={release} accesiones={len(gtdb)}")

    records = []
    for pilot in pilots:
        genome = pilot["genome"]
        acc = pilot_acc[genome]
        row = meta.get(acc, {})

        ncbi_info = {"found": bool(row)}
        species_ncbi = None
        lineage_list: list[tuple[str, str]] = []
        if row:
            taxid = int(row["taxid"])
            lineage_list = L.lineage_for(taxid, parent, rank, names)
            _, species_ncbi = split_binomial(
                next((n for r, n in lineage_list if r == "species"), None)
            )
            species_ncbi_full = next((n for r, n in lineage_list if r == "species"), None)
            ncbi_info.update(
                {
                    "taxid": taxid,
                    "species_taxid": int(row["species_taxid"]),
                    "organism_name": row["organism_name"],
                    "assembly_level": row["assembly_level"],
                    "relation_to_type_material": row["relation_to_type_material"],
                    "refseq_category": row["refseq_category"],
                    "release_date": row["seq_rel_date"],
                    "lineage": ncbi_lineage_string(lineage_list),
                    "lineage_ranks": [
                        {"rank": r, "name": n} for r, n in lineage_list if r in
                        ("domain", "superkingdom", "phylum", "class", "order", "family",
                         "genus", "species", "strain")
                    ],
                }
            )
        else:
            species_ncbi_full = None
            ncbi_info["lineage"] = "unknown"

        # GTDB
        gt_lineage = L.gtdb_lineage_for(acc, gtdb)
        gs = gtdb_species(gt_lineage) if gt_lineage else None
        gg = gtdb_genus(gt_lineage) if gt_lineage else None
        gtdb_info = {
            "release": release,
            "accession": f"RS_{acc}" if gt_lineage else None,
            "lineage": gt_lineage or "unknown",
            "species": gs,
            "genus": gg,
        }

        # especie final y consenso
        species_gtdb = None if is_placeholder_species(gs) else gs
        g_ncbi_genus, s_ncbi = split_binomial(species_ncbi_full)
        g_gtdb_genus, s_gtdb = split_binomial(species_gtdb)

        species = None
        species_source = "unknown"
        if species_ncbi_full and species_gtdb:
            species = species_ncbi_full
            species_source = "ncbi+gtdb (coinciden en epiteto)" if s_ncbi == s_gtdb else "ncbi (conflicto con gtdb)"
        elif species_ncbi_full:
            species = species_ncbi_full
            species_source = "ncbi"
        elif species_gtdb:
            species = species_gtdb
            species_source = "gtdb"
        elif row:
            species = "unknown"
        else:
            species = "unknown"

        conflicts = []
        if species_ncbi_full and species_gtdb and s_ncbi != s_gtdb:
            if epithet_key(s_ncbi) == epithet_key(s_gtdb):
                conflicts.append(
                    f"mismo epiteto con distinta terminacion latina: "
                    f"NCBI='{species_ncbi_full}' vs GTDB='{species_gtdb}'"
                )
            else:
                conflicts.append(
                    f"NCBI='{species_ncbi_full}' vs GTDB='{species_gtdb}' (epiteto distinto)"
                )
        if species_ncbi_full and g_gtdb_genus and g_ncbi_genus:
            if strip_gtdb_suffix(g_gtdb_genus) != strip_gtdb_suffix(g_ncbi_genus):
                conflicts.append(
                    f"genero NCBI='{g_ncbi_genus}' vs GTDB='{g_gtdb_genus}'"
                )
            elif g_gtdb_genus != g_ncbi_genus:
                conflicts.append(
                    f"genero renombrado en GTDB: NCBI='{g_ncbi_genus}' -> GTDB='{g_gtdb_genus}' (mismo genero base)"
                )

        # ANI
        genome_hits = hits.get(genome, [])
        conspecific_hits, hetero_hits = [], []
        for h in genome_hits:
            ref_meta = panel.get(h["ref"], {})
            ref_taxid = ref_meta.get("species_taxid") or ""
            if row and ref_meta.get("role") == "local_panel" and ref_meta.get("species_taxid") == "":
                ref_species = None
            else:
                ref_species = ref_meta.get("organism_name") or None
            h["ref_organism"] = ref_species
            h["ref_role"] = ref_meta.get("role", "")
            # conspecifico SOLO si comparte el species_taxid de ESTE piloto
            if row and ref_taxid and str(row["species_taxid"]) == str(ref_taxid):
                conspecific_hits.append(h)
            else:
                hetero_hits.append(h)

        best = genome_hits[0] if genome_hits else None
        best_cons = conspecific_hits[0] if conspecific_hits else None
        best_hetero = hetero_hits[0] if hetero_hits else None

        ani_value = None
        if best:
            ani_value = round(best["ani_f"], 2)
            # no redondear a 100.00 algo que no es identico
            if ani_value >= 100.0 and best["ani_f"] < 100.0:
                ani_value = 99.99
        closest_ref = best["ref"] if best else None

        ani_evidence = {
            "method": "fastANI",
            "fastani_version": "1.33",
            "min_fraction": float(L.FASTANI_MIN_FRACTION),
            "species_threshold": SPECIES_ANI,
            "n_hits": len(genome_hits),
            "n_refs_searched": sum(1 for g in panel if g != genome),
            "best_hit": None,
            "best_conspecific": None,
            "best_heterospecific": None,
        }
        for label, h in (("best_hit", best), ("best_conspecific", best_cons),
                         ("best_heterospecific", best_hetero)):
            if h:
                ani_evidence[label] = {
                    "ref": h["ref"],
                    "ref_organism": h.get("ref_organism"),
                    "ref_role": h.get("ref_role"),
                    "ani": round(h["ani_f"], 4),
                    "orthologous_matches": int(h["orthologous_matches"]),
                    "total_fragments": int(h["total_fragments"]),
                    "af": float(h["af"]),
                }

        # criterios de confianza
        cons_ok = bool(best_cons and best_cons["ani_f"] >= SPECIES_ANI)
        hetero_conflict = bool(best_hetero and best_hetero["ani_f"] >= SPECIES_ANI)
        if hetero_conflict:
            conflicts.append(
                f"ANI {fmt_ani(best_hetero['ani_f'])}% >= {SPECIES_ANI}% contra "
                f"{best_hetero['ref']} ({best_hetero.get('ref_organism')})"
            )

        has_ncbi = bool(species_ncbi_full)
        has_gtdb = bool(gt_lineage)
        same_epithet = bool(
            s_ncbi and s_gtdb
            and (s_ncbi == s_gtdb or epithet_key(s_ncbi) == epithet_key(s_gtdb))
        )
        genus_renamed = bool(
            g_ncbi_genus and g_gtdb_genus and g_gtdb_genus != g_ncbi_genus
        )

        if not has_ncbi and not has_gtdb:
            confidence = "unknown"
        elif hetero_conflict or (has_ncbi and has_gtdb and not same_epithet):
            confidence = "low"
        elif (same_epithet and not genus_renamed
              and (cons_ok or (not best_cons and best_hetero
                               and best_hetero["ani_f"] <= SPECIES_ANI - 2.0
                               and float(best_hetero["af"]) >= 0.2))):
            confidence = "high"
        elif same_epithet:
            confidence = "medium"
        elif has_ncbi and not has_gtdb:
            confidence = "medium"
        else:
            confidence = "medium"

        status = "assigned" if confidence in ("high", "medium") else (
            "provisional" if confidence == "low" else "unknown")

        notes = []
        if not best:
            notes.append(
                "Sin hit de ANI por encima de minFraction="
                f"{L.FASTANI_MIN_FRACTION}: no hay referencia cercana en el panel; "
                "la especie NO se valida por ANI.")
        if best and not best_cons:
            notes.append(
                "No hay genoma conspecifico independiente en el panel: el ANI no "
                "puede confirmar la frontera de especie (posible especie singlete).")
        if cons_ok:
            notes.append(
                f"ANI {fmt_ani(best_cons['ani_f'])}% contra {best_cons['ref']} "
                f"({best_cons.get('ref_organism')}) >= {SPECIES_ANI}%: respalda la especie.")
        if row.get("relation_to_type_material") not in (None, "", "na"):
            notes.append(
                f"El propio ensamblado es material tipo ({row['relation_to_type_material']}): "
                "el nombre de especie es el ancla nomenclatural.")
        if genus_renamed:
            notes.append(
                f"GTDB renombra el genero ('{g_ncbi_genus}' -> '{g_gtdb_genus}'): "
                "confianza limitada a 'medium' por nomenclatura, no por ANI.")
        if best and float(best["af"]) < 0.2:
            notes.append(
                f"ANI {fmt_ani(best['ani_f'])}% calculado sobre solo "
                f"{round(float(best['af'])*100, 1)}% de los fragmentos de la query "
                "(af < 0.2): evidencia de ANI debil, sirve como cota de genero/familia, "
                "no como ANI de especie.")
        if best and best["ani_f"] >= 99.99:
            notes.append(
                f"ANI {round(best['ani_f'], 4)}% contra {best['ref']}: genomas practicamente "
                "indistinguibles (posible reenvio/duplicado del mismo aislado).")
        if species_gtdb is None and gt_lineage:
            notes.append("GTDB no asigna especie nombrada a este accession (placeholder o sin s__).")

        record = {
            "genome": genome,
            "closest_ref": closest_ref,
            "ani": ani_value,
            "gtdb_lineage": gt_lineage if gt_lineage else "unknown",
            "species": species,
            "confidence": confidence,
            "evidence": {
                "ani": ani_evidence,
                "ncbi": ncbi_info,
                "gtdb": gtdb_info,
                "species_source": species_source,
                "conflicts": conflicts,
                "notes": notes,
                "checks": {
                    "pilot_phylum_declared": pilot.get("phylum"),
                    "contigs_declared": pilot.get("contigs"),
                    "genome_bp_declared": pilot.get("genome_bp"),
                },
            },
            "species_ncbi": species_ncbi_full,
            "species_gtdb": species_gtdb,
            "ncbi_lineage": ncbi_info.get("lineage", "unknown"),
            "ncbi_taxid": ncbi_info.get("taxid"),
            "status": status,
            "taxonomy_version": {"pipeline": VERSION, "gtdb": release,
                                 "ncbi_taxdump": "local data/taxdump",
                                 "fastani": "1.33"},
        }
        records.append(record)
        L.write_json(L.OUT / f"{genome}.json", record)
        L.log(f"[json] {genome}: {species} | ani={ani_value} | {confidence}")

    # summary
    from collections import Counter

    conf = Counter(r["confidence"] for r in records)
    species_counter = Counter(r["species"] for r in records)
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pipeline_version": VERSION,
        "n_genomes": len(records),
        "n_ani_computed": sum(1 for r in records if r["ani"] is not None),
        "n_ani_null": sum(1 for r in records if r["ani"] is None),
        "n_gtdb_lineage": sum(1 for r in records if r["gtdb_lineage"] != "unknown"),
        "species_threshold": SPECIES_ANI,
        "confidence_counts": dict(conf),
        "status_counts": dict(Counter(r["status"] for r in records)),
        "n_with_conflicts": sum(1 for r in records if r["evidence"]["conflicts"]),
        "species_counts": dict(species_counter),
        "sources": {
            "ncbi_assembly_summary": str(L.ASSEMBLY_SUMMARY),
            "ncbi_taxdump": str(L.TAXDUMP),
            "gtdb_taxonomy_release": release,
            "gtdb_taxonomy_cache": str(L.GTDB_DIR / "gtdb_taxonomy.tsv"),
            "ani": "fastANI 1.33 (binario estatico, release GitHub v1.34)",
            "ani_hits": str(hits_path),
            "refpanel": str(L.OUT / "refpanel.tsv"),
        },
        "genomes": [
            {
                "genome": r["genome"],
                "species": r["species"],
                "species_ncbi": r["species_ncbi"],
                "species_gtdb": r["species_gtdb"],
                "closest_ref": r["closest_ref"],
                "ani": r["ani"],
                "confidence": r["confidence"],
                "status": r["status"],
                "conflicts": r["evidence"]["conflicts"],
            }
            for r in records
        ],
    }
    L.write_json(L.OUT / "summary.json", summary)
    L.log(f"[json] summary -> {L.OUT / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
