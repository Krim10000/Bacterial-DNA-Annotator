#!/usr/bin/env python3
"""07 — Asigna especie/linaje a los 400 genomas y escribe `full/<genome>.json`.

Mismas fuentes y mismas reglas que la capa piloto (03), pero:
  * la metadata NCBI/GTDB viene ya resuelta en `full/genome_meta.tsv`;
  * el ANI viene de `full/ani_hits_full.tsv` y cada genoma declara su **alcance**
    (`scope_rank`/`scope_name`), porque no se comparó contra todo el panel;
  * no se toca `data/taxonomy/<genome>.json` del piloto.

Salidas: data/taxonomy/full/<genome>.json (400) y full/summary_full.json
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"
SPECIES_ANI = L.ANI_SPECIES_THRESHOLD
VERSION = "2.0-full"


# ------------------------------------------------------------------ helpers
def strip_gtdb_suffix(name: str | None) -> str | None:
    return re.sub(r"_[A-Z]+$", "", name) if name else None


def is_placeholder_species(name: str | None) -> bool:
    return True if not name else bool(re.search(r"\bsp\d{5,}\b|\bsp\.", name))


def split_binomial(name: str | None) -> tuple[str | None, str | None]:
    """`Candidatus Nitrospira allomarina` -> (Nitrospira, allomarina)."""
    if not name:
        return None, None
    parts = name.split()
    if parts and parts[0].lower() == "candidatus" and len(parts) >= 3:
        parts = parts[1:]
    return (parts[0], parts[1]) if len(parts) >= 2 else (parts[0], None)


def strip_gtdb_token(name: str | None) -> str | None:
    """Quita el sufijo de desambiguacion GTDB: `parvum_A` -> `parvum`."""
    return re.sub(r"_[A-Z]$", "", name) if name else None


def epithet_key(epithet: str | None) -> str | None:
    if not epithet:
        return None
    e = epithet.lower()
    for suffix in ("ensis", "ense", "um", "us", "a", "e", "i", "o", "u"):
        if e.endswith(suffix) and len(e) > len(suffix) + 2:
            return e[: -len(suffix)]
    return e


def fmt_ani(value: float) -> float:
    v = round(value, 2)
    if v >= 100.0 and value < 100.0:
        v = 99.99
    return v


def gtdb_parts(lineage: str) -> tuple[str | None, str | None]:
    if not lineage or lineage == "unknown":
        return None, None
    genus = species = None
    for part in lineage.split(";"):
        if part.startswith("g__") and len(part) > 3:
            genus = part[3:]
        if part.startswith("s__") and len(part) > 3:
            species = part[3:]
    return genus, species


def ncbi_lineage_string(ranks: dict[str, str]) -> str:
    order = [("d__", "domain"), ("p__", "phylum"), ("c__", "class"), ("o__", "order"),
             ("f__", "family"), ("g__", "genus"), ("s__", "species")]
    return ";".join(prefix + (ranks.get(key) or "") for prefix, key in order)


# ------------------------------------------------------------------ main
def main() -> int:
    FULL.mkdir(parents=True, exist_ok=True)
    meta = []
    with open(FULL / "genome_meta.tsv", encoding="utf-8") as fh:
        meta = list(csv.DictReader(fh, delimiter="\t"))

    plan = {}
    with open(FULL / "scope_plan.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            plan[r["genome"]] = r

    panel = {}
    with open(FULL / "panel_full.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panel[r["genome"]] = r

    # numero de referencias realmente usadas por grupo (para declararlo)
    group_refs: dict[tuple[str, str], int] = {}
    gp = FULL / "ani_groups.json"
    if gp.exists():
        for g in json.loads(gp.read_text(encoding="utf-8")).get("groups", []):
            group_refs[(g.get("rank", ""), g.get("name", ""))] = g.get("n_refs", 0)

    hits: dict[str, list[dict]] = {}
    with open(FULL / "ani_hits_full.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["ani_f"] = float(r["ani"])
            hits.setdefault(r["genome"], []).append(r)
    for v in hits.values():
        v.sort(key=lambda r: -r["ani_f"])

    records = []
    for row in meta:
        genome = row["genome"]
        sc = plan.get(genome, {"scope_rank": "none", "scope_name": "", "n_refs": "0"})
        ranks = {k: row.get(k, "") for k in
                 ("domain", "phylum", "class", "order", "family", "genus", "species_ncbi")}
        species_ncbi_full = row.get("species_ncbi") or None
        g_ncbi, s_ncbi = split_binomial(species_ncbi_full)
        gt_lineage = row.get("gtdb_lineage") or "unknown"
        g_gtdb, s_gtdb_raw = gtdb_parts(gt_lineage)
        species_gtdb = None if is_placeholder_species(s_gtdb_raw) else s_gtdb_raw
        g_gtdb_base, s_gtdb = split_binomial(species_gtdb)

        # especie final
        s_ncbi_cmp = strip_gtdb_token(s_ncbi)
        s_gtdb_cmp = strip_gtdb_token(s_gtdb)
        same_epithet = bool(s_ncbi_cmp and s_gtdb_cmp and
                            (s_ncbi_cmp == s_gtdb_cmp
                             or epithet_key(s_ncbi_cmp) == epithet_key(s_gtdb_cmp)))
        if species_ncbi_full and species_gtdb:
            species = species_ncbi_full
            species_source = ("ncbi+gtdb (coinciden en epiteto)" if same_epithet
                              else "ncbi (epiteto GTDB distinto)")
        elif species_ncbi_full:
            species = species_ncbi_full
            species_source = "ncbi (sin GTDB r232)" if gt_lineage == "unknown" else "ncbi"
        elif species_gtdb:
            species = species_gtdb
            species_source = "gtdb (sin NCBI)"
        else:
            species = "unknown"
            species_source = "unknown"

        g_ncbi_cmp = strip_gtdb_token(g_ncbi)
        g_gtdb_cmp = strip_gtdb_token(g_gtdb)
        genus_renamed = bool(g_ncbi_cmp and g_gtdb_cmp and g_gtdb_cmp != g_ncbi_cmp)
        conflicts = []
        if species_ncbi_full and species_gtdb and s_ncbi_cmp != s_gtdb_cmp:
            if epithet_key(s_ncbi_cmp) == epithet_key(s_gtdb_cmp):
                conflicts.append(f"mismo epiteto con distinta terminacion: "
                                 f"NCBI='{species_ncbi_full}' vs GTDB='{species_gtdb}'")
            else:
                conflicts.append(f"NCBI='{species_ncbi_full}' vs GTDB='{species_gtdb}' "
                                 "(epiteto distinto)")
        if genus_renamed:
            conflicts.append(f"genero NCBI='{g_ncbi}' vs GTDB='{g_gtdb}'")
        if gt_lineage == "unknown":
            conflicts.append("sin linea GTDB r232 para este accession")

        # ANI
        ghits = hits.get(genome, [])
        cons, hetero = [], []
        for h in ghits:
            ref = panel.get(h["ref"], {})
            h["ref_organism"] = ref.get("organism_name") or None
            h["ref_role"] = ref.get("role", "")
            if row.get("species_taxid") and ref.get("species_taxid") and \
                    str(row["species_taxid"]) == str(ref["species_taxid"]):
                cons.append(h)
            else:
                hetero.append(h)
        best = ghits[0] if ghits else None
        bc = cons[0] if cons else None
        bh = hetero[0] if hetero else None

        ani_value = None
        if best:
            ani_value = fmt_ani(best["ani_f"])

        # el alcance declarado es el del hit que manda (puede venir de la pasada
        # dirigida por familia/conspecíficos, no del barrido por filo)
        if best and best.get("scope_rank"):
            scope_rank = best["scope_rank"]
            scope_name = best.get("scope_name", "")
        else:
            scope_rank = sc["scope_rank"]
            scope_name = sc["scope_name"]
        n_refs = group_refs.get((scope_rank, scope_name), int(sc["n_refs"] or 0))

        ani_ev = {
            "method": "fastANI",
            "fastani_version": "1.33",
            "min_fraction": float(L.FASTANI_MIN_FRACTION),
            "species_threshold": SPECIES_ANI,
            "scope_rank": scope_rank,
            "scope_name": scope_name,
            "n_refs": n_refs,
            "n_hits": len(ghits),
            "best_hit": None,
            "best_conspecific": None,
            "best_heterospecific": None,
        }
        for label, h in (("best_hit", best), ("best_conspecific", bc),
                         ("best_heterospecific", bh)):
            if h:
                ani_ev[label] = {"ref": h["ref"], "ref_organism": h.get("ref_organism"),
                                 "ref_role": h.get("ref_role"), "ani": round(h["ani_f"], 4),
                                 "orthologous_matches": int(h["orthologous_matches"]),
                                 "total_fragments": int(h["total_fragments"]),
                                 "af": float(h["af"])}

        cons_ok = bool(bc and bc["ani_f"] >= SPECIES_ANI)
        hetero_conflict = bool(bh and bh["ani_f"] >= SPECIES_ANI)
        if hetero_conflict:
            conflicts.append(f"ANI {fmt_ani(bh['ani_f'])}% >= {SPECIES_ANI}% contra "
                             f"{bh['ref']} ({bh.get('ref_organism')})")

        # confianza
        if not species_ncbi_full and gt_lineage == "unknown":
            confidence = "unknown"
        elif hetero_conflict or (species_ncbi_full and species_gtdb and not same_epithet):
            confidence = "low"
        elif (same_epithet and not genus_renamed
              and (cons_ok or (not bc and bh and bh["ani_f"] <= SPECIES_ANI - 2.0
                               and float(bh["af"]) >= 0.2))):
            confidence = "high"
        else:
            confidence = "medium"
        status = ("assigned" if confidence in ("high", "medium")
                  else "provisional" if confidence == "low" else "unknown")

        notes = []
        if not best:
            notes.append("Sin hit de ANI (alcance sin referencias o sin coincidencias): "
                         "la especie no se valida por ANI.")
        if scope_rank in ("order", "class", "phylum"):
            notes.append(f"El ANI se calculó en alcance '{scope_rank}' "
                         f"({scope_name}) con {n_refs} referencias, no a nivel de "
                         "género: es una cota inferior de la similitud disponible.")
        if best and float(best["af"]) < 0.2:
            notes.append(f"ANI {fmt_ani(best['ani_f'])}% sobre solo "
                         f"{round(float(best['af'])*100,1)}% de la query (af < 0.2): "
                         "evidencia débil.")
        if cons_ok:
            notes.append(f"ANI {fmt_ani(bc['ani_f'])}% contra {bc['ref']} "
                         f"({bc.get('ref_organism')}) >= {SPECIES_ANI}%: respalda la especie.")
        if not bc and best:
            notes.append("Sin conspecífico independiente en el alcance: la frontera de "
                         "especie no se valida por ANI.")
        if genus_renamed:
            notes.append(f"GTDB renombra el género ('{g_ncbi}' -> '{g_gtdb}'): confianza "
                         "limitada a 'medium' por nomenclatura.")
        if gt_lineage == "unknown":
            notes.append("Sin línea GTDB r232: especie anclada solo en NCBI (taxdump).")

        rec = {
            "genome": genome,
            "closest_ref": best["ref"] if best else None,
            "ani": ani_value,
            "gtdb_lineage": gt_lineage,
            "species": species,
            "confidence": confidence,
            "evidence": {
                "ani": ani_ev,
                "ncbi": {
                    "taxid": int(row["taxid"]) if row.get("taxid") else None,
                    "species_taxid": int(row["species_taxid"]) if row.get("species_taxid") else None,
                    "organism_name": row.get("organism_name") or "",
                    "assembly_level": row.get("assembly_level") or "",
                    "relation_to_type_material": row.get("relation_to_type_material") or "",
                    "lineage": ncbi_lineage_string(ranks),
                    "ranks": {k: v for k, v in ranks.items() if v},
                },
                "gtdb": {"release": "r232", "lineage": gt_lineage,
                         "species": species_gtdb, "genus": g_gtdb},
                "species_source": species_source,
                "conflicts": conflicts,
                "notes": notes,
                "annotation_dir": row.get("annotation_dir") or "",
            },
            "species_ncbi": species_ncbi_full,
            "species_gtdb": species_gtdb,
            "ncbi_lineage": ncbi_lineage_string(ranks),
            "ncbi_taxid": int(row["taxid"]) if row.get("taxid") else None,
            "status": status,
            "layer": "full",
            "taxonomy_version": {"pipeline": VERSION, "gtdb": "r232",
                                 "ncbi_taxdump": "local data/taxdump",
                                 "fastani": "1.33"},
        }
        records.append(rec)
        L.write_json(FULL / f"{genome}.json", rec)

    # ---------------- summary ----------------
    conf = Counter(r["confidence"] for r in records)
    status_c = Counter(r["status"] for r in records)
    scope_c = Counter(r["evidence"]["ani"]["scope_rank"] for r in records)
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pipeline_version": VERSION,
        "layer": "full",
        "n_genomes": len(records),
        "n_ani_computed": sum(1 for r in records if r["ani"] is not None),
        "n_ani_null": sum(1 for r in records if r["ani"] is None),
        "n_gtdb_lineage": sum(1 for r in records if r["gtdb_lineage"] != "unknown"),
        "n_species_named": sum(1 for r in records if r["species"] != "unknown"),
        "species_threshold": SPECIES_ANI,
        "confidence_counts": dict(conf),
        "status_counts": dict(status_c),
        "scope_counts": dict(scope_c),
        "n_with_conflicts": sum(1 for r in records if r["evidence"]["conflicts"]),
        "n_high_confidence_ani": sum(1 for r in records if r["confidence"] == "high"),
        "n_hits": sum(r["evidence"]["ani"]["n_hits"] for r in records),
        "sources": {
            "ncbi_assembly_summary": str(L.ASSEMBLY_SUMMARY),
            "ncbi_taxdump": str(L.TAXDUMP),
            "gtdb_taxonomy_release": "r232",
            "fastani": "1.33",
            "panel": str(FULL / "panel_full.tsv"),
            "ani_hits": str(FULL / "ani_hits_full.tsv"),
            "scope_plan": str(FULL / "scope_plan.tsv"),
        },
        "genomes": [
            {"genome": r["genome"], "species": r["species"],
             "species_ncbi": r["species_ncbi"], "species_gtdb": r["species_gtdb"],
             "closest_ref": r["closest_ref"], "ani": r["ani"],
             "confidence": r["confidence"], "status": r["status"],
             "scope_rank": r["evidence"]["ani"]["scope_rank"],
             "n_hits": r["evidence"]["ani"]["n_hits"],
             "conflicts": r["evidence"]["conflicts"]}
            for r in records
        ],
    }
    L.write_json(FULL / "summary_full.json", summary)
    L.log(f"[full] {len(records)} JSON | ANI {summary['n_ani_computed']}/"
          f"{len(records)} | GTDB {summary['n_gtdb_lineage']} | "
          f"conf {dict(conf)} | status {dict(status_c)}")
    L.log(f"[full] summary -> {FULL / 'summary_full.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
