#!/usr/bin/env python3
"""05 — Prepara la capa taxonómica de los 400 genomas (metadata, no ANI).

Universo: los 400 genomas de `data/annotations_ab_initio/` (que coinciden 1:1 con
los 400 FASTA de `data/raw/` + `data/raw_excluded/` ya cacheados).

Reutiliza todo lo del piloto:
  * `data/taxonomy/refpanel.tsv` (695 genomas: 400 base + 285 material tipo + 10
    conspecíficos) y sus FASTA en `$TAXONOMY_WORK/refs/`.
  * `data/assembly_summary.txt` y `data/taxdump/` (linaje NCBI).
  * taxonomía GTDB r232 ya cacheada (por accession).

Salidas en `data/taxonomy/full/`:
  genome_meta.tsv   una fila por genoma base (linaje NCBI + GTDB + phylum...)
  panel_full.tsv    panel de referencia etiquetado por rango (filo/clase/orden/familia)
  scope_plan.tsv    qué referencias se usarán para cada query y por qué rango
"""

from __future__ import annotations

import csv

import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"
RANKS = ["family", "order", "class", "phylum"]  # filo solo como último recurso
MIN_REFS = 1  # basta una referencia del mismo rango (el alcance se declara)


def base_genomes() -> list[str]:
    root = L.DATA / "annotations_ab_initio"
    return sorted(p.name for p in root.iterdir() if p.is_dir())


def fasta_for(genome: str) -> Path | None:
    p = L.REF_DIR / f"{genome}.fna"
    return p if p.exists() else None


def main() -> int:
    FULL.mkdir(parents=True, exist_ok=True)
    genomes = base_genomes()
    L.log(f"[full] genomas base: {len(genomes)}")
    acc_of = {g: L.pilot_accession(g) for g in genomes}

    # --- metadata NCBI de los 400 -------------------------------------------
    meta = L.load_assembly_summary(set(acc_of.values()))
    L.log(f"[full] filas en assembly_summary: {len(meta)}")

    # --- linaje NCBI (taxdump) ----------------------------------------------
    L.log("[full] leyendo taxdump")
    parent, rank, names = L.taxdump_lineage()

    def ranks_of(taxid: int) -> dict[str, str]:
        lin = L.lineage_for(taxid, parent, rank, names)
        out: dict[str, str] = {}
        for r, n in lin:
            if r in ("domain", "superkingdom"):
                out.setdefault("domain", n)
            elif r in ("phylum", "class", "order", "family", "genus", "species"):
                out.setdefault(r, n)
        return out

    # --- GTDB ---------------------------------------------------------------
    gtdb = L.load_gtdb_taxonomy()
    L.log(f"[full] GTDB accesiones: {len(gtdb)}")

    # --- panel --------------------------------------------------------------
    panel: list[dict] = []
    with open(L.OUT / "refpanel.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panel.append(r)

    panel_ranks: dict[str, dict[str, str]] = {}
    for r in panel:
        tid = r.get("taxid") or ""
        if not tid and r["role"] == "local_panel":
            # los 10 pilotos quedaron sin metadata en refpanel; se resuelve aquí
            m = meta.get(L.pilot_accession(r["genome"]))
            tid = m["taxid"] if m else ""
        if tid:
            panel_ranks[r["genome"]] = ranks_of(int(tid))
        else:
            panel_ranks[r["genome"]] = {}
    L.log("[full] rangos del panel resueltos: "
          f"{sum(1 for v in panel_ranks.values() if v)}/{len(panel)}")

    # --- filas de metadata por genoma base ----------------------------------
    rows = []
    missing_fasta, missing_meta, missing_gtdb = [], [], []
    for g in genomes:
        fa = fasta_for(g)
        if fa is None:
            missing_fasta.append(g)
        m = meta.get(acc_of[g])
        if not m:
            missing_meta.append(g)
        rk = ranks_of(int(m["taxid"])) if m else {}
        gl = L.gtdb_lineage_for(acc_of[g], gtdb)
        if not gl:
            missing_gtdb.append(g)
        gs = None
        if gl:
            for part in gl.split(";"):
                if part.startswith("s__") and len(part) > 3:
                    gs = part[3:]
                    break
        rows.append(
            {
                "genome": g,
                "accession": acc_of[g],
                "taxid": m["taxid"] if m else "",
                "organism_name": m["organism_name"] if m else "",
                "species_taxid": m["species_taxid"] if m else "",
                "assembly_level": m["assembly_level"] if m else "",
                "relation_to_type_material": m["relation_to_type_material"] if m else "",
                "domain": rk.get("domain", ""),
                "phylum": rk.get("phylum", ""),
                "class": rk.get("class", ""),
                "order": rk.get("order", ""),
                "family": rk.get("family", ""),
                "genus": rk.get("genus", ""),
                "species_ncbi": rk.get("species", ""),
                "gtdb_lineage": gl or "unknown",
                "species_gtdb": gs or "",
                "fasta_path": str(fa) if fa else "",
                "annotation_dir": str(L.DATA / "annotations_ab_initio" / g),
            }
        )

    fields = list(rows[0].keys())
    with open(FULL / "genome_meta.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # --- panel etiquetado ---------------------------------------------------
    for r in panel:
        rk = panel_ranks.get(r["genome"], {})
        for k in ("domain", "phylum", "class", "order", "family"):
            r[f"rk_{k}"] = rk.get(k, "")
    pfields = list(panel[0].keys())
    with open(FULL / "panel_full.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=pfields, delimiter="\t")
        w.writeheader()
        w.writerows(panel)

    # --- plan de alcance por query -----------------------------------------
    by_rank: dict[str, dict[str, list[str]]] = {
        rk: defaultdict(list) for rk in RANKS
    }
    panel_ids = {r["genome"] for r in panel}
    for r in panel:
        for rk in RANKS:
            val = r[f"rk_{rk}"]
            if val:
                by_rank[rk][val].append(r["genome"])

    plan = []
    scope_counter = Counter()
    pair_total = 0
    for row in rows:
        g = row["genome"]
        chosen = None
        for rk in RANKS:
            val = row.get(rk, "")
            if not val:
                continue
            refs = [x for x in by_rank[rk].get(val, []) if x != g]
            if len(refs) >= MIN_REFS:
                chosen = (rk, val, refs)
                break
        if chosen is None:
            plan.append({"genome": g, "scope_rank": "none", "scope_name": "",
                         "n_refs": 0, "refs": ""})
            scope_counter["none"] += 1
            continue
        rk, val, refs = chosen
        plan.append({"genome": g, "scope_rank": rk, "scope_name": val,
                     "n_refs": len(refs), "refs": ",".join(sorted(refs))})
        scope_counter[rk] += 1
        pair_total += len(refs)

    with open(FULL / "scope_plan.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "scope_rank", "scope_name",
                                           "n_refs", "refs"], delimiter="\t")
        w.writeheader()
        w.writerows(plan)

    # --- informe de preparación --------------------------------------------
    phyla = Counter(r["phylum"] for r in rows)
    fams = Counter(r["family"] for r in rows)
    info = {
        "n_genomes": len(genomes),
        "n_panel": len(panel),
        "missing_fasta": missing_fasta,
        "missing_assembly_summary": missing_meta,
        "missing_gtdb": missing_gtdb,
        "phyla": dict(phyla),
        "n_phyla": len(phyla),
        "top_families": fams.most_common(10),
        "scope_counts": dict(scope_counter),
        "estimated_pairs": pair_total,
        "min_refs": MIN_REFS,
    }
    L.write_json(FULL / "prep_stats.json", info)
    print(json.dumps(info, indent=1, ensure_ascii=False)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
