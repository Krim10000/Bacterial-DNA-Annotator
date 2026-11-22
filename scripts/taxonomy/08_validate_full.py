#!/usr/bin/env python3
"""08 — Valida la capa completa (400) contra las fuentes de evidencia.

  * 400 JSON, todos parseables y con los campos del contrato.
  * `ani`/`closest_ref` deben existir como par en `full/ani_hits_full.tsv` y la
    referencia debe estar en `full/panel_full.tsv`.
  * `gtdb_lineage` debe coincidir literalmente con GTDB r232 o ser "unknown"
    (y "unknown" solo se acepta si el accession no está en GTDB).
  * `ncbi_taxid` debe existir en `data/taxdump/nodes.dmp`.
  * `species` debe estar respaldada por NCBI o GTDB.
  * `summary_full.json` debe cuadrar con los 400 JSON.
  * **Muestra**: 20 genomas al azar (semilla fija) re-verificados uno a uno.
  * **Separación de capas**: los 10 JSON del piloto no deben haber cambiado de
    pipeline (`1.0` vs `2.0-full`).

Salida: data/taxonomy/full/validation_full.json (exit 1 si hay errores)
"""

from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"
REQUIRED = ["genome", "closest_ref", "ani", "gtdb_lineage", "species",
            "confidence", "evidence"]
CONF = {"high", "medium", "low", "unknown"}
STATUS = {"assigned", "provisional", "unknown"}
SAMPLE_N = 20


def main() -> int:
    errors: list[str] = []
    warnings: list[str] = []
    per: list[dict] = []

    panel = set()
    with open(FULL / "panel_full.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panel.add(r["genome"])

    hits: dict[tuple[str, str], float] = {}
    with open(FULL / "ani_hits_full.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hits[(r["genome"], r["ref"])] = float(r["ani"])

    gtdb = L.load_gtdb_taxonomy()
    meta = {}
    with open(FULL / "genome_meta.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["genome"]] = r

    docs: dict[str, dict] = {}
    for genome in meta:
        path = FULL / f"{genome}.json"
        rec = {"genome": genome, "errors": [], "warnings": []}
        if not path.exists():
            errors.append(f"{genome}: JSON ausente")
            rec["errors"].append("JSON ausente")
            per.append(rec)
            continue
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            errors.append(f"{genome}: JSON invalido ({exc})")
            rec["errors"].append("JSON invalido")
            per.append(rec)
            continue
        docs[genome] = doc

        for f in REQUIRED:
            if f not in doc:
                rec["errors"].append(f"falta {f}")
        if doc.get("genome") != genome:
            rec["errors"].append("campo genome no coincide con el fichero")
        if doc.get("layer") != "full":
            rec["errors"].append("layer != 'full'")

        ani = doc.get("ani")
        ref = doc.get("closest_ref")
        if ani is not None:
            if not isinstance(ani, (int, float)) or not (70.0 <= float(ani) <= 100.0):
                rec["errors"].append(f"ani fuera de rango: {ani}")
            if not ref:
                rec["errors"].append("ani numerico sin closest_ref")
            elif ref not in panel:
                rec["errors"].append(f"closest_ref {ref} fuera del panel")
            else:
                h = hits.get((genome, ref))
                if h is None:
                    rec["errors"].append("par ANI ausente en ani_hits_full.tsv")
                elif abs(h - float(ani)) > 0.011:
                    rec["errors"].append(f"ani {ani} != hits {h}")
        elif ref:
            rec["warnings"].append("ani=null con closest_ref no nulo")

        gl = doc.get("gtdb_lineage")
        acc = L.pilot_accession(genome)
        expected = gtdb.get("RS_" + acc) or gtdb.get("GB_" + acc)
        if gl == "unknown":
            if expected:
                rec["errors"].append("gtdb_lineage=unknown pero el accession SI esta en GTDB")
        elif not isinstance(gl, str) or not gl.startswith(("d__Bacteria", "d__Archaea")):
            rec["errors"].append("gtdb_lineage no empieza por d__")
        elif expected is None:
            rec["errors"].append("gtdb_lineage presente pero accession ausente de GTDB")
        elif expected != gl:
            rec["errors"].append("gtdb_lineage no coincide con GTDB r232")

        if doc.get("confidence") not in CONF:
            rec["errors"].append(f"confidence invalida: {doc.get('confidence')}")
        if doc.get("status") not in STATUS:
            rec["errors"].append(f"status invalido: {doc.get('status')}")

        species = doc.get("species")
        if species != "unknown":
            if not isinstance(species, str) or len(species.split()) < 2:
                rec["errors"].append(f"species no binomial: {species}")
            else:
                backed = species in (doc.get("species_ncbi") or "", doc.get("species_gtdb") or "")
                ranks = doc.get("evidence", {}).get("ncbi", {}).get("ranks", {})
                if not backed and ranks.get("species") != species:
                    rec["errors"].append(f"species '{species}' sin respaldo NCBI/GTDB")

        if doc.get("ncbi_taxid") is not None and not isinstance(doc["ncbi_taxid"], int):
            rec["errors"].append("ncbi_taxid no entero")

        for e in rec["errors"]:
            errors.append(f"{genome}: {e}")
        for w in rec["warnings"]:
            warnings.append(f"{genome}: {w}")
        per.append(rec)

    # taxids reales
    wanted = {d["ncbi_taxid"] for d in docs.values() if isinstance(d.get("ncbi_taxid"), int)}
    existing: set[int] = set()
    with open(L.TAXDUMP / "nodes.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.split("\t|\t")
            if len(p) >= 3:
                existing.add(int(p[0]))
    for g, d in docs.items():
        t = d.get("ncbi_taxid")
        if isinstance(t, int) and t not in existing:
            errors.append(f"{g}: taxid {t} no existe en nodes.dmp")

    # summary
    summary = None
    sp = FULL / "summary_full.json"
    if not sp.exists():
        errors.append("summary_full.json ausente")
    else:
        summary = json.loads(sp.read_text(encoding="utf-8"))
        if summary.get("n_genomes") != len(docs):
            errors.append(f"summary n_genomes={summary.get('n_genomes')} != {len(docs)}")
        listed = {g["genome"] for g in summary.get("genomes", [])}
        if listed != set(docs):
            errors.append("summary.genomes != JSON en disco")
        n_ani = sum(1 for d in docs.values() if d.get("ani") is not None)
        if summary.get("n_ani_computed") != n_ani:
            errors.append(f"summary n_ani_computed={summary.get('n_ani_computed')} != {n_ani}")
        for g in summary.get("genomes", []):
            d = docs.get(g["genome"], {})
            if d and (g.get("ani") != d.get("ani") or g.get("species") != d.get("species")):
                errors.append(f"summary difiere para {g['genome']}")

    # separación de capas: piloto intacto
    pilot_files = sorted(L.OUT.glob("GCF_*.json"))
    pilot_docs = [json.loads(p.read_text(encoding="utf-8")) for p in pilot_files]
    pilot_versions = {d.get("taxonomy_version", {}).get("pipeline") for d in pilot_docs}
    if pilot_versions != {"1.0"}:
        errors.append(f"los JSON del piloto cambiaron de version: {pilot_versions}")
    if len(pilot_files) != 10:
        errors.append(f"se esperaban 10 JSON de piloto, hay {len(pilot_files)}")

    # muestra re-verificada
    random.seed(20261120)
    sample = sorted(random.sample(sorted(docs), min(SAMPLE_N, len(docs))))
    sample_checks = []
    for g in sample:
        d = docs[g]
        ref = d.get("closest_ref")
        row = meta[g]
        chk = {
            "genome": g,
            "ani": d.get("ani"),
            "closest_ref": ref,
            "ref_in_panel": ref in panel if ref else None,
            "ref_hit_matches": (abs(hits.get((g, ref), -1) - float(d["ani"])) <= 0.011)
                               if ref and d.get("ani") is not None else None,
            "gtdb_literal_ok": (d.get("gtdb_lineage") == "unknown"
                                and (gtdb.get("RS_" + L.pilot_accession(g))
                                     or gtdb.get("GB_" + L.pilot_accession(g))) is None)
                               or (d.get("gtdb_lineage") == (gtdb.get("RS_" + L.pilot_accession(g))
                                                             or gtdb.get("GB_" + L.pilot_accession(g)))),
            "species_backed": d.get("species") in (d.get("species_ncbi"),
                                                   d.get("species_gtdb")),
            "ncbi_organism": row.get("organism_name"),
            "scope_rank": d.get("evidence", {}).get("ani", {}).get("scope_rank"),
            "n_hits": d.get("evidence", {}).get("ani", {}).get("n_hits"),
        }
        sample_checks.append(chk)
    bad_sample = [c["genome"] for c in sample_checks
                  if c["ref_in_panel"] is False or c["ref_hit_matches"] is False
                  or c["gtdb_literal_ok"] is False]
    if bad_sample:
        errors.append(f"muestra con fallos: {bad_sample}")

    result = {
        "ok": not errors,
        "layer": "full",
        "n_genomes": len(docs),
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "errors": errors,
        "warnings": warnings,
        "sample_size": len(sample_checks),
        "sample_checks": sample_checks,
        "coverage": {
            "ani": sum(1 for d in docs.values() if d.get("ani") is not None),
            "ani_null": sum(1 for d in docs.values() if d.get("ani") is None),
            "gtdb": sum(1 for d in docs.values() if d.get("gtdb_lineage") != "unknown"),
            "gtdb_unknown": sum(1 for d in docs.values() if d.get("gtdb_lineage") == "unknown"),
            "species_named": sum(1 for d in docs.values() if d.get("species") != "unknown"),
        },
        "pilot_layer": {"n_files": len(pilot_files),
                        "pipelines": sorted(str(v) for v in pilot_versions)},
        "per_genome_errors": [p for p in per if p["errors"]],
    }
    L.write_json(FULL / "validation_full.json", result)
    print(json.dumps({k: result[k] for k in
                      ("ok", "n_genomes", "n_errors", "n_warnings", "coverage",
                       "sample_size", "errors", "warnings", "pilot_layer")},
                     indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
