#!/usr/bin/env python3
"""04 — Valida los JSON de `data/taxonomy/`.

No valida "que suene bien": valida que cada valor tenga respaldo en un fichero
de evidencia y que nada se haya inventado.

  * `ani` y `closest_ref` deben existir como fila en `ani_hits.tsv`.
  * `closest_ref` debe existir en `refpanel.tsv`.
  * `gtdb_lineage` debe coincidir literalmente con la taxonomia GTDB cacheada.
  * `ncbi_taxid` debe existir en `data/taxdump/nodes.dmp`.
  * `species` debe aparecer en el linaje NCBI o GTDB del propio genoma.
  * `summary.json` debe cuadrar con los 10 JSON.

Salida: data/taxonomy/validation.json  (exit 1 si hay errores)
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

REQUIRED = ["genome", "closest_ref", "ani", "gtdb_lineage", "species",
            "confidence", "evidence"]
CONFIDENCE = {"high", "medium", "low", "unknown"}
STATUS = {"assigned", "provisional", "unknown"}


def main() -> int:
    L.ensure_dirs()
    pilots = L.read_pilot_tsv()
    pilot_keys = [p["genome"] for p in pilots]
    errors: list[str] = []
    warnings: list[str] = []
    checks: list[dict] = []

    panel = set()
    with open(L.OUT / "refpanel.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            panel.add(r["genome"])

    hits: dict[tuple[str, str], float] = {}
    with open(L.OUT / "ani_hits.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hits[(r["genome"], r["ref"])] = float(r["ani"])

    gtdb = L.load_gtdb_taxonomy()
    taxids_in_taxdump: set[int] = set()
    wanted: set[int] = set()
    docs: dict[str, dict] = {}

    for key in pilot_keys:
        path = L.OUT / f"{key}.json"
        rec = {"genome": key, "file": path.name, "errors": [], "warnings": []}
        if not path.exists():
            rec["errors"].append("JSON ausente")
            errors.append(f"{key}: JSON ausente")
            checks.append(rec)
            continue
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            rec["errors"].append(f"JSON invalido: {exc}")
            errors.append(f"{key}: JSON invalido ({exc})")
            checks.append(rec)
            continue
        docs[key] = doc

        for field in REQUIRED:
            if field not in doc:
                rec["errors"].append(f"falta campo {field}")
        if doc.get("genome") != key:
            rec["errors"].append(f"campo genome='{doc.get('genome')}' != nombre de archivo")

        ani = doc.get("ani")
        if ani is not None:
            if not isinstance(ani, (int, float)):
                rec["errors"].append("ani no es numerico")
            elif not (70.0 <= float(ani) <= 100.0):
                rec["errors"].append(f"ani fuera de rango razonable: {ani}")
            else:
                ref = doc.get("closest_ref")
                if not ref:
                    rec["errors"].append("ani numerico sin closest_ref")
                else:
                    if ref not in panel:
                        rec["errors"].append(f"closest_ref {ref} no esta en refpanel.tsv")
                    hit = hits.get((key, ref))
                    if hit is None:
                        rec["errors"].append(f"par ({key},{ref}) no esta en ani_hits.tsv")
                    elif abs(hit - float(ani)) > 0.011:
                        rec["errors"].append(
                            f"ani {ani} != ani_hits.tsv {hit} para {ref}")
        else:
            if doc.get("closest_ref"):
                rec["warnings"].append("ani=null pero closest_ref no es null")
                if doc["closest_ref"] not in panel:
                    rec["errors"].append("closest_ref no esta en refpanel.tsv")

        gl = doc.get("gtdb_lineage")
        if gl == "unknown":
            rec["warnings"].append("gtdb_lineage=unknown")
        elif not isinstance(gl, str):
            rec["errors"].append("gtdb_lineage no es string")
        else:
            if not gl.startswith(("d__Bacteria", "d__Archaea")):
                rec["errors"].append(f"gtdb_lineage no empieza por d__: {gl[:30]}")
            acc = L.pilot_accession(key)
            expected = gtdb.get("RS_" + acc) or gtdb.get("GB_" + acc)
            if expected is None:
                rec["errors"].append(f"accesion {acc} ausente de la taxonomia GTDB")
            elif expected != gl:
                rec["errors"].append("gtdb_lineage no coincide con la taxonomia GTDB")

        conf = doc.get("confidence")
        if conf not in CONFIDENCE:
            rec["errors"].append(f"confidence invalida: {conf}")
        if doc.get("status") not in STATUS:
            rec["errors"].append(f"status invalido: {doc.get('status')}")

        species = doc.get("species")
        if species != "unknown":
            if not isinstance(species, str) or len(species.split()) < 2:
                rec["errors"].append(f"species no es binomial: {species}")
            else:
                in_ncbi = species == (doc.get("species_ncbi") or "")
                in_gtdb = species == (doc.get("species_gtdb") or "")
                lineage_species = [
                    r for r in doc.get("evidence", {}).get("ncbi", {}).get("lineage_ranks", [])
                    if r.get("rank") == "species"
                ]
                in_lineage = any(r.get("name") == species for r in lineage_species)
                if not (in_ncbi or in_gtdb or in_lineage):
                    rec["errors"].append(
                        f"species '{species}' no respaldada por NCBI ni GTDB")
            # chequeo de phylum declarado vs linaje NCBI
            phylum_declared = doc.get("evidence", {}).get("checks", {}).get("pilot_phylum_declared")
            ranks = {r["rank"]: r["name"] for r in
                     doc.get("evidence", {}).get("ncbi", {}).get("lineage_ranks", [])}
            phylum_ncbi = ranks.get("phylum")
            if phylum_declared and phylum_ncbi and phylum_declared != phylum_ncbi:
                rec["warnings"].append(
                    f"phylum declarado '{phylum_declared}' != NCBI '{phylum_ncbi}'")

        tid = doc.get("ncbi_taxid")
        if tid is not None:
            if not isinstance(tid, int):
                rec["errors"].append("ncbi_taxid no es entero")
            else:
                wanted.add(tid)

        for msg in rec["errors"]:
            errors.append(f"{key}: {msg}")
        checks.append(rec)

    # taxids reales en el taxdump
    if wanted:
        with open(L.TAXDUMP / "nodes.dmp", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                p = line.split("\t|\t")
                if len(p) >= 3:
                    taxids_in_taxdump.add(int(p[0]))
        for key, doc in docs.items():
            tid = doc.get("ncbi_taxid")
            if isinstance(tid, int) and tid not in taxids_in_taxdump:
                errors.append(f"{key}: taxid {tid} no existe en nodes.dmp")

    # summary
    summary = None
    sp = L.OUT / "summary.json"
    if not sp.exists():
        errors.append("summary.json ausente")
    else:
        summary = json.loads(sp.read_text(encoding="utf-8"))
        if summary.get("n_genomes") != len(docs):
            errors.append(
                f"summary.n_genomes={summary.get('n_genomes')} != {len(docs)} JSON")
        listed = {g["genome"] for g in summary.get("genomes", [])}
        if listed != set(docs):
            errors.append("summary.genomes no coincide con los JSON en disco")
        n_ani = sum(1 for d in docs.values() if d.get("ani") is not None)
        if summary.get("n_ani_computed") != n_ani:
            errors.append(
                f"summary.n_ani_computed={summary.get('n_ani_computed')} != {n_ani}")
        for g in summary.get("genomes", []):
            d = docs.get(g["genome"], {})
            if d and g.get("ani") != d.get("ani"):
                errors.append(f"summary: ani distinto para {g['genome']}")
            if d and g.get("species") != d.get("species"):
                errors.append(f"summary: species distinto para {g['genome']}")

    result = {
        "ok": not errors,
        "n_genomes": len(docs),
        "n_errors": len(errors),
        "n_warnings": sum(len(c["warnings"]) for c in checks),
        "errors": errors,
        "per_genome": checks,
    }
    L.write_json(L.OUT / "validation.json", result)

    print(json.dumps({"ok": result["ok"], "errors": errors,
                      "warnings": [w for c in checks for w in c["warnings"]]},
                     indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
