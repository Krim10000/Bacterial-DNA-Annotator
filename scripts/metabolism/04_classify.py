#!/usr/bin/env python3
"""04 - Clasifica cada gap del piloto en T1 / T1_review / T2 / T3.

Regla (ver README):
  T1_bacterial  hit a un acc Swiss-Prot **bacteriano** ligado por
                `rhea2uniprot_sprot.tsv` a una reacción Rhea de transporte
                exacta del metabolito (ChEBI a ambos lados), con
                pident >= 35,0, qcov >= 0,70 y evalue <= 1e-5.
  T1_review     hay hit a un candidato, pero es eucariota/arquea o no alcanza
                los umbrales. **No es T1 y no habilita importación.**
  T2            `transport_reaction_only`: existe reacción Rhea de transporte
                exacta para el metabolito, pero ningún hit T1.
  T3            `no_transport`: no existe reacción de transporte exacta (el
                metabolito no aparece, o aparece solo en un lado de la reacción).

Entradas: ref/ + pilot_genomes.tsv + diamond/<genome>.hits.tsv +
          data/metabolism/transport_audit.tsv (solo lectura)
Salidas:  data/metabolism/transport_t1/<genome>.tsv (10)
          data/metabolism/transport_t1/t1_hits.tsv
          data/metabolism/transport_t1/summary.json
"""
from __future__ import annotations

import importlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
common = importlib.import_module("00_common")


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def load_hits(path: Path) -> dict[str, list[tuple[float, float, float, float, str]]]:
    """acc -> [(pident, qcov_fraccion, evalue, bitscore, query), ...]."""
    hits: dict[str, list] = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 10:
                continue
            try:
                acc = f[1].split("|")[1]
            except IndexError:
                continue
            # DIAMOND qcovhsp viene en porcentaje (0-100); se normaliza a fracción.
            hits[acc].append(
                (float(f[2]), float(f[8]) / 100.0, float(f[6]), float(f[7]), f[0])
            )
    return dict(hits)


def rank_key(hit: tuple[float, float, float, float, str]) -> tuple:
    """Mejor hit: bitscore alto, evalue bajo, pident alto."""
    pident, _qcov, evalue, bitscore, _query = hit
    return (-bitscore, evalue, -pident)


def first_failure(meta: dict, hit: tuple) -> str:
    pident, qcov, evalue, _bits, _query = hit
    reasons = []
    if meta["superkingdom"] != "Bacteria":
        reasons.append(f"taxid_{meta['superkingdom'].lower()}")
    if pident < common.MIN_PIDENT:
        reasons.append(f"pident<{common.MIN_PIDENT:g}")
    if qcov < common.MIN_QCOV:
        reasons.append(f"qcov<{common.MIN_QCOV:g}")
    if evalue > common.MAX_EVALUE:
        reasons.append(f"evalue>{common.MAX_EVALUE:g}")
    return "|".join(reasons) if reasons else ""


def main() -> int:
    db, chebi2rhea, rhea2acc, acc2rhea = common.load_ref()
    _, _, one_side = common.load_rhea()

    pilot = list(common.read_tsv(common.OUT / "pilot_genomes.tsv"))
    pilot_genomes = [row["genome"] for row in pilot]
    if len(pilot_genomes) != 10:
        raise SystemExit(f"ERROR: se esperaban 10 pilotos, hay {len(pilot_genomes)}")

    gaps_by_genome: dict[str, list[dict]] = defaultdict(list)
    for row in common.read_tsv(common.TRANSPORT_AUDIT):
        if row["genome"] in pilot_genomes:
            gaps_by_genome[row["genome"]].append(row)
    missing = [g for g in pilot_genomes if not gaps_by_genome.get(g)]
    if missing:
        raise SystemExit(f"ERROR: pilotos sin gaps en transport_audit: {missing}")

    status_counts: dict[str, Counter] = {}
    t1_rows = []
    all_rows: dict[str, list] = {}

    for prow in pilot:
        genome = prow["genome"]
        hits_path = common.DIAMOND_DIR / f"{genome}.hits.tsv"
        if not hits_path.exists():
            raise SystemExit(f"ERROR: falta {hits_path}; ejecutar 03_run_diamond.py")
        hits = load_hits(hits_path)
        counts = Counter()
        rows = []

        for gap in gaps_by_genome[genome]:
            name = common.clean(gap["gap"])
            chebi = common.clean(gap["chebi"])
            rheas = chebi2rhea.get(chebi, []) if chebi else []
            row = {
                "genome": genome,
                "gap": name,
                "chebi": chebi,
                "rhea": "",
                "acc": "",
                "pident": "",
                "qcov": "",
                "evalue": "",
                "status": "",
                "evidence": "",
            }

            if not rheas:
                row["status"] = "no_transport"
                if not chebi:
                    row["evidence"] = "missing_chebi"
                elif chebi in one_side:
                    row["evidence"] = "no_exact_transport_rhea;one_side_only_reaction"
                else:
                    row["evidence"] = "no_exact_transport_rhea"
                counts["no_transport"] += 1
                rows.append(row)
                continue

            candidates = set()
            for rid in rheas:
                candidates |= rhea2acc.get(rid, set())
            row["rhea"] = ";".join(f"RHEA:{r}" for r in rheas[:5])

            candidate_hits = []  # (acc, hit)
            for acc in candidates:
                for hit in hits.get(acc, ()):
                    candidate_hits.append((acc, hit))

            passing = []
            for acc, hit in candidate_hits:
                meta = db.get(acc)
                if meta is None:
                    continue
                pident, qcov, evalue = hit[0], hit[1], hit[2]
                if (
                    meta["superkingdom"] == "Bacteria"
                    and pident >= common.MIN_PIDENT
                    and qcov >= common.MIN_QCOV
                    and evalue <= common.MAX_EVALUE
                ):
                    passing.append((acc, hit))
            if passing:
                passing.sort(key=lambda x: rank_key(x[1]))
                acc, hit = passing[0]
                confidence = (
                    "high"
                    if hit[0] >= 40.0 and hit[1] >= 0.80 and hit[2] <= 1e-30
                    else "medium"
                )
                linked = sorted(r for r in rheas if acc in rhea2acc.get(r, set()))
                row.update(
                    {
                        "rhea": ";".join(f"RHEA:{r}" for r in (linked or rheas)[:5]),
                        "acc": acc,
                        "pident": common.fmt_pident(hit[0]),
                        "qcov": common.fmt_qcov(hit[1]),
                        "evalue": common.fmt_evalue(hit[2]),
                        "status": "t1_bacterial",
                        "evidence": (
                            "diamond_allhits;rhea2uniprot_exact;chebi_both_sides;"
                            f"taxid_bacteria;thresholds_pass;confidence={confidence};"
                            f"query={hit[4]}"
                        ),
                    }
                )
                counts["t1_bacterial"] += 1
                t1_rows.append(
                    [
                        genome,
                        name,
                        chebi,
                        row["rhea"],
                        acc,
                        row["pident"],
                        row["qcov"],
                        row["evalue"],
                        db[acc]["organism"],
                        prow["phylum"],
                        confidence,
                    ]
                )
            elif candidate_hits:
                candidate_hits.sort(key=lambda x: rank_key(x[1]))
                acc, hit = candidate_hits[0]
                meta = db.get(acc, {"superkingdom": "unknown"})
                fail = first_failure(meta, hit) or "below_threshold"
                linked = sorted(r for r in rheas if acc in rhea2acc.get(r, set()))
                row.update(
                    {
                        "rhea": ";".join(f"RHEA:{r}" for r in (linked or rheas)[:5]),
                        "acc": acc,
                        "pident": common.fmt_pident(hit[0]),
                        "qcov": common.fmt_qcov(hit[1]),
                        "evalue": common.fmt_evalue(hit[2]),
                        "status": "t1_review",
                        "evidence": (
                            "diamond_allhits;rhea2uniprot_exact;chebi_both_sides;"
                            f"taxid={meta['superkingdom'].lower()};fail={fail};"
                            f"confidence=low;query={hit[4]}"
                        ),
                    }
                )
                counts["t1_review"] += 1
            else:
                row["status"] = "transport_reaction_only"
                if candidates:
                    row["evidence"] = (
                        "rhea_transport_exact;chebi_both_sides;"
                        f"n_rhea={len(rheas)};n_candidates={len(candidates)};"
                        "no_genome_hit"
                    )
                else:
                    row["evidence"] = (
                        "rhea_transport_exact;chebi_both_sides;"
                        f"n_rhea={len(rheas)};no_curated_swissprot"
                    )
                counts["transport_reaction_only"] += 1
            rows.append(row)

        rows.sort(key=lambda r: (common.STATUS_RANK[r["status"]], r["gap"]))
        all_rows[genome] = rows
        status_counts[genome] = counts
        common.write_tsv(
            common.OUT / f"{genome}.tsv",
            common.TSV_COLUMNS,
            ([r[c] for c in common.TSV_COLUMNS] for r in rows),
        )
        print(
            f"[04] {genome}: {len(rows)} gaps -> T1={counts['t1_bacterial']} "
            f"review={counts['t1_review']} T2={counts['transport_reaction_only']} "
            f"T3={counts['no_transport']}",
            flush=True,
        )

    common.write_tsv(
        common.OUT / "t1_hits.tsv",
        [
            "genome",
            "gap",
            "chebi",
            "rhea",
            "acc",
            "pident",
            "qcov",
            "evalue",
            "organism",
            "phylum",
            "confidence",
        ],
        t1_rows,
    )

    # --- summary ---------------------------------------------------------------
    db_stats = {}
    ref_stats_path = common.REF / "ref_stats.json"
    if ref_stats_path.exists():
        db_stats = json.load(open(ref_stats_path, encoding="utf-8"))
    run_stats_path = common.DIAMOND_DIR / "run_stats.json"
    run_stats = json.load(open(run_stats_path, encoding="utf-8")) if run_stats_path.exists() else {}

    totals = Counter()
    for counts in status_counts.values():
        totals.update(counts)
    total_gaps = sum(totals.values())
    if total_gaps != totals["t1_bacterial"] + totals["t1_review"] + totals["transport_reaction_only"] + totals["no_transport"]:
        raise SystemExit("ERROR: los conteos de estado no suman el total de gaps")

    pilot_out = []
    for prow in pilot:
        g = prow["genome"]
        c = status_counts[g]
        pilot_out.append(
            {
                "genome": g,
                "phylum": prow["phylum"],
                "organism": prow["organism"],
                "proteins": int(prow["proteins"]),
                "gaps_total": len(all_rows[g]),
                "t1_bacterial": c["t1_bacterial"],
                "t1_review": c["t1_review"],
                "transport_reaction_only": c["transport_reaction_only"],
                "no_transport": c["no_transport"],
            }
        )

    t1_by_phylum = Counter(r[9] for r in t1_rows)
    inputs = [
        common.TRANSPORT_AUDIT,
        common.GAP_TABLE,
        common.RHEA2ACC,
        common.RHEA_PKL,
        common.TRANSPORT_FASTA,
        common.REF / "transport_db.tsv",
        common.REF / "chebi_transport_rhea.tsv",
        common.REF / "rhea2acc_transport.tsv",
    ]
    summary = {
        "schema_version": 1,
        "module": "transport_t1",
        "agent": "transport-builder",
        "generated_at": now(),
        "scope": (
            "10 genomas piloto; capa T1 bacteriana de transporte con sustrato, "
            "separada de la capa T2"
        ),
        "definitions": {
            "t1_bacterial": (
                "transportador Swiss-Prot bacteriano ligado a la reacción Rhea de "
                "transporte exacta del metabolito; pident>=35, qcov>=0.70, evalue<=1e-5"
            ),
            "t1_review": (
                "hay hit a un candidato, pero es no bacteriano o no alcanza los "
                "umbrales; NO es T1 y no habilita importación"
            ),
            "transport_reaction_only": (
                "T2: existe reacción Rhea de transporte exacta del metabolito, "
                "pero ningún hit T1"
            ),
            "no_transport": (
                "T3: no existe reacción Rhea de transporte exacta (el metabolito "
                "no aparece, o aparece solo en un lado de la reacción)"
            ),
        },
        "transport_definition": {
            "source": "data/metabolism/rhea_reactions.pkl + data/metabolism/rhea2uniprot_sprot.tsv",
            "exact_substrate_rule": "chebi presente en substrates y products de la misma reacción Rhea",
            "why": (
                "el flag transport de Rhea también marca reacciones que comparten "
                "un ion (H+/Na+); exigir el metabolito a ambos lados evita T2 falsos "
                "como acil-CoA sintetasas y decarboxilasas"
            ),
        },
        "inputs": {
            str(p.relative_to(common.ROOT)): {
                "sha256": common.sha256_file(p),
                "rows": common.count_rows(p) if p.suffix == ".tsv" else None,
            }
            for p in inputs
        },
        "database": {
            "fasta": "data/db/transport.fasta",
            "dmnd": "data/db/transport.dmnd",
            "entries": db_stats.get("transport_db_entries"),
            "superkingdom": db_stats.get("superkingdom"),
            "rhea_transport_flag": db_stats.get("rhea_transport_flag"),
            "rhea_transport_exact_substrate": db_stats.get("rhea_transport_exact_substrate"),
            "metabolites_with_exact_transport": db_stats.get("metabolites_with_exact_transport"),
            "rhea_transport_with_swissprot": db_stats.get("rhea_transport_with_swissprot"),
            "rhea_transport_without_swissprot": db_stats.get("rhea_transport_without_swissprot"),
            "parse_check": db_stats.get("parse_check"),
        },
        "method": {
            "diamond": run_stats.get("parameters", {}),
            "candidate_filter": (
                "acc ligado por rhea2uniprot_sprot a la misma Rhea de transporte "
                "exacta del ChEBI del gap"
            ),
            "bacterial_filter": "superkingdom Bacteria resuelto con data/taxdump/nodes.dmp",
            "thresholds": {
                "pident_min": common.MIN_PIDENT,
                "qcov_min": common.MIN_QCOV,
                "evalue_max": common.MAX_EVALUE,
                "search_evalue_max": common.SEARCH_EVALUE,
            },
            "hit_ranking": "bitscore desc, evalue asc, pident desc (HSP de mayor bitscore)",
            "confidence": (
                "high: pident>=40, qcov>=0.80 y evalue<=1e-30; medium: cumple "
                "umbrales T1 pero no 'high'; low: t1_review"
            ),
            "row_order": "T1, T1_review, T2, T3; luego gap alfabético",
            "one_row_per_gap": True,
        },
        "pilot": pilot_out,
        "totals": {
            "genomes": len(pilot_out),
            "gaps": total_gaps,
            "t1_bacterial": totals["t1_bacterial"],
            "t1_review": totals["t1_review"],
            "transport_reaction_only": totals["transport_reaction_only"],
            "no_transport": totals["no_transport"],
            "t1_by_phylum": dict(sorted(t1_by_phylum.items())),
        },
        "t1": [
            {
                "genome": r[0],
                "gap": r[1],
                "chebi": r[2],
                "rhea": r[3],
                "acc": r[4],
                "pident": float(r[5]),
                "qcov": float(r[6]),
                "evalue": float(r[7]),
                "organism": r[8],
                "phylum": r[9],
                "confidence": r[10],
            }
            for r in t1_rows
        ],
        "tools": {
            "diamond": run_stats.get("diamond", "desconocida"),
            "python": sys.version.split()[0],
        },
        "limits": [
            "solo 306 de los 400 genomas base tienen gaps en transport_audit.tsv",
            "solo las Rhea de transporte con proteína Swiss-Prot curada pueden dar T1; las demás quedan T2 no buscables",
            "el transporte por grupo (PTS) y las translocaciones con modificación química no se detectan porque el ChEBI cambia entre lados",
            "no se usan dominios Pfam/TM ni contexto de operón; un hit es evidencia de homología, no de función medida",
            "qcov es la cobertura del HSP de mayor bitscore, no la cobertura combinada de varios HSPs",
            "las Rhea de transporte son de dirección no especificada (UN): T1 prueba transporte del sustrato, no importación neta",
            "en 2 de los 4 T1 (D-metionina) el mejor hit Swiss-Prot de la consulta es una ATPasa ABC no específica (YknY/MceG), no MetN; el T1 se apoya en un hit secundario a MetN y requiere revisión antes de marcar importación",
            "T1_review no es T1: no debe usarse para marcar importación",
        ],
        "validation": "ver validation.json (05_validate.py)",
    }
    common.write_json(common.OUT / "summary.json", summary)

    print(
        f"[04] total: T1={totals['t1_bacterial']} review={totals['t1_review']} "
        f"T2={totals['transport_reaction_only']} T3={totals['no_transport']} "
        f"de {total_gaps} gaps"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
