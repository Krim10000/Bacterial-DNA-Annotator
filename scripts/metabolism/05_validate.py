#!/usr/bin/env python3
"""05 - Validación independiente de la capa T1 de transporte.

Comprueba estructura, semántica y conteos:
  - 10 TSV con las 10 columnas exactas y una fila por gap del piloto;
  - T1: acc bacteriano, métricas sobre umbral, Rhea de transporte exacta del
    ChEBI y acc ligado a esa Rhea;
  - T1_review: tiene hit y **no** cumple el filtro T1 (no se cuenta como T1);
  - T2: Rhea de transporte exacta y sin acc/métricas;
  - T3: sin Rhea de transporte exacta y sin acc/métricas;
  - conteos por genoma y globales iguales a summary.json;
  - t1_hits.tsv coincide con las filas T1;
  - hashes de las entradas iguales a los registrados en summary.json
    (evidencia de que las tablas base no se modificaron).

Salida: data/metabolism/transport_t1/validation.json
"""
from __future__ import annotations

import importlib
import json
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
common = importlib.import_module("00_common")


def main() -> int:
    errors: list[str] = []
    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if not ok:
            errors.append(f"{name}: {detail}")

    summary_path = common.OUT / "summary.json"
    if not summary_path.exists():
        raise SystemExit("ERROR: falta summary.json; ejecutar 04_classify.py")
    summary = json.load(open(summary_path, encoding="utf-8"))

    db, chebi2rhea, rhea2acc, _acc2rhea = common.load_ref()
    _, _, one_side = common.load_rhea()

    # --- 1. integridad de la referencia ---------------------------------------
    ref_stats_path = common.REF / "ref_stats.json"
    check("ref_stats_existe", ref_stats_path.exists())
    ref_stats = json.load(open(ref_stats_path, encoding="utf-8")) if ref_stats_path.exists() else {}
    check(
        "db_entries_coinciden",
        ref_stats.get("transport_db_entries") == len(db),
        f"ref={ref_stats.get('transport_db_entries')} db={len(db)}",
    )
    check(
        "taxonomia_resuelta",
        all(m["superkingdom"] in {"Bacteria", "Eukaryota", "Archaea", "other"} for m in db.values()),
        "superkingdom resuelto para todas las entradas de la base",
    )
    check(
        "parse_check_ok",
        ref_stats.get("parse_check", {}).get("db_accs_equal_rhea2uniprot_accs") is True,
    )

    # --- 2. pilotos y archivos -------------------------------------------------
    pilot = list(common.read_tsv(common.OUT / "pilot_genomes.tsv"))
    pilot_genomes = [p["genome"] for p in pilot]
    check("pilot_10", len(pilot) == 10, f"n={len(pilot)}")
    check("pilot_genomas_unicos", len(set(pilot_genomes)) == len(pilot_genomes))

    gaps_by_genome: dict[str, list[dict]] = {}
    for row in common.read_tsv(common.TRANSPORT_AUDIT):
        if row["genome"] in pilot_genomes:
            gaps_by_genome.setdefault(row["genome"], []).append(row)

    expected_header = common.TSV_COLUMNS
    recomputed_totals = Counter()
    genome_results = {}
    t1_from_files = []

    for genome in pilot_genomes:
        path = common.OUT / f"{genome}.tsv"
        if not path.exists():
            check(f"tsv_existe:{genome}", False, str(path))
            continue
        with open(path, encoding="utf-8") as fh:
            lines = [line.rstrip("\n") for line in fh if line.strip()]
        header = lines[0].split("\t") if lines else []
        check(f"columnas:{genome}", header == expected_header, f"header={header}")
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]

        audit = gaps_by_genome.get(genome, [])
        check(
            f"filas_vs_audit:{genome}",
            len(rows) == len(audit),
            f"tsv={len(rows)} audit={len(audit)}",
        )
        expected_keys = {(r["genome"], r["gap"], r["chebi"]) for r in audit}
        got_keys = {(r["genome"], r["gap"], r["chebi"]) for r in rows}
        check(f"claves:{genome}", expected_keys == got_keys, "claves (genome,gap,chebi) iguales")
        check(f"sin_duplicados:{genome}", len(got_keys) == len(rows))

        counts = Counter(r["status"] for r in rows)
        check(
            f"estados_validos:{genome}",
            set(counts) <= set(common.STATUSES),
            f"estados={sorted(counts)}",
        )
        check(
            f"suma_estados:{genome}",
            sum(counts.values()) == len(rows),
            f"suma={sum(counts.values())} filas={len(rows)}",
        )

        for r in rows:
            st = r["status"]
            acc, chebi, rhea = r["acc"], r["chebi"], r["rhea"]
            rheas = [int(x.split("RHEA:")[-1]) for x in rhea.split(";") if x]
            has_metrics = bool(r["pident"] and r["qcov"] and r["evalue"])

            if st == "t1_bacterial":
                ok = True
                if not acc:
                    errors.append(f"{genome}/{r['gap']}: T1 sin acc")
                    ok = False
                else:
                    meta = db.get(acc)
                    if not meta:
                        errors.append(f"{genome}/{r['gap']}: acc {acc} fuera de la base")
                        ok = False
                    elif meta["superkingdom"] != "Bacteria":
                        errors.append(
                            f"{genome}/{r['gap']}: T1 con taxón {meta['superkingdom']} ({acc})"
                        )
                        ok = False
                if not has_metrics:
                    errors.append(f"{genome}/{r['gap']}: T1 sin métricas")
                    ok = False
                else:
                    try:
                        pident, qcov, evalue = float(r["pident"]), float(r["qcov"]), float(r["evalue"])
                    except ValueError:
                        errors.append(f"{genome}/{r['gap']}: métricas no numéricas")
                        ok = False
                    else:
                        if pident < common.MIN_PIDENT or qcov < common.MIN_QCOV or evalue > common.MAX_EVALUE:
                            errors.append(
                                f"{genome}/{r['gap']}: T1 bajo umbral pident={pident} qcov={qcov} e={evalue}"
                            )
                            ok = False
                if not rheas:
                    errors.append(f"{genome}/{r['gap']}: T1 sin Rhea")
                    ok = False
                else:
                    for rid in rheas:
                        if rid not in chebi2rhea.get(chebi, []):
                            errors.append(
                                f"{genome}/{r['gap']}: Rhea {rid} no es transporte exacto de {chebi}"
                            )
                            ok = False
                        if acc and acc not in rhea2acc.get(rid, set()):
                            errors.append(
                                f"{genome}/{r['gap']}: acc {acc} no está ligado a Rhea {rid}"
                            )
                            ok = False
                if "thresholds_pass" not in r["evidence"]:
                    errors.append(f"{genome}/{r['gap']}: evidencia T1 sin thresholds_pass")
                    ok = False
                conf = ""
                for token in r["evidence"].split(";"):
                    if token.startswith("confidence="):
                        conf = token.split("=", 1)[1]
                if conf not in {"high", "medium"}:
                    errors.append(f"{genome}/{r['gap']}: T1 con confidence inválida {conf!r}")
                    ok = False
                t1_from_files.append(
                    [
                        genome,
                        r["gap"],
                        chebi,
                        rhea,
                        acc,
                        r["pident"],
                        r["qcov"],
                        r["evalue"],
                        db.get(acc, {}).get("organism", ""),
                        next(p["phylum"] for p in pilot if p["genome"] == genome),
                        conf,
                    ]
                )
            elif st == "t1_review":
                if not acc or not has_metrics:
                    errors.append(f"{genome}/{r['gap']}: review sin acc/métricas")
                else:
                    meta = db.get(acc)
                    if meta is None:
                        errors.append(f"{genome}/{r['gap']}: review con acc {acc} fuera de la base")
                        meta = {"superkingdom": "unknown"}
                    try:
                        pident, qcov, evalue = float(r["pident"]), float(r["qcov"]), float(r["evalue"])
                    except ValueError:
                        errors.append(f"{genome}/{r['gap']}: review con métricas no numéricas")
                    else:
                        passes = (
                            meta["superkingdom"] == "Bacteria"
                            and pident >= common.MIN_PIDENT
                            and qcov >= common.MIN_QCOV
                            and evalue <= common.MAX_EVALUE
                        )
                        if passes:
                            errors.append(
                                f"{genome}/{r['gap']}: review que en realidad cumple T1"
                            )
                    if not rheas:
                        errors.append(f"{genome}/{r['gap']}: review sin Rhea")
            elif st == "transport_reaction_only":
                if acc or has_metrics:
                    errors.append(f"{genome}/{r['gap']}: T2 no debe tener acc/métricas")
                if not rheas:
                    errors.append(f"{genome}/{r['gap']}: T2 sin Rhea")
                for rid in rheas:
                    if rid not in chebi2rhea.get(chebi, []):
                        errors.append(
                            f"{genome}/{r['gap']}: T2 con Rhea {rid} no exacta para {chebi}"
                        )
            elif st == "no_transport":
                if acc or has_metrics or rheas:
                    errors.append(f"{genome}/{r['gap']}: T3 no debe tener acc/métricas/Rhea")
                if chebi and chebi in chebi2rhea:
                    errors.append(
                        f"{genome}/{r['gap']}: T3 pero {chebi} tiene transporte exacto"
                    )
            else:
                errors.append(f"{genome}/{r['gap']}: estado desconocido {st}")

        for status in common.STATUSES:
            recomputed_totals[status] += counts[status]
        genome_results[genome] = {
            "rows": len(rows),
            "t1_bacterial": counts["t1_bacterial"],
            "t1_review": counts["t1_review"],
            "transport_reaction_only": counts["transport_reaction_only"],
            "no_transport": counts["no_transport"],
        }

    # --- 3. t1_hits.tsv ---------------------------------------------------------
    t1_hits_path = common.OUT / "t1_hits.tsv"
    check("t1_hits_existe", t1_hits_path.exists())
    if t1_hits_path.exists():
        t1_file_rows = []
        for r in common.read_tsv(t1_hits_path):
            t1_file_rows.append(
                [
                    r["genome"],
                    r["gap"],
                    r["chebi"],
                    r["rhea"],
                    r["acc"],
                    r["pident"],
                    r["qcov"],
                    r["evalue"],
                    r["organism"],
                    r["phylum"],
                    r["confidence"],
                ]
            )
        check(
            "t1_hits_coincide",
            sorted(t1_file_rows) == sorted(t1_from_files),
            f"t1_hits={len(t1_file_rows)} filas_T1={len(t1_from_files)}",
        )
        check(
            "t1_hits_solo_piloto",
            all(r[0] in pilot_genomes for r in t1_file_rows),
            "todos los genomas de t1_hits están en el piloto",
        )

    # --- 4. summary -------------------------------------------------------------
    summary_totals = summary.get("totals", {})
    for status in common.STATUSES:
        check(
            f"summary_{status}",
            summary_totals.get(status) == recomputed_totals[status],
            f"summary={summary_totals.get(status)} recomputado={recomputed_totals[status]}",
        )
    check(
        "summary_gaps",
        summary_totals.get("gaps") == sum(recomputed_totals.values()),
        f"summary={summary_totals.get('gaps')} recomputado={sum(recomputed_totals.values())}",
    )
    check("summary_genomes", summary_totals.get("genomes") == 10)
    for p in summary.get("pilot", []):
        g = p["genome"]
        res = genome_results.get(g, {})
        check(
            f"summary_piloto:{g}",
            p.get("t1_bacterial") == res.get("t1_bacterial")
            and p.get("t1_review") == res.get("t1_review")
            and p.get("transport_reaction_only") == res.get("transport_reaction_only")
            and p.get("no_transport") == res.get("no_transport"),
            f"summary={p} tsv={res}",
        )

    # --- 4b. lista T1 del summary ----------------------------------------------
    summary_t1 = summary.get("t1", [])
    summary_t1_norm = [
        [
            r.get("genome", ""),
            r.get("gap", ""),
            r.get("chebi", ""),
            r.get("rhea", ""),
            r.get("acc", ""),
            f"{float(r.get('pident', 0)):.1f}",
            f"{float(r.get('qcov', 0)):.3f}",
            f"{float(r.get('evalue', 0)):.2e}",
            r.get("organism", ""),
            r.get("phylum", ""),
            r.get("confidence", ""),
        ]
        for r in summary_t1
    ]
    check(
        "summary_t1_coincide",
        sorted(summary_t1_norm) == sorted(t1_from_files),
        f"summary={len(summary_t1_norm)} filas_T1={len(t1_from_files)}",
    )

    # --- 5. entradas sin modificar ---------------------------------------------
    inputs_ok = True
    for rel, meta in summary.get("inputs", {}).items():
        path = common.ROOT / rel
        if not path.exists():
            inputs_ok = False
            errors.append(f"input ausente: {rel}")
            continue
        actual = common.sha256_file(path)
        if actual != meta.get("sha256"):
            inputs_ok = False
            errors.append(f"input modificado desde la clasificación: {rel}")
    check("inputs_sin_modificar", inputs_ok, "hashes de entradas iguales a summary.json")

    ok = not errors
    result = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "ok": ok,
        "errors_total": len(errors),
        "errors": errors[:100],
        "checks_total": len(checks),
        "checks_failed": sum(1 for c in checks if not c["ok"]),
        "genomes_ok": sum(
            1
            for g, r in genome_results.items()
            if not any(e.startswith(g + "/") for e in errors)
        ),
        "genomes": genome_results,
        "totals": dict(recomputed_totals),
        "checks": checks,
    }
    common.write_json(common.OUT / "validation.json", result)

    print(f"[05] checks={len(checks)} fallos={result['checks_failed']} errores={len(errors)}")
    for e in errors[:20]:
        print(f"  ERROR {e}")
    print(f"[05] ok={ok} -> {common.OUT / 'validation.json'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
