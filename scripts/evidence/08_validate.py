#!/usr/bin/env python3
"""08_validate.py — verifica el contrato y que los estados no inventen evidencia.

Comprobaciones:
  * esquema exacto y 500 filas únicas;
  * `status` = `homolog+dominio` con vocabulario permitido;
  * `homolog_present`/`homolog_weak` implican métricas numéricas y umbrales;
  * `best_acc` es una accesión productora no eucariota de ESA brecha;
  * `domain_present` implica un dominio esperado del gap realmente detectado;
  * `no_domain` implica panel completo y ausencia real (no falta de dato);
  * `not_checked` no lleva métricas ni dominios inventados;
  * support/family/productores coinciden con las tablas base;
  * las 77 accionables no están;
  * `summary.json` parsea, sus conteos coinciden y el hash de la tabla es válido.
Devuelve código 1 si algo falla.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ACTIONABLE, BATCH, EVID, GENERIC, KEY, MAX_EVALUE, R2A,
                    MIN_PIDENT_STRONG, MIN_PIDENT_WEAK, MIN_QCOV_STRONG,
                    MIN_QCOV_WEAK, README_OUT, REAL_GAPS, REFINED, REVIEW,
                    SUMMARY_OUT, TABLE_OUT, WORK, read_tsv)

FIELDS = ["genome", "gap", "family", "support", "producers", "best_acc",
          "pident", "qcov", "evalue", "domains", "status", "evidence"]
HOMOLOG = {"homolog_present", "homolog_weak", "homolog_low",
           "no_homolog", "not_checked"}
DOMAIN = {"domain_present", "no_domain", "not_checked"}

ERRORS: list[str] = []
CHECKS = 0


def check(cond: bool, msg: str) -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        ERRORS.append(msg)


def main() -> int:
    table = TABLE_OUT
    check(table.exists(), f"falta {table.name}")
    if not table.exists():
        print(f"FALLA: falta {table.name}")
        return 1
    with open(table, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        raw = list(reader)
    check(header == FIELDS, f"columnas: {header}")
    rows = [dict(zip(FIELDS, r)) for r in raw]
    selected = {(r["genome"], r["gap"]) for r in read_tsv(WORK / "selected_gaps.tsv")}
    check(len(rows) == len(selected), f"filas: {len(rows)} vs lote {len(selected)}")
    check({(r["genome"], r["gap"]) for r in rows} == selected,
          "la tabla no cubre exactamente el lote seleccionado")
    check(len({(r["genome"], r["gap"]) for r in rows}) == len(rows),
          "claves duplicadas")
    check(all(len(r) == len(FIELDS) for r in raw), "filas con columnas de más/menos")

    review = {(r["genome"], r["gap"]): r for r in read_tsv(REVIEW)}
    real = {(r["genome"], r["gap"]): r for r in read_tsv(REAL_GAPS)}
    refined = {(r["genome"], r["gap"]): r for r in read_tsv(REFINED)}
    acc_info = {r["acc"]: r for r in read_tsv(WORK / "acc_info.tsv")}
    gap_dom = {(r["genome"], r["gap"]): r for r in read_tsv(WORK / "gap_domains.tsv")}
    hom = {(r["genome"], r["gap"]): r for r in read_tsv(WORK / "homology_hits.tsv")}
    det = {(r["genome"], r["pfam_acc"]) for r in read_tsv(WORK / "domain_hits.tsv")}
    panel = {r["requested_acc"] for r in read_tsv(WORK / "panel_meta.tsv")}
    name_by_acc = {r["requested_acc"]: r["name"]
                   for r in read_tsv(WORK / "panel_meta.tsv")}
    generic_accs = {a for a, n in name_by_acc.items() if n in GENERIC}

    scanned_dom = {x.name[:-7] for x in (WORK / "domtbl").glob("*.tsv.gz")}
    scanned_hom = {x.name[:-7] for x in (WORK / "diamond").glob("*.tsv.gz")}
    actionable = {k for k, r in review.items() if r["classification"] in ACTIONABLE}
    check(not ({k for k in [(r["genome"], r["gap"]) for r in rows]} & actionable),
          "hay filas accionables en la evidencia")

    status_c: Counter = Counter()
    for r in rows:
        key = (r["genome"], r["gap"])
        tag = f"{key}"
        status_c[r["status"]] += 1
        check(key in review, f"{tag}: clave fuera de gap_review")
        if key not in review:
            continue
        check(review[key]["classification"] == "unknown",
              f"{tag}: no era unknown")
        if key in real:
            check(int(r["support"]) == int(real[key]["support_400"]),
                  f"{tag}: support no coincide")
        check(r["family"] == review[key]["family"], f"{tag}: family no coincide")

        # productores exactamente los del refined
        if key in refined:
            want = {p for p in refined[key]["producers"].split(";") if p}
            got = {p for p in r["producers"].split(";") if p}
            check(got == want, f"{tag}: productores no coinciden con refined")

        parts = r["status"].split("+")
        check(len(parts) == 2, f"{tag}: status mal formado {r['status']}")
        h_state, d_state = (parts + ["", ""])[:2]
        check(h_state in HOMOLOG, f"{tag}: homolog status inválido {h_state}")
        check(d_state in DOMAIN, f"{tag}: domain status inválido {d_state}")

        # --- homología ------------------------------------------------------
        if h_state != "not_checked":
            check(r["genome"] in scanned_hom, f"{tag}: genoma sin escaneo de homología")
        if h_state in {"homolog_present", "homolog_weak", "homolog_low"}:
            check(bool(r["best_acc"]), f"{tag}: homolog sin best_acc")
            try:
                pid, qcov, ev = float(r["pident"]), float(r["qcov"]), float(r["evalue"])
            except ValueError:
                check(False, f"{tag}: métricas no numéricas")
                pid = qcov = ev = -1
            if h_state == "homolog_low":
                # hay hit, pero no alcanza el umbral débil
                check(not (pid >= MIN_PIDENT_WEAK and qcov >= MIN_QCOV_WEAK
                           and ev <= MAX_EVALUE),
                      f"{tag}: homolog_low que sí pasa el umbral débil")
            else:
                check(ev <= MAX_EVALUE, f"{tag}: evalue {ev} > {MAX_EVALUE}")
                check(pid >= MIN_PIDENT_WEAK, f"{tag}: pident {pid} < débil")
                check(qcov >= MIN_QCOV_WEAK, f"{tag}: qcov {qcov} < débil")
            if h_state == "homolog_present":
                check(pid >= MIN_PIDENT_STRONG, f"{tag}: pident {pid} < fuerte")
                check(qcov >= MIN_QCOV_STRONG, f"{tag}: qcov {qcov} < fuerte")
            # best_acc debe ser productor de la brecha y no eucariota
            accs = {a for a in gap_accs(r, acc_info) if a}
            check(r["best_acc"] in accs,
                  f"{tag}: best_acc {r['best_acc']} no es productor")
            check(acc_info.get(r["best_acc"], {}).get("euk") == "0",
                  f"{tag}: best_acc eucariota")
            # coherencia con la tabla de homología
            check(r["best_acc"] == hom[key]["best_acc"],
                  f"{tag}: best_acc difiere de homology_hits")
        else:
            check(not r["best_acc"], f"{tag}: {h_state} con best_acc")
            check(not r["pident"] and not r["qcov"] and not r["evalue"],
                  f"{tag}: {h_state} con métricas")

        # --- dominio --------------------------------------------------------
        expected = [x for x in gap_dom[key]["domains"].split(";") if x]
        specific = [x for x in expected if x not in generic_accs]
        if d_state in {"domain_present", "no_domain"}:
            check(r["genome"] in scanned_dom, f"{tag}: genoma sin escaneo de dominios")
        if d_state == "domain_present":
            check(bool(r["domains"]), f"{tag}: domain_present sin dominios")
            for token in r["domains"].split(";"):
                check("(" in token, f"{tag}: dominio mal formado {token}")
                acc = token.split("(")[1].split(",")[0] if "(" in token else ""
                check(acc in expected, f"{tag}: dominio {acc} no esperado")
                check(acc in panel, f"{tag}: dominio {acc} fuera del panel")
                check((r["genome"], acc) in det,
                      f"{tag}: dominio {acc} no está en domain_hits")
        elif d_state == "no_domain":
            check(r["domains"] == "none", f"{tag}: no_domain sin 'none'")
            check(gap_dom[key]["status"] == "ok",
                  f"{tag}: no_domain sin base de dominios")
            check(bool(specific), f"{tag}: no_domain sin dominios específicos")
            check(set(expected) <= panel, f"{tag}: no_domain con panel incompleto")
            check(not any((r["genome"], a) in det for a in specific),
                  f"{tag}: no_domain pero hay dominio detectado")
        else:  # not_checked
            check(not r["domains"], f"{tag}: not_checked con dominios")
            base = gap_dom[key]
            ok_not_checked = (base["status"] != "ok" or not specific
                              or bool(set(expected) - panel))
            check(ok_not_checked, f"{tag}: not_checked sin justificación")

        # evidence no vacío y sin tabs
        check(bool(r["evidence"]), f"{tag}: evidence vacío")
        check("\t" not in r["evidence"] and "\n" not in r["evidence"],
              f"{tag}: evidence con tab/salto")

    # --- summary ---------------------------------------------------------------
    sp = SUMMARY_OUT
    check(sp.exists(), f"falta {sp.name}")
    summary = json.loads(sp.read_text())
    ec = summary["evidence_counts"]
    check(ec["status"] == dict(status_c), "summary.status no coincide")
    check(sum(ec["status"].values()) == len(rows),
          "summary.status no suma el total de filas")
    check(summary["outputs"]["table_sha256"] == hashlib.sha256(table.read_bytes()).hexdigest(),
          "hash de la tabla no coincide con summary")

    # --- no_bacterial_reference (lote 03) ------------------------------------------------------
    if BATCH in ("03", "04"):
        ep = EVID / f"no_bacterial_reference_batch_{BATCH}.tsv"
        check(ep.exists(), f"falta no_bacterial_reference_batch_{BATCH}.tsv")
        if ep.exists():
            er = read_tsv(ep)
            expected_euk = set()
            for r in rows:
                key = (r["genome"], r["gap"])
                accs = REFINED_ROWS.get(key, [])
                prok = [a for a in accs if acc_info.get(a, {}).get("euk") == "0"]
                if not prok:
                    expected_euk.add(key)
            got_euk = {(x["genome"], x["gap"]) for x in er}
            check(got_euk == expected_euk,
                  f"no_bacterial_reference: {len(got_euk)} vs esperadas {len(expected_euk)}")
            check(len(er) == len(expected_euk), "no_bacterial_reference con duplicados")
            for x in er:
                check(x["status"] == "not_checked+not_checked",
                      f"{x['genome']},{x['gap']}: no_bacterial_reference con evidencia")
                check("no es ausencia" in x["note"],
                      f"{x['genome']},{x['gap']}: no_bacterial_reference sin nota anti-ausencia")

    # --- README ----------------------------------------------------------------
    rd = README_OUT.read_text().lower()
    for needle in ["lote", "genomas", "método", "límites", "3 hallazgos"]:
        check(needle in rd, f"README sin '{needle}'")

    print(f"comprobaciones: {CHECKS}")
    if ERRORS:
        print(f"FALLOS: {len(ERRORS)}")
        for e in ERRORS[:40]:
            print("  -", e)
        return 1
    print("OK: contrato y estados verificados")
    return 0


def gap_accs(row: dict, acc_info: dict) -> list[str]:
    """Accesiones no eucariotas del gap (reconstruidas del refined)."""
    key = (row["genome"], row["gap"])
    ref = REFINED_ROWS[key]
    out = []
    for acc in ref:
        if acc_info.get(acc, {}).get("euk") == "0":
            out.append(acc)
    return out


REFINED_ROWS: dict[tuple, list[str]] = {}


def load_refined_accs() -> None:
    r2a = read_tsv(R2A)
    rid: dict[str, set[str]] = {}
    for r in r2a:
        if r["RHEA_ID"].isdigit() and r["ID"]:
            rid.setdefault("RHEA:" + r["RHEA_ID"], set()).add(r["ID"])
    for r in read_tsv(REFINED):
        key = (r["genome"], r["gap"])
        accs: set[str] = set()
        for p in r["producers"].split(";"):
            accs |= rid.get(p.strip(), set())
        REFINED_ROWS[key] = sorted(accs)


if __name__ == "__main__":
    load_refined_accs()
    raise SystemExit(main())
