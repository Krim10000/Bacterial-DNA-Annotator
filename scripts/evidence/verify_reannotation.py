#!/usr/bin/env python3
"""verify_reannotation.py -- verificación de las salidas de reannotation-builder.

Comprueba, sin escribir nada:
  1. `reannotation.tsv`: 10 columnas, 77 filas, claves únicas e iguales a los 77
     accionables de `gap_review.tsv`.
  2. decisiones válidas, campos obligatorios completos, sin tabs/saltos.
  3. `candidate` no vacío => el ORF existe en el proteoma del genoma y tiene un
     hit del dominio del caso en `domain_hits.tsv`.
  4. `summary.json` coherente con el TSV (conteos y clases) y con rutas citadas.
  5. `README.md` presente y con las secciones de método, límites y hallazgos.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "metabolism" / "reannotation"
TSV = OUT / "reannotation.tsv"
SUM = OUT / "summary.json"
README = OUT / "README.md"
REVIEW = ROOT / "data" / "metabolism" / "review" / "gap_review.tsv"
NOHOM = ROOT / "data" / "metabolism" / "NO_HOMOLOG_DOMAINS.tsv"
HITS = OUT / "domain_hits.tsv"
ANNOT = ROOT / "data" / "annotations_ab_initio"
COLS = ["genome", "gap", "family", "support", "classification", "domain",
        "candidate", "decision", "reason", "evidence"]
DECISIONS = {"reannotate", "validate_domain", "keep_missing"}
ACTIONABLE = {"annotation_gap", "domain_present"}

errors: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        errors.append(name)


def main() -> int:
    for p in (TSV, SUM, README, HITS):
        check(f"existe {p.name}", p.exists())

    raw = TSV.read_text().splitlines()
    rows = list(csv.DictReader(open(TSV), delimiter="\t"))
    header = raw[0].split("\t") if raw else []
    check("columnas exactas", header == COLS, str(header))
    check("77 filas", len(rows) == 77, f"{len(rows)}")
    check("10 columnas por fila", all(len(l.split("\t")) == len(COLS) for l in raw[1:]),
          f"malformadas={sum(1 for l in raw[1:] if len(l.split(chr(9))) != len(COLS))}")

    review = {(r["genome"], r["gap"]) for r in csv.DictReader(open(REVIEW), delimiter="\t")
              if r["classification"] in ACTIONABLE}
    keys = [(r["genome"], r["gap"]) for r in rows]
    check("claves únicas", len(keys) == len(set(keys)),
          f"dup={len(keys) - len(set(keys))}")
    check("claves == 77 accionables", set(keys) == review,
          f"faltan={len(review - set(keys))} sobran={len(set(keys) - review)}")
    check("clasificaciones válidas", all(r["classification"] in ACTIONABLE for r in rows))
    check("decisiones válidas", all(r["decision"] in DECISIONS for r in rows),
          str(Counter(r["decision"] for r in rows)))
    check("sin campos vacíos obligatorios",
          all(str(r[k]).strip() for r in rows for k in
              ("genome", "gap", "family", "support", "classification", "domain",
               "decision", "reason", "evidence")))
    check("sin tabs/saltos en reason/evidence",
          all("\t" not in r["reason"] and "\n" not in r["reason"]
              and "\t" not in r["evidence"] and "\n" not in r["evidence"] for r in rows))
    check("support numérico", all(r["support"].isdigit() for r in rows))

    nohom = {(r["genome"], r["gap"]): r for r in
             csv.DictReader(open(NOHOM), delimiter="\t")}
    check("domain == arquitectura de NO_HOMOLOG_DOMAINS",
          all(r["domain"] == nohom[(r["genome"], r["gap"])]["domains"] for r in rows))

    hits = {(r["genome"], r["orf"], r["domain"]) for r in
            csv.DictReader(open(HITS), delimiter="\t")}
    bad_cand, bad_hit = [], []
    for r in rows:
        if not r["candidate"]:
            if r["decision"] in ("reannotate", "validate_domain"):
                bad_cand.append((r["genome"], r["gap"]))
            continue
        faa = ANNOT / r["genome"] / "proteins.faa"
        found = False
        if faa.exists():
            for line in open(faa):
                if line.startswith(">") and line[1:].strip() == r["candidate"]:
                    found = True
                    break
        if not found:
            bad_cand.append((r["genome"], r["gap"], r["candidate"]))
        exp = set(nohom[(r["genome"], r["gap"])]["domains"].split(";"))
        exp |= {"AdoMet_dc"}   # modelo extra del productor para dcAdoMet
        if not any((r["genome"], r["candidate"], d) in hits for d in exp):
            bad_hit.append((r["genome"], r["gap"], r["candidate"]))
    check("candidato existe en el proteoma", not bad_cand, str(bad_cand[:3]))
    check("candidato tiene hit de dominio del caso", not bad_hit, str(bad_hit[:3]))
    check("keep_missing sin candidato",
          all(not r["candidate"] for r in rows if r["decision"] == "keep_missing"))

    summary = json.load(open(SUM))
    dec = Counter(r["decision"] for r in rows)
    check("summary total == 77", summary["counts"]["total"] == 77)
    check("summary by_decision == TSV",
          summary["counts"]["by_decision"] == dict(dec), str(dec))
    cls = Counter(r["classification"] for r in rows)
    check("summary by_classification == TSV",
          summary["counts"]["by_classification"] == dict(cls), str(cls))
    check("summary validation sin fallos",
          summary["validation"]["expected_77"]
          and summary["validation"]["duplicates"] == 0
          and summary["validation"]["invalid_decisions"] == 0
          and summary["validation"]["empty_fields"] == 0)
    check("summary 3 decisiones documentadas",
          set(summary["decision_policy"]) == DECISIONS)
    check("summary no_confirmed", summary.get("status") == "proposed_not_confirmed")
    check("summary cita base intacta",
          summary["source_of_truth"]["base_table_not_modified"]
          == "data/metabolism/gap_interpretation.tsv")

    text = README.read_text()
    for sec in ("## Método", "## Límites", "## Hallazgos"):
        check(f"README sección {sec}", sec in text)
    check("README menciona 77", "77" in text)

    print()
    if errors:
        print(f"FALLARON {len(errors)} checks: {errors}")
        return 1
    print("OK: 77 casos, 10 columnas, decisiones válidas, candidatos verificados, "
          "summary y README consistentes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
