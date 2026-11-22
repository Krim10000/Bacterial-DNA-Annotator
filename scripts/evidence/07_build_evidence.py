#!/usr/bin/env python3
"""07_build_evidence.py — cruza homología + dominios y escribe la evidencia.

Por cada una de las 500 brechas del lote produce una fila con:

  * `best_acc`, `pident`, `qcov`, `evalue` del mejor hit a cualquier enzima
    productora no eucariota (DIAMOND all-hits);
  * `domains` con los dominios Pfam detectados que estaban entre los esperados
    de la brecha (hmmscan --cut_ga);
  * `status` de dos ejes, para no mezclar evidencia ausente con no evaluada:
        homolog_present | homolog_weak | no_homolog | not_checked
      + domain_present | no_domain | not_checked
    (p. ej. `homolog_present+domain_present`, `no_homolog+no_domain`);
  * `evidence` con el detalle legible por máquina.

Nada se promueve sin base: si faltan accesiones, dominios o proteoma, el eje
correspondiente queda `not_checked` y el `evidence` explica por qué.

Salidas:
    data/metabolism/evidence/gap_evidence.tsv
    data/metabolism/evidence/summary.json
    data/metabolism/evidence/README.md
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ACTIONABLE, BATCH, EVID, GENERIC, MAX_EVALUE,
                    MIN_PIDENT_STRONG, MIN_PIDENT_WEAK, MIN_QCOV_STRONG,
                    MIN_QCOV_WEAK, README_OUT, REVIEW, SUMMARY_OUT, TABLE_OUT,
                    WORK, clean, read_tsv, write_json, write_tsv)

FIELDS = ["genome", "gap", "family", "support", "producers", "best_acc",
          "pident", "qcov", "evalue", "domains", "status", "evidence"]
ALLOWED_HOMOLOG = {"homolog_present", "homolog_weak", "homolog_low",
                   "no_homolog", "not_checked"}
ALLOWED_DOMAIN = {"domain_present", "no_domain", "not_checked"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        lines = (out.stdout or out.stderr).strip().splitlines()
        for line in lines:
            line = line.lstrip("# ").strip()
            if line.startswith("HMMER") or "diamond version" in line:
                return line
        return lines[0] if lines else ""
    except Exception:
        return ""


def fmt_evalue(value: float) -> str:
    return f"{value:.2e}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", default="gap_evidence_batch_01")
    a = ap.parse_args()

    sel = read_tsv(WORK / "selected_gaps.tsv")
    dom_expect = {(r["genome"], r["gap"]): r for r in read_tsv(WORK / "gap_domains.tsv")}
    hom = {(r["genome"], r["gap"]): r for r in read_tsv(WORK / "homology_hits.tsv")}
    dom_hits = read_tsv(WORK / "domain_hits.tsv")
    panel = read_tsv(WORK / "panel_meta.tsv")
    acc_info = {r["acc"]: r for r in read_tsv(WORK / "acc_info.tsv")}
    universe = (json.loads((WORK / "universe.json").read_text())
                if (WORK / "universe.json").exists() else None)

    scanned_dom = {x.name[:-7] for x in (WORK / "domtbl").glob("*.tsv.gz")}
    scanned_hom = {x.name[:-7] for x in (WORK / "diamond").glob("*.tsv.gz")}
    panel_accs = {r["requested_acc"] for r in panel}
    name_by_acc = {r["requested_acc"]: r["name"] for r in panel}
    generic_accs = {acc for acc, name in name_by_acc.items() if name in GENERIC}

    # dominios detectados por genoma -> acc -> mejor hit
    detected: dict[tuple[str, str], dict] = {}
    for r in dom_hits:
        key = (r["genome"], r["pfam_acc"])
        ev = float(r["i_evalue"])
        prev = detected.get(key)
        if prev is None or ev < prev["i_evalue"]:
            detected[key] = {"i_evalue": ev, "source": r["source"],
                             "protein_id": r["protein_id"],
                             "pfam_name": r["pfam_name"],
                             "qcov_dom": float(r["qcov_dom"])}

    rows = []
    for s in sel:
        key = (s["genome"], s["gap"])
        h = hom.get(key, {})
        d = dom_expect.get(key, {})
        producers = [p for p in s["producers"].split(";") if p]
        accs = [x for x in s["producer_accs"].split(";") if x]
        prok_accs = [x for x in accs if acc_info.get(x, {}).get("euk") == "0"]
        euk_accs = [x for x in accs if x not in prok_accs]

        # --- eje homología ---------------------------------------------------
        homolog = h.get("status", "not_checked")
        if homolog not in ALLOWED_HOMOLOG:
            homolog = "not_checked"
        if s["genome"] not in scanned_hom:
            homolog, h = "not_checked", {}
        elif not prok_accs:
            homolog = "not_checked"
            h = {}
        # verificación anti-invención del acc reportado
        if h.get("best_acc") and h["best_acc"] not in prok_accs:
            raise AssertionError(f"best_acc {h['best_acc']} no es productor de {key}")

        # --- eje dominio -----------------------------------------------------
        expected = [x for x in d.get("domains", "").split(";") if x]
        panel_missing = [x for x in expected if x not in panel_accs]
        specific = [x for x in expected if x not in generic_accs]
        found = []
        for acc in specific:
            hit = detected.get((s["genome"], acc))
            if hit:
                found.append((acc, hit))
        found.sort(key=lambda x: x[1]["i_evalue"])
        generic_found = [x for x in expected
                         if x in generic_accs and (s["genome"], x) in detected]
        if s["genome"] not in scanned_dom:
            domain = "not_checked"
            why_domain = "proteoma no escaneado con el panel"
        elif not expected or d.get("status") != "ok":
            domain = "not_checked"
            why_domain = d.get("reason", "sin dominio Pfam esperado")
        elif not specific:
            domain = "not_checked"
            why_domain = "solo dominios genéricos; no sostienen dominio presente"
        elif found:
            domain = "domain_present"
            why_domain = f"{len(found)} dominio(s) específico(s) detectado(s)"
        elif panel_missing:
            domain = "not_checked"
            why_domain = f"{len(panel_missing)} dominio(s) sin HMM en el panel"
        else:
            domain = "no_domain"
            why_domain = "esperados evaluados y ausentes"

        # --- columnas --------------------------------------------------------
        if found:
            dom_txt = ";".join(
                f"{hit['pfam_name']}({acc},i={fmt_evalue(hit['i_evalue'])},"
                f"{hit['source']})" for acc, hit in found)
        elif domain == "no_domain":
            dom_txt = "none"
        else:
            dom_txt = ""

        if homolog == "not_checked":
            info_h = f"homologo=not_checked ({h.get('note', 'sin base')})"
        elif h.get("best_acc"):
            org = acc_info.get(h["best_acc"], {}).get("organism", "")
            info_h = (f"homologo={homolog} (acc={h['best_acc']}; "
                      f"pident={h['pident']}; qcov={h['qcov']}; "
                      f"evalue={h['evalue']}; fuente={h['source']}; "
                      f"accs_hit={h['n_accs_hit']}/{len(prok_accs)}; "
                      f"n_hits={h['n_hits']}; organismo={org[:60]})")
        else:
            info_h = f"homologo={homolog} (sin hit; accs={len(prok_accs)})"
        info_d = (f"dominio={domain} ({why_domain}; "
                  f"esperados={'yes' if expected else 'none'}"
                  f"={','.join(expected) if expected else '-'}; "
                  f"panel={len(expected) - len(panel_missing)}/{len(expected)}; "
                  f"politica={d.get('domain_policy', '-')}; "
                  f"genericos_detectados={len(generic_found)})")
        info_p = (f"productores={len(producers)} RHEA; "
                  f"accs={len(accs)} ({len(prok_accs)} no eucariotas)")
        evidence = " | ".join([info_h, info_d, info_p, f"lote={a.batch}"])

        rows.append({
            "genome": s["genome"], "gap": s["gap"], "family": s["family"],
            "support": s["support"],
            "producers": ";".join(producers),
            "best_acc": h.get("best_acc", ""),
            "pident": h.get("pident", ""), "qcov": h.get("qcov", ""),
            "evalue": h.get("evalue", ""),
            "domains": dom_txt,
            "status": f"{homolog}+{domain}",
            "evidence": clean(evidence),
            "_homolog": homolog, "_domain": domain,
            "_prok": bool(prok_accs), "_euk_accs": euk_accs,
        })

    rows.sort(key=lambda r: (-int(r["support"]), r["genome"], r["gap"]))
    out = TABLE_OUT
    write_tsv(out, rows, FIELDS)
    euk_rows = [{
        "genome": r["genome"], "gap": r["gap"], "family": r["family"],
        "support": r["support"], "producers": r["producers"],
        "accs_total": len(r["_euk_accs"]),
        "accs_reference": len(r["_euk_accs"]),
        "reference_sample": ";".join(r["_euk_accs"][:5]),
        "status": r["status"],
        "note": ("sin productor bacteriano curado en Rhea; no es ausencia; "
                 "no es genoma eucariota"),
    } for r in rows if not r["_prok"]]
    if BATCH in ("03", "04"):
        write_tsv(EVID / f"no_bacterial_reference_batch_{BATCH}.tsv", euk_rows,
                  ["genome", "gap", "family", "support", "producers",
                   "accs_total", "accs_reference", "reference_sample",
                   "status", "note"])

    # --- summary -------------------------------------------------------------
    hom_c = Counter(r["_homolog"] for r in rows)
    dom_c = Counter(r["_domain"] for r in rows)
    status_c = Counter(r["status"] for r in rows)
    strong = [r for r in rows if r["_homolog"] == "homolog_present"]
    weak = [r for r in rows if r["_homolog"] == "homolog_weak"]
    low = [r for r in rows if r["_homolog"] == "homolog_low"]
    both = [r for r in rows if r["_homolog"] == "homolog_present"
            and r["_domain"] == "domain_present"]
    neither = [r for r in rows if r["_homolog"] == "no_homolog"
               and r["_domain"] == "no_domain"]
    any_positive = [r for r in rows
                    if r["best_acc"] or r["domains"] not in {"", "none"}]

    review = {(r["genome"], r["gap"]): r for r in read_tsv(REVIEW)}
    prev_homolog = Counter(review[(r["genome"], r["gap"])]["homolog_status"]
                           for r in rows)
    rescued = [r for r in rows
               if review[(r["genome"], r["gap"])]["homolog_status"] == "not_checked"
               and r["_homolog"] in {"homolog_present", "homolog_weak"}]
    src_hom, src_hom_strong, src_dom = Counter(), Counter(), Counter()
    src_dom_rows: Counter = Counter()
    for r in rows:
        m = re.search(r"fuente=(\w+)", r["evidence"])
        if r["best_acc"] and m:
            src_hom[m.group(1)] += 1
            if r["_homolog"] == "homolog_present":
                src_hom_strong[m.group(1)] += 1
        combo = set()
        for tok in r["domains"].split(";"):
            mm = re.search(r",(ncbi|abinitio)\)$", tok)
            if mm:
                src_dom[mm.group(1)] += 1
                combo.add(mm.group(1))
        for x in combo:
            src_dom_rows[x] += 1
    sel_json = json.loads((WORK / "selection.json").read_text())
    hom_json = json.loads((WORK / "homology.json").read_text())
    dom_json = json.loads((WORK / "domains.json").read_text())
    panel_json = json.loads((WORK / "panel_meta.json").read_text())

    def med(vals: list[float]) -> float | None:
        if not vals:
            return None
        vals = sorted(vals)
        return round(vals[len(vals) // 2], 1)

    sup = sorted(int(r["support"]) for r in rows)
    db_json = json.loads((WORK / "producer_db.json").read_text())
    euk_excluded = db_json["euk_accs_excluded"]
    panel_models_n = panel_json["panel_models"]
    prio = dict(Counter(s["priority"] for s in sel))
    n_metab = len({r["gap"] for r in rows})
    n_fam = len({r["family"] for r in rows})
    summary = {
        "batch": a.batch,
        "universe": universe,
        "no_bacterial_reference": len(euk_rows),
        "generated_at": None,  # se rellena abajo con la fecha del sistema
        "batch_scope": {
            "rows": len(rows),
            "genomes": len({r["genome"] for r in rows}),
            "metabolites": len({r["gap"] for r in rows}),
            "families": len({r["family"] for r in rows}),
            "priorities": prio,
            "support_min": sup[0],
            "support_median": sup[len(sup) // 2],
            "support_max": sup[-1],
            "excluded_actionable": sum(
                1 for r in review.values() if r["classification"] in ACTIONABLE),
            "selection_rule": sel_json.get("rule", ""),
            "unknown_total": sel_json.get("unknown_total"),
            "excluded_previous_batch": sel_json.get("excluded_previous_batch"),
        },
        "methods": {
            "homology": {
                "tool": tool_version(["diamond", "version"]),
                "db": "productores no eucariotas del lote (work/producer_db.dmnd)",
                "db_sequences": db_json["prok_accs_in_db"],
                "proteomes": ["ncbi (data/annotations/<genome>.faa.gz)",
                              "abinitio (data/annotations_ab_initio/<genome>/proteins.faa)"],
                "mode": "blastp all-hits, evalue<=1e-3, max-target-seqs=500",
                "thresholds": hom_json["thresholds"],
            },
            "domains": {
                "tool": tool_version(["hmmscan", "-h"]),
                "cut": "hmmscan --cut_ga",
                "panel_models": panel_json["panel_models"],
                "panel_families_requested": panel_json["families_requested"],
                "panel_families_failed": panel_json["families_failed"],
                "domain_source": "UniProt REST xref_pfam por accesión productora",
                "policy": "core > majority (>=60%) > any, como el resto del proyecto",
            },
            "taxonomy": {
                "source": "data/taxdump (ancestro 2759 Eukaryota)",
                "euk_accs_excluded": euk_excluded,
            },
        },
        "evidence_counts": {
            "homolog_axis": dict(hom_c),
            "domain_axis": dict(dom_c),
            "status": dict(status_c),
            "rows_with_any_positive": len(any_positive),
            "homolog_present": len(strong),
            "homolog_weak": len(weak),
            "homolog_low": len(low),
            "homolog_and_domain": len(both),
            "neither_homolog_nor_domain": len(neither),
            "not_checked_domain": dom_c.get("not_checked", 0),
            "not_checked_homolog": hom_c.get("not_checked", 0),
            "pident_median_strong": med([float(r["pident"]) for r in strong]),
            "pident_median_weak": med([float(r["pident"]) for r in weak]),
        },
        "gaps_still_unknown_after": summary_gaps_unknown(rows),
        "breakdowns": {
            "by_family": _by_family(rows),
            "top_metabolites_with_evidence": _by_metabolite(rows)[:12],
            "metabolites_with_any_row_without_positive_evidence": _no_evidence_metabolites(rows),
        },
        "source_breakdown": {
            "homolog_best_by_source": dict(src_hom),
            "homolog_strong_by_source": dict(src_hom_strong),
            "domain_tokens_by_source": dict(src_dom),
            "domain_rows_by_source": dict(src_dom_rows),
            "note": "abinitio = ORFs predichos sin anotar; fuente declarada en evidence",
        },
        "base_alignment": {
            "review_rows": sel_json["input_rows"]["review"],
            "queue_rows": sel_json["input_rows"]["queue"],
            "real_rows": sel_json["input_rows"]["real"],
            "common_keys": sel_json["input_rows"]["common_keys"],
            "review_keys_missing_in_base": sel_json["input_rows"]["excluded_missing_base"],
            "note": (
                "universo alineado: review = queue = real"
                if (sel_json["input_rows"]["review"]
                    == sel_json["input_rows"]["common_keys"])
                else "las tablas base difieren; el lote usa la intersección "
                     "de claves de review/queue/real"),
        },
        "comparison_with_gap_review": {
            "previous_homolog_status_in_batch": dict(prev_homolog),
            "rescued_from_not_checked": len(rescued),
        },
        "findings": [
            {
                "id": 1,
                "title": "La homología dirigida rescata brechas que gap_review dejaba en unknown",
                "metric": len(rescued),
                "statement": (
                    f"{len(rescued)} filas que gap_review marcaba homolog_status="
                    f"not_checked ahora tienen homólogo fuerte o débil; "
                    f"{len(strong)} filas del lote tienen homología fuerte. "
                    "No es ausencia de enzima: es reacción no mapeada."),
            },
            {
                "id": 2,
                "title": "Los dominios Pfam añaden señal donde la identidad no llega",
                "metric": dom_c.get("domain_present", 0),
                "statement": (
                    f"{dom_c.get('domain_present', 0)} filas tienen un dominio Pfam "
                    "específico del productor; "
                    f"{len(both)} combinan homología fuerte y dominio. "
                    "Son candidatas a homología divergente o anotación incompleta, "
                    "no a pérdida de la enzima."),
            },
            {
                "id": 3,
                "title": "El bloque sin evidencia positiva es el duro restante",
                "metric": len(neither),
                "statement": (
                    f"{len(neither)} filas no muestran ni homólogo ni dominio; "
                    f"{hom_c.get('not_checked', 0)} no son evaluables (sin productor "
                    f"no eucariota). En total, {summary_gaps_unknown(rows)} de las "
                    f"{len(rows)} quedan sin evidencia positiva. Ahí cabe hablar de "
                    "posible ausencia/auxotrofía, pero el lote no la confirma."),
            },
        ],
        "limits": [
            f"Lote de {len(rows)} filas (no las "
            f"{sel_json.get('unknown_total', '?')} unknown del universo): "
            f"{prio.get('alta', 0)} alta + {prio.get('media', 0)} media.",
            f"Las {len(rows)} filas corresponden a {n_metab} metabolitos y "
            f"{n_fam} familias; no es una muestra aleatoria.",
            "Homología y dominio son evidencia de presencia, no de función: no cierran la brecha.",
            f"Productores eucariotas excluidos de la homología y de los dominios "
            f"esperados ({euk_excluded} accesiones).",
            "Los ORFs ab_initio pueden ser espurios; la fuente del hit queda declarada.",
            f"El panel cubre solo las familias Pfam pedidas por el lote "
            f"({panel_models_n}); no es Pfam-A completo.",
            "Sin validación experimental ni contexto genómico (operón, dirección, regulación).",
            "Homólogo y dominio pueden venir de proteínas distintas: no se exige co-localización.",
            "`not_checked` no es ausencia: es falta de base para ese eje.",
            (f"{hom_c.get('not_checked', 0)} filas no tienen ningún productor no "
             "eucariota en Rhea (solo eucariotas): quedan fuera de la homología y "
             "del dominio, concentradas en metabolismo lipídico."),
            (f"Universo: {universe['enzymatic_missing']} enzymatic_missing; "
             f"{universe['excluded_union']} claves excluidas (cubiertas 01+02, "
             f"accionables, promovidas, reannotated); pool {universe['pool']}."
             if universe else
             f"Universo vigente: {sel_json['input_rows']['review']} filas de revisión, "
             f"{sel_json['input_rows']['common_keys']} claves comunes con "
             f"real_gaps/REVIEW_QUEUE; se excluyen "
             f"{sel_json['excluded_previous_batch']} claves del lote previo y "
             f"{sel_json['excluded_actionable']} accionables."),
        ],
        "outputs": {
            "table": str(out.relative_to(out.parents[2])),
            "table_sha256": sha256(out),
            "work_dir": "data/metabolism/evidence/work/",
        },
    }
    if universe:
        summary["batch_scope"]["excluded_actionable"] = universe["actionable_total"]
        summary["batch_scope"]["excluded_keys"] = universe["excluded_union"]
        summary["limits"].append(
            f"Universo lote {BATCH}: {universe['enzymatic_missing']} "
            f"enzymatic_missing; pool restante {universe['pool']} tras "
            f"excluir {universe['excluded_union']} claves (no ~1.927: las "
            "promovidas y reannotated ya estaban en las otras listas).")
    from datetime import datetime, timezone
    summary["generated_at"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    write_json(SUMMARY_OUT, summary)

    write_readme(summary, rows)
    print(json.dumps(summary["evidence_counts"], indent=2, ensure_ascii=False))
    print(f"-> {out}")
    return 0


def _positive(r: dict) -> bool:
    return bool(r["best_acc"]) or r["domains"] not in {"", "none"}


def _by_family(rows: list[dict]) -> list[dict]:
    out = {}
    for r in rows:
        f = out.setdefault(r["family"], {"family": r["family"], "rows": 0,
                                         "homolog_present": 0,
                                         "domain_present": 0,
                                         "no_positive_evidence": 0})
        f["rows"] += 1
        f["homolog_present"] += r["_homolog"] == "homolog_present"
        f["domain_present"] += r["_domain"] == "domain_present"
        f["no_positive_evidence"] += not _positive(r)
    return sorted(out.values(), key=lambda x: (-x["rows"], x["family"]))


def _by_metabolite(rows: list[dict]) -> list[dict]:
    out = {}
    for r in rows:
        m = out.setdefault(r["gap"], {"gap": r["gap"], "rows": 0,
                                      "with_positive_evidence": 0,
                                      "homolog_present": 0,
                                      "domain_present": 0})
        m["rows"] += 1
        m["with_positive_evidence"] += _positive(r)
        m["homolog_present"] += r["_homolog"] == "homolog_present"
        m["domain_present"] += r["_domain"] == "domain_present"
    return sorted(out.values(), key=lambda x: (-x["with_positive_evidence"],
                                               -x["rows"], x["gap"]))


def _no_evidence_metabolites(rows: list[dict]) -> list[str]:
    return sorted({r["gap"] for r in rows if not _positive(r)})


def summary_gaps_unknown(rows: list[dict]) -> int:
    return len([r for r in rows
                if not r["best_acc"] and r["domains"] in {"", "none"}])


def write_readme(summary: dict, rows: list[dict]) -> None:
    c = summary["evidence_counts"]
    b = summary["batch_scope"]
    m = summary["methods"]
    u = summary.get("universe")
    if u:
        excl_line = (
            f"- Exclusiones: {u['covered_total']} cubiertas por lotes previos "
            f"({u['promoted_total']} ya promovidas), {u['actionable_total']} "
            f"accionables y {u['reannotated_total']} reannotated; unión "
            f"{u['excluded_union']} claves.")
    else:
        excl_line = (
            f"- Excluidas las {b['excluded_actionable']} accionables vigentes "
            "(`annotation_gap` + `domain_present`, de `reannotation-builder`) y "
            f"las {b.get('excluded_previous_batch', 0)} del lote previo.")

    def examples(cond, n=6):
        out = []
        for r in rows:
            if cond(r):
                out.append((r["genome"], r["gap"], r["best_acc"],
                            r["pident"], r["domains"], r["status"]))
            if len(out) >= n:
                break
        return out

    lines = [
        f"# Evidencia por brecha (lote {BATCH})",
        "",
        f"**Generado:** {summary['generated_at']} · **lote:** {summary['batch']}",
        "",
        "## Qué es",
        "",
        f"Evidencia dirigida para {b['rows']} brechas `unknown` de "
        f"{b.get('unknown_total', '?')} del universo:",
        f"{b['priorities'].get('alta', 0)} alta, "
        f"{b['priorities'].get('media', 0)} media y "
        f"{b['priorities'].get('baja', 0)} baja "
        f"(las del lote previo quedan excluidas).",
        "Cada fila busca dos evidencias independientes de que la enzima productora",
        "existe en el genoma aunque la reacción no esté mapeada:",
        "",
        "1. **Homología all-hits** (DIAMOND) de las proteínas del genoma contra las",
        "   enzimas productoras Rhea de esa brecha (todas las accesiones, sin recorte).",
        "2. **Dominios Pfam** (hmmscan `--cut_ga`) contra el panel de familias que",
        "   espera esa brecha, con política `core` > `majority` > `any`.",
        "",
        "No se modifica la tabla base, ni `webapp/`, ni `annotations_ab_initio/`.",
        "No hay commits.",
        "",
    ] + ([
        "**Cambio de universo respecto al lote 01:** 23 brechas `reannotated`",
        "salieron de la base y 12 T1 inválidos entraron como `unknown`; el",
        "universo pasó de 4.021 a 4.010 filas y las 42 altas quedaron cubiertas",
        "en el lote 01. Este lote es, por tanto, 100% prioridad media.",
    ] if BATCH == "02" else [
        "**Universo recompuesto (lote 03):** `gap_interpretation.tsv` tiene",
        "3.667 `enzymatic_missing`. Se excluyen las 1.500 cubiertas por los",
        "lotes 01+02 (343 ya promovidas a `homology_present`), 77 accionables",
        "y 23 `reannotated`; la unión es 1.577 claves y el pool real 2.456.",
        "Este lote toma 1.000: 25 alta y 975 media.",
    ] if BATCH == "03" else [
        "**Universo de cierre (lote 04):** 3.433 `enzymatic_missing`; se",
        "excluyen las 2.500 cubiertas por 01-03, 77 accionables, 23",
        "reannotated y las promovidas a `homology_present`. Pool real 1.456",
        "(1.510 menos 54 accionables). Se procesan todas; las brechas sin",
        "productor bacteriano curado se separan en",
        "`no_bacterial_reference_batch_04.tsv` y no se tratan como ausencia.",
    ] if BATCH == "04" else []) + [
        "",
        "## Lote",
        "",
        f"- {b['rows']} filas · {b['genomes']} genomas · {b['metabolites']} metabolitos · {b['families']} familias.",
        f"- Prioridad: {b['priorities'].get('alta', 0)} alta / "
        f"{b['priorities'].get('media', 0)} media / "
        f"{b['priorities'].get('baja', 0)} baja.",
        f"- Support: min {b['support_min']} · mediana {b['support_median']} · max {b['support_max']}.",
        *([f"- no_bacterial_reference: {summary['no_bacterial_reference']} filas sin productor bacteriano curado, "
           f"separadas en `no_bacterial_reference_batch_{BATCH}.tsv` (no son ausencia)."]
          if BATCH in ("03", "04") else []),
        excl_line,
        f"- Universo vigente: {summary['base_alignment']['review_rows']} filas, "
        f"{b.get('unknown_total', '?')} `unknown`; "
        + ("sin altas nuevas en el pool (cubiertas o promovidas)."
           if b['priorities'].get('alta', 0) == 0 else
           "las de prioridad alta entran completas."),
        "",
        "## Métodos",
        "",
        f"- DIAMOND: `{m['homology']['tool']}`; DB de {m['homology']['db_sequences']} productores no eucariotas.",
        f"- Búsqueda: {m['homology']['mode']}.",
        f"  Umbral fuerte: {m['homology']['thresholds']['strong']}; "
        f"débil: {m['homology']['thresholds']['weak']}.",
        "- Dos proteomas: NCBI anotada y ORFs ab_initio; la fuente del hit queda en `evidence`.",
        f"- Dominios: `{m['domains']['tool']}` con {m['domains']['panel_models']} modelos Pfam",
        f"  ({m['domains']['panel_families_requested']} pedidos, {len(m['domains']['panel_families_failed'])} fallidos).",
        f"- Taxonomía: {m['taxonomy']['source']}; {m['taxonomy']['euk_accs_excluded']} accesiones eucariotas excluidas.",
        "",
        "## Estados",
        "",
        "`status` tiene dos ejes para no confundir *ausencia* con *no evaluado*:",
        "",
        "| eje | valores |",
        "|---|---|",
        "| homología | `homolog_present` (fuerte), `homolog_weak`, `homolog_low` (hay hit pero no pasa umbrales), `no_homolog`, `not_checked` |",
        "| dominio | `domain_present`, `no_domain`, `not_checked` |",
        "",
        "El valor final es `homolog+dominio` (p. ej. `homolog_present+domain_present`).",
        "`no_domain` solo se emite si TODOS los dominios esperados estaban en el panel",
        "y ninguno apareció; si falta un HMM o no hay dominios Pfam, es `not_checked`.",
        "",
        "## Resultados",
        "",
        "| estado | n |",
        "|---|---:|",
    ]
    for k, v in sorted(c["status"].items(), key=lambda kv: -kv[1]):
        lines.append(f"| `{k}` | {v} |")
    lines += [
        "",
        f"- Con homología fuerte: **{c['homolog_present']}** · débil: **{c['homolog_weak']}**"
        f" · baja (hit bajo umbral): **{c['homolog_low']}**.",
        f"- Con dominio detectado: **{c['domain_axis'].get('domain_present', 0)}**.",
        f"- Ambos (homología fuerte + dominio): **{c['homolog_and_domain']}**.",
        f"- Sin homólogo ni dominio: **{c['neither_homolog_nor_domain']}**"
        f" ({summary['gaps_still_unknown_after']} sin evidencia positiva al sumar las no evaluadas).",
        f"- `not_checked` de dominio: **{c['not_checked_domain']}** · de homología: **{c['not_checked_homolog']}**.",
        "",
        f"- Fuente de los homólogos: NCBI {summary['source_breakdown']['homolog_best_by_source'].get('ncbi', 0)}"
        f" / ab initio {summary['source_breakdown']['homolog_best_by_source'].get('abinitio', 0)}"
        f" (fuertes: {summary['source_breakdown']['homolog_strong_by_source'].get('ncbi', 0)}"
        f" / {summary['source_breakdown']['homolog_strong_by_source'].get('abinitio', 0)}).",
        f"- Dominios detectados: filas con NCBI {summary['source_breakdown']['domain_rows_by_source'].get('ncbi', 0)}"
        f", ab initio {summary['source_breakdown']['domain_rows_by_source'].get('abinitio', 0)}"
        " (una fila puede tener ambos).",
        "",
        "## 3 hallazgos",
        "",
        "1. **La homología dirigida rescata brechas que el panel previo dejaba en `unknown`.**",
        f"   De las {summary['batch_scope']['rows']} filas, "
        f"{summary['comparison_with_gap_review']['previous_homolog_status_in_batch']['not_checked']}"
        " venían como `homolog_status=not_checked`;"
        f" {summary['comparison_with_gap_review']['rescued_from_not_checked']} ahora tienen homólogo"
        f" fuerte o débil ({c['homolog_present']} fuertes). No es ausencia de enzima: es reacción no mapeada.",
        "2. **Los dominios añaden señal donde la identidad no llega.**",
        f"   {c['domain_axis'].get('domain_present', 0)} filas tienen un dominio Pfam específico del productor,"
        f" y {c['domain_axis'].get('domain_present', 0) - c['homolog_and_domain'] - c['status'].get('homolog_weak+domain_present', 0)}"
        " de ellas no tienen ni siquiera homología débil; son candidatas a"
        " homología divergente o anotación incompleta,",
        "   no a pérdida de la enzima.",
        "3. **El bloque sin evidencia positiva es el duro restante.**",
        f"   {c['neither_homolog_nor_domain']} filas no muestran ni homólogo ni dominio;"
        f" {c['not_checked_homolog']} no son evaluables (sin productor no eucariota);"
        f" en total, {summary['gaps_still_unknown_after']} de las {summary['batch_scope']['rows']}"
        f" ({round(100 * summary['gaps_still_unknown_after'] / summary['batch_scope']['rows'])}%)"
        " quedan sin evidencia positiva.",
        "   Ahí sí cabe hablar de posible ausencia/auxotrofía, pero el lote no la confirma.",
        "",
        "## Límites",
        "",
    ]
    lines += [f"- {x}" for x in summary["limits"]]
    lines += [
        "",
        "## Reproducir",
        "",
        "```bash",
        "cd <raíz-del-repositorio>",
        f"bash scripts/evidence/{'run_all.sh' if BATCH == '01' else f'run_batch_{BATCH}.sh'}",
        "python3 scripts/evidence/08_validate.py",
        "```",
        "",
        "Pasos: `01_select_worklist.py` → `02_build_producer_db.py` →",
        "`03_resolve_domains.py` → `04_build_domain_panel.py` →",
        "`05_run_homology.py` → `06_run_domains.py` → `07_build_evidence.py`.",
        "Los intermedios (selección, caché Pfam, crudos DIAMOND/domtbl, panel) quedan",
        f"en `{WORK.name}/`; `03` y `04` reanudan desde caché y `05`/`06` desde los crudos.",
        "",
        "## Archivos",
        "",
        "| archivo | contenido |",
        "|---|---|",
        f"| `{TABLE_OUT.name}` | {b['rows']} filas × 12 columnas (contrato) |",
        f"| `{SUMMARY_OUT.name}` | lote, métodos, conteos, límites y hash de la tabla |",
        f"| `{README_OUT.name}` | este documento |",
        *([f"| `no_bacterial_reference_batch_{BATCH}.tsv` | "
           f"{summary['no_bacterial_reference']} filas sin productor "
           "bacteriano curado (no ausencia) |"] if BATCH == "04" else []),
        f"| `{WORK.name}/` | intermedios reproducibles (no integrables) |",
        "",
    ]
    if u:
        idx = lines.index("## Lote")
        lines[idx:idx] = [
            "## Universo",
            "",
            f"- `enzymatic_missing` vigente: **{u['enzymatic_missing']}** filas.",
            f"- Cubiertas por lotes previos: **{u['covered_total']}** "
            f"({u['covered_in_universe']} siguen en el universo; "
            f"{u['promoted_total']} ya promovidas fuera).",
            f"- Accionables: **{u['actionable_total']}** "
            f"({u['actionable_in_universe']} en el universo; "
            f"{u['reannotated_in_actionable']} reannotated).",
            f"- Pool restante: **{u['pool']}** claves. La resta simple no "
            "coincide porque las promovidas ya estaban entre las cubiertas y "
            "las reannotated entre las accionables (las listas se solapan).",
            "",
        ]
    README_OUT.write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
