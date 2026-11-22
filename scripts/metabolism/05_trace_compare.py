#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_trace_compare.py — Compara los reportes de 04_network_trace.py entre genomas. (v1.0)

Una tarea: leer los {genome}_network_trace.json que ya escribió 04 (NO
re-analiza redes ni Rhea) y separar las brechas por frecuencia entre genomas:

  🎯 SISTEMÁTICAS — presentes en >= --min-freq de los genomas (por defecto, en
     TODOS). Casi siempre artefactos estructurales (clases genéricas de Rhea,
     cofactores tipo NTP/acil-CoA) o biología central que el anotador no
     cubre: aparecen igual en cualquier genoma, así que no discriminan.
  🤝 COMPARTIDAS — en más de un genoma pero no en todos.
  🔬 ESPECÍFICAS — solo en un genoma: candidatas reales de enzima faltante
     (o de anotación perdida en ese genoma).

Además agrega la tría de 04 >= v3.3 (campo 'triage' de gaps_missing_enzyme);
con reportes antiguos la tría sale como 'unknown'.

Salida:
  - consola: resumen por genoma + top de brechas sistemáticas/específicas
  - --out (por defecto {reports-dir}/trace_compare.json): JSON completo con la
    matriz de frecuencias, la tría agregada y el resumen por genoma.

Uso:
    python scripts/metabolism/05_trace_compare.py
    python scripts/metabolism/05_trace_compare.py \
        --reports-dir data/metabolism --out data/metabolism/trace_compare.json
"""
import argparse
import json
import logging
import math
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

SCRIPT_VERSION = "1.0"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-7s | %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("trace_compare")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_REPORTS_DIR = PROJECT_ROOT / "data" / "metabolism"

# Misma tría que 04 v3.3 (copia, para no importar el otro script)
TRIAGE_LABELS = {
    "closable_with_genome_enzyme": "closable (enzima ya en la red)",
    "generic_class": "clase genérica (artefacto probable)",
    "generic_producer": "productor genérico (is-a de ChEBI pendiente)",
    "genuine_missing_enzyme": "genuina (sin productor específico en Rhea)",
}


def load_reports(reports_dir):
    """Lee *_network_trace.json y devuelve {genome_id: reporte}."""
    reports = {}
    for p in sorted(reports_dir.glob("*_network_trace.json")):
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            log.warning(f"{p.name}: no se pudo leer ({exc}); se omite")
            continue
        if not isinstance(data, dict) or "gaps_missing_enzyme" not in data:
            log.warning(f"{p.name}: no parece un reporte de 04; se omite")
            continue
        gid = data.get("genome") or p.name[:-len("_network_trace.json")]
        reports[gid] = data
    return reports


def gap_frequency(reports):
    """chebi -> {name, triage (mayoritario), genomes (ordenada)}.

    La tría se agrega por voto mayoritario entre los reportes que mencionan
    la brecha (determinista en empates: orden alfabético de la cubeta).
    """
    freq = {}
    for gid, rep in reports.items():
        for e in rep.get("gaps_missing_enzyme", []):
            c = e.get("chebi")
            if c is None:
                continue
            d = freq.setdefault(c, {"name": e.get("name", c),
                                    "votes": defaultdict(int),
                                    "genomes": set()})
            d["genomes"].add(gid)
            d["votes"][e.get("triage", "unknown")] += 1
    out = {}
    for c, d in freq.items():
        tri = max(sorted(d["votes"].items()), key=lambda kv: kv[1])[0]
        out[c] = {"name": d["name"], "triage": tri,
                  "genomes": sorted(d["genomes"])}
    return out


def per_genome_summary(reports):
    """Resumen de counts por genoma, robusto a reportes viejos."""
    out = {}
    for gid, rep in reports.items():
        st = rep.get("stats", {})
        out[gid] = {
            "reactions_matched": st.get("reactions_matched"),
            "metabolites": st.get("metabolites"),
            "gaps_total": len(rep.get("gaps_missing_enzyme", [])),
            "gaps_by_class": rep.get("gaps_by_class", {}),
            "sources": len(rep.get("sources_importable", [])),
            "sinks": len(rep.get("sinks_final_products", [])),
            "downstream_gaps": len(rep.get("downstream_gaps", [])),
        }
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Compara reportes de 04_network_trace.py entre genomas "
                    "(brechas sistemáticas vs específicas).")
    ap.add_argument("--reports-dir", type=Path, default=DEFAULT_REPORTS_DIR,
                    help="Directorio con los {genome}_network_trace.json "
                         "(por defecto data/metabolism)")
    ap.add_argument("--out", type=Path, default=None,
                    help="JSON de salida (por defecto "
                         "{reports-dir}/trace_compare.json)")
    ap.add_argument("--min-freq", type=float, default=1.0,
                    help="Fracción de genomas para considerar una brecha "
                         "sistemática (1.0 = en todos; 0.9 = en el 90%%)")
    ap.add_argument("--top", type=int, default=15,
                    help="Cuántas brechas mostrar por sección")
    args = ap.parse_args()

    reports = load_reports(args.reports_dir)
    if not reports:
        log.error(f"No hay reportes *_network_trace.json en {args.reports_dir}")
        return 1
    n = len(reports)
    threshold = max(2, math.ceil(args.min_freq * n))
    if n < 2:
        log.warning("solo hay 1 reporte: sin comparación real, todas las "
                    "brechas quedan como 'específicas'")

    freq = gap_frequency(reports)
    systematic, shared, specific = [], [], []
    for c, d in freq.items():
        k = len(d["genomes"])
        if k >= threshold:
            systematic.append(c)
        elif k == 1:
            specific.append(c)
        else:
            shared.append(c)
    systematic.sort(key=lambda c: (-len(freq[c]["genomes"]), freq[c]["name"]))
    shared.sort(key=lambda c: (-len(freq[c]["genomes"]), freq[c]["name"]))
    specific.sort(key=lambda c: freq[c]["name"].lower())

    log.info(f"== Comparación de trazas: {n} genomas "
             f"(umbral sistemática = {threshold}/{n}) ==")
    for gid in sorted(reports):
        rep = reports[gid]
        st = rep.get("stats", {})
        gc = rep.get("gaps_by_class") or {}
        tri = "; ".join(
            f"{TRIAGE_LABELS.get(k, k).split(' (')[0]} {v}"
            for k, v in sorted(gc.items()) if v) or "n/d"
        log.info(f"  {gid}: {st.get('reactions_matched', '?')} reacciones | "
                 f"brechas {len(rep.get('gaps_missing_enzyme', []))} ({tri}) | "
                 f"sinks {len(rep.get('sinks_final_products', []))} | "
                 f"fuentes {len(rep.get('sources_importable', []))}")

    log.info(f"\n🎯 Brechas SISTEMÁTICAS (en >= {threshold}/{n} genomas — "
             f"artefacto estructural probable: no discrimina entre genomas): "
             f"{len(systematic)}")
    for c in systematic[:args.top]:
        d = freq[c]
        lab = TRIAGE_LABELS.get(d["triage"], d["triage"])
        log.info(f"  {d['name']} ({c}) [{lab}] — {len(d['genomes'])}/{n} genomas")

    log.info(f"\n🔬 Brechas ESPECÍFICAS de un genoma (candidatas reales de "
             f"enzima faltante o anotación perdida): {len(specific)}")
    for c in specific[:args.top]:
        d = freq[c]
        lab = TRIAGE_LABELS.get(d["triage"], d["triage"])
        log.info(f"  {d['name']} ({c}) [{lab}] — solo en {d['genomes'][0]}")

    log.info(f"\n🤝 Brechas COMPARTIDAS por algunos genomas "
             f"(2..{threshold - 1}): {len(shared)}")
    for c in shared[:args.top]:
        d = freq[c]
        lab = TRIAGE_LABELS.get(d["triage"], d["triage"])
        gs = ", ".join(d["genomes"][:6]) + ("..." if len(d["genomes"]) > 6
                                            else "")
        log.info(f"  {d['name']} ({c}) [{lab}] — {len(d['genomes'])}/{n}: {gs}")

    out_path = args.out or (args.reports_dir / "trace_compare.json")
    payload = {
        "script": "05_trace_compare",
        "script_version": SCRIPT_VERSION,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "reports_dir": str(args.reports_dir),
        "genomes": sorted(reports),
        "systematic_threshold": threshold,
        "min_freq": args.min_freq,
        "per_genome": per_genome_summary(reports),
        "gap_frequency": [
            {"chebi": c, "name": freq[c]["name"], "triage": freq[c]["triage"],
             "freq": len(freq[c]["genomes"]),
             "genomes": freq[c]["genomes"]}
            for c in sorted(freq, key=lambda x: (-len(freq[x]["genomes"]),
                                                 freq[x]["name"]))],
        "systematic_gaps": systematic,
        "shared_gaps": shared,
        "specific_gaps": specific,
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    log.info(f"\nReporte escrito: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
