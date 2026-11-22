#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_trace_effect.py — Gestionar y medir una RE-CORRIDA de 04 (antes/después).

El flujo completo de mejora es:

    1. Congelar el estado actual de las trazas  (snapshot 'base')
    2. Re-correr 04 en batch con nuevos flags   (p. ej. --resolve-generics)
    3. Congelar el estado nuevo                 (snapshot '<tag>')
    4. Comparar: qué brechas se resolvieron, cuáles persisten, cuáles son nuevas

04 escribe {genome}_network_trace.json EN el mismo fichero (lo sobreescribe),
así que SIN snapshot del 'antes' la comparación es imposible. Este script
automatiza los 4 pasos y produce un delta por genoma y por brecha.

Subcomandos
-----------
  snapshot --tag base
      Copia los *_network_trace.json actuales de data/metabolism/ a
      data/metabolism/trace_runs/base/ (no sobreescribe un tag existente).

  run --tag rg [--04-args "--resolve-generics"] [--with-05 --min-freq 0.85]
      Ejecuta los 4 pasos: snapshot 'base' (solo si no existe), lanza 04 en
      batch con los flags dados (hereda la consola), snapshot '<tag>' y
      compare. Con --with-05 termina invocando 05 sobre el estado nuevo.

  compare --a base --b rg [--top 20]
      Delta entre dos snapshots: por genoma (brechas/sinks/fuentes/cubetas) y
      por brecha (resueltas / nuevas / cambio de cubeta). Escribe
      data/metabolism/trace_effect_<a>_vs_<b>.json.

Uso típico (todo en uno):
    python scripts/metabolism/06_trace_effect.py run --tag rg \
        --04-args "--resolve-generics" --with-05 --min-freq 0.85

O paso a paso:
    python scripts/metabolism/06_trace_effect.py snapshot --tag base
    python scripts/metabolism/04_network_trace.py --resolve-generics
    python scripts/metabolism/06_trace_effect.py snapshot --tag rg
    python scripts/metabolism/06_trace_effect.py compare --a base --b rg
"""
import argparse
import datetime
import json
import logging
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_VERSION = "1.0"
log = logging.getLogger("trace_effect")
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-7s | %(message)s",
                    datefmt="%H:%M:%S")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT_04 = Path(__file__).resolve().parent / "04_network_trace.py"
SCRIPT_05 = Path(__file__).resolve().parent / "05_trace_compare.py"

# Mismas claves/orden que TRIAGE_LABELS de 04 v3.3 (campo 'triage' del reporte)
TRIAGE_KEYS = [
    "closable_with_genome_enzyme",
    "generic_class",
    "generic_producer",
    "genuine_missing_enzyme",
]
TRIAGE_SHORT = {
    "closable_with_genome_enzyme": "closable",
    "generic_class": "clase genérica",
    "generic_producer": "productor genérico",
    "genuine_missing_enzyme": "genuina",
}


def net_dir_default():
    return PROJECT_ROOT / "data" / "metabolism"


def runs_dir(net_dir: Path) -> Path:
    return net_dir / "trace_runs"


def snapshot(net_dir: Path, tag: str, overwrite: bool = False) -> int:
    """Copia los trazas actuales ({genome}_network_trace.json de primer nivel)
    a trace_runs/<tag>/. Devuelve el número de ficheros copiados."""
    dst = runs_dir(net_dir) / tag
    if dst.exists() and not overwrite:
        log.error(f"El snapshot '{tag}' ya existe: {dst}")
        log.error("Elige otro --tag o usa --overwrite si de verdad quieres "
                  "reemplazarlo.")
        return -1
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for p in sorted(net_dir.glob("*_network_trace.json")):
        shutil.copy2(p, dst / p.name)
        n += 1
    log.info(f"Snapshot '{tag}': {n} trazas -> {dst}")
    return n


def load_dir(d: Path):
    """Lee todos los reportes de traza de un directorio -> {genome: report}."""
    out = {}
    for p in sorted(d.glob("*_network_trace.json")):
        try:
            rep = json.loads(p.read_text(encoding="utf-8"))
            out[p.name[:-len("_network_trace.json")]] = rep
        except (json.JSONDecodeError, OSError) as e:
            log.warning(f"No se pudo leer {p.name}: {e}")
    return out


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def compare(a_dir: Path, b_dir: Path, top: int, out_path: Path):
    """Delta entre dos snapshots. Escribe JSON y resume en consola."""
    A, B = load_dir(a_dir), load_dir(b_dir)
    if not A or not B:
        log.error(f"Snapshots vacíos: {a_dir} ({len(A)}), {b_dir} ({len(B)})")
        return None
    common = sorted(set(A) & set(B))
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))
    log.info(f"Comparando {a_dir.name} vs {b_dir.name}: {len(common)} genomas "
             f"comunes (solo en A: {len(only_a)}, solo en B: {len(only_b)})")

    # ---- por genoma ----
    per_genome = []
    tri_a = {k: 0 for k in TRIAGE_KEYS}
    tri_b = {k: 0 for k in TRIAGE_KEYS}
    for g in common:
        ga = A[g].get("gaps_missing_enzyme", [])
        gb = B[g].get("gaps_missing_enzyme", [])
        for r in ga:
            tri_a[r.get("triage", "genuine_missing_enzyme")] += 1
        for r in gb:
            tri_b[r.get("triage", "genuine_missing_enzyme")] += 1
        per_genome.append({
            "genome": g,
            "gaps_before": len(ga), "gaps_after": len(gb),
            "delta": len(gb) - len(ga),
            "sinks_before": len(A[g].get("sinks_final_products", [])),
            "sinks_after": len(B[g].get("sinks_final_products", [])),
            "sources_before": len(A[g].get("sources_importable", [])),
            "sources_after": len(B[g].get("sources_importable", [])),
        })

    n = len(common)
    gaps_a = [p["gaps_before"] for p in per_genome]
    gaps_b = [p["gaps_after"] for p in per_genome]
    resolved_n = sum(1 for p in per_genome if p["delta"] < 0)
    worsened_n = sum(1 for p in per_genome if p["delta"] > 0)

    # ---- por brecha (identidad = chebi) ----
    # Índice chebi -> nombre, construido una vez por lado (evita O(genomas)
    # por brecha al etiquetar las listas resueltas/nuevas).
    names_a, names_b = {}, {}

    def gap_map(rep):
        return {r["chebi"]: r for r in rep.get("gaps_missing_enzyme", [])
                if "chebi" in r}

    for g in common:
        for c, r in gap_map(A[g]).items():
            names_a.setdefault(c, r.get("name", c))
        for c, r in gap_map(B[g]).items():
            names_b.setdefault(c, r.get("name", c))

    freq_a, freq_b = {}, {}          # chebi -> en cuántos genomas del set
    tri_changes = {}                 # chebi -> {from,to,genomes:[...]}
    resolved, new = {}, {}
    for g in common:
        ma, mb = gap_map(A[g]), gap_map(B[g])
        for c, r in ma.items():
            freq_a[c] = freq_a.get(c, 0) + 1
        for c, r in mb.items():
            freq_b[c] = freq_b.get(c, 0) + 1
        for c in set(ma) - set(mb):
            resolved[c] = resolved.get(c, 0) + 1
        for c in set(mb) - set(ma):
            new[c] = new.get(c, 0) + 1
        for c in set(ma) & set(mb):
            ta, tb = ma[c].get("triage"), mb[c].get("triage")
            if ta != tb:
                e = tri_changes.setdefault(c, {"from": ta, "to": tb,
                                               "genomes": []})
                e["genomes"].append(g)

    def entry(chebi, freq):
        return {"chebi": chebi,
                "name": names_a.get(chebi) or names_b.get(chebi) or chebi,
                "freq": freq}

    resolved_l = [entry(c, f) for c, f in
                  sorted(resolved.items(), key=lambda kv: -kv[1])]
    new_l = [entry(c, f) for c, f in
             sorted(new.items(), key=lambda kv: -kv[1])]
    tri_changes_l = [{"chebi": c,
                      "name": names_a.get(c) or names_b.get(c) or c,
                      "from": e["from"], "to": e["to"],
                      "genomes": len(e["genomes"])}
                     for c, e in sorted(tri_changes.items(),
                                        key=lambda kv: -len(kv[1]["genomes"]))]
    persistent = len(set(freq_a) & set(freq_b))

    summary = {
        "script_version": SCRIPT_VERSION,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "before": str(a_dir), "after": str(b_dir),
        "genomes_compared": n,
        "genomes_only_in_before": only_a, "genomes_only_in_after": only_b,
        "cohort": {
            "gaps_total_before": sum(gaps_a), "gaps_total_after": sum(gaps_b),
            "gaps_mean_before": round(sum(gaps_a) / n, 1) if n else 0,
            "gaps_mean_after": round(sum(gaps_b) / n, 1) if n else 0,
            "gaps_median_before": _median(gaps_a) if n else 0,
            "gaps_median_after": _median(gaps_b) if n else 0,
            "genomes_improved": resolved_n, "genomes_worsened": worsened_n,
            "triage_before": tri_a, "triage_after": tri_b,
        },
        "gaps": {
            "resolved_count": len(resolved_l),
            "new_count": len(new_l),
            "persistent_count": persistent,
            "resolved_top": resolved_l[:100],
            "new_top": new_l[:100],
            "triage_changes": tri_changes_l[:100],
        },
        "per_genome": sorted(per_genome, key=lambda p: p["delta"]),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                        encoding="utf-8")

    # ---- consola ----
    c = summary["cohort"]
    log.info(f"\n== EFECTO {a_dir.name} -> {b_dir.name} ({n} genomas) ==")
    log.info(f"Brechas totales: {c['gaps_total_before']} -> "
             f"{c['gaps_total_after']}  "
             f"({c['gaps_total_after'] - c['gaps_total_before']:+d})")
    log.info(f"Media/genoma: {c['gaps_mean_before']} -> "
             f"{c['gaps_mean_after']}   mediana: {c['gaps_median_before']} -> "
             f"{c['gaps_median_after']}")
    log.info(f"Genomas que mejoran: {c['genomes_improved']}/{n}  "
             f"empeoran: {c['genomes_worsened']}/{n}")
    log.info("Cubetas de tría (antes -> después):")
    for k in TRIAGE_KEYS:
        log.info(f"  {TRIAGE_SHORT[k]:20s} {tri_a[k]:6d} -> {tri_b[k]:6d}  "
                 f"({tri_b[k] - tri_a[k]:+d})")
    log.info(f"\n🎯 Brechas RESUELTAS (desaparecieron): {len(resolved_l)}; "
             f"nuevas: {len(new_l)}; persistentes: {persistent}")
    for r in resolved_l[:top]:
        log.info(f"  ✔ {r['name']} ({r['chebi']}) — resuelta en "
                 f"{r['freq']}/{n} genomas")
    if new_l:
        log.info(f"  ⚠ NUEVAS brechas (top):")
        for r in new_l[:top]:
            log.info(f"  ✚ {r['name']} ({r['chebi']}) — aparece en "
                     f"{r['freq']}/{n} genomas")
    if tri_changes_l:
        log.info(f"  ↻ Cambio de cubeta (top):")
        for r in tri_changes_l[:top]:
            log.info(f"  ~ {r['name']} ({r['chebi']}): "
                     f"{TRIAGE_SHORT.get(r['from'], r['from'])} -> "
                     f"{TRIAGE_SHORT.get(r['to'], r['to'])} "
                     f"en {r['genomes']} genomas")
    log.info(f"\nPer-genoma (mayores mejoras):")
    for p in sorted(per_genome, key=lambda x: x["delta"])[:top]:
        log.info(f"  {p['genome']}: {p['gaps_before']} -> {p['gaps_after']} "
                 f"({p['delta']:+d})")
    log.info(f"Reporte escrito: {out_path}")
    return summary


def run_full(net_dir: Path, tag: str, args_04: str, with_05: bool,
             min_freq: float):
    """Snapshot 'base' -> 04 batch -> snapshot '<tag>' -> compare."""
    base = runs_dir(net_dir) / "base"
    if base.exists():
        log.info(f"Snapshot 'base' ya existe, se reutiliza: {base}")
    else:
        n = snapshot(net_dir, "base")
        if n <= 0:
            log.error("Sin snapshot 'base': no hay nada que comparar. ¿Hay "
                      f"*_network_trace.json en {net_dir}?")
            return
    cmd = [sys.executable, str(SCRIPT_04)] + shlex.split(args_04)
    log.info(f"Ejecutando 04: {' '.join(cmd)}")
    rc = subprocess.run(cmd, cwd=str(PROJECT_ROOT)).returncode
    if rc != 0:
        log.error(f"04 terminó con código {rc}: no se compara ni se hace "
                  "snapshot del estado nuevo.")
        sys.exit(rc)
    if snapshot(net_dir, tag) < 0:
        return
    compare(runs_dir(net_dir) / "base", runs_dir(net_dir) / tag,
            top=15, out_path=net_dir / f"trace_effect_base_vs_{tag}.json")
    if with_05:
        cmd5 = [sys.executable, str(SCRIPT_05), "--min-freq", str(min_freq)]
        log.info(f"Ejecutando 05: {' '.join(cmd5)}")
        subprocess.run(cmd5, cwd=str(PROJECT_ROOT)).returncode


def main():
    ap = argparse.ArgumentParser(
        description="Gestión y medición de re-corridas de 04 (antes/después).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--net-dir", type=Path, default=None,
                       help="Directorio de redes/trazas "
                            "(por defecto data/metabolism)")

    p_snap = sub.add_parser("snapshot",
                            help="Congelar trazas actuales bajo un tag")
    common(p_snap)
    p_snap.add_argument("--tag", required=True)
    p_snap.add_argument("--overwrite", action="store_true")

    p_run = sub.add_parser(
        "run", help="snapshot base -> 04 batch -> snapshot tag -> compare")
    common(p_run)
    p_run.add_argument("--tag", required=True,
                       help="Tag para el estado posterior (p. ej. 'rg')")
    p_run.add_argument("--04-args", dest="args_04", default="--resolve-generics",
                       help="Flags extra para 04 (por defecto "
                            "'--resolve-generics'; añade p. ej. "
                            "'--rhea-dir X --net-dir Y' si no usas los "
                            "defaults)")
    p_run.add_argument("--with-05", action="store_true",
                       help="Al terminar, correr 05 sobre el estado nuevo")
    p_run.add_argument("--min-freq", type=float, default=0.85,
                       help="--min-freq para 05 (por defecto 0.85)")

    p_cmp = sub.add_parser("compare", help="Delta entre dos snapshots")
    common(p_cmp)
    p_cmp.add_argument("--a", required=True, help="Tag 'antes'")
    p_cmp.add_argument("--b", required=True, help="Tag 'después'")
    p_cmp.add_argument("--top", type=int, default=15,
                       help="Cuántos ejemplos listar en consola")
    p_cmp.add_argument("--out", type=Path, default=None,
                       help="JSON de salida (por defecto "
                            "trace_effect_<a>_vs_<b>.json en net-dir)")

    args = ap.parse_args()
    net_dir = args.net_dir or net_dir_default()

    if args.cmd == "snapshot":
        snapshot(net_dir, args.tag, overwrite=args.overwrite)
    elif args.cmd == "run":
        run_full(net_dir, args.tag, args.args_04, args.with_05,
                 args.min_freq)
    elif args.cmd == "compare":
        a, b = runs_dir(net_dir) / args.a, runs_dir(net_dir) / args.b
        out = args.out or (net_dir /
                           f"trace_effect_{args.a}_vs_{args.b}.json")
        compare(a, b, args.top, out)


if __name__ == "__main__":
    main()
