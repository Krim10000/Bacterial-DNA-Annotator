#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_network_view.py — Portal interactivo de redes metabólicas (todos los genomas).

UN único HTML (portal) con:

  * Menú desplegable de genomas — no se pide --genome: el usuario elige en
    la interfaz y el grafo se recarga al vuelo.
  * Vista RUTAS (por defecto): diagrama de flujo por brechas — un bloque
    compacto por metabolito-brecha (enzima ausente → brecha → consumidores);
    pocos elementos, inspección directa. La red completa de cientos de nodos
    queda en la pestaña RED, perezosa (solo se construye al abrirla).
  * Interruptor ANTES / DESPUÉS: "antes" muestra la red sin los resultados
    del 08; "después" resalta en la red los candidatos a enzima faltante
    validados por plegado (halo amarillo pulsante) y añade al panel los
    veredictos (pLDDT, dominio, best hit).
  * Botón "Ejecutar 08": lanza 08_fold_candidates.py como subproceso desde
    la propia interfaz, con opciones (ORFs, herramienta de similitud,
    presupuesto de plegado, proveedor ESM y ESM_API_KEY/ESM_API_URL
    opcionales) y log en vivo; al terminar recarga los candidatos y pasa a
    "después" automáticamente.

Dos modos de uso:

  1) Servidor (POR DEFECTO):  python 07_network_view.py
       -> arranca el backend en http://127.0.0.1:8765 y ABRE EL NAVEGADOR
          (evítalo con --no-browser). Backend /api/*: payloads bajo demanda,
          candidatos, run08, status. Las rutas (--gene-map, --proteins-dir,
          --genomes-dir-fna, --genome-meta) que pases al servidor se
          reenvían a 08. (--serve se acepta por compatibilidad.)
  2) Estático (offline):   python 07_network_view.py --static
       -> escribe network_portal.html con todos los genomas embebidos; el
          botón de 08 muestra la pista de cómo lanzar el servidor.

Colores clave: naranja = reacción codificada solo por proteína(s)
hipotética(s); rojo intenso discontinuo = nodo FANTASMA (enzima ausente
de Rhea que produciría la brecha); halos amarillos = candidatos 08. El
tamaño del ADN es explícito en la cabecera (bp/Mb + aviso GENOMA REDUCIDO
<1 Mb), en cada panel (gen en bp y % del genoma) y en la tira de contexto
del contig (±6 kb, ancho ∝ bp).

Entradas: {genome}_network.json + {genome}_network_trace.json (04), cache
Rhea, {genome}_gene_map.tsv|json (opcional), --genome-meta, --cohort (05),
--candidates dir de 08. Salida: servidor en vivo o network_portal.html en
data/metabolism/network_views/ (--static).
"""
import argparse
import glob
import gzip
import json
import os
import subprocess
import sys
import threading
import time
import urllib.parse
import webbrowser
from collections import defaultdict
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
m04 = importlib.import_module("04_network_trace")

SCRIPT_VERSION = "2.2"
DEFAULT_PORT = 8765

COLORS = {
    "chem_default": "#94a3b8",
    "chem_source": "#0284c7",
    "chem_gap": "#dc2626",
    "chem_sink": "#16a34a",
    "rxn_ok": "#0f766e",
    "rxn_hypothetical": "#f97316",
    "rxn_ghost": "#be123c",
    "candidate": "#facc15",
}
TRIAGE_SHORT = {
    "closable_with_genome_enzyme": "closable",
    "generic_class": "clase genérica",
    "generic_producer": "productor genérico",
    "genuine_missing_enzyme": "genuina",
}


# ----------------------------------------------------------------- entradas
def load_gene_map(path: Path):
    """{rhea_int: [filas]} desde {genome}_gene_map.tsv|json (también la usa 08)."""
    rows = []
    if path.suffix == ".json":
        rows = json.loads(path.read_text(encoding="utf-8"))
    else:
        lines = path.read_text(encoding="utf-8").splitlines()
        header = None
        for line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip().lower() for h in parts]
                continue
            rows.append(dict(zip(header, parts)))
    by_rhea = {}
    for r in rows:
        norm = m04.normalize_rhea_id(r.get("rhea_id", ""))
        if not norm:
            continue
        row = {
            "gene_id": r.get("gene_id", ""),
            "product": r.get("product", ""),
            "locus_tag": r.get("locus_tag", ""),
            "contig": r.get("contig", ""),
            "start": int(r.get("start") or 0),
            "end": int(r.get("end") or 0),
            "strand": r.get("strand", ""),
            "length_bp": int(r.get("gene_length_bp") or 0) or
                         (abs(int(r.get("end") or 0) - int(r.get("start") or 0)) + 1),
            "hypothetical": str(r.get("hypothetical", "")).strip().lower()
                            in ("true", "1", "yes", "y", "sí", "si", "hipothetical"),
            "sequence": (r.get("sequence") or "").strip(),  # 08 la usa; 07 no
        }
        by_rhea.setdefault(norm[1], []).append(row)
    return by_rhea


def load_genome_meta(path: Path, genome: str):
    """(total_length_bp, gc_percent, n_genes, n_hypothetical) o None.

    Acepta JSON {genome: {...}} o lista de dicts con 'genome_id', o TSV con
    cabecera: genome_id  total_length_bp  [gc_percent] [n_genes]
    [n_hypothetical]."""
    if not path or not path.exists():
        return None
    m = None
    if path.suffix == ".json":
        all_m = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(all_m, dict):
            m = all_m.get(genome)
        elif isinstance(all_m, list):
            m = next((x for x in all_m if x.get("genome_id") == genome), None)
    else:
        hdr, m = None, None
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [h.strip().lower() for h in parts]
                if "genome_id" not in hdr:  # sin cabecera: posición fija
                    hdr = ["genome_id", "total_length_bp", "gc_percent",
                           "n_genes", "n_hypothetical"]
                continue
            if parts and parts[0].strip() == genome:
                m = dict(zip(hdr, parts))
                break
    if not m:
        return None
    return {
        "total_length_bp": int(m.get("total_length_bp") or 0) or None,
        "gc_percent": float(m["gc_percent"]) if m.get("gc_percent") else None,
        "n_genes": int(m["n_genes"]) if m.get("n_genes") else None,
        "n_hypothetical": int(m["n_hypothetical"]) if m.get("n_hypothetical") else None,
    }


def load_cohort(path: Path):
    """Set de CHEBI sistemáticos desde el trace_compare.json del 05
    ('systematic_gaps' es una lista de chebi; se aceptan dicts también)."""
    if not path or not path.exists():
        return set()
    data = json.loads(path.read_text(encoding="utf-8"))
    out = set()
    for e in data.get("systematic_gaps", []):
        out.add(e.get("chebi") if isinstance(e, dict) else e)
    return {c for c in out if c}


def load_candidates_list(cand_dir: Path, genome: str):
    """Lista completa de candidatos desde el output de 08 (o [])."""
    p = cand_dir / f"{genome}_candidates.json" if cand_dir else None
    if not p or not p.exists():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return [c for c in data.get("candidates", []) if c.get("id")]
    except (json.JSONDecodeError, OSError):
        return []


# ------------------------------------------------------------- grafo payload
BRANCH_HUB_DEGREE = 10


def compute_branches(nodes, edges, gaps, hub_degree=BRANCH_HUB_DEGREE):
    """Ramas del metabolismo para la vista RAMAS (cajas por ruta).

    Cada rama = componente conexo de la red tras CORTAR en los metabolitos
    'hub' (grado >= hub_degree: ATP, agua y compañía ya van filtrados como
    moneda, pero nucleótidos-azúcar y otros hubs siguen cosiendo ramas
    distintas). Una rama queda INCOMPLETA si alguna de sus reacciones toca
    un metabolito en brecha del 04 (dead end: se consume y nada lo produce).
    Devuelve lista de dicts {id, name, complete, reactions, gaps, n_genes,
    n_rxn} ordenada por accionabilidad (incompletas primero, tamaño)."""
    chem = {n["id"]: n for n in nodes if n["kind"] == "chem"}
    rxn = {n["id"]: n for n in nodes if n["kind"] == "rxn"}
    # la moneda (ATP, agua...) cose TODAS las ramas: se corta siempre para
    # esta vista, independientemente de --keep-currency
    cut = {cid for cid, n in chem.items() if n.get("currency")}
    deg = defaultdict(int)
    for e in edges:
        if e[0] in chem:
            deg[e[0]] += 1
        if e[1] in chem:
            deg[e[1]] += 1
    cut |= {c for c, k in deg.items() if k >= hub_degree}

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    for e in edges:
        if len(e) > 2:      # arista fantasma: no conecta ramas reales
            continue
        s, d = e
        if (s in chem and s in cut) or (d in chem and d in cut):
            continue
        if s not in chem and d not in chem:
            continue
        ra, rb = find(s), find(d)
        if ra != rb:
            parent[ra] = rb

    comps = defaultdict(lambda: {"rxn": [], "chem": [], "gaps": [],
                                 "deg": defaultdict(int)})
    for nid in list(rxn) + [c for c in chem if c not in cut]:
        c = comps[find(nid)]
        c["rxn" if nid in rxn else "chem"].append(nid)

    role = {n["id"]: n.get("role", "none") for n in nodes}
    out = []
    for c in comps.values():
        if len(c["rxn"]) < 2:
            continue     # caja trivial (una reacción suelta): queda en RED
        gaps_b = sorted(ch for ch in c["chem"] if ch in gaps)
        sinks = [ch for ch in c["chem"] if role.get(ch) == "sink"]
        sources = [ch for ch in c["chem"] if role.get(ch) == "source"]
        if sinks:
            nm = " → ".join(chem[ch]["label"] for ch in sinks[:2])
        elif sources:
            nm = "de " + " / ".join(chem[ch]["label"] for ch in sources[:2])
        else:
            top = sorted(c["chem"], key=lambda ch: -deg[ch])[:2]
            nm = " / ".join(chem[ch]["label"] for ch in top) if top else "?"
        n_genes = sum(len(rxn[r].get("genes") or []) for r in c["rxn"])
        out.append({"id": f"b{len(out) + 1}", "name": nm,
                    "complete": not gaps_b, "reactions": sorted(c["rxn"]),
                    "gaps": gaps_b, "n_gaps": len(gaps_b),
                    "n_genes": n_genes, "n_rxn": len(c["rxn"])})
    out.sort(key=lambda b: (b["complete"], -b["n_rxn"]))
    return out


def meta_status(args):
    """'ok' | 'file_missing' | 'not_given' — para el badge ADN del portal."""
    p = getattr(args, "genome_meta", None)
    if not p:
        return "not_given"
    return "ok" if p.exists() else "file_missing"


def _fna_facts(path: Path):
    """(total_bp, gc_percent) de un FASTA nucleotídico .fna[.gz] — para el
    badge ADN cuando no hay genome_meta.tsv."""
    op = (gzip.open(path, "rt", encoding="utf-8", errors="replace")
          if path.suffix == ".gz"
          else open(path, "rt", encoding="utf-8", errors="replace"))
    tot = gc = at = 0
    with op as f:
        for line in f:
            if line.startswith(">"):
                continue
            s = line.strip().upper()
            tot += len(s)
            gc += s.count("G") + s.count("C")
            at += s.count("A") + s.count("T")
    gc_pct = round(100.0 * gc / (gc + at), 2) if (gc + at) else None
    return tot, gc_pct


def build_payload(genome, net_dir, rhea_dir, gene_map_by_rhea, genome_meta,
                  systematic_ch, keep_currency, top_refs=2,
                  hub_degree=BRANCH_HUB_DEGREE, genome_meta_status=None):
    with open(net_dir / f"{genome}_network.json", encoding="utf-8") as f:
        net = json.load(f)
    genome_rheas = set()
    for n in net.get("nodes", []):
        if n.get("type") == "reaction":
            norm = m04.normalize_rhea_id(n.get("id", ""))
            if norm:
                genome_rheas.add(norm[1])

    reactions, equations, _master = m04.load_reaction_data(rhea_dir)
    names = m04.load_chebi_names(rhea_dir)

    trace_p = net_dir / f"{genome}_network_trace.json"
    trace = json.loads(trace_p.read_text(encoding="utf-8")) if trace_p.exists() else {}
    gaps = {g["chebi"]: g for g in trace.get("gaps_missing_enzyme", [])}
    sources = {s["chebi"] for s in trace.get("sources_importable", [])}
    sinks = {s["chebi"] for s in trace.get("sinks_final_products", [])}

    total_bp = (genome_meta or {}).get("total_length_bp")

    nodes, edges = [], []
    seen_chem = set()

    def chem_node(chebi):
        if chebi in seen_chem:
            return
        seen_chem.add(chebi)
        role = "gap" if chebi in gaps else ("source" if chebi in sources
                                            else ("sink" if chebi in sinks else "none"))
        tri = gaps[chebi].get("triage") if chebi in gaps else None
        nodes.append({
            "id": chebi, "kind": "chem",
            "label": names.get(chebi, chebi),
            "role": role, "triage": TRIAGE_SHORT.get(tri, tri),
            "systematic": chebi in systematic_ch,
            "currency": (not keep_currency and chebi in m04.CURRENCY_CHEBI),
        })

    def rxn_node(rid, ghost=False, gap_chebi=None):
        rid_str = f"RHEA:{rid}"
        nid = f"ghost:{rid_str}" if ghost else rid_str
        genes = gene_map_by_rhea.get(rid, []) if not ghost else []
        hyp = bool(genes) and all(g["hypothetical"] for g in genes)
        node = {
            "id": nid, "kind": "ghost" if ghost else "rxn",
            "rhea": rid_str,
            "equation": equations.get(rid, ""),
            "label": (f"{rid_str} (ausente)" if ghost else rid_str),
            "hypothetical": hyp,
            "gap_for": gap_chebi,
            "gap_name": names.get(gap_chebi, gap_chebi) if gap_chebi else None,
            "systematic": bool(gap_chebi and gap_chebi in systematic_ch),
            "genes": [{k: g[k] for k in ("gene_id", "product", "locus_tag",
                                         "contig", "start", "end", "strand",
                                         "length_bp")} for g in genes],
        }
        if total_bp:
            for g in node["genes"]:
                g["pct_genome"] = round(100.0 * g["length_bp"] / total_bp, 3)
        nodes.append(node)
        return node

    def add_edges(rid, node_id, skip_prod=None):
        e = reactions.get(rid)
        if not e:
            return
        if e.get("transport"):
            return  # mover entre compartimentos no aporta topología visible
        for s in e["substrates"]:
            chem_node(s)
            edges.append([s, node_id])
        for p in e["products"]:
            if skip_prod and p == skip_prod:
                continue  # el enlace al gap lo aporta la arista fantasma explícita
            chem_node(p)
            edges.append([node_id, p])

    for rid in sorted(genome_rheas):
        if rid not in reactions:
            continue
        rxn_node(rid)
        add_edges(rid, f"RHEA:{rid}")

    # nodos fantasma: enzima ausente que produce cada brecha (top productores)
    for chebi, g in sorted(gaps.items()):
        chem_node(chebi)
        for r in g.get("producing_rhea", [])[:top_refs]:
            norm = m04.normalize_rhea_id(r)
            if not norm:
                continue
            rid = norm[1]
            if f"RHEA:{rid}" in {n["id"] for n in nodes}:
                add_edges(rid, f"RHEA:{rid}")
            else:
                rxn_node(rid, ghost=True, gap_chebi=chebi)
                add_edges(rid, f"ghost:RHEA:{rid}", skip_prod=chebi)
                edges.append([f"ghost:RHEA:{rid}", chebi, "ghost"])

    stats = {
        "reactions_in_genome": len(genome_rheas),
        "metabolites": len(seen_chem),
        "gaps": len(gaps), "sources": len(sources), "sinks": len(sinks),
        "hypothetical_rxns": sum(1 for n in nodes
                                 if n["kind"] == "rxn" and n["hypothetical"]),
        "ghost_rxns": sum(1 for n in nodes if n["kind"] == "ghost"),
    }
    branches = compute_branches(nodes, edges, gaps, hub_degree=hub_degree)
    stats["branches"] = len(branches)
    stats["branches_incomplete"] = sum(1 for b in branches
                                       if not b["complete"])
    meta = {
        "genome": genome, "script_version": SCRIPT_VERSION,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "dna_total_length_bp": total_bp,
        "gc_percent": (genome_meta or {}).get("gc_percent"),
        "n_genes": (genome_meta or {}).get("n_genes"),
        "n_hypothetical": (genome_meta or {}).get("n_hypothetical"),
        "reduced_genome": bool(total_bp and total_bp < 1_000_000),
        "has_gene_map": bool(gene_map_by_rhea),
        "genome_meta_status": genome_meta_status or "not_given",
    }
    return {"meta": meta, "stats": stats, "colors": COLORS,
            "branches": branches, "nodes": nodes, "edges": edges}


# ------------------------------------------------------------------- HTML
HTML_HEAD = """<!DOCTYPE html>
<html lang="es"><head><meta charset="utf-8">
<title>Redes metabólicas — portal</title>
<style>
:root { --bg:#0b1220; --panel:#111a2e; --line:#1e2a44; --txt:#dbe4f3; --dim:#7c8aa5; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--txt);
       font:13px/1.45 -apple-system, "Segoe UI", Roboto, sans-serif; }
header { padding:10px 16px; border-bottom:1px solid var(--line);
         display:flex; flex-wrap:wrap; gap:8px 18px; align-items:center; }
header h1 { font-size:15px; margin:0; color:var(--dim); font-weight:600; }
select#genomeSel { background:#16233d; color:var(--txt); border:1px solid #2c3f66;
         border-radius:5px; padding:5px 10px; font-size:14px; font-weight:600; }
.badge { background:#16233d; border:1px solid var(--line); border-radius:4px;
         padding:2px 8px; font-size:12px; color:var(--dim); }
.badge.dna { color:#7dd3fc; }
.badge.reduced { color:#fda4af; border-color:#be123c; font-weight:600; }
.badge.cand { color:#fde68a; }
.modes { display:flex; border:1px solid var(--line); border-radius:5px;
         overflow:hidden; }
.modes button { background:var(--panel); color:var(--dim); border:0;
         padding:5px 12px; font-size:13px; cursor:pointer; }
.modes button.active { background:#1d2f52; color:#fff; font-weight:600; }
.modes button.ready { color:#fde68a; font-weight:600; }
.modes button:disabled { opacity:.4; cursor:not-allowed; }
#btn-run { background:#14532d; color:#dcfce7; border:1px solid #166534;
         border-radius:5px; padding:5px 12px; font-size:13px; cursor:pointer; }
#btn-run:hover { background:#166534; }
.viewtabs { display:flex; border:1px solid var(--line); border-radius:5px;
         overflow:hidden; }
.viewtabs button { background:var(--panel); color:var(--dim); border:0;
         padding:5px 12px; font-size:13px; cursor:pointer; }
.viewtabs button.active { background:#1d2f52; color:#fff; font-weight:600; }
#main { display:flex; height:calc(100vh - 50px); }
#graph { flex:1; position:relative; overflow:hidden; cursor:grab; }
#graph:active { cursor:grabbing; }
svg#svg { width:100%; height:100%; }
#flow { position:absolute; inset:0; overflow-y:auto; padding:12px 16px 40px;
        display:none; }
.fsec { border:1px solid var(--line); border-radius:8px; background:var(--panel);
        margin:0 0 10px; overflow:hidden; }
.fsec > .fhead { display:flex; gap:8px; align-items:center; flex-wrap:wrap;
        padding:8px 12px; cursor:pointer; }
.fsec > .fhead:hover { background:#16233d; }
.fsec .fbody { padding:2px 12px 12px; border-top:1px solid var(--line); }
.frow { display:flex; align-items:center; flex-wrap:wrap; gap:6px; margin:7px 0; }
.flabel { color:#7c8aa5; font-size:11px; width:118px; flex-shrink:0;
        text-transform:uppercase; letter-spacing:.4px; }
.fbox { border:1px solid #2c3f66; border-radius:5px; padding:4px 9px;
        font-size:12px; background:#16233d; cursor:pointer; }
.fbox:hover { border-color:#5b7bb0; }
.fbox.ok { border-color:#15803d; background:#0f2a1a; }
.fbox.hyp { border-color:#c2410c; background:#2a1708; color:#fdba74; }
.fbox.ghost { border-style:dashed; border-color:#dc2626; background:#2a0a0a;
        color:#fca5a5; }
.fbox.halo { box-shadow:0 0 0 2px #fbbf24, 0 0 10px #fbbf24; }
.fbox .rh { font-weight:600; }
.fmet { background:#7f1d1d; border:1px solid #b91c1c; color:#fecaca;
        font-weight:600; padding:4px 10px; border-radius:5px; font-size:12px;
        cursor:pointer; }
.fmet.sys { box-shadow:0 0 0 2px #fbbf24; }
.farrow { color:#7c8aa5; font-size:13px; }
.fnone { color:#7c8aa5; font-style:italic; font-size:12px; }
.fnote { color:#7c8aa5; font-size:12px; }
#branches { position:absolute; inset:0; overflow-y:auto; padding:0 0 40px;
            display:none; }
.bgrid { display:grid; grid-template-columns:repeat(auto-fill, minmax(250px, 1fr));
         gap:10px; padding:10px 16px; }
.bcard { background:var(--panel); border:1px solid var(--line);
         border-radius:8px; padding:10px 12px; cursor:pointer; }
.bcard:hover { border-color:#3b82f6; }
.bcard.ok { border-left:4px solid #4ade80; }
.bcard.in { border-left:4px solid #fbbf24; }
.bname { font-weight:600; margin-bottom:6px; line-height:1.3; }
.bmeta { font-size:12px; color:var(--dim); }
.bmiss { margin-top:5px; color:#f87171; font-size:12px; }
.tag.ok { background:#14532d; color:#bbf7d0; }
.tag.in { background:#713f12; color:#fde68a; }
#panel { width:380px; border-left:1px solid var(--line); background:var(--panel);
         padding:14px 16px; overflow-y:auto; }
#panel h2 { font-size:14px; margin:0 0 6px; }
.tag { display:inline-block; padding:1px 7px; border-radius:10px;
       font-size:11px; margin:1px 3px 1px 0; }
.tag.gap { background:#7f1d1d; color:#fecaca; }
.tag.source { background:#0c4a6e; color:#bae6fd; }
.tag.sink { background:#14532d; color:#bbf7d0; }
.tag.sys { background:#111; color:#fbbf24; border:1px solid #fbbf24; }
.tag.hyp { background:#7c2d12; color:#fed7aa; }
.tag.fuerte { background:#14532d; color:#bbf7d0; }
.tag.validar { background:#713f12; color:#fde68a; }
.tag.descartado { background:#1f2937; color:#9ca3af; }
.gene { border:1px solid var(--line); border-radius:6px; margin:8px 0;
        padding:8px 10px; background:#0d1729; }
.gene .len { color:#7dd3fc; font-weight:600; }
.dna-strip { margin-top:8px; }
#legend { position:absolute; left:10px; bottom:10px; background:rgba(13,23,41,.92);
          border:1px solid var(--line); border-radius:8px; padding:8px 12px;
          font-size:12px; z-index:4; }
#legend .row { display:flex; align-items:center; gap:8px; margin:2px 0; }
.dot { width:11px; height:11px; border-radius:50%; display:inline-block; }
.sq { width:11px; height:11px; transform:rotate(45deg); display:inline-block; }
#controls { position:absolute; top:8px; left:10px; display:flex; gap:6px;
            flex-wrap:wrap; align-items:center; z-index:4; }
#controls input, #controls button { background:var(--panel); color:var(--txt);
     border:1px solid var(--line); border-radius:4px; padding:4px 8px; font-size:12px; }
#runbox { position:absolute; top:8px; right:10px; width:330px; max-height:calc(100% - 24px);
          background:rgba(13,23,41,.97); border:1px solid var(--line);
          border-radius:8px; padding:8px 12px; font-size:12px; z-index:6;
          overflow-y:auto; }
#runbox .rhead { display:flex; justify-content:space-between; align-items:center;
          gap:8px; cursor:pointer; padding:2px 0 6px; }
#runbox .rhead b { color:#fde68a; font-size:13px; }
#runbox .rhead .min { color:var(--dim); border:1px solid var(--line);
          border-radius:4px; padding:1px 7px; background:var(--panel); }
#runbox.min #runform, #runbox.min #runlog { display:none !important; }
#runbox label { display:block; margin:6px 0 2px; color:var(--dim); }
#runbox input, #runbox select { width:100%; background:#16233d; color:var(--txt);
     border:1px solid var(--line); border-radius:4px; padding:4px 6px; font-size:12px; }
#runbox .row2 { display:flex; gap:6px; margin-top:8px; }
#runlog { max-height:200px; overflow-y:auto; background:#050a14;
          border:1px solid var(--line); border-radius:4px; margin-top:8px;
          padding:6px; font:11px/1.35 ui-monospace, monospace;
          white-space:pre-wrap; color:#86efac; }
.flash { animation: pulse 1.2s ease-in-out infinite; }
@keyframes pulse { 0%,100% { stroke-opacity:.15; } 50% { stroke-opacity:1; } }
.eq { font-family: ui-monospace, monospace; font-size:12px; color:#a5b4fc; }
.dimnote { color:var(--dim); font-size:12px; }
</style></head><body>
"""

HTML_BODY = """<header>
  <h1>Red metabólica:</h1>
  <select id="genomeSel"></select>
  <span class="badge dna" id="dna-total">ADN: n/d</span>
  <span class="badge reduced" id="dna-flag" style="display:none">
    GENOMA REDUCIDO (&lt;1 Mb): ruta faltante puede ser pérdida génica</span>
  <div class="modes">
    <button id="btn-before">ANTES</button><button id="btn-after">DESPUÉS</button>
  </div>
  <div class="viewtabs">
    <button id="btn-branches" class="active" title="Ramas del metabolismo en cajas: verde = ruta completa, ámbar = incompleta (toca una brecha). Click en una caja para sus pasos">RAMAS</button>
    <button id="btn-flow" title="Diagrama de flujo por brechas: pocos elementos, fáciles de inspeccionar">RUTAS</button>
    <button id="btn-net" title="Red metabólica completa (cientos de nodos)">RED</button>
  </div>
  <span class="badge" id="b-rxn"></span><span class="badge" id="b-met"></span>
  <span class="badge" id="b-gap"></span><span class="badge" id="b-hyp"></span>
  <span class="badge" id="b-ghost"></span>
  <span class="badge cand" id="b-cand" style="display:none"></span>
  <button id="btn-run">⚙ Ejecutar 08…</button>
</header>
<div id="main"><div id="graph">
  <div id="flow"></div>
  <div id="branches"></div>
  <div id="controls">
    <input id="search" placeholder="buscar metabolito / RHEA / gen" size="26">
    <label><input type="checkbox" id="f-gap"> solo brechas</label>
    <label><input type="checkbox" id="f-hyp"> solo no identificadas</label>
    <label><input type="checkbox" id="f-cur"> mostrar moneda</label>
  </div>
  <svg id="svg"></svg>
  <div id="legend"></div>
  <div id="runbox" style="display:none"></div>
</div><div id="panel"><h2>Panel</h2><div class="dimnote">
Elige un genoma y haz click en un nodo: reacción → genes → segmento de ADN (bp).
El interruptor ANTES/DESPUÉS muestra el efecto del plegado (08).
</div></div></div>
"""

HTML_JS = r"""
// ================== portal ==================
const PORTAL = __PORTAL__;
const C = PORTAL.colors;
const kb = x => x >= 1e6 ? (x/1e6).toFixed(2)+' Mb' : (x/1e3).toFixed(1)+' kb';
const esc = s => String(s??'').replace(/</g,'&lt;');
const $ = id => document.getElementById(id);

// ---------- estado ----------
let S = { genome: null, P: null, cands: [], mode: 'before', view: 'branches',
          built: false,
          candById: {}, candByRxn: {}, loadToken: 0 };
const byId = {};
let nbrs = {}, nodesPos = {};
let alpha = 1, tickN = 0, runToken = 0;
let flowShowAll = false;
let view = { k: 1, x: 0, y: 0 };
let selected = null, searchQ = '', fOnlyGaps = false, fOnlyHyp = false,
    showCurrency = false;

// ---------- leyenda fija ----------
const legendRows = [
  ['metabolito', C.chem_default, 'dot'],
  ['fuente (importable)', C.chem_source, 'dot'],
  ['brecha (sin productor)', C.chem_gap, 'dot'],
  ['sink (producto final)', C.chem_sink, 'dot'],
  ['enzima identificada', C.rxn_ok, 'sq'],
  ['ENZIMA NO IDENTIFICADA (prot. hipotética)', C.rxn_hypothetical, 'sq'],
  ['ENZIMA AUSENTE (fantasma, falta en el genoma)', C.rxn_ghost, 'sq'],
];
$('legend').innerHTML =
  '<b>Leyenda</b><div class="row"><span class="dot" style="background:' +
  C.candidate + ';box-shadow:0 0 6px ' + C.candidate + '"></span>candidato 08 (plegado)</div>'
  + legendRows.map(r =>
      '<div class="row"><span class="' + r[2] + '" style="background:' + r[1] +
      '"></span>' + r[0] + '</div>').join('');

// ---------- desplegable de genomas ----------
(function fillGenomes() {
  const sel = $('genomeSel');
  sel.innerHTML = PORTAL.genomes.map(g =>
    '<option value="' + esc(g) + '">' + esc(g) + '</option>').join('');
  sel.addEventListener('change', () => loadGenome(sel.value));
})();

// ---------- carga de un genoma ----------
async function loadGenome(g) {
  const token = ++S.loadToken;
  let P, cands;
  try {
    if (PORTAL.server) {
      P = await (await fetch('/api/payload/' + encodeURIComponent(g))).json();
      cands = await (await fetch('/api/cands/' + encodeURIComponent(g))).json();
    } else {
      P = PORTAL.data[g];
      cands = PORTAL.cands[g] || { candidates: [] };
    }
  } catch (e) {
    $('panel').innerHTML = '<h2>Error</h2><p class="dimnote">' + esc(e) + '</p>';
    return;
  }
  if (token !== S.loadToken) return;  // ya se pidió otro genoma
  S.genome = g; S.P = P; S.cands = (cands && cands.candidates) || [];
  indexCands();
  indexNodes();          // la vista RUTAS también necesita el índice de nodos
  S.built = false;
  flowShowAll = false;
  updateHeader();
  renderBranches();
  if (S.view === 'flow') renderFlow();
  if (S.view === 'net') rebuildGraph();
  $('panel').innerHTML = '<h2>' + esc(g) + '</h2>' + defaultPanelHtml();
  setMode('before', true);   // siempre aterriza en ANTES; DESPUÉS queda "encendido"
  refreshRunBox();
}

function indexCands() {
  S.candById = {}; S.candByRxn = {};
  S.cands.forEach(c => {
    S.candById[c.id] = c;
    (c.rxn_ids || []).forEach(r => (S.candByRxn[r] ||= []).push(c));
  });
}

function defaultPanelHtml() {
  return '<div class="dimnote">Vista RAMAS: una caja por rama del metabolismo ' +
    '— verde = ruta completa, ámbar = incompleta (toca una brecha). Click en ' +
    'una caja para ver sus pasos (reacción → genes → ADN). Las otras vistas: ' +
    'RUTAS (un bloque por brecha) y RED (grafo completo). El interruptor ' +
    'ANTES/DESPUÉS muestra el efecto del plegado (08).</div>';
}

// ---------- cabecera ----------
function updateHeader() {
  const M = S.P.meta, St = S.P.stats;
  $('dna-total').textContent =
    'ADN: ' + (M.dna_total_length_bp ? kb(M.dna_total_length_bp) + ' (' +
    M.dna_total_length_bp.toLocaleString('es') + ' bp)'
    : (M.genome_meta_status === 'file_missing'
        ? 'n/d (--genome-meta: fichero no encontrado)'
        : 'n/d (falta --genome-meta)'));
  $('dna-flag').style.display = M.reduced_genome ? '' : 'none';
  $('b-rxn').textContent = 'reacciones: ' + St.reactions_in_genome;
  $('b-met').textContent = 'metabolitos: ' + St.metabolites;
  $('b-gap').textContent = 'brechas: ' + St.gaps;
  $('b-hyp').textContent = 'hipotéticas: ' + St.hypothetical_rxns;
  $('b-ghost').textContent = 'enzimas ausentes: ' + St.ghost_rxns;
  updateCandBadge();
}
function candCounts() {
  const k = { fuerte: 0, validar: 0, descartado: 0, otro: 0 };
  S.cands.forEach(c => { k[c.verdict in k ? c.verdict : 'otro']++; });
  return k;
}
function updateCandBadge() {
  const b = $('b-cand');
  const k = candCounts();
  const has = S.cands.length > 0;
  $('btn-after').disabled = !has;
  $('btn-after').classList.toggle('ready', has);
  if (!has || S.mode !== 'after') { b.style.display = 'none'; return; }
  b.style.display = '';
  b.textContent = 'candidatos 08 — fuertes: ' + k.fuerte +
    ' · a validar: ' + k.validar + ' · descartados: ' + k.descartado;
}

// ---------- interruptor ANTES/DESPUÉS ----------
function setMode(m, force) {
  if (!force && S.mode === m) return;
  if (m === 'after' && !S.cands.length) return;   // sin resultados aún
  S.mode = m;
  $('btn-before').classList.toggle('active', m === 'before');
  $('btn-after').classList.toggle('active', m === 'after');
  if (S.view === 'net') draw();
  else if (S.view === 'branches') renderBranches();
  else renderFlow();
  if (selected) select(selected);
  updateCandBadge();
}
$('btn-before').addEventListener('click', () => setMode('before'));
$('btn-after').addEventListener('click', () => setMode('after'));

// ---------- vistas (RAMAS por defecto; RUTAS = flujo por brechas; RED) ----------
function setView(v) {
  if (S.view === v) return;
  S.view = v;
  ['branches', 'flow', 'net'].forEach(x =>
    $('btn-' + x).classList.toggle('active', v === x));
  $('flow').style.display = v === 'flow' ? 'block' : 'none';
  $('branches').style.display = v === 'branches' ? 'block' : 'none';
  $('svg').style.display = v === 'net' ? '' : 'none';
  $('controls').style.display = v === 'net' ? '' : 'none';
  $('legend').style.display = v === 'net' ? '' : 'none';
  if (v === 'net') { if (!S.built) rebuildGraph(); }
  else if (v === 'branches') renderBranches();
  else renderFlow();
}
$('btn-branches').addEventListener('click', () => setView('branches'));
$('btn-flow').addEventListener('click', () => setView('flow'));
$('btn-net').addEventListener('click', () => setView('net'));

// ---------- vista RAMAS (cajas por rama del metabolismo) ----------
function renderBranches() {
  const host = $('branches');
  if (!host) return;
  const B = S.P.branches || [];
  const q = searchQ.toLowerCase();
  let h = '<div class="fnote" style="padding:12px 16px 0">RAMAS del metabolismo: ' +
    'una caja por rama (componente de la red tras separar los metabolitos hub ' +
    'que cosen rutas distintas). <b style="color:#4ade80">Verde</b> = ruta ' +
    'completa (ninguna de sus reacciones toca una brecha), <b style="color:#fbbf24">' +
    'ámbar</b> = incompleta. Click en una caja: pasos, genes y ADN.' +
    (S.mode === 'after' ? ' Con DESPUÉS activo, el halo amarillo marca pasos con candidato del 08.' : '') +
    '</div>';
  if (!B.length) {
    h += '<div class="fnote" style="padding:12px 16px">Sin ramas: la red ' +
      'queda en un solo componente o todo son hubs. Baja el corte con ' +
      '--branch-hub-degree o mira la pestaña RED.</div>';
    host.innerHTML = h;
    return;
  }
  const cards = B.filter(b => !q ||
    (b.name + ' ' + b.reactions.join(' ') + ' ' + b.gaps.join(' '))
      .toLowerCase().includes(q));
  h += '<div class="bgrid">' + cards.map(b => {
    const st = b.complete ? 'ok' : 'in';
    const miss = b.gaps && b.gaps.length
      ? '<div class="bmiss">falta: ' + b.gaps.map(g =>
          byId[g] ? esc(byId[g].label) : esc(g)).join(', ') + '</div>'
      : '';
    return '<div class="bcard ' + st + '" data-b="' + esc(b.id) + '">' +
      '<div class="bname">' + esc(b.name) + '</div>' +
      '<div class="bmeta"><span class="tag ' + st + '">' +
      (b.complete ? 'COMPLETA' : 'INCOMPLETA') + '</span> ' +
      b.reactions.length + ' pasos' +
      (b.n_genes ? ' · ' + b.n_genes + ' gen' + (b.n_genes > 1 ? 'es' : '') : '') +
      (b.n_gaps ? ' · ' + b.n_gaps + ' brecha' + (b.n_gaps > 1 ? 's' : '') : '') +
      miss + '</div></div>';
  }).join('') + '</div>';
  if (!cards.length) h += '<div class="fnote" style="padding:0 16px">Ninguna ' +
    'rama coincide con la búsqueda.</div>';
  host.innerHTML = h;
  host.querySelectorAll('.bcard').forEach(el =>
    el.addEventListener('click', () => branchDetail(el.dataset.b)));
}

function branchDetail(bid) {
  const b = (S.P.branches || []).find(x => x.id === bid);
  if (!b) return;
  let h = '<h2>' + esc(b.name) + '</h2>' +
    '<div style="margin:6px 0"><span class="tag ' + (b.complete ? 'ok' : 'in') +
    '">' + (b.complete ? 'RUTA COMPLETA' : 'RUTA INCOMPLETA') + '</span> ' +
    b.reactions.length + ' pasos' +
    (b.n_genes ? ' · ' + b.n_genes + ' gen' + (b.n_genes > 1 ? 'es' : '') : '') +
    '</div>';
  if (b.gaps && b.gaps.length) {
    h += '<div class="gene" style="border-color:#b91c1c"><div>Metabolitos en ' +
      'brecha dentro de la rama (se consumen y nada del genoma los produce): <b>' +
      b.gaps.map(g => byId[g] ? esc(byId[g].label) : esc(g)).join(', ') +
      '</b>. Falta la enzima que los produce — críala con el 08 (⚙).</div></div>';
  }
  h += '<h2 style="margin-top:10px">Pasos de la ruta</h2><div>' +
    b.reactions.map(r => fboxHtml(r)).join('<div class="farrow">↓</div>') +
    '</div><p class="dimnote" style="margin-top:8px">Click en un paso para ' +
    'reacción → genes → ADN.</p>';
  const panel = $('panel');
  panel.innerHTML = h;
  panel.scrollTop = 0;
  panel.querySelectorAll('.fbox').forEach(el =>
    el.addEventListener('click', ev => {
      ev.stopPropagation(); select(el.dataset.id);
    }));
}

function fboxHtml(id) {
  const n = byId[id];
  if (!n) return '';
  const cls = n.kind === 'ghost' ? 'ghost' : (n.hypothetical ? 'hyp' : 'ok');
  const halo = (S.mode === 'after' && haloFor(n)) ? ' halo' : '';
  const label = n.kind === 'ghost'
    ? '<span class="rh">' + esc(n.rhea || '') + '</span> <span class="fnote">' +
      'ausente — produciría la brecha</span>'
    : '<span class="rh">' + esc(n.rhea || id) + '</span>' +
      (n.equation ? ' <span class="fnote">' + esc(n.equation) + '</span>' : '');
  return '<span class="fbox ' + cls + halo + '" data-id="' + esc(id) + '">' +
    label + '</span>';
}

function flowRow(label, inner) {
  return '<div class="frow"><span class="flabel">' + label + '</span>' +
    inner + '</div>';
}

function renderFlow() {
  const host = $('flow');
  const gaps = S.P.nodes.filter(n => n.kind === 'chem' && n.role === 'gap');
  gaps.sort((a, b) => ((b.systematic ? 1 : 0) - (a.systematic ? 1 : 0)) ||
    String(a.label).localeCompare(String(b.label)));
  const inE = {}, outE = {};
  S.P.edges.forEach(e => {
    (inE[e[1]] ||= []).push(e[0]);
    (outE[e[0]] ||= []).push(e[1]);
  });
  const feeders = id => (inE[id] || []).filter(x => byId[x] && byId[x].kind !== 'chem');
  const consumers = id => (outE[id] || []).filter(x => byId[x] && byId[x].kind === 'rxn');
  const q = searchQ.toLowerCase();
  const FLOW_MAX = 30;   // con genomas reales hay cientos de brechas
  const lim = (q || flowShowAll) ? gaps.length : Math.min(FLOW_MAX, gaps.length);
  let h = '<div class="fnote">Vista RUTAS — un bloque por brecha (metabolito ' +
    'sin enzima que lo produzca en el genoma): pocos elementos, inspección ' +
    'directa. Click en una caja para reacción → genes → ADN; la red completa ' +
    'está en la pestaña RED.' +
    (S.mode === 'after' ? ' Halo amarillo = reacción con candidato del 08.' :
     ' Tras ejecutar el 08, pasa a DESPUÉS para ver los candidatos.') + '</div>';
  if (gaps.length > FLOW_MAX && !q) {
    h += '<div class="fnote" style="margin:6px 0">' + gaps.length +
      ' brechas: se muestran las ' + lim + ' más relevantes ' +
      '(sistemáticas primero). Filtra con el buscador o ' +
      (flowShowAll ? '' : '<button class="fbox" id="flow-all">mostrar todas</button>') +
      '</div>';
  }
  let shown = 0;
  gaps.slice(0, lim).forEach((g, i) => {
    const ghosts = feeders(g.id).filter(id => byId[id].kind === 'ghost');
    const prods = feeders(g.id).filter(id => byId[id].kind === 'rxn');
    const cons = consumers(g.id);
    const hay = [g.label, g.id, g.triage || '', ...ghosts, ...prods, ...cons]
      .join(' ').toLowerCase();
    if (q && !hay.includes(q)) return;
    shown++;
    const chip = '<span class="fmet' + (g.systematic ? ' sys' : '') +
      '" data-id="' + esc(g.id) + '">' + esc(g.label) + ' · ' + esc(g.id) +
      '</span>';
    h += '<div class="fsec"><div class="fhead"><b>Brecha ' + (i + 1) + '/' +
      gaps.length + '</b>' + chip +
      (g.triage ? '<span class="tag">' + esc(g.triage) + '</span>' : '') +
      (g.systematic ? '<span class="tag sys">SISTEMÁTICA (05)</span>' : '') +
      '<span class="fnote" style="margin-left:auto">' + cons.length +
      ' consumidor(es)</span></div><div class="fbody">' +
      flowRow('enzima ausente',
        ghosts.length ? ghosts.map(fboxHtml).join('<span class="farrow">+</span>')
        : '<span class="fnone">ninguna reacción de Rhea lo produciría</span>') +
      flowRow('brecha', chip) +
      (prods.length ? flowRow('productor (en el genoma)',
        prods.map(fboxHtml).join('<span class="farrow">+</span>')) : '') +
      flowRow('consumidores',
        cons.length ? cons.map(fboxHtml).join('<span class="farrow">+</span>')
        : '<span class="fnone">sin reacción consumidora en el genoma</span>') +
      '</div></div>';
  });
  if (!shown) h += '<div class="fnote">Sin brechas' +
    (q ? ' que coincidan con la búsqueda' : '') + '.</div>';
  host.innerHTML = h;
  host.querySelectorAll('.fbox, .fmet').forEach(el =>
    el.addEventListener('click', ev => {
      ev.stopPropagation(); select(el.dataset.id);
    }));
  const allBtn = host.querySelector('#flow-all');
  if (allBtn) allBtn.addEventListener('click', () => {
    flowShowAll = true; renderFlow();
  });
  host.querySelectorAll('.fhead').forEach(el =>
    el.addEventListener('click', () => {
      const b = el.nextElementSibling;
      if (b) b.style.display = b.style.display === 'none' ? '' : 'none';
    }));
}

// ---------- grafo ----------
function indexNodes() {
  for (const k in byId) delete byId[k];
  S.P.nodes.forEach(n => byId[n.id] = n);
}

function rebuildGraph() {
  S.built = true;
  indexNodes();
  nbrs = {}; nodesPos = {}; selected = null;
  S.P.edges.forEach(e => {
    (nbrs[e[0]] ||= new Set()).add(e[1]);
    (nbrs[e[1]] ||= new Set()).add(e[0]);
  });
  const svg = $('svg');
  const W = svg.clientWidth, H = svg.clientHeight;
  const rxs = S.P.nodes.filter(n => n.kind !== 'chem');
  const chs = S.P.nodes.filter(n => n.kind === 'chem');
  rxs.forEach((n, i) => {
    const a = 2*Math.PI*i/Math.max(rxs.length,1);
    nodesPos[n.id] = { x: W/2 + 120*Math.cos(a), y: H/2 + 120*Math.sin(a), vx:0, vy:0 };
  });
  chs.forEach((n, i) => {
    const a = 2*Math.PI*i/Math.max(chs.length,1) + 0.7;
    nodesPos[n.id] = { x: W/2 + 300*Math.cos(a), y: H/2 + 300*Math.sin(a), vx:0, vy:0 };
  });
  alpha = 1; tickN = 0;
  draw();
  const token = ++runToken;
  (function step() {
    if (token !== runToken || alpha < 0.005 || tickN > 320) return;
    tickBody(); draw(); requestAnimationFrame(step);
  })();
}

function tickBody() {
  const ids = Object.keys(nodesPos);
  const step = ids.length > 1200 ? 2 : 1;
  for (let i = 0; i < ids.length; i++) {
    const a = nodesPos[ids[i]];
    for (let j = i+1; j < ids.length; j += step) {
      const b = nodesPos[ids[j]];
      let dx = a.x-b.x, dy = a.y-b.y, d2 = dx*dx+dy*dy || 1;
      if (d2 > 250000) continue;
      const f = 900*alpha/d2, d = Math.sqrt(d2);
      dx /= d; dy /= d;
      a.vx += dx*f; a.vy += dy*f; b.vx -= dx*f; b.vy -= dy*f;
    }
  }
  S.P.edges.forEach(e => {
    const a = nodesPos[e[0]], b = nodesPos[e[1]];
    if (!a || !b) return;
    let dx = b.x-a.x, dy = b.y-a.y, d = Math.sqrt(dx*dx+dy*dy) || 1;
    const f = (d-70)*0.02*alpha;
    a.vx += dx/d*f; a.vy += dy/d*f; b.vx -= dx/d*f; b.vy -= dy/d*f;
  });
  ids.forEach(id => {
    const p = nodesPos[id];
    p.vx += (window.innerWidth/2 - p.x) * 0.0015 * alpha;
    p.vy += (window.innerHeight/2 - p.y) * 0.0015 * alpha;
    p.x += Math.max(-8, Math.min(8, p.vx*8));
    p.y += Math.max(-8, Math.min(8, p.vy*8));
    p.vx *= 0.6; p.vy *= 0.6;
  });
  alpha *= 0.985; tickN++;
}

// ---------- visibilidad y render ----------
function nodeColor(n) {
  if (n.kind === 'chem')
    return n.role === 'gap' ? C.chem_gap : n.role === 'source' ? C.chem_source :
           n.role === 'sink' ? C.chem_sink : C.chem_default;
  if (n.kind === 'ghost') return C.rxn_ghost;
  return n.hypothetical ? C.rxn_hypothetical : C.rxn_ok;
}
function haloFor(n) {
  if (S.mode !== 'after') return false;
  if (n.rhea && S.candByRxn[n.rhea]) return true;
  return (n.genes || []).some(g => S.candById[g.gene_id]);
}
function visible(n) {
  if (!showCurrency && n.currency) return false;
  if (fOnlyGaps && !(n.role === 'gap' || n.gap_for || n.kind === 'ghost')) return false;
  if (fOnlyHyp && !(n.hypothetical || n.kind === 'ghost')) return false;
  if (searchQ) {
    const q = searchQ.toLowerCase();
    const hay = [n.label, n.rhea, n.id,
      (n.genes||[]).map(g => g.gene_id+' '+g.locus_tag+' '+g.product).join(' ')]
      .join(' ').toLowerCase();
    if (!hay.includes(q)) {
      let near = false;
      (nbrs[n.id]||new Set()).forEach(o => {
        const on = byId[o];
        if (on && [on.label, on.rhea, on.id].join(' ').toLowerCase().includes(q))
          near = true;
      });
      if (!near) return false;
    }
  }
  return true;
}

function draw() {
  const svg = $('svg');
  let parts = [];
  S.P.edges.forEach(e => {
    const na = byId[e[0]], nb = byId[e[1]];
    if (!na || !nb || !visible(na) || !visible(nb)) return;
    const a = nodesPos[e[0]], b = nodesPos[e[1]];
    if (!a || !b) return;
    const ghost = e[2] === 'ghost';
    parts.push('<line x1="'+a.x+'" y1="'+a.y+'" x2="'+b.x+'" y2="'+b.y+
      '" stroke="'+(ghost ? C.rxn_ghost : '#223252')+
      '" stroke-width="'+(ghost ? 1.6 : 1)+'"'+
      (ghost ? ' stroke-dasharray="5,4"' : '')+' opacity="0.8"/>');
  });
  S.P.nodes.forEach(n => {
    if (!visible(n)) return;
    const p = nodesPos[n.id]; if (!p) return;
    const col = nodeColor(n);
    const halo = haloFor(n);
    if (n.kind === 'chem') {
      parts.push('<circle cx="'+p.x+'" cy="'+p.y+'" r="'+(n.role==='gap'?8:6)+
        '" fill="'+col+'" stroke="'+(n.systematic?'#fbbf24':'#0b1220')+
        '" stroke-width="'+(n.systematic?2.5:1)+'" data-id="'+n.id+'" class="nd"/>');
    } else {
      const s = n.kind === 'ghost' ? 10 : 8;
      parts.push('<g data-id="'+n.id+'" class="nd"><rect x="'+(p.x-s)+'" y="'+(p.y-s)+
        '" width="'+2*s+'" height="'+2*s+'" transform="rotate(45 '+p.x+' '+p.y+
        ')" fill="'+col+'" stroke="'+(halo ? C.candidate : '#0b1220')+
        '" stroke-width="'+(halo ? 3 : 1)+'"'+(n.kind==='ghost'?' stroke-dasharray="4,3"':'')+
        (halo ? ' class="flash"' : '')+'/><text x="'+p.x+'" y="'+(p.y+s+11)+
        '" fill="'+col+'" font-size="9" text-anchor="middle">'+
        (n.rhea||n.id).replace('RHEA:','R')+'</text></g>');
    }
  });
  svg.innerHTML = parts.join('');
  svg.style.transform = 'translate(' + view.x + 'px,' + view.y + 'px) scale(' + view.k + ')';
  svg.style.transformOrigin = '0 0';
  svg.querySelectorAll('.nd').forEach(el => {
    el.addEventListener('click', ev => { ev.stopPropagation(); select(el.dataset.id); });
  });
}

// ---------- panel ----------
function fmtBP(bp, total) {
  const t = bp >= 1000 ? (bp/1000).toFixed(2)+' kb' : bp+' bp';
  const pct = total ? ' (' + (100*bp/total).toFixed(2) + '% del genoma)' : '';
  return '<span class="len">' + t + '</span>' + pct;
}
function candTag(c) {
  return '<span class="tag ' + esc(c.verdict) + '">' + esc(c.verdict) +
    ' · pLDDT ' + (c.mean_plddt != null ? Math.round(c.mean_plddt) : '—') +
    ' · dominio ' + (c.structured_run != null ? c.structured_run + ' aa' : '—') +
    '</span>';
}
function candCard(c) {
  const bh = c.best_hit || {};
  return '<div class="gene">' +
    '<div><b>' + esc(c.id) + '</b> <span class="dimnote">(' + esc(c.kind) + ')</span></div>' +
    candTag(c) +
    '<div class="dimnote">' + esc(c.reason || '') + '</div>' +
    '<div>🧬 ADN: ' + esc(c.contig || 'n/d') + ':' + (c.start ?? '—') + '-' + (c.end ?? '—') +
    ' (' + esc(c.strand || '?') + ') · ' + (c.length_bp || '—') + ' bp</div>' +
    (bh.uniprot ? '<div>🎯 best hit: ' + esc(bh.uniprot) + ' (' + bh.pident + '% id, cob. ' +
      Math.round(100*(bh.coverage||0)) + '%) → brecha ' + esc(bh.gap_chebi || 'n/d') + '</div>' : '') +
    '</div>';
}
function geneCard(g) {
  const cands = S.mode === 'after' && S.candById[g.gene_id]
    ? [S.candById[g.gene_id]] : [];
  return '<div class="gene">' +
    '<div><b>' + esc(g.gene_id) + '</b>' + (g.locus_tag ? ' · ' + esc(g.locus_tag) : '') + '</div>' +
    '<div class="dimnote">' + esc(g.product || 'producto no anotado') + '</div>' +
    '<div>📍 ' + esc(g.contig || 'contig n/d') + ':' + g.start + '-' + g.end +
    ' (' + esc(g.strand || '?') + ')</div>' +
    '<div>🧬 tamaño del segmento: ' + fmtBP(g.length_bp, S.P.meta.dna_total_length_bp) + '</div>' +
    cands.map(candCard).join('') +
    '</div>';
}
function dnaStrip(g0) {
  const same = [];
  S.P.nodes.forEach(n => (n.genes||[]).forEach(g => {
    if (g.contig && g0.contig && g.contig === g0.contig && g.start <= g0.end+6000
        && g.end >= g0.start-6000) same.push(g);
  }));
  if (same.length < 2) return '<div class="dna-strip dimnote">Sin vecinos mapeados ' +
    'en la misma región de ADN.</div>';
  const lo = Math.min(...same.map(g => g.start)), hi = Math.max(...same.map(g => g.end));
  const span = hi-lo || 1, wpx = 340;
  const sel = g => g.gene_id === g0.gene_id;
  const strip = same.sort((a,b) => a.start-b.start).map(g =>
    '<rect x="' + (20+(g.start-lo)/span*wpx) + '" y="14" width="' +
    Math.max(3,(g.end-g.start)/span*wpx) + '" height="16" rx="2" fill="' +
    (sel(g) ? C.rxn_hypothetical : '#334766') + '"/>' +
    '<text x="' + (20+(g.start-lo)/span*wpx) + '" y="42" font-size="8" fill="' +
    (sel(g) ? '#fdba74' : '#7c8aa5') + '" text-anchor="middle">' +
    (g.end-g.start+1) + ' bp</text>').join('');
  return '<div class="dna-strip"><div class="dimnote">Contexto del contig ' +
    esc(g0.contig) + ' (±6 kb; ancho ∝ bp):</div><svg width="' + (wpx+40) +
    '" height="48">' + strip + '</svg></div>';
}
function candidatesForRxn(n) {
  if (S.mode !== 'after') return [];
  let out = (n.rhea && S.candByRxn[n.rhea]) ? S.candByRxn[n.rhea].slice() : [];
  if (n.gap_for) {
    S.cands.forEach(c => {
      if (c.best_hit && c.best_hit.gap_chebi === n.gap_for && !out.includes(c))
        out.push(c);
    });
  }
  return out;
}
function select(id) {
  const n = byId[id];
  if (!n) return;
  selected = id;
  let h = '';
  if (n.kind === 'chem') {
    h += '<h2>' + esc(n.label) + '</h2><div class="dimnote">' + esc(n.id) + '</div>';
    const tags = [];
    if (n.role !== 'none') tags.push('<span class="tag ' + n.role + '">' +
      ({gap:'BRECHA (sin enzima que lo produzca)', source:'FUENTE (importable)',
        sink:'SINK (producto final)'}[n.role]) + '</span>');
    if (n.triage) tags.push('<span class="tag">' + esc(n.triage) + '</span>');
    if (n.systematic) tags.push('<span class="tag sys">SISTEMÁTICA (cohorte 05)</span>');
    h += '<div style="margin:8px 0">' + tags.join('') + '</div>';
    const prod = S.P.edges.filter(e => e[1] === id && e[2] !== 'ghost').map(e => e[0]);
    const prodG = S.P.edges.filter(e => e[1] === id && e[2] === 'ghost').map(e => e[0]);
    if (prodG.length) h += '<p>⚡ <b>Enzimas ausentes</b> que lo producirían: ' +
      prodG.map(r => byId[r] ? esc(byId[r].rhea) : r).join(', ') + '</p>';
    if (prod.length) h += '<p>Producido en el genoma por: ' +
      prod.filter(r => byId[r] && byId[r].kind === 'rxn').map(r => esc(r)).join(', ') + '</p>';
    const cands = candidatesForRxn({gap_for: n.id, rhea: null});
    if (cands.length) h += '<h2 style="margin-top:10px">Candidatos 08 (plegado)</h2>' +
      cands.map(candCard).join('');
  } else {
    const ghost = n.kind === 'ghost';
    h += '<h2>' + (ghost ? '👻 ' : '⚙ ') + esc(n.rhea) +
      (ghost ? ' — ENZIMA AUSENTE' : '') + '</h2>';
    if (n.gap_for) h += '<div>Resuelve la brecha: <b>' + esc(n.gap_name) +
      '</b> (' + esc(n.gap_for) + ')' +
      (n.systematic ? ' <span class="tag sys">SISTEMÁTICA</span>' : '') + '</div>';
    if (n.equation) h += '<p class="eq">' + esc(n.equation) + '</p>';
    if (ghost) {
      h += '<div class="gene"><div>Esta reacción <b>no está en el genoma</b>: ' +
        'es la enzima faltante que produciría el metabolito-brecha. ' +
        (S.mode === 'after'
          ? 'Abajo, las candidatas del 08 cuyo best hit explica esta brecha.'
          : 'Ejecuta el 08 (botón ⚙) para cribar hipotéticas y ORFs.') + '</div></div>';
    } else if (n.hypothetical) {
      h += '<div class="gene" style="border-color:#f97316"><div>' +
        '<span class="tag hyp">ENZIMA NO IDENTIFICADA</span> codificada solo por ' +
        'proteína(s) hipotética(s) — verificar función por plegado (08).</div></div>';
    }
    const cands = candidatesForRxn(n);
    if (cands.length) h += '<h2 style="margin-top:10px">Candidatos 08 (plegado)</h2>' +
      cands.map(candCard).join('');
    if (n.genes && n.genes.length) {
      h += '<h2 style="margin-top:10px">Genes → ADN</h2>' +
        n.genes.map(geneCard).join('');
      h += dnaStrip(n.genes[0]);
    } else if (!ghost) {
      h += '<p class="dimnote">Sin procedencia génica: falta ' +
        esc(S.genome) + '_gene_map.tsv en --gene-map.</p>';
    }
  }
  const panel = $('panel');
  panel.innerHTML = h;
  panel.scrollTop = 0;
}

// ================== ejecutar 08 desde la interfaz ==================
let pollTimer = null;
let runOpen = false, runMin = false;   // cajón cerrado por defecto: no estorba

function refreshRunBox() {
  const box = $('runbox');
  if (!PORTAL.server || !S.genome) { box.style.display = 'none'; return; }
  box.style.display = runOpen ? '' : 'none';
  box.classList.toggle('min', runMin);
  if (box.dataset.built !== '1') {
    box.dataset.built = '1';
    box.innerHTML =
      '<div class="rhead" id="run-head"><b>Ejecutar 08 en <span id="run-title"></span></b>' +
      '<span class="min" id="run-min">▾</span></div>' +
      '<div id="runform">' +
      '<label><input type="checkbox" id="run-orfs" style="width:auto"> incluir ORFs de los contigs</label>' +
      '<label>similitud</label><select id="run-tool">' +
      '<option value="auto">auto (mmseqs/blast si hay)</option>' +
      '<option value="none">sin similitud (solo estructura)</option>' +
      '<option value="mmseqs">mmseqs2</option><option value="blastp">blastp</option></select>' +
      '<label>máx. secuencias a plegar</label>' +
      '<input id="run-maxfold" type="number" value="50" min="1" max="500">' +
      '<label>servicio de plegado ESM</label><select id="run-provider">' +
      '<option value="auto">auto (según entorno/.env: si hay key -> Neurosnap)</option>' +
      '<option value="neurosnap">Neurosnap ESMFold (requiere API key)</option>' +
      '<option value="direct">endpoint directo (ESM_API_URL)</option></select>' +
      '<label>ESM_API_URL (opcional)</label>' +
      '<input id="run-url" type="text" placeholder="https://tu-endpoint/fold">' +
      '<label>ESM_API_KEY (opcional; solo para el subproceso)</label>' +
      '<input id="run-key" type="password" placeholder="••••••••">' +
      '<div class="fnote">Consejo: escribe la clave una vez en el archivo ' +
      '<b>.env</b> de la raíz del proyecto (ESM_API_KEY=...) y déjalo vacío aquí.</div>' +
      '<div class="row2"><button id="run-go" style="flex:1;background:#14532d;color:#dcfce7;' +
      'border:1px solid #166534;border-radius:4px;padding:5px;font-size:12px;cursor:pointer">' +
      'Lanzar</button><button id="run-hide" style="background:var(--panel);color:var(--dim);' +
      'border:1px solid var(--line);border-radius:4px;padding:5px;cursor:pointer">✕</button></div>' +
      '</div>' +
      '<div id="runlog" style="display:none"></div>';
    $('run-hide').addEventListener('click', () => { runOpen = false; box.style.display = 'none'; });
    $('run-go').addEventListener('click', launch08);
    $('run-head').addEventListener('click', () => {
      runMin = !runMin;
      box.classList.toggle('min', runMin);
      $('run-min').textContent = runMin ? '▸' : '▾';
    });
  }
  $('run-title').textContent = S.genome;
  if (runMin) $('run-min').textContent = '▸';
}

$('btn-run').addEventListener('click', () => {
  if (!PORTAL.server) {
    alert('Modo estático: para ejecutar el 08 desde la interfaz lanza\n\n' +
          '  python 07_network_view.py\n\ny abre http://127.0.0.1:8765');
    return;
  }
  runOpen = !runOpen;
  refreshRunBox();
});

async function launch08() {
  const body = {
    genome: S.genome,
    orfs: $('run-orfs').checked,
    tool: $('run-tool').value,
    max_fold: parseInt($('run-maxfold').value || '50', 10),
    esm_provider: $('run-provider').value,
  };
  const k = $('run-key').value.trim(), u = $('run-url').value.trim();
  if (k) body.esm_api_key = k;
  if (u) body.esm_api_url = u;
  const log = $('runlog');
  log.style.display = '';
  log.textContent = 'Lanzando 08 para ' + S.genome + '…\n';
  const r = await fetch('/api/run08', { method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body) });
  if (!r.ok) {
    const j = await r.json().catch(() => ({}));
    log.textContent += 'ERROR: ' + (j.error || r.status) + '\n';
    return;
  }
  $('run-go').textContent = 'En marcha…';
  pollTimer = setInterval(poll08, 1200);
  poll08();
}

async function poll08() {
  const st = await (await fetch('/api/status')).json();
  const log = $('runlog');
  if (st.log && st.log.length) log.textContent = st.log.slice(-60).join('');
  log.scrollTop = log.scrollHeight;
  if (st.running) return;
  clearInterval(pollTimer); pollTimer = null;
  $('run-go').textContent = 'Lanzar';
  if (!st.started) return;
  const g = st.genome || S.genome;
  const cands = await (await fetch('/api/cands/' + encodeURIComponent(g))).json();
  if (g === S.genome) {   // solo si el usuario sigue viendo ese genoma
    S.cands = (cands && cands.candidates) || [];
    indexCands();
    setMode(S.cands.length ? 'after' : 'before', true);
    const k = candCounts();
    log.textContent += '\n== Completado (rc=' + st.rc + '). Candidatos: fuertes ' +
      k.fuerte + ' · a validar ' + k.validar + ' · descartados ' + k.descartado +
      ' — vista DESPUÉS activada ==';
    if (S.cands.length) draw();
  } else {
    log.textContent += '\n== Completado (rc=' + st.rc + ') para ' + esc(g) +
      ': vuelve a seleccionarlo para ver los candidatos ==';
  }
}

// ================== interacción del grafo ==================
const graphEl = document.getElementById('graph');
let drag = null;
graphEl.addEventListener('mousedown', e => {
  if (S.view !== 'net' || e.target.closest('#flow') || e.target.closest('#runbox'))
    return;   // en RUTAS el ratón es scroll/selección, no pan del grafo
  drag = { x: e.clientX, y: e.clientY, vx: view.x, vy: view.y };
});
window.addEventListener('mousemove', e => {
  if (drag) { view.x = drag.vx + (e.clientX-drag.x); view.y = drag.vy + (e.clientY-drag.y); draw(); }
});
window.addEventListener('mouseup', () => drag = null);
graphEl.addEventListener('wheel', e => {
  if (S.view !== 'net' || e.target.closest('#flow') || e.target.closest('#runbox'))
    return;   // deja que la vista RUTAS haga scroll nativo
  e.preventDefault();
  const k2 = view.k * (e.deltaY < 0 ? 1.15 : 0.87);
  view.x = e.offsetX - (e.offsetX-view.x)*k2/view.k;
  view.y = e.offsetY - (e.offsetY-view.y)*k2/view.k;
  view.k = k2; draw();
}, { passive: false });

$('search').addEventListener('input', e => {
  searchQ = e.target.value.trim();
  if (S.view === 'flow') renderFlow(); else draw();
});
$('f-gap').addEventListener('change', e => { fOnlyGaps = e.target.checked; draw(); });
$('f-hyp').addEventListener('change', e => { fOnlyHyp = e.target.checked; draw(); });
$('f-cur').addEventListener('change', e => { showCurrency = e.target.checked; draw(); });

// ---------- arranque ----------
// aplica el estado de la vista por defecto (RAMAS) al DOM; sin esto el
// contenedor queda con el display:none del CSS y la pantalla en blanco
(function initView() {
  $('branches').style.display = 'block';   // el CSS base lo trae oculto
  $('flow').style.display = 'none';
  $('svg').style.display = 'none';
  $('controls').style.display = 'none';
  $('legend').style.display = 'none';
})();
loadGenome($('genomeSel').value);
"""

HTML_TAIL = "</script></body></html>\n"


def _jsafe(obj):
    """JSON seguro para incrustar en <script>."""
    return json.dumps(obj, ensure_ascii=False).replace("</", "<\\/")


def render_portal(payloads, cands_by_genome, server, genomes=None):
    portal = {
        "server": bool(server),
        "genomes": sorted(genomes if genomes is not None
                          else payloads.keys()),
        "colors": COLORS,
    }
    if not server:
        portal["data"] = payloads
        portal["cands"] = cands_by_genome
    body = HTML_HEAD + HTML_BODY + "<script>" + HTML_JS + HTML_TAIL
    return body.replace("__PORTAL__", _jsafe(portal))


# -------------------------------------------------------------- servidor
class _FoldJob:
    """Estado del subproceso 08 (uno a la vez)."""

    def __init__(self):
        self.lock = threading.Lock()
        self.proc = None
        self.genome = None
        self.log = []
        self.rc = None
        self.started = False

    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, cmd, env, genome):
        with self.lock:
            if self.proc is not None and self.proc.poll() is None:
                return False
            self.genome = genome
            self.log = []
            self.rc = None
            self.started = True
            self.proc = subprocess.Popen(
                cmd, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1)
        threading.Thread(target=self._pump, daemon=True).start()
        return True

    def _pump(self):
        assert self.proc and self.proc.stdout
        for line in self.proc.stdout:
            with self.lock:
                self.log.append(line)
                if len(self.log) > 2000:
                    self.log = self.log[-2000:]
        with self.lock:
            self.rc = self.proc.wait()

    def stop(self):
        with self.lock:
            if self.proc and self.proc.poll() is None:
                self.proc.terminate()
                return True
            return False

    def status(self):
        with self.lock:
            return {
                "running": self.proc is not None and self.proc.poll() is None,
                "started": self.started,
                "genome": self.genome,
                "rc": self.rc,
                "log": self.log[-60:],
            }


class Portal:
    """Estado compartido del servidor (payloads cacheados + job de 08)."""

    def __init__(self, args):
        self.args = args
        self.net_dir = args.net_dir
        self.rhea_dir = args.rhea_dir
        self._payloads = {}
        self._lock = threading.Lock()
        self._gm_tried = set()          # intentos de auto-genemap por genoma
        self._fna_memo = {}             # (path, mtime, size) -> (bp, gc%)
        self.job = _FoldJob()
        if getattr(args, "genome_meta", None) and not args.genome_meta.exists():
            m04.log.warning(f"--genome-meta: no existe {args.genome_meta} — "
                            "el badge ADN quedará n/d (no afecta al gene map)")
        self.systematic = load_cohort(args.cohort or args.net_dir / "trace_compare.json")

    def genomes(self):
        return sorted(p.name[:-len("_network.json")]
                      for p in self.net_dir.glob("*_network.json"))

    def _gene_map_path(self, genome):
        if not self.args.gene_map:
            return None
        for suffix in (f"{genome}_gene_map.tsv", f"{genome}_gene_map.json"):
            p = self.args.gene_map / suffix
            if p.exists():
                return p
        if self.args.gene_map.suffix in (".tsv", ".json"):
            return self.args.gene_map
        return None

    def _ensure_gene_map(self, genome):
        """Integración con el 09: si falta {genome}_gene_map.tsv|json (o
        existe pero está VACÍO: salió sin ECs y luego llegó la anotación
        buena, p. ej. tras descargar los .gbff con 10_fetch_gbff.py) y hay
        anotación descubrible, lo genera al vuelo una única vez por genoma."""
        gm = self._gene_map_path(genome)
        if gm and load_gene_map(gm):
            return gm                     # existe y tiene filas
        if not self.args.gene_map or self.args.no_auto_genemap \
                or genome in self._gm_tried:
            return gm
        if gm:
            m04.log.info(f"[{genome}] gene map vacío presente: regenerando "
                         "con el 09")
        self._gm_tried.add(genome)
        try:
            m09 = importlib.import_module("09_make_gene_map")
            out = m09.build_gene_map(
                genome, annotations_dir=self.args.annotations_dir,
                rhea_dir=self.rhea_dir, net_dir=self.net_dir,
                out_dir=self.args.gene_map)
            m04.log.info(f"[{genome}] gene map autogenerado por 09: {out}")
        except SystemExit as e:
            m04.log.warning(f"[{genome}] sin gene map autogenerado: {e}")
        except Exception:  # noqa: BLE001 — el portal no debe caer por esto
            m04.log.exception(f"[{genome}] auto-genemap falló")
        return self._gene_map_path(genome)

    def _fna_meta(self, genome):
        """Badge ADN sin genome_meta.tsv: total bp + GC desde el
        .fna[.gz] de --genomes-dir-fna (memoizado por path/mtime/size)."""
        d = getattr(self.args, "genomes_dir_fna", None)
        if not d:
            return None
        d = Path(d)
        p = None
        # cualquier fichero del genoma con sufijo FASTA (nombres NCBI
        # _genomic.fna.gz, renombrados .fna/.fasta, con o sin .gz)
        fasta_suf = (".fna", ".fna.gz", ".fasta", ".fasta.gz",
                     ".ffn", ".ffn.gz", ".fnn", ".fnn.gz")
        for cand in sorted(d.glob(glob.escape(genome) + "*")):
            if cand.name.lower().endswith(fasta_suf):
                p = cand
                break
        if not p:
            m04.log.warning(
                f"[{genome}] badge ADN n/d: no hay {genome}*.fna[.gz] en "
                f"{d} (ejemplos ahí: "
                f"{[x.name for x in sorted(d.glob('*'))[:3]]})")
            return None
        st = p.stat()
        key = (str(p), st.st_mtime_ns, st.st_size)
        hit = self._fna_memo.get(key)
        if hit is None:
            hit = _fna_facts(p)
            self._fna_memo[key] = hit
            m04.log.info(f"[{genome}] badge ADN desde {p.name}: "
                         f"{hit[0]} bp, GC {hit[1]}%")
        tot, gc = hit
        if not tot:
            return None
        return {"total_length_bp": tot, "gc_percent": gc}, "computed_fna"

    def payload(self, genome):
        with self._lock:
            if genome in self._payloads:
                return self._payloads[genome]
        gm_path = self._ensure_gene_map(genome)
        gene_map = load_gene_map(gm_path) if gm_path else {}
        meta = load_genome_meta(self.args.genome_meta, genome)
        gm_status = meta_status(self.args)
        if not meta:
            fna = self._fna_meta(genome)
            if fna:
                meta, gm_status = fna
        t0 = time.perf_counter()
        pl = build_payload(genome, self.net_dir, self.rhea_dir, gene_map,
                           meta, self.systematic, self.args.keep_currency,
                           hub_degree=self.args.branch_hub_degree,
                           genome_meta_status=gm_status)
        with self._lock:
            self._payloads[genome] = pl
        st = pl["stats"]
        m04.log.info(f"[{genome}] payload en {time.perf_counter() - t0:.2f}s:"
                     f" {st['reactions_in_genome']} reacciones,"
                     f" {st['branches']} ramas"
                     f" ({st['branches_incomplete']} incompletas)")
        return pl

    def cands(self, genome):
        return load_candidates_list(self.args.candidates, genome)

    def render_page(self):
        # servidor: el HTML va SIN payloads (el JS los pide por genoma a
        # /api/payload). Embeberlos aquí era lo que congelaba la carga
        # inicial minutos con muchos genomas.
        return render_portal({}, {}, server=True, genomes=self.genomes())

    def run08(self, body):
        genome = str(body.get("genome") or "")
        if genome not in self.genomes():
            return False, f"genoma desconocido: {genome}"
        cmd = [sys.executable,
               str(Path(__file__).resolve().parent / "08_fold_candidates.py"),
               "--genome", genome,
               "--net-dir", str(self.net_dir),
               "--rhea-dir", str(self.rhea_dir),
               "--out-dir", str(self.args.candidates)]
        if self.args.gene_map:
            cmd += ["--gene-map", str(self.args.gene_map)]
        if self.args.genome_meta:
            cmd += ["--genome-meta", str(self.args.genome_meta)]
        if self.args.proteins_dir:
            cmd += ["--proteins-dir", str(self.args.proteins_dir)]
        if body.get("orfs"):
            if not self.args.genomes_dir_fna:
                return False, ("ORFs requiere lanzar el servidor con "
                               "--genomes-dir-fna")
            cmd += ["--orfs", "--genomes-dir-fna",
                    str(self.args.genomes_dir_fna)]
        cmd += ["--tool", str(body.get("tool") or "auto")]
        if body.get("max_fold"):
            cmd += ["--max-fold", str(int(body["max_fold"]))]
        env = dict(os.environ)
        prov = str(body.get("esm_provider") or "auto")
        if prov in ("neurosnap", "direct"):
            env["ESM_API_PROVIDER"] = prov
        if body.get("esm_api_url"):
            env["ESM_API_URL"] = str(body["esm_api_url"])
        if body.get("esm_api_key"):
            env["ESM_API_KEY"] = str(body["esm_api_key"])
        if prov == "neurosnap" and not env.get("ESM_API_KEY", "").strip():
            return False, ("Neurosnap requiere ESM_API_KEY: ponla en el "
                           "archivo .env de la raíz del proyecto "
                           "(ESM_API_KEY=...) o en este formulario")
        started = self.job.start(cmd, env, genome)
        if not started:
            return False, "ya hay una ejecución de 08 en marcha"
        return True, None


def make_handler(portal):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):  # silenciar access log
            pass

        def _send(self, code, body, ctype="application/json; charset=utf-8"):
            data = body if isinstance(body, bytes) else body.encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _json(self, obj, code=200):
            self._send(code, json.dumps(obj, ensure_ascii=False))

        def do_GET(self):
            path = urllib.parse.urlparse(self.path).path
            if path in ("/", "/index.html"):
                self._send(200, portal.render_page(), "text/html; charset=utf-8")
                return
            parts = path.strip("/").split("/")
            if len(parts) == 3 and parts[0] == "api":
                kind, g = parts[1], urllib.parse.unquote(parts[2])
                if kind == "payload":
                    try:
                        self._json(portal.payload(g))
                    except Exception as e:  # noqa: BLE001
                        self._json({"error": str(e)}, 500)
                    return
                if kind == "cands":
                    self._json({"genome": g, "candidates": portal.cands(g)})
                    return
            if path == "/api/genomes":
                self._json({"genomes": portal.genomes(),
                            "running": portal.job.running()})
                return
            if path == "/api/status":
                self._json(portal.job.status())
                return
            self._send(404, '{"error": "not found"}')

        def do_POST(self):
            path = urllib.parse.urlparse(self.path).path
            if path == "/api/run08":
                try:
                    n = int(self.headers.get("Content-Length") or 0)
                    body = json.loads(self.rfile.read(n).decode("utf-8") or "{}")
                except json.JSONDecodeError:
                    self._json({"error": "JSON inválido"}, 400)
                    return
                ok, err = portal.run08(body)
                self._json({"ok": ok, **({"error": err} if err else {})},
                           200 if ok else 409)
                return
            if path == "/api/stop08":
                stopped = portal.job.stop()
                self._json({"ok": stopped,
                            **({} if stopped else {"error": "no hay proceso"})})
                return
            self._send(404, '{"error": "not found"}')

    return Handler


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="Portal interactivo de redes metabólicas "
                    "(desplegable de genomas, antes/después, ejecuta 08).")
    ap.add_argument("--net-dir", type=Path, default=None)
    ap.add_argument("--rhea-dir", type=Path, default=None)
    ap.add_argument("--gene-map", type=Path, default=None,
                    help="dir con {genome}_gene_map.tsv|.json (o fichero único); "
                         "si falta el de un genoma se AUTO-GENERA con 09 desde "
                         "data/annotations (ver --annotations-dir)")
    ap.add_argument("--annotations-dir", type=Path, default=None,
                    help="dir de anotación para el auto-genemap del 09 "
                         "(por defecto data/annotations, "
                         "data/annotations_ab_initio y data/proteomes)")
    ap.add_argument("--no-auto-genemap", action="store_true",
                    help="no generar gene maps al vuelo con el 09")
    ap.add_argument("--branch-hub-degree", type=int, default=BRANCH_HUB_DEGREE,
                    help="grado mínimo de metabolito-hub que separa ramas en "
                         "la vista RAMAS (por defecto %d)" % BRANCH_HUB_DEGREE)
    ap.add_argument("--genome-meta", type=Path, default=None)
    ap.add_argument("--cohort", type=Path, default=None,
                    help="trace_compare.json del 05 (brechas sistemáticas)")
    ap.add_argument("--candidates", type=Path, default=None,
                    help="dir de salida de 08 (estado DESPUÉS)")
    ap.add_argument("--proteins-dir", type=Path, default=None,
                    help="dir con {genome}.faa[.gz] (reenviado a 08)")
    ap.add_argument("--genomes-dir-fna", type=Path, default=None,
                    help="dir con {genome}.fna[.gz] (reenviado a 08 para --orfs)")
    ap.add_argument("--keep-currency", action="store_true")
    ap.add_argument("--env-file", type=Path, default=None,
                    help=".env con claves (ESM_API_KEY...); por defecto "
                         ".env en la raíz del proyecto")
    ap.add_argument("--serve", action="store_true",
                    help="modo servidor (ya es el modo por defecto; "
                         "se mantiene por compatibilidad)")
    ap.add_argument("--static", action="store_true",
                    help="solo generar network_portal.html, sin servidor")
    ap.add_argument("--no-browser", action="store_true",
                    help="no abrir el navegador al arrancar el servidor")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    m04.load_env_file(args.env_file)   # ESM_API_KEY etc. (.env no sobreescribe)

    net_dir = args.net_dir or m04.NET_DIR
    args.net_dir = net_dir
    args.rhea_dir = args.rhea_dir or m04.RHEA_DIR
    args.candidates = args.candidates or net_dir / "fold_candidates"

    if not args.static:  # servidor: el modo por defecto
        portal = Portal(args)
        httpd = ThreadingHTTPServer(("127.0.0.1", args.port),
                                    make_handler(portal))
        url = f"http://127.0.0.1:{args.port}"
        print(f"Portal en {url}  (Ctrl+C para salir)")
        print(f"  07_network_view v{SCRIPT_VERSION}: payloads perezosos + "
              f"badge ADN desde --genomes-dir-fna (si no hay genome_meta)")
        print(f"  genomas: {len(portal.genomes())} "
              f"(los payloads se construyen al abrir cada uno)")
        print("  vista RUTAS por defecto (flujo por brechas); RED completa "
              "en su pestaña")
        if os.environ.get("ESM_API_KEY", "").strip():
            print("  ESM_API_KEY: presente (.env/entorno) — el 08 usará "
                  "Neurosnap con auto")
        else:
            print("  ESM_API_KEY: ausente — ponla en .env (raíz del "
                  "proyecto) o en el formulario ⚙")
        print("  lanza el 08 desde la interfaz con el botón ⚙")
        if not args.no_browser:
            threading.Timer(0.8, lambda: webbrowser.open(url)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nAdiós")
        return

    # modo estático: un portal con todos los genomas embebidos
    systematic = load_cohort(args.cohort or net_dir / "trace_compare.json")
    genomes = sorted(p.name[:-len("_network.json")]
                     for p in net_dir.glob("*_network.json"))
    if not genomes:
        print(f"No hay redes en {net_dir}")
        return
    payloads, cands = {}, {}
    for genome in genomes:
        try:
            gm = None
            if args.gene_map:
                gm = (args.gene_map / f"{genome}_gene_map.tsv")
                if not gm.exists():
                    gm = (args.gene_map / f"{genome}_gene_map.json")
                if not gm.exists() and args.gene_map.suffix in (".tsv", ".json"):
                    gm = args.gene_map
            gene_map = load_gene_map(gm) if gm and gm.exists() else {}
            meta = load_genome_meta(args.genome_meta, genome)
            payloads[genome] = build_payload(genome, net_dir, args.rhea_dir,
                                             gene_map, meta, systematic,
                                             args.keep_currency,
                                             hub_degree=args.branch_hub_degree,
                                             genome_meta_status=meta_status(args))
            cands[genome] = {"candidates": load_candidates_list(
                args.candidates, genome)}
            print(f"[{genome}] {payloads[genome]['stats']} "
                  f"(candidatos: {len(cands[genome]['candidates'])})")
        except Exception as e:  # noqa: BLE001 — un genoma no debe tumbar el portal
            print(f"[{genome}] ERROR, se omite: {e}")
    out_dir = args.out_dir or net_dir / "network_views"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "network_portal.html"
    out.write_text(render_portal(payloads, cands, server=False),
                   encoding="utf-8")
    print(f"Portal estático ({len(payloads)} genomas) -> {out}\n"
          f"Para ejecutar el 08 desde la interfaz: "
          f"python {Path(__file__).name}  # abre http://127.0.0.1:{DEFAULT_PORT}")


if __name__ == "__main__":
    main()
