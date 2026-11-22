#!/usr/bin/env python3
"""Construye el grafo Rhea y busca cadenas alternativas para 5 familias piloto.

Entradas (solo lectura)
-----------------------
- data/metabolism/rhea_reactions.pkl      Rhea: 4 variantes (UN/LR/RL/BI) por reacción.
- data/kegg/pathway_reactions.tsv         reacciones KEGG -> Rhea maestro.
- out_fl/refined/refined_gaps.tsv         brechas priorizadas (producers, label).
- data/metabolism/real_gaps.tsv           brechas con soporte en 400 genomas.
- data/metabolism/thermo_results.tsv      ΔG' (opcional; casi vacío).
- data/reactions/tsv/chebiId_name.tsv     nombres ChEBI usados por Rhea.

Salidas (solo `--out-dir`, por defecto data/chains/)
---------------------------------------------------
- families.tsv          configuración de familias y precursores (auditable).
- targets.tsv           todas las brechas con su familia asignada y ranking.
- <family>.tsv          cadenas: source, target, path, length, reversible, evidence.
- summary.json          estadísticas, método, límites y hallazgos.
- validation.json       escrito por validate_chains.py (no por este script).

Política de dirección
---------------------
Rhea documenta que en la reacción maestra (`UN`) *los dos lados son
arbitrarios* (ver data/metabolism/DIRECTION_THERMO.md). Por eso este módulo
construye rutas **topológicas**: cada reacción se puede recorrer en ambos
sentidos y la columna `reversible` declara `unknown` salvo que exista ΔG' en
`thermo_results.tsv`. No se inventa dirección ni termodinámica.

Uso
---
    python scripts/regulation/build_chains.py
    python scripts/regulation/validate_chains.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

SCHEMA_VERSION = "1.0.0"
TOOL = "scripts/regulation/build_chains.py"

# ---------------------------------------------------------------------------
# Configuración de familias
# ---------------------------------------------------------------------------
# `pathways`: mapas KEGG que definen el universo de reacciones de la familia.
# `sources`: precursores externos declarados (ChEBI + nombre Rhea + motivo).
# No se derivan automáticamente: son una decisión de modelado explícita.
FAMILIES: dict[str, dict] = {
    "cisteina_metionina": {
        "label": "Cisteína y metionina",
        "pathways": ["map00270", "map00920"],
        "sources": [
            ("CHEBI:33384", "L-serine", "esqueleto de cisteína"),
            ("CHEBI:58340", "O-acetyl-L-serine", "precursor directo de cisteína"),
            ("CHEBI:57476", "L-homoserine", "precursor de metionina vía aspartato"),
            ("CHEBI:29991", "L-aspartate", "precursor de homoserina/treonina"),
            ("CHEBI:16189", "sulfate", "asimilación de azufre"),
            ("CHEBI:35235", "L-cysteine", "donante de azufre para metionina"),
            ("CHEBI:57844", "L-methionine", "aminoácido central de la familia"),
        ],
    },
    "glutation": {
        "label": "Glutatión",
        "pathways": ["map00480"],
        "sources": [
            ("CHEBI:29985", "L-glutamate", "componente de glutatión"),
            ("CHEBI:35235", "L-cysteine", "componente de glutatión"),
            ("CHEBI:57305", "glycine", "componente de glutatión"),
            ("CHEBI:35491", "L-cystine", "fuente de cisteína exógena"),
        ],
    },
    "folato": {
        "label": "Folato",
        "pathways": ["map00790", "map00670"],
        "sources": [
            ("CHEBI:37565", "GTP", "origen de la pteridina (GTP ciclohidrolasa)"),
            ("CHEBI:17836", "4-aminobenzoate", "PABA, componente del folato"),
            ("CHEBI:29985", "L-glutamate", "cola de glutamato"),
            ("CHEBI:17001", "7,8-dihydroneopterin", "intermedio de pteridina"),
        ],
    },
    "purinas": {
        "label": "Purinas",
        "pathways": ["map00230"],
        "sources": [
            ("CHEBI:58017", "5-phospho-alpha-D-ribose 1-diphosphate", "PRPP de novo"),
            ("CHEBI:58359", "L-glutamine", "donante de nitrógeno"),
            ("CHEBI:57305", "glycine", "C4/C5 del anillo de purina"),
            ("CHEBI:29991", "L-aspartate", "donante de N en adenilosuccinato"),
            ("CHEBI:15740", "formate", "C2/C8 del anillo"),
            ("CHEBI:17596", "inosine", "salvage de nucleósido"),
            ("CHEBI:58115", "GMP", "salvage de nucleótido"),
        ],
    },
    "pared_celular": {
        "label": "Pared celular",
        "pathways": ["map00550", "map00540", "map00541", "map00520"],
        "sources": [
            ("CHEBI:57705", "UDP-N-acetyl-alpha-D-glucosamine", "precursor de peptidoglicano"),
            ("CHEBI:70757", "UDP-N-acetyl-alpha-D-muramate", "precursor de peptidoglicano"),
            ("CHEBI:57972", "L-alanine", "enlace peptídico del tallo"),
            ("CHEBI:57416", "D-alanine", "dipéptido D-Ala-D-Ala"),
            ("CHEBI:29986", "D-glutamate", "segundo aminoácido del tallo"),
            ("CHEBI:57791", "meso-2,6-diaminopimelate", "tercer aminoácido del tallo (Gram-)"),
            ("CHEBI:58885", "UDP-alpha-D-glucose", "precursor de UDP-glucuronato"),
            ("CHEBI:57477", "dTDP-alpha-D-glucose", "precursor de dTDP-ramnosa"),
        ],
    },
}

# Metabolismo "moneda": se permite como origen/fin, pero no como nodo interno.
# Evita atajos artificiales del tipo X -> ATP -> Z entre reacciones no
# relacionadas. La lista es explícita para poder auditarla.
CURRENCY: frozenset[str] = frozenset({
    "CHEBI:15378",  # H+
    "CHEBI:15377",  # H2O
    "CHEBI:43474",  # phosphate
    "CHEBI:33019",  # diphosphate
    "CHEBI:16526",  # CO2
    "CHEBI:28938",  # NH4+
    "CHEBI:16134",  # NH3
    "CHEBI:16240",  # H2O2
    "CHEBI:15379",  # O2
    # nucleótidos / energía
    "CHEBI:30616", "CHEBI:456216", "CHEBI:456215",  # ATP, ADP, AMP
    "CHEBI:37565", "CHEBI:58189", "CHEBI:58115",   # GTP, GDP, GMP
    "CHEBI:46398", "CHEBI:58223", "CHEBI:60377",   # UTP, UDP, CMP
    "CHEBI:37568",                                  # dTTP
    # redox
    "CHEBI:57540", "CHEBI:57945", "CHEBI:58349", "CHEBI:57783",
    "CHEBI:57692", "CHEBI:58307", "CHEBI:17499", "CHEBI:13193",
    # acarreadores
    "CHEBI:57287", "CHEBI:57288", "CHEBI:58343",
})

# Precursores/cofactores centrales que se asumen disponibles en el citoplasma.
# NO son orígenes de cadena: solo pueden satisfacer el papel de co-sustrato.
# Igual que CURRENCY, es una decisión de modelado explícita y auditable.
COFACTOR_POOL: frozenset[str] = frozenset({
    "CHEBI:58702",  # phosphoenolpyruvate
    "CHEBI:16810",  # 2-oxoglutarate
    "CHEBI:15361",  # pyruvate
    "CHEBI:16452",  # oxaloacetate
    "CHEBI:58017",  # PRPP
    "CHEBI:29985",  # L-glutamate
    "CHEBI:58359",  # L-glutamine
    "CHEBI:57305",  # glycine
    "CHEBI:33384",  # L-serine
    "CHEBI:29991",  # L-aspartate
    "CHEBI:57972",  # L-alanine
    "CHEBI:15740",  # formate
    "CHEBI:58228",  # carbamoyl phosphate
})

MAX_LEN = 6
MAX_PATHS_PER_FAMILY = 50
MAX_PATHS_PER_PAIR = 25
MAX_CO_SUBSTRATES = 2
EXPANSION_CAP = 120_000
TARGETS_PER_FAMILY = 8
MIN_FAMILY_SCORE = 2  # 2*productor + consumidor en reacciones KEGG de la familia
THERMO_THRESHOLD = 10.0  # kJ/mol, ver data/metabolism/DIRECTION_THERMO.md


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: Iterable[dict]) -> int:
    n = 0
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})
            n += 1
    return n


def load_chebi_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names[parts[0]] = parts[1].strip()
    return names


def load_reactions(pkl: Path) -> dict[int, dict]:
    import pickle

    with pkl.open("rb") as fh:
        raw = pickle.load(fh)
    reactions = {}
    for rid, r in raw.items():
        if r.get("type") != "UN":  # solo maestros (citados por el resto del repo)
            continue
        if r.get("transport"):
            continue
        reactions[int(rid)] = r
    return reactions


def load_kegg(path: Path) -> tuple[dict[str, set[int]], dict[str, str]]:
    by_family_pathway: dict[str, set[int]] = collections.defaultdict(set)
    names: dict[str, str] = {}
    for row in read_tsv(path):
        pid = row["pathway_id"]
        rhea = row["rhea"].strip()
        if not rhea.startswith("RHEA:"):
            continue
        by_family_pathway[pid].add(int(rhea.split(":")[1]))
        names[pid] = row["pathway_name"]
    return by_family_pathway, names


def load_gaps(refined: Path, real: Path) -> dict[str, dict]:
    """Une refined_gaps + real_gaps por ChEBI."""
    gaps: dict[str, dict] = {}

    def get(chebi: str, name: str) -> dict:
        return gaps.setdefault(chebi, {
            "chebi": chebi, "name": name, "n_rows": 0,
            "support_400": 0, "n_enz": 0,
            "producers": set(), "consumers": set(),
            "pathway_counts": collections.Counter(),
            "in_refined": False, "in_real": False,
        })

    for row in read_tsv(refined):
        g = get(row["chebi"], row["gap"])
        g["in_refined"] = True
        g["n_rows"] += 1
        g["n_enz"] = max(g["n_enz"], int(row.get("n_enz_producers") or 0))
        for p in (row.get("producers") or "").split(";"):
            if p.strip():
                g["producers"].add(p.strip())

    for row in read_tsv(real):
        g = get(row["chebi"], row["gap"])
        g["in_real"] = True
        g["n_rows"] += 1
        g["support_400"] = max(g["support_400"], int(row.get("support_400") or 0))
        for p in (row.get("pathways") or "").split(";"):
            p = p.strip().split()[0] if p.strip() else ""
            if p:
                g["pathway_counts"][p] += 1
    return gaps


def load_thermo(path: Path) -> dict[int, dict]:
    """rid -> {status, dg} usando la fila con |dg| menor para no exagerar."""
    out: dict[int, dict] = {}
    if not path.exists():
        return out
    by_rid: dict[int, list[dict]] = collections.defaultdict(list)
    for row in read_tsv(path):
        rhea = (row.get("reaction") or "").strip()
        if not rhea.startswith("RHEA:"):
            continue
        try:
            dg = float(row["dg"]) if row.get("dg") not in (None, "", "nan") else None
        except ValueError:
            dg = None
        by_rid[int(rhea.split(":")[1])].append({"status": row.get("status", ""), "dg": dg})
    for rid, rows in by_rid.items():
        rows.sort(key=lambda r: abs(r["dg"]) if r["dg"] is not None else 9e9)
        out[rid] = rows[0]
    return out


def step_co_substrates(reactions: dict[int, dict], prev: str, rid: int,
                       nxt: str) -> set[str]:
    """Co-sustratos no-moneda que la reacción exige al recorrer prev -> nxt.

    Si el sentido es inverso (prev es producto), los co-sustratos son los otros
    productos de la reacción. No se comprueba disponibilidad: se declara.
    """
    r = reactions.get(rid)
    if not r:
        return set()
    if prev in r["substrates"] and nxt in r["products"]:
        others = r["substrates"] - {prev}
    elif prev in r["products"] and nxt in r["substrates"]:
        others = r["products"] - {nxt}
    else:
        return set()
    return {m for m in others if m not in CURRENCY}


def missing_co_substrates(reactions: dict[int, dict], chain: list[str],
                          rids: list[int], sources: set[str],
                          pool: frozenset[str] = COFACTOR_POOL) -> list[str]:
    available = set(sources) | set(chain) | set(pool) | set(CURRENCY)
    missing: set[str] = set()
    for i, rid in enumerate(rids):
        for m in step_co_substrates(reactions, chain[i], rid, chain[i + 1]):
            if m not in available:
                missing.add(m)
    return sorted(missing)


# ---------------------------------------------------------------------------
# Grafo
# ---------------------------------------------------------------------------
def build_graph(reactions: dict[int, dict]) -> dict[str, list[tuple[int, str]]]:
    """Grafo bipartito colapsado metabolito-metabolito.

    Cada reacción aporta aristas en ambos sentidos porque los lados del maestro
    son arbitrarios. Las aristas guardan (rhea_id, vecino).
    """
    adj: dict[str, list[tuple[int, str]]] = collections.defaultdict(list)
    for rid, r in reactions.items():
        # ordenar para que el grafo (y por tanto la búsqueda) sea determinista
        for s in sorted(r["substrates"]):
            for p in sorted(r["products"]):
                adj[s].append((rid, p))
                adj[p].append((rid, s))
    for node in adj:
        adj[node].sort()
    return adj


def bfs_from(target: str, adj: dict[str, list[tuple[int, str]]],
             allowed_currency: set[str], limit: int = MAX_LEN) -> dict[str, int]:
    dist = {target: 0}
    q = collections.deque([target])
    while q:
        u = q.popleft()
        if dist[u] >= limit:
            continue
        for _rid, v in adj.get(u, ()):  # type: ignore[arg-type]
            if v in CURRENCY and v not in allowed_currency:
                continue
            if v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def enumerate_paths(adj, source: str, target: str, dist: dict[str, int],
                    reactions: dict[int, dict], available_static: set[str] | None,
                    max_len: int = MAX_LEN, max_paths: int = MAX_PATHS_PER_PAIR,
                    cap: int = EXPANSION_CAP) -> tuple[list[list[tuple[int, str]]], int]:
    """Caminos simples más cortos primero (iterative deepening).

    Un paso solo es válido si todos los co-sustratos no-moneda que exige están
    en `available_static` (fuentes + pool) o ya aparecieron en la cadena. Eso
    evita rutas que dependen de metabolitos exóticos no disponibles.

    Devuelve (caminos, expansiones). Cada camino es [(rhea_id, metabolito), ...].
    """
    if source == target or source not in dist:
        return [], 0
    out: list[list[tuple[int, str]]] = []
    expansions = 0

    for length in range(1, max_len + 1):
        if len(out) >= max_paths:
            break

        def dfs(u: str, path: list[tuple[int, str]], used: set[int],
                visited: set[str], depth: int) -> None:
            nonlocal expansions
            expansions += 1
            if expansions > cap or len(out) >= max_paths:
                return
            if depth == length:
                if u == target:
                    out.append(list(path))
                return
            remaining = length - depth
            for rid, v in adj.get(u, ()):
                if rid in used or v in visited:
                    continue
                if v in CURRENCY and v != target:
                    continue
                if dist.get(v, 10**9) > remaining - 1:
                    continue
                if available_static is not None:
                    req = step_co_substrates(reactions, u, rid, v)
                    if req - available_static - visited:
                        continue
                path.append((rid, v))
                used.add(rid)
                visited.add(v)
                dfs(v, path, used, visited, depth + 1)
                path.pop()
                used.discard(rid)
                visited.discard(v)
                if expansions > cap or len(out) >= max_paths:
                    return

        dfs(source, [], set(), {source}, 0)

    return out, expansions


# ---------------------------------------------------------------------------
# Familias y objetivos
# ---------------------------------------------------------------------------
def assign_families(gaps: dict[str, dict], reactions: dict[int, dict],
                    family_rx: dict[str, set[int]]):
    """Asigna cada brecha a una familia y devuelve (asignación, productores, consumidores).

    Puntaje = 2*(reacciones de la familia que la producen) + (que la consumen).
    La producción pesa más porque la cadena termina en el target. Empates:
    orden de declaración en FAMILIES.
    """
    prod: dict[str, collections.Counter] = {f: collections.Counter() for f in FAMILIES}
    cons: dict[str, collections.Counter] = {f: collections.Counter() for f in FAMILIES}
    for fam, rids in family_rx.items():
        for rid in rids:
            r = reactions.get(rid)
            if not r:
                continue
            for m in r["products"]:
                prod[fam][m] += 1
            for m in r["substrates"]:
                cons[fam][m] += 1

    order = list(FAMILIES)
    assigned: dict[str, str] = {}
    for chebi in gaps:
        scores = {
            fam: (2 * prod[fam][chebi] + cons[fam][chebi], prod[fam][chebi], cons[fam][chebi])
            for fam in FAMILIES
        }
        best = max(scores.items(), key=lambda kv: (kv[1][0], kv[1][1], -order.index(kv[0])))
        if best[1][0] == 0:
            continue
        assigned[chebi] = best[0]
    return assigned, prod, cons


def reversible_for_step(rid: int, thermo: dict[int, dict]) -> str:
    row = thermo.get(rid)
    if not row:
        return "unknown"
    status, dg = row["status"], row["dg"]
    if status == "thermo_reversible":
        return "reversible"
    if status in {"thermo_forward", "thermo_reverse"} and dg is not None and abs(dg) >= THERMO_THRESHOLD:
        return "direction_specific"
    return "unknown"


def classify_path(steps: list[tuple[int, str]], thermo: dict[int, dict]) -> tuple[str, str]:
    kinds = [reversible_for_step(rid, thermo) for rid, _ in steps]
    n_thermo = sum(1 for rid, _ in steps if rid in thermo)
    if any(k == "direction_specific" for k in kinds):
        return "direction_specific", f"{n_thermo}/{len(steps)}"
    if kinds and all(k == "reversible" for k in kinds):
        return "reversible", f"{n_thermo}/{len(steps)}"
    if n_thermo:
        return "mixed", f"{n_thermo}/{len(steps)}"
    return "unknown", f"0/{len(steps)}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("."), help="raíz del proyecto")
    ap.add_argument("--out-dir", type=Path, default=None, help="por defecto <root>/data/chains")
    ap.add_argument("--include-unmapped-rhea", action="store_true",
                    help="usar todo Rhea (por defecto solo reacciones mapeadas a KEGG)")
    ap.add_argument("--strict-cofactors", action="store_true",
                    help="exigir que todo co-sustrato esté en fuentes/pool/cadena "
                         "(por defecto se permite hasta MAX_CO_SUBSTRATES externos)")
    args = ap.parse_args()

    root = args.root.resolve()
    out_dir = (args.out_dir or root / "data" / "chains").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs = {
        "rhea_reactions": root / "data/metabolism/rhea_reactions.pkl",
        "kegg_pathway_reactions": root / "data/kegg/pathway_reactions.tsv",
        "refined_gaps": root / "out_fl/refined/refined_gaps.tsv",
        "real_gaps": root / "data/metabolism/real_gaps.tsv",
        "thermo_results": root / "data/metabolism/thermo_results.tsv",
        "chebi_names": root / "data/reactions/tsv/chebiId_name.tsv",
    }
    for name, path in inputs.items():
        if not path.exists() and name != "thermo_results":
            raise SystemExit(f"falta entrada {name}: {path}")

    names = load_chebi_names(inputs["chebi_names"])
    reactions = load_reactions(inputs["rhea_reactions"])
    kegg_by_pathway, kegg_names = load_kegg(inputs["kegg_pathway_reactions"])
    gaps = load_gaps(inputs["refined_gaps"], inputs["real_gaps"])
    thermo = load_thermo(inputs["thermo_results"])

    family_rx = {
        fam: {rid for pid in cfg["pathways"] for rid in kegg_by_pathway.get(pid, set())}
        for fam, cfg in FAMILIES.items()
    }
    kegg_reactions = {rid for rids in kegg_by_pathway.values() for rid in rids}
    if args.include_unmapped_rhea:
        graph_reactions = reactions
    else:
        graph_reactions = {rid: r for rid, r in reactions.items() if rid in kegg_reactions}

    nodes = set()
    for r in reactions.values():
        nodes |= r["substrates"] | r["products"]

    adj = build_graph(graph_reactions)
    assigned, prod_counter, cons_counter = assign_families(gaps, reactions, family_rx)

    # --- configuración auditable -------------------------------------------
    conf_rows = []
    for fam, cfg in FAMILIES.items():
        for pid in cfg["pathways"]:
            conf_rows.append({
                "family": fam, "pathway_id": pid,
                "pathway_name": kegg_names.get(pid, ""),
                "n_family_reactions": len(kegg_by_pathway.get(pid, ())),
            })
    write_tsv(out_dir / "families.tsv",
              ["family", "pathway_id", "pathway_name", "n_family_reactions"], conf_rows)

    sources_rows = []
    for fam, cfg in FAMILIES.items():
        for chebi, name, reason in cfg["sources"]:
            sources_rows.append({
                "family": fam, "source": name, "source_chebi": chebi,
                "in_rhea_pickle": chebi in nodes, "rationale": reason,
            })
    write_tsv(out_dir / "sources.tsv",
              ["family", "source", "source_chebi", "in_rhea_pickle", "rationale"], sources_rows)
    missing_sources = [r for r in sources_rows if not r["in_rhea_pickle"]]

    # --- targets ------------------------------------------------------------
    target_rows = []
    per_family_targets: dict[str, list[dict]] = collections.defaultdict(list)
    for chebi, g in gaps.items():
        fam = assigned.get(chebi)
        row = {
            "family": fam or "",
            "target": g["name"],
            "target_chebi": chebi,
            "support_400": g["support_400"],
            "n_enz_producers": g["n_enz"],
            "n_gap_rows": g["n_rows"],
            "in_refined_gaps": g["in_refined"],
            "in_real_gaps": g["in_real"],
            "family_producer_reactions": prod_counter[fam][chebi] if fam else 0,
            "family_consumer_reactions": cons_counter[fam][chebi] if fam else 0,
            "in_rhea_pickle": chebi in nodes,
            "chosen": False,
        }
        target_rows.append(row)
        if fam and chebi in nodes:
            per_family_targets[fam].append(row)

    chosen: dict[str, list[dict]] = {}
    for fam, rows in per_family_targets.items():
        rows = [r for r in rows
                if 2 * r["family_producer_reactions"] + r["family_consumer_reactions"] >= MIN_FAMILY_SCORE]
        per_family_targets[fam] = rows
        rows.sort(key=lambda r: (-r["support_400"], -r["n_enz_producers"], r["target"]))
        chosen[fam] = rows[:TARGETS_PER_FAMILY]
        for r in chosen[fam]:
            r["chosen"] = True
    target_rows.sort(key=lambda r: (r["family"], not r["chosen"], -r["support_400"], r["target"]))
    write_tsv(out_dir / "targets.tsv",
              ["family", "target", "target_chebi", "support_400", "n_enz_producers",
               "n_gap_rows", "in_refined_gaps", "in_real_gaps",
               "family_producer_reactions", "family_consumer_reactions",
               "in_rhea_pickle", "chosen"], target_rows)

    # --- búsqueda de cadenas ------------------------------------------------
    out_header = [
        "source", "target", "path", "length", "reversible", "evidence",
        "family", "source_chebi", "target_chebi", "path_metabolites", "path_chebi",
        "kegg_family_steps", "thermo_steps", "co_substrates", "n_co_substrates",
        "final_step_is_known_producer", "alternative_route", "target_rank",
        "target_support_400", "score",
    ]
    summary_families: dict[str, dict] = {}
    findings: list[str] = []
    total_paths = 0

    for fam, cfg in FAMILIES.items():
        rows = chosen.get(fam, [])
        candidates: list[dict] = []
        pair_stats = {"pairs": 0, "with_paths": 0, "expansions": 0}
        for rank, tgt in enumerate(rows, start=1):
            tchebi = tgt["target_chebi"]
            known_producers = {p.split(":")[1] for p in gaps[tchebi]["producers"] if ":" in p}
            allowed_currency = {s[0] for s in cfg["sources"]}
            available_static = allowed_currency | set(COFACTOR_POOL) | set(CURRENCY)
            dist = bfs_from(tchebi, adj, allowed_currency)
            for schebi, sname, _reason in cfg["sources"]:
                if schebi == tchebi:
                    continue
                pair_stats["pairs"] += 1
                paths, expansions = enumerate_paths(
                    adj, schebi, tchebi, dist, reactions,
                    available_static if args.strict_cofactors else None)
                pair_stats["expansions"] += expansions
                if paths:
                    pair_stats["with_paths"] += 1
                for steps in paths:
                    rids = [rid for rid, _ in steps]
                    chain = [schebi] + [m for _, m in steps]
                    fam_steps = sum(1 for rid in rids if rid in family_rx[fam])
                    final_rid = str(rids[-1])
                    is_known = final_rid in known_producers
                    rev, thermo_steps = classify_path(steps, thermo)
                    missing = missing_co_substrates(
                        reactions, chain, rids, {s[0] for s in cfg["sources"]})
                    frac = fam_steps / len(rids)
                    score = round((1.0 / len(rids)) * (0.5 + 0.5 * frac), 4)
                    candidates.append({
                        "source": sname,
                        "target": tgt["target"],
                        "path": ";".join(f"RHEA:{rid}" for rid in rids),
                        "length": len(rids),
                        "reversible": rev,
                        "evidence": (f"rhea_master_verified={len(rids)}/{len(rids)}; "
                                     f"kegg_family_steps={fam_steps}/{len(rids)}; "
                                     f"thermo={thermo_steps}; "
                                     f"co_substrates={len(missing)}; "
                                     f"known_producers={len(known_producers)}"),
                        "family": fam,
                        "source_chebi": schebi,
                        "target_chebi": tchebi,
                        "path_metabolites": " -> ".join(
                            names.get(m, m) for m in chain),
                        "path_chebi": ";".join(chain),
                        "kegg_family_steps": fam_steps,
                        "thermo_steps": thermo_steps,
                        "co_substrates": ";".join(names.get(m, m) for m in missing),
                        "n_co_substrates": len(missing),
                        "final_step_is_known_producer": "yes" if is_known else "no",
                        "alternative_route": "no" if is_known else "yes",
                        "target_rank": rank,
                        "target_support_400": tgt["support_400"],
                        "score": score,
                        "_rids": rids,
                        "_family_frac": frac,
                    })

        candidates = [c for c in candidates
                      if c["n_co_substrates"] <= MAX_CO_SUBSTRATES]
        # dedupe por (source, target, reacciones)
        uniq: dict[tuple, dict] = {}
        for c in candidates:
            key = (c["source_chebi"], c["target_chebi"], tuple(c["_rids"]))
            uniq.setdefault(key, c)
        candidates = list(uniq.values())

        def rank_key(c: dict) -> tuple:
            # preferimos rutas que NO terminen en un productor conocido (alternativas),
            # luego cortas, luego con más pasos dentro de la familia.
            return (c["n_co_substrates"],
                    0 if c["alternative_route"] == "yes" else 1,
                    c["length"], -c["kegg_family_steps"], c["source"])

        # Reparto justo: cada target recibe hasta `per_target` rutas, priorizando
        # fuentes distintas para que las alternativas no sean variantes triviales.
        groups: dict[int, list[dict]] = collections.defaultdict(list)
        for c in candidates:
            groups[c["target_rank"]].append(c)
        per_target = max(1, MAX_PATHS_PER_FAMILY // max(1, len(rows)))
        selected: list[dict] = []
        for target_rank in sorted(groups):
            ranked = sorted(groups[target_rank], key=rank_key)
            picked: list[dict] = []
            seen_sources: set[str] = set()
            for c in ranked:
                if c["source_chebi"] not in seen_sources:
                    picked.append(c)
                    seen_sources.add(c["source_chebi"])
                if len(picked) >= per_target:
                    break
            for c in ranked:  # rellenar si no hay fuentes suficientes
                if len(picked) >= per_target:
                    break
                if c not in picked:
                    picked.append(c)
            selected.extend(picked)

        used = {(c["source_chebi"], c["target_chebi"], tuple(c["_rids"])) for c in selected}
        leftovers = sorted((c for c in candidates
                            if (c["source_chebi"], c["target_chebi"], tuple(c["_rids"])) not in used),
                           key=lambda c: (c["length"], rank_key(c)))
        selected.extend(leftovers[: max(0, MAX_PATHS_PER_FAMILY - len(selected))])
        selected = selected[:MAX_PATHS_PER_FAMILY]
        selected.sort(key=lambda c: (c["target_rank"], c["n_co_substrates"],
                                     c["length"], rank_key(c)[1],
                                     -c["kegg_family_steps"], c["source"]))

        fam_path = out_dir / f"{fam}.tsv"
        n = write_tsv(fam_path, out_header, selected)
        total_paths += n

        lengths = collections.Counter(c["length"] for c in selected)
        revs = collections.Counter(c["reversible"] for c in selected)
        alts = sum(1 for c in selected if c["alternative_route"] == "yes")
        clean = sum(1 for c in selected if c["n_co_substrates"] == 0)
        with_paths = {c["target_chebi"] for c in selected}
        summary_families[fam] = {
            "label": cfg["label"],
            "pathways": {pid: kegg_names.get(pid, "") for pid in cfg["pathways"]},
            "n_family_reactions": len(family_rx[fam]),
            "n_sources": len(cfg["sources"]),
            "n_targets": len(rows),
            "targets": [{"name": t["target"], "chebi": t["target_chebi"],
                         "support_400": t["support_400"], "rank": i}
                        for i, t in enumerate(rows, start=1)],
            "n_candidate_paths": len(candidates),
            "n_paths": n,
            "paths_by_length": {str(k): v for k, v in sorted(lengths.items())},
            "reversible": dict(revs),
            "alternative_routes": alts,
            "paths_without_extra_co_substrates": clean,
            "targets_with_path": len(with_paths),
            "targets_without_path": sorted(
                t["target"] for t in rows if t["target_chebi"] not in with_paths),
            "pairs_explored": pair_stats["pairs"],
            "pairs_with_path": pair_stats["with_paths"],
            "search_expansions": pair_stats["expansions"],
            "file": fam_path.name,
        }

    # --- hallazgos (con números reales) ------------------------------------
    all_rows = []
    for fam in FAMILIES:
        all_rows.extend(read_tsv(out_dir / f"{fam}.tsv"))
    n_alt = sum(1 for r in all_rows if r["alternative_route"] == "yes")
    n_known = len(all_rows) - n_alt
    short = sum(1 for r in all_rows if int(r["length"]) <= 3)
    clean = sum(1 for r in all_rows if r["n_co_substrates"] == "0")
    fam_best = sorted(summary_families.items(), key=lambda kv: -kv[1]["n_paths"])
    unreachable = {f: v["targets_without_path"] for f, v in summary_families.items() if v["targets_without_path"]}

    clean_by_family = {f: v["paths_without_extra_co_substrates"] for f, v in summary_families.items()}
    fam_clean = sorted(clean_by_family.items(), key=lambda kv: -kv[1])
    with_thermo = sum(1 for r in all_rows if r["reversible"] != "unknown")
    target_names = {t["name"] for v in summary_families.values() for t in v["targets"]}
    coupled: list[dict] = []
    for r in all_rows:
        req = {x for x in (r["co_substrates"] or "").split(";") if x}
        inter = sorted(req & target_names)
        if inter:
            coupled.append({"family": r["family"], "target": r["target"],
                            "requires": inter, "source": r["source"]})

    findings.append(
        f"{len(all_rows)} cadenas topológicas en {len(FAMILIES)} familias: "
        f"{short} con longitud <=3, {clean} sin co-sustratos externos y {n_alt} "
        f"evitan terminar en un productor conocido del gap (ruta alternativa real, "
        f"no el paso ausente)."
    )
    findings.append(
        f"Rutas más limpias: {fam_clean[0][0]}={fam_clean[0][1]}/50, "
        f"{fam_clean[1][0]}={fam_clean[1][1]}/50 y {fam_clean[2][0]}={fam_clean[2][1]}/50 "
        f"cadenas sin co-sustratos externos. La más dependiente es "
        f"{fam_clean[-1][0]}={fam_clean[-1][1]}/50 (SAM/THF fuera de la cadena); "
        "esos requerimientos se declaran en `co_substrates`."
    )
    findings.append(
        f"La dirección es el límite principal: {len(all_rows) - with_thermo}/{len(all_rows)} "
        f"cadenas no tienen ningún ΔG' (reversible=unknown) y solo {with_thermo} tienen "
        f"termodinámica asignada; Rhea maestro no define sentido, así que son rutas "
        f"topológicas, no flujos."
    )
    if coupled:
        findings.append(
            f"{len(coupled)} cadenas exigen como co-sustrato otra brecha del mismo piloto "
            f"({len({c['target'] for c in coupled})} targets afectados, p. ej. "
            f"{coupled[0]['target']} requiere {', '.join(coupled[0]['requires'])}): "
            "no cierran la brecha por sí solas."
        )

    summary = {
        "schema_version": SCHEMA_VERSION,
        "tool": TOOL,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": {
            name: {"path": str(path.relative_to(root)), "sha256": sha256(path),
                   "bytes": path.stat().st_size}
            for name, path in inputs.items() if path.exists()
        },
        "rhea": {
            "n_reactions_master_nontransport": len(reactions),
            "n_reactions_in_graph": len(graph_reactions),
            "universe": "all_rhea" if args.include_unmapped_rhea else "kegg_mapped_only",
            "n_metabolites_graph": len({m for r in graph_reactions.values()
                                        for m in r["substrates"] | r["products"]}),
            "n_metabolites_all_rhea": len(nodes),
            "n_edges_metabolite_pairs": sum(len(v) for v in adj.values()),
            "currency_excluded": sorted(CURRENCY),
        },
        "method": {
            "graph": "bipartito colapsado metabolito-reacción, ambos sentidos",
            "reaction_universe": "reacciones Rhea mapeadas a alguna ruta KEGG del TSV"
                                 if not args.include_unmapped_rhea else "todo Rhea",
            "direction_policy": "Rhea UN no define sentido; reversible=unknown salvo ΔG'",
            "max_length": MAX_LEN,
            "max_paths_per_family": MAX_PATHS_PER_FAMILY,
            "max_paths_per_pair": MAX_PATHS_PER_PAIR,
            "targets_per_family": TARGETS_PER_FAMILY,
            "min_family_score": MIN_FAMILY_SCORE,
            "target_selection": "asignación por rol en reacciones KEGG de la familia "
                                "(2*productor + consumidor); top por support_400",
            "thermo_threshold_kj": THERMO_THRESHOLD,
            "max_external_co_substrates": (0 if args.strict_cofactors else MAX_CO_SUBSTRATES),
            "strict_cofactors": args.strict_cofactors,
            "score": "1/length * (0.5 + 0.5*fracción_pasos_de_la_familia); no es confianza",
        },
        "families": summary_families,
        "totals": {
            "families": len(FAMILIES),
            "targets": sum(v["n_targets"] for v in summary_families.values()),
            "paths": total_paths,
            "alternative_routes": n_alt,
            "paths_with_known_final_producer": n_known,
            "paths_without_extra_co_substrates": clean,
            "targets_without_path": sum(len(v) for v in unreachable.values()),
        },
        "findings": findings,
        "coupled_gaps": coupled[:30],
        "limitations": [
            "Rhea maestro (UN): los lados son arbitrarios; no se declara dirección fisiológica.",
            "Sin transporte: no se modela captación ni exportación (fuera de alcance).",
            "Sin estequiometría: no se balancea la reacción; los co-sustratos no-moneda "
            "que cada paso exige se declaran en `co_substrates` (máximo "
            f"{MAX_CO_SUBSTRATES} externos por cadena).",
            "Los pasos se citan solo si existen en rhea_reactions.pkl (tipo UN, no transporte).",
            f"Caminos <={MAX_LEN} pasos: rutas de novo largas quedan fuera o aparecen "
            "como variantes parciales.",
            "El universo de reacciones se limita a Rhea mapeado a KEGG "
            "(4215 de 17177 reacciones maestras no-transporte); reacciones huérfanas "
            "de Rhea no se exploran.",
            "Cobertura de ΔG' mínima (thermo_results.tsv cubre pocas reacciones).",
            "La familia se asigna por pertenencia KEGG: incluye vecinos de mapa "
            "(p. ej. ornithine en map00480, molybdopterin en map00790).",
        ],
        "warnings": {
            "sources_missing_in_pickle": [r["source_chebi"] for r in missing_sources],
            "targets_without_path": unreachable,
        },
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[chains] familias={len(FAMILIES)} targets="
          f"{sum(v['n_targets'] for v in summary_families.values())} paths={total_paths} "
          f"alternativas={n_alt} sin_ruta={sum(len(v) for v in unreachable.values())}")
    for fam, v in summary_families.items():
        print(f"  {fam:22s} paths={v['n_paths']:2d} target={v['targets_with_path']}/{v['n_targets']} "
              f"lens={v['paths_by_length']} rev={v['reversible']}")
    if missing_sources:
        print(f"  [warn] fuentes ausentes en Rhea: {[r['source_chebi'] for r in missing_sources]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
