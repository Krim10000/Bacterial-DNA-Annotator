#!/usr/bin/env python3
"""C0 — evaluador retrospectivo de guilt-by-association (GBA) sobre represores con efector.

Pregunta que responde
---------------------
Dado un represor con efector conocido, ¿la ruta metabólica que regula (vecindad /
regulón) rankea al efector verdadero mejor que el azar?

Cadena de datos (todo local, sin GPU):
    TF -> genes regulados (NetworkRegulatorGene)
       -> EC (regulation/ecoli_k12/genes.tsv)
       -> Rhea master (rhea2ec.tsv)
       -> compuestos participantes (rhea-reaction-smiles.tsv + rhea-chebi-smiles.tsv)
       -> expansión a k saltos en el grafo compuesto-reacción

Nulos obligatorios (los cuatro de PLAN_ALOSTERIA §8.1):
    N1 azar global        compuesto al azar del universo de reacciones
    N2 prior de familia   solo "efectores que usan TFs de la misma familia" (LOO)
    N3 prior de frecuencia solo "efectores más frecuentes en general" (LOO)
    N4 regulón aleatorio  conjunto de genes del mismo tamaño, sin la biología

Salidas (idempotentes):
    data/allostery/c0_gba.tsv     una fila por (TF, efector)
    data/allostery/c0_gba.json    métricas + nulos
    data/allostery/c0_gba.md      reporte legible

Uso:
    python scripts/allostery/00_c0_gba.py [--maxdist 3] [--reps 300]
"""
from __future__ import annotations

import argparse
import collections
import csv
import difflib
import json
import random
import re
import statistics
import unicodedata
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem.inchi import MolToInchiKey
except Exception:  # modo degradado: comparación exacta de SMILES
    Chem = None
    MolToInchiKey = None
else:
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")   # InChI falla en R-groups y polímeros: no es error
    except Exception:
        pass

ROOT = Path(__file__).resolve().parents[2]
RD = ROOT / "data" / "regulondb"
RXN = ROOT / "data" / "reactions" / "tsv"
ECO = ROOT / "data" / "regulation" / "ecoli_k12" / "genes.tsv"
OUT = ROOT / "data" / "allostery"

_KEY_CACHE: dict[str, str] = {}


def chem_key(smiles: str) -> str:
    """Clave insensible a tautómeros, protonación y estereoquímica: bloque esqueleto
    del InChIKey. Sin RDKit devuelve el SMILES tal cual (modo degradado)."""
    if smiles in _KEY_CACHE:
        return _KEY_CACHE[smiles]
    key = smiles
    if Chem is not None:
        m = Chem.MolFromSmiles(smiles)
        if m is not None:
            try:
                ik = MolToInchiKey(m)
            except Exception:
                ik = None
            key = ik.split("-")[0] if ik else Chem.MolToSmiles(m)
    _KEY_CACHE[smiles] = key
    return key


HTML = re.compile(r"<[^>]+>")
ROLE_TAIL = re.compile(r"\s+(transcriptional\s+)?(activator|repressor|dual\s+regulator)$", re.I)

# Alias manuales: nombre en RegulonDB -> nombre en ChEBI. Solo lo que el fuzzy no cierra.
ALIAS = {
    "cyclic amp": "camp",
    "3 5 cyclic amp": "camp",
    "ppgpp": "guanosine 3 5 bis diphosphate",
    "alpha d xylopyranose": "alpha d xylose",
    "d trehalose 6 phosphate": "trehalose 6 phosphate",
    "alpha d glucuronate": "alpha d glucuronic acid",
    "s 2 acetolactate": "2 acetolactate",
    "l carnitine": "carnitine",
}


def norm(s: str) -> str:
    s = HTML.sub("", s)
    s = s.replace("&alpha;", "alpha").replace("&beta;", "beta").replace("&amp;", "and")
    s = s.replace("\u03b1", "alpha").replace("\u03b2", "beta").replace("\u00d7", "x")
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower().replace("_", " ")
    s = re.sub(r"[^a-z0-9+ ]", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def clean_effector(tf: str, raw: str) -> str:
    """Quita HTML, prefijo 'TF-' y cola de rol. 'MetJ-<i>S</i>-adenosil...' -> 'S-adenosil...'."""
    s = HTML.sub("", raw).strip()
    m = re.match(r"^([A-Za-z0-9]+)\s*-\s*(.+)$", s)
    if m and m.group(1).lower() == tf.lower():
        s = m.group(2)
    s = ROLE_TAIL.sub("", s).strip()
    return s


def load_chebi():
    name2chebi, chebi2name = {}, {}
    for line in (RXN / "chebiId_name.tsv").read_text(errors="ignore").splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        chebi, name = parts[0].strip(), parts[1].strip()
        chebi2name[chebi] = name
        name2chebi.setdefault(norm(name), chebi)
    chebi2smiles, chebi2key = {}, {}
    for line in (RXN / "rhea-chebi-smiles.tsv").read_text(errors="ignore").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0].startswith("CHEBI:"):
            chebi, smi = parts[0].strip(), parts[1].strip()
            chebi2smiles[chebi] = smi
            chebi2key[chebi] = chem_key(smi)
    return name2chebi, chebi2name, chebi2smiles, chebi2key


def resolve_effector(tf, raw, name2chebi):
    """Devuelve (chebi, nombre_normalizado, tipo_de_match) o (None, motivo)."""
    n = norm(clean_effector(tf, raw))
    if n in ALIAS:
        n = ALIAS[n]
    if n in name2chebi:
        return name2chebi[n], n, "exact"
    cand = difflib.get_close_matches(n, name2chebi.keys(), n=1, cutoff=0.90)
    if cand:
        return name2chebi[cand[0]], cand[0], f"fuzzy:{n}"
    return None, n, "unresolved"


def load_reactions():
    ec2master = collections.defaultdict(set)
    for r in csv.DictReader(open(RXN / "rhea2ec.tsv"), delimiter="\t"):
        if r["ID"].strip():
            ec2master[r["ID"].strip()].add(r["MASTER_ID"].strip())
    master2dirs = {}
    for r in csv.DictReader(open(RXN / "rhea-directions.tsv"), delimiter="\t"):
        master2dirs[r["RHEA_ID_MASTER"]] = [
            r["RHEA_ID_LR"], r["RHEA_ID_RL"], r["RHEA_ID_BI"]
        ]
    dir2tokens = {}
    for line in (RXN / "rhea-reaction-smiles.tsv").read_text(errors="ignore").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            dir2tokens[parts[0].strip()] = [chem_key(t) for t in re.split(r"[.>]", parts[1]) if t]
    master2comp = {}
    for m, dirs in master2dirs.items():
        comp = set()
        for d in dirs:
            comp.update(dir2tokens.get(d, ()))
        if comp:
            master2comp[m] = comp
    return ec2master, master2comp


def load_genes():
    gene2ec, gene2fam = {}, {}
    for r in csv.DictReader(open(ECO), delimiter="\t"):
        g = (r.get("gene_name") or "").strip() or (r.get("locus_tag") or "").strip()
        if not g:
            continue
        gene2ec[g] = [x.strip() for x in re.split(r"[|,;]", r.get("ec") or "") if x.strip()]
        if r.get("fam"):
            gene2fam[g] = r["fam"].strip()
    return gene2ec, gene2fam


def load_tf_effectors():
    """Represores con efector en UNA sola conformación (los 36 decisivos) + la evidencia."""
    roles = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(RD / "NetworkRegulatorGene.tsv"), delimiter="\t"):
        roles[r["2)regulatorName"]][r["6)function"]] += 1

    def role(t):
        c = roles.get(t, collections.Counter())
        pos, neg = c.get("+", 0), c.get("-", 0)
        if pos and neg:
            return "dual"
        if pos:
            return "activador"
        if neg:
            return "represor"
        return "sin_aristas"

    out = []
    for r in csv.DictReader(open(RD / "TFSet.tsv"), delimiter="\t"):
        tf = r["2)name"]
        ea = [x for x in (r["10)tfActiveConformationsEffector"] or "").split("|") if x.strip()]
        ei = [x for x in (r["11)tfInactiveConformationsEffector"] or "").split("|") if x.strip()]
        if not (ea or ei) or role(tf) != "represor":
            continue
        if ea and ei:
            modo, effs = "ambiguo", ea + ei
        elif ei:
            modo, effs = "inductor", ei
        else:
            modo, effs = "correpresor", ea
        out.append({
            "tf": tf, "modo": modo, "conf": r["17)confidenceLevel"],
            "pmid": r["18)tfConformationPMID"].strip(), "efectores": effs,
        })
    return out


def load_regulon():
    reg = collections.defaultdict(set)
    for r in csv.DictReader(open(RD / "NetworkRegulatorGene.tsv"), delimiter="\t"):
        g = (r["5)regulatedGeneName"] or "").strip()
        if g:
            reg[r["2)regulatorName"]].add(g)
    return reg


def build_adj(comp2rxn, rxn2comp):
    """Adyacencia compuesto -> compuestos que comparten reacción (una sola vez)."""
    adj = {}
    for c, rxns in comp2rxn.items():
        nb = set()
        for rx in rxns:
            nb |= rxn2comp.get(rx, set())
        nb.discard(c)
        adj[c] = nb
    return adj


def expand(seeds, adj, maxdist):
    """BFS compuesto -> compuesto, devolviendo {compuesto: distancia}."""
    dist = {c: 0 for c in seeds}
    frontier = set(seeds)
    for d in range(1, maxdist + 1):
        nxt = set()
        for c in frontier:
            for c2 in adj.get(c, ()):
                if c2 not in dist:
                    nxt.add(c2)
        if not nxt:
            break
        for c in nxt:
            dist[c] = d
        frontier = nxt
    return dist


def rank_of(target, dist, support, prior):
    """Orden: distancia, luego soporte en el regulón, luego prior de familia, luego SMILES."""
    ordered = sorted(dist, key=lambda c: (dist[c], -support.get(c, 0),
                                          -prior.get(c, 0), c))
    try:
        return ordered.index(target) + 1, len(ordered)
    except ValueError:
        return None, len(ordered)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxdist", type=int, default=3)
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    name2chebi, chebi2name, chebi2smiles, chebi2key = load_chebi()
    ec2master, master2comp = load_reactions()
    gene2ec, gene2fam = load_genes()
    tfs = load_tf_effectors()
    regulon = load_regulon()

    # grafo compuesto <-> reacción
    rxn2comp_raw = master2comp
    comp2rxn = collections.defaultdict(set)
    for rx, comps in rxn2comp_raw.items():
        for c in comps:
            comp2rxn[c].add(rx)
    universe = set(comp2rxn)
    adj = build_adj(comp2rxn, rxn2comp_raw)
    pool = [g for g in gene2ec if gene2ec[g]]

    # prior de familia y de frecuencia global (LOO): sobre efectores resueltos
    resolved, attempts = [], []
    for t in tfs:
        for eff in t["efectores"]:
            chebi, n, how = resolve_effector(t["tf"], eff, name2chebi)
            attempts.append({"tf": t["tf"], "norm": n, "match": how})
            if chebi and chebi in chebi2key:
                resolved.append({**t, "efector": clean_effector(t["tf"], eff),
                                 "chebi": chebi, "norm": n, "match": how,
                                 "key": chebi2key[chebi]})
    fam_of = {t: gene2fam.get(t) for t in regulon}
    by_fam = collections.defaultdict(collections.Counter)
    by_eff = collections.Counter()
    for r in resolved:
        by_fam[fam_of.get(r["tf"])][r["chebi"]] += 1
        by_eff[r["chebi"]] += 1

    rows, nulls_rand, null_dists = [], [], []
    for r in resolved:
        tf, chebi, key = r["tf"], r["chebi"], r["key"]
        rxn_seeds = set()
        n_genes = 0
        for g in regulon.get(tf, ()):
            ecs = gene2ec.get(g, ())
            if not ecs:
                continue
            n_genes += 1
            for ec in ecs:
                for m in ec2master.get(ec, ()):
                    if m in rxn2comp_raw:
                        rxn_seeds.add(m)
        if not rxn_seeds:
            rows.append({**r, "estado": "sin_ruta", "rank": None, "n_cand": 0})
            continue
        seeds = set()
        for rx in rxn_seeds:
            seeds |= rxn2comp_raw[rx]
        if key not in universe:
            # El SMILES de ChEBI no coincide con el de Rhea (tautómero/protonación):
            # es un artefacto de medida, no biología. Se cuenta aparte.
            rows.append({**r, "estado": "efector_ausente_en_rhea", "rank": None,
                         "n_cand": 0, "dist": None, "n_genes": n_genes,
                         "n_rxn": len(rxn_seeds)})
            continue
        d = expand(seeds, adj, args.maxdist)
        ball = expand({key}, adj, args.maxdist)   # comp -> distancia al efector
        # prior LOO
        fam = fam_of.get(tf)
        prior = collections.Counter()
        for c, k in by_fam.get(fam, {}).items():
            if c != chebi:
                prior[c] = k
        support = collections.Counter()
        for rx in rxn_seeds:
            for c in rxn2comp_raw.get(rx, ()):
                support[c] += 1
        rank, n_cand = rank_of(key, d, support, prior)
        dist_real = d.get(key)
        estado = "ok" if dist_real is not None else "fuera_de_alcance"
        rows.append({**r, "estado": estado,
                     "rank": rank, "n_cand": n_cand, "dist": dist_real,
                     "n_genes": n_genes, "n_rxn": len(rxn_seeds)})

        # N4 comparable: mismo estadístico (distancia) para un regulón aleatorio del
        # mismo tamaño. Por simetría, distancia(regulón falso, efector) = min distancia
        # del efector a los compuestos semilla del regulón falso.
        for _ in range(args.reps):
            pick = random.sample(pool, min(n_genes, len(pool)))
            s0 = set()
            for g in pick:
                for ec in gene2ec[g]:
                    for m in ec2master.get(ec, ()):
                        if m in rxn2comp_raw:
                            s0 |= rxn2comp_raw[m]
            dm = min((ball[c] for c in s0 if c in ball), default=None)
            null_dists.append(dm)
        # N1: el conjunto de candidatos como fracción del universo (sanity de cobertura)
        nulls_rand.append(len(d) / len(universe) if universe else 0.0)

    # N2/N3: prior solo, sobre los candidatos con prior>0
    prior_rows = []
    for r in resolved:
        fam = fam_of.get(r["tf"])
        cand = {c: 1 for c in by_fam.get(fam, {}) if c != r["chebi"]}
        rank, n_cand = rank_of(r["key"], cand, {}, {}) if r["key"] in cand else (None, len(cand))
        cand3 = {c: 1 for c in by_eff if c != r["chebi"]}
        rank3, n3 = rank_of(r["key"], cand3, {}, {}) if r["key"] in cand3 else (None, len(cand3))
        prior_rows.append({"tf": r["tf"], "rank_fam": rank, "n_fam": n_cand,
                           "rank_freq": rank3, "n_freq": n3})

    ok = [r for r in rows if r["estado"] == "ok" and r["rank"]]
    fuera = [r for r in rows if r["estado"] == "fuera_de_alcance"]
    def rec(k):
        return sum(1 for r in ok if r["rank"] <= k) / len(ok) if ok else None

    def frac_dists(ds, thr):
        return sum(1 for x in ds if x is not None and x <= thr) / len(ds) if ds else None

    real_dists = [r["dist"] for r in ok] + [None for r in fuera]
    summary = {
        "n_tf_con_efector": len(tfs),
        "n_resueltos_chebi": len(resolved),
        "n_evaluados": len(ok),
        "n_fuera_de_alcance": len(fuera),
        "n_efector_ausente_en_rhea": sum(1 for r in rows if r["estado"] == "efector_ausente_en_rhea"),
        "ausentes": [f"{r['tf']}:{r['efector']}" for r in rows if r["estado"] == "efector_ausente_en_rhea"],
        "unresolved": sorted({a["norm"] for a in attempts if a["match"] == "unresolved"}),
        "sin_ruta": [r["tf"] for r in rows if r["estado"] == "sin_ruta"],
        "maxdist": args.maxdist,
        "recall": {f"@{k}": rec(k) for k in (1, 3, 5, 10)},
        "mrr": statistics.mean(1 / r["rank"] for r in ok) if ok else None,
        "median_rank": statistics.median(r["rank"] for r in ok) if ok else None,
        "median_n_cand": statistics.median(r["n_cand"] for r in ok) if ok else None,
        "azar_en_el_conjunto@5": (5 / statistics.median([r["n_cand"] for r in ok])) if ok else None,
        "cobertura_conjunto_fraccion_universo": statistics.mean(nulls_rand) if nulls_rand else None,
        "dist_real": {f"<={t}": frac_dists(real_dists, t) for t in (0, 1, 2)},
        "dist_azar_regulon": {f"<={t}": frac_dists(null_dists, t) for t in (0, 1, 2)},
        "null_prior_familia": {
            "n_con_prior": sum(1 for p in prior_rows if p["rank_fam"]),
            "recall": {f"@{k}": (sum(1 for p in prior_rows if p["rank_fam"] and p["rank_fam"] <= k)
                                 / max(1, sum(1 for p in prior_rows if p["rank_fam"]))) for k in (1, 3, 5)},
        },
        "null_prior_frecuencia": {
            "n_con_prior": sum(1 for p in prior_rows if p["rank_freq"]),
            "recall": {f"@{k}": (sum(1 for p in prior_rows if p["rank_freq"] and p["rank_freq"] <= k)
                                 / max(1, sum(1 for p in prior_rows if p["rank_freq"]))) for k in (1, 3, 5)},
        },
        "por_modo": {},
    }
    for modo in ("inductor", "correpresor"):
        sub = [r for r in ok if r["modo"] == modo]
        if sub:
            summary["por_modo"][modo] = {
                "n": len(sub),
                "recall@5": sum(1 for r in sub if r["rank"] <= 5) / len(sub),
                "median_rank": statistics.median(r["rank"] for r in sub),
            }

    with open(OUT / "c0_gba.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "tf", "modo", "conf", "pmid", "efector", "chebi", "key", "match", "estado",
            "rank", "n_cand", "dist", "n_genes", "n_rxn"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in w.fieldnames})
    (OUT / "c0_gba.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    write_md(summary, ok, prior_rows)
    print(json.dumps(summary, indent=1, ensure_ascii=False))


def write_md(summary, ok, prior_rows):
    L = ["# C0 — guilt-by-association, evaluación retrospectiva", ""]
    L.append(f"Represores con efector: {summary['n_tf_con_efector']} · "
             f"resueltos a ChEBI: {summary['n_resueltos_chebi']} · "
             f"evaluados: {summary['n_evaluados']} · "
             f"fuera de alcance: {summary['n_fuera_de_alcance']} · "
             f"SMILES ausente en Rhea (artefacto): {summary['n_efector_ausente_en_rhea']} · "
             f"maxdist={summary['maxdist']}")
    L += ["", "## Ranking dentro del conjunto", "", "| métrica | valor |", "|---|---:|"]
    for k, v in summary["recall"].items():
        L.append(f"| recall{k} | {v:.3f} |" if v is not None else f"| recall{k} | — |")
    L.append(f"| MRR | {summary['mrr']:.3f} |" if summary["mrr"] else "| MRR | — |")
    L.append(f"| rango mediano | {summary['median_rank']} |")
    L.append(f"| candidatos medianos | {summary['median_n_cand']} |")
    a = summary.get("azar_en_el_conjunto@5")
    L.append(f"| azar dentro del conjunto @5 | {a:.5f} |" if a else "| azar dentro del conjunto @5 | — |")

    L += ["", "## Nulo comparable (mismo estadístico: distancia)", "",
          "| distancia | regulón real | regulón aleatorio (mismo tamaño) |", "|---|---:|---:|"]
    for t in ("<=0", "<=1", "<=2"):
        r = summary["dist_real"].get(t)
        n = summary["dist_azar_regulon"].get(t)
        L.append(f"| {t} | {r:.3f} | {n:.3f} |")

    L += ["", "## Priors solos (sin ruta)", "",
          f"- prior de familia: {summary['null_prior_familia']['n_con_prior']} casos con prior; "
          f"recall@1/3/5 = "
          + "/".join(f"{v:.2f}" for v in summary['null_prior_familia']['recall'].values()),
          f"- prior de frecuencia: {summary['null_prior_frecuencia']['n_con_prior']} casos; "
          f"recall@1/3/5 = "
          + "/".join(f"{v:.2f}" for v in summary['null_prior_frecuencia']['recall'].values()),
          ]
    if summary["por_modo"]:
        L += ["", "## Por modo", "", "| modo | n | recall@5 | rango mediano |", "|---|---:|---:|---:|"]
        for m, d in summary["por_modo"].items():
            L.append(f"| {m} | {d['n']} | {d['recall@5']:.3f} | {d['median_rank']} |")

    L += ["", "## Casos", "", "| TF | modo | efector | rank | n_cand | dist |", "|---|---|---|---:|---:|---:|"]
    for r in sorted(ok, key=lambda x: x["rank"]):
        L.append(f"| {r['tf']} | {r['modo']} | {r['efector']} | {r['rank']} | "
                 f"{r['n_cand']} | {r['dist']} |")
    (OUT / "c0_gba.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
