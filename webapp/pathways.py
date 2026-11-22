"""Vista acotada de rutas metabólicas para la GUI.

No dibuja cientos de nodos: muestra una ruta por vez, con sus pasos en orden.
Cada paso sabe si la enzima está presente en el genoma (verde), falta (rojo) o
no hay datos (gris), y guarda la reacción y el locus del candidato.
"""
from __future__ import annotations

import csv
import json
import pickle
from functools import lru_cache
from pathlib import Path

# Rutas curadas: nombre -> pasos (enzima, EC, etiqueta de reacción).
PATHWAYS: dict[str, list[dict]] = {
    "Riboflavina": [
        {"enzyme": "RibA", "ec": "3.5.4.25", "label": "GTP → pirimidina"},
        {"enzyme": "RibB", "ec": "4.1.99.12", "label": "pirimidina → DHBP"},
        {"enzyme": "RibD", "ec": "1.1.1.193", "label": "reducción del ribosilo"},
        {"enzyme": "RibE", "ec": "2.5.1.9", "label": "riboflavina sintasa"},
        {"enzyme": "RibF", "ec": "2.7.1.26", "label": "flavina quinasa"},
    ],
    "Folato": [
        {"enzyme": "FolE", "ec": "3.5.4.16", "label": "GTP → neopterina"},
        {"enzyme": "FolB", "ec": "4.1.2.25", "label": "neopterina → pterina"},
        {"enzyme": "FolK", "ec": "2.7.6.3", "label": "fosforilación"},
        {"enzyme": "FolP", "ec": "2.5.1.15", "label": "DHPS"},
        {"enzyme": "FolC", "ec": "6.3.2.17", "label": "DHFS"},
    ],
    "Hemo": [
        {"enzyme": "HemB", "ec": "4.2.1.24", "label": "ALA → PBG"},
        {"enzyme": "HemC", "ec": "2.5.1.61", "label": "PBG → HMB"},
        {"enzyme": "HemD", "ec": "4.2.1.75", "label": "HMB → uroporfirinógeno III"},
        {"enzyme": "HemE", "ec": "4.1.1.37", "label": "URO-D"},
        {"enzyme": "HemF", "ec": "1.3.3.3", "label": "coproporfirinógeno oxidasa"},
        {"enzyme": "HemH", "ec": "4.99.1.1", "label": "ferroquelatasa"},
    ],
    "Cofactor de molibdeno": [
        {"enzyme": "MoaA", "ec": "4.1.99.18", "label": "GTP → cPMP"},
        {"enzyme": "MoaC", "ec": "4.6.1.10", "label": "cPMP sintasa"},
        {"enzyme": "MoaD/E", "ec": "2.7.7.80", "label": "molibdopterina sintasa"},
        {"enzyme": "MoeA/MobB", "ec": "2.7.7.75", "label": "activación del cofactor"},
    ],
    "MTA y poliaminas": [
        {"enzyme": "MtnN", "ec": "3.2.2.9", "label": "MTA nucleosidasa"},
        {"enzyme": "MtnA", "ec": "3.1.1.77", "label": "MTA → metiltioribosa"},
        {"enzyme": "MTA fosforilasa", "ec": "2.4.2.28", "label": "MTA → adenina"},
    ],
    "Ácidos grasos": [
        {"enzyme": "FabD", "ec": "2.3.1.39", "label": "malonil-CoA → ACP"},
        {"enzyme": "FabH", "ec": "2.3.1.180", "label": "inicio de cadena"},
        {"enzyme": "FabG", "ec": "1.1.1.100", "label": "cetoacil reductasa"},
        {"enzyme": "FabZ", "ec": "4.2.1.59", "label": "deshidratasa"},
        {"enzyme": "FabI", "ec": "1.3.1.9", "label": "enoyl-ACP reductasa"},
    ],
    "Glucólisis": [
        {"enzyme": "HK", "ec": "2.7.1.1", "label": "glucosa → G6P"},
        {"enzyme": "PGI", "ec": "5.3.1.9", "label": "G6P → F6P"},
        {"enzyme": "PFK", "ec": "2.7.1.11", "label": "F6P → F1,6BP"},
        {"enzyme": "FBA", "ec": "4.1.2.13", "label": "F1,6BP → DHAP + GAP"},
        {"enzyme": "TPI", "ec": "5.3.1.1", "label": "DHAP → GAP"},
        {"enzyme": "GAPDH", "ec": "1.2.1.12", "label": "GAP → 1,3BPG"},
        {"enzyme": "PGK", "ec": "2.7.2.3", "label": "1,3BPG → 3PG"},
        {"enzyme": "ENO", "ec": "4.2.1.11", "label": "2PG → PEP"},
        {"enzyme": "PK", "ec": "2.7.1.40", "label": "PEP → piruvato"},
    ],
    "TCA": [
        {"enzyme": "CS", "ec": "2.3.3.1", "label": "acetil-CoA + OAA → citrato"},
        {"enzyme": "ACON", "ec": "4.2.1.3", "label": "citrato → isocitrato"},
        {"enzyme": "IDH", "ec": "1.1.1.42", "label": "isocitrato → 2-OG"},
        {"enzyme": "KGDH", "ec": "1.2.4.2", "label": "2-OG → succinil-CoA"},
        {"enzyme": "SCS", "ec": "6.2.1.5", "label": "succinil-CoA → succinato"},
        {"enzyme": "SDH", "ec": "1.3.5.1", "label": "succinato → fumarato"},
        {"enzyme": "FH", "ec": "4.2.1.2", "label": "fumarato → malato"},
        {"enzyme": "MDH", "ec": "1.1.1.37", "label": "malato → OAA"},
    ],
    "Carotenoides": [
        {"enzyme": "PSY", "ec": "2.5.1.32", "label": "GGPP → fitoeno"},
        {"enzyme": "PDS", "ec": "1.3.5.5", "label": "fitoeno → fitoflueno"},
        {"enzyme": "ZDS", "ec": "1.3.5.6", "label": "ζ-caroteno → licopeno"},
        {"enzyme": "LCY", "ec": "5.5.1.19", "label": "licopeno → β-caroteno"},
    ],
    "Salvamento de pirimidinas": [
        {"enzyme": "Tdk", "ec": "2.7.1.21", "label": "timidina → dTMP"},
        {"enzyme": "DeoA", "ec": "2.4.2.4", "label": "timidina → timina"},
        {"enzyme": "CDA", "ec": "3.5.4.5", "label": "citidina → uridina"},
        {"enzyme": "ThyA", "ec": "2.1.1.45", "label": "dUMP → dTMP"},
    ],
}


DOMAIN_EC = {
    "RibD_C": "1.1.1.193", "GTP_cyclohydro2": "3.5.4.25", "GTP_cyclohydroI": "3.5.4.16",
    "Glyoxalase": "4.4.1.5", "Glyoxalase_2": "4.4.1.5", "Glyoxalase_3": "4.4.1.5",
    "Glyoxalase_4": "4.4.1.5", "DHFR_1": "1.5.1.3", "URO-D": "4.1.1.37",
    "Coprogen_oxidas": "1.3.3.3", "FolB": "4.1.2.25", "FolK": "2.7.6.3",
    "FolP": "2.5.1.15", "FolE": "3.5.4.16", "ThiE": "2.5.1.3", "TPK": "2.7.6.2",
    "Transketolase": "2.2.1.1", "MobB": "2.7.7.75",
    "MoaC": "4.6.1.10", "MoaA": "4.1.99.18",
}
EC_DOMAINS: dict[str, set[str]] = {}
for _d, _e in DOMAIN_EC.items():
    EC_DOMAINS.setdefault(_e, set()).add(_d)


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))


@lru_cache(maxsize=4)
def _chebi_names(path: str) -> dict[str, str]:
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            out[p[0].strip()] = p[1].strip()
    return out


@lru_cache(maxsize=2)
def _reactions(cache: str) -> dict[str, dict]:
    return pickle.load(open(cache, "rb"))


@lru_cache(maxsize=2)
def _ec_to_rhea(path: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for r in _read_tsv(Path(path)):
        if not r.get("ID"):
            continue
        out.setdefault(r["ID"], set()).add(f"RHEA:{r['RHEA_ID']}")
    return out


@lru_cache(maxsize=2)
def _rhea_acc(path: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[3]:
            out.setdefault(f"RHEA:{p[0]}", set()).add(p[3])
    return out


@lru_cache(maxsize=8)
def _candidates(path: str) -> list[dict]:
    try:
        return json.load(open(path)).get("candidates", [])
    except Exception:
        return []


def pathway_view(genome: str, pathway: str, data_root: Path, out_root: Path) -> list[dict]:
    """Devuelve los pasos con estado, reacción y locus del candidato."""
    net_path = data_root / "metabolism" / f"{genome}_network.json"
    net = {n["id"] for n in json.load(open(net_path)).get("nodes", [])
           if n.get("type") == "reaction"} if net_path.exists() else set()

    refined = [r for r in _read_tsv(out_root / "refined" / "refined_gaps.tsv")
               if r["genome"] == genome]
    gap_producers = set()
    for r in refined:
        gap_producers |= set(r.get("producers", "").split(";"))

    closures = [r for r in _read_tsv(out_root / "closures_labeled.tsv")
                if r["genome"] == genome]
    r2a = _rhea_acc(str(data_root / "metabolism" / "rhea2uniprot_sprot.tsv")) \
        if (data_root / "metabolism" / "rhea2uniprot_sprot.tsv").exists() else {}
    cands = _candidates(str(out_root / f"{genome}_ghost_candidates.json"))
    by_acc: dict[str, dict] = {}
    for c in cands:
        acc = c.get("hit")
        if not acc or c.get("self_hit"):
            continue
        if acc not in by_acc or c.get("rank", 0) > by_acc[acc].get("rank", 0):
            by_acc[acc] = c

    # dominio Pfam (panel) -> locus, vía cluster y miembros
    mem = {}
    mem_path = out_root / "dedup" / "dedup_members.tsv"
    if mem_path.exists():
        for r in _read_tsv(mem_path):
            mem.setdefault(r["cluster"], []).append((r["genome"], r["locus"]))
    domain_loci: dict[str, dict[str, str]] = {}
    panel = out_root / "dedup" / "panel_domtbl.txt"
    if panel.exists():
        for line in open(panel):
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 23:
                continue
            seq, dom = p[0], p[3]
            for g, loc in mem.get(seq, []):
                domain_loci.setdefault(g, {}).setdefault(dom, loc)

    rx_path = None
    for cand in (data_root / "metabolism" / "rhea_reactions.pkl",
                 out_root.parent / "cache" / "rhea_reactions.pkl",
                 Path("cache/rhea_reactions.pkl")):
        if cand.exists():
            rx_path = cand
            break
    rhea = _reactions(str(rx_path)) if rx_path else {}
    names = _chebi_names(str(data_root / "reactions" / "tsv" / "chebiId_name.tsv"))
    ec2r = _ec_to_rhea(str(data_root / "reactions" / "tsv" / "rhea2ec.tsv"))

    steps = []
    for step in PATHWAYS[pathway]:
        rheas = ec2r.get(step["ec"], set())
        present = sorted(rheas & net)
        missing = sorted(rheas & gap_producers)
        status = "presente" if present else ("falta" if missing else "sin datos")
        rhea_id = (present or missing or [""])[0]
        eq = ""
        if rhea_id:
            rid = int(rhea_id.split(":")[1])
            r = rhea.get(rid, {})
            subs = " + ".join(names.get(c, c) for c in sorted(r.get("substrates", [])))[:80]
            prods = " + ".join(names.get(c, c) for c in sorted(r.get("products", [])))[:80]
            eq = f"{subs} → {prods}" if subs or prods else ""
        locus = ""
        evidence = ""
        for c in closures:
            if c.get("rhea") == rhea_id:
                locus = f"{c['genome']} · {c['locus']}"
                evidence = f"cierre etiquetado ({c.get('label_confidence', '')})"
                break
        if not locus:
            targets: set[str] = set()
            for r in rheas:
                targets |= r2a.get(r, set())
            hits = targets & set(by_acc)
            if hits:
                best = max((by_acc[a] for a in hits), key=lambda c: c.get("rank", 0))
                locus = f"{best['contig']}:{best['start']}-{best['end']}{best['strand']}"
                evidence = (f"ORF de novo; hit {best.get('hit')} "
                            f"({best.get('hit_pident')}% id, cov {best.get('hit_qcov')}); "
                            f"rank {best.get('rank')}")
        if not locus:
            for dom in EC_DOMAINS.get(step["ec"], set()):
                loc = domain_loci.get(genome, {}).get(dom)
                if loc:
                    locus = loc
                    evidence = f"dominio Pfam {dom} en ORF no anotado"
                    break
        steps.append({
            "enzyme": step["enzyme"], "ec": step["ec"], "label": step["label"],
            "status": status, "rhea": rhea_id, "equation": eq,
            "locus": locus, "evidence": evidence, "n_rhea": len(rheas),
        })
    return steps


# --------------------------------------------------------------------------
# Catálogo por categorías (escala a cientos de rutas)
# --------------------------------------------------------------------------

CATEGORIES: list[tuple[str, str]] = [
    ("Pigmentos", r"caroten|lycopen|melanin|violacein|prodigiosin|chlorophyll|bacteriochlorophyll|phycobilin|pigment"),
    ("Cofactores y vitaminas", r"riboflav|folate|pterin|thiamin|biotin|heme|porphyr|molybd|molibd|cobalamin|menaquinone|ubiquinone|\bNAD\b|\bFAD\b|\bFMN\b|\bCoA\b|lipoat|pyridox|siroheme"),
    ("Carbohidratos", r"sugar|glucose|fructose|galactose|ribose|xylose|mannose|sucrose|lactose|glucon|glucar|glycol|glyoxyl|pyruvate|glyoxylate|sedoheptulose|erythrose"),
    ("Lípidos", r"acyl-CoA|fatty|lipid|phosphatid|glycerol|ceramide|mycolic|isopren|terpen|geranyl|farnesyl|palmit|stear|oleoyl"),
    ("Nucleótidos", r"purine|pyrimidine|thymidine|cytidine|uridine|adenine|guanine|dUTP|dATP|dGTP|dCTP|ADP-ribose|nucleoside|nucleotide|inosine|xanthine|hypoxanthine"),
    ("Aminoácidos", r"^L-|amino acid|cystathion|prephenate|anthranilate|chorismate|lysine|arginine|proline|glutamate|aspartate|serine|threonine|methionine|cysteine|tryptophan|tyrosine|phenylalanine|valine|leucine|isoleucine|alanine|glycine|ornithine|citrulline"),
    ("Metabolismo central", r"\bTCA\b|citrate|isocitrate|succinate|fumarate|malate|oxaloacetate|acetyl-CoA|glyoxylate|pyruvate|2-oxoglutarate"),
]


def categorize(name: str) -> str:
    for cat, pat in CATEGORIES:
        import re
        if re.search(pat, name or "", re.I):
            return cat
    return "Otros"


def catalog(genome: str, out_root: Path) -> dict[str, list[dict]]:
    rows = [r for r in _read_tsv(out_root / "refined" / "refined_gaps.tsv")
            if r["genome"] == genome]
    out: dict[str, list[dict]] = {}
    for r in rows:
        out.setdefault(categorize(r["gap"]), []).append(r)
    return dict(sorted(out.items(), key=lambda kv: [c for c, _ in CATEGORIES].index(kv[0])
                       if kv[0] in [c for c, _ in CATEGORIES] else 99))


def _eq(reactions: dict, names: dict, rhea_id: str) -> str:
    if not rhea_id:
        return ""
    try:
        rid = int(rhea_id.split(":")[1])
    except ValueError:
        return ""
    r = reactions.get(rid, {})
    subs = " + ".join(names.get(c, c) for c in sorted(r.get("substrates", [])))[:70]
    prods = " + ".join(names.get(c, c) for c in sorted(r.get("products", [])))[:70]
    return f"{subs} → {prods}" if subs or prods else ""


def local_context(genome: str, gap_row: dict, data_root: Path, out_root: Path,
                  max_steps: int = 6) -> list[dict]:
    """Contexto local acotado de una brecha: aguas arriba, productores y consumidores."""
    net_path = data_root / "metabolism" / f"{genome}_network.json"
    net = {n["id"] for n in json.load(open(net_path)).get("nodes", [])
           if n.get("type") == "reaction"} if net_path.exists() else set()
    rx_path = data_root / "metabolism" / "rhea_reactions.pkl"
    if not rx_path.exists():
        rx_path = out_root.parent / "cache" / "rhea_reactions.pkl"
    reactions = _reactions(str(rx_path)) if rx_path.exists() else {}
    names = _chebi_names(str(data_root / "reactions" / "tsv" / "chebiId_name.tsv"))
    chebi = gap_row["chebi"]
    producers = [x for x in gap_row.get("producers", "").split(";") if x]
    consumers = [r for r in net if chebi in reactions.get(int(r.split(":")[1]), {}).get("substrates", set())] \
        if reactions else []

    cards: list[dict] = []
    # aguas arriba: reacciones presentes que producen los sustratos de los productores
    GENERIC = {"CHEBI:15377", "CHEBI:15378", "CHEBI:15379", "CHEBI:16526",
               "CHEBI:30616", "CHEBI:456216", "CHEBI:456215", "CHEBI:57945",
               "CHEBI:57783", "CHEBI:57287", "CHEBI:58349", "CHEBI:57540"}
    upstream: list[str] = []
    for r in producers[:2]:
        subs = sorted(reactions.get(int(r.split(":")[1]), {}).get("substrates", set()))[:2]
        for sub in subs:
            if sub in GENERIC:
                continue
            if names.get(sub, "") in {"A", "AH2", "H(+)", "H+", "H2O",
                                      "NADH", "NADPH", "NAD(+)", "NADP(+)", "CoA"}:
                continue
            for nr in net:
                if sub in reactions.get(int(nr.split(":")[1]), {}).get("products", set()):
                    if nr not in upstream:
                        upstream.append(nr)
                        break
            if len(upstream) >= 2:
                break
    for r in upstream[:2]:
        cards.append({"role": "aguas arriba", "status": "presente", "rhea": r,
                      "label": "reacción previa", "equation": _eq(reactions, names, r),
                      "locus": "", "evidence": "reacción presente en la red del genoma"})
    for r in producers[:2]:
        cards.append({"role": "productor", "status": "falta", "rhea": r,
                      "label": "enzima productora", "equation": _eq(reactions, names, r),
                      "locus": "", "evidence": "reacción ausente: define la brecha"})
    cards.append({"role": "brecha", "status": "falta", "rhea": "",
                  "label": gap_row["gap"], "equation": f"metabolito faltante: {gap_row['gap']}",
                  "locus": "", "evidence": f"producers={len(producers)}; consumers={len(consumers)}"})
    for r in consumers[:1]:
        cards.append({"role": "consumidor", "status": "presente", "rhea": r,
                      "label": "reacción posterior", "equation": _eq(reactions, names, r),
                      "locus": "", "evidence": "reacción presente en la red del genoma"})
    return cards[:max_steps]


# --------------------------------------------------------------------------
# Catálogo KEGG (cientos de rutas con nombre)
# --------------------------------------------------------------------------

@lru_cache(maxsize=2)
def _kegg_rows(path: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not Path(path).exists():
        return out
    for r in _read_tsv(Path(path)):
        d = out.setdefault(r["pathway_id"], {"name": r["pathway_name"], "rheas": set()})
        d["rheas"].add(r["rhea"])
    return out


def kegg_catalog(genome: str, data_root: Path, out_root: Path) -> list[dict]:
    net_path = data_root / "metabolism" / f"{genome}_network.json"
    net = {n["id"] for n in json.load(open(net_path)).get("nodes", [])
           if n.get("type") == "reaction"} if net_path.exists() else set()
    refined = [r for r in _read_tsv(out_root / "refined" / "refined_gaps.tsv")
               if r["genome"] == genome]
    gap_producers = set()
    for r in refined:
        gap_producers |= set(r.get("producers", "").split(";"))
    rows = []
    kegg_path = data_root / "kegg" / "pathway_reactions.tsv"
    if not kegg_path.exists():
        kegg_path = out_root.parent / "kegg" / "pathway_reactions.tsv"
    for pid, d in _kegg_rows(str(kegg_path)).items():
        present = d["rheas"] & net
        missing = d["rheas"] & gap_producers
        if present or missing:
            rows.append({"id": pid, "name": d["name"], "present": len(present),
                         "missing": len(missing), "total": len(d["rheas"])})
    rows.sort(key=lambda r: (r["missing"], r["present"]), reverse=True)
    return rows


def kegg_pathway_steps(genome: str, pathway_id: str, data_root: Path,
                       out_root: Path, max_steps: int = 12) -> list[dict]:
    net_path = data_root / "metabolism" / f"{genome}_network.json"
    net = {n["id"] for n in json.load(open(net_path)).get("nodes", [])
           if n.get("type") == "reaction"} if net_path.exists() else set()
    refined = [r for r in _read_tsv(out_root / "refined" / "refined_gaps.tsv")
               if r["genome"] == genome]
    gap_producers = set()
    for r in refined:
        gap_producers |= set(r.get("producers", "").split(";"))
    rx_path = data_root / "metabolism" / "rhea_reactions.pkl"
    if not rx_path.exists():
        rx_path = out_root.parent / "cache" / "rhea_reactions.pkl"
    reactions = _reactions(str(rx_path)) if rx_path.exists() else {}
    names = _chebi_names(str(data_root / "reactions" / "tsv" / "chebiId_name.tsv"))
    kegg_path = data_root / "kegg" / "pathway_reactions.tsv"
    if not kegg_path.exists():
        kegg_path = out_root.parent / "kegg" / "pathway_reactions.tsv"
    d = _kegg_rows(str(kegg_path)).get(pathway_id, {})
    steps = []
    for rhea in sorted(d.get("rheas", set())):
        status = "presente" if rhea in net else ("falta" if rhea in gap_producers else "sin datos")
        steps.append({
            "rhea": rhea, "status": status,
            "equation": _eq(reactions, names, rhea),
        })
    order = {"falta": 0, "presente": 1, "sin datos": 2}
    steps.sort(key=lambda s: order[s["status"]])
    return steps[:max_steps]


def kegg_catalog_from_net(net: set[str], gap_producers: set[str],
                          data_root: Path, out_root: Path) -> list[dict]:
    """Catálogo KEGG para una red construida dinámicamente (genoma subido)."""
    kegg_path = data_root / "kegg" / "pathway_reactions.tsv"
    if not kegg_path.exists():
        kegg_path = out_root.parent / "kegg" / "pathway_reactions.tsv"
    rows = []
    for pid, d in _kegg_rows(str(kegg_path)).items():
        present = d["rheas"] & net
        missing = d["rheas"] & gap_producers
        if present or missing:
            rows.append({"id": pid, "name": d["name"], "present": len(present),
                         "missing": len(missing), "total": len(d["rheas"])})
    rows.sort(key=lambda r: (r["missing"], r["present"]), reverse=True)
    return rows
