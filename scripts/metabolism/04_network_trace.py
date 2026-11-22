#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_network_trace.py — Auto-descubrimiento de rutas y detección de brechas. (v3.3)

Encuentra fuentes importables, brechas (falta de enzima) y productos finales
(sinks) de una red metabólica RHEA por genoma, cruzando contra todo Rhea, y
escribe data/metabolism/{genome}_network_trace.json.

================================================================================
SEMÁNTICA DE DIRECCIÓN (v3.1) — por qué el cableado ya no es "todo o nada"
================================================================================
Rhea es deliberadamente NEUTRA en dirección: cada familia se publica como
maestra (UN) + variante LR + variante RL + variante BI (4 miembros SIEMPRE,
verificado en rhea.rdf.gz: 19005 maestros x 4 = 76020). La presencia de la
variante BI NO indica reversibilidad fisiológica, y ni el RDF ni
rhea-directions.tsv codifican qué dirección es la válida (transformableTo solo
enlaza el lado L con el R de la misma familia).

La dirección curada existe en rhea2uniprot.tsv: UniProt añade filas LR o RL
(para el maestro) SOLO cuando la actividad catalítica de esa proteína tiene
dirección definida; si es reversible, solo hay fila UN. Reglas de cableado:

  1. Id del genoma LR o RL -> UN sentido, tal cual su ecuación (el id ya es la
     forma direccional curada).
  2. Id UN o BI -> VOTO UniProt para su maestro:
       - Si la red del genoma contiene ids UniProt (auto-detección en nodes y
         edges), se usa el voto de las proteínas DEL GENOMA.
       - Si no, se usa el voto GLOBAL de UniProt (todas las proteínas).
       votos {LR}  -> un sentido natural (L->R);
       votos {RL}  -> un sentido invertido (R->L);
       votos vacíos o LR+RL mixtos -> AMBOS sentidos (reversible).
  3. Sin rhea2uniprot.tsv y sin proteínas en la red -> UN/BI se cablean en
     ambos sentidos (red químicamente cerrada: brechas/sinks ~ 0; se avisa).

Modos: --all-reversible (todo bidireccional: debe dar 0/0/0 en un genoma
completo; diagnóstico) y --no-reversible (v2: todo un sentido canónico).

Modo batch: si se omite --genome se procesan TODOS los *_network.json de
data/metabolism (el RDF se parsea una sola vez; el voto UniProt también).

Historial:
  v1 leía rhea-relationships.tsv con un formato inexistente -> grafo vacío.
  v2 cacheó participantes/ecuaciones del RDF pero cableó TODO en un sentido
     canónico -> en E. coli K-12 generó 210 brechas y 211 sinks falsos.
  v3 intentó corregirlo asumiendo "familia con BI = reversible" (falso: BI
     existe para todas) -> 819/819 bidireccionales y 0/0/0/0.
  v3.1 usa el voto UniProt (verificado: glucokinasa->LR, PEP carboxilasa->RL,
     MDH/fumarasa/transaldolasa->mixto, catalasa->LR, hisB fosfatasa->LR).
  v3.2 corrige el ÍNDICE Rhea-wide: producers_of/consumers_of contaban las
     variantes RL de familias que CONSUMEN un metabolito como "productoras"
     (v3.1 real: UTP "lo producen 42 familias" — la mayoría eran galU y otras
     consumidoras al revés; glycine betaine con betA/betB invertidas; 72/72
     sinks marcados como brecha aguas abajo). Ahora el índice es FISIOLÓGICO
     por familia (voto UniProt), un sink no se marca como brecha aguas abajo
     por sus propias familias, y --resolve-generics expande también el índice
     (la familia ndk aparece como productor real de UTP vía su clase).
  v3.2.1 anota en consola las familias reversibles/sin voto cuya ecuación
     canónica muestra el metabolito al revés (UTP ya no parece "producido" por
     sus consumidores: la nota dice que se cuenta el sentido inverso).
  v3.3 TRÍA de brechas en 4 cubetas accionables (consola + reporte): closable
     (enzima ya en la red: cableado/voto), clase genérica (artefacto probable),
     productor genérico (resoluble con inferencia is-a de ChEBI) y genuina
     (sin productor específico en Rhea). gen_map se calcula SIEMPRE (barato)
     para que UTP caiga en 'productor genérico' aunque no se use
     --resolve-generics; el cableado sigue expandiendo clases solo con el flag.

Descarga (una sola vez):
    curl -o data/metabolism/tsv/rhea.rdf.gz \
        https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz
    curl -o data/metabolism/tsv/rhea2uniprot.tsv \
        https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv
"""
import argparse
import gzip
import json
import logging
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime
from pathlib import Path

try:
    from networkx import DiGraph
except ImportError:  # shim mínimo para no depender de networkx
    class DiGraph:
        def __init__(self):
            self._succ, self._pred = {}, {}

        def add_edge(self, u, v, **attr):
            self._succ.setdefault(u, {})[v] = attr
            self._pred.setdefault(v, {})[u] = attr
            self._succ.setdefault(v, {})
            self._pred.setdefault(u, {})

        def has_edge(self, u, v):
            return v in self._succ.get(u, {})

        def number_of_nodes(self):
            return len(self._succ)

        def number_of_edges(self):
            return sum(len(s) for s in self._succ.values())

        def nodes(self):
            return list(self._succ)

        def in_degree(self, n):
            return len(self._pred.get(n, {}))

        def out_degree(self, n):
            return len(self._succ.get(n, {}))

        @property
        def pred(self):
            return self._pred


SCRIPT_VERSION = "3.3"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-7s | %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("network_trace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RHEA_DIR = PROJECT_ROOT / "data" / "metabolism" / "tsv"
NET_DIR = PROJECT_ROOT / "data" / "metabolism"


def load_env_file(path=None):
    """Carga KEY=VALUE de un .env (por defecto PROJECT_ROOT/.env) SIN
    sobreescribir lo que ya viene del entorno. Los valores con comillas se
    limpian; '#' inicia comentario. Nunca registra ni devuelve los valores."""
    path = Path(path) if path else PROJECT_ROOT / ".env"
    if not path.exists():
        return []
    loaded = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key, val = key.strip(), val.strip()
        if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
            val = val[1:-1]
        if key and key not in os.environ:
            os.environ[key] = val
            loaded.append(key)
    if loaded:
        log.info(f".env cargado ({path.name}): claves "
                 f"{', '.join(sorted(loaded))}")
    return loaded

RDF_GZ = RHEA_DIR / "rhea.rdf.gz"
PARTICIPANTS_TSV = RHEA_DIR / "rhea_participants.tsv"   # cache: RHEA_ID, TYPE, ROLE, CHEBI_ID, TRANSPORT[, MASTER_ID]
EQUATIONS_TSV = RHEA_DIR / "rhea_equations.tsv"         # cache: RHEA_ID, TYPE, EQUATION
RDF_DOWNLOAD = "https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz"

NS_RDF = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
NS_RH = "{http://rdf.rhea-db.org/}"

SIDE_RE = re.compile(r"(\d+)_(L|R)")

# Metabolitos de moneda (co-sustratos ubícuos): saturan el análisis topológico
# y generan aristas espurias. Rhea usa microespecies pH 7.3 (ATP=30616...).
CURRENCY_CHEBI = {
    "CHEBI:15377",   # H2O
    "CHEBI:15378",   # H(+)
    "CHEBI:30616",   # ATP
    "CHEBI:456216",  # ADP
    "CHEBI:456215",  # AMP
    "CHEBI:57540",   # NAD(+)
    "CHEBI:57945",   # NADH
    "CHEBI:58349",   # NADP(+)
    "CHEBI:57783",   # NADPH
    "CHEBI:57287",   # CoA
    "CHEBI:57288",   # acetyl-CoA
    "CHEBI:57692",   # FAD
    "CHEBI:58307",   # FADH2
    "CHEBI:58210",   # FMN
    "CHEBI:57618",   # FMNH2
    "CHEBI:16526",   # CO2
    "CHEBI:15379",   # O2
    "CHEBI:43474",   # phosphate
    "CHEBI:33019",   # diphosphate
    "CHEBI:28938",   # NH4(+)
    "CHEBI:17996",   # chloride
    "CHEBI:18420",   # Mg(2+)
    "CHEBI:29103",   # K(+)
    "CHEBI:29101",   # Na(+)
    "CHEBI:29108",   # Ca(2+)
    "CHEBI:16189",   # sulfate
}

RHEA_ID_RE = re.compile(r"(\d+)\s*$")

# Clases genéricas de Rhea: 'a carboxylate', 'an alcohol', '[protein]-X', ...
GENERIC_NAME_RE = re.compile(r"^(an?\s+[a-z(]|\[protein)")

# Cubetas de la tría de brechas (v3.3): clave = valor del campo 'triage' del
# reporte; orden = prioridad de clasificación en classify_gap.
TRIAGE_LABELS = {
    "closable_with_genome_enzyme": "closable (enzima ya en la red)",
    "generic_class": "clase genérica (artefacto probable)",
    "generic_producer": "productor genérico (is-a de ChEBI pendiente)",
    "genuine_missing_enzyme": "genuina (sin productor específico en Rhea)",
}

# Clases genéricas de nucleótidos -> instancias (para --resolve-generics).
GENERIC_NUCLEOTIDE_PATTERNS = [
    (re.compile(r"2['’′-]?deoxy.*(tri|di|mono)phosphate", re.I), {
        "triphosphate": ["dATP", "dGTP", "dCTP", "dTTP", "dITP"],
        "diphosphate": ["dADP", "dGDP", "dCDP", "dTDP"],
        "monophosphate": ["dAMP", "dGMP", "dCMP", "dTMP"]}),
    (re.compile(r"(ribonucleoside|nucleoside).*triphosphate", re.I),
     {"triphosphate": ["ATP", "GTP", "CTP", "UTP", "ITP", "XTP"]}),
    (re.compile(r"(ribonucleoside|nucleoside).*diphosphate", re.I),
     {"diphosphate": ["ADP", "GDP", "CDP", "UDP", "IDP", "XDP"]}),
    (re.compile(r"(ribonucleoside|nucleoside).*monophosphate", re.I),
     {"monophosphate": ["AMP", "GMP", "CMP", "UMP", "IMP", "XMP"]}),
]

# Accesiones UniProtKB (6 o 10 caracteres, formato estándar)
UNIPROT_ACC_RE = re.compile(
    r"\b(?:[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9])"
    r"(?:[A-Z][A-Z0-9]{3}[0-9])?\b")


def normalize_rhea_id(raw):
    """Acepta 'RHEA:10001', 'rhea:10001', '10001' -> ('RHEA:10001', 10001)."""
    m = RHEA_ID_RE.search(str(raw).strip())
    if not m:
        return None
    n = int(m.group(1))
    return f"RHEA:{n}", n


_MEMO_CACHE = {}          # (fichero, mtime, size) -> objeto parseado


def _memo_parse(path, loader):
    """Parsea `path` una sola vez por (mtime, size); el 07 construye payloads
    por genoma y sin esto re-leería el TSV de Rhea en cada uno (minutos)."""
    st = path.stat()
    key = (str(path), st.st_mtime_ns, st.st_size)
    hit = _MEMO_CACHE.get(key)
    if hit is None:
        hit = loader()
        _MEMO_CACHE[key] = hit
    return hit


def load_chebi_names(rhea_dir):
    file = rhea_dir / "chebiId_name.tsv"
    if not file.exists():
        return {}

    def _parse():
        names = {}
        with open(file, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and (parts[0].startswith("CHEBI:") or
                                        parts[0].strip().isdigit()):
                    names[parts[0].strip()] = parts[1].strip()
        return names

    return _memo_parse(file, _parse)


def parse_rhea_rdf(rdf_path):
    """
    Parseo en streaming de rhea.rdf.gz, en DOS fases:
      Reaction (maestro UN)       --rh:side------------------->  ReactionSide ({master}_L / {master}_R)
      DirectionalReaction (LR/RL) --rh:substrates / rh:products-->  ReactionSide
      BidirectionalReaction (BI)  --rh:substratesOrProducts-->  ReactionSide
      ReactionSide                --rh:contains[N]----------->  ReactionParticipant
      ReactionParticipant         --rh:compound-------------->  Compound
      Compound                    --rh:accession = "CHEBI:xxxxx"
    El maestro de cada variante se deduce del prefijo de sus ReactionSide.
    """
    sides = {}
    participant_comp = {}
    compounds = {}
    eq_of = {}
    refs_of = {}
    obsolete_ids = set()
    transport_ids = set()

    with gzip.open(rdf_path, "rb") as fh:
        it = ET.iterparse(fh, events=("start", "end"))
        _, root = next(it)
        for event, elem in it:
            if event != "end" or elem.tag != f"{NS_RDF}Description":
                continue
            key = elem.get(f"{NS_RDF}about", "").rsplit("/", 1)[-1]

            if SIDE_RE.fullmatch(key):
                pids = {c.get(f"{NS_RDF}resource", "").rsplit("/", 1)[-1]
                        for c in elem
                        if c.tag.startswith(NS_RH) and
                        c.tag[len(NS_RH):].startswith("contains")}
                pids.discard("")
                sides.setdefault(key, set()).update(pids)

            elif key.startswith("Participant_"):
                for c in elem:
                    if c.tag == f"{NS_RH}compound":
                        cid = c.get(f"{NS_RDF}resource", "").rsplit("/", 1)[-1]
                        if cid:
                            participant_comp[key] = cid

            elif key.startswith("Compound_"):
                for c in elem:
                    if c.tag == f"{NS_RH}accession":
                        acc = (c.text or "").strip()
                        if acc.startswith("CHEBI:"):
                            compounds[key] = acc

            elif key.isdigit():
                rid = int(key)
                refs = refs_of.setdefault(rid, {"subs": set(), "prods": set(),
                                                "both": set(), "own": set()})
                for c in elem:
                    tag = c.tag[len(NS_RH):] if c.tag.startswith(NS_RH) else ""
                    res = c.get(f"{NS_RDF}resource")
                    if tag == "equation":
                        eq_of[rid] = (c.text or "").strip()
                    elif tag == "status" and res and \
                            res.rsplit("/", 1)[-1] == "Obsolete":
                        obsolete_ids.add(rid)
                    elif tag == "isTransport" and \
                            (c.text or "").strip().lower() == "true":
                        transport_ids.add(rid)
                    elif res:
                        sk = res.rsplit("/", 1)[-1]
                        if SIDE_RE.fullmatch(sk):
                            if tag == "substrates":
                                refs["subs"].add(sk)
                            elif tag == "products":
                                refs["prods"].add(sk)
                            elif tag == "substratesOrProducts":
                                refs["both"].add(sk)
                            elif tag == "side":
                                refs["own"].add(sk)
            root.clear()

    def side_chebis(side_keys):
        out = set()
        for sk in side_keys:
            for pid in sides.get(sk, ()):
                cid = participant_comp.get(pid)
                if cid and cid in compounds:
                    out.add(compounds[cid])
        return out

    reactions, master_of = {}, {}
    for rid, refs in refs_of.items():
        subs_sides, prods_sides = refs["subs"], refs["prods"]
        if refs["own"] and not subs_sides and not refs["both"]:
            rtype = "UN"
            subs_sides = {s for s in refs["own"] if s.endswith("_L")}
            prods_sides = {s for s in refs["own"] if s.endswith("_R")}
            master_of[rid] = rid
        elif refs["both"] and not subs_sides and not prods_sides:
            rtype = "BI"
            subs_sides = {s for s in refs["both"] if s.endswith("_L")}
            prods_sides = {s for s in refs["both"] if s.endswith("_R")}
        else:
            sample = next(iter(subs_sides)) if subs_sides else \
                next(iter(prods_sides), "")
            rtype = "LR" if sample.endswith("_L") else "RL"
        if rtype != "UN":
            m = SIDE_RE.match(sample) if sample else None
            if m:
                master_of[rid] = int(m.group(1))
        if rid in obsolete_ids:
            continue
        subs, prods = side_chebis(subs_sides), side_chebis(prods_sides)
        if not subs or not prods:
            continue
        reactions[rid] = {
            "type": rtype,
            "substrates": subs,
            "products": prods,
            "transport": rid in transport_ids,
        }
    return reactions, eq_of, master_of


def load_reaction_data(rhea_dir, force=False):
    """Carga participantes+ecuaciones. Cache TSV si existe (v2 o v3); si no, RDF.
    El TSV se parsea UNA sola vez por (mtime, size): los dict se comparten
    entre genomas — solo lectura (nobody muta substrates/products)."""
    part_tsv = rhea_dir / PARTICIPANTS_TSV.name
    eq_tsv = rhea_dir / EQUATIONS_TSV.name
    rdf_gz = rhea_dir / RDF_GZ.name

    if part_tsv.exists() and eq_tsv.exists() and not force:
        def _parse():
            log.info(f"Leyendo cache de participantes: {part_tsv.name}")
            reactions, equations, master_of = {}, {}, {}
            with open(part_tsv, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 4:
                        rid, rtype, role, chebi = (int(parts[0]), parts[1],
                                                   parts[2], parts[3])
                        is_transport = len(parts) >= 5 and \
                            parts[4].strip().lower() == "true"
                        entry = reactions.setdefault(
                            rid, {"type": rtype, "substrates": set(),
                                  "products": set(), "transport": is_transport})
                        entry[role].add(chebi)
                        if is_transport:
                            entry["transport"] = True
                        if len(parts) >= 6 and parts[5].strip().isdigit():
                            master_of[rid] = int(parts[5])
            with open(eq_tsv, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3:
                        equations[int(parts[0])] = parts[2]
            return reactions, equations, master_of

        return _memo_parse(part_tsv, _parse)

    if not rdf_gz.exists():
        log.error(f"No se encontró {rdf_gz}")
        log.error("Descárgalo (una sola vez) con:")
        log.error(f"    curl -o {rdf_gz} {RDF_DOWNLOAD}")
        return {}, {}, {}

    log.info(f"Parseando {rdf_gz.name} (solo la primera vez; luego se cachea)...")
    reactions, equations, master_of = parse_rhea_rdf(rdf_gz)
    log.info(f"RDF parseado: {len(reactions)} reacciones con participantes.")

    with open(part_tsv, "w", encoding="utf-8") as f:
        f.write("# RHEA_ID\tDIRECTION_TYPE\tROLE\tCHEBI_ID\tTRANSPORT\tMASTER_ID\n")
        for rid in sorted(reactions):
            e = reactions[rid]
            for role in ("substrates", "products"):
                for chebi in sorted(e[role]):
                    f.write(f"{rid}\t{e['type']}\t{role}\t{chebi}\t"
                            f"{'true' if e.get('transport') else 'false'}\t"
                            f"{master_of.get(rid, '')}\n")
    with open(eq_tsv, "w", encoding="utf-8") as f:
        f.write("# RHEA_ID\tDIRECTION_TYPE\tEQUATION\n")
        for rid in sorted(equations):
            f.write(f"{rid}\t{reactions.get(rid, {}).get('type', 'UN')}\t"
                    f"{equations[rid]}\n")
    log.info(f"Cache escrito: {part_tsv.name}, {eq_tsv.name}")
    return reactions, equations, master_of


# ------------------------------------------------------------------ familias
def compute_families(reactions, master_of):
    """Devuelve (family_of, family_members): familia = maestro RHEA.

    Con MASTER_ID (cache v3 o RDF) se usa directamente; si no, se reconstruye
    por firma de participantes (la variante LR comparte lados con el maestro;
    la RL los tiene intercambiados).
    """
    family_of = {}
    if master_of:
        for rid in reactions:
            family_of[rid] = master_of.get(rid, rid)
    else:
        parent = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        sig2ids = defaultdict(list)
        for rid, e in reactions.items():
            sig2ids[(frozenset(e["substrates"]),
                     frozenset(e["products"]))].append(rid)
        for ids in sig2ids.values():
            for rid in ids:
                union(ids[0], rid)
        for (left, right), ids in sig2ids.items():
            swapped = (right, left)
            if swapped in sig2ids:
                union(ids[0], sig2ids[swapped][0])
        reps = {}
        for rid in reactions:
            reps.setdefault(find(rid), []).append(rid)
        for members in reps.values():
            rep = next((r for r in sorted(members)
                        if reactions[r]["type"] == "UN"), min(members))
            for r in members:
                family_of[r] = rep

    members = defaultdict(list)
    for rid, fam in family_of.items():
        members[fam].append(rid)
    return family_of, members


# ------------------------------------------------------------- voto UniProt
def find_uniprot_file(rhea_dir):
    for name in ("rhea2uniprot.tsv", "rhea2uniprot.tsv.gz",
                 "rhea2uniprot_sprot.tsv", "rhea2uniprot_sprot.tsv.gz"):
        p = rhea_dir / name
        if p.exists():
            return p
    return None


def read_uniprot_votes(path, want_accessions=None):
    """rhea2uniprot -> {master_id: {dirs LR/RL}}.

    UniProt solo registra LR/RL cuando la actividad catalítica tiene dirección
    definida; filas UN se ignoran. Con want_accessions se limita a las
    proteínas del genoma (voto organismo); sin él, voto global.
    """
    votes = defaultdict(set)
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            header = None
            i_dir, i_master, i_id = 1, 2, 3
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    if line.startswith("#"):
                        continue
                    header = [h.strip().upper() for h in parts]
                    if all(x in header for x in ("DIRECTION", "MASTER_ID", "ID")):
                        i_dir = header.index("DIRECTION")
                        i_master = header.index("MASTER_ID")
                        i_id = header.index("ID")
                    continue
                if len(parts) <= max(i_dir, i_master, i_id):
                    continue
                if want_accessions is not None and \
                        parts[i_id].strip() not in want_accessions:
                    continue
                d = parts[i_dir].strip().upper()
                if d in ("LR", "RL") and parts[i_master].strip().isdigit():
                    votes[int(parts[i_master])].add(d)
    except (OSError, EOFError) as exc:
        log.warning(f"No se pudo leer {path.name}: {exc}")
        return {}
    return dict(votes)


def extract_uniprot_ids(obj, acc=None):
    """Busca accesiones UniProtKB en toda la estructura de la red (nodes/edges)."""
    if acc is None:
        acc = set()
    if isinstance(obj, dict):
        for v in obj.values():
            extract_uniprot_ids(v, acc)
    elif isinstance(obj, (list, tuple, set)):
        for v in obj:
            extract_uniprot_ids(v, acc)
    elif isinstance(obj, str):
        acc.update(UNIPROT_ACC_RE.findall(obj))
    return acc


def vote_direction(family, org_votes, global_votes):
    """-> 'both' | 'one_lr' (sentido natural) | 'one_rl' (invertido)."""
    dirs = set()
    if family in org_votes:
        dirs = set(org_votes[family])
    elif family in global_votes:
        dirs = set(global_votes[family])
    dirs.discard("UN")
    if not dirs or dirs == {"LR", "RL"}:
        return "both"
    if dirs == {"RL"}:
        return "one_rl"
    return "one_lr"


def resolve_direction(entry, rid, family_of, org_votes, global_votes,
                      reversible_mode, all_reversible):
    """-> False (un sentido, lados propios), 'reverse' (solo invertido) o True."""
    if all_reversible:
        return True
    if not reversible_mode:
        return False                       # modo v2: todo un sentido canónico
    if entry["type"] in ("LR", "RL"):
        return False                       # id direccional: un sentido propio
    v = vote_direction(family_of.get(rid, rid), org_votes, global_votes)
    if v == "both":
        return True
    if v == "one_rl":
        return "reverse"
    return False                           # one_lr: sentido natural (L->R)


def resolve_nucleotide_generics(reactions, names):
    """--resolve-generics: clases genéricas de nucleótidos -> instancias ChEBI."""
    name_index = defaultdict(list)
    for cid, nm in names.items():
        name_index[nm.lower()].append(cid)
    result = {}
    seen = set()
    for e in reactions.values():
        for c in e["substrates"] | e["products"]:
            if c in seen:
                continue
            seen.add(c)
            nm = names.get(c, "")
            if not GENERIC_NAME_RE.match(nm):
                continue
            for pat, table in GENERIC_NUCLEOTIDE_PATTERNS:
                if not pat.search(nm):
                    continue
                kind = re.search(r"(triphosphate|diphosphate|monophosphate)",
                                 nm, re.I)
                key = kind.group(1).lower() if kind else None
                instances = set()
                for iname in table.get(key, []):
                    instances.update(name_index.get(iname.lower(), []))
                instances.discard(c)
                if instances:
                    result[c] = sorted(instances)
                break
    return result


def add_reaction_edges(G, entry, rid, filter_currency, currency_dropped,
                       direction, gen_map=None, gen_stats=None):
    """Añade aristas Sustrato->Producto según la dirección resuelta.

    Las reacciones de transporte no generan aristas (mueven el MISMO CHEBI
    entre compartimentos). Las de moneda se excluyen de ambos extremos.
    direction: False = un sentido (lados propios); True = ambos; 'reverse' = solo invertido.
    """
    if entry.get("transport"):
        return 0
    subs = entry["substrates"]
    prods = entry["products"]

    def keep(c):
        return not filter_currency or c not in CURRENCY_CHEBI

    def note_currency(c):
        if filter_currency and c in CURRENCY_CHEBI and currency_dropped is not None:
            currency_dropped.add(c)

    pairs = [(s, p) for s in subs for p in prods if s != p]
    both = False
    if direction == "reverse":
        pairs = [(p, s) for (s, p) in pairs]
    elif direction is True:
        both = True
    n = 0
    for s, p in pairs:
        note_currency(s)
        note_currency(p)
        if not keep(s) or not keep(p):
            continue
        if not G.has_edge(s, p):
            G.add_edge(s, p, enzyme=f"RHEA:{rid}")
        n += 1
        if both and not G.has_edge(p, s):
            G.add_edge(p, s, enzyme=f"RHEA:{rid}")
    # --resolve-generics: conectar instancias específicas con el otro lado
    if gen_map:
        sk = [s for s in subs if keep(s)]
        pk = [p for p in prods if keep(p)]
        for s in list(sk) + list(pk):
            if s in gen_map:
                others = pk if s in subs else sk
                for inst in gen_map[s]:
                    if not keep(inst):
                        continue
                    for o in others:
                        if inst == o or not keep(o):
                            continue
                        a, b = (inst, o) if s in subs else (o, inst)
                        if not G.has_edge(a, b):
                            G.add_edge(a, b, enzyme=f"RHEA:{rid}·gen")
                            if gen_stats is not None:
                                gen_stats.add(a)
                            n += 1
                        if both and not G.has_edge(b, a):
                            G.add_edge(b, a, enzyme=f"RHEA:{rid}·gen")
    return n


def master_orientation(entry):
    """Lados en orientación maestro: (sustratos_maestro, productos_maestro)."""
    if entry["type"] == "RL":
        return entry["products"], entry["substrates"]
    return entry["substrates"], entry["products"]


def family_direction_sides(fam, entry, global_votes):
    """(productos, sustratos) de la familia según su voto (orientación maestro).

    Sin voto o con votos LR+RL mixtos cuenta ambos lados (reversible); voto RL
    intercambia los lados del maestro; voto LR deja el sentido natural.
    """
    msubs, mprods = master_orientation(entry)
    dirs = set(global_votes.get(fam, ()))
    dirs.discard("UN")
    if not dirs or dirs == {"LR", "RL"}:
        return mprods | msubs, msubs | mprods
    if dirs == {"RL"}:
        return msubs, mprods
    return mprods, msubs


def build_physiological_index(reactions, family_of, global_votes, gen_map=None):
    """Índice FISIOLÓGICO por familia: productores/consumidores solo en las
    direcciones válidas (voto UniProt).

    En el índice en bruto, el hijo RL de una familia que consume X aparece con
    X en 'products' y contaminaba el índice (v3.1 real: UTP con 42
    'productores', glycine betaine 'consumible' por betA/betB al revés). Aquí
    cada familia registra únicamente sus direcciones curadas; las familias sin
    voto se tratan como reversibles (ambos lados).
    """
    fam_members = defaultdict(list)
    for rid, fam in family_of.items():
        fam_members[fam].append(rid)
    fam_products, fam_substrates = defaultdict(set), defaultdict(set)
    for fam, mem in fam_members.items():
        rep = fam if fam in reactions else mem[0]
        prods, subs = family_direction_sides(fam, reactions[rep], global_votes)
        for c in prods:
            fam_products[c].add(fam)
        for c in subs:
            fam_substrates[c].add(fam)
        if gen_map:                             # instancias de clases genéricas
            for c in prods:
                for inst in gen_map.get(c, ()):
                    fam_products[inst].add(fam)
            for c in subs:
                for inst in gen_map.get(c, ()):
                    fam_substrates[inst].add(fam)
    producers_of = {c: sorted(v) for c, v in fam_products.items()}
    consumers_of = {c: sorted(v) for c, v in fam_substrates.items()}
    return producers_of, consumers_of


def build_raw_index(reactions, gen_map=None):
    """Índice en bruto (modos v2 y --all-reversible): todas las variantes."""
    producers_of, consumers_of = {}, {}
    for rid, e in reactions.items():
        if e.get("transport"):
            continue
        for c in e["products"]:
            producers_of.setdefault(c, []).append(rid)
        for c in e["substrates"]:
            consumers_of.setdefault(c, []).append(rid)
    if gen_map:
        for c, insts in gen_map.items():
            for rid in producers_of.get(c, []):
                for inst in insts:
                    producers_of.setdefault(inst, []).append(rid)
            for rid in consumers_of.get(c, []):
                for inst in insts:
                    consumers_of.setdefault(inst, []).append(rid)
    return producers_of, consumers_of


def direction_note(fam, chem, reactions, global_votes, as_producer):
    """Nota de dirección para el display: aclara por qué la familia cuenta.

    La ecuación canónica del maestro puede mostrar el metabolito en el lado
    contrario al que se contó en el índice (familias sin voto = reversibles,
    voto RL = sentido inverso). Devuelve una nota corta o "".
    """
    entry = reactions.get(fam)
    if entry is None:
        return ""
    msubs, mprods = master_orientation(entry)
    opposite = chem in (msubs if as_producer else mprods)
    if not opposite:
        return ""
    dirs = set(global_votes.get(fam, ()))
    dirs.discard("UN")
    if dirs == {"RL"}:
        return " [voto UniProt RL: cuenta el sentido inverso]"
    if not dirs or dirs == {"LR", "RL"}:
        return " [familia reversible: cuenta también el sentido inverso]"
    return ""


def classify_gap(c, fam_prod, genome_families, family_of, reactions,
                 global_votes, names, generic_flag, gen_instances=()):
    """Tría de una brecha -> cubeta accionable (ver TRIAGE_LABELS).

    Prioridad: la enzima ya está en la red (cableado/voto) > el metabolito es
    él mismo una clase genérica > es instancia de una clase genérica conocida
    o algún productor lo da vía una clase genérica (ambos resolubles con
    inferencia is-a de ChEBI / --resolve-generics) > genuina.
    """
    if any(family_of.get(x, x) in genome_families for x in fam_prod):
        return "closable_with_genome_enzyme"
    if generic_flag:
        return "generic_class"
    if c in gen_instances:
        return "generic_producer"
    for fam in fam_prod:
        entry = reactions.get(fam)
        if entry is None:
            continue
        prods, _ = family_direction_sides(fam, entry, global_votes)
        if any(p != c and GENERIC_NAME_RE.match(names.get(p, ""))
               for p in prods):
            return "generic_producer"
    return "genuine_missing_enzyme"


def family_list(ids, family_of):
    """Dedup por familia: devuelve los REPRESENTANTES (maestro preferido)."""
    fams = {}
    for rid in ids:
        fam = family_of.get(rid, rid)
        fams.setdefault(fam, rid)
    return sorted(fams)


# ------------------------------------------------------------------ análisis
def process_genome(genome, net_dir, ctx, args):
    """Analiza una red {genome}_network.json y escribe su reporte."""
    reactions, equations, names = ctx["reactions"], ctx["equations"], ctx["names"]
    family_of = ctx["family_of"]

    net_file = net_dir / f"{genome}_network.json"
    if not net_file.exists():
        log.warning(f"[{genome}] no existe {net_file.name}; se omite")
        return
    with open(net_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    # voto de dirección con las proteínas del propio genoma (si están en la red)
    up_ids = extract_uniprot_ids(data.get("nodes", []))
    extract_uniprot_ids(data.get("edges", []), up_ids)
    org_votes = {}
    if up_ids and ctx["up_file"]:
        org_votes = read_uniprot_votes(ctx["up_file"], want_accessions=up_ids)
        log.info(f"[{genome}] {len(up_ids)} ids UniProt en la red -> voto de "
                 f"dirección por proteína ({len(org_votes)} familias con voto)")
    elif ctx["up_file"]:
        log.info(f"[{genome}] sin ids UniProt en la red -> voto de dirección "
                 f"GLOBAL de UniProt ({len(ctx['global_votes'])} familias con voto)")

    genome_rheas, skipped = set(), []
    for n in data.get("nodes", []):
        if n.get("type") != "reaction":
            continue
        norm = normalize_rhea_id(n.get("id", ""))
        if norm:
            genome_rheas.add(norm[1])
        else:
            skipped.append(str(n.get("id")))

    G = DiGraph()
    currency_dropped = set()
    missing_in_rhea, matched = [], []
    n_both = n_one = n_rev = 0
    type_counts = defaultdict(int)
    gen_stats = set()
    for rid in sorted(genome_rheas):
        entry = reactions.get(rid)
        if entry is None:
            missing_in_rhea.append(rid)
            continue
        matched.append(rid)
        type_counts[entry["type"]] += 1
        direction = resolve_direction(entry, rid, family_of,
                                      org_votes, ctx["global_votes"],
                                      ctx["reversible_mode"],
                                      args.all_reversible)
        if direction is True:
            n_both += 1
        elif direction == "reverse":
            n_rev += 1
        else:
            n_one += 1
        add_reaction_edges(G, entry, rid, not args.keep_currency,
                           currency_dropped, direction, ctx["wiring_gen_map"],
                           gen_stats)

    log.info(f"[{genome}] Analizando ({len(genome_rheas)} reacciones "
             f"enzimáticas; tipos: {dict(sorted(type_counts.items()))})...")
    if skipped:
        log.warning(f"[{genome}] {len(skipped)} ids sin formato RHEA numérico "
                    f"fueron ignorados: {skipped[:5]}")
    if missing_in_rhea:
        log.warning(f"[{genome}] {len(missing_in_rhea)} RHEA ids sin "
                    f"participantes en el RDF (¿obsoletos o sin estequiometría?): "
                    f"{['RHEA:%d' % r for r in missing_in_rhea[:5]]}")
    if not type_counts.get("LR") and not type_counts.get("RL") and not up_ids \
            and not ctx["global_votes"]:
        log.warning(f"[{genome}] la red solo contiene ids UN/BI y no hay fuente "
                    f"de dirección (rhea2uniprot ausente) -> red químicamente "
                    f"cerrada; brechas/sinks ~ 0. Descarga rhea2uniprot.tsv.")
    log.info(f"[{genome}] Grafo: {G.number_of_nodes()} metabolitos, "
             f"{G.number_of_edges()} conexiones"
             + (f" ({len(currency_dropped)} metabolitos de moneda excluidos)"
                if currency_dropped else ""))
    log.info(f"[{genome}] Cableado: {n_both} bidireccionales / {n_one} un "
             f"sentido / {n_rev} invertidos por voto UniProt "
             f"(modo {ctx['mode_name']})")

    # ---- Clasificación de metabolitos del genoma ----
    sources, gaps, sinks = [], [], []
    for node in G.nodes():
        inn, out = G.in_degree(node), G.out_degree(node)
        if out > 0 and inn == 0:
            (gaps if node in ctx["producers_of"] else sources).append(node)
        elif inn > 0 and out == 0:
            sinks.append(node)

    genome_families = {family_of.get(r, r) for r in matched}

    def nm(ch):
        return names.get(ch, ch)

    def gen_flag(ch):
        return bool(GENERIC_NAME_RE.match(nm(ch)))

    fam_prod_of = {c: family_list(v, family_of)
                   for c, v in ctx["producers_of"].items()}
    fam_cons_of = {c: family_list(v, family_of)
                   for c, v in ctx["consumers_of"].items()}

    # tría de brechas (v3.3): cubeta accionable por brecha
    gen_instances = {i for lst in ctx["gen_map"].values() for i in lst}
    triage_of = {c: classify_gap(c, fam_prod_of.get(c, []), genome_families,
                                 family_of, reactions, ctx["global_votes"],
                                 names, gen_flag(c), gen_instances)
                 for c in gaps}
    tri_counts = defaultdict(int)
    for c in gaps:
        tri_counts[triage_of[c]] += 1

    # brecha aguas abajo: un sink que OTRA familia consume fisiológicamente.
    # La reversa de la propia familia que lo produjo no cuenta (v3.2: glycine
    # betaine ya no sale como 'consumible' por betA/betB al revés).
    downstream_gaps = []
    for node in sinks:
        prod_fams = set()
        for _, attrs in G.pred[node].items():
            m = re.match(r"RHEA:(\d+)", (attrs or {}).get("enzyme", ""))
            if m:
                prod_fams.add(family_of.get(int(m.group(1)), int(m.group(1))))
        if any(f not in prod_fams for f in fam_cons_of.get(node, [])):
            downstream_gaps.append(node)

    log.info(f"\n[{genome}] === SELF-DISCOVERY: Perfil Metabólico ===")

    log.info(f"\n📥 Fuentes (consumidos sin productor conocido en Rhea -> "
             f"importables): {len(sources)}")
    for s in sorted(sources, key=lambda x: G.out_degree(x),
                    reverse=True)[:args.top]:
        log.info(f"  📥 {nm(s)} ({s})")

    log.info(f"\n🧩 Brechas reales (consumidos sin productor en el genoma, pero "
             f"SÍ con productor en Rhea): {len(gaps)}")
    for gnode in sorted(gaps, key=lambda x: G.out_degree(x),
                        reverse=True)[:args.top]:
        fam_prod = family_list(ctx["producers_of"].get(gnode, []), family_of)
        closable = [r for r in fam_prod if family_of.get(r, r) in genome_families]
        extra = ""
        if gen_flag(gnode):
            extra += " [CLASE GENÉRICA]"
        if closable:
            extra += f" — cerrable con enzima YA en el genoma ({len(closable)})"
        else:
            tri = triage_of.get(gnode, "genuine_missing_enzyme")
            extra += f" — tría: {TRIAGE_LABELS[tri]}"
        log.info(f"  ❌ {nm(gnode)} ({gnode}){extra} — falta enzima; en Rhea lo "
                 f"producen ({len(fam_prod)} familias):")
        for rid in fam_prod[:3]:
            eq = equations.get(rid, "")
            note = direction_note(rid, gnode, reactions, ctx["global_votes"],
                                  as_producer=True)
            log.info(f"      · RHEA:{rid} [{reactions[rid]['type']}] {eq}{note}")

    log.info(
        f"\n🧮 Tría de las {len(gaps)} brechas: "
        f"{tri_counts['genuine_missing_enzyme']} genuinas (sin enzima "
        f"conocida), {tri_counts['generic_producer']} vía clase genérica de "
        f"productor (inferencia is-a de ChEBI pendiente), "
        f"{tri_counts['generic_class']} clases genéricas (artefacto probable), "
        f"{tri_counts['closable_with_genome_enzyme']} closables (enzima ya en "
        f"la red)")

    log.info(f"\n📤 Productos finales (sinks: producidos y nunca consumidos): "
             f"{len(sinks)}")
    for s in sorted(sinks, key=lambda x: G.in_degree(x), reverse=True)[:args.top]:
        gen = " [CLASE GENÉRICA]" if gen_flag(s) else ""
        log.info(f"  ✅ {nm(s)} ({s}){gen}")

    log.info(f"\n⚠️ Sinks con consumidores conocidos en Rhea (posible brecha aguas "
             f"abajo): {len(downstream_gaps)}")
    for s in sorted(downstream_gaps, key=lambda x: G.in_degree(x),
                    reverse=True)[:args.top]:
        parts = []
        for r in family_list(ctx["consumers_of"].get(s, []), family_of)[:3]:
            eq = equations.get(r, f"RHEA:{r}")
            shown = f"{eq[:70]}..." if len(eq) > 73 else eq
            note = direction_note(r, s, reactions, ctx["global_votes"],
                                  as_producer=False)
            parts.append(shown + note)
        log.info(f"  ⚠️ {nm(s)} ({s}) — en Rhea lo consumen: {'; '.join(parts)}")

    if not args.no_report:
        report = {
            "genome": genome,
            "script_version": SCRIPT_VERSION,
            "generated": datetime.now().isoformat(timespec="seconds"),
            "stats": {
                "reactions_in_genome": len(genome_rheas),
                "reactions_matched": len(matched),
                "reactions_obsolete": [f"RHEA:{r}" for r in missing_in_rhea],
                "genome_reaction_types": dict(sorted(type_counts.items())),
                "metabolites": G.number_of_nodes(),
                "edges": G.number_of_edges(),
                "currency_filtered": sorted(currency_dropped),
                "currency_filter_disabled": bool(args.keep_currency),
                "reversible_mode": ctx["mode_name"],
                "wiring": {"bidirectional": n_both, "one_way": n_one,
                           "reversed_by_vote": n_rev},
                "uniprot_ids_in_network": len(up_ids),
                "uniprot_votes_organism": len(org_votes),
                "uniprot_votes_global": len(ctx["global_votes"]),
                "master_ids_from_cache_or_rdf": bool(ctx["master_of"]),
                "generic_instances_resolved": sorted(gen_stats)
                    if ctx["wiring_gen_map"] else [],
            },
            "sources_importable": [{"chebi": c, "name": nm(c)}
                                   for c in sorted(sources)],
            "gaps_missing_enzyme": [
                {"chebi": c, "name": nm(c),
                 "generic": gen_flag(c),
                 "triage": triage_of.get(c, "genuine_missing_enzyme"),
                 "closable_in_genome": [f"RHEA:{r}" for r in
                                        [x for x in fam_prod_of.get(c, [])
                                         if family_of.get(x, x)
                                         in genome_families]],
                 "producers_total": len(fam_prod_of.get(c, [])),
                 "producing_rhea": [f"RHEA:{r}" for r in
                                    fam_prod_of.get(c, [])[:20]]}
                for c in sorted(gaps)],
            "gaps_by_class": {k: tri_counts[k] for k in TRIAGE_LABELS},
            "sinks_final_products": [{"chebi": c, "name": nm(c),
                                      "generic": gen_flag(c)}
                                     for c in sorted(sinks)],
            "downstream_gaps": [{"chebi": c, "name": nm(c),
                                 "generic": gen_flag(c),
                                 "consumers_total": len(fam_cons_of.get(c, [])),
                                 "consuming_rhea": [f"RHEA:{r}" for r in
                                                    fam_cons_of.get(c, [])[:20]]}
                                for c in sorted(downstream_gaps)],
        }
        out_path = net_dir / f"{genome}_network_trace.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
        log.info(f"\n[{genome}] Reporte escrito: {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Auto-descubrimiento de rutas y detección de brechas "
                    "metabólicas (v3.3, voto UniProt + tría de brechas).")
    ap.add_argument("--genome", nargs="?", default=None,
                    help="Ej. GCF_000006985.1_ASM698v1; si se omite, procesa "
                         "TODOS los *_network.json de --net-dir")
    ap.add_argument("--top", type=int, default=10,
                    help="Cuántos ejemplos mostrar por categoría")
    ap.add_argument("--keep-currency", action="store_true",
                    help="No filtrar ATP/H2O/NAD+... del grafo")
    ap.add_argument("--no-reversible", action="store_true",
                    help="Comportamiento v2: todo en un solo sentido "
                         "(canonical L->R)")
    ap.add_argument("--reversible", action="store_true",
                    help="(obsoleto en v3: se ignora)")
    ap.add_argument("--all-reversible", action="store_true",
                    help="Cablear TODAS las reacciones en ambos sentidos "
                         "(diagnóstico: debe dar 0 brechas en un genoma completo)")
    ap.add_argument("--resolve-generics", action="store_true",
                    help="Expandir clases genéricas de nucleótidos "
                         "(a nucleoside diphosphate -> ADP/UDP/...)")
    ap.add_argument("--force-parse", action="store_true",
                    help="Re-parsear rhea.rdf.gz aunque exista cache")
    ap.add_argument("--no-report", action="store_true",
                    help="No escribir {genome}_network_trace.json")
    ap.add_argument("--rhea-dir", type=Path, default=None,
                    help="Directorio TSV de Rhea (por defecto data/metabolism/tsv)")
    ap.add_argument("--net-dir", type=Path, default=None,
                    help="Directorio de redes (por defecto data/metabolism)")
    args = ap.parse_args()

    net_dir = args.net_dir or NET_DIR
    rhea_dir = args.rhea_dir or RHEA_DIR
    reversible_mode = not args.no_reversible

    # ---- fase 1: datos de Rhea (una sola vez, compartidos por todos los genomas)
    names = load_chebi_names(rhea_dir)
    reactions, equations, master_of = load_reaction_data(rhea_dir,
                                                         force=args.force_parse)
    if not reactions:
        return
    family_of, members = compute_families(reactions, master_of)

    # gen_map se calcula SIEMPRE (es barato): la tría de brechas lo usa para
    # detectar 'productor genérico' aunque el cableado no expanda clases.
    gen_map = resolve_nucleotide_generics(reactions, names)
    wiring_gen_map = gen_map if args.resolve_generics else {}
    if wiring_gen_map:
        log.info(f"Clases genéricas de nucleótidos resolubles: {len(gen_map)} "
                 f"({sum(len(v) for v in gen_map.values())} instancias ChEBI)")

    up_file = find_uniprot_file(rhea_dir)
    global_votes = {}
    if up_file:
        global_votes = read_uniprot_votes(up_file)
        n_lr = sum(1 for v in global_votes.values() if v == {"LR"})
        n_rl = sum(1 for v in global_votes.values() if v == {"RL"})
        n_mix = sum(1 for v in global_votes.values() if v == {"LR", "RL"})
        log.info(f"Voto de dirección UniProt ({up_file.name}): "
                 f"{n_lr} familias solo LR, {n_rl} solo RL, {n_mix} mixtas "
                 f"(resto reversible/sin dato)")
    else:
        log.warning("rhea2uniprot(.tsv/.gz) no encontrado en "
                    f"{rhea_dir}: sin voto de dirección; UN/BI se cablean en "
                    "ambos sentidos. Descárgalo con:")
        log.warning(f"    curl -o {rhea_dir / 'rhea2uniprot.tsv'} "
                    "https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv")

    # ---- índice Rhea-wide de productores/consumidores ----
    # Fisiológico en el modo v3.1 (cada familia solo aporta sus direcciones
    # válidas según el voto UniProt); en bruto en los modos v2 y
    # --all-reversible, como en las versiones anteriores.
    if reversible_mode and not args.all_reversible:
        producers_of, consumers_of = build_physiological_index(
            reactions, family_of, global_votes, wiring_gen_map)
        log.info("Índice Rhea-wide: direcciones fisiológicas por familia "
                 "(voto UniProt)")
    else:
        producers_of, consumers_of = build_raw_index(reactions, wiring_gen_map)
        log.info(f"Índice Rhea-wide: todas las variantes "
                 f"({'v2' if not reversible_mode else 'all-reversible'})")

    ctx = {
        "reactions": reactions,
        "equations": equations,
        "names": names,
        "family_of": family_of,
        "master_of": master_of,
        "producers_of": producers_of,
        "consumers_of": consumers_of,
        "gen_map": gen_map,
        "wiring_gen_map": wiring_gen_map,
        "global_votes": global_votes,
        "up_file": up_file,
        "reversible_mode": reversible_mode,
        "mode_name": ("off (v2, una dirección)" if not reversible_mode else
                      ("all-reversible" if args.all_reversible else
                       "v3.1 (voto UniProt)")),
    }

    # ---- fase 2: una o todas las redes
    if args.genome:
        genomes = [args.genome]
    else:
        genomes = sorted(p.name[:-len("_network.json")]
                         for p in net_dir.glob("*_network.json"))
        if not genomes:
            log.error(f"No se encontraron redes *_network.json en {net_dir}")
            return
        log.info(f"Modo batch: {len(genomes)} redes -> {', '.join(genomes)}")
    for genome in genomes:
        process_genome(genome, net_dir, ctx, args)


if __name__ == "__main__":
    main()
