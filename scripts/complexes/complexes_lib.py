#!/usr/bin/env python3
"""complexes_lib.py — núcleo del buscador de complejos proteicos.

Idea
----
El código genético da una cadena de aminoácidos, no la estructura cuaternaria.
Este módulo **no predice interfaces**: ordena hipótesis de "esta proteína va con
esta otra" combinando evidencia genómica y funcional, y las agrupa en
candidatos a complejo. Cada arista declara de dónde sale.

Señales por par
---------------
- genómicas: mismo operón, adyacencia, distancia, misma EC, sintenia conservada;
- STRING: experimentos, vecindario, fusión de dominios, coocurrencia,
  coexpresión y (opcional) bases de datos curadas.

La vía `database` de STRING incluye complejos curados, así que **se excluye por
defecto** al ajustar el modelo: si se usa, la validación contra ComplexPortal
sería parcialmente circular. `--include-string-db` permite medir ese techo.

Contrato
--------
- determinista: sin aleatoriedad en la predicción (la semilla solo se usa en CV);
- sin red: la descarga vive en `fetch_gold_k12.py`;
- solo escribe bajo `data/complexes/`.
"""
from __future__ import annotations

import csv
import gzip
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # dna_poc/
DATA = ROOT / "data"
CX = DATA / "complexes"

# Orden fijo de los canales de STRING en el fichero protein.links.full.
# El header trae "cooccurence" (sic).
CHANNEL_ORDER = ["experiments", "neighborhood", "fusion",
                 "cooccurence", "coexpression", "database"]
CHANNEL_FEATURE = {
    "experiments": "string_exp",
    "neighborhood": "string_nbr",
    "fusion": "string_fus",
    "cooccurence": "string_cooc",
    "coexpression": "string_coex",
    "database": "string_db",
}
FEATURES_GENOMIC = ["same_operon", "adjacent", "log_gap", "log_dist",
                    "synteny", "same_ec"]
FEATURES_STRING = [CHANNEL_FEATURE[c] for c in CHANNEL_ORDER]
FEATURES_ALL = FEATURES_GENOMIC + FEATURES_STRING

# Evidencia legible para el TSV de salida.
CHANNEL_LABEL = {
    "string_exp": "experimentos",
    "string_nbr": "vecindario",
    "string_fus": "fusión",
    "string_cooc": "coocurrencia",
    "string_coex": "coexpresión",
    "string_db": "base de datos",
}


# --------------------------------------------------------------------------- #
# genes
# --------------------------------------------------------------------------- #
@dataclass
class Gene:
    locus_tag: str
    gene_name: str = ""
    ec: str = ""
    contig: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    uniprot: str = ""


def _strand(raw: str) -> str:
    raw = (raw or "").strip()
    if raw in ("-", "-1", "reverse", "R", "complement"):
        return "-"
    return "+"


def _norm_ec(raw: str) -> str:
    """ECs separadas por ';' o '|', ordenadas y sin repetir."""
    parts = re.split(r"[;|,]", raw or "")
    seen: list[str] = []
    for p in parts:
        p = p.strip()
        if p and p not in seen:
            seen.append(p)
    return ";".join(seen)


def read_genes_tsv(path: Path) -> dict[str, Gene]:
    """Lee el `genes.tsv` del pipeline (10_build_gene_table.py)."""
    out: dict[str, Gene] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tag = (r.get("locus_tag") or r.get("gene_id") or "").strip()
            if not tag:
                continue
            try:
                start = int(float(r.get("start") or 0))
                end = int(float(r.get("end") or 0))
            except ValueError:
                continue
            out[tag] = Gene(
                locus_tag=tag,
                gene_name=(r.get("gene_name") or r.get("gene") or "").strip(),
                ec=_norm_ec(r.get("ec") or ""),
                contig=(r.get("contig") or "").strip(),
                start=start, end=end,
                strand=_strand(r.get("strand", "1")),
                uniprot=(r.get("uniprot") or "").strip(),
            )
    return out


def read_genes_gbff(path: Path) -> dict[str, Gene]:
    """CDS de un GenBank con locus_tag, nombre, EC y UniProtKB/Swiss-Prot.

    No usa Biopython. El mapeo UniProtKB solo se usa para el gold de K-12.
    """
    out: dict[str, Gene] = {}
    accession = ""
    cur: dict | None = None

    def flush() -> None:
        nonlocal cur
        if not cur:
            return
        tag = cur.get("locus_tag") or ""
        loc = cur.get("loc") or ""
        m = re.search(r"(\d+)\.\.(\d+)", loc.replace(" ", ""))
        if tag and m and int(m.group(1)) <= int(m.group(2)):
            g = Gene(
                locus_tag=tag,
                gene_name=cur.get("gene_name", ""),
                ec=_norm_ec(cur.get("ec", "")),
                contig=cur.get("contig", ""),
                start=int(m.group(1)), end=int(m.group(2)),
                strand="-" if "complement" in loc else "+",
                uniprot=cur.get("uniprot", ""),
            )
            # un mismo b-number no debería repetirse; si pasa, gana el primero
            out.setdefault(tag, g)
        cur = None

    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("ACCESSION"):
            parts = line.split()
            if len(parts) > 1:
                accession = parts[1]
        if line.startswith("//"):
            flush()
            continue
        if line.startswith("     ") and len(line) > 5 and line[5] != " ":
            flush()
            key = line[5:21].strip()
            if key == "CDS":
                cur = {"locus_tag": "", "gene_name": "", "ec": "",
                       "uniprot": "", "contig": accession,
                       "loc": line[21:].strip(), "in_loc": True}
            continue
        if cur is None:
            continue
        body = line[21:].strip()
        if body.startswith("/"):
            cur["in_loc"] = False
            q = body[1:]
            key, _, val = q.partition("=")
            val = val.strip().strip('"')
            if key == "locus_tag":
                cur["locus_tag"] = val
            elif key == "gene" and not cur["gene_name"]:
                cur["gene_name"] = val
            elif key == "EC_number":
                cur["ec"] = (cur["ec"] + ";" + val) if cur["ec"] else val
            elif key == "db_xref" and val.startswith("UniProtKB/Swiss-Prot:"):
                cur["uniprot"] = val.split(":")[-1]
        elif body and cur.get("in_loc"):
            # solo la location; las continuaciones de /product o /note no van aquí
            cur["loc"] += body
    flush()
    return out


def write_genes_tsv(genes: dict[str, Gene], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["contig", "start", "end", "strand", "locus_tag",
            "gene_name", "ec", "uniprot"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for g in sorted(genes.values(), key=lambda x: (x.contig, x.start)):
            w.writerow([g.contig, g.start, g.end, g.strand, g.locus_tag,
                        g.gene_name, g.ec, g.uniprot])


# --------------------------------------------------------------------------- #
# operones, sintenia y STRING
# --------------------------------------------------------------------------- #
def read_operons(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Devuelve (gen -> operón, operón -> genes)."""
    op_of: dict[str, str] = {}
    operons: dict[str, list[str]] = {}
    if not path or not Path(path).is_file():
        return op_of, operons
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            oid = r.get("operon_id") or ""
            tags = [g for g in (r.get("genes") or "").split(",") if g]
            if not oid or not tags:
                continue
            operons[oid] = tags
            for g in tags:
                op_of[g] = oid
    return op_of, operons


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_string_full(path: Path, species: str = "511145") -> dict[tuple[str, str], tuple[float, ...]]:
    """protein.links.full: 6 canales + combined, normalizados a [0,1].

    Devuelve {(gen_a, gen_b) ordenado: (exp, nbr, fus, cooc, coex, db, combined)}.
    """
    out: dict[tuple[str, str], tuple[float, ...]] = {}
    if not path or not Path(path).is_file():
        return out
    with _open_text(path) as fh:
        header = fh.readline().rstrip("\n").split()
        if not header:
            return out
        idx = {h: i for i, h in enumerate(header)}
        # el header de STRING trae "cooccurence" (sic); toleramos ambas
        cooc = next((h for h in ("cooccurence", "cooccurrence") if h in idx), None)
        if cooc is None or "database" not in idx:
            return out
        cols = [idx["experiments"], idx["neighborhood"], idx["fusion"],
                idx[cooc], idx["coexpression"], idx["database"],
                idx["combined_score"]]
        prefix = f"{species}."
        for line in fh:
            p = line.split()
            if len(p) < len(header):
                continue
            a = p[0][len(prefix):] if p[0].startswith(prefix) else p[0]
            b = p[1][len(prefix):] if p[1].startswith(prefix) else p[1]
            key = (a, b) if a < b else (b, a)
            out[key] = tuple(int(float(p[i])) / 1000.0 for i in cols)
    return out


def read_string_links(path: Path, species: str = "511145") -> dict[tuple[str, str], tuple[float, ...]]:
    """protein.links (solo combined). Rellena los canales con 0."""
    out: dict[tuple[str, str], tuple[float, ...]] = {}
    if not path or not Path(path).is_file():
        return out
    with _open_text(path) as fh:
        header = fh.readline().split()
        if not header:
            return out
        col = header.index("combined_score")
        prefix = f"{species}."
        for line in fh:
            p = line.split()
            if len(p) <= col:
                continue
            a = p[0][len(prefix):] if p[0].startswith(prefix) else p[0]
            b = p[1][len(prefix):] if p[1].startswith(prefix) else p[1]
            key = (a, b) if a < b else (b, a)
            v = int(float(p[col])) / 1000.0
            out[key] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, v)
    return out


def read_synteny(path: Path) -> dict[tuple[str, str], float]:
    """`k12_synteny.tsv`: fracción de genomas que conservan la adyacencia."""
    out: dict[tuple[str, str], float] = {}
    if not path or not Path(path).is_file():
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a, b = r.get("left"), r.get("right")
            if not a or not b:
                continue
            try:
                out[(a, b) if a < b else (b, a)] = float(r.get("frac") or 0)
            except ValueError:
                continue
    return out


# --------------------------------------------------------------------------- #
# ComplexPortal -> gold
# --------------------------------------------------------------------------- #
_PART = re.compile(r"(?:^|\|)([A-Z][A-Z0-9]{4,9})\((\d+)\)")

# Clases de evidencia de ComplexPortal (columna "Evidence Code").
_EVIDENCE_CLASS = {
    "ECO:0000353": "experimental",       # interacción física, aserción manual
    "ECO:0005544": "inferida",           # por ortología
    "ECO:0005546": "inferida",           # reconstrucción por parálogos
    "ECO:0005547": "inferida",           # reconstrucción por conocimiento previo
}


def classify_evidence(code: str) -> str:
    """Traduce el ECO de ComplexPortal a una clase legible.

    No todo lo "documentado" es experimental: ComplexPortal también cura
    complejos por ortología, parálogos o conocimiento previo.
    """
    for eco, cls in _EVIDENCE_CLASS.items():
        if eco in (code or ""):
            return cls
    return "otra"


def read_complexportal(path: Path, acc2loc: dict[str, str],
                       drop_incomplete: bool = True) -> tuple[list[dict], dict]:
    """Complejos curados de ComplexPortal mapeados a locus_tag.

    `drop_incomplete=True` descarta complejos con subunidades sin mapear: un
    subconjunto parcial inventaría pares que no son el complejo completo.
    """
    stats = {"total": 0, "unmapped": 0, "homomers": 0, "kept": 0}
    out: list[dict] = []
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return out, stats
    for r in rows[1:]:
        if len(r) < 12:
            continue
        stats["total"] += 1
        members = _PART.findall(r[4])
        if not members:
            continue
        accs = {a for a, _ in members}
        missing = {a for a in accs if a not in acc2loc}
        locs = sorted({acc2loc[a] for a in accs if a in acc2loc})
        if len(locs) < 2:
            stats["homomers"] += 1
            continue
        if drop_incomplete and missing:
            stats["unmapped"] += 1
            continue
        out.append({
            "complex_id": r[0],
            "complex_name": r[1],
            "subunits": len(members),
            "locus_tags": locs,
            "properties": r[11],
            "evidence_code": r[5],
            "evidence_class": classify_evidence(r[5]),
            "source": "ComplexPortal",
        })
    stats["kept"] = len(out)
    return out, stats


# --------------------------------------------------------------------------- #
# contexto y features
# --------------------------------------------------------------------------- #
@dataclass
class Context:
    genes: dict[str, Gene]
    order: list[str] = field(default_factory=list)
    pos: dict[str, int] = field(default_factory=dict)
    op_of: dict[str, str] = field(default_factory=dict)
    synteny: dict[tuple[str, str], float] = field(default_factory=dict)
    string: dict[tuple[str, str], tuple[float, ...]] = field(default_factory=dict)


def pair_key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a < b else (b, a)


def build_context(genes: dict[str, Gene], op_of: dict[str, str] | None = None,
                  synteny: dict[tuple[str, str], float] | None = None,
                  string: dict[tuple[str, str], tuple[float, ...]] | None = None) -> Context:
    order = sorted(genes, key=lambda t: (genes[t].contig, genes[t].start, genes[t].end))
    return Context(genes=genes, order=order,
                   pos={g: i for i, g in enumerate(order)},
                   op_of=dict(op_of or {}),
                   synteny=dict(synteny or {}),
                   string=dict(string or {}))


def pair_geometry(a: str, b: str, ctx: Context) -> dict[str, float | int]:
    """Gap intergénico y distancia entre inicios, con `a` aguas arriba."""
    ga, gb = ctx.genes[a], ctx.genes[b]
    if (ga.contig, ga.start) > (gb.contig, gb.start):
        ga, gb = gb, ga
    same_contig = ga.contig == gb.contig and ga.contig != ""
    gap = gb.start - ga.end - 1 if same_contig else 0
    adjacent = 0
    if same_contig and ga.strand == gb.strand:
        ia, ib = ctx.pos.get(ga.locus_tag, -1), ctx.pos.get(gb.locus_tag, -1)
        adjacent = 1 if abs(ib - ia) == 1 else 0
    return {"gap": int(gap), "dist": int(abs(ga.start - gb.start)),
            "adjacent": adjacent, "upstream": ga.locus_tag,
            "downstream": gb.locus_tag}


def features_row(a: str, b: str, ctx: Context) -> dict[str, float]:
    """Todas las señales de un par. `a` y `b` se ordenan por coordenada."""
    ga, gb = ctx.genes[a], ctx.genes[b]
    if (ga.contig, ga.start) > (gb.contig, gb.start):
        ga, gb = gb, ga
    geo = pair_geometry(a, b, ctx)
    same_op = 1.0 if (ctx.op_of.get(ga.locus_tag)
                      and ctx.op_of.get(ga.locus_tag) == ctx.op_of.get(gb.locus_tag)) else 0.0
    eca = set((ga.ec or "").split(";")) - {""}
    ecb = set((gb.ec or "").split(";")) - {""}
    same_ec = 1.0 if eca & ecb else 0.0
    key = pair_key(a, b)
    ch = ctx.string.get(key, (0.0,) * 7)
    row = {
        "same_operon": same_op,
        "adjacent": float(geo["adjacent"]),
        "log_gap": math.log1p(max(0, geo["gap"])),
        "log_dist": math.log1p(geo["dist"]),
        "synteny": float(ctx.synteny.get(key, 0.0)),
        "same_ec": same_ec,
    }
    for i, channel in enumerate(CHANNEL_ORDER):
        row[CHANNEL_FEATURE[channel]] = ch[i]
    return row


def feature_vector(a: str, b: str, ctx: Context, names: list[str]) -> list[float]:
    row = features_row(a, b, ctx)
    return [row[n] for n in names]


def candidate_pairs(ctx: Context, min_string: float = 0.2,
                    min_combined: float = 0.3) -> list[tuple[str, str]]:
    """Pares que merecen puntuarse. Evita O(n²) en genomas grandes.

    Incluye: todos los pares intra-operón, genes adyacentes, aristas de STRING
    con señal (sin contar `database`) y adyacencias conservadas.
    """
    out: set[tuple[str, str]] = set()

    by_op: dict[str, list[str]] = defaultdict(list)
    for g, op in ctx.op_of.items():
        if g in ctx.genes:
            by_op[op].append(g)
    for tags in by_op.values():
        tags = sorted(tags)
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                out.add(pair_key(tags[i], tags[j]))

    for a, b in zip(ctx.order, ctx.order[1:]):
        ga, gb = ctx.genes[a], ctx.genes[b]
        if ga.contig == gb.contig and ga.strand == gb.strand:
            out.add(pair_key(a, b))

    for key, ch in ctx.string.items():
        if key[0] not in ctx.genes or key[1] not in ctx.genes:
            continue
        if max(ch[:5]) >= min_string or ch[6] >= min_combined:
            out.add(key)

    for key, frac in ctx.synteny.items():
        if frac > 0 and key[0] in ctx.genes and key[1] in ctx.genes:
            out.add(key)
    return sorted(out)


# --------------------------------------------------------------------------- #
# modelo
# --------------------------------------------------------------------------- #
def _sigmoid(z: float) -> float:
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)


def fit_model(X: list[list[float]], y: list[int], names: list[str],
              C: float = 1.0, seed: int = 0) -> dict:
    """Regresión logística balanceada con estandarización guardada."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression

    arr = np.asarray(X, dtype=float)
    yv = np.asarray(y, dtype=int)
    mean = arr.mean(axis=0)
    scale = arr.std(axis=0)
    scale[scale == 0] = 1.0
    clf = LogisticRegression(max_iter=3000, class_weight="balanced",
                             C=C, random_state=seed)
    clf.fit((arr - mean) / scale, yv)
    return {
        "version": 1,
        "features": list(names),
        "mean": [float(v) for v in mean],
        "scale": [float(v) for v in scale],
        "coef": [float(v) for v in clf.coef_[0]],
        "intercept": float(clf.intercept_[0]),
        "n_pos": int(yv.sum()),
        "n_neg": int(len(yv) - yv.sum()),
    }


def predict(model: dict, x: list[float]) -> float:
    z = model["intercept"]
    for i, v in enumerate(x):
        z += model["coef"][i] * (v - model["mean"][i]) / model["scale"][i]
    return _sigmoid(z)


def predict_many(model: dict, X: list[list[float]]) -> list[float]:
    return [predict(model, x) for x in X]


def save_model(model: dict, path: Path) -> None:
    import json
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(model, indent=2, sort_keys=True) + "\n")


def load_model(path: Path) -> dict:
    import json
    return json.loads(Path(path).read_text())


# --------------------------------------------------------------------------- #
# clustering
# --------------------------------------------------------------------------- #
class _UF:
    def __init__(self):
        self.parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def cluster_edges(edges: dict[tuple[str, str], float], top_k: int = 2,
                  min_p: float = 0.1) -> list[list[str]]:
    """Componentes de la gráfica de **top-k mutuo**.

    Sin el filtro mutuo, un umbral bajo produce una componente gigante por
    encadenamiento (en K-12: 748 genes → 14 clusters, el mayor con cientos).
    Con top-2 mutuo: ~210 clusters y Jaccard medio 0.68 contra ComplexPortal.
    """
    edges = {k: v for k, v in edges.items() if v >= min_p}
    top: dict[str, list[tuple[float, str]]] = defaultdict(list)
    for (a, b), p in edges.items():
        top[a].append((p, b))
        top[b].append((p, a))
    topset = {g: {b for _, b in sorted(lst, reverse=True)[:top_k]}
              for g, lst in top.items()}
    uf = _UF()
    for (a, b), _ in edges.items():
        if b in topset.get(a, ()) and a in topset.get(b, ()):
            uf.union(a, b)
    comps: dict[str, list[str]] = defaultdict(list)
    for g in list(uf.parent):
        comps[uf.find(g)].append(g)
    return sorted((sorted(c) for c in comps.values() if len(c) > 1),
                  key=lambda c: (-len(c), c[0]))


def evidence_tags(a: str, b: str, ctx: Context, min_string: float = 0.2) -> list[str]:
    row = features_row(a, b, ctx)
    tags: list[str] = []
    if row["same_operon"]:
        tags.append("operón")
    if row["adjacent"]:
        tags.append("adyacente")
    if row["synteny"] > 0:
        tags.append(f"sintenia:{row['synteny']:.2f}")
    if row["same_ec"]:
        tags.append("misma EC")
    for name in FEATURES_STRING:
        if row[name] >= min_string:
            tags.append(f"STRING:{CHANNEL_LABEL[name]}:{row[name]:.2f}")
    return tags


def confidence(n_genes: int, min_p: float, mean_p: float) -> str:
    """Tres niveles; los cortes se calibraron en `validate_k12.py`."""
    if min_p >= 0.90 and n_genes <= 4:
        return "alta"
    if mean_p >= 0.70:
        return "media"
    return "baja"


def cluster_summary(cluster: list[str], edges: dict[tuple[str, str], float],
                    ctx: Context) -> dict:
    ps = [edges.get(pair_key(a, b), 0.0) for i, a in enumerate(cluster)
          for b in cluster[i + 1:]]
    ps = ps or [0.0]
    tags: list[str] = []
    seen: set[str] = set()
    for i, a in enumerate(cluster):
        for b in cluster[i + 1:]:
            for t in evidence_tags(a, b, ctx):
                key = t.split(":")[0]
                if key not in seen:
                    seen.add(key)
                    tags.append(t)
    return {
        "n_genes": len(cluster),
        "mean_p": round(sum(ps) / len(ps), 4),
        "min_p": round(min(ps), 4),
        "evidence": ";".join(tags),
        "confidence": confidence(len(cluster), min(ps), sum(ps) / len(ps)),
    }
