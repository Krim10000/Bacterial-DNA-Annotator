#!/usr/bin/env python3
"""regulondb_gt.py — verdad de referencia de RegulonDB (E. coli K-12) para L2/L3.

RegulonDB v13+ ya no publica ficheros planos: los datasets se sirven por
GraphQL (`getDataOfFile`). Este script los descarga una vez, los normaliza a
TSVs limpios y deja un manifiesto con versión, fecha, hashes y recuentos.

Por qué existe
--------------
El predictor de operones (L2) y el regulón (L3) se validan contra RegulonDB,
no contra una lista de operones escrita a mano. Cada operón y cada promotor
trae evidencia y PMIDs, y los genes vienen en b-numbers (`b0001`...), que son
exactamente los `locus_tag` de NC_000913.3 — el mismo genoma que
`data/test_ecoli.gbff`. No hace falta mapear identificadores.

Salidas (`data/regulondb/gt/`)
------------------------------
  gt_operons.tsv    operón + coordenadas + hebra + genes + confianza
  gt_tus.tsv        unidad de transcripción (subconjunto del operón, con su
                    propio promotor y factor sigma)
  gt_promoters.tsv  TSS, sigma, cajas -10/-35 con posición y secuencia
  gt_tfs.tsv        factor de transcripción, gen, b-number, efector (metabolito)
  manifest.json     versión, fecha, recuentos, sha256 de cada raw

Detalle importante: un operón puede tener varias TUs (promotores internos). La
comparación correcta para L2 es el **operón** (`gt_operons.tsv`, 2.605 filas),
no la TU. Para L3, la arista TF→TU está en `gt_tus.tsv`.

Uso
---
    python scripts/validation/regulondb_gt.py            # usa raw si existe
    python scripts/validation/regulondb_gt.py --fetch    # fuerza descarga
    python scripts/validation/regulondb_gt.py --offline  # prohíbe red

Cita: Salgado H, Gama-Castro S, et al. RegulonDB, Nucleic Acids Research.
RegulonDB es de uso académico/no comercial.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from collections import defaultdict
from pathlib import Path

URL = "https://regulondb.ccg.unam.mx/graphql"
QUERY = '{getDataOfFile(fileName:"%s"){content rdbVersion creationDate}}'
DATASETS = ("TuPromOperonTFSet", "OperonSet", "PromoterSet", "TFSet",
            "TerminatorSet", "RegulatorSet")
ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "regulondb"
GT_DIR = RAW_DIR / "gt"
STRAND = {"forward": "+", "reverse": "-", "fwd": "+", "rev": "-"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str) -> tuple[str, dict]:
    """Descarga un dataset por GraphQL. Devuelve (contenido, metadatos)."""
    req = urllib.request.Request(
        URL, method="POST",
        data=json.dumps({"query": QUERY % name}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        payload = json.loads(resp.read())
    data = payload["data"]["getDataOfFile"]
    return data["content"], {"rdb_version": data["rdbVersion"],
                             "creation_date": data["creationDate"]}


def ensure_raw(names, force: bool, offline: bool) -> dict:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    meta = {}
    for name in names:
        path = RAW_DIR / f"{name}.tsv"
        if path.exists() and not force:
            meta[name] = {"cached": True}
            continue
        if offline:
            raise SystemExit(f"falta {path} y --offline está activo")
        content, info = fetch(name)
        path.write_text(content)
        meta[name] = {**info, "cached": False}
        print(f"  descargado {name}: {len(content) / 1e6:.2f} MB")
    return meta


def read_tsv(path: Path) -> list[dict]:
    """Lee un TSV de RegulonDB. La cabecera viene como '1)col 2)col ...'."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = [c.split(")")[-1].strip() for c in rows[0]]
    out = []
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        out.append({h: (row[i].strip() if i < len(row) else "")
                    for i, h in enumerate(header)})
    return out


def _genes(field: str) -> list[str]:
    return [g.strip() for g in field.split("|") if g.strip()]


def build_operons(tus: list[dict], operons: list[dict]) -> list[dict]:
    """Une los b-numbers de las TUs de cada operón (unión = genes del operón).

    Un operón con promotores internos aparece partido en varias TUs; la unión
    de sus genes reconstruye el operón completo. Se marca `mapping_ok` para que
    el benchmark pueda excluir los 24 casos de elementos de inserción donde los
    nombres no resuelven (insAB4 vs insA4|insB4).
    """
    union: dict[str, set[str]] = defaultdict(set)
    for tu in tus:
        union[tu["operonId"]].update(_genes(tu["tuGenesBnumber"]))
    out = []
    for op in operons:
        try:
            n = int(op["numberOfGenes"])
        except (ValueError, KeyError):
            n = 0
        genes = union.get(op["operonId"], set())
        left, right = op.get("firstGeneLeftPos", ""), op.get("lastGeneRightPos", "")
        try:
            start, end = sorted((int(left), int(right)))
        except ValueError:
            start = end = 0
        out.append({
            "operon_id": op["operonId"],
            "operon_name": op.get("operonName", ""),
            "strand": STRAND.get(op.get("strand", "").lower(), "?"),
            "start": start, "end": end, "n_genes": n,
            "genes": ",".join(sorted(genes)),
            "gene_names": op.get("operonGenes", ""),
            "n_genes_tus": len(genes),
            "mapping_ok": "yes" if n and len(genes) == n else "no",
            "confidence": op.get("confidenceLevel", ""),
            "evidence": op.get("operonEvidence", ""),
        })
    return out


def build_tus(tus: list[dict]) -> list[dict]:
    return [{
        "tu_id": t["tuId"], "tu_name": t.get("tuName", ""),
        "operon_id": t.get("operonId", ""), "operon_name": t.get("operonName", ""),
        "genes": ",".join(sorted(_genes(t.get("tuGenesBnumber", "")))),
        "first_gene_pos": t.get("firstGenePos", ""),
        "last_gene_pos": t.get("lastGenePos", ""),
        "promoter_pos": t.get("promoterPos", ""),
        "sigma_factor": t.get("sigmaFactor", ""),
        "tfs": t.get("tfsNames", ""),
        "terminator_positions": t.get("terminatorPositions", ""),
    } for t in tus]


def build_promoters(prom: list[dict]) -> list[dict]:
    return [{
        "promoter_id": p["id"], "name": p.get("name", ""),
        "strand": STRAND.get(p.get("strand", "").lower(), "?"),
        "tss": p.get("posTSS", ""), "sigma_factor": p.get("sigmaFactor", ""),
        "first_gene": p.get("firstGeneName", ""),
        "dist_to_first_gene": p.get("distToFirstGene", ""),
        "minus10_pos": p.get("boxMinus10pos", ""),
        "minus10_seq": p.get("boxMinus10seq", ""),
        "minus35_pos": p.get("boxMinus35pos", ""),
        "minus35_seq": p.get("boxMinus35seq", ""),
        "confidence": p.get("confidenceLevel", ""),
    } for p in prom]


def build_tfs(tfs: list[dict]) -> list[dict]:
    return [{
        "tf_id": t["id"], "name": t.get("name", ""),
        "gene": t.get("geneCodingForTF", ""),
        "bnumber": t.get("geneBnumberCodingForTF", ""),
        "active_conformation": t.get("tfActiveConformations", ""),
        "effector": t.get("tfActiveConformationsEffector", ""),
    } for t in tfs]


def build_aliases(tus: list[dict]) -> list[dict]:
    """b-number <-> nombre de gen, alineados por posición dentro de cada TU.

    Es la única fuente fiable del alias: en `OperonSet` los nombres vienen sin
    su b-number, y emparejarlos por conjunto daría asignaciones arbitrarias.
    """
    alias: dict[str, str] = {}
    for tu in tus:
        names = _genes(tu.get("tuGenesNames", ""))
        bnums = _genes(tu.get("tuGenesBnumber", ""))
        for name, bnum in zip(names, bnums):
            alias.setdefault(name, bnum)
    return [{"name": n, "bnumber": b} for n, b in sorted(alias.items())]


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fetch", action="store_true", help="fuerza descarga")
    ap.add_argument("--offline", action="store_true", help="prohíbe red")
    ap.add_argument("--out", type=Path, default=GT_DIR)
    args = ap.parse_args()

    print(f"RegulonDB: {len(DATASETS)} datasets")
    meta = ensure_raw(DATASETS, args.fetch, args.offline)

    raw = {n: read_tsv(RAW_DIR / f"{n}.tsv") for n in DATASETS}
    operons = build_operons(raw["TuPromOperonTFSet"], raw["OperonSet"])
    tus = build_tus(raw["TuPromOperonTFSet"])
    promoters = build_promoters(raw["PromoterSet"])
    tfs = build_tfs(raw["TFSet"])

    write_tsv(args.out / "gt_operons.tsv", operons)
    write_tsv(args.out / "gt_tus.tsv", tus)
    write_tsv(args.out / "gt_promoters.tsv", promoters)
    write_tsv(args.out / "gt_tfs.tsv", tfs)
    write_tsv(args.out / "gt_gene_aliases.tsv", build_aliases(raw["TuPromOperonTFSet"]))

    versions = {m.get("rdb_version") for m in meta.values() if m.get("rdb_version")}
    manifest = {
        "source": URL, "datasets": DATASETS,
        "rdb_version": sorted(versions)[-1] if versions else None,
        "creation_date": next((m.get("creation_date") for m in meta.values()
                               if m.get("creation_date")), None),
        "genome": "Escherichia coli K-12 MG1655 (NC_000913.3)",
        "id_space": "b-numbers (b0001...) == locus_tag de NC_000913.3",
        "counts": {
            "operons": len(operons),
            "operons_multi": sum(1 for o in operons if o["n_genes"] >= 2),
            "operons_mapping_ok": sum(1 for o in operons if o["mapping_ok"] == "yes"),
            "tus": len(tus),
            "tus_with_sigma": sum(1 for t in tus if t["sigma_factor"]),
            "tus_with_tf": sum(1 for t in tus if t["tfs"]),
            "promoters": len(promoters),
            "promoters_with_minus10": sum(1 for p in promoters if p["minus10_seq"]),
            "promoters_with_minus35": sum(1 for p in promoters if p["minus35_seq"]),
            "tfs": len(tfs),
            "terminators": len(raw["TerminatorSet"]),
        },
        "sha256": {n: _sha256(RAW_DIR / f"{n}.tsv") for n in DATASETS},
        "citation": "Salgado H, Gama-Castro S, et al. RegulonDB, Nucleic Acids Res.",
        "license": "académico / no comercial",
    }
    (RAW_DIR / "manifest.json").write_text(json.dumps(manifest, indent=1))

    c = manifest["counts"]
    print(f"\nGT {args.out}")
    print(f"  operones {c['operons']} (>=2 genes: {c['operons_multi']}, "
          f"mapeo ok: {c['operons_mapping_ok']})")
    print(f"  TUs {c['tus']} | con sigma {c['tus_with_sigma']} | con TF {c['tus_with_tf']}")
    print(f"  promotores {c['promoters']} (-10: {c['promoters_with_minus10']}, "
          f"-35: {c['promoters_with_minus35']}) | TFs {c['tfs']}")
    if manifest["rdb_version"]:
        print(f"  RegulonDB {manifest['rdb_version']} ({manifest['creation_date']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
