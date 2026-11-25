#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Agrupa los factores sigma de los 400 genomas por ortólogo.

Lee `data/sigma/<genoma>/sigma_factors.tsv`, extrae las proteínas de
`data/annotations_ab_initio/<genoma>/proteins.faa` (o de
`data/regulation/<genoma>/proteins.faa`), corre DIAMOND all-vs-all y agrupa
por identidad/cobertura. Añade a cada factor:

    countable       yes/no
    quality_flag    validated | consensus | partial_domain | no_product | fragment
    ortholog_group  OG0001...
    copy_number     copias en el genoma
    is_expansion    yes/no

También escribe `data/sigma/sigma_orthologs.tsv` y actualiza
`data/sigma/summary_all.tsv`.

Uso:
    python scripts/sigma/20_group_sigma.py --all
    python scripts/sigma/20_group_sigma.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import (  # noqa: E402
    OUTPUT_FACTOR_COLS, _read_tsv, _sha256, _write_tsv, factor_quality,
)

SEP = "@@"
FACTOR_COLS_EXTRA = ["countable", "quality_flag", "ortholog_group",
                     "copy_number", "is_expansion"]
SUMMARY_EXTRA = ["n_factors_countable", "n_factors_unresolved",
                 "n_ortholog_groups", "n_ortholog_groups_countable",
                 "n_expansions"]


# --------------------------------------------------------------------------
# Proteínas
# --------------------------------------------------------------------------
def _load_faa(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def _protein_paths(data_root: Path, genome: str) -> list[Path]:
    return [
        data_root / "annotations_ab_initio" / genome / "proteins.faa",
        data_root / "regulation" / genome / "proteins.faa",
        data_root / "annotations" / f"{genome}.faa.gz",
    ]


def _load_genome_proteins(data_root: Path, genome: str) -> dict[str, str]:
    for path in _protein_paths(data_root, genome):
        if path.exists():
            got = _load_faa(path)
            if got:
                return got
    return {}


# --------------------------------------------------------------------------
# Clustering
# --------------------------------------------------------------------------
def _find(parent: dict[str, str], x: str) -> str:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _union(parent: dict[str, str], a: str, b: str) -> None:
    ra, rb = _find(parent, a), _find(parent, b)
    if ra != rb:
        parent[rb] = ra


def _cluster_from_hits(ids: list[str], hits: list[tuple[str, str, float, float, float]],
                       min_pident: float = 50.0, min_cov: float = 70.0,
                       max_evalue: float = 1e-5) -> dict[str, str]:
    """Union-find sobre hits; devuelve id -> representante de cluster."""
    parent = {i: i for i in ids}
    for q, s, pident, qcov, evalue in hits:
        if q not in parent or s not in parent or q == s:
            continue
        if pident >= min_pident and qcov >= min_cov and evalue <= max_evalue:
            _union(parent, q, s)
    groups: dict[str, list[str]] = defaultdict(list)
    for i in ids:
        groups[_find(parent, i)].append(i)
    # Representante estable: el miembro lexicográficamente menor.
    out = {}
    for members in groups.values():
        rep = sorted(members)[0]
        for m in members:
            out[m] = rep
    return out


def _run_diamond(faa: Path, work: Path, threads: int) -> list[tuple[str, str, float, float, float]]:
    db = work / "sigma_proteins"
    hits = work / "sigma_hits.tsv"
    subprocess.run(["diamond", "makedb", "--in", str(faa), "-d", str(db), "--quiet"],
                   check=True)
    subprocess.run([
        "diamond", "blastp", "-q", str(faa), "-d", str(db), "-o", str(hits),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "evalue",
        "--max-target-seqs", "40", "--threads", str(max(1, threads)),
        "--quiet", "--masking", "0",
    ], check=True)
    out = []
    for line in hits.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        try:
            out.append((f[0], f[1], float(f[2]), float(f[3]), float(f[4])))
        except ValueError:
            continue
    return out


# --------------------------------------------------------------------------
# Anotación de grupos
# --------------------------------------------------------------------------
def _read_factors(path: Path) -> list[dict[str, str]]:
    return _read_tsv(path)


def _factor_key(row: dict[str, str]) -> str:
    return f"{row.get('genome', '')}{SEP}{row.get('gene_id', '')}"


def _quality(row: dict[str, str]) -> tuple[str, str, bool]:
    if (row.get("countable") or "").strip() in {"yes", "no"}:
        return (row.get("quality_flag", ""), row.get("factor_note", ""),
                row["countable"] == "yes")
    return factor_quality(row)


def group_all(data_root: Path, limit: int = 0, threads: int = 8) -> dict:
    sigma_root = data_root / "sigma"
    summaries = sorted(sigma_root.glob("*/sigma_factors.tsv"))
    # excluir el genoma de prueba
    summaries = [p for p in summaries if "test" not in p.parent.name]
    if limit:
        summaries = summaries[:limit]
    if not summaries:
        raise FileNotFoundError("no hay data/sigma/*/sigma_factors.tsv")

    rows_by_genome: dict[str, list[dict[str, str]]] = {}
    seqs: dict[str, str] = {}
    meta: dict[str, dict[str, str]] = {}
    missing_seq = 0
    for path in summaries:
        genome = path.parent.name
        rows = _read_factors(path)
        rows_by_genome[genome] = rows
        proteins = _load_genome_proteins(data_root, genome)
        for row in rows:
            gid = row.get("gene_id", "")
            seq = proteins.get(gid) or proteins.get(row.get("locus_tag", ""))
            if not seq:
                missing_seq += 1
                continue
            key = _factor_key(row)
            seqs[key] = seq
            meta[key] = row

    work = sigma_root / "_work"
    work.mkdir(parents=True, exist_ok=True)
    faa = work / "sigma_proteins.faa"
    with faa.open("w", encoding="utf-8") as fh:
        for key in sorted(seqs):
            fh.write(f">{key}\n{seqs[key]}\n")

    hits = _run_diamond(faa, work, threads) if seqs else []
    clusters = _cluster_from_hits(sorted(seqs), hits)

    # Agrupar por representante y ordenar por tamaño.
    by_rep: dict[str, list[str]] = defaultdict(list)
    for key, rep in clusters.items():
        by_rep[rep].append(key)
    ordered = sorted(by_rep.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    group_of: dict[str, str] = {}
    group_meta: list[dict] = []
    for idx, (rep, members) in enumerate(ordered, 1):
        og = f"OG{idx:04d}"
        for m in members:
            group_of[m] = og
        genomes = {meta[m]["genome"] for m in members if m in meta}
        fams = Counter(meta[m].get("sigma_type", "") for m in members if m in meta)
        countable = 0
        for m in members:
            row = meta.get(m)
            if row is not None:
                _, _, ok = _quality(row)
                countable += int(ok)
        group_meta.append({
            "ortholog_group": og, "representative": rep, "n_members": len(members),
            "n_genomes": len(genomes), "n_copies": len(members),
            "n_countable": countable,
            "families": dict(fams.most_common()),
            "genomes": ",".join(sorted(genomes)),
        })

    # Asignar grupo/copias a cada factor y reescribir.
    for genome, rows in rows_by_genome.items():
        # calidad y contable
        for row in rows:
            qflag, qnote, countable = _quality(row)
            row["countable"] = "yes" if countable else "no"
            row["quality_flag"] = qflag
            if qnote and not row.get("factor_note"):
                row["factor_note"] = qnote
            key = _factor_key(row)
            row["ortholog_group"] = group_of.get(key, "")
            row["copy_number"] = ""
            row["is_expansion"] = ""
        # copias por grupo dentro del genoma
        counts = Counter(r["ortholog_group"] for r in rows if r["ortholog_group"])
        for row in rows:
            og = row["ortholog_group"]
            if og:
                row["copy_number"] = str(counts[og])
                row["is_expansion"] = "yes" if counts[og] > 1 else "no"
        out = sigma_root / genome / "sigma_factors.tsv"
        _write_tsv(out, rows, OUTPUT_FACTOR_COLS)
        # resumen
        sp = sigma_root / genome / "sigma_summary.json"
        try:
            summ = json.loads(sp.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            summ = {"genome": genome}
        summ["n_factors"] = len(rows)
        summ["n_factors_countable"] = sum(1 for r in rows if r["countable"] == "yes")
        summ["n_factors_unresolved"] = sum(1 for r in rows if r.get("sigma_type") == "sigma_unresolved")
        countable_rows = [r for r in rows if r.get("countable") == "yes"]
        c_groups = Counter(r["ortholog_group"] for r in countable_rows if r.get("ortholog_group"))
        summ["n_ortholog_groups"] = len({r["ortholog_group"] for r in rows if r["ortholog_group"]})
        summ["n_ortholog_groups_countable"] = len(c_groups)
        summ["n_expansions"] = len({g for g, n in c_groups.items() if n > 1})
        summ["sha256"] = _sha256([out, sigma_root / genome / "sigma_promoters.tsv"])
        sp.write_text(json.dumps(summ, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")

    # Tabla de grupos
    og_path = sigma_root / "sigma_orthologs.tsv"
    _write_tsv(og_path, group_meta, ["ortholog_group", "representative", "n_members",
                                     "n_genomes", "n_copies", "n_countable",
                                     "families", "genomes"])

    # Actualizar summary_all
    sp = sigma_root / "summary_all.tsv"
    if sp.exists():
        rows = _read_tsv(sp)
        for row in rows:
            try:
                summ = json.loads((sigma_root / row["genome"] / "sigma_summary.json")
                                  .read_text(encoding="utf-8"))
                for col in SUMMARY_EXTRA:
                    row[col] = str(summ.get(col, ""))
            except (OSError, ValueError):
                pass
        cols = list(rows[0].keys()) if rows else []
        for col in SUMMARY_EXTRA:
            if col not in cols:
                cols.append(col)
        _write_tsv(sp, rows, cols)

    return {
        "genomes": len(rows_by_genome), "factors": sum(len(v) for v in rows_by_genome.values()),
        "with_sequence": len(seqs), "missing_sequence": missing_seq,
        "groups": len(group_meta), "hits": len(hits),
        "orthologs": str(og_path),
    }


def selftest() -> str:
    ids = ["a", "b", "c", "d"]
    hits = [("a", "b", 90.0, 95.0, 1e-30), ("c", "d", 30.0, 95.0, 1e-30)]
    cl = _cluster_from_hits(ids, hits)
    assert cl["a"] == cl["b"]
    assert cl["c"] != cl["d"]
    assert factor_quality({"sigma_type": "sigma_unresolved", "source_status": "novel",
                           "product": "", "source_evidence": "domains:Sigma70_r2"})[2] is False
    return "clustering y calidad sigma OK"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0
    if not args.all:
        ap.error("indica --all o --selftest")
    data_root = ROOT / "data"
    t0 = time.time()
    res = group_all(data_root, limit=args.limit, threads=args.threads)
    print(json.dumps(res, ensure_ascii=False, indent=2))
    print(f"tiempo: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
