#!/usr/bin/env python3
"""Experimento para elegir el umbral de deduplicación entre genomas.

Toma proteínas YA anotadas con acceso Swiss-Prot conocido, las agrupa con
`diamond cluster` a varios umbrales de identidad y mide:

  - pureza: ¿el cluster mezcla genes distintos?
  - recall de pares: ¿junta copias del mismo gen?
  - F1 y compresión (nº de clusters).

La verdad de referencia es el nombre de gen (`GN=`) del hit Swiss-Prot; si no
hay, el acceso. Así, dos ortólogos que caen en accesos distintos pero comparten
`GN=` cuentan como el mismo gen.

Uso:
    python3 -m ghost.dedup_eval --data-root <data> --out ./out/dedup \
        --sample 20000 --thresholds 80,90,95,97,100
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import random
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

GEN_RE = re.compile(r"GN=(\S+)")


def collect_annotations(data_root: Path) -> list[dict]:
    """(proteína, acceso, clave de gen, genoma) de functional.tsv."""
    rows: list[dict] = []
    for tsv in sorted(glob.glob(str(data_root / "annotations_ab_initio" / "*" / "functional.tsv"))):
        genome = os.path.basename(os.path.dirname(tsv))
        with open(tsv) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 13:
                    continue
                sseqid = p[1]
                acc = sseqid.split("|")[1] if "|" in sseqid else sseqid
                # Solo hits de longitud completa: evita que un fragmento
                # herede el nombre de una proteína grande.
                try:
                    aln, qlen, slen = int(p[3]), int(p[4]), int(p[5])
                except ValueError:
                    continue
                if aln / max(1, qlen) < 0.7 or aln / max(1, slen) < 0.7:
                    continue
                title = " ".join(p[12:])
                m = GEN_RE.search(title)
                key = m.group(1) if m else acc
                rows.append({"pid": p[0], "acc": acc, "key": key, "genome": genome})
    return rows


def sample_rows(rows: list[dict], n: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    by_key: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_key[r["key"]].append(r)
    # mitad de claves repetidas (garantiza pares), mitad aleatoria
    repeated = [r for k, rs in by_key.items() if len(rs) >= 2 for r in rs]
    rng.shuffle(repeated)
    half = min(n // 2, len(repeated))
    sample = repeated[:half]
    seen = {r["pid"] for r in sample}
    rest = [r for r in rows if r["pid"] not in seen]
    rng.shuffle(rest)
    sample += rest[: n - len(sample)]
    return sample


def write_fasta(data_root: Path, sample: list[dict], path: Path) -> None:
    wanted = {r["pid"] for r in sample}
    genomes = {r["genome"] for r in sample}
    found: dict[str, str] = {}
    for g in genomes:
        faa = data_root / "annotations_ab_initio" / g / "proteins.faa"
        if not faa.exists():
            continue
        name = None
        buf: list[str] = []
        with open(faa) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in wanted:
                        found[name] = "".join(buf)
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if name in wanted:
            found[name] = "".join(buf)
    with open(path, "w") as fh:
        for r in sample:
            seq = found.get(r["pid"])
            if not seq:
                continue
            fh.write(f">{r['pid']}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def run_cluster(db: Path, out: Path, identity: float, threads: int) -> None:
    subprocess.run(
        ["diamond", "cluster", "--db", str(db), "--out", str(out),
         "--approx-id", f"{identity:g}", "--threads", str(threads),
         "--memory-limit", "4G", "--quiet"],
        check=True,
    )


def parse_clusters(path: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                clusters[p[0]].append(p[1])
    return clusters


def metrics(clusters: dict[str, list[str]], key_of: dict[str, str]) -> dict:
    n = sum(len(m) for m in clusters.values())
    sizes = [len(m) for m in clusters.values()]
    mixed = 0
    same_in_cluster = 0
    all_in_cluster = 0
    weighted_pure = 0.0
    for members in clusters.values():
        if len(members) < 2:
            continue
        keys = Counter(key_of.get(x, x) for x in members)
        if len(keys) > 1:
            mixed += 1
        c = len(members)
        all_in_cluster += c * (c - 1) // 2
        same_in_cluster += sum(v * (v - 1) // 2 for v in keys.values())
        weighted_pure += max(keys.values()) / c * c
    total_same = 0
    by_key: Counter = Counter(key_of.get(x, x) for members in clusters.values() for x in members)
    for v in by_key.values():
        total_same += v * (v - 1) // 2
    precision = same_in_cluster / all_in_cluster if all_in_cluster else 1.0
    recall = same_in_cluster / total_same if total_same else 1.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "sequences": n,
        "clusters": len(sizes),
        "singletons": sum(1 for s in sizes if s == 1),
        "largest_cluster": max(sizes) if sizes else 0,
        "reduction": round(1 - len(sizes) / n, 4) if n else 0,
        "mixed_clusters": mixed,
        "pair_precision": round(precision, 4),
        "pair_recall": round(recall, 4),
        "f1": round(f1, 4),
        "weighted_purity": round(weighted_pure / n, 4) if n else 0,
        "same_key_pairs": total_same,
    }


def choose(rows: list[dict], min_precision: float = 0.95) -> dict:
    ok = [r for r in rows if r["pair_precision"] >= min_precision]
    if ok:
        return min(ok, key=lambda r: r["clusters"])
    return max(rows, key=lambda r: r["f1"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample", type=int, default=20000)
    ap.add_argument("--thresholds", default="80,90,95,97,100")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    rows = collect_annotations(a.data_root)
    print(f"anotaciones con Swiss-Prot: {len(rows):,}")
    sample = sample_rows(rows, a.sample, a.seed)
    key_of = {r["pid"]: r["key"] for r in sample}
    acc_of = {r["pid"]: r["acc"] for r in sample}
    fasta = a.out / "sample.faa"
    write_fasta(a.data_root, sample, fasta)
    print(f"muestra: {len(sample):,} proteínas -> {fasta}")

    db = a.out / "sample.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(fasta), "-d", str(db.with_suffix(""))],
                   check=True, capture_output=True)

    results = []
    for t in [float(x) for x in a.thresholds.split(",")]:
        out = a.out / f"clusters_{int(t)}.tsv"
        run_cluster(db, out, t, a.threads)
        clusters = parse_clusters(out)
        m = metrics(clusters, key_of)
        m["threshold"] = t
        results.append(m)
        print(f"  {t:>3.0f}%: clusters={m['clusters']:,} mezclados={m['mixed_clusters']:,} "
              f"precisión={m['pair_precision']:.3f} recall={m['pair_recall']:.3f} F1={m['f1']:.3f}")

    best = choose(results)
    with open(a.out / "dedup_thresholds.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    (a.out / "dedup_choice.json").write_text(json.dumps(
        {"chosen": best, "criterion": "menor nº de clusters con precisión >= 0.95",
         "results": results}, indent=1))
    print(f"\numbral elegido: {best['threshold']:.0f}%  "
          f"(precisión={best['pair_precision']}, recall={best['pair_recall']}, F1={best['f1']})")
    print(f"  -> {a.out / 'dedup_thresholds.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
