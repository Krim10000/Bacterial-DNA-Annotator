#!/usr/bin/env python3
"""pipeline.py — TSS y TTS sobre regiones candidatas, para cualquier especie.

Idea
----
Ni un promotor ni un terminador viven dentro de un gen, así que analizar el
genoma entero desperdicia ~80% del trabajo. Pero recortar el espacio intergénico
en seco pierde los terminadores que se solapan con el final del gen. La solución
es partir del intergénico y **entrar 50 pb en los genes vecinos**.

Qué hace
--------
1. Lee un GBFF y construye regiones candidatas (intergénico + margen).
2. **TTS**: escanea cada región buscando horquilla + poli-U (Rho-independiente).
3. **TSS**: si se le pasa un TSV de TSS (pos, strand), los asigna al gen que
   empiezan; si no, se omite (no se inventan).
4. **Clasifica anomalías** de cada detección: intergénica, antisentido, dentro
   de un gen en el mismo sentido, o en el borde de un contig.
5. Si se da un conjunto de referencia, mide cobertura a ±N pb, estratificada.

Por qué importa la clasificación de anomalías
---------------------------------------------
De los 375 terminadores curados de K-12, 11 son **antisentido** (terminan
transcritos de la hebra contraria, no operones) y 19 caen **dentro de un gen**,
probablemente por desajuste de coordenadas. Mezclarlos con los demás infla la
cobertura y contamina cualquier benchmark. Aquí salen etiquetados, no escondidos.

Salidas
-------
    <out>.tsv      una fila por detección, con su clase y su ubicación
    <out>.json     resumen por tipo y por clase de anomalía

Uso
---
    python scripts/regulation/pipeline.py \
        --gbff data/test_ecoli.gbff \
        --tss data/regulondb/tss_k12.tsv \
        --ref-tts data/regulondb/TerminatorSet.tsv \
        --out data/operon_finder/k12_tss_tts
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import annotate_hairpins as H  # noqa: E402
import terminator_classes as T  # noqa: E402

FEATURES = ("CDS", "tRNA", "rRNA", "ncRNA")
EDGE = 100          # pb desde el borde del contig para marcar "borde"
NEAR_GENE = 600     # pb máximo para asignar un TSS a un gen


def load_genome(gbff: Path):
    from Bio import SeqIO
    contigs = {}
    for rec in SeqIO.parse(str(gbff), "genbank"):
        feats = []
        for f in rec.features:
            if f.type not in FEATURES:
                continue
            strand = "+" if f.location.strand == 1 else "-"
            name = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag")
                    or f.qualifiers.get("product") or [f.type])[0]
            for p in f.location.parts:
                feats.append({"start": int(p.start), "end": int(p.end),
                              "strand": strand, "type": f.type, "name": name})
        feats.sort(key=lambda x: x["start"])
        contigs[rec.id] = {"seq": str(rec.seq), "feats": feats}
    return contigs


def build_regions(feats: list[dict], L: int, margin: int, min_len: int):
    merged = []
    for f in feats:
        a, b = f["start"], f["end"]
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    inter, prev = [], 0
    for a, b in merged:
        if a > prev:
            inter.append((prev, a))
        prev = max(prev, b)
    if prev < L:
        inter.append((prev, L))
    out = []
    for a, b in inter:
        if b - a < min_len:
            continue
        out.append((max(0, a - margin), min(L, b + margin)))
    out.sort()
    merged2 = []
    for a, b in out:
        if merged2 and a <= merged2[-1][1]:
            merged2[-1] = (merged2[-1][0], max(merged2[-1][1], b))
        else:
            merged2.append((a, b))
    return merged2


def location(feats: list[dict], pos: int, strand: str, L: int,
             margin: int = 50) -> tuple[str, dict | None]:
    """Dónde cae una detección, separando borde de gen de interior.

    El margen entra a propósito 50 pb en los genes vecinos, así que una detección
    allí NO es una anomalía. Solo se marca `inside_gene` si está más adentro que
    el margen: eso sí es un caso raro (antisentido, atenuador interno, o
    desajuste de coordenadas).
    """
    if pos < EDGE or pos > L - EDGE:
        return "contig_edge", None
    hit = None
    for f in feats:
        if f["start"] <= pos <= f["end"]:
            hit = f
            if f["strand"] != strand:
                return "antisense", f
    if hit:
        borde = min(pos - hit["start"], hit["end"] - pos)
        if borde <= margin + 10:
            return "gene_edge", hit
        return "inside_gene", hit
    return "intergenic", None


def gene_downstream(feats: list[dict], pos: int, strand: str, maxd: int = NEAR_GENE):
    best = None
    for f in feats:
        if f["type"] != "CDS" or f["strand"] != strand:
            continue
        if strand == "+":
            d = f["start"] - pos
            if 0 <= d <= maxd and (best is None or d < best[0]):
                best = (d, f)
        else:
            d = pos - f["end"]
            if 0 <= d <= maxd and (best is None or d < best[0]):
                best = (d, f)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--tss", default=None, help="TSV con pos,strand (opcional)")
    ap.add_argument("--ref-tts", default=None, help="TSS/TTS de referencia para medir cobertura")
    ap.add_argument("--batter", default=None,
                    help="BED de terminadores de BATTER-TPE; se añaden como TTS "
                         "de clase 'batter' (detector externo, mayor recall)")
    ap.add_argument("--margin", type=int, default=50)
    ap.add_argument("--min-region", type=int, default=20)
    ap.add_argument("--max-region", type=int, default=400,
                    help="las regiones enormes se truncan para el escaneo")
    ap.add_argument("--tol", type=int, default=25, help="tolerancia de comparación (pb)")
    ap.add_argument("--out", required=True, help="prefijo de salida")
    a = ap.parse_args()

    contigs = load_genome(Path(a.gbff))
    print(f"contigs: {len(contigs)} | "
          f"pb totales: {sum(len(c['seq']) for c in contigs.values()):,}")

    rows = []
    n_regions = 0
    for contig, info in contigs.items():
        L = len(info["seq"])
        seq, feats = info["seq"], info["feats"]
        regions = build_regions(feats, L, a.margin, a.min_region)
        n_regions += len(regions)
        for start, end in regions:
            chunk = seq[start:min(end, start + a.max_region)].upper()
            hp = H.best_hairpin(chunk, use_dg=False)
            if not hp:
                continue
            pos = start + int(hp.get("start", 0))
            left = [f for f in feats if f["end"] <= pos and f["type"] == "CDS"]
            strand = max(left, key=lambda f: f["end"])["strand"] if left else "+"
            loc, feat = location(feats, pos, strand, L, a.margin)
            rows.append({"type": "TTS", "contig": contig, "pos": pos,
                         "strand": strand, "class": hp.get("cls", ""),
                         "stem": hp.get("stem", ""), "polyU": hp.get("polyU", ""),
                         "gc": round(hp.get("gc", 0), 3),
                         "location": loc,
                         "feature": feat["name"] if feat else "",
                         "gene": "", "dist_gene": ""})

        if a.batter:
            pass  # se procesa tras recorrer los contigs (necesita el mapa completo)

        if a.tss:
            for r in csv.DictReader(open(a.tss), delimiter="\t"):
                p = (r.get("pos") or "").strip()
                if not p.isdigit():
                    continue
                pos = int(p)
                st = (r.get("strand") or "+").strip()
                st = "-" if st in ("-", "-1", "reverse") else "+"
                if not (0 <= pos <= L):
                    continue
                loc, feat = location(feats, pos, st, L, a.margin)
                g = gene_downstream(feats, pos, st)
                rows.append({"type": "TSS", "contig": contig, "pos": pos,
                             "strand": st, "class": (r.get("sigma") or "").strip(),
                             "stem": "", "polyU": "", "gc": "",
                             "location": loc,
                             "feature": feat["name"] if feat else "",
                             "gene": g[1]["name"] if g else "",
                             "dist_gene": g[0] if g else ""})

    if a.batter and Path(a.batter).exists():
        n_b = 0
        for line in open(a.batter):
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].lstrip("-").isdigit():
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            pos = (start + end) // 2
            contig = f[0]
            info = contigs.get(contig) or next(iter(contigs.values()))
            feats, L = info["feats"], len(info["seq"])
            left = [x for x in feats if x["end"] <= pos and x["type"] == "CDS"]
            strand = max(left, key=lambda x: x["end"])["strand"] if left else "+"
            loc, feat = location(feats, pos, strand, L, a.margin)
            rows.append({"type": "TTS", "contig": contig, "pos": pos,
                         "strand": strand, "class": "batter", "stem": "",
                         "polyU": "", "gc": "", "location": loc,
                         "feature": feat["name"] if feat else "",
                         "gene": "", "dist_gene": ""})
            n_b += 1
        print(f"terminadores BATTER añadidos: {n_b}")

    rows.sort(key=lambda x: (x["contig"], x["pos"]))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tsv = out.with_suffix(".tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["type", "contig", "pos", "strand",
                                           "class", "stem", "polyU", "gc",
                                           "location", "feature", "gene",
                                           "dist_gene"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # resumen
    resumen = {"regiones": n_regions,
               "total": len(rows),
               "por_tipo": dict(Counter(r["type"] for r in rows)),
               "por_clase": dict(Counter(f"{r['type']}:{r['location']}" for r in rows)),
               "por_clase_tts": dict(Counter(r["class"] for r in rows if r["type"] == "TTS"))}

    if a.ref_tts:
        ref = []
        for r in csv.DictReader(open(a.ref_tts), delimiter="\t"):
            v = (r.get("2)leftEndPos") or r.get("start") or r.get("pos") or "").strip()
            if v.isdigit():
                ref.append((int(v), (r.get("7)type") or r.get("type") or "").strip()))
        tts = sorted(r["pos"] for r in rows if r["type"] == "TTS")
        import bisect
        ok = 0
        for p, _ in ref:
            i = bisect.bisect_left(tts, p - a.tol)
            if i < len(tts) and tts[i] <= p + a.tol:
                ok += 1
        resumen["cobertura_ref"] = {"n": len(ref), "encontrados": ok,
                                    "pct": round(100 * ok / max(1, len(ref)), 1),
                                    "tol": a.tol}

    json_path = out.with_suffix(".json")
    json_path.write_text(json.dumps(resumen, indent=1, ensure_ascii=False))

    print(f"regiones escaneadas: {n_regions}")
    print(f"detecciones: {len(rows)}  {resumen['por_tipo']}")
    for k, v in sorted(resumen["por_clase"].items()):
        print(f"  {k:26s} {v}")
    if "cobertura_ref" in resumen:
        c = resumen["cobertura_ref"]
        print(f"referencia: {c['encontrados']}/{c['n']} ({c['pct']}%) a ±{c['tol']} pb")
    print(f"-> {tsv}\n-> {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
