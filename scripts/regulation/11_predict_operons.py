#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""11_predict_operons.py — unidades de transcripción (L2) por genoma.

Método (combinación de evidencia, nunca una sola señal):
  * distancia intergénica, con umbral calibrado POR GENOMA a partir del valle
    de la distribución bimodal de huecos (no un 50 pb fijo);
  * acoplamiento traduccional (solapamiento / hueco <= 3 pb);
  * RBS del gen de abajo (ausente o débil => compatible con cotranscripción);
  * sintenía: mismo par de ortólogos adyacente en otros genomas del set
    (`data/regulation/synteny_pairs.tsv`, construido con --build-synteny);
  * coherencia funcional (texto de anotación: prefijo génico o palabra de ruta).

Un gen pertenece a un único operón; las coordenadas son monótonas; un ARN o un
cambio de hebra corta la unidad. No hay predictor de terminadores: no se usa su
ausencia como evidencia.

Salidas: `data/regulation/<G>/operons.tsv`, `operon_genes.tsv`, `meta.json`.

Uso:
  python scripts/regulation/11_predict_operons.py --build-synteny
  python scripts/regulation/11_predict_operons.py --genome GCF_...
  python scripts/regulation/11_predict_operons.py --all
  python scripts/regulation/11_predict_operons.py --selftest
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

import numpy as np

LO, HI, BIN = -20, 300, 5
T_HI_CAP = 400          # techo de la zona ambigua (huecos > esto nunca se unen)
MERGE_SCORE = 0.50      # banda lejana: exige 2 evidencias
REG_SCORE = 0.75        # un regulador suele tener promotor propio
REG_RE = re.compile(r"transcriptional (?:regulator|repressor|activator)|"
                    r"response regulator|two[- ]component|sigma factor", re.I)
LOOSE_SCORE = 0.35      # banda del bulto secundario: basta 1

OP_COLS = ["operon_id", "label", "contig", "start", "end", "strand", "n_genes",
           "genes", "gene_names", "gene_ids", "n_cds", "span_bp", "min_gap", "max_gap",
           "n_links_evidence", "status", "evidence", "attempted", "confidence"]
OG_COLS = ["operon_id", "position", "gene_id", "gene_name", "contig", "start",
           "end", "strand", "gap_to_prev", "link_rule", "link_score",
           "link_evidence"]


# --------------------------------------------------------------------------- #
# calibración del umbral
# --------------------------------------------------------------------------- #
def _hist(gaps):
    xs = np.clip(np.asarray(gaps, dtype=float), LO, HI + 0.01)
    idx = np.clip(((xs - LO) / BIN).astype(int), 0, int((HI - LO) / BIN))
    return np.bincount(idx, minlength=int((HI - LO) / BIN)).astype(float)


def _smooth(h, w=3):
    x = np.arange(-w, w + 1) / max(1.0, w / 2.0)
    k = np.exp(-0.5 * x * x)
    return np.convolve(h, k / k.sum(), mode="same")


def calibrate(gaps: list[int]) -> dict:
    """Valle de la distribución de huecos -> umbral t_lo (pb)."""
    if not gaps:
        return {"t_lo": 40, "t_hi": T_HI_CAP, "method": "fallback_empty",
                "hist": []}
    h = _hist(gaps)
    s = _smooth(h)
    c = LO + BIN * np.arange(len(h)) + BIN / 2.0
    lo_i = int(np.searchsorted(c, 20.0))
    hi_i = int(np.searchsorted(c, 140.0))
    seg = s[lo_i:hi_i]
    method = "histogram_valley"
    if len(seg) < 3:
        t_lo, method = 40, "fallback_range"
    else:
        # PRIMER mínimo local tras el pico pequeño (el valle entre los dos
        # modos), no el mínimo global de la cola larga.
        peak = float(s[:lo_i + 1].max()) if lo_i else 1.0
        cand = [j for j in range(1, len(seg) - 1)
                if seg[j] <= seg[j - 1] and seg[j] <= seg[j + 1]]
        cand = [j for j in cand if c[lo_i + j] >= 25.0
                and seg[j] <= 0.6 * max(peak, 1e-9)]
        if cand:
            t_lo = int(round(c[lo_i + cand[0]]))
        else:
            t_lo, method = 40, "fallback_no_valley"
    t_lo = int(min(80, max(15, t_lo)))
    # segundo nivel: final del "bulto" secundario de operones laxos
    t_mid = int(round(1.8 * t_lo))
    cand2 = [j for j in range(1, len(seg) - 1)
             if seg[j] <= seg[j - 1] and seg[j] <= seg[j + 1]
             and c[lo_i + j] >= t_lo + 15]
    if cand2:
        t_mid = int(round(c[lo_i + cand2[0]]))
    t_mid = int(min(140, max(t_lo + 10, t_mid)))
    # techo duro de la zona ambigua
    t_hi = int(min(T_HI_CAP, max(300, 4 * t_lo)))
    keep = [(int(LO + BIN * i), int(h[i])) for i in range(len(h)) if h[i] > 0]
    return {"t_lo": t_lo, "t_mid": t_mid, "t_hi": t_hi, "method": method,
            "n_gaps": len(gaps), "hist": keep}


# --------------------------------------------------------------------------- #
# sintenía global
# --------------------------------------------------------------------------- #
def build_synteny(limit: int = 0) -> dict:
    """Cuenta pares de ortólogos adyacentes (mismo contig/hebra) en el set."""
    dirs = sorted(p for p in C.AB.glob("GCF_*") if p.is_dir())
    if limit:
        dirs = dirs[:limit]
    counts: Counter = Counter()
    for i, d in enumerate(dirs, 1):
        g = C.load_genome(d.name)
        if g is None:
            continue
        for unit in C.transcription_units(g):
            for a, b in zip(unit, unit[1:]):
                if C.gap_between(a, b) > 400:
                    continue
                if a.acc and b.acc:
                    counts[(a.acc, b.acc)] += 1
        if i % 50 == 0:
            print(f"  [{i}/{len(dirs)}] pares={len(counts)}", flush=True)
    rows = [{"acc_a": a, "acc_b": b, "n_genomes": n}
            for (a, b), n in sorted(counts.items()) if n >= 2]
    C.write_tsv(C.REG / "synteny_pairs.tsv", rows,
                ["acc_a", "acc_b", "n_genomes"])
    meta = {"genomes": len(dirs), "pairs_all": len(counts),
            "pairs_kept": len(rows), "min_genomes": 2}
    C.dump_json(C.REG / "synteny_meta.json", meta)
    print(f"sintenía: {len(rows)} pares (>=2 genomas) de {len(counts)} en "
          f"{len(dirs)} genomas")
    return meta


def load_synteny(path: Path) -> dict:
    if not path.exists():
        return {}
    return {(r["acc_a"], r["acc_b"]): int(r["n_genomes"])
            for r in C.read_tsv(path)}


# --------------------------------------------------------------------------- #
# predicción
# --------------------------------------------------------------------------- #
def link_evidence(up: C.Gene, down: C.Gene, gap: int, syn: int,
                  stats: dict | None = None) -> dict:
    strong, why, spec = C.coherent(up, down, stats)
    rs = C.rbs_strength(down.rbs_motif)
    return {
        "gap": gap,
        "coupled": gap <= 3,
        "overlap": gap < 0,
        "rbs": rs,
        "syn": syn,
        "coherent": strong,
        "coh_why": why,
        "coh_spec": spec,
        "regulator": bool(REG_RE.search(down.product or "") or
                              REG_RE.search(up.product or "")),
    }


def link_score(ev: dict) -> float:
    s = 0.0
    if ev["syn"] >= 3:
        s += 0.45
    elif ev["syn"] == 2:
        s += 0.30
    elif ev["syn"] == 1:
        s += 0.15
    if ev["coherent"]:
        s += 0.30 + 0.25 * ev.get("coh_spec", 0.0)
    if ev["coupled"]:
        s += 0.25
    if ev["rbs"] == 0:
        s += 0.20
    elif ev["rbs"] == 1:
        s += 0.10
    return round(min(1.0, s), 4)


def decide_link(ev: dict, t_lo: int, t_mid: int, t_hi: int) -> tuple[bool, str, float]:
    """Dos niveles: la banda del bulto secundario exige menos evidencia."""
    gap = ev["gap"]
    score = link_score(ev)
    if gap <= t_lo:
        return True, "dist_close", score
    if gap > t_hi:
        return False, "dist_far", score
    if gap <= t_mid:
        return (score >= LOOSE_SCORE), "ambiguous_loose", score
    if ev.get("regulator"):
        # no separar por defecto: un regulador adyacente suele ser su propia TU
        return (score >= REG_SCORE), "ambiguous_regulator", score
    return (score >= MERGE_SCORE), "ambiguous_score", score


def _evidence_str(ev: dict, rule: str, merged: bool) -> str:
    parts = [f"gap={ev['gap']}", f"rule={rule}"]
    if ev["overlap"]:
        parts.append("overlap")
    if ev["coupled"]:
        parts.append("coupled")
    parts.append(f"rbs={ev['rbs']}")
    if ev["syn"]:
        parts.append(f"synteny={ev['syn']}")
    if ev["coherent"]:
        parts.append(ev["coh_why"])
    return ";".join(parts)


def operon_label(genes: list[C.Gene]) -> str:
    names = [g.gene_name or "" for g in genes]
    if len(genes) >= 2 and all(names):
        k = 0
        for chars in zip(*names):
            if len(set(chars)) != 1:
                break
            k += 1
        if k >= 3:
            return names[0][:k] + "".join(n[k:] for n in names)
    if names and names[0]:
        return names[0] if len(genes) == 1 else f"{names[0]}-{names[-1]}"
    return f"op_{genes[0].start}"


def predict(genome_id: str, synteny: dict, out_dir: Path) -> dict:
    t0 = time.time()
    g = C.load_input(genome_id)
    if g is None:
        raise SystemExit(f"genoma no encontrado: {genome_id}")
    units = C.transcription_units(g)
    all_gaps = [C.gap_between(a, b) for u in units for a, b in zip(u, u[1:])]
    stats = C.cooccurrence_stats(g.genes)
    ncrna: dict[tuple, list] = {}
    for x in g.genes:
        if x.gene_type == "ncRNA":
            ncrna.setdefault((x.contig, x.strand), []).append(x)
    cal = calibrate(all_gaps)
    t_lo, t_mid, t_hi = cal["t_lo"], cal["t_mid"], cal["t_hi"]

    ops, og = [], []
    for unit in units:
        groups: list[list[C.Gene]] = [[unit[0]]] if unit else []
        links: list[dict] = []
        for a, b in zip(unit, unit[1:]):
            gap = C.gap_between(a, b)
            syn = synteny.get((a.acc, b.acc), 0) if (a.acc and b.acc) else 0
            ev = link_evidence(a, b, gap, syn, stats)
            merged, rule, score = decide_link(ev, t_lo, t_mid, t_hi)
            links.append({"ev": ev, "rule": rule, "score": score,
                          "merged": merged})
            if merged:
                groups[-1].append(b)
            else:
                groups.append([b])
        for gi, genes in enumerate(groups):
            # sRNA/ncRNA dentro del intervalo del operón = miembro del transcrito
            lo, hi = min(x.start for x in genes), max(x.end for x in genes)
            inner = [r for r in ncrna.get((genes[0].contig, genes[0].strand), [])
                     if r.start <= hi and r.end >= lo]   # solape, no contención
            if inner:
                genes = sorted(genes + inner,
                               key=lambda x: (-x.end if x.strand == -1 else x.start))
            start = min(x.start for x in genes)
            end = max(x.end for x in genes)
            strand = genes[0].strand
            oid = f"{genome_id}:{genes[0].contig}:{start}:{'+' if strand == 1 else '-'}"
            ln = links[max(0, gi - 1):gi] if gi > 0 else []
            idx0 = 0
            # vínculos internos de este grupo
            base = sum(len(x) for x in groups[:gi])
            n_cds = sum(1 for x in genes if not x.is_rna)
            internal = links[base:base + n_cds - 1] if n_cds > 1 else []
            cds_i = 0
            for j, gene in enumerate(genes):
                lk = None
                if not gene.is_rna and cds_i > 0 and cds_i - 1 < len(internal):
                    lk = internal[cds_i - 1]
                if not gene.is_rna:
                    cds_i += 1
                og.append({
                    "operon_id": oid, "position": j + 1,
                    "gene_id": gene.gene_id, "gene_name": gene.gene_name,
                    "contig": gene.contig, "start": gene.start,
                    "end": gene.end, "strand": gene.strand,
                    "gap_to_prev": lk["ev"]["gap"] if lk else "",
                    "link_rule": lk["rule"] if lk else "",
                    "link_score": C.fnum(lk["score"]) if lk else "",
                    "link_evidence": _evidence_str(lk["ev"], lk["rule"],
                                                   True) if lk else "",
                })
            n_ev = sum(1 for lk in internal
                       if lk["ev"]["syn"] or lk["ev"]["coherent"])
            amb = any(lk["rule"] == "ambiguous_score" for lk in internal)
            if n_cds == 1 and len(genes) == 1:
                st, conf = "consensus", "low"
                why = "unidad monocistrónica; no hubo evidencia para unir vecinos"
            elif amb:
                st, conf = "anomalous", "low"
                why = "al menos un vínculo se unió en la zona ambigua"
            elif n_ev == len(internal) and internal:
                st, conf = "validated", "high"
                why = "todos los vínculos con evidencia independiente de distancia"
            else:
                st, conf = "consensus", "medium"
                why = "vínculos unidos solo por distancia"
            d = C.decide(st, validated_by="11_predict_operons",
                         reason=why, attempted="distancia;RBS;sintenía;función",
                         confidence=conf,
                         evidence=f"t_lo={t_lo};t_mid={t_mid};t_hi={t_hi};n={len(genes)}")
            ops.append(C.status_row({
                "operon_id": oid, "label": operon_label([x for x in genes if not x.is_rna] or genes),
                "contig": genes[0].contig, "start": start, "end": end,
                "strand": strand, "n_genes": len(genes), "n_cds": n_cds,
                # b-number/locus_tag primero: es el ID que resuelve RegulonDB
                "genes": ",".join(x.locus_tag or x.gene_name or x.gene_id
                                  for x in genes),
                "gene_names": ",".join(x.gene_name or x.gene_id for x in genes),
                "gene_ids": ",".join(x.gene_id for x in genes),
                "span_bp": end - start + 1,
                "min_gap": min((lk["ev"]["gap"] for lk in internal),
                               default=""),
                "max_gap": max((lk["ev"]["gap"] for lk in internal),
                               default=""),
                "n_links_evidence": n_ev,
            }, d))

    ops.sort(key=lambda r: (r["contig"], int(r["start"])))
    oid_order = {r["operon_id"]: i for i, r in enumerate(ops)}
    og.sort(key=lambda r: (oid_order[r["operon_id"]], r["position"]))
    C.write_tsv(out_dir / "operons.tsv", ops, OP_COLS)
    C.write_tsv(out_dir / "operon_genes.tsv", og, OG_COLS)

    multi = [o for o in ops if int(o["n_cds"]) >= 2]
    n_cds = len(g.cds())
    meta = {
        "genome": genome_id, "n_units": len(ops), "n_multigene": len(multi),
        "n_singleton": len(ops) - len(multi),
        "n_cds": n_cds,
        "pct_cds_in_multigene": C.fnum(
            100.0 * sum(int(o["n_cds"]) for o in multi) / max(1, n_cds), 2),
        "mean_genes_multigene": C.fnum(
            sum(int(o["n_cds"]) for o in multi) / max(1, len(multi)), 3),
        "max_genes": max((o["n_genes"] for o in ops), default=0),
        "threshold": cal,
        "merge_score": MERGE_SCORE,
        "status_counts": dict(sorted(Counter(o["status"] for o in ops).items())),
        "elapsed_s": round(time.time() - t0, 2),
        "sha256_operons": C.sha256_files([out_dir / "operons.tsv",
                                          out_dir / "operon_genes.tsv"]),
    }
    C.dump_json(out_dir / "operons_meta.json", meta)
    return meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--build-synteny", action="store_true")
    ap.add_argument("--synteny", type=Path, default=C.REG / "synteny_pairs.tsv")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        return selftest()
    if a.build_synteny:
        build_synteny(a.limit)
        return 0
    syn = load_synteny(a.synteny)
    if not syn:
        print("aviso: sin tabla de sintenía; se sigue sin esa evidencia",
              file=sys.stderr)

    if a.all:
        dirs = sorted(p for p in C.AB.glob("GCF_*") if p.is_dir())
        if a.limit:
            dirs = dirs[:a.limit]
        t0 = time.time()
        ok = 0
        for i, d in enumerate(dirs, 1):
            try:
                m = predict(d.name, syn, C.REG / d.name)
                ok += 1
                if i % 50 == 0:
                    print(f"  [{i}/{len(dirs)}] operones={m['n_units']} "
                          f"t_lo={m['threshold']['t_lo']}", flush=True)
            except Exception as e:                        # noqa: BLE001
                print(f"  FALLO {d.name}: {e}", file=sys.stderr)
        print(f"operones: {ok}/{len(dirs)} genomas en {time.time()-t0:.1f}s")
        return 0 if ok else 1

    if not a.genome:
        ap.error("indica --genome, --all o --build-synteny")
    m = predict(a.genome, syn, C.REG / a.genome)
    print(f"{a.genome}: {m['n_units']} unidades, {m['n_multigene']} "
          f"multigénicas, {m['pct_cds_in_multigene']}% CDS, "
          f"media {m['mean_genes_multigene']}, t_lo={m['threshold']['t_lo']} "
          f"({m['threshold']['method']})")
    return 0


def selftest() -> int:
    import tempfile
    # calibración: distribución bimodal sintética
    gaps = [-5, 0, 2, 5, 8, 10, 12, 15, 18, 20, 25, 30, 35] * 4 + \
           [120, 150, 180, 200, 250, 300] * 10
    cal = calibrate(gaps)
    assert 15 <= cal["t_lo"] <= 80, cal
    assert cal["t_hi"] >= cal["t_lo"], cal
    # unión por distancia
    ev = {"gap": 5, "coupled": False, "overlap": False, "rbs": 2, "syn": 0,
          "coherent": False, "coh_why": ""}
    assert decide_link(ev, 40, 70, 200)[0] is True
    ev_far = dict(ev, gap=500)
    assert decide_link(ev_far, 40, 70, 200)[0] is False
    # zona ambigua: sin evidencia no se une, con sintenía+coherencia sí
    ev_amb = {"gap": 80, "coupled": False, "overlap": False, "rbs": 2,
              "syn": 0, "coherent": False, "coh_why": ""}
    assert decide_link(ev_amb, 40, 70, 200)[0] is False
    ev_amb2 = dict(ev_amb, syn=4, coherent=True, coh_why="prefix:lac", coh_spec=1.0)
    assert decide_link(ev_amb2, 40, 70, 200)[0] is True
    ev_amb3 = dict(ev_amb, coherent=True, coh_why="prefix:lac", coh_spec=1.0)
    assert decide_link(ev_amb3, 40, 90, 200)[0] is True
    assert decide_link(ev_amb3, 40, 70, 200)[0] is True   # prefijo raro basta lejos
    ev_amb4 = dict(ev_amb, coherent=True, coh_why="kw:protein", coh_spec=0.0)
    assert decide_link(ev_amb4, 40, 70, 200)[0] is False  # token común no
    # etiqueta de operón
    gs = [C.Gene("c", 1, 3, 1, "1", gene_name="lacZ"),
          C.Gene("c", 5, 7, 1, "2", gene_name="lacY"),
          C.Gene("c", 9, 11, 1, "3", gene_name="lacA")]
    assert operon_label(gs) == "lacZYA", operon_label(gs)
    # predicción sobre un genoma sintético
    g = C.Genome("syn", "gbff")
    g.contigs = {"c": 2000}
    g.genes = [
        C.Gene("c", 100, 400, 1, "a", gene_name="lacZ", acc="A"),
        C.Gene("c", 405, 700, 1, "b", gene_name="lacY", acc="B"),
        C.Gene("c", 1500, 1800, 1, "c2", gene_name="foo", acc="C"),
        C.Gene("c", 250, 500, -1, "d", gene_name="bar", acc="D"),
    ]
    with tempfile.TemporaryDirectory(dir=str(C.REG)) as td:
        C.REG.joinpath(".tmp").mkdir(parents=True, exist_ok=True)
        old = C.load_input
        C.load_input = lambda *x, **k: g
        try:
            m = predict("syn", {("A", "B"): 3}, Path(td))
        finally:
            C.load_input = old
        ops = C.read_tsv(Path(td) / "operons.tsv")
        assert m["n_units"] == 3, (m["n_units"], ops)
        assert any(o["n_genes"] == "2" for o in ops), ops
        # determinismo: mismo hash en dos corridas
        h1 = m["sha256_operons"]
        C.load_input = lambda *x, **k: g
        try:
            m2 = predict("syn", {("A", "B"): 3}, Path(td))
        finally:
            C.load_input = old
        assert m2["sha256_operons"] == h1
    print("selftest 11_predict_operons: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
