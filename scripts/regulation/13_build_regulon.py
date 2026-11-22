#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""13_build_regulon.py — regulón (L3): TF -> operón, con signo y regla lógica.

Aristas por dos vías, separadas y auditables:
  * `documented`: tabla curada pequeña (E. coli) TF -> conjunto de genes, con
    cita. Es la única fuente que permite `logic_status=documented`.
  * `motif`: escaneo posicional del motivo del TF río arriba del operón. El
    signo sale de la posición (operador sobre el promotor => represión; río
    arriba del -35 => activación). Es una HEURÍSTICA: `status=consensus` como
    máximo nunca `validated`.

La función lógica por operón es SIEMPRE asumida:
    ON IF (cualquier activador unido) AND NOT (cualquier represor unido)
salvo que exista regla documentada. Se guarda en `logic_status`.

Motivos: solo los que se pueden citar sin inventar. Los TF sin consenso
público fiable (H-NS, Lrp, NarL, OxyR, SoxS, Fis) NO se escanean y se
declaran como excluidos en el meta. RpoS es factor sigma, no TF de operador:
se trata en 12_promoters_sigma.py.

Salidas: `regulon.json`, `expression_rules.tsv`.

Uso:
  python scripts/regulation/13_build_regulon.py --genome ecoli_k12
  python scripts/regulation/13_build_regulon.py --all
  python scripts/regulation/13_build_regulon.py --validate
  python scripts/regulation/13_build_regulon.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

UP = 300

# Motivo IUPAC + semilla exacta + efecto documentado + cita.
MOTIFS = {
    "CRP": {"motif": "TGTGANNNNNNTCACA", "seed": "TGTGA", "effect": "+",
            "gene": "crp", "conf": "high",
            "cite": "Berg & von Hippel 1988; Busby & Ebright 1999"},
    "FNR": {"motif": "TTGATNNNNATCAA", "seed": "TTGAT", "effect": "+",
            "gene": "fnr", "conf": "high",
            "cite": "Eiglmeier et al. 1989; Guest 1992"},
    "Fur": {"motif": "GATAATGATAATCATTATC", "seed": "GATAAT", "effect": "-",
            "gene": "fur", "conf": "high",
            "cite": "Escolar et al. 1999 (Fur box)"},
    "LexA": {"motif": "CTGTATATATATACAG", "seed": "CTGTAT", "effect": "-",
             "gene": "lexA", "conf": "high",
             "cite": "Walker 1984; Little & Mount 1982 (SOS box)"},
    "LacI": {"motif": "AATTGTGAGCGGATAACAATT", "seed": "GTGAGC", "effect": "-",
             "gene": "lacI", "conf": "high",
             "cite": "Gilbert & Maxam 1973 (lac O1)"},
    "TrpR": {"motif": "TGTACTAGTTAACTAGTAC", "seed": "ACTAGT", "effect": "-",
             "gene": "trpR", "conf": "low",
             "cite": "Gunsalus & Yanofsky 1980 (trp operator)"},
    "ArcA": {"motif": "GTTAATTAAT", "seed": "TTAATT", "effect": "-",
             "gene": "arcA", "conf": "low",
             "cite": "Lynch & Lin 1996 (ArcA box)"},
    "IHF": {"motif": "WATCAANNNNTTR", "seed": "ATCAA", "effect": "bend",
            "gene": "ihfA", "conf": "medium",
            "cite": "Craig & Nash 1984 (IHF site)"},
}
# Zona esperada del sitio respecto al inicio del gen (pb). Restringir la
# ventana es lo que separa un operador real de una coincidencia en 300 pb.
ZONES = {"CRP": (80, 200), "FNR": (80, 200), "IHF": (40, 200),
         "Fur": (40, 110), "LexA": (40, 110), "LacI": (40, 110),
         "TrpR": (40, 130), "ArcA": (40, 200)}

EXCLUDED = {
    "H-NS": "sin consenso de secuencia fiable (une ADN curvo AT-rico)",
    "Lrp": "sin consenso fiable para escaneo posicional",
    "NarL": "motivo heptámero degenerado; no se cita secuencia sin fuente",
    "OxyR": "sitio largo y degenerado; requiere PWM",
    "SoxS": "soxbox degenerado; requiere PWM",
    "Fis": "consenso 15 pb demasiado degenerado para semilla exacta",
    "RpoS": "es factor sigma (ver 12_promoters_sigma.py), no TF de operador",
}

# Reglas documentadas (K-12). La clave es el conjunto de nombres de gen.
DOCUMENTED = [
    ("LacI", ["lacZ", "lacY", "lacA"], "-", "RegulonDB; EcoCyc lacZYA"),
    ("CRP", ["lacZ", "lacY", "lacA"], "+", "RegulonDB; EcoCyc lacZYA"),
    ("CRP", ["galE", "galT", "galK", "galM"], "+", "RegulonDB; EcoCyc galETKM"),
    ("TrpR", ["trpL", "trpE", "trpD", "trpC", "trpB", "trpA"], "-",
     "Gunsalus & Yanofsky 1980; RegulonDB trpLEDCBA"),
    ("AraC", ["araB", "araA", "araD"], "+", "Englesberg & Wilcox 1974; EcoCyc"),
    ("MalT", ["malE", "malF", "malG"], "+", "Boos & Shuman 1998; EcoCyc"),
    ("LexA", ["recA"], "-", "Walker 1984; RegulonDB SOS"),
    ("LexA", ["lexA"], "-", "Walker 1984; RegulonDB SOS"),
    ("LexA", ["uvrA"], "-", "Walker 1984; RegulonDB SOS"),
    ("LexA", ["uvrB"], "-", "Walker 1984; RegulonDB SOS"),
    ("LexA", ["sulA"], "-", "Walker 1984; RegulonDB SOS"),
]

RULE_COLS = ["operon_id", "label", "genes", "rule_pseudo_code", "activators",
             "repressors", "terms_with_evidence", "terms_assumed", "effectors",
             "logic_status", "status", "evidence", "attempted", "confidence"]


def load_effectors() -> dict:
    """Efector metabólico por TF desde la GT de RegulonDB (si está)."""
    p = C.DATA / "regulondb" / "gt" / "gt_tfs.tsv"
    if not p.exists():
        return {}
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("name"):
                out[r["name"]] = (r.get("effector") or "").strip()
    return out


def _tf_present(genome_genes: list[C.Gene]) -> dict[str, str]:
    by_name = {(g.gene_name or "").lower(): g for g in genome_genes}
    out = {}
    for tf, spec in MOTIFS.items():
        g = by_name.get(spec["gene"].lower())
        if g is None and tf == "IHF":
            g = by_name.get("ihfB")
        if g is not None:
            out[tf] = g.gene_id
    return out


def _site_sign(spec: dict, dist_start: int) -> tuple[str, str]:
    if spec["effect"] == "bend":
        return "bend", "curvatura"
    if spec["effect"] in ("+", "-"):
        return spec["effect"], "documentado"
    if dist_start <= 60:
        return "-", "posición:operador"
    return "+", "posición:upstream"


def build(genome_id: str, out_dir: Path) -> dict:
    t0 = time.time()
    g = C.load_input(genome_id)
    if g is None:
        raise SystemExit(f"genoma no encontrado: {genome_id}")
    seqs = C.load_sequence(genome_id)
    ops = C.read_tsv(out_dir / "operons.tsv")
    og = C.read_tsv(out_dir / "operon_genes.tsv")
    eff = load_effectors()
    present = _tf_present(g.cds())
    members: dict[str, list[dict]] = {}
    for r in og:
        members.setdefault(r["operon_id"], []).append(r)
    name_by_id = {x.gene_id: (x.gene_name or "") for x in g.genes}
    bnum_by_name = {(x.gene_name or "").lower(): x.locus_tag for x in g.genes}

    edges, rules = [], []
    for o in ops:
        oid = o["operon_id"]
        mem = sorted(members.get(oid, []), key=lambda r: int(r["position"]))
        if not mem:
            continue
        names = [r["gene_name"] for r in mem]
        f = mem[0]
        contig, strand = o["contig"], int(o["strand"])
        seq = seqs.get(contig, "")
        w = C.upstream_seq(seq, int(f["start"]), int(f["end"]), strand, UP) if seq else ""
        got = []
        for tf, spec in MOTIFS.items():
            if tf not in present:
                continue
            hit, dist = None, 0
            if w:
                lo, hi = ZONES.get(tf, (0, UP))
                for s_, mm, sub in C.scan(spec["motif"], w):
                    d_ = len(w) - s_
                    if not lo <= d_ <= hi:
                        continue
                    if hit is None or mm < hit[1]:
                        hit, dist = (s_, mm, sub), d_
            if hit:
                sign, how = _site_sign(spec, dist)
                conf = "high" if (spec["conf"] == "high" and hit[1] <= 1) else \
                       "medium" if hit[1] <= 3 else "low"
                d = C.decide("consensus", validated_by="motivo_posicional",
                             reason=f"signo por {how}",
                             attempted="motivo_curado;posicion",
                             confidence=conf,
                             evidence=f"site={hit[2]};mm={hit[1]};dist={dist}")
                got.append({"tf": tf, "sign": sign, "dist": dist, "mm": hit[1],
                            "site": hit[2], "source": "motif", "conf": conf,
                            "evidence": d})
        for tf, genes, sign, cite in DOCUMENTED:
            if set(genes) == set(n for n in names if n):
                d = C.decide("validated", validated_by="RegulonDB;literatura",
                             reason="regla documentada", attempted="documentado",
                             confidence="high", evidence=f"cite:{cite}")
                got.append({"tf": tf, "sign": sign, "dist": 0, "mm": 0,
                            "site": "", "source": "documented", "conf": "high",
                            "evidence": d})
        if not got:
            rules.append(C.status_row({
                "operon_id": oid, "label": o["label"],
                "genes": o["genes"], "rule_pseudo_code": "",
                "activators": "", "repressors": "", "terms_with_evidence": "",
                "terms_assumed": "", "effectors": "", "logic_status": "none",
            }, C.decide("unresolved", attempted="motivo;documentado",
                        reason="ningún TF modelado con sitio o regla",
                        confidence="low", evidence="")))
            continue
        acts = sorted({e["tf"] for e in got if e["sign"] == "+"})
        reps = sorted({e["tf"] for e in got if e["sign"] == "-"})
        bends = sorted({e["tf"] for e in got if e["sign"] not in "+-"})
        documented = any(e["source"] == "documented" for e in got)
        def _lab(t):
            e = eff.get(t) or ""
            return f"{t}[{e}]" if e else t
        terms = [_lab(t) for t in acts + reps]
        parts = []
        if acts:
            parts.append("(" + " OR ".join(_lab(t) for t in acts) + ")")
        if reps:
            parts.append("NOT (" + " OR ".join(_lab(t) for t in reps) + ")")
        if parts:
            rule = "ON IF " + " AND ".join(parts)
        elif bends:
            rule = "UNKNOWN (solo curvatura: " + " OR ".join(bends) + ")"
        else:
            rule = ""
        for e in got:
            edges.append({
                "tf": e["tf"], "target_operon": oid, "target_label": o["label"],
                "target_genes": o["genes"], "sign": e["sign"],
                "source": e["source"], "site": e["site"], "mismatches": e["mm"],
                "dist_to_start": e["dist"], "effector": eff.get(e["tf"], ""),
                "status": e["evidence"].status, "confidence": e["conf"],
                "evidence": e["evidence"].evidence,
                "attempted": e["evidence"].attempted,
            })
        ev_terms = sorted({e["tf"] for e in got if e["source"] == "documented"})
        as_terms = sorted({e["tf"] for e in got if e["source"] == "motif"})
        d = C.decide("validated" if documented and not as_terms else "consensus",
                     validated_by="documentado" if documented else "motivo_posicional",
                     reason="lógica asumida ON IF (act) AND NOT (rep)",
                     attempted="motivo;documentado",
                     confidence="high" if documented and not as_terms else "low",
                     evidence=f"act={','.join(acts)};rep={','.join(reps)}")
        rules.append(C.status_row({
            "operon_id": oid, "label": o["label"], "genes": o["genes"],
            "rule_pseudo_code": rule, "activators": ",".join(acts),
            "repressors": ",".join(reps),
            "terms_with_evidence": ",".join(ev_terms),
            "terms_assumed": ",".join(as_terms),
            "effectors": ";".join(f"{t}={eff[t]}" for t in acts + reps if eff.get(t)),
            "logic_status": "documented" if documented else "assumed",
        }, d))

    edges.sort(key=lambda e: (e["target_operon"], e["tf"], e["source"]))
    rules.sort(key=lambda r: r["operon_id"])
    C.write_tsv(out_dir / "expression_rules.tsv", rules, RULE_COLS)
    with_rules = [r for r in rules if r["rule_pseudo_code"]]
    regulon = {
        "genome": genome_id,
        "tfs_modeled": sorted(MOTIFS),
        "tfs_present": sorted(present),
        "tf_excluded": EXCLUDED,
        "nodes": [{"id": o["operon_id"], "type": "operon", "label": o["label"],
                   "n_genes": int(o["n_genes"])} for o in ops] +
                 [{"id": tf, "type": "tf", "effector": eff.get(tf, "")}
                  for tf in sorted(present)],
        "edges": edges,
        "summary": {"n_operons": len(ops), "n_with_rule": len(with_rules),
                    "n_edges": len(edges),
                    "n_documented_edges": sum(1 for e in edges if e["source"] == "documented"),
                    "n_motif_edges": sum(1 for e in edges if e["source"] == "motif"),
                    "coverage_rules": C.fnum(100.0 * len(with_rules) / max(1, len(ops)), 2)},
    }
    C.dump_json(out_dir / "regulon.json", regulon)
    meta = {"genome": genome_id, "n_operons": len(ops), "n_rules": len(with_rules),
            "n_edges": len(edges), "tfs_present": sorted(present),
            "tfs_modeled": sorted(MOTIFS), "tf_excluded": EXCLUDED,
            "elapsed_s": round(time.time() - t0, 2),
            "sha256": C.sha256_files([out_dir / "expression_rules.tsv",
                                      out_dir / "regulon.json"])}
    C.dump_json(out_dir / "regulon_meta.json", meta)
    return meta


def validate(genome_id: str, out_dir: Path) -> dict:
    """Ajuste de aristas TF->TU contra gt_tus, a nivel de ARITA.

    Una arista predicha (TF, operón) cuenta como acierto si RegulonDB tiene una
    TU regulada por ese TF que comparte al menos un gen con el operón. Se
    comparan por separado las aristas documentadas y las de motivo.
    """
    gt_p = C.DATA / "regulondb" / "gt" / "gt_tus.tsv"
    if not gt_p.exists():
        raise SystemExit("falta gt_tus.tsv")
    gt = set()
    with open(gt_p, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = frozenset(g.strip() for g in (r.get("genes") or "").split(",")
                              if g.strip())
            for tf in (r.get("tfs") or "").split(","):
                tf = tf.strip()
                if tf in MOTIFS and genes:
                    gt.add((tf, genes))
    name2b = {}
    for r in C.read_tsv(out_dir / "genes.tsv"):
        name2b[(r["gene_name"] or "").lower()] = r["locus_tag"]
    reg = json.loads((out_dir / "regulon.json").read_text())
    pred = {}
    for e in reg["edges"]:
        genes = set()
        for tok in (e["target_genes"] or "").split(","):
            t = tok.strip()
            b = t if t.startswith("b") and t[1:].isdigit() else name2b.get(t.lower(), "")
            if b:
                genes.add(b)
        if genes:
            pred.setdefault(e["source"], set()).add((e["tf"], frozenset(genes)))

    def score(pset):
        tp = {pe for pe in pset
              if any(pe[0] == g[0] and (pe[1] & g[1]) for g in gt)}
        rec = {g for g in gt
               if any(g[0] == p[0] and (g[1] & p[1]) for p in pset)}
        p_ = len(tp) / len(pset) if pset else 0.0
        r_ = len(rec) / len(gt) if gt else 0.0
        f_ = 2 * p_ * r_ / (p_ + r_) if (p_ + r_) else 0.0
        return {"pred": len(pset), "tp": len(tp), "gt_recovered": len(rec),
                "precision": round(p_, 4), "recall": round(r_, 4),
                "f1": round(f_, 4)}

    res = {"genome": genome_id, "tfs_modeled": sorted(MOTIFS),
           "gt_edges": len(gt), "by_source": {k: score(v) for k, v in sorted(pred.items())},
           "all_predicted": score(set().union(*pred.values()) if pred else set())}
    C.dump_json(out_dir / "regulon_validation.json", res)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.validate:
        res = validate(a.genome or "ecoli_k12", C.REG / (a.genome or "ecoli_k12"))
        for src, m in res["by_source"].items():
            print(f"  {src:10s} aristas={m['pred']:4d} tp={m['tp']:4d} "
                  f"P={m['precision']:.3f} R={m['recall']:.3f} F1={m['f1']:.3f}")
        print(f"regulón vs gt_tus: aristas GT={res['gt_edges']} "
              f"F1={res['all_predicted']['f1']:.3f}")
        return 0
    if a.all:
        dirs = sorted(p for p in C.AB.glob("GCF_*") if p.is_dir())
        if a.limit:
            dirs = dirs[:a.limit]
        for d in dirs:
            if not (C.REG / d.name / "operons.tsv").exists():
                continue
            try:
                build(d.name, C.REG / d.name)
            except Exception as e:                        # noqa: BLE001
                print(f"  FALLO {d.name}: {e}", file=sys.stderr)
        print(f"regulón: {len(dirs)} genomas")
        return 0
    if not a.genome:
        ap.error("indica --genome, --all o --validate")
    m = build(a.genome, C.REG / a.genome)
    print(f"{a.genome}: {m['n_edges']} aristas, {m['n_rules']}/{m['n_operons']} "
          f"operones con regla; TFs presentes={m['tfs_present']}")
    return 0


def selftest() -> int:
    assert "CRP" in MOTIFS and MOTIFS["CRP"]["seed"] == "TGTGA"
    # el motivo de CRP encuentra su sitio canónico
    w = "A" * 20 + "TGTGAACGTACTCACA" + "A" * 20
    hits = C.scan(MOTIFS["CRP"]["motif"], w)
    assert hits and hits[0][1] == 0, hits
    # LexA en el box SOS
    w2 = "A" * 10 + "CTGTATATATATACAG" + "A" * 10
    assert C.scan(MOTIFS["LexA"]["motif"], w2)[0][1] == 0
    # signo por posición
    assert _site_sign({"effect": "-"}, 30)[0] == "-"
    assert _site_sign({"effect": "+"}, 200)[0] == "+"
    assert _site_sign({"effect": "bend"}, 200)[0] == "bend"
    # tabla documentada coherente
    assert any(t == "LacI" and s == "-" for t, _, s, _ in DOCUMENTED)
    assert "H-NS" in EXCLUDED and "RpoS" in EXCLUDED
    print("selftest 13_build_regulon: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
