#!/usr/bin/env python3
"""build_reannotation.py -- propuestas de re-anotación/validación para los 77 accionables.

Une, por caso (`genome`,`gap`):
  1. la clasificación del revisor (`data/metabolism/review/gap_review.tsv`),
  2. la arquitectura de dominios esperada (`NO_HOMOLOG_DOMAINS.tsv`),
  3. los hits de dominio del proteoma (`domain_hits.tsv`, ver `scan_domains.py`),
  4. la anotación del ORF candidato (`data/annotations_ab_initio/<g>/functional.tsv`
     + coordenadas de `annotation.gff3`),
  5. la reacción/reacción productora (`out_fl/refined/refined_gaps.tsv` y
     `rhea2ec.tsv`), verificando con `rhea_reactions.pkl` que el metabolito esté
     del lado producto de la reacción citada,
  6. la ausencia de homólogo dirigido (`homology_closure_top1000.tsv`).

Regla de decisión (documentada en README.md):
  - `reannotate`     : un dominio **específico del productor** aparece a --cut_ga
                       en un ORF cuya anotación actual es hipotética/débil o ya
                       nombra la enzima productora (reconectar mapeo).
  - `validate_domain`: hay señal de dominio productor (específico a umbral
                       relajado, o dominio de familia amplia a --cut_ga, o
                       arquitectura incompleta) pero no alcanza para re-anotar.
  - `keep_missing`   : sin dominio del productor; lo hallado es dominio de la
                       enzima consumidora, genérico, o de una familia no
                       diagnóstica. La brecha queda como está.

Salidas (solo en `data/metabolism/reannotation/`):
  - `reannotation.tsv` (11 columnas)
  - `summary.json`
No se modifica la tabla base ni ningún otro directorio.
"""
from __future__ import annotations

import csv
import hashlib
import json
import pickle
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / "data" / "metabolism" / "review" / "gap_review.tsv"
NOHOM = ROOT / "data" / "metabolism" / "NO_HOMOLOG_DOMAINS.tsv"
HITS = ROOT / "data" / "metabolism" / "reannotation" / "domain_hits.tsv"
REFINED = ROOT / "out_fl" / "refined" / "refined_gaps.tsv"
RHEA2EC = ROOT / "data" / "metabolism" / "rhea2ec.tsv"
RHEA_PKL = ROOT / "data" / "metabolism" / "rhea_reactions.pkl"
HOMCLOS = ROOT / "data" / "metabolism" / "homology_closure_top1000.tsv"
MAPPING_EC = ROOT / "data" / "metabolism" / "mapping_ec.tsv"
ANNOT = ROOT / "data" / "annotations_ab_initio"
OUTDIR = ROOT / "data" / "metabolism" / "reannotation"
OUT_TSV = OUTDIR / "reannotation.tsv"
OUT_SUM = OUTDIR / "summary.json"

ACTIONABLE = ("annotation_gap", "domain_present")
COLS = ["genome", "gap", "family", "support", "classification", "domain",
        "candidate", "decision", "reason", "evidence"]
DECISIONS = ("reannotate", "validate_domain", "keep_missing")

# --- Configuración por gap -------------------------------------------------
# `key_domains` son los dominios Pfam que identifican a la enzima PRODUCTORA
# (verificados contra las accesiones de referencia de los productores y contra
# Pfam). `broad=True` marca familias amplias que solas no justifican re-anotar.
# `consumer` son dominios de pasos posteriores/consumidores, que el revisor
# contó como "arquitectura" pero que no prueban que el productor esté presente.
GAP_CFG = {
    "6-hydroxymethyl-7,8-dihydropterin": {
        "enzyme": "dihidroneopterin aldolasa (FolB)",
        "ec": "4.1.2.25",
        "key": {"FolB": "specific"},
        "consumer": ["HPPK", "Pterin_bind"],
        "producer_rhea": ["RHEA:10540", "RHEA:51016"],
        "terms": ["dihydroneopterin", "folb", "aldolase"],
        "note": "HPPK (FolK) y Pterin_bind (FolP) son pasos posteriores, no el productor",
    },
    "S-adenosyl 3-(methylsulfanyl)propylamine": {
        "enzyme": "S-adenosilmetionina descarboxilasa (SpeD)",
        "ec": "4.1.1.50",
        "key": {"AdoMet_dc": "specific"},
        "consumer": ["Spermine_synth", "Spermine_synt_N"],
        "producer_rhea": ["RHEA:15981"],
        "terms": ["adenosylmethionine decarboxylase", "sped", "decarboxylase"],
        "note": "Spermine_synth es la espermidina sintasa (SpeE), que CONSUME el metabolito",
    },
    "(4R,5S)-dethiobiotin": {
        "enzyme": "dethiobiotina sintasa (BioD)",
        "ec": "6.3.3.3",
        "key": {"AAA_26": "broad"},
        "consumer": ["Aminotran_3"],
        "producer_rhea": ["RHEA:15805", "RHEA:63684"],
        "terms": ["dethiobiotin", "biod"],
        "note": "AAA_26 (PF13500) es familia AAA amplia; Aminotran_3 es genérico",
    },
    "(R)-4'-phosphopantetheine": {
        "enzyme": "fosfopantotenoilcisteína descarboxilasa (CoaC/PPCDC)",
        "ec": "4.1.1.36",
        "key": {"Flavoprotein": "shared", "DFP": "specific", "ACP_PD": "specific"},
        "consumer": ["NUDIX"],
        "producer_rhea": ["RHEA:16793", "RHEA:20537", "RHEA:64988", "RHEA:67592"],
        "terms": ["phosphopantothenoyl", "coac", "phosphopantetheine"],
        "note": "NUDIX es superfamilia amplia; Flavoprotein/DFP (PF02441/PF04127) son la PPCDC",
    },
    "L-cysteinylglycine": {
        "enzyme": "gama-glutamiltransferasa (EC 3.4.19.13) / ChaC (EC 4.3.2.7) / EC 3.5.1.128",
        "ec": "3.4.19.13/4.3.2.7/3.5.1.128",
        "key": {"G_glu_transpept": "specific", "ChaC": "specific",
                "CN_hydrolase": "broad"},
        "consumer": [],
        "producer_rhea": ["RHEA:28807", "RHEA:47724", "RHEA:54532"],
        "terms": ["cysteinylglycine", "gamma-glutamyl", "glutathione", "chac"],
        "note": "CN_hydrolase es superfamilia carbon-nitrogen hidrolasa",
    },
    "5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate": {
        "enzyme": "AIR carboxilasa / PurE (EC 4.1.1.21 / 5.4.99.18)",
        "ec": "4.1.1.21/5.4.99.18",
        "key": {"AIRC": "specific"},
        "consumer": ["ATP-grasp", "PurK_C", "RS_preATP-grasp-like", "SAICAR_synt"],
        "producer_rhea": ["RHEA:10792", "RHEA:13193"],
        "terms": ["phosphoribosylaminoimidazole carboxylase", "pure", "air carboxylase"],
        "note": "PurE necesita el paso previo PurK (ATP-grasp/PurK_C); arquitectura incompleta",
    },
    "L-ornithine": {
        "enzyme": "productor no identificable (55 reacciones; sin Pfam consenso)",
        "ec": "n/a",
        "key": {},
        "consumer": ["Ala_racemase_N", "Ala_racemase_C", "ArgJ", "Asp_Glu_race"],
        "producer_rhea": ["RHEA:11940", "RHEA:13201", "RHEA:13265", "RHEA:15349"],
        "terms": ["ornithine"],
        "note": "Ala_racemase/Asp_Glu_race son racemasas/aminotransferasas no diagnósticas; "
                "RHEA:11584 (ornithine racemase) consume L-ornitina en dirección master",
    },
    "(2S)-2-[5-amino-1-(5-phospho-beta-D-ribosyl)imidazole-4-carboxamido]succinate": {
        "enzyme": "SAICAR sintasa (PurC, EC 6.3.2.6)",
        "ec": "6.3.2.6",
        "key": {"SAICAR_synt": "specific"},
        "consumer": ["AIRC", "ATP-grasp", "PurK_C", "RS_preATP-grasp-like"],
        "producer_rhea": ["RHEA:22628"],
        "terms": ["saicar", "phosphoribosylaminoimidazole-succinocarboxamide", "purc"],
        "note": "AIRC/ATP-grasp/PurK_C son pasos previos del mismo brazo de purinas",
    },
    "UDP-N-acetyl-alpha-D-muramate": {
        "enzyme": "UDP-N-acetilenolpiruvoilglucosamina reductasa (MurB, EC 1.3.1.98)",
        "ec": "1.3.1.98",
        "key": {"MurB_C": "specific"},
        "consumer": ["FAD_binding_4"],
        "producer_rhea": ["RHEA:12248", "RHEA:29563", "RHEA:53716"],
        "terms": ["udp-n-acetylenolpyruvoylglucosamine reductase", "murb"],
        "note": "FAD_binding_4 es dominio de unión a FAD genérico",
    },
}

WEAK_WORDS = ("uncharacterized", "uncharacterised", "hypothetical",
              "putative", "probable", "duf", "upf", "unknown protein",
              "protein of unknown")


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm_rhea(x: str) -> int:
    return int(x.split(":")[-1])


def load_rhea() -> dict:
    d = pickle.load(open(RHEA_PKL, "rb"))
    return {int(k): v for k, v in d.items()}


def load_annotation(genome: str) -> dict[str, dict]:
    """ORF -> {acc, name, pident, desc} desde functional.tsv (DIAMOND vs Swiss-Prot)."""
    path = ANNOT / genome / "functional.tsv"
    out = {}
    if not path.exists():
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        acc = f[1]
        name = acc.split("|")[1] if acc.count("|") >= 2 else acc
        out[f[0]] = {"acc": name, "pident": f[2], "desc": f[-1]}
    return out


def load_coords(genome: str) -> dict[str, str]:
    """ORF -> 'contig:start-end(strand)' desde annotation.gff3."""
    path = ANNOT / genome / "annotation.gff3"
    out = {}
    if not path.exists():
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        m = re.search(r"ID=([^;]+)", f[8])
        if m:
            out[m.group(1)] = f"{f[0]}:{f[3]}-{f[4]}({f[6]})"
    return out


def classify_annotation(ann: dict | None, terms: list[str]) -> tuple[str, str]:
    """Devuelve (estado, detalle): unannotated/weak/match/conflict."""
    if not ann or not ann.get("desc"):
        return "unannotated", "sin hit Swiss-Prot"
    desc = ann["desc"]
    low = desc.lower()
    detail = f"{ann['acc']} pident={ann['pident']}"
    if any(w in low for w in WEAK_WORDS):
        return "weak", detail
    if any(t in low for t in terms):
        return "match", detail
    return "conflict", detail


def main() -> int:
    cases = [r for r in read_tsv(REVIEW) if r["classification"] in ACTIONABLE]
    expected = {(r["genome"], r["gap"]): r
                for r in read_tsv(NOHOM)}
    refined = {(r["genome"], r["gap"]): r for r in read_tsv(REFINED)}
    homclos = {(r["genome"], r["gap"]): r for r in read_tsv(HOMCLOS)}
    mapped = {(r["genome"], r["gap"]): r for r in read_tsv(MAPPING_EC)}
    rhea2ec = {}
    for r in read_tsv(RHEA2EC):
        rhea2ec[int(r["RHEA_ID"])] = r["ID"]
    rhea = load_rhea()
    hits_by_genome = defaultdict(list)
    for r in read_tsv(HITS):
        hits_by_genome[r["genome"]].append(r)

    ann_cache: dict[str, dict] = {}
    coord_cache: dict[str, dict] = {}
    out_rows = []
    per_case_debug = []

    for c in cases:
        key = (c["genome"], c["gap"])
        genome, gap = key
        cfg = GAP_CFG[gap]
        exp = expected.get(key, {})
        dom_expected = exp.get("domains", "")
        ref = refined.get(key, {})
        transport = c.get("transport", "")
        support = c.get("support", "")
        # Verificación química de las reacciones productoras citadas.
        chebi = ref.get("chebi", "")
        ver_prod, ver_sub = [], []
        for rid in (ref.get("producers", "") or "").split(";"):
            if not rid:
                continue
            rr = rhea.get(norm_rhea(rid), {})
            if chebi and chebi in rr.get("products", set()):
                ver_prod.append(rid)
            elif chebi and chebi in rr.get("substrates", set()):
                ver_sub.append(rid)

        ec_str = ",".join(sorted({rhea2ec.get(norm_rhea(r), "?")
                                  for r in ver_prod})) or "?"

        if genome not in ann_cache:
            ann_cache[genome] = load_annotation(genome)
            coord_cache[genome] = load_coords(genome)
        ann = ann_cache[genome]

        keyd = cfg["key"]
        ga = [h for h in hits_by_genome.get(genome, [])
              if h["mode"] == "ga" and h["domain"] in keyd]
        rel = [h for h in hits_by_genome.get(genome, [])
               if h["mode"] == "relaxed" and h["domain"] in keyd]
        cons_hits = [h for h in hits_by_genome.get(genome, [])
                     if h["mode"] == "ga" and h["domain"] in cfg["consumer"]]

        # Mejor hit por (ORF, dominio).
        def best(rows):
            d = {}
            for h in rows:
                k = (h["orf"], h["domain"])
                if k not in d or float(h["dom_ievalue"]) < float(d[k]["dom_ievalue"]):
                    d[k] = h
            return d
        best_ga, best_rel = best(ga), best(rel)
        by_orf = defaultdict(dict)
        for (orf, dom), h in best_ga.items():
            by_orf[orf][dom] = h
        # ORF candidato: más dominios clave, luego mayor score, luego más corto.
        cand_orf = ""
        if by_orf:
            ranked = sorted(
                by_orf.items(),
                key=lambda kv: (-len(kv[1]),
                                -max(float(x["dom_score"]) for x in kv[1].values()),
                                int(next(iter(kv[1].values()))["orf_len"])))
            cand_orf = ranked[0][0]
        cand_hits = by_orf.get(cand_orf, {})
        ann_state, ann_detail = classify_annotation(ann.get(cand_orf), cfg["terms"])
        cand_desc = ann.get(cand_orf, {}).get("desc", "")
        cand_coord = coord_cache[genome].get(cand_orf, "")
        cand_len = next(iter(cand_hits.values()))["orf_len"] if cand_hits else ""

        specific = [d for d, s in keyd.items() if s == "specific"]
        shared = [d for d, s in keyd.items() if s in ("shared", "broad")]
        ga_specific = sorted(d for d in cand_hits if d in specific)
        ga_shared = sorted(d for d in cand_hits if d in shared)

        def fmt_hit(h):
            return (f"{h['domain']}@{h['orf']} ie={h['dom_ievalue']} "
                    f"cov={h['ali_cov']} len={h['orf_len']}")

        # ---------------- decisión ----------------
        rel_best = sorted(best_rel.items(),
                          key=lambda kv: float(kv[1]["dom_ievalue"]))[0] \
            if best_rel else None

        def qualifies_relaxed(orf: str, h: dict) -> bool:
            need = 0.5 if keyd.get(h["domain"]) == "specific" else 0.8
            state = classify_annotation(ann.get(orf), cfg["terms"])[0]
            return float(h["ali_cov"]) >= need and state in ("unannotated", "weak")

        if ga_specific:
            doms = "+".join(ga_specific)
            if ann_state in ("unannotated", "weak"):
                decision = "reannotate"
                reason = (f"dominio productor {doms} a --cut_ga en ORF {cand_orf} "
                          f"({cand_len} aa, {ann_state}); proponer "
                          f"{cfg['enzyme']} (EC {cfg['ec']}) y conectar la reacción "
                          f"productora")
            elif ann_state == "match":
                decision = "reannotate"
                reason = (f"enzima ya nombrada ({cand_desc[:70]}) y dominio productor "
                          f"{doms} a --cut_ga en el mismo ORF {cand_orf}; re-anotar = "
                          f"reconectar el mapeo a la reacción productora (EC {cfg['ec']})")
            elif ann_state == "conflict" and _pident(ann_detail) < 55:
                decision = "reannotate"
                reason = (f"dominio productor {doms} a --cut_ga en ORF {cand_orf} con "
                          f"anotación distinta y débil ({cand_desc[:70]}, "
                          f"pident={_pident(ann_detail)}); proponer {cfg['enzyme']} "
                          f"(EC {cfg['ec']})")
            else:
                decision = "validate_domain"
                reason = (f"dominio productor {doms} a --cut_ga en ORF {cand_orf}, "
                          f"pero anotación fuerte distinta ({cand_desc[:70]}, "
                          f"pident={_pident(ann_detail)}); validar dominio/estructura "
                          f"antes de re-anotar")
            if gap.startswith("5-amino-1-(5-phospho-D-ribosyl)imidazole"):
                same_orf = [h for h in cons_hits if h["orf"] == cand_orf]
                purk = any(h["domain"] in ("ATP-grasp", "PurK_C") for h in cons_hits)
                if same_orf:
                    decision = "reannotate"
                    reason = (f"AIRC (PurE) y {same_orf[0]['domain']} en el mismo ORF "
                              f"{cand_orf} (proteína bifuncional); la anotación actual "
                              f"({cand_desc[:60]}) solo cubre el paso posterior: añadir "
                              f"la actividad PurE (EC {cfg['ec']}) para producir CAIR")
                elif not purk:
                    decision = "validate_domain"
                    reason = (f"AIRC (PurE) a --cut_ga en ORF {cand_orf}, pero falta el "
                              f"paso previo PurK (ATP-grasp/PurK_C) y no hay dominio "
                              f"posterior en el mismo ORF; validar la ruta antes de "
                              f"re-anotar")
        elif ga_shared:
            doms = "+".join(ga_shared)
            if ann_state in ("unannotated", "weak") or (
                    ann_state == "conflict" and _pident(ann_detail) < 55):
                decision = "validate_domain"
                reason = (f"dominio de familia amplia/compartida {doms} a --cut_ga en "
                          f"ORF {cand_orf} (anotación {ann_state}, "
                          f"pident={_pident(ann_detail) or 'n/d'}); no es diagnóstico "
                          f"por sí solo: validar dominio/estructura y EC antes de "
                          f"re-anotar")
            else:
                decision = "keep_missing"
                reason = (f"dominio de familia amplia {doms} a --cut_ga en ORF "
                          f"{cand_orf} pero la anotación es fuerte y distinta "
                          f"({cand_desc[:60]}, pident={_pident(ann_detail)}): no "
                          f"alcanza para re-anotar, la brecha sigue")
        elif rel_best and qualifies_relaxed(rel_best[0][0], rel_best[1]):
            (orf, dom), h = rel_best
            cand_orf = orf
            cand_len = h["orf_len"]
            cand_coord = coord_cache[genome].get(orf, "")
            ann_state, ann_detail = classify_annotation(ann.get(orf), cfg["terms"])
            cand_desc = ann.get(orf, {}).get("desc", "")
            decision = "validate_domain"
            reason = (f"señal del dominio productor {dom} solo a umbral relajado "
                      f"(ie={h['dom_ievalue']}, cov={h['ali_cov']}) en ORF {orf} "
                      f"({h['orf_len']} aa, {ann_state}); sin hit a --cut_ga: validar "
                      f"HMM/estructura antes de re-anotar")
        elif rel_best:
            (orf, dom), h = rel_best
            cand_orf = orf
            cand_len = h["orf_len"]
            cand_coord = coord_cache[genome].get(orf, "")
            ann_state, ann_detail = classify_annotation(ann.get(orf), cfg["terms"])
            cand_desc = ann.get(orf, {}).get("desc", "")
            decision = "keep_missing"
            reason = (f"hit parcial de {dom} solo a umbral relajado "
                      f"(ie={h['dom_ievalue']}, cov={h['ali_cov']}) en ORF {orf}, "
                      f"anotado como {cand_desc[:50] or 'sin hit'}: no alcanza para "
                      f"proponer re-anotación, la brecha sigue")
        else:
            decision = "keep_missing"
            if not keyd:
                reason = (f"sin modelo Pfam específico del productor para este gap "
                          f"({cfg['note']}); los dominios hallados no permiten elegir "
                          f"una enzima productora: la brecha sigue")
            elif cons_hits:
                cons = "+".join(sorted({h["domain"] for h in cons_hits}))
                reason = (f"solo dominios de paso posterior/genérico ({cons}); el "
                          f"dominio productor ({'/'.join(keyd) or 'sin modelo'}) no "
                          f"aparece: no hay base para re-anotar, la brecha sigue")
            else:
                reason = (f"sin hit del dominio productor ({'/'.join(keyd) or 'sin modelo'}) "
                          f"ni en --cut_ga ni a umbral relajado; la brecha sigue")

        # ---------------- evidencia ----------------
        ev = [f"revision={c['classification']}",
              f"dominio_esperado={dom_expected or 'sin_dominios'}",
              f"clave_productor={'/'.join(keyd) or 'sin_modelo'}"]
        if cand_hits:
            ev.append("hits_origen=" + ";".join(fmt_hit(h) for h in cand_hits.values()))
        if best_rel and not best_ga:
            ev.append("hits_relajados=" + ";".join(fmt_hit(h) for h in best_rel.values()))
        if cons_hits:
            ev.append("hits_paso_posterior=" + ";".join(
                sorted({f"{h['domain']}@{h['orf']}" for h in cons_hits})))
        if cand_orf:
            ev.append(f"orf={cand_orf} coord={cand_coord} len={cand_len}")
            ev.append(f"anotacion={cand_desc[:160] or 'sin_hit'} ({ann_detail})")
        ev.append("productor_verificado=" + (
            f"{';'.join(ver_prod)} EC=" + ",".join(
                sorted({rhea2ec.get(norm_rhea(r), '?') for r in ver_prod}))
            if ver_prod else "ninguna reacción cita el gap como producto en dirección master"))
        if ver_sub:
            ev.append("reacciones_que_lo_consumen=" + ";".join(ver_sub))
        hc = homclos.get(key, {}).get("status", "sin_fila_sin_busqueda_dirigida")
        ev.append(f"homologo={c.get('homolog_status', '')} homolog_closure={hc}")
        ev.append(f"soporte={support}/400 transporte={transport}")
        if ref:
            ev.append(f"label={ref.get('label', '')}")
        ev.append("mapeo_ec=" + ("fila_presente" if key in mapped else "sin_fila"))
        ev.append("nota=" + cfg["note"])
        evidence = " | ".join(x.replace("\t", " ").replace("\n", " ") for x in ev)

        row = {
            "genome": genome, "gap": gap, "family": c["family"],
            "support": support, "classification": c["classification"],
            "domain": dom_expected,
            "candidate": cand_orf if decision != "keep_missing" else "",
            "decision": decision,
            "reason": re.sub(r"\s+", " ", reason).strip(), "evidence": evidence,
        }
        out_rows.append(row)
        per_case_debug.append({
            "genome": genome, "gap": gap,
            "decision": decision, "cand": cand_orf,
            "ga_specific": ga_specific, "ga_shared": ga_shared,
            "n_ga": len(ga), "n_rel": len(rel), "n_cons": len(cons_hits),
            "ann_state": ann_state, "ann": cand_desc[:80],
        })

    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # ---------------- validación y summary ----------------
    keys = [(r["genome"], r["gap"]) for r in out_rows]
    dec = Counter(r["decision"] for r in out_rows)
    cross = Counter((r["classification"], r["decision"]) for r in out_rows)
    by_family = Counter(r["family"] for r in out_rows)
    fam_dec = defaultdict(Counter)
    gap_dec = defaultdict(Counter)
    for r in out_rows:
        fam_dec[r["family"]][r["decision"]] += 1
        gap_dec[r["gap"]][r["decision"]] += 1
    checks = {
        "rows": len(out_rows),
        "unique_keys": len(set(keys)),
        "duplicates": len(keys) - len(set(keys)),
        "expected_77": len(keys) == 77,
        "empty_fields": sum(1 for r in out_rows for k in COLS
                            if k not in ("candidate",) and not str(r[k]).strip()),
        "invalid_decisions": sum(1 for r in out_rows if r["decision"] not in DECISIONS),
        "decisions_sum_to_rows": sum(dec.values()) == len(out_rows),
        "decision_split_by_class": {
            f"{k[0]}|{k[1]}": v for k, v in sorted(cross.items())},
        "evidence_has_no_tabs_or_newlines": all(
            "\t" not in r["evidence"] and "\n" not in r["evidence"] for r in out_rows),
        "reason_has_no_tabs_or_newlines": all(
            "\t" not in r["reason"] and "\n" not in r["reason"] for r in out_rows),
        "candidate_exists_or_empty": all(
            (not r["candidate"]) or _candidate_in_proteome(r["genome"], r["candidate"])
            for r in out_rows),
        "class_counts": dict(Counter(r["classification"] for r in out_rows)),
    }
    summary = {
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(
            timespec="seconds"),
        "agent": "reannotation-builder",
        "status": "proposed_not_confirmed",
        "scope": "77 accionables de gap_review.tsv (17 annotation_gap + 60 domain_present)",
        "source_of_truth": {
            "base_table_not_modified": "data/metabolism/gap_interpretation.tsv",
            "review": "data/metabolism/review/gap_review.tsv",
            "domain_architecture": "data/metabolism/NO_HOMOLOG_DOMAINS.tsv",
            "domain_hits": "data/metabolism/reannotation/domain_hits.tsv (hmmscan --cut_ga y relajado)",
            "proteomes": "data/annotations_ab_initio/<genome>/proteins.faa (pyrodigal)",
            "orf_annotation": "data/annotations_ab_initio/<genome>/functional.tsv (DIAMOND vs Swiss-Prot)",
            "reactions": [
                "out_fl/refined/refined_gaps.tsv",
                "data/metabolism/rhea2ec.tsv",
                "data/metabolism/rhea_reactions.pkl",
            ],
            "no_homolog": "data/metabolism/homology_closure_top1000.tsv",
        },
        "decision_policy": {
            "reannotate": "dominio específico del productor a --cut_ga en ORF hipotético/débil o ya anotado como la enzima; propuesta de anotación/mapeo",
            "validate_domain": "señal de dominio productor insuficiente: familia amplia, umbral relajado o arquitectura incompleta",
            "keep_missing": "sin dominio del productor; lo hallado es consumidor/genérico; la brecha no se cierra",
        },
        "counts": {
            "total": len(out_rows),
            "genomes": len({r["genome"] for r in out_rows}),
            "by_classification": dict(Counter(r["classification"] for r in out_rows)),
            "by_decision": dict(dec),
            "class_x_decision": {f"{k[0]}|{k[1]}": v for k, v in sorted(cross.items())},
            "candidate_orfs": sum(1 for r in out_rows if r["candidate"]),
            "by_gap": dict(Counter(r["gap"] for r in out_rows)),
            "by_family": dict(by_family),
            "family_x_decision": {f: dict(c) for f, c in sorted(fam_dec.items())},
            "gap_x_decision": {g: dict(c) for g, c in sorted(gap_dec.items())},
        },
        "validation": checks,
        "caveats": [
            "Solo se usó evidencia de dominio HMM (Pfam) y anotación existente; no se hizo búsqueda de homología nueva ni validación EC/experimental.",
            "Los hits 'relaxed' (E<=1e-3) son señal divergente, no prueba de función.",
            "El panel cubre los 26 modelos Pfam que usan los 77 casos (NO_HOMOLOG_DOMAINS.domains + AdoMet_dc); otros dominios del productor no fueron probados.",
            "Para dcAdoMet el revisor contó Spermine_synth (SpeE, consumidora); el productor SpeD (AdoMet_dc) se escaneó aparte: ver decisión por caso.",
            "La anotación de ORF es el mejor hit Swiss-Prot (DIAMOND), no una anotación curada; 'putative' se trata como débil.",
            "4 reacciones listadas como productoras de dcAdoMet por refined_gaps lo consumen en dirección master (RHEA:23244/33387/36487/42864); solo RHEA:15981 lo produce.",
            "No se modificó data/metabolism/gap_interpretation.tsv, webapp/, data/annotations_ab_initio/ ni data/regulation/. Sin commits.",
        ],
        "how_to_reproduce": [
            "python3 scripts/evidence/scan_domains.py --force",
            "python3 scripts/evidence/build_reannotation.py",
            "python3 scripts/evidence/verify_reannotation.py",
        ],
    }
    with open(OUT_SUM, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")

    print(json.dumps({"rows": len(out_rows), "decisions": dict(dec),
                      "class_x_decision": summary["counts"]["class_x_decision"],
                      "validation": checks}, indent=2, ensure_ascii=False))
    return 0


def _pident(detail: str) -> float:
    m = re.search(r"pident=([0-9.]+)", detail or "")
    return float(m.group(1)) if m else 0.0


def _ec(cfg: dict, rhea2ec: dict) -> str:
    ecs = sorted({rhea2ec.get(norm_rhea(r), "?") for r in cfg["producer_rhea"]})
    return ",".join(ecs) if ecs else "?"


def _candidate_in_proteome(genome: str, orf: str) -> bool:
    path = ANNOT / genome / "proteins.faa"
    if not path.exists():
        return False
    for line in open(path):
        if line.startswith(">") and line[1:].strip() == orf:
            return True
    return False


if __name__ == "__main__":
    raise SystemExit(main())
