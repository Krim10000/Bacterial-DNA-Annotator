#!/usr/bin/env python3
"""Anota la capa de membrana/transporte de un genoma bacteriano.

Entradas (solo lectura):
  data/annotations_ab_initio/<genome>/annotation.gff3
  data/annotations_ab_initio/<genome>/proteins.faa
  data/annotations_ab_initio/<genome>/functional.tsv  (DIAMOND vs Swiss-Prot)

Salidas:
  data/membrane/<genome>.tsv
  data/membrane/gff/<genome>.gff3

Metodo:
  - hmmscan contra un panel Pfam curado (transporte, secrecion, membrana,
    lipoproteinas, senal).
  - Heuristica de hidropatia Kyte-Doolittle para segmentos TM (no es TMHMM).
  - Heuristica de signal peptide N-terminal (no es SignalP).
  - Motivo lipobox [LVI][ASTM][AG]C en los primeros 40 aa (no es LipoP).
  - Sustrato SOLO si el hit Swiss-Prot lo nombra explicitamente.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
import tempfile
import urllib.parse
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2, "X": 0.0, "B": 0.0, "Z": 0.0, "U": 0.0, "O": 0.0,
}

SUBSTRATES = [
    # orden: mas especifico primero
    (r"ferric\s+citrate", "ferric citrate"),
    (r"ferric\s+enterobactin", "ferric enterobactin"),
    (r"ferric\s+hydroxamate", "ferric hydroxamate"),
    (r"ferric\s+siderophore", "ferric siderophore"),
    (r"sulfate", "sulfate"),
    (r"thiosulfate", "thiosulfate"),
    (r"sulfonate", "sulfonate"),
    (r"phosphonate", "phosphonate"),
    (r"phosphate", "phosphate"),
    (r"molybdate", "molybdate"),
    (r"tungstate", "tungstate"),
    (r"nitrate", "nitrate"),
    (r"nitrite", "nitrite"),
    (r"ammonium|ammonia", "ammonium"),
    (r"urea", "urea"),
    (r"\bzinc\b", "zinc"),
    (r"manganese", "manganese"),
    (r"magnesium", "magnesium"),
    (r"\bnickel\b", "nickel"),
    (r"cobalt", "cobalt"),
    (r"copper", "copper"),
    (r"calcium", "calcium"),
    (r"potassium|\bK\+", "potassium"),
    (r"sodium|\bNa\+", "sodium"),
    (r"ferrous|\biron\b(?![- ]sulfur)|\bFe2\+|\bFe3\+", "iron"),
    (r"haem|heme", "heme"),
    (r"siderophore", "siderophore"),
    (r"lactose", "lactose"),
    (r"maltose", "maltose"),
    (r"maltodextrin", "maltodextrin"),
    (r"sucrose", "sucrose"),
    (r"trehalose", "trehalose"),
    (r"cellobiose", "cellobiose"),
    (r"chitobiose", "chitobiose"),
    (r"glucose", "glucose"),
    (r"galactose", "galactose"),
    (r"fructose", "fructose"),
    (r"mannose", "mannose"),
    (r"rhamnose", "rhamnose"),
    (r"fucose", "fucose"),
    (r"xylose", "xylose"),
    (r"arabinose", "arabinose"),
    (r"ribose", "ribose"),
    (r"glucuronate|glucuronic", "glucuronate"),
    (r"galacturonate|galacturonic", "galacturonate"),
    (r"glycerol", "glycerol"),
    (r"mannitol", "mannitol"),
    (r"sorbitol", "sorbitol"),
    (r"N-acetylglucosamine|GlcNAc", "N-acetylglucosamine"),
    (r"oligopeptide", "oligopeptide"),
    (r"dipeptide", "dipeptide"),
    (r"tripeptide", "tripeptide"),
    (r"peptide", "peptide"),
    (r"amino acid", "amino acid"),
    (r"arginine", "arginine"),
    (r"glutamate|glutamic", "glutamate"),
    (r"glutamine", "glutamine"),
    (r"aspartate|aspartic", "aspartate"),
    (r"asparagine", "asparagine"),
    (r"lysine", "lysine"),
    (r"serine", "serine"),
    (r"threonine", "threonine"),
    (r"methionine", "methionine"),
    (r"cysteine", "cysteine"),
    (r"glycine", "glycine"),
    (r"proline", "proline"),
    (r"valine", "valine"),
    (r"leucine", "leucine"),
    (r"isoleucine", "isoleucine"),
    (r"histidine", "histidine"),
    (r"phenylalanine", "phenylalanine"),
    (r"tyrosine", "tyrosine"),
    (r"tryptophan", "tryptophan"),
    (r"alanine", "alanine"),
    (r"nucleoside", "nucleoside"),
    (r"nucleotide", "nucleotide"),
    (r"purine", "purine"),
    (r"pyrimidine", "pyrimidine"),
    (r"adenine", "adenine"),
    (r"guanine", "guanine"),
    (r"cytosine", "cytosine"),
    (r"uracil", "uracil"),
    (r"thymine", "thymine"),
    (r"xanthine", "xanthine"),
    (r"hypoxanthine", "hypoxanthine"),
    (r"adenosine", "adenosine"),
    (r"guanosine", "guanosine"),
    (r"uridine", "uridine"),
    (r"cytidine", "cytidine"),
    (r"biotin", "biotin"),
    (r"thiamine|thiamin", "thiamine"),
    (r"riboflavin", "riboflavin"),
    (r"cobalamin|vitamin B12", "cobalamin"),
    (r"folate", "folate"),
    (r"pantothenate", "pantothenate"),
    (r"nicotinamide|nicotinic", "nicotinamide"),
    (r"lactate", "lactate"),
    (r"acetate", "acetate"),
    (r"citrate", "citrate"),
    (r"succinate", "succinate"),
    (r"malate", "malate"),
    (r"fumarate", "fumarate"),
    (r"formate", "formate"),
    (r"oxalate", "oxalate"),
    (r"pyruvate", "pyruvate"),
    (r"2-oxoglutarate|alpha-ketoglutarate", "2-oxoglutarate"),
    (r"dicarboxylate", "dicarboxylate"),
    (r"tricarboxylate", "tricarboxylate"),
    (r"spermidine", "spermidine"),
    (r"putrescine", "putrescine"),
    (r"polyamine", "polyamine"),
    (r"choline", "choline"),
    (r"betaine", "betaine"),
    (r"carnitine", "carnitine"),
    (r"taurine", "taurine"),
    (r"multidrug", "multidrug"),
    (r"tetracycline", "tetracycline"),
    (r"macrolide", "macrolide"),
    (r"aminoglycoside", "aminoglycoside"),
    (r"beta[- ]lactam|β[- ]lactam", "beta-lactam"),
    (r"chloramphenicol", "chloramphenicol"),
    (r"fluoroquinolone", "fluoroquinolone"),
    (r"antibiotic", "antibiotic"),
    (r"\bdrug\b", "drug"),
    (r"arsenite", "arsenite"),
    (r"arsenate", "arsenate"),
    (r"chromate", "chromate"),
    (r"mercury", "mercury"),
    (r"cadmium", "cadmium"),
    (r"tellurite", "tellurite"),
    (r"colicin", "colicin"),
    (r"bacteriocin", "bacteriocin"),
    (r"phospholipid", "phospholipid"),
    (r"fatty acid", "fatty acid"),
    (r"lipopolysaccharide|\bLPS\b", "lipopolysaccharide"),
    (r"teichoic acid", "teichoic acid"),
    (r"\blipid\b", "lipid"),
]
SUBSTRATE_RE = [(re.compile(p, re.I), name) for p, name in SUBSTRATES]
TRANSPORT_TEXT_RE = re.compile(
    r"transport|permease|symport|antiport|uptake|ABC|channel|pump|binding protein|receptor",
    re.I)


def clean_seq(s: str) -> str:
    return re.sub(r"[^A-Za-z]", "", s).upper().rstrip("*")


def load_faa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = clean_seq("".join(chunks))
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if cur is not None:
        seqs[cur] = clean_seq("".join(chunks))
    return seqs


def load_gff(path: Path):
    contigs: dict[str, int] = {}
    feats = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    contigs[p[1]] = int(p[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            if f[2] != "CDS":
                continue
            attrs = dict()
            for part in f[8].split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    attrs[k] = urllib.parse.unquote(v)
            pid = attrs.get("ID") or attrs.get("protein_id")
            if not pid:
                continue
            feats.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]),
                              strand=f[6], protein_id=pid))
    return contigs, feats


def load_functional(path: Path) -> dict[str, dict]:
    best: dict[str, dict] = {}
    if not path.exists():
        return best
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 13:
                continue
            if row[0] in {"qseqid", "protein_id"}:
                continue
            qseqid, sseqid, pident, evalue, stitle = row[0], row[1], row[2], row[10], row[12]
            try:
                pident_f = float(pident)
                evalue_f = float(evalue)
            except ValueError:
                continue
            prev = best.get(qseqid)
            if prev is None or evalue_f < prev["evalue"]:
                best[qseqid] = dict(sseqid=sseqid, pident=pident_f, evalue=evalue_f, stitle=stitle)
    return best


def load_panel_models(path: Path) -> dict[str, dict]:
    models = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            models[row["pfam_name"]] = row
    return models


def parse_domtbl(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 23:
                continue
            try:
                rec = dict(
                    protein=f[3], pfam=f[0], pfam_acc=f[1],
                    full_e=float(f[6]), full_score=float(f[7]),
                    i_e=float(f[12]), dom_score=float(f[13]),
                    ali_from=int(f[17]), ali_to=int(f[18]), acc=float(f[21]),
                )
            except ValueError:
                continue
            hits[rec["protein"]].append(rec)
    return hits


def run_hmmscan(panel: Path, faa: Path, out: Path, cpu: int = 4) -> None:
    cmd = ["hmmscan", "--cpu", str(cpu), "--noali", "--cut_ga",
           "--domtblout", str(out), "-o", "/dev/null", str(panel), str(faa)]
    subprocess.run(cmd, check=True)


def mean_kd(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(KD.get(a, 0.0) for a in seq) / len(seq)


def tm_segments(seq: str, window: int = 19, threshold: float = 2.0,
                min_len: int = 18) -> list[dict]:
    n = len(seq)
    if n < window:
        return []
    kd = [KD.get(a, 0.0) for a in seq]
    # media movil centrada
    prefix = [0.0]
    for v in kd:
        prefix.append(prefix[-1] + v)
    hydro_centers = []
    half = window // 2
    for i in range(half, n - half):
        m = (prefix[i + half + 1] - prefix[i - half]) / window
        if m >= threshold:
            hydro_centers.append((i, m))
    if not hydro_centers:
        return []
    # fusionar centros contiguos (o separados por <=2)
    groups = []
    cur = [hydro_centers[0]]
    for item in hydro_centers[1:]:
        if item[0] - cur[-1][0] <= 2:
            cur.append(item)
        else:
            groups.append(cur)
            cur = [item]
    groups.append(cur)
    raw = []
    for g in groups:
        start = max(0, g[0][0] - half)
        end = min(n, g[-1][0] + half + 1)
        # recortar bordes poco hidrofobicos
        while start < end and kd[start] < 0:
            start += 1
        while end > start and kd[end - 1] < 0:
            end -= 1
        if end - start < min_len:
            continue
        m = mean_kd(seq[start:end])
        if m >= threshold:
            raw.append([start, end])
    # fusionar segmentos solapados o casi contiguos (ventanas centradas)
    merged: list[list[int]] = []
    for start, end in sorted(raw):
        if merged and start <= merged[-1][1] + 2:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    segs = []
    for start, end in merged:
        m = mean_kd(seq[start:end])
        if end - start < min_len or m < threshold:
            continue
        segs.append(dict(start=start + 1, end=end, length=end - start, mean=round(m, 2)))
    return segs


def sp_signal(seq: str) -> dict | None:
    """Heuristica de signal peptide (n-region + h-region + sitio de corte A-X-A).

    No sustituye SignalP. Solo se reporta como evidencia debil.
    """
    if len(seq) < 25:
        return None
    best = None
    for length in (9, 11, 13):
        for start in range(0, 21 - length + 1):
            w = seq[start:start + length]
            m = mean_kd(w)
            if m < 2.2:
                continue
            if not any(a in "LIVMFWAC" for a in w):
                continue
            if best is None or m > best[0]:
                best = (m, start, start + length)
    if best is None:
        return None
    m, start, end = best
    if end >= len(seq) - 1:
        return None
    # residuos -1 y -3 pequenos (posiciones end-1 y end-3, 0-based)
    small = set("ASGCTV")
    a = seq[end - 1] if end - 1 < len(seq) else ""
    b = seq[end - 3] if end - 3 >= 0 else ""
    if a not in small and b not in small:
        return None
    # n-region con al menos 1 residuo y sin demasiadas cargas en h-region
    if start < 1:
        return None
    return dict(start=start + 1, end=end, h_len=end - start, mean=round(m, 2),
                minus1=a, minus3=b)


def lipobox(seq: str) -> dict | None:
    if len(seq) < 8:
        return None
    for m in re.finditer(r"(?=([LVI][ASTM][AG]C))", seq[:40]):
        i = m.start(1)
        cpos = i + 4  # 1-based posicion de la C
        if i < 3 or cpos > 40:
            continue
        upstream = seq[max(0, i - 12):i]
        # la h-region puede ser corta: maxima media de una ventana de 7 aa
        best = max((mean_kd(upstream[j:j + 7])
                    for j in range(0, max(1, len(upstream) - 6))), default=0.0)
        if best < 1.8:
            continue
        return dict(cys=cpos, motif=seq[i:i + 4], upstream_mean=round(best, 2))
    return None


def find_substrate(stitle: str) -> list[str]:
    # quitar metadatos OS=/OX=/GN=... para no confundir nombres de organismo
    # (p. ej. "Icterohaemorrhagiae" contiene "haem") con el sustrato
    desc = re.split(r"\s+OS=", stitle, maxsplit=1)[0]
    found = []
    for rx, name in SUBSTRATE_RE:
        if rx.search(desc):
            if name not in found:
                found.append(name)
    return found


def choose_primary(hits: list[dict], models: dict[str, dict]) -> dict | None:
    if not hits:
        return None
    # prioridad por rol y luego por E-value
    prio = {"secretion": 0, "transporter": 1, "membrane": 2, "lipoprotein": 3, "sbp": 4, "signal": 5}
    def key(h):
        role = models.get(h["pfam"], {}).get("role", "zz")
        return (prio.get(role, 9), h["full_e"], -h["full_score"])
    return sorted(hits, key=key)[0]


def classify(seq: str, hits: list[dict], models: dict[str, dict]):
    tm = tm_segments(seq)
    sp = sp_signal(seq)
    lb = lipobox(seq)

    # descartar el primer segmento si es la h-region de un SP y no hay mas
    tm_eff = list(tm)
    if sp and tm_eff:
        s = tm_eff[0]
        if s["start"] <= sp["end"] + 2 and s["end"] <= 35 and len(tm_eff) == 1:
            tm_eff = []

    primary = choose_primary(hits, models)
    role = models.get(primary["pfam"], {}).get("role") if primary else None
    lip_hit = any(models.get(h["pfam"], {}).get("role") == "lipoprotein" for h in hits)
    sbp_hit = any(models.get(h["pfam"], {}).get("role") == "sbp" for h in hits)
    sec_hit = any(models.get(h["pfam"], {}).get("role") == "secretion" for h in hits)
    tr_hit = any(models.get(h["pfam"], {}).get("role") == "transporter" for h in hits)
    mem_hit = any(models.get(h["pfam"], {}).get("role") == "membrane" for h in hits)
    sig_hit = any(models.get(h["pfam"], {}).get("role") == "signal" for h in hits)

    ftype = None
    if sec_hit:
        ftype = "secretion"
        primary = choose_primary([h for h in hits if models.get(h["pfam"], {}).get("role") == "secretion"], models)
    elif lb and lip_hit:
        ftype = "lipoprotein"
    elif tr_hit:
        ftype = "transporter"
        primary = choose_primary([h for h in hits if models.get(h["pfam"], {}).get("role") == "transporter"], models)
    elif lb:
        # lipobox con region hidrofobica aguas arriba: senal fuerte de lipoproteina,
        # aunque la heuristica de SP no alcance el umbral.
        ftype = "lipoprotein"
    elif len(tm_eff) >= 2:
        ftype = "tm"
    elif sp and len(tm_eff) == 1:
        # SP + helice TM fuera de la h-region: ancla de membrana, no proteina exportada soluble
        ftype = "tm"
    elif sp:
        ftype = "signal_peptide"
    elif mem_hit:
        ftype = "tm"
    elif sbp_hit:
        # SBP sin senal: no es una proteina de membrana; se omite
        return None
    else:
        return None

    # familia
    family = ""
    if primary:
        m = models.get(primary["pfam"], {})
        family = f"{primary['pfam_acc'].split('.')[0]} {primary['pfam']}"
        if m.get("subgroup") and m["subgroup"] not in {"lipoprotein", "SBP", "lipobox"}:
            family += f" ({m['subgroup']}"
            if m.get("tc") and m["tc"] != ".":
                family += f"; TC {m['tc']}"
            family += ")"
    elif lb and ftype == "lipoprotein":
        family = "lipobox (heuristica)"
    elif ftype == "tm":
        family = "hidropatia KD (sin Pfam)"
    elif ftype == "signal_peptide":
        family = "sin Pfam (SP heuristica)"

    # evidencia
    ev = []
    if primary:
        ev.append(f"hmmscan:Pfam={primary['pfam_acc']}/{primary['pfam']},fullE={primary['full_e']:.1e},"
                  f"iE={primary['i_e']:.1e},cov={primary['acc']:.2f},role={models[primary['pfam']]['role']}")
    elif hits:
        ev.append(f"hmmscan:Pfam={hits[0]['pfam_acc']}/{hits[0]['pfam']},fullE={hits[0]['full_e']:.1e}")
    if len(tm_eff) > 0:
        ev.append("kd_tm=" + ",".join(f"{s['start']}-{s['end']}({s['mean']})" for s in tm_eff))
    if sp:
        ev.append(f"sp_heur=nterm{sp['start']}-{sp['end']}(h={sp['h_len']},kd={sp['mean']},-1={sp['minus1']},-3={sp['minus3']})")
    if lb:
        ev.append(f"lipobox={lb['motif']}@C{lb['cys']}(up_kd={lb['upstream_mean']})")
    if len(tm) > len(tm_eff):
        ev.append("tm_hidrofobico_1_descartado_por_sp")
    return dict(type=ftype, family=family, evidence=ev, tm=tm_eff, sp=sp, lb=lb,
                primary=primary, sbp=sbp_hit, lip_hit=lip_hit, sec_hit=sec_hit)


def confidence(rec: dict, models: dict[str, dict]) -> str:
    p = rec.get("primary")
    ftype = rec["type"]
    e = p["full_e"] if p else 1.0
    cov = p["acc"] if p else 0.0
    role = models.get(p["pfam"], {}).get("role") if p else None

    if ftype == "secretion":
        if e <= 1e-20 and cov >= 0.5:
            return "high"
        if e <= 1e-5:
            return "medium"
        return "low"
    if ftype == "transporter":
        if e <= 1e-20 and (len(rec["tm"]) >= 1 or role == "transporter"):
            return "high"
        if e <= 1e-5:
            return "medium"
        return "low"
    if ftype == "lipoprotein":
        if rec["lb"] and rec["lip_hit"] and e <= 1e-5:
            return "high" if e <= 1e-20 else "medium"
        if rec["lb"] and rec["sp"]:
            return "medium"
        return "low"
    if ftype == "tm":
        if len(rec["tm"]) >= 4:
            return "medium"
        return "low"
    if ftype == "signal_peptide":
        if rec["sp"] and rec["sbp"] and e <= 1e-10:
            return "medium"
        if rec["sp"] and e <= 1e-10:
            return "medium"
        return "low"
    return "low"


def build_evidence(rec: dict, sub_info: dict | None, extra: list[str]) -> str:
    ev = list(rec["evidence"])
    if sub_info:
        ev.append(f"swissprot={sub_info['sseqid']},pident={sub_info['pident']:.1f},"
                  f"E={sub_info['evalue']:.1e}")
        ev.append(f"substrate_text={sub_info['substrates']}")
    ev.extend(extra)
    return ";".join(ev)


def gff_escape(v: str) -> str:
    safe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:_-+/"
    return urllib.parse.quote(v, safe=safe)


def annotate(genome: str, annot_dir: Path, panel: Path, models: dict[str, dict],
             out_tsv: Path, out_gff: Path, cpu: int = 4) -> dict:
    contigs, feats = load_gff(annot_dir / "annotation.gff3")
    seqs = load_faa(annot_dir / "proteins.faa")
    func = load_functional(annot_dir / "functional.tsv")

    with tempfile.TemporaryDirectory() as td:
        domtbl = Path(td) / "hits.domtbl"
        run_hmmscan(panel, annot_dir / "proteins.faa", domtbl, cpu=cpu)
        hits_by_prot = parse_domtbl(domtbl)

    rows = []
    for f in feats:
        pid = f["protein_id"]
        seq = seqs.get(pid, "")
        if not seq:
            continue
        hits = hits_by_prot.get(pid, [])
        rec = classify(seq, hits, models)
        if rec is None:
            continue
        ft = rec["type"]

        # sustrato: solo transporters y SBPs (signal_peptide), solo Swiss-Prot explicito
        substrate = ""
        sub_info = None
        fh = func.get(pid)
        if fh and fh["pident"] >= 50.0 and fh["evalue"] <= 1e-5:
            subs = find_substrate(fh["stitle"])
            # solo transporte/SBP y solo si el texto tiene contexto de transporte;
            # evita enzimas con "phosphate"/"heme" en el nombre
            transport_context = bool(TRANSPORT_TEXT_RE.search(fh["stitle"]))
            allowed = ((ft == "transporter" and transport_context)
                       or (ft == "signal_peptide" and rec["sbp"]))
            if subs and allowed:
                substrate = "|".join(subs[:2])
                sub_info = dict(fh)
                sub_info["substrates"] = substrate
        if substrate and sub_info is None:
            substrate = ""

        # refuerzo de lipoproteina por texto Swiss-Prot
        extra = []
        if fh and fh["pident"] >= 50 and re.search(r"lipoprotein", fh["stitle"], re.I):
            extra.append(f"swissprot_lipoprotein_text={fh['sseqid']},pident={fh['pident']:.1f}")
            if ft == "lipoprotein" and not rec["lip_hit"]:
                rec["type"] = "lipoprotein"
        # familia con TC para SEC/SBP cuando aplica
        if not rec["family"] and rec["type"] == "lipoprotein":
            rec["family"] = "lipobox (heuristica)"

        ev = build_evidence(rec, sub_info, extra)
        conf = confidence(rec, models)

        rows.append(dict(
            genome=genome, contig=f["contig"], start=f["start"], end=f["end"],
            strand=f["strand"], protein_id=pid, type=ft,
            family=rec["family"], substrate=substrate, evidence=ev, confidence=conf,
        ))

    # ordenar por contig (orden de aparicion) y coordenada
    order = {c: i for i, c in enumerate(contigs)}
    rows.sort(key=lambda r: (order.get(r["contig"], 9999), r["start"], r["end"], r["protein_id"], r["type"]))

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_gff.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["genome", "contig", "start", "end", "strand",
                                       "protein_id", "type", "family", "substrate",
                                       "evidence", "confidence"])
        w.writeheader()
        w.writerows(rows)

    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for c, ln in contigs.items():
            fh.write(f"##sequence-region {c} 1 {ln}\n")
        for i, r in enumerate(rows, 1):
            rid = f"{genome}_{r['protein_id']}_{r['type']}_{i}"
            attrs = [
                f"ID={gff_escape(rid)}",
                f"protein_id={gff_escape(r['protein_id'])}",
                f"membrane_type={gff_escape(r['type'])}",
                f"family={gff_escape(r['family'])}",
                f"substrate={gff_escape(r['substrate'])}",
                f"evidence={gff_escape(r['evidence'])}",
                f"confidence={gff_escape(r['confidence'])}",
            ]
            fh.write("\t".join([
                r["contig"], "membrane-builder", "polypeptide",
                str(r["start"]), str(r["end"]), ".", r["strand"], ".", ";".join(attrs)
            ]) + "\n")

    return dict(rows=rows, contigs=contigs, n_proteins=len(seqs), n_features=len(feats))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, help="nombre del directorio GCF_...")
    ap.add_argument("--annot-root", default=str(ROOT / "data/annotations_ab_initio"))
    ap.add_argument("--panel-dir", default=str(ROOT / "data/membrane/pfam_panel"))
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane"))
    ap.add_argument("--gff-dir", default=None,
                    help="directorio de GFF3 (por defecto <outdir>/gff)")
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()

    annot_dir = Path(args.annot_root) / args.genome
    if not annot_dir.is_dir():
        sys.exit(f"[annot] no existe {annot_dir}")
    panel = Path(args.panel_dir) / "membrane_panel.hmm"
    models = load_panel_models(Path(args.panel_dir) / "panel_models.tsv")
    out_tsv = Path(args.outdir) / f"{args.genome}.tsv"
    if args.gff_dir:
        out_gff = Path(args.gff_dir) / f"{args.genome}.gff3"
    else:
        out_gff = Path(args.outdir) / "gff" / f"{args.genome}.gff3"

    res = annotate(args.genome, annot_dir, panel, models, out_tsv, out_gff, cpu=args.cpu)
    counts = defaultdict(int)
    for r in res["rows"]:
        counts[r["type"]] += 1
    print(json.dumps(dict(genome=args.genome, contigs=len(res["contigs"]),
                          proteins=res["n_proteins"], rows=len(res["rows"]),
                          counts=counts), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
