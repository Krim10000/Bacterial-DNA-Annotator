#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_validate_markers.py — FASE A del gate de ORCA: validación de marcadores dnaA/dnaN.

Separación de responsabilidades (corrige el diseño acoplado del 35):

  A (este script): busca y VALIDA los marcadores dnaA/dnaN por genoma usando la
      anotación RefSeq YA DESCARGADA (data/annotations/*.gff.gz + *.faa.gz) y el
      genoma crudo (data/raw/*_genomic.fna.gz). Tolerante a numeración, guiones,
      guiones bajos y mayúsculas (dnaA, dnaa1, dnaA_2, dnaa-3...). Verifica:
        - unicidad del marcador (paralogs -> ambiguous, se listan todos)
        - evidencia gene= (curada) vs product= (solo texto -> weak)
        - posición interior (no pegada a bordes de contig: artefacto ya diagnosticado)
        - coherencia de operón: distancia dnaA<->dnaN
        - pseudogenes / CDS parciales
      Recupera NOMBRE y SECUENCIA (nt del genoma crudo con revcomp; aa del faa
      RefSeq por protein_id) y lista TODAS las alternativas encontradas.
      NO predice nada: solo emite veredictos.

  B (35_orca_from_gbk.py): consume los veredictos y corre ORCA SOLO sobre
      genomas con status=validated.

Estatus por genoma:
  validated   ambos marcadores únicos por gene=, interiores, no pseudogenes
  weak        marcador solo por product= (o pseudogénero/parcial) — revisar a mano
  ambiguous   >1 candidato gene= para algún marcador (se listan todos)
  missing     falta dnaA o dnaN
  error       problema de lectura del GFF (se reporta y continúa)

Todo offline: no usa red. Fuentes:
  data/annotations/{GCF}*.gff.gz   CDS con gene=/product=/Dbxref (EC:...)
  data/annotations/{GCF}*.faa.gz   proteínas RefSeq (secuencia aa)
  data/raw/{GCF}*_genomic.fna.gz   genoma crudo (secuencia nt por replicón)

Salidas:
  data/marker_validation/{GCF}.json             detalle: marcadores, secuencias, alternativas
  data/marker_validation/validation_summary.json
  data/marker_validation/validated_genomes.txt  (entradas para la fase B)

Uso:
  python scripts/oric/01_validate_markers.py --only GCF_000006685.1 -v
  python scripts/oric/01_validate_markers.py                 # los 400 (~10-20 min)
  python scripts/oric/01_validate_markers.py --selftest      # prueba autónoma
"""

import argparse
import gzip
import json
import logging
import re
import time
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

LOG = logging.getLogger("35a")

# EC en Dbxref (EC:3.6.1.-) o dentro del producto ("(EC 2.7.1.1)").
# Acepta clases con guion (clase/subclase indefinida) y 1-3 subclases.
EC_RE = re.compile(r"(?<![A-Za-z])EC[:=]\s*([0-9\-]+(?:\.[0-9\-]+){1,3})", re.I)

MARKERS = {
    "dnaA": {
        "gene_re": r"dnaa[-_]?\d*",   # dnaA, dnaa, dnaA1, dnaa_2, dnaa-3
        "product_re": r"(?i)(chromosomal\s+replication\s+initiator|\bdnaa\b)",
    },
    "dnaN": {
        "gene_re": r"dnan[-_]?\d*",   # dnaN, dnan, dnan1, dnan_2
        "product_re": r"(?i)(sliding\s+clamp|dna\s+polymerase\s+iii.*\bbeta\b)",
    },
}

_COMP_DNA = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "K": "M", "M": "K",
    "S": "S", "W": "W", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}
_REVCOMP_TABLE = str.maketrans(
    {**_COMP_DNA, **{k.lower(): v.lower() for k, v in _COMP_DNA.items()}}
)


def revcomp(seq):
    return seq.translate(_REVCOMP_TABLE)[::-1]


# ----------------------------------------------------------------------------- IO

def _read_text(path):
    """Lee un archivo (plain o .gz) como texto decodificado manualmente.

    Evita gzip.open(..., 'rt', errors=...): algunas builds rechazan
    encoding/errors en el modo texto de gzip; binario + decode funciona igual
    en todas partes.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        return fh.read().decode("utf-8", errors="replace")


def parse_attrs(s):
    """Parsea atributos GFF; NCBI usa percent-encoding (%20 -> espacio)."""
    d = {}
    for part in s.strip().rstrip(";").split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        d[k.strip()] = unquote(v.strip()).strip('"').strip()
    return d


def parse_gff_cds(path):
    """CDS del GFF RefSeq. Deduplica por ID (CDS multi-segmento, raro en procariotas)."""
    cds, seen = [], set()
    for line in _read_text(path).splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        a = parse_attrs(p[8])
        fid = a.get("ID") or "noid-%d" % len(cds)
        if fid in seen:
            continue
        seen.add(fid)
        try:
            start, end = int(p[3]), int(p[4])
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        haystack = (a.get("Dbxref", "") + ";" + a.get("product", "") + ";" +
                    a.get("Note", ""))
        cds.append({
            "replicon": p[0],
            "start": start, "end": end, "strand": p[6],
            "gene": a.get("gene"), "product": a.get("product"),
            "protein_id": a.get("protein_id") or a.get("Name"),
            "locus_tag": a.get("locus_tag"),
            "ec": sorted(set(EC_RE.findall(haystack))),
            "pseudogene": a.get("pseudogene") == "true" or a.get("pseudo") == "true",
            "partial": a.get("partial") == "true",
        })
    return cds


def index_fasta(path):
    """FASTA (faa o fna, plain o .gz) -> {acc_versionado: seq}, {acc_sin_version: acc}."""
    seqs, acc, buf = {}, None, []
    for line in _read_text(path).splitlines():
        line = line.strip()
        if line.startswith(">"):
            if acc is not None:
                seqs[acc] = "".join(buf)
            acc = line[1:].split()[0]
            buf = []
        elif acc is not None and line:
            buf.append(line)
    if acc is not None:
        seqs[acc] = "".join(buf)
    base = {}
    for k in seqs:
        b = k.split(".")[0]
        base.setdefault(b, k)
    return seqs, base


def find_fna(raw_dir, stem):
    """Localiza el genoma crudo: prefiere {stem}_genomic.fna.gz."""
    raw_dir = Path(raw_dir)
    if not raw_dir.is_dir():
        return None
    hits = sorted(raw_dir.glob("%s*.fna.gz" % stem)) or sorted(raw_dir.glob("%s*.gz" % stem))
    if not hits:
        return None
    pref = [h for h in hits if h.name.endswith("_genomic.fna.gz")]
    return (pref or hits)[0]


# ------------------------------------------------------------------- marcadores

def find_candidates(cds, marker):
    """Todos los CDS que casan el marcador. Evidencia 'gene' (curada) > 'product'."""
    spec = MARKERS[marker]
    out = []
    for f in cds:
        gene = (f.get("gene") or "").strip()
        prod = (f.get("product") or "").strip()
        ev = None
        if gene and re.fullmatch(spec["gene_re"], gene, re.I):
            ev = "gene"
        if ev is None and prod and re.search(spec["product_re"], prod):
            ev = "product"
        if ev:
            out.append(dict(f, evidence=ev))
    return out


def _pick_primary(cands, gfa):
    """Entre varios candidatos: replicón mayor primero, luego el más 5'."""
    def key(c):
        seq = gfa.get(c["replicon"]) or gfa.get(c["replicon"].split(".")[0]) or ""
        return (-len(seq), c["start"])
    return sorted(cands, key=key)[0]


def _compact(c):
    return {
        "gene": c.get("gene"), "product": c.get("product"),
        "protein_id": c.get("protein_id"), "replicon": c["replicon"],
        "start": c["start"], "end": c["end"], "strand": c["strand"],
        "evidence": c["evidence"], "ec": c.get("ec") or [],
        "pseudogene": c.get("pseudogene", False), "partial": c.get("partial", False),
    }


def process_genome(stem, paths, edge_margin=1000, max_distance=20000):
    rec = {
        "genome": stem, "status": "error", "reasons": [], "notes": [],
        "markers": {}, "alternatives": {},
        "distance_dnaA_dnaN_bp": None, "operon_consistent": None,
        "sources": {k: str(v) for k, v in paths.items() if v},
        "schema": "marker_validation_v1",
    }
    try:
        cds = parse_gff_cds(paths["gff"])
    except Exception as e:
        rec["reasons"].append("gff_error:%s" % e)
        return rec
    if not cds:
        rec["reasons"].append("gff_sin_CDS")
        return rec

    fa, faa_base, gfa = {}, {}, {}
    if paths.get("faa"):
        try:
            fa, faa_base = index_fasta(paths["faa"])
        except Exception as e:
            rec["notes"].append("faa_error:%s" % e)
    if paths.get("fna"):
        try:
            gfa, _ = index_fasta(paths["fna"])
        except Exception as e:
            rec["notes"].append("fna_error:%s" % e)

    miss_any = amb_any = weak_any = False
    primaries = {}

    for marker in ("dnaA", "dnaN"):
        cands = find_candidates(cds, marker)
        gene_c = [c for c in cands if c["evidence"] == "gene"]
        prod_c = [c for c in cands if c["evidence"] == "product"]
        primary_raw, weak = None, False

        if len(gene_c) == 1:
            primary_raw = gene_c[0]
        elif len(gene_c) > 1:
            primary_raw = _pick_primary(gene_c, gfa)
        elif len(prod_c) == 1:
            primary_raw, weak = prod_c[0], True
        elif len(prod_c) > 1:
            primary_raw, weak = _pick_primary(prod_c, gfa), True

        if primary_raw is None:
            miss_any = True
            rec["reasons"].append("missing:%s" % marker)
            rec["alternatives"][marker] = []
            continue

        pid = primary_raw.get("protein_id") or ""
        aa = fa.get(pid) or (faa_base.get(pid.split(".")[0]) if pid else None)
        rep = primary_raw["replicon"]
        gseq = gfa.get(rep) or gfa.get(rep.split(".")[0])
        rep_len = len(gseq) if gseq else 0
        nt = None
        if gseq:
            nt = gseq[primary_raw["start"] - 1:primary_raw["end"]]
            if primary_raw["strand"] == "-":
                nt = revcomp(nt)
        interior = None
        if rep_len:
            interior = (primary_raw["start"] > edge_margin and
                        (rep_len - primary_raw["end"]) > edge_margin)

        primary = dict(primary_raw, aa_seq=aa,
                       aa_len=(len(aa) if aa else None),
                       nt_seq=nt, replicon_len=(rep_len or None),
                       interior=interior)
        primaries[marker] = primary

        rec["alternatives"][marker] = [
            _compact(c) for c in (gene_c + prod_c) if c is not primary_raw
        ]
        rec["markers"][marker] = {
            "gene": primary_raw.get("gene"),
            "product": primary_raw.get("product"),
            "protein_id": pid or None,
            "locus_tag": primary_raw.get("locus_tag"),
            "replicon": rep,
            "start": primary_raw["start"], "end": primary_raw["end"],
            "strand": primary_raw["strand"],
            "length_bp": primary_raw["end"] - primary_raw["start"] + 1,
            "aa_len": primary["aa_len"],
            "ec": primary_raw.get("ec") or [],
            "evidence": primary_raw["evidence"],
            "pseudogene": primary_raw.get("pseudogene", False),
            "partial": primary_raw.get("partial", False),
            "interior": interior,
            "replicon_len": primary["replicon_len"],
            "nt_seq": nt, "aa_seq": aa,
        }

        if primary_raw.get("pseudogene"):
            weak = True
            rec["reasons"].append("pseudogene:%s" % marker)
        if primary_raw.get("partial"):
            weak = True
            rec["reasons"].append("partial_cds:%s" % marker)
        if interior is False:
            rec["reasons"].append("edge_proximity:%s" % marker)
        if weak:
            weak_any = True
            if primary_raw["evidence"] != "gene":
                rec["reasons"].append("weak_evidence:%s" % marker)
        if len(gene_c) > 1:
            amb_any = True
            rec["reasons"].append("ambiguous:%s(%d)" % (marker, len(gene_c)))

    a, n = primaries.get("dnaA"), primaries.get("dnaN")
    if a and n:
        if a["replicon"] == n["replicon"]:
            d = abs((a["start"] + a["end"]) // 2 - (n["start"] + n["end"]) // 2)
            rec["distance_dnaA_dnaN_bp"] = d
            rec["operon_consistent"] = d <= max_distance
            if d > max_distance:
                rec["reasons"].append("distance_warn:%dkb" % (d // 1000))
        else:
            rec["reasons"].append("different_replicons")
            weak_any = True

    if miss_any:
        rec["status"] = "missing"
    elif amb_any:
        rec["status"] = "ambiguous"
    elif weak_any:
        rec["status"] = "weak"
    else:
        rec["status"] = "validated"
    return rec


# ------------------------------------------------------------------ descubrir

def discover(ann_dir, only=None, limit=None):
    ann_dir = Path(ann_dir)
    out = []
    for gff in sorted(ann_dir.glob("*.gff.gz")):
        stem = gff.name[:-len(".gff.gz")]
        if only and only not in stem:
            continue
        faa = ann_dir / (stem + ".faa.gz")
        out.append({
            "stem": stem,
            "gff": gff,
            "faa": faa if faa.exists() else None,
        })
    if limit:
        out = out[:int(limit)]
    return out


def run_validation(ann_dir, raw_dir, out_dir, only=None, limit=None,
                   edge_margin=1000, max_distance=20000, verbose=False):
    items = discover(ann_dir, only, limit)
    if not items:
        LOG.warning("Sin GFFs en %s (¿--only mal escrito? ¿directorio vacío?)", ann_dir)
        return []
    LOG.info("=== 35A: validación de marcadores dnaA/dnaN (%d genomas) ===", len(items))
    LOG.info("GFF+faa: %s | fasta crudo: %s | salida: %s", ann_dir, raw_dir, out_dir)
    records, t0 = [], time.time()
    for i, it in enumerate(items, 1):
        fna = find_fna(raw_dir, it["stem"])
        rec = process_genome(it["stem"], {"gff": it["gff"], "faa": it["faa"],
                                          "fna": fna}, edge_margin, max_distance)
        records.append(rec)
        m = rec["markers"]

        def _fmt(mk):
            if mk not in m:
                return "%s:—" % mk
            d = m[mk]
            label = d["gene"] or (d["product"] or "")[:32]
            return "%s %s:%d-%d(%s) %s aa=%s" % (mk, d["replicon"], d["start"],
                                                 d["end"], d["strand"], label,
                                                 d["aa_len"] if d["aa_len"] else "?")

        dist = rec["distance_dnaA_dnaN_bp"]
        dist_s = ("%.1f kb" % (dist / 1000)) if dist is not None else "—"
        line = ("[%d/%d] %s: %s | %s | %s | dist=%s | alt dnaA=%d dnaN=%d" %
                (i, len(items), rec["genome"], rec["status"], _fmt("dnaA"),
                 _fmt("dnaN"), dist_s,
                 len(rec["alternatives"].get("dnaA", [])),
                 len(rec["alternatives"].get("dnaN", []))))
        (LOG.info if rec["status"] == "validated" else LOG.warning)(line)
        if verbose and rec["status"] != "validated":
            LOG.info("    razones: %s", "; ".join(rec["reasons"]) or "—")
        if verbose:
            for mk, d in m.items():
                if d.get("nt_seq"):
                    LOG.info("    %s nt[%d]: %s...", mk, len(d["nt_seq"]),
                             d["nt_seq"][:60])
                if d.get("aa_seq"):
                    LOG.info("    %s aa[%d]: %s...", mk, len(d["aa_seq"]),
                             d["aa_seq"][:60])
    LOG.info("=== %d genomas en %.1f s ===", len(records), time.time() - t0)
    write_outputs(out_dir, records)
    counts = Counter(r["status"] for r in records)
    LOG.info("=== Resumen 35A ===")
    for st in ("validated", "weak", "ambiguous", "missing", "error"):
        LOG.info("  %-10s %d", st, counts.get(st, 0))
    LOG.info("  Salidas: %s/{{genoma}.json, validation_summary.json, validated_genomes.txt}",
             out_dir)
    return records


def write_outputs(out_dir, records):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for r in records:
        (out_dir / (r["genome"] + ".json")).write_text(
            json.dumps(r, ensure_ascii=False, indent=2))
    counts = Counter(r["status"] for r in records)
    summary = {
        "schema": "marker_validation_summary_v1",
        "counts": {k: counts.get(k, 0) for k in
                   ("validated", "weak", "ambiguous", "missing", "error")},
        "genomes": [{"genome": r["genome"], "status": r["status"],
                     "distance_dnaA_dnaN_bp": r["distance_dnaA_dnaN_bp"],
                     "reasons": r["reasons"]} for r in records],
    }
    (out_dir / "validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    (out_dir / "validated_genomes.txt").write_text(
        "".join(r["genome"] + "\n" for r in records if r["status"] == "validated"))


# ------------------------------------------------------------------- selftest

def _selftest():
    import random
    import shutil
    root = Path("data/_selftest_35a")
    shutil.rmtree(root, ignore_errors=True)
    ann, raw, out = root / "ann", root / "raw", root / "out"
    for d in (ann, raw):
        d.mkdir(parents=True)
    rng = random.Random(42)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def w_gz(path, text):
        with gzip.open(path, "wt") as fh:
            fh.write(text)

    def cds_line(seqid, s, e, strand, attrs):
        return "%s\tRefSeq\tCDS\t%d\t%d\t.\t%s\t0\t%s\n" % (seqid, s, e, strand, attrs)

    def w_fasta_gz(path, records):
        with gzip.open(path, "wt") as fh:
            for h, s in records:
                fh.write(">" + h + "\n")
                for i in range(0, len(s), 70):
                    fh.write(s[i:i + 70] + "\n")

    # --- genoma A: limpio -> validated (dnaN en strand '-' para probar revcomp)
    accA, seqidA = "GCF_999999999.1_TESTA", "NC_999001.2"
    contig = list(rnd(120000))
    dnaa_nt = "ATG" + rnd(1359) + "TAA"          # 1365 nt -> 455 aa
    contig[10000:11365] = list(dnaa_nt)
    dnan_fwd = "ATG" + rnd(1095) + "TAA"         # 1101 nt -> 366 aa, instalado en '-'
    contig[12000:13101] = list(revcomp(dnan_fwd))
    contigA = "".join(contig)
    aaA, aaN = "M" + "K" * 454, "M" + "A" * 365
    gffA = "#!genome-build TESTA\n"
    gffA += cds_line(seqidA, 10001, 11365, "+",
                     "ID=cds-1;Parent=gene-1;Dbxref=Genbank:WP_TEST001.1,EC:3.6.1.-;"
                     "Name=WP_TEST001.1;gene=dnaA;"
                     "product=Chromosomal%20replication%20initiator%20protein%20DnaA")
    gffA += cds_line(seqidA, 50001, 50900, "+",
                     "ID=cds-9;Name=WP_TEST009.1;"
                     "product=DnaA-binding%20transcriptional%20regulator")
    gffA += cds_line(seqidA, 12001, 13101, "-",
                     "ID=cds-2;Dbxref=Genbank:WP_TEST002.1;Name=WP_TEST002.1;gene=dnaN;"
                     "product=DNA%20polymerase%20sliding%20clamp")
    w_gz(ann / (accA + ".gff.gz"), gffA)
    w_gz(ann / (accA + ".faa.gz"),
         ">WP_TEST001.1 chromosomal replication initiator protein DnaA [Test]\n" +
         "\n".join(aaA[i:i + 70] for i in range(0, len(aaA), 70)) + "\n" +
         ">WP_TEST002.1 DNA polymerase sliding clamp [Test]\n" +
         "\n".join(aaN[i:i + 70] for i in range(0, len(aaN), 70)) + "\n")
    w_gz(raw / (accA + "_genomic.fna.gz"),
         ">" + seqidA + " Test bacterium chromosome\n" +
         "\n".join(contigA[i:i + 70] for i in range(0, len(contigA), 70)) + "\n")

    # --- genoma B: dnaA duplicado (dnaa + dnaA2) y sin dnaN -> ambiguous
    accB, seqidB = "GCF_999999998.1_TESTB", "NC_999002.1"
    contigB = list(rnd(100000))
    dnaab = "ATG" + rnd(894) + "TAA"             # 900 nt
    contigB[7999:8999] = list(dnaab)
    dnaa2 = "ATG" + rnd(994) + "TAA"             # 1000 nt
    contigB[90000:91000] = list(dnaa2)
    contigB = "".join(contigB)
    gffB = "#!genome-build TESTB\n"
    gffB += cds_line(seqidB, 8001, 9000, "+",
                     "ID=cds-1;Name=WP_TEST003.1;gene=dnaA;"
                     "product=Chromosomal%20replication%20initiator%20protein%20DnaA")
    gffB += cds_line(seqidB, 90001, 91000, "+",
                     "ID=cds-2;Name=WP_TEST004.1;gene=dnaA2;"
                     "product=Chromosomal%20replication%20initiator%20protein%20DnaA")
    w_gz(ann / (accB + ".gff.gz"), gffB)
    w_gz(ann / (accB + ".faa.gz"),
         ">WP_TEST003.1 chromosomal replication initiator protein DnaA [Test]\n" +
         "\n".join(("M" + "K" * 299)[i:i + 70] for i in range(0, 300, 70)) + "\n" +
         ">WP_TEST004.1 chromosomal replication initiator protein DnaA [Test]\n" +
         "\n".join(("M" + "K" * 299)[i:i + 70] for i in range(0, 300, 70)) + "\n")
    w_gz(raw / (accB + "_genomic.fna.gz"),
         ">" + seqidB + " Test bacterium chromosome\n" +
         "\n".join(contigB[i:i + 70] for i in range(0, len(contigB), 70)) + "\n")

    try:
        recs = run_validation(ann, raw, out, verbose=True)
        A = next(r for r in recs if r["genome"] == accA)
        B = next(r for r in recs if r["genome"] == accB)

        assert A["status"] == "validated", A["status"]
        mA, mN = A["markers"]["dnaA"], A["markers"]["dnaN"]
        assert len(mA["nt_seq"]) == 1365 and mA["nt_seq"].startswith("ATG")
        assert mA["nt_seq"].endswith("TAA")
        assert len(mA["aa_seq"]) == 455
        assert mA["ec"] == ["3.6.1.-"], mA["ec"]          # EC con guion: parsea bien
        assert len(mN["nt_seq"]) == 1101 and mN["nt_seq"].startswith("ATG")
        assert mN["nt_seq"].endswith("TAA")               # revcomp correcto en strand -
        assert len(mN["aa_seq"]) == 366 and mN["strand"] == "-"
        assert A["distance_dnaA_dnaN_bp"] == 1868, A["distance_dnaA_dnaN_bp"]
        assert A["operon_consistent"] is True
        assert mA["interior"] is True and mN["interior"] is True
        altA = A["alternatives"]["dnaA"]
        assert len(altA) == 1 and altA[0]["evidence"] == "product", altA
        assert "DnaA-binding" in (altA[0]["product"] or "")

        assert B["status"] == "missing", B["status"]  # missing prevalece sobre ambiguous
        assert "ambiguous:dnaA(2)" in B["reasons"], B["reasons"]
        assert "missing:dnaN" in B["reasons"]
        altsB = B["alternatives"]["dnaA"]
        assert any(c["gene"] == "dnaA2" for c in altsB), altsB

        txt = (out / "validated_genomes.txt").read_text().split()
        assert txt == [accA], txt
        ssum = json.loads((out / "validation_summary.json").read_text())
        assert ssum["counts"]["validated"] == 1 and ssum["counts"]["missing"] == 1
        print("SELFTEST 35A: PASS (validated, revcomp, EC con guion, alternativas, "
              "ambiguous+dnaA2, missing:dnaN, salidas)")
        return 0
    finally:
        shutil.rmtree(root, ignore_errors=True)


# ----------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="FASE A: validación de marcadores dnaA/dnaN "
                    "(entradas de 35B/ORCA; no predice nada)")
    ap.add_argument("--ann-dir", default="data/annotations",
                    help="GFF+faa RefSeq descargados (default: data/annotations)")
    ap.add_argument("--raw-dir", default="data/raw",
                    help="Genomas crudos *_genomic.fna.gz (default: data/raw)")
    ap.add_argument("--out-dir", default="data/marker_validation")
    ap.add_argument("--only", help="subcadena del accession (p.ej. GCF_000006685.1)")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--edge-margin", type=int, default=1000,
                    help="bp mínimas al borde del contig (default 1000)")
    ap.add_argument("--max-distance", type=int, default=20000,
                    help="distancia máx dnaA<->dnaN para operón consistente (bp)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="muestra secuencias recuperadas y razones detalladas")
    ap.add_argument("--selftest", action="store_true",
                    help="prueba autónoma con genomas sintéticos y sale")
    args = ap.parse_args(argv)
    logging.basicConfig(format="%(asctime)s | %(levelname)-7s | %(message)s",
                        datefmt="%H:%M:%S", level=logging.INFO)
    if args.selftest:
        return _selftest()
    run_validation(args.ann_dir, args.raw_dir, args.out_dir,
                   only=args.only, limit=args.limit,
                   edge_margin=args.edge_margin,
                   max_distance=args.max_distance, verbose=args.verbose)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
