#!/usr/bin/env python3
"""Prototipo end-to-end: genoma crudo → proteínas → represores → sitios → efector.

Flujo:
  1. Llamada de genes con pyrodigal (sin dependencias externas).
  2. Detección de dominios de TF con hmmsearch sobre el panel de 299 HMM.
  3. Familia asignada por el dominio; homólogo más cercano por alineamiento contra los
     243 TFs de referencia de RegulonDB.
  4. Sitios de unión: se escanea la región aguas arriba (−250 pb) con los PWM de la
     familia; percentil frente a ventanas barajadas del propio genoma.
  5. Efector: se **hereda** del homólogo si la identidad ≥ 40 % (`inherited`), o se usa el
     prior de familia (`family_prior`), o queda `unknown`.

Todas las salidas llevan `status` y `evidence`: nada se declara confirmado.

Uso:
    python scripts/allostery/18_genome_pipeline.py genoma.fna --out data/allostery/genome_run
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "scripts" / "allostery"
sys.path.insert(0, str(X))
from _common import FAMILY_MAP  # noqa: E402

RD = ROOT / "data" / "regulondb"
HMM = ROOT / "data" / "allostery" / "work" / "tf_domains.hmm"
PWM = RD / "gt" / "tf_pwms.tsv"
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = str.maketrans("ACGT", "TGCA")
CORE_W = 14
UPSTREAM = 250

# Códigos de familia de genes.tsv (columna `fam`) -> familia canónica
FAM_CODE_MAP = {
    "LACI": "LacI/GalR", "GALS": "LacI/GalR", "PURR": "LacI/GalR", "CYTR": "LacI/GalR",
    "LYSR": "LysR", "GNTR": "GntR", "ARAC": "AraC", "MARR": "MarR", "TETR": "TetR",
    "DEOR": "DeoR", "DEORC": "DeoR", "ARSR": "ArsR/SmtB", "MERR": "MerR", "FUR": "Fur",
    "ICLR": "IclR", "RPIR": "RpiR", "ASNC": "AsnC", "PADR": "PadR", "RRF2": "Rrf2",
    "NIKR": "NikR", "CSOR": "CsoR", "METJ": "MetJ", "NRDR": "NrdR", "HXLR": "HxlR",
    "TRPR": "TrpR", "DtxR": "DtxR", "FE_DEP_REPRESS": "DtxR", "COPR": "CopY",
}


def fam_code_to_family(code: str) -> str:
    return FAM_CODE_MAP.get(code.upper(), code)


DBD_OK = ("HTH", "LacI", "TetR_N", "MerR", "DtxR", "Fur", "PadR", "Rrf2", "NrdR",
          "MetJ", "HxlR", "NikR", "CopY", "CsoR", "MarR", "AsnC", "IclR",
          "RpiR", "SorC", "Trp_repressor", "LysR", "DeoR", "UTRA", "AraC", "GntR",
          "AbrB", "Crp", "Fnr", "OmpR", "NarL", "LuxR", "Xre", "Cro", "CI", "LexA")


# Familias predominantemente represoras (LysR y AraC son mayoritariamente activadoras)
REPRESSOR_FAMILIES = {
    "LacI/GalR", "TetR", "GntR", "DeoR", "MarR", "ArsR/SmtB", "IclR", "RpiR",
    "AsnC", "PadR", "Rrf2", "HxlR", "CopY", "CsoR", "NrdR", "MetJ", "DtxR", "Fur",
}


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- 1. genes (reuso)
def predict_genes(fasta: Path, outdir: Path) -> tuple[list[dict], dict]:
    """Reutiliza `webapp.annotator.annotate` (llamada de genes del repo) en vez de
    duplicar pyrodigal. Devuelve genes normalizados + info del genoma."""
    sys.path.insert(0, str(ROOT / "webapp"))
    from annotator import annotate  # noqa: E402

    outdir.mkdir(parents=True, exist_ok=True)
    data = fasta.read_bytes()
    res = annotate(fasta.stem, data, {"cds": True})
    (outdir / "proteins.faa").write_text(res["proteins_fasta"])
    (outdir / "annotation.gff3").write_text(res["gff"])

    aa = {p["id"]: p["aa"] for p in res["proteins"]}
    genes = []
    for f in res["features"]:
        if f.get("type") != "CDS":
            continue
        genes.append({"contig": f["contig"], "gene": f["id"],
                      "start": int(f["start"]) - 1, "end": int(f["end"]),
                      "strand": f["strand"], "protein": aa.get(f["id"], "")})
    info = {"bp": res["bp"], "gc": res["gc"], "contigs": res["contigs"]}
    return genes, info


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


# ---------------------------------------------------------------- 1b. homología (DIAMOND)
def build_ref_db(faa: Path, dmnd: Path) -> bool:
    """Construye la DB de TFs de referencia con diamond makedb (cacheada)."""
    if dmnd.exists() and dmnd.stat().st_mtime >= faa.stat().st_mtime:
        return True
    r = subprocess.run(["diamond", "makedb", "--in", str(faa), "-d", str(dmnd), "--quiet"],
                       capture_output=True, text=True)
    return r.returncode == 0 and dmnd.exists()


def load_sprot_titles(faa: Path) -> dict[str, dict]:
    out = {}
    if not faa.exists():
        return out
    for line in faa.read_text(errors="ignore").splitlines():
        if not line.startswith(">"):
            continue
        head = line[1:]
        sid = head.split()[0]
        m = re.search(r"OS=([^=]+?)(?: OX=| PE=| SV=| GN=|$)", head)
        desc = re.sub(r"^\S+\s+", "", head)
        out[sid] = {"desc": desc[:120], "os": (m.group(1).strip() if m else "")}
    return out


def run_diamond_db(cand_faa: Path, dmnd: Path, outdir: Path, tag: str) -> dict[str, dict]:
    if not dmnd.exists():
        return {}
    tsv = outdir / f"homologs_{tag}.tsv"
    r = subprocess.run(
        ["diamond", "blastp", "-q", str(cand_faa), "-d", str(dmnd), "-o", str(tsv),
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "evalue",
         "bitscore", "qcovhsp", "--threads", "2", "--quiet", "--max-target-seqs", "3"],
        capture_output=True, text=True)
    if r.returncode != 0 or not tsv.exists():
        return {}
    best = {}
    for line in tsv.read_text(errors="ignore").splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        q, sid = f[0], f[1]
        try:
            rec = {"gene": sid, "pident": float(f[2]), "evalue": float(f[5]),
                   "bitscore": float(f[6]), "qcov": float(f[7])}
        except ValueError:
            continue
        if q not in best or rec["bitscore"] > best[q]["bitscore"]:
            best[q] = rec
    return best


def run_diamond(cand_faa: Path, db_faa: Path, outdir: Path) -> dict[str, dict]:
    """diamond blastp de los candidatos contra la DB de TFs de referencia.

    Devuelve {proteína: {gene, pident, qcov, evalue}} con el mejor hit por bitscore.
    Si diamond no está o falla, devuelve {} y el llamador cae a k-mers.
    """
    # DB compartida entre corridas (no por salida)
    dmnd = ROOT / "data" / "allostery" / "work" / "ref_tfs.dmnd"
    if not build_ref_db(db_faa, dmnd):
        return {}
    tsv = outdir / "homologs.tsv"
    r = subprocess.run(
        ["diamond", "blastp", "-q", str(cand_faa), "-d", str(dmnd), "-o", str(tsv),
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "evalue",
         "bitscore", "qcovhsp", "--threads", "2", "--quiet", "--max-target-seqs", "5"],
        capture_output=True, text=True)
    if r.returncode != 0 or not tsv.exists():
        return {}
    best = {}
    for line in tsv.read_text(errors="ignore").splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        q, s = f[0], f[1]
        try:
            rec = {"gene": s, "pident": float(f[2]), "evalue": float(f[5]),
                   "bitscore": float(f[6]), "qcov": float(f[7])}
        except ValueError:
            continue
        if q not in best or rec["bitscore"] > best[q]["bitscore"]:
            best[q] = rec
    return best


# ---------------------------------------------------------------- 2. dominios TF
def find_tf_domains(faa: Path, outdir: Path) -> dict[str, dict]:
    dom = outdir / "tf.domtbl"
    subprocess.run(["hmmsearch", "--domtblout", str(dom), "--cpu", "2", "-E", "1e-3",
                    "--noali", str(HMM), str(faa)],
                   capture_output=True, text=True, check=True)
    best = {}
    for line in dom.read_text(errors="ignore").splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        prot, hmm = f[0], f[3]
        if not hmm.startswith(DBD_OK):     # el panel trae modelos no-TF
            continue
        score, ev = float(f[13]), float(f[12])
        if ev > 1e-5:
            continue
        if prot not in best or score > best[prot]["score"]:
            best[prot] = {"hmm": hmm, "score": score, "evalue": ev,
                          "from": int(f[19]), "to": int(f[20])}
    return best


# ---------------------------------------------------------------- 3. referencia
def reference_families(prots: dict[str, str]) -> dict[str, str]:
    """Familia de cada TF de referencia por hmmsearch (cacheado)."""
    cache = ROOT / "data" / "allostery" / "work" / "ref_tf_families.json"
    if cache.exists():
        return json.loads(cache.read_text())
    tmp = ROOT / "data" / "allostery" / "work" / "ref_tfs.faa"
    with open(tmp, "w") as fh:
        for k, v in prots.items():
            fh.write(f">{k}\n{v}\n")
    dom = ROOT / "data" / "allostery" / "work" / "ref_tfs.domtbl"
    subprocess.run(["hmmsearch", "--domtblout", str(dom), "--cpu", "2", "-E", "1e-3",
                    "--noali", str(HMM), str(tmp)], capture_output=True, text=True)
    best = {}
    for line in dom.read_text(errors="ignore").splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        if float(f[12]) > 1e-5:
            continue
        if f[0] not in best or float(f[13]) > best[f[0]][0]:
            best[f[0]] = (float(f[13]), f[3])
    out = {k: FAMILY_MAP.get(h, h) for k, (_, h) in best.items()}
    cache.write_text(json.dumps(out, indent=1))
    return out


def load_reference():
    from Bio import SeqIO

    recs = []
    for rec in SeqIO.parse(str(ROOT / "data" / "test_ecoli.gbff"), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            gene = feat.qualifiers.get("gene", [None])[0]
            tr = feat.qualifiers.get("translation", [None])[0]
            if gene and tr:
                recs.append({"gene": gene, "protein": tr})
    eff = {}
    for r in csv.DictReader(open(RD / "TFSet.tsv"), delimiter="\t"):
        tf = r["2)name"]
        ea = r["10)tfActiveConformationsEffector"].strip()
        ei = r["11)tfInactiveConformationsEffector"].strip()
        if not (ea or ei):
            continue
        eff[tf.lower()] = {"modo": ("correpresor" if ea and not ei else
                            "inductor" if ei and not ea else "ambiguo"),
                   "efector": (ea or ei).split("|")[0],
                   "conf": r["17)confidenceLevel"]}
    tfnames = {r["2)name"].lower() for r in csv.DictReader(open(RD / "TFSet.tsv"), delimiter="\t")}
    prots = {r["gene"]: r["protein"] for r in recs}
    fam = reference_families({k: v for k, v in prots.items() if k.lower() in tfnames})
    return recs, eff, fam


def load_multiorg_effectors() -> dict[str, dict]:
    """TF (por nombre) -> efector, desde RegPrecise y minado de UniProt.

    Se indexa por nombre de TF, no por organismo: sirve para transferir a un candidato
    cuyo homólogo se llama igual, aunque venga de otra especie.
    """
    out = {}
    rp = ROOT / "data" / "allostery" / "work" / "regprecise_effectors.tsv"
    if rp.exists():
        for r in csv.DictReader(open(rp), delimiter="\t"):
            name = r["tf"].strip().lower()
            if name and name not in out:
                out[name] = {"efector": r["efector"], "fuente": "regprecise",
                             "organismo": r["organismo"]}
    up = ROOT / "data" / "allostery" / "work" / "uniprot_effectors.tsv"
    if up.exists():
        for r in csv.DictReader(open(up), delimiter="\t"):
            name = (r.get("tf") or "").strip().lower()
            if name and name not in out:
                out[name] = {"efector": r["efector"], "fuente": "uniprot_text",
                             "organismo": r.get("organismo", "")}
    return out


def load_pwms():
    import ast

    out = {}
    for r in csv.DictReader(open(PWM), delimiter="\t"):
        try:
            m = ast.literal_eval(r["pwm"])
        except Exception:
            continue
        if not m or len(m) < CORE_W:
            continue
        info = [max(row) - min(row) for row in m]
        bi = max(range(len(m) - CORE_W + 1), key=lambda i: sum(info[i:i + CORE_W]))
        core = np.array(m[bi:bi + CORE_W])
        cons = "".join("ACGT"[int(np.argmax(row))] for row in core)
        self_s = float(sum(core[i][IDX[cons[i]]] for i in range(CORE_W)))
        if self_s > 0:
            out[r["tf"]] = {"core": core, "self": self_s}
    return out


def best_score(core, self_s, seq: str) -> float:
    if len(seq) < CORE_W:
        return -1e9
    best = -1e9
    for s in range(len(seq) - CORE_W + 1):
        acc = 0.0
        ok = True
        for i, ch in enumerate(seq[s:s + CORE_W]):
            b = IDX.get(ch, -1)
            if b < 0:
                ok = False
                break
            acc += core[i][b]
        if ok:
            best = max(best, acc)
    return best / self_s


def similarity(a: str, b: str, k: int = 4) -> float:
    """Jaccard de k-mers como proxy de identidad (rápido y suficiente para el prototipo).

    Calibración: proteínas idénticas → 1,0; homólogas de la misma familia → 0,2-0,5;
    no relacionadas → <0,05.
    """
    A = {a[i:i + k] for i in range(max(0, len(a) - k))}
    B = {b[i:i + k] for i in range(max(0, len(b) - k))}
    if not A or not B:
        return 0.0
    return len(A & B) / len(A | B)


# ---------------------------------------------------------------- 4. sitios
def upstream(seq: str, gene: dict, length: int = UPSTREAM) -> str:
    if gene["strand"] == "+":
        s = max(0, gene["start"] - length)
        return seq[s:gene["start"]]
    e = min(len(seq), gene["end"] + length)
    return rc(seq[gene["end"]:e])


def scan_sites(region: str, pwms: dict, bg_get, top: int = 3):
    hits = []
    for tf, p in pwms.items():
        s = best_score(p["core"], p["self"], region)
        if s <= -1e8:
            continue
        bgv = bg_get(tf, p)
        pct = float((bgv < s).mean()) if bgv is not None and len(bgv) else None
        hits.append({"pwm_de": tf, "score": round(s, 3),
                     "percentil": None if pct is None else round(pct, 4)})
    hits.sort(key=lambda h: -(h["percentil"] or 0))
    return hits[:top]


def make_bg_get(genome_seqs: dict[str, str], n: int = 600, L: int = 20):
    """Fondo genómico perezoso: solo se calcula para los PWM que se usan."""
    rng = np.random.default_rng(0)
    wins = []
    for seq in genome_seqs.values():
        if len(seq) <= L:
            continue
        starts = rng.integers(0, len(seq) - L, max(1, n // max(1, len(genome_seqs))))
        wins += [seq[s:s + L] for s in starts]
    cache = {}

    def get(tf: str, p: dict):
        if tf not in cache:
            cache[tf] = np.array([best_score(p["core"], p["self"], w) for w in wins])
        return cache[tf]

    return get


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("--out", default=str(ROOT / "data" / "allostery" / "genome_run"))
    ap.add_argument("--max-tfs", type=int, default=40)
    ap.add_argument("--find-repressors", action="store_true",
                    help="mostrar solo familias predominantemente represoras")
    ap.add_argument("--tier", choices=["rapida", "estandar"], default="estandar",
                    help="rapida: sin escaneo de sitios (segundos) · "
                         "estandar: incluye sitios PWM (minutos)")
    args = ap.parse_args()

    out = Path(args.out)
    fasta = Path(args.fasta)
    data = fasta.read_bytes()
    key = hashlib.sha256(
        data + f"{args.max_tfs}|{args.find_repressors}|{args.tier}".encode()).hexdigest()[:16]
    cache = ROOT / "data" / "allostery" / "work" / "cache" / key
    if cache.exists() and (cache / "summary.json").exists():
        out.mkdir(parents=True, exist_ok=True)
        for f in cache.iterdir():
            if f.is_file():
                (out / f.name).write_bytes(f.read_bytes())
        print(f"cache hit {key}: {out}", file=sys.stderr)
        print((cache / "summary.json").read_text())
        return
    out.mkdir(parents=True, exist_ok=True)

    timings = {}
    t0 = time.perf_counter()
    print("1/5 genes+anotación (reuso webapp.annotator)…", file=sys.stderr)
    genes, ginfo = predict_genes(fasta, out)
    seqs = read_fasta(fasta)
    timings["genes_s"] = round(time.perf_counter() - t0, 1)
    print(f"    {len(genes)} genes · {timings['genes_s']} s", file=sys.stderr)

    t1 = time.perf_counter()
    print("2/5 dominios de TF (hmmsearch)…", file=sys.stderr)
    doms = find_tf_domains(out / "proteins.faa", out)
    timings["dominios_s"] = round(time.perf_counter() - t1, 1)
    print(f"    {len(doms)} candidatos · {timings['dominios_s']} s", file=sys.stderr)

    t2 = time.perf_counter()
    print("3/5 referencia, DB de TFs y PWM…", file=sys.stderr)
    refs, eff, fam_of_ref = load_reference()
    multi = load_multiorg_effectors()
    pwms = load_pwms()
    print(f"    efectores multi-organismo: {len(multi)} TFs", file=sys.stderr)
    tf_names = {r["2)name"].lower() for r in csv.DictReader(open(RD / "TFSet.tsv"), delimiter="\t")}
    ref_tfs = {r["gene"]: r["protein"] for r in refs if r["gene"].lower() in tf_names}
    ref_faa = ROOT / "data" / "allostery" / "work" / "ref_tfs.faa"
    if not ref_faa.exists() or ref_faa.stat().st_size < 1000:
        with open(ref_faa, "w") as fh:
            for k, v in ref_tfs.items():
                fh.write(f">{k}\n{v}\n")
    hits = run_diamond(out / "proteins.faa", ref_faa, out)

    # Unión de evidencias: un candidato puede venir del dominio HMM o de un homólogo
    # claro de un TF de referencia (el panel HMM no cubre todas las familias DBD).
    cand = set(doms)
    for prot, h in hits.items():
        if h["pident"] >= 40.0 and h["qcov"] >= 70.0 and h["evalue"] <= 1e-5:
            cand.add(prot)
    print(f"    candidatos HMM={len(doms)} ∪ homólogos={len(cand) - len(doms)} "
          f"= {len(cand)}", file=sys.stderr)
    sprot_dmnd = ROOT / "data" / "allostery" / "work" / "sprot_tf.dmnd"
    sprot_hits = run_diamond_db(out / "proteins.faa", sprot_dmnd, out, "sprot") if sprot_dmnd.exists() else {}
    sprot_titles = load_sprot_titles(ROOT / "data" / "allostery" / "work" / "sprot_tf.faa")
    homolog_metodo = "diamond" if hits else "kmer"
    print(f"    Swiss-Prot TFs: {len(sprot_hits)} hits", file=sys.stderr)
    timings["referencia_s"] = round(time.perf_counter() - t2, 1)
    print(f"    referencia {len(ref_tfs)} TFs · homología: {homolog_metodo}", file=sys.stderr)

    t3 = time.perf_counter()
    bg_get = make_bg_get(seqs)
    timings["fondo_s"] = 0.0  # perezoso: se contabiliza dentro de la asignación

    t4 = time.perf_counter()
    print("5/5 homólogo, familia, efector y sitios…", file=sys.stderr)
    by_name = {g["gene"]: g for g in genes}
    rows = []
    order = sorted(cand, key=lambda p: -(doms.get(p, {}).get("score", 0)))
    # con filtro de represores no se recorta antes de filtrar (sesgaría a los mejor
    # puntuados, que suelen ser activadores de familias grandes como LysR)
    limite = len(order) if args.find_repressors else args.max_tfs
    for prot in order[:limite]:
        organismo = descripcion = ""
        g = by_name.get(prot)
        if not g:
            continue
        d = doms.get(prot) or {"hmm": "sin_dominio_HMM", "score": 0.0, "evalue": None}
        family = FAMILY_MAP.get(d["hmm"], d["hmm"])
        if d["hmm"] == "sin_dominio_HMM":
            fam_ref = fam_code_to_family(fam_of_ref.get((hits.get(prot) or {}).get("gene", "").lower(), ""))
            family = fam_ref or "sin_dominio_HMM"
        fam_refs = [r for r in refs if fam_of_ref.get(r["gene"]) and
                    fam_code_to_family(fam_of_ref[r["gene"]]) == family]
        pool = fam_refs or refs[:80]
        best = max(pool, key=lambda r: similarity(g["protein"], r["protein"]))
        ident = similarity(g["protein"], best["protein"])
        dh = hits.get(prot)
        if dh and dh["pident"] >= 40.0 and dh["qcov"] >= 70.0 and dh["evalue"] <= 1e-5:
            best = {"gene": dh["gene"], "protein": ref_tfs.get(dh["gene"], "")}
            ident = dh["pident"] / 100.0
            homolog_metodo_gen = "diamond"
        elif sprot_hits.get(prot) and sprot_hits[prot]["pident"] >= 40.0 \
                and sprot_hits[prot]["qcov"] >= 70.0 and sprot_hits[prot]["evalue"] <= 1e-5:
            sh = sprot_hits[prot]
            info_s = sprot_titles.get(sh["gene"], {})
            best = {"gene": sh["gene"], "protein": ""}
            ident = sh["pident"] / 100.0
            homolog_metodo_gen = "diamond_sprot"
            organismo = info_s.get("os", "")
            descripcion = info_s.get("desc", "")
        else:
            homolog_metodo_gen = "kmer"
            if ident < 0.25:               # por debajo, no es un homólogo: no se reporta
                best = {"gene": "—", "protein": ""}
                ident = 0.0
                homolog_metodo_gen = "sin_homologo"
        e = eff.get(best["gene"].lower())
        mnem = ""
        if best["gene"].startswith("sp|"):
            mnem = best["gene"].split("|")[-1].split("_")[0].lower()
        m = multi.get(mnem) or multi.get(best["gene"].lower())
        if e and ident >= (0.35 if homolog_metodo_gen == "kmer" else 0.40):
            modo, efector, status = e["modo"], e["efector"], "inherited_homolog"
        elif m:
            modo, efector, status = "desconocido", m["efector"], m["fuente"] + "_tf_name"
        else:
            cnt = collections.Counter(v["efector"] for k, v in eff.items()
                                      if fam_code_to_family(fam_of_ref.get(k, "")) == family)
            modo = "desconocido"
            efector = cnt.most_common(1)[0][0] if cnt else None
            status = "family_prior" if efector else "unknown"
        sites = []
        if args.tier == "estandar":
            region = upstream(seqs.get(g["contig"], ""), g)
        # PWM restringidos a la familia del candidato; si no hay, solo el del homólogo
        pwm_lc = {k.lower(): (k, v) for k, v in pwms.items()}
        fam_tfs = [r["gene"] for r in refs
                   if fam_code_to_family(fam_of_ref.get(r["gene"].lower(), "")) == family]
        sub_pwms = {}
        for t in fam_tfs:
            hit = pwm_lc.get(t.lower())
            if hit:
                sub_pwms[hit[0]] = hit[1]
        if not sub_pwms:
            hit = pwm_lc.get(best["gene"].lower())
            if hit:
                sub_pwms[hit[0]] = hit[1]
        if len(sub_pwms) > 8:                      # acota coste: 8 PWM por familia
            sub_pwms = dict(list(sub_pwms.items())[:8])
            sites = scan_sites(region, sub_pwms or pwms, bg_get)
        rows.append({"gen": g["gene"], "contig": g["contig"], "inicio": g["start"],
                     "fin": g["end"], "hebra": g["strand"], "aa": len(g["protein"]),
                     "dominio": d["hmm"], "familia": family, "hmm_evalue": d["evalue"],
                     "homologo": best["gene"], "identidad": round(ident, 3),
                     "homolog_metodo": homolog_metodo_gen,
                     "organismo": organismo, "descripcion": descripcion,
                     "modo": modo, "efector": efector, "status": status,
                     "sitios": sites})

    if args.find_repressors:
        before = len(rows)
        rows = [r for r in rows if r["familia"] in REPRESSOR_FAMILIES or
                r["modo"] == "correpresor"]
        print(f"    filtro represores: {before} -> {len(rows)}", file=sys.stderr)

    with open(out / "represores.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "gen", "contig", "inicio", "fin", "hebra", "aa", "dominio", "familia",
            "hmm_evalue", "homologo", "identidad", "homolog_metodo", "organismo",
            "descripcion", "modo", "efector", "status", "sitios"])
        w.writeheader()
        for r in rows:
            w.writerow({**r, "sitios": json.dumps(r["sitios"], ensure_ascii=False)})

    timings["asignacion_s"] = round(time.perf_counter() - t4, 1)
    timings["total_s"] = round(time.perf_counter() - t0, 1)
    peak_mb = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)

    summary = {"genoma": fasta.name, "nivel": args.tier,
               "profunda_bajo_demanda": ["conservacion_entre_genomas (E2)",
                                         "cofolding_interfaz (E4, GPU)",
                                         "MM/GBSA y ABFE del efector (GPU, pago)"],
               "bp": ginfo["bp"], "gc": ginfo["gc"],
               "genes": len(genes), "tf_candidatos": len(doms),
               "homolog_metodo": homolog_metodo,
               "recursos": {"tiempos_s": timings, "pico_ram_mb": peak_mb,
                            "cpus": __import__("os").cpu_count()},
               "solo_represores": bool(args.find_repressors),
               "reportados": len(rows),
               "con_efector_heredado": sum(1 for r in rows if r["status"] == "inherited_homolog"),
               "con_prior_de_familia": sum(1 for r in rows if r["status"] == "family_prior"),
               "sin_efector": sum(1 for r in rows if r["status"] == "unknown"),
               "familias": dict(collections.Counter(r["familia"] for r in rows))}
    (out / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    try:
        cache.mkdir(parents=True, exist_ok=True)
        for f in out.iterdir():
            if f.is_file():
                (cache / f.name).write_bytes(f.read_bytes())
    except Exception:
        pass

    L = [f"# Informe de reguladores — {fasta.name}", "",
         f"Genes: {len(genes)} · candidatos TF: {len(doms)} · reportados: {len(rows)}",
         "", "| gen | familia | homólogo | ident. | modo | efector | status |", "|---|---|---|---:|---|---|---|"]
    for r in rows:
        L.append(f"| {r['gen']} | {r['familia']} | {r['homologo']} | {r['identidad']:.2f} | "
                 f"{r['modo']} | {r['efector'] or '—'} | {r['status']} |")
    (out / "informe.md").write_text("\n".join(L) + "\n")
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
