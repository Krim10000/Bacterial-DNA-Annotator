#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_fold_candidates.py — Criba de enzimas faltantes por plegado (ESMFold + Foldseek).

v1.3 (2026-11): hipotéticas SIN anotar + Foldseek SwissProt (AF2).

Embudo (orden de aplicación):
  1) TRIAGE del 04 (ON por defecto): solo pasan gaps con
     triage == 'genuine_missing_enzyme'.
  2) COHORTE del 05 (--use-cohort, OFF): excluye gaps en
     'systematic_gaps' (presentes en TODOS los genomas = artefacto).
     --cohort-specific-only reproduce el filtro antiguo (solo specific).
  3) CANDIDATAS, dos fuentes:
     a) Hipotéticas YA MAPPED a producing_rhea de gaps genuinos (pequeño,
        valida anotación existente).
     b) --include-unannotated-hyp: hypothetical protein del .faa que NO
        están en el gene map (el caso real de "enzima faltante"):
        proteínas sin EC asignado.
     c) --orfs: ORFs de contigs sin gen anotado.
  4) SIMILITUD PRIMERO (--strategy similarity_first): no pliega si hay hit
     mmseqs/blastp >=30% id y >=70% cob.
  5) PLEGADO con ESMFold (solo si el paso 4 no resuelve).
  6) FOLDSEEK tras plegar (--no-foldseek lo desactiva): easy-search del PDB
     contra data/db/foldseek_swissprot (AF2 SwissProt). Hit con TM > 0.5
     contra una proteína con Rhea relevante -> veredicto fuerte.

Veredicto:
  fuerte      seq >=30/70  O  dominio estructurado + hit  O  Foldseek TM>0.5
              contra Rhea relevante  O  (transitively) hit Foldseek con Rhea
              en acc_to_rheas ∩ action_rheas
  validar     dominio estructurado sin evidencia de similitud
  descartado  pLDDT<50 o sin dominio
  sin plegado API falló

Nuevos flags v1.3:
  --include-unannotated-hyp     pliega hypotheticals sin RHEA
  --hyp-no-coords               (aviso) unannotated no tienen coordenadas
  --foldseek-db PATH            default data/db/foldseek_swissprot
  --no-foldseek                 desactiva Foldseek
  --foldseek-top N              hits a conservar por candidata (default 3)
  --foldseek-tm FLOAT           TM-score mínimo (default 0.5)
  --cohort-specific-only        filtro cohort antiguo
  --max-unannotated-hyp N       tope de hipotéticas sin anotar (default 200)
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
m04 = importlib.import_module("04_network_trace")
m07 = importlib.import_module("07_network_view")

log = logging.getLogger("fold_candidates")
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-7s | %(message)s",
                    datefmt="%H:%M:%S")

SCRIPT_VERSION = "1.3"
ESM_URL_DEFAULT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
NEUROSNAP_BASE_DEFAULT = "https://neurosnap.ai/api"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
REDUCED_GENOME_BP = 1_000_000
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_ALT = {"B": "X", "Z": "X", "J": "X", "U": "C", "O": "K"}
IUPAC_COMP = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")

CODONS = {}
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for b1 in _BASES:
    for b2 in _BASES:
        for b3 in _BASES:
            CODONS[b1 + b2 + b3] = _AAS[_i]
            _i += 1
STARTS = {"ATG", "GTG", "TTG"}


# ============================================================ resultado vacío
def _empty_result(genome, gaps, total_bp, reduced, context_note, funnel=None,
                  diag="sin candidatas"):
    return {
        "genome": genome,
        "candidates": [],
        "candidates_folded": 0,
        "candidates_with_pdb": 0,
        "candidates_skipped_fold": 0,
        "gaps": len(gaps) if gaps else 0,
        "gap_refs": 0,
        "context": context_note,
        "reduced_genome": reduced,
        "dna": {"genome_size_bp": total_bp, "reduced_genome": reduced},
        "funnel": funnel or {},
        "diag": diag,
    }


# ============================================================ filtros gaps
def _apply_gap_filters(gaps, use_triage, cohort_cfg, genome):
    """Devuelve (gaps_filtrados, funnel).

    cohort_cfg: None | {'mode': 'exclude_systematic'|'specific_only',
                        'set': set(...)}
    """
    funnel = {
        "gaps_total": len(gaps),
        "by_triage": {},
        "after_triage": len(gaps),
        "after_cohort": len(gaps),
        "cohort_mode": None,
        "candidate_rheas": 0,
    }
    out = []
    warned_old = False
    if use_triage:
        by_triage = defaultdict(int)
        for g in gaps:
            tri = g.get("triage")
            if tri is None:
                if not warned_old:
                    log.warning(f"[{genome}] trace 04 sin 'triage' — se "
                                "mantienen todos los gaps.")
                    warned_old = True
                out.append(g)
            else:
                by_triage[tri] += 1
                if tri == "genuine_missing_enzyme":
                    out.append(g)
        funnel["by_triage"] = dict(by_triage)
        funnel["after_triage"] = len(out)
    else:
        out = list(gaps)

    if cohort_cfg:
        s = cohort_cfg["set"]
        mode = cohort_cfg["mode"]
        funnel["cohort_mode"] = mode
        if mode == "exclude_systematic":
            out = [g for g in out if g.get("chebi") not in s]
        elif mode == "specific_only":
            out = [g for g in out if g.get("chebi") in s]
        funnel["after_cohort"] = len(out)

    cand_rheas = set()
    for g in out:
        for r in g.get("producing_rhea", []):
            norm = m04.normalize_rhea_id(r)
            if norm:
                cand_rheas.add(norm[1])
    funnel["candidate_rheas"] = len(cand_rheas)
    return out, funnel


# ============================================================ referencias
def rhea_to_uniprot(rhea_dir: Path):
    up = rhea_dir / "rhea2uniprot.tsv"
    if not up.exists():
        up = rhea_dir / "rhea2uniprot.tsv.gz"
    if not up.exists():
        log.warning(f"Sin {rhea_dir}/rhea2uniprot.tsv: sin referencias por "
                    "brecha. Descárgalo de https://ftp.expasy.org/databases/"
                    "rhea/tsv/rhea2uniprot.tsv")
        return {}
    opener = gzip.open if up.suffix == ".gz" else open
    rhea2acc, master2acc = {}, {}
    header = None
    with opener(up, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip().upper() for h in parts]
                continue
            if len(parts) < 4:
                continue
            d = dict(zip(header, parts))
            acc = (d.get("ID") or d.get("UNIPROT") or "").strip()
            rid_raw = (d.get("RHEA_ID_DIR") or d.get("RHEA_ID_LR")
                       or d.get("RHEA_ID") or "")
            mid_raw = (d.get("RHEA_ID_MASTER") or d.get("MASTER_ID") or "")
            rid = m04.normalize_rhea_id(rid_raw)
            mid = m04.normalize_rhea_id(mid_raw)
            if not acc or acc.upper() == "ID" or len(acc) < 4:
                continue
            if rid:
                rhea2acc.setdefault(rid[1], set()).add(acc)
            if mid:
                master2acc.setdefault(mid[1], set()).add(acc)
    for rid, accs in master2acc.items():
        rhea2acc.setdefault(rid, set()).update(accs)
    return {k: sorted(v) for k, v in rhea2acc.items()}


def fetch_uniprot_fasta(acc, cache_dir: Path, timeout=60) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = cache_dir / f"{acc}.fasta"
    if p.exists():
        return read_fasta_first(p)
    try:
        with urllib.request.urlopen(UNIPROT_FASTA.format(acc=acc),
                                    timeout=timeout) as r:
            text = r.read().decode("utf-8", errors="replace")
        p.write_text(text, encoding="utf-8")
        time.sleep(0.3)
        return read_fasta_first(p)
    except Exception as e:
        log.warning(f"UniProt {acc}: {e}")
        return ""


def open_text(path):
    raw = open(path, "rb")
    if hasattr(raw, "peek"):
        gz = raw.peek(2)[:2] == b"\x1f\x8b"
    else:
        gz = raw.read(2) == b"\x1f\x8b"
        raw.seek(0)
    if gz:
        raw.close()
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")


def read_fasta_first(path: Path) -> str:
    seq = []
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def read_fasta_dict(path: Path):
    out, cur_id, desc, seq = {}, None, "", []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id:
                    out[cur_id] = (desc, "".join(seq))
                head = line[1:].split()
                cur_id = head[0] if head else ""
                desc = line[1:]
                seq = []
            elif line:
                seq.append(line.strip())
    if cur_id:
        out[cur_id] = (desc, "".join(seq))
    return out


# ============================================================ ORFs
def find_orfs(dna: str, min_aa: int, max_orfs: int):
    dna = "".join(dna.split()).upper()
    rc = dna.translate(IUPAC_COMP)[::-1]
    orfs = []
    for strand, s in (("+", dna), ("-", rc)):
        for frame in range(3):
            for i in range(frame, len(s) - 2, 3):
                if s[i:i + 3] not in STARTS:
                    continue
                j = i + 3
                prot = ["M"]
                stop = False
                while j + 3 <= len(s):
                    aa = CODONS.get(s[j:j + 3], "X")
                    if aa == "*":
                        stop = True
                        break
                    prot.append(aa)
                    j += 3
                if not stop or len(prot) < min_aa:
                    continue
                nt_len = (j + 3) - i
                if strand == "+":
                    gstart, gend = i + 1, i + nt_len
                else:
                    gstart = len(s) - (i + nt_len) + 1
                    gend = len(s) - i
                orfs.append({
                    "id": "", "seq": "".join(prot),
                    "start": gstart, "end": gend,
                    "strand": strand, "length_bp": nt_len,
                })
    orfs.sort(key=lambda o: -o["length_bp"])
    return orfs[:max_orfs]


def sanitize_protein(seq: str):
    seq = "".join(seq.split()).upper().replace("*", "").replace("-", "")
    out, subs = [], {}
    for ch in seq:
        if ch in AA:
            out.append(ch)
        elif ch in AA_ALT:
            out.append(AA_ALT[ch])
            subs[ch] = AA_ALT[ch]
        else:
            out.append("X")
            subs[ch] = "X"
    n_x = out.count("X")
    frac_x = n_x / len(out) if out else 1.0
    note = ""
    if subs:
        note = "residuos no estándar: " + ", ".join(
            f"{k}-> {v} x{seq.count(k)}" for k, v in sorted(subs.items()))
    return "".join(out), note, frac_x


# ============================================================ ESMFold API
def pick_provider(explicit: str) -> str:
    if explicit and explicit != "auto":
        return explicit
    env = os.environ.get("ESM_API_PROVIDER", "").strip().lower()
    if env in ("neurosnap", "direct"):
        return env
    if os.environ.get("ESM_API_KEY", "").strip():
        return "neurosnap"
    return "direct"


def _http(req: urllib.request.Request, timeout: int) -> bytes:
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _multipart_field(field: str, value: str):
    boundary = "----dna_poc_" + hashlib.md5(value.encode()).hexdigest()[:16]
    body = (f"--{boundary}\r\n"
            f"Content-Disposition: form-data; name=\"{field}\"\r\n\r\n"
            f"{value}\r\n--{boundary}--\r\n").encode()
    return body, f"multipart/form-data; boundary={boundary}"


def esm_fold_neurosnap(seq, seq_id, genome, base, api_key, timeout, poll,
                       job_timeout, retries):
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", seq_id)[:40] or "seq"
    note = urllib.parse.quote(f"{genome}-{seq_id}"[:60], safe="")
    payload = json.dumps({"aa": {name: seq}})
    hdr = {"X-API-KEY": api_key, "User-Agent": "dna_poc/08_fold_candidates"}
    last = None
    for attempt in range(1, retries + 1):
        try:
            body, ctype = _multipart_field("Input Sequences", payload)
            req = urllib.request.Request(
                f"{base}/job/submit/ESMFold?note={note}", data=body,
                method="POST", headers={**hdr, "Content-Type": ctype})
            resp = _http(req, timeout).decode("utf-8", "replace").strip()
            try:
                job = json.loads(resp)
            except json.JSONDecodeError:
                job = resp.strip('" ')
            if isinstance(job, dict):
                job = job.get("job_id") or job.get("id") or ""
            if not job:
                raise ValueError(f"submit sin job id: {resp[:120]}")
            log.info(f"    Neurosnap job {job} ({len(seq)} aa, intento {attempt})")
            deadline = time.time() + job_timeout
            while True:
                st = _http(urllib.request.Request(
                    f"{base}/job/status/{job}", headers=hdr),
                    timeout).decode("utf-8", "replace").strip()
                try:
                    st = json.loads(st)
                except json.JSONDecodeError:
                    pass
                if isinstance(st, dict):
                    st = st.get("status", "")
                st = str(st).strip().lower()
                if st == "completed":
                    break
                if st in ("failed", "cancelled", "deleted", "error"):
                    raise ValueError(f"job terminó con estado '{st}'")
                if time.time() > deadline:
                    raise TimeoutError(f"job {job} sin terminar en "
                                       f"{job_timeout}s (estado: {st or 'n/d'})")
                time.sleep(poll)
            data = json.loads(_http(urllib.request.Request(
                f"{base}/job/data/{job}", headers=hdr),
                timeout).decode("utf-8", "replace"))
            pdb_name = None
            for fi in data.get("out", []):
                fname = (fi[0] if isinstance(fi, (list, tuple))
                         else fi.get("name", "") if isinstance(fi, dict)
                         else str(fi))
                if str(fname).endswith(".pdb"):
                    pdb_name = str(fname)
                    break
            if not pdb_name:
                raise ValueError(f"sin .pdb en la salida del job {job}")
            pdb = _http(urllib.request.Request(
                f"{base}/job/file/{job}/out/{urllib.parse.quote(pdb_name)}",
                headers=hdr), timeout).decode("utf-8", "replace")
            if "ATOM" not in pdb:
                raise ValueError("PDB descargado sin átomos")
            return pdb
        except Exception as e:
            last = e
            log.warning(f"Neurosnap intento {attempt}/{retries}: {e}")
            time.sleep(min(5 * attempt, 15))
    log.error(f"Neurosnap falló para {seq_id} ({len(seq)} aa): {last}")
    return ""


def esm_fold_direct(seq, url, api_key, timeout, retries):
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url, data=seq.encode(), method="POST",
                headers={"Content-Type": "text/plain",
                         "User-Agent": "dna_poc/08_fold_candidates"})
            if api_key:
                req.add_header("X-API-Key", api_key)
                req.add_header("Authorization", f"Bearer {api_key}")
            pdb = _http(req, timeout).decode("utf-8", "replace")
            if "ATOM" not in pdb:
                raise ValueError("respuesta sin átomos")
            return pdb
        except Exception as e:
            last_err = e
            log.warning(f"ESMFold intento {attempt}/{retries}: {e}")
            time.sleep(2 * attempt)
    log.error(f"ESMFold falló para seq de {len(seq)} aa: {last_err}")
    return ""


def esm_fold(seq, cache_dir, timeout, retries, pause, provider="auto",
             seq_id="seq", genome="", poll=5.0, job_timeout=900):
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(seq.encode()).hexdigest()
    p = cache_dir / f"{key}.pdb"
    if p.exists():
        return p.read_text(encoding="utf-8", errors="replace"), p
    api_key = os.environ.get("ESM_API_KEY", "").strip()
    prov = pick_provider(provider)
    if prov == "neurosnap":
        if not api_key:
            log.error("proveedor neurosnap sin ESM_API_KEY")
            return "", None
        base = os.environ.get("ESM_API_URL", NEUROSNAP_BASE_DEFAULT).rstrip("/")
        pdb = esm_fold_neurosnap(seq, seq_id, genome, base, api_key, timeout,
                                 poll, job_timeout, retries)
    else:
        url = os.environ.get("ESM_API_URL", ESM_URL_DEFAULT)
        pdb = esm_fold_direct(seq, url, api_key, timeout, retries)
    if pdb:
        p.write_text(pdb, encoding="utf-8")
        time.sleep(pause)
        return pdb, p
    return "", None


def plddt_from_pdb(pdb_text: str):
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:76]) * 100.0)
            except ValueError:
                continue
    if not vals:
        return None, None, None
    good = [v >= 70 for v in vals]
    best = run = 0
    for g in good:
        run = run + 1 if g else 0
        best = max(best, run)
    return (sum(vals) / len(vals),
            sum(1 for g in good if g) / len(good),
            best)


# ============================================================ Foldseek
def _parse_foldseek_acc(target: str) -> str:
    """'AF-P12345-F1-model_v4' -> 'P12345'; 'P12345_...' -> 'P12345'."""
    t = target.strip()
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", t)
    if m:
        return m.group(1)
    for sep in ("_", "|", "."):
        if sep in t:
            parts = t.split(sep)
            for p in parts:
                if 4 <= len(p) <= 10 and p[0].isalpha():
                    return p
    return t


def foldseek_search(pdb_path: Path, db_prefix: Path, out_dir: Path,
                    top: int, timeout: int):
    """Ejecuta foldseek easy-search y devuelve lista de hits:
    [{target, acc, qtmscore, ttmscore, alntmscore, seqid, prob}]."""
    exe = shutil.which("foldseek") or str(Path.cwd() / "foldseek" / "bin" / "foldseek")
    if not Path(exe).exists():
        log.warning(f"foldseek no encontrado en PATH ni en {exe}")
        return []
    if not Path(str(db_prefix) + ".dbtype").exists():
        log.warning(f"DB Foldseek no encontrada: {db_prefix}.dbtype")
        return []
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = hashlib.md5(str(pdb_path).encode()).hexdigest()[:10]
    aln = out_dir / f"{tag}.m8"
    tmp = out_dir / f"{tag}_tmp"
    tmp.mkdir(exist_ok=True)
    try:
        subprocess.run(
            [exe, "easy-search", str(pdb_path), str(db_prefix), str(aln),
             str(tmp),
             "--format-output",
             "query,target,qtmscore,ttmscore,alntmscore,seqid,prob",
             "--max-seqs", str(max(top * 2, 5)),
             "-v", "1"],
            check=True, capture_output=True, timeout=timeout)
    except subprocess.CalledProcessError as e:
        log.warning(f"foldseek falló: {e.stderr.decode()[:200] if e.stderr else e}")
        return []
    except subprocess.TimeoutExpired:
        log.warning(f"foldseek timeout ({timeout}s)")
        return []
    hits = []
    if aln.exists():
        for line in aln.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                hits.append({
                    "target": parts[1],
                    "acc": _parse_foldseek_acc(parts[1]),
                    "qtmscore": round(float(parts[2]), 3),
                    "ttmscore": round(float(parts[3]), 3),
                    "alntmscore": round(float(parts[4]), 3),
                    "seqid": round(float(parts[5]), 1),
                    "prob": round(float(parts[6]), 3),
                })
            except ValueError:
                continue
    hits.sort(key=lambda h: -h["qtmscore"])
    return hits[:top]


# ============================================================ similitud
def similarity_search(query_faa, ref_faa, tool, workdir):
    if tool == "auto":
        tool = ("mmseqs" if shutil.which("mmseqs") else
                ("blastp" if shutil.which("blastp") else "none"))
    if tool == "none":
        log.info("Sin mmseqs2/blastp: veredicto por estructura sola.")
        return []
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        if tool == "mmseqs":
            aln = workdir / "res.tab"
            subprocess.run(
                ["mmseqs", "easy-search", str(query_faa), str(ref_faa),
                 str(aln), str(workdir / "tmp"),
                 "--format-output", "query,target,pident,alnlen,qlen,tlen",
                 "-s", "7.5", "--threads", "4"],
                check=True, capture_output=True)
            out = []
            for line in aln.read_text().splitlines():
                q, t, pid, alen, qlen, tlen = line.split("\t")[:6]
                cov = int(alen) / max(1, int(qlen))
                out.append((q, t, float(pid), cov))
            return out
        if tool == "blastp":
            dbdir = workdir / "db"
            dbdir.mkdir(exist_ok=True)
            subprocess.run(["makeblastdb", "-dbtype", "prot", "-in",
                            str(ref_faa), "-out", str(dbdir / "ref")],
                           check=True, capture_output=True)
            res = workdir / "blast.tab"
            subprocess.run(["blastp", "-query", str(query_faa), "-db",
                            str(dbdir / "ref"), "-outfmt",
                            "6 qseqid sseqid pident length qlen",
                            "-evalue", "1e-3", "-out", str(res)],
                           check=True, capture_output=True)
            out = []
            for line in res.read_text().splitlines():
                q, t, pid, alen, qlen = line.split("\t")[:5]
                cov = int(alen) / max(1, int(qlen))
                out.append((q, t, float(pid), cov))
            return out
    except Exception as e:
        log.warning(f"Similitud ({tool}) falló: {e}")
        return []
    return []


def find_seq_file(d, genome, ext):
    if not d:
        return None
    d = Path(d)
    for name in (f"{genome}{ext}", f"{genome}{ext}.gz"):
        p = d / name
        if p.exists():
            return p
    for pat in (f"{genome}*{ext}.gz", f"{genome}*{ext}"):
        hits = sorted(d.glob(pat))
        if hits:
            return hits[0]
    return None


# ============================================================ unannotated hyp
def find_unannotated_hypotheticals(faa_path: Path, gene_map_ids: set,
                                   max_hyps: int, min_aa: int):
    """Hypothetical proteins del .faa que NO están en el gene map.
    Devuelve lista [{id, seq, desc}] ordenada por longitud desc."""
    out = []
    for pid, (desc, seq) in read_fasta_dict(faa_path).items():
        # NCBI protein FASTA: id suele ser WP_/YP_/NP_ + descripción textual.
        # Algunos gene maps usan el locus_tag como gene_id.
        base_id = pid.split("|")[-1] if "|" in pid else pid
        if pid in gene_map_ids or base_id in gene_map_ids:
            continue
        if "hypothetical" not in desc.lower():
            continue
        if len(seq) < min_aa:
            continue
        out.append({"id": pid, "seq": seq, "desc": desc,
                    "length_bp": len(seq) * 3})
    out.sort(key=lambda x: -len(x["seq"]))
    return out[:max_hyps]


# ============================================================ núcleo
def process_genome(genome, net_dir, rhea_dir, gene_map_path, proteins_faa,
                   genome_meta_path, contigs_fna, out_dir, up_cache_dir,
                   esm_cache_dir, work_root, args, rhea2acc, acc_to_rheas,
                   cohort_cfg=None):
    trace_p = net_dir / f"{genome}_network_trace.json"
    if not trace_p.exists():
        log.warning(f"[{genome}] sin {trace_p.name}: corre 04 primero")
        return None
    trace = json.loads(trace_p.read_text(encoding="utf-8"))
    all_gaps = trace.get("gaps_missing_enzyme", [])

    genes = m07.load_gene_map(gene_map_path) if gene_map_path else {}
    meta = m07.load_genome_meta(genome_meta_path, genome)
    total_bp = (meta or {}).get("total_length_bp")
    reduced = bool(total_bp and total_bp < REDUCED_GENOME_BP)
    context_note = (
        f"genoma reducido ({total_bp} bp)" if reduced else
        (f"tamaño del genoma: {total_bp} bp" if total_bp
         else "tamaño del genoma desconocido (falta --genome-meta)"))

    gaps, gap_funnel = _apply_gap_filters(
        all_gaps, not args.no_use_triage, cohort_cfg, genome)
    funnel = dict(gap_funnel)

    if not gaps:
        msg = (f"[{genome}] 0 gaps accionables tras filtros "
               f"(triage: {gap_funnel.get('by_triage')}, "
               f"cohort_mode: {gap_funnel.get('cohort_mode')}). Nada que plegar.")
        log.info(msg)
        return _empty_result(genome, gaps, total_bp, reduced,
                             context_note, funnel=funnel,
                             diag="filtros eliminaron todos los gaps")

    # secuencias disponibles
    seqs = {}
    if proteins_faa and Path(proteins_faa).exists():
        seqs = read_fasta_dict(Path(proteins_faa))
    elif any(row.get("sequence") for rows in genes.values() for row in rows):
        for rows in genes.values():
            for row in rows:
                if row.get("sequence"):
                    seqs[row["gene_id"]] = (row["gene_id"], row["sequence"])

    candidates = []
    seen_ids = set()

    # A) hipotéticas ya asociadas a RHEA de algún gap genuino
    action_rheas = set()
    for g in gaps:
        for r in g.get("producing_rhea", []):
            norm = m04.normalize_rhea_id(r)
            if norm:
                action_rheas.add(norm[1])
    hyp_all = {rid for rid, rows in genes.items()
               if rows and all(r["hypothetical"] for r in rows)}
    hyp_rxn_ids = hyp_all & action_rheas

    funnel.update({
        "gene_map_empty": not genes,
        "hyp_rxn_total": len(hyp_all),
        "hyp_rxn_relevant": len(hyp_rxn_ids),
        "hyp_annotated": 0,
        "hyp_unannotated": 0,
        "orfs_added": 0,
        "total_candidates": 0,
        "strong_seq_skipped": 0,
        "folded": 0,
        "foldseek_hits": 0,
    })

    def add_candidate(cid, seq, kind, contig, start, end, strand, rxn_ids,
                      pool=None):
        pool = candidates if pool is None else pool
        if cid in seen_ids:
            return
        seen_ids.add(cid)
        seq, seq_note, frac_x = sanitize_protein(seq or "")
        if not seq:
            return
        if frac_x > 0.05:
            log.warning(f"[{genome}] {cid}: {frac_x:.0%} X — no se pliega")
            return
        pool.append({
            "id": cid, "kind": kind, "seq": seq, "seq_note": seq_note,
            "contig": contig, "start": start, "end": end,
            "strand": strand, "length_bp": (end - start + 1) if end else 0,
            "rxn_ids": sorted(rxn_ids),
        })

    for rid in sorted(hyp_rxn_ids):
        for row in genes[rid]:
            seq = row.get("sequence") or seqs.get(row["gene_id"], ("", ""))[1]
            add_candidate(row["gene_id"], seq, "hypothetical_annotated",
                          row["contig"], row["start"], row["end"],
                          row["strand"], [f"RHEA:{rid}"])
    funnel["hyp_annotated"] = len(candidates)

    # B) hypotheticals SIN anotar (no en el gene map)
    if args.include_unannotated_hyp:
        if proteins_faa and Path(proteins_faa).exists():
            gm_ids = {row["gene_id"] for rows in genes.values() for row in rows}
            unann = find_unannotated_hypotheticals(
                Path(proteins_faa), gm_ids,
                args.max_unannotated_hyp, args.orf_min_aa)
            added = 0
            for h in unann:
                add_candidate(h["id"], h["seq"], "hypothetical_unannotated",
                              None, 0, 0, None, [])
                added += 1
            funnel["hyp_unannotated"] = added
            log.info(f"[{genome}] hypotheticals SIN anotar (no en gene map): "
                     f"{added} (de {len(unann)} encontradas en {proteins_faa})")
        else:
            log.warning("--include-unannotated-hyp requiere --proteins-dir")
    else:
        log.info(f"[{genome}] {funnel['hyp_annotated']} hipotéticas "
                 f"ANOTADAS (ya en gene map); para incluir sin anotar usa "
                 f"--include-unannotated-hyp")

    # C) ORFs
    if args.orfs:
        if contigs_fna and Path(contigs_fna).exists():
            contigs = read_fasta_dict(Path(contigs_fna))
            orf_pool = []
            for cid, (_desc, dna) in contigs.items():
                unk = sum(1 for ch in dna.upper() if ch not in "ACGT")
                if unk:
                    log.info(f"    contig {cid}: {unk} no-ACGT "
                             f"({100.0*unk/max(1,len(dna)):.1f}%)")
                for orf in find_orfs(dna, args.orf_min_aa, args.max_orfs):
                    oid = f"ORF_{cid}_{orf['start']}_{orf['end']}{orf['strand']}"
                    add_candidate(oid, orf["seq"], "orf", cid, orf["start"],
                                  orf["end"], orf["strand"], [], pool=orf_pool)
            orf_pool.sort(key=lambda c: -len(c["seq"]))
            budget = max(0, args.max_fold - len(candidates))
            added = orf_pool[:budget]
            candidates.extend(added)
            funnel["orfs_added"] = len(added)
        else:
            log.warning("--orfs sin --genome-fna")

    if not candidates:
        # diagnóstico específico
        if funnel["gene_map_empty"]:
            diag = "gene_map vacío — corre 09 o pasa --gene-map"
        elif funnel["hyp_rxn_total"] == 0:
            diag = ("gene_map sin filas 'hypothetical' — ninguna proteína del "
                    "genoma está marcada como hipotética")
        elif funnel["hyp_rxn_relevant"] == 0:
            diag = (f"{funnel['hyp_rxn_total']} hipotéticas anotadas pero "
                    "NINGUNA produce un RHEA de los gaps genuinos — "
                    "considera --include-unannotated-hyp")
        elif not args.include_unannotated_hyp and not args.orfs:
            diag = "usa --include-unannotated-hyp o --orfs para ampliar"
        else:
            diag = "sin candidatas tras aplicar todas las fuentes"
        log.info(f"[{genome}] {diag}")
        log.info(f"[{genome}] Funnel: {funnel}")
        return _empty_result(genome, gaps, total_bp, reduced,
                             context_note, funnel=funnel, diag=diag)

    if len(candidates) > args.max_fold:
        candidates.sort(key=lambda c: -len(c["seq"]))
        log.info(f"[{genome}] {len(candidates)} candidatas -> top "
                 f"{args.max_fold} por longitud")
        candidates = candidates[:args.max_fold]
    funnel["total_candidates"] = len(candidates)

    # referencias por brecha
    gap_refs = {}
    for g in gaps:
        accs = []
        for r in g.get("producing_rhea", [])[:args.ref_rhea]:
            norm = m04.normalize_rhea_id(r)
            if norm:
                accs.extend(rhea2acc.get(norm[1], []))
        if accs:
            gap_refs[g["chebi"]] = accs[:args.ref_per_gap]

    query_faa = work_root / f"{genome}.cand.faa"
    ref_faa = work_root / f"{genome}.ref.faa"
    with open(query_faa, "w") as f:
        for c in candidates:
            f.write(f">{c['id']}\n{c['seq']}\n")
    ref_seqs, ref_of_gap = {}, {}
    for chebi, accs in gap_refs.items():
        for acc in accs:
            if acc in ref_seqs:
                continue
            s = fetch_uniprot_fasta(acc, Path(up_cache_dir))
            if s:
                ref_seqs[acc] = s
                ref_of_gap[acc] = chebi
    with open(ref_faa, "w") as f:
        for acc, s in ref_seqs.items():
            f.write(f">{acc}\n{s}\n")

    sims = similarity_search(query_faa, ref_faa, args.tool,
                             Path(work_root) / f"{genome}_sim") \
        if ref_seqs else []
    best_sim = {}
    for q, t, pid, cov in sims:
        cur = best_sim.get(q)
        if cur is None or (pid, cov) > (cur[1], cur[2]):
            best_sim[q] = (t, pid, cov)

    # log embudo
    by_tri = funnel.get("by_triage", {}) or {}
    tri_str = ", ".join(f"{k}:{v}" for k, v in sorted(by_tri.items())) or "n/a"
    log.info(
        f"[{genome}] Embudo: {funnel['gaps_total']} gaps 04 → "
        f"{funnel.get('after_triage', '?')} genuinos ({tri_str}) → "
        f"{funnel.get('after_cohort', '?')} post-cohort → "
        f"{funnel['hyp_annotated']} hip-anot + "
        f"{funnel['hyp_unannotated']} hip-sin-anot + "
        f"{funnel['orfs_added']} ORFs → {funnel['total_candidates']} cand.")
    log.info(f"[{genome}] Similitud: {len(ref_seqs)} refs, "
             f"{len(best_sim)} candidatas con hit")

    # veredicto + foldseek
    folded = skipped_fold = foldseek_hit_count = 0
    fk_cfg = {
        "enabled": not args.no_foldseek,
        "db": Path(args.foldseek_db) if args.foldseek_db else None,
        "top": args.foldseek_top,
        "tm": args.foldseek_tm,
        "timeout": args.foldseek_timeout,
    }
    if fk_cfg["enabled"] and fk_cfg["db"]:
        if not Path(str(fk_cfg["db"]) + ".dbtype").exists():
            log.warning(f"[{genome}] foldseek DB no encontrada "
                        f"({fk_cfg['db']}.dbtype) — Foldseek desactivado")
            fk_cfg["enabled"] = False

    for c in candidates:
        hit = best_sim.get(c["id"])
        best = None
        if hit:
            t, pid, cov = hit
            best = {"uniprot": t, "pident": round(pid, 1),
                    "coverage": round(cov, 2),
                    "gap_chebi": ref_of_gap.get(t)}
        strong_seq = bool(best and best["pident"] >= 30.0
                          and best["coverage"] >= 0.70)

        if args.strategy == "similarity_first" and strong_seq:
            c.update(verdict="fuerte", folded=False,
                     reason="secuencia fuerte (>=30% id, >=70% cob.); "
                            "plegado omitido",
                     reason_category="sequence_strong",
                     mean_plddt=None, frac_structured=None,
                     structured_run=None, best_hit=best,
                     foldseek_hits=[])
            skipped_fold += 1
            continue

        if len(c["seq"]) > args.max_aa:
            c.update(verdict="descartado", folded=False,
                     reason=f"secuencia > {args.max_aa} aa",
                     reason_category="too_long",
                     mean_plddt=None, frac_structured=None,
                     structured_run=None, best_hit=best,
                     foldseek_hits=[])
            continue

        pdb, pdb_path = esm_fold(c["seq"], Path(esm_cache_dir), args.timeout,
                                 args.retries, args.pause,
                                 provider=args.esm_provider,
                                 seq_id=c["id"], genome=genome,
                                 poll=args.esm_poll,
                                 job_timeout=args.esm_job_timeout)
        if not pdb:
            c.update(verdict="sin plegado", folded=False,
                     reason="API falló", reason_category="api_failed",
                     mean_plddt=None, frac_structured=None,
                     structured_run=None, best_hit=best,
                     foldseek_hits=[])
            continue
        folded += 1
        mean_p, frac, run = plddt_from_pdb(pdb)
        structured = (mean_p is not None and mean_p >= 50 and run >= 80)

        # ---- foldseek (solo si plegado) ----
        fs_hits = []
        if fk_cfg["enabled"] and pdb_path:
            fs_hits = foldseek_search(
                pdb_path, fk_cfg["db"],
                Path(work_root) / f"{genome}_fs",
                fk_cfg["top"], fk_cfg["timeout"])
            if fs_hits:
                foldseek_hit_count += 1
                # enriquece: rheas del acc + si toca un gap
                for h in fs_hits:
                    rheas = acc_to_rheas.get(h["acc"], set())
                    h["rhea_hits"] = [f"RHEA:{r}" for r in sorted(rheas)][:5]
                    h["relevant"] = bool(rheas & action_rheas)
        relevant_fs = [h for h in fs_hits if h["relevant"]
                       and h["qtmscore"] >= fk_cfg["tm"]]
        best_fs = (max(relevant_fs, key=lambda h: h["qtmscore"])
                   if relevant_fs else (fs_hits[0] if fs_hits else None))

        # ---- veredicto ----
        if structured and (strong_seq or relevant_fs):
            verdict, cat = "fuerte", ("structured_and_similar"
                                      if strong_seq else
                                      "structural_hit_experimental")
            reason = ("dominio estructurado + similitud seq"
                      if strong_seq else
                      f"dominio estructurado + Foldseek TM "
                      f"{best_fs['qtmscore']:.2f} vs {best_fs['acc']}")
        elif structured and best_fs and best_fs["qtmscore"] >= fk_cfg["tm"]:
            verdict, cat = "validar", "structural_match"
            reason = f"dominio estructurado + Foldseek TM " \
                     f"{best_fs['qtmscore']:.2f} vs {best_fs['acc']}"
        elif structured:
            verdict, cat = "validar", "structured_no_similar"
            reason = "dominio estructurado, sin similitud fuerte"
        else:
            verdict, cat = "descartado", "low_plddt"
            reason = (f"desorden (pLDDT {mean_p and round(mean_p)}, "
                      f"dominio {run} aa)")

        c.update(verdict=verdict, folded=True, reason=reason,
                 reason_category=cat, mean_plddt=mean_p,
                 frac_structured=frac, structured_run=run, best_hit=best,
                 foldseek_hits=fs_hits)
        fs_log = (f" · Fs {best_fs['acc']} TM {best_fs['qtmscore']:.2f}"
                  if best_fs else "")
        log.info(f"  [{c['id']}] {c['length_bp']} bp · pLDDT "
                 f"{mean_p and round(mean_p)} · dom {run} aa{fs_log} -> "
                 f"{verdict}")

    funnel["folded"] = folded
    funnel["strong_seq_skipped"] = skipped_fold
    funnel["foldseek_hits"] = foldseek_hit_count
    log.info(f"[{genome}] Plegado: {folded} API, {skipped_fold} saltadas "
             f"por secuencia fuerte. Foldseek: {foldseek_hit_count} con hit.")

    # salida
    out = {
        "script": "08_fold_candidates", "script_version": SCRIPT_VERSION,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "genome": genome, "strategy": args.strategy,
        "dna": {"genome_size_bp": total_bp, "reduced_genome": reduced,
                "context": context_note},
        "gaps_total": len(all_gaps),
        "gaps_after_filters": len(gaps),
        "funnel": funnel,
        "candidates_folded": folded,
        "candidates_skipped_fold": skipped_fold,
        "candidates": candidates,
        "gap_reference_uniprot": gap_refs,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    jp = out_dir / f"{genome}_candidates.json"
    jp.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                  encoding="utf-8")
    csv_path = out_dir / f"{genome}_candidates.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "kind", "verdict", "reason", "reason_category",
                    "folded", "seq_notes", "mean_plddt", "structured_run_aa",
                    "best_hit", "pident", "coverage", "gap_chebi",
                    "fs_best_acc", "fs_best_tm", "fs_relevant",
                    "contig", "start", "end", "strand", "length_bp",
                    "pct_of_genome", "rxn_ids", "genome_size_bp",
                    "reduced_genome"])
        for c in candidates:
            bh = c.get("best_hit") or {}
            fs = c.get("foldseek_hits") or []
            fs_best = fs[0] if fs else {}
            w.writerow([c["id"], c["kind"], c["verdict"], c["reason"],
                        c.get("reason_category", ""), c.get("folded", False),
                        c.get("seq_note", ""), c.get("mean_plddt"),
                        c.get("structured_run"),
                        bh.get("uniprot", ""), bh.get("pident", ""),
                        bh.get("coverage", ""), bh.get("gap_chebi", ""),
                        fs_best.get("acc", ""), fs_best.get("qtmscore", ""),
                        fs_best.get("relevant", ""),
                        c.get("contig") or "", c.get("start") or "",
                        c.get("end") or "", c.get("strand") or "",
                        c["length_bp"],
                        round(100.0 * c["length_bp"] / total_bp, 3)
                        if total_bp and c["length_bp"] else "",
                        ";".join(c["rxn_ids"]),
                        total_bp or "", reduced])
    log.info(f"[{genome}] {folded}/{len(candidates)} plegadas -> {jp.name}, "
             f"{csv_path.name}")
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Criba de enzimas faltantes: hipotéticas + ORFs + ESMFold "
                    "+ Foldseek (v1.3).")
    ap.add_argument("--genome", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--net-dir", type=Path, default=None)
    ap.add_argument("--rhea-dir", type=Path, default=None)
    ap.add_argument("--gene-map", type=Path, default=None)
    ap.add_argument("--proteins-dir", type=Path, default=None)
    ap.add_argument("--genome-meta", type=Path, default=None)
    ap.add_argument("--genomes-dir-fna", type=Path, default=None)
    ap.add_argument("--orfs", action="store_true")
    ap.add_argument("--orf-min-aa", type=int, default=80)
    ap.add_argument("--max-orfs", type=int, default=20)
    ap.add_argument("--max-fold", type=int, default=50)
    ap.add_argument("--max-aa", type=int, default=1200)
    ap.add_argument("--ref-rhea", type=int, default=3)
    ap.add_argument("--ref-per-gap", type=int, default=5)
    ap.add_argument("--tool", default="auto",
                    choices=["auto", "mmseqs", "blastp", "none"])
    ap.add_argument("--esm-provider", default="auto",
                    choices=["auto", "neurosnap", "direct"])
    ap.add_argument("--esm-poll", type=float, default=5.0)
    ap.add_argument("--esm-job-timeout", type=int, default=900)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--pause", type=float, default=1.0)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--env-file", type=Path, default=None)
    ap.add_argument("--strategy", default="similarity_first",
                    choices=["similarity_first", "fold_all"])
    ap.add_argument("--no-use-triage", action="store_true")
    ap.add_argument("--use-cohort", action="store_true",
                    help="excluir gaps systematic_gaps (artefactos en TODOS "
                         "los genomas)")
    ap.add_argument("--cohort-specific-only", action="store_true",
                    help="filtro antiguo: solo gaps en specific_gaps (1 genoma)")
    ap.add_argument("--cohort", type=Path, default=None)
    # v1.3
    ap.add_argument("--include-unannotated-hyp", action="store_true",
                    help="pliega hypothetical proteins del .faa que NO están "
                         "en el gene map (sin EC asignada)")
    ap.add_argument("--max-unannotated-hyp", type=int, default=200,
                    help="tope de hypotheticals sin anotar (default 200)")
    ap.add_argument("--foldseek-db", default=None,
                    help="prefijo de la DB Foldseek "
                         "(default data/db/foldseek_swissprot)")
    ap.add_argument("--no-foldseek", action="store_true")
    ap.add_argument("--foldseek-top", type=int, default=3)
    ap.add_argument("--foldseek-tm", type=float, default=0.5)
    ap.add_argument("--foldseek-timeout", type=int, default=60)
    args = ap.parse_args()

    m04.load_env_file(args.env_file)

    net_dir = args.net_dir or m04.NET_DIR
    rhea_dir = args.rhea_dir or m04.RHEA_DIR
    out_dir = args.out_dir or net_dir / "fold_candidates"
    up_cache = rhea_dir / "uniprot_cache"
    esm_cache = net_dir / "esm_cache"
    work_root = net_dir / "fold_work"
    work_root.mkdir(parents=True, exist_ok=True)

    # foldseek DB default
    if args.foldseek_db is None and not args.no_foldseek:
        default_db = m04.PROJECT_ROOT / "data" / "db" / "foldseek_swissprot"
        if Path(str(default_db) + ".dbtype").exists():
            args.foldseek_db = str(default_db)
            log.info(f"Foldseek DB: {default_db}")

    if args.all:
        genomes = sorted(p.name[:-len("_network_trace.json")] for p in
                         net_dir.glob("*_network_trace.json"))
    elif args.genome:
        genomes = [args.genome]
    else:
        ap.error("usa --genome o --all")

    # cohort (05)
    cohort_cfg = None
    if args.use_cohort or args.cohort_specific_only:
        cohort_path = args.cohort or (net_dir / "trace_compare.json")
        if cohort_path.exists():
            try:
                cohort = json.loads(cohort_path.read_text(encoding="utf-8"))
                if args.cohort_specific_only:
                    s = set(cohort.get("specific_gaps", []))
                    mode = "specific_only"
                    log.info(f"Cohorte 05 ({cohort_path.name}): modo "
                             f"specific_only, {len(s)} gaps")
                else:
                    s = set(cohort.get("systematic_gaps", []))
                    mode = "exclude_systematic"
                    log.info(f"Cohorte 05 ({cohort_path.name}): modo "
                             f"exclude_systematic, {len(s)} gaps a excluir")
                cohort_cfg = {"mode": mode, "set": s}
            except Exception as e:
                log.warning(f"No se pudo leer {cohort_path}: {e}")
        else:
            log.warning(f"--use-cohort sin {cohort_path}; se ignora")

    rhea2acc = rhea_to_uniprot(rhea_dir)
    log.info(f"rhea2uniprot: {len(rhea2acc)} RHEA con UniProt")
    acc_to_rheas = defaultdict(set)
    for rid, accs in rhea2acc.items():
        for a in accs:
            acc_to_rheas[a].add(rid)

    summary = []
    for genome in genomes:
        gm = None
        if args.gene_map:
            gm = (args.gene_map / f"{genome}_gene_map.tsv")
            if not gm.exists():
                gm = (args.gene_map / f"{genome}_gene_map.json")
            if not gm.exists() and args.gene_map.suffix in (".tsv", ".json"):
                gm = args.gene_map
        faa = find_seq_file(args.proteins_dir, genome, ".faa")
        fna = find_seq_file(args.genomes_dir_fna, genome, ".fna")
        try:
            res = process_genome(genome, net_dir, rhea_dir, gm, faa,
                                 args.genome_meta, fna, out_dir, up_cache,
                                 esm_cache, work_root, args, rhea2acc,
                                 dict(acc_to_rheas), cohort_cfg=cohort_cfg)
            if res:
                strong = sum(1 for c in res["candidates"]
                             if c["verdict"] == "fuerte")
                valid = sum(1 for c in res["candidates"]
                            if c["verdict"] == "validar")
                dna = res.get("dna") or {}
                summary.append([
                    genome, len(res["candidates"]),
                    res.get("candidates_folded", 0),
                    res.get("candidates_skipped_fold", 0),
                    strong, valid,
                    dna.get("genome_size_bp") or "",
                    dna.get("reduced_genome", False)])
        except Exception as e:
            log.exception(f"[{genome}] falló: {e}")

    if summary:
        sp = out_dir / "all_candidates_summary.csv"
        with open(sp, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["genome", "candidates", "folded", "skipped_fold",
                        "fuerte", "validar", "genome_size_bp",
                        "reduced_genome"])
            w.writerows(summary)
        log.info(f"Resumen del lote -> {sp}")


if __name__ == "__main__":
    main()
