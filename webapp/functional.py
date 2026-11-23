"""Anotación funcional automática y clasificación en tres caminos.

Camino A: hit Swiss-Prot con función nombrada.
Camino B: hit a proteína uncharacterized/DUF/UPF.
Camino C: sin hit Swiss-Prot.

No inventa función: conserva el hit, la identidad, el e-value y la fuente.
Si falta DIAMOND o la base, devuelve un aviso y no rompe la anotación.
"""
from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

RE_OS = re.compile(r"\bOS=(.*?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)")
RE_GN = re.compile(r"\bGN=([^\s]+)")
RE_EC_PAREN = re.compile(r"\(EC[ :]*([\d\.]+\.-?[\d\-\.]*)\)")
RE_EC_EQ = re.compile(r"\bEC[=:]\s*([\d\.]+\.-?[\d\-\.]*)")

SUBUNIT_NAME = re.compile(r"\bsubunit\b|\bcomponent\b", re.I)

UNCHARACTERIZED = re.compile(
    r"uncharacterized|unknown|duf\d+|upf\d+|putative|hypothetical|"
    r"protein of unknown function|y[a-z]{2,4}_",
    re.I,
)


@dataclass
class Hit:
    cds_id: str
    uniprot: str = ""
    entry_name: str = ""
    gene: str = ""
    product: str = ""
    organism: str = ""
    ec: str = ""
    pident: float = 0.0
    aln_len: int = 0
    evalue: float = 1.0
    bitscore: float = 0.0
    path: str = "C"
    next_step: str = ""
    label: str = ""
    evidence_strength: str = "none"
    evidence_toward: str = "unknown"
    subunit_hint: str = "none"
    subunit_evidence: str = ""
    interpretation: str = ""


def parse_hit(qseqid: str, sseqid: str, stitle: str, pident: float,
              aln_len: int, evalue: float, bitscore: float) -> Hit:
    hit = Hit(cds_id=qseqid, pident=pident, aln_len=aln_len,
              evalue=evalue, bitscore=bitscore)
    parts = sseqid.split("|")
    if len(parts) >= 3:
        hit.uniprot, hit.entry_name = parts[1], parts[2]
    else:
        hit.entry_name = sseqid
    m = RE_OS.search(stitle)
    if m:
        hit.organism = m.group(1).strip()
    m = RE_GN.search(stitle)
    if m:
        hit.gene = m.group(1).split(",")[0].strip()
    if not hit.gene and "_" in hit.entry_name:
        hit.gene = hit.entry_name.split("_")[0].lower()
    cut = stitle.find(" OS=")
    hit.product = stitle[:cut].strip() if cut > 0 else stitle.strip()
    m = RE_EC_PAREN.search(stitle) or RE_EC_EQ.search(stitle)
    if m:
        hit.ec = m.group(1)
    if SUBUNIT_NAME.search(hit.product):
        hit.subunit_hint = "named_subunit"
        hit.subunit_evidence = "product_name_only"
    if UNCHARACTERIZED.search(hit.product) or UNCHARACTERIZED.search(hit.entry_name):
        hit.path = "B"
        hit.label = "B_uncharacterized_homolog"
        hit.evidence_toward = "uncharacterized_family"
        hit.next_step = "dominios/TM/contexto + estructura selectiva"
        hit.interpretation = ("Homología a proteína no caracterizada/DUF: "
                              "evidencia de familia, no de función.")
    else:
        hit.path = "A"
        hit.label = "A_named_function"
        hit.evidence_toward = "named_function"
        hit.next_step = "función nombrada; no requiere estructura"
        hit.interpretation = ("Función nombrada por homología Swiss-Prot; "
                              "no confirmada experimentalmente.")
    return hit


def _evidence(hit: Hit, query_len: int) -> tuple[str, str, str]:
    """Devuelve (label, fuerza, interpretación) de forma conservadora."""
    cov = min(1.0, hit.aln_len / max(1, query_len))
    if hit.path == "A":
        strength = "strong" if hit.pident >= 70 and cov >= 0.80 else "moderate"
        label = "A_named_function"
        toward = "named_function"
        text = ("Función nombrada por homología; fuerte solo como homología, "
                "no como función confirmada.")
    elif hit.path == "B":
        strength = "moderate" if hit.pident >= 60 and cov >= 0.70 else "weak"
        label = "B_uncharacterized_homolog"
        toward = "uncharacterized_family"
        text = ("Homología a proteína no caracterizada/DUF; la función sigue "
                "sin resolverse.")
    else:
        strength = "none"
        label = "C_no_swissprot_hit"
        toward = "unknown"
        text = ("Sin hit Swiss-Prot; no prueba ausencia de función. "
                "Requiere dominio/contexto/estructura.")
    return label, strength, f"{text} (pident={hit.pident:.0f}%, cov={cov:.0%})"


def _row_for(protein: dict[str, Any], hit: Hit | None) -> dict[str, Any]:
    seq = str(protein.get("aa", ""))
    if hit is None:
        return {
            "protein_id": protein["id"], "path": "C",
            "label": "C_no_swissprot_hit", "evidence_strength": "none",
            "evidence_toward": "unknown", "subunit_hint": "none",
            "subunit_evidence": "",
            "interpretation": ("Sin hit Swiss-Prot; no prueba ausencia de "
                               "función. Requiere dominio/contexto/estructura."),
            "uniprot": "", "gene": "", "product": "", "ec": "", "organism": "",
            "pident": "", "coverage": "", "evalue": "",
            "next_step": "dominios/TM/contexto/estructura selectiva",
        }
    label, strength, text = _evidence(hit, len(seq))
    return {
        "protein_id": protein["id"], "path": hit.path, "label": label,
        "evidence_strength": strength, "evidence_toward": hit.evidence_toward,
        "subunit_hint": hit.subunit_hint, "subunit_evidence": hit.subunit_evidence,
        "interpretation": text, "uniprot": hit.uniprot, "gene": hit.gene,
        "product": hit.product, "ec": hit.ec, "organism": hit.organism,
        "pident": hit.pident,
        "coverage": round(min(1.0, hit.aln_len / max(1, len(seq))), 3),
        "evalue": hit.evalue, "next_step": hit.next_step,
    }


def _diamond_db(db: Path) -> Path:
    if db.suffix == ".dmnd":
        return db
    return Path(str(db) + ".dmnd")


def run_diamond(proteins_faa: Path, db: Path, out_tsv: Path,
                threads: int = 4, evalue: str = "1e-5") -> None:
    exe = shutil.which("diamond")
    if exe is None:
        raise RuntimeError("diamond no está en PATH")
    db_path = _diamond_db(db)
    if not db_path.exists():
        raise RuntimeError(f"no existe la base DIAMOND: {db_path}")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        exe, "blastp", "-d", str(db_path), "-q", str(proteins_faa),
        "-o", str(out_tsv), "--outfmt", "6", "qseqid", "sseqid",
        "pident", "length", "evalue", "bitscore", "stitle",
        "--evalue", evalue, "--max-target-seqs", "1",
        "--threads", str(max(1, threads)), "--quiet",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if res.returncode != 0:
        raise RuntimeError((res.stderr or res.stdout or "DIAMOND falló")[-500:])


def parse_m8(path: Path) -> dict[str, Hit]:
    hits: dict[str, Hit] = {}
    if not path.exists():
        return hits
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 7:
            continue
        try:
            hit = parse_hit(cols[0], cols[1], cols[6], float(cols[2]),
                            int(cols[3]), float(cols[4]), float(cols[5]))
        except ValueError:
            continue
        if hit.cds_id not in hits or hit.bitscore > hits[hit.cds_id].bitscore:
            hits[hit.cds_id] = hit
    return hits


def annotate(proteins: list[dict[str, Any]], data_root: Path, outdir: Path,
             genome: str, threads: int = 4) -> dict[str, Any]:
    """Anota proteínas y escribe functional.tsv/json + protein_paths.tsv."""
    outdir.mkdir(parents=True, exist_ok=True)
    faa = outdir / f"{genome}.proteins.faa"
    faa.write_text("".join(f">{p['id']}\n{p['aa']}\n" for p in proteins),
                   encoding="utf-8")
    db = data_root / "db" / "swissprot.dmnd"
    raw = outdir / f"{genome}.functional.tsv"
    result: dict[str, Any] = {
        "ok": False, "n_proteins": len(proteins), "n_hits": 0,
        "warnings": [], "hits": {}, "paths": {"A": 0, "B": 0, "C": 0},
    }
    try:
        run_diamond(faa, db, raw, threads=threads)
    except Exception as exc:  # noqa: BLE001
        result["warnings"].append(f"anotación funcional omitida: {exc}")
        # sin hit no significa camino C: es "no evaluado"
        result["status"] = "not_evaluated"
        (outdir / f"{genome}.functional_annotation.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        return result

    hits = parse_m8(raw)
    by_id = {p["id"]: p for p in proteins}
    rows = []
    for p in proteins:
        row = _row_for(p, hits.get(p["id"]))
        rows.append(row)
        result["paths"][row["path"]] += 1

    result["rows"] = rows
    result.update(
        ok=True, status="ok", n_hits=len(hits),
        hits={k: asdict(v) for k, v in hits.items()},
        paths=dict(result["paths"]),
        match_rate=round(len(hits) / max(1, len(by_id)), 4),
    )
    (outdir / f"{genome}.functional_annotation.json").write_text(
        json.dumps({"report": {k: v for k, v in result.items()
                               if k not in ("hits", "rows")},
                    "hits": result["hits"]}, ensure_ascii=False, indent=2),
        encoding="utf-8")
    write_paths(rows, outdir, genome)
    return result


def write_paths(rows: list[dict[str, Any]], outdir: Path, genome: str) -> Path:
    """Escribe protein_paths.tsv con todas las columnas presentes."""
    path = outdir / f"{genome}.protein_paths.tsv"
    cols: list[str] = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    cols = cols or ["protein_id", "path"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(rows)
    return path


def gene_positions(features: list[dict[str, Any]],
                   hits: dict[str, dict[str, Any]]) -> dict[str, int]:
    """Devuelve posiciones de genes dnaA/dnaN a partir de la anotación funcional."""
    out: dict[str, int] = {}
    for f in features:
        h = hits.get(f.get("id", ""))
        if not h:
            continue
        gene = str(h.get("gene", "")).lower()
        if gene in ("dnaa", "dnan") and gene not in out:
            out[gene] = int(f.get("start", 0))
    return out
