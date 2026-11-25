#!/usr/bin/env python3
"""Descarga FASTA/GFF de los genomas de reference_genomes.json desde NCBI.

- Busca un ensamblado RefSeq completo; prefiere "reference genome" y la cepa.
- No descarga dos veces el mismo accession.
- Guarda en data/reference_genomes/<id>/.
- Actualiza el manifiesto con accession, rutas y estado downloaded.
- Reanudable, --limit, --level, --dry-run.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data" / "reference_genomes.json"
OUT = ROOT / "data" / "reference_genomes"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get_json(url: str, timeout: int = 30) -> dict:
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(min(20, 4 * (attempt + 1)))
    raise last  # type: ignore[misc]


def _taxid(species: str) -> str | None:
    term = urllib.parse.quote(f'"{species}"[Scientific Name]')
    data = _get_json(f"{EUTILS}/esearch.fcgi?db=taxonomy&retmode=json&term={term}")
    ids = data.get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None


def _assembly_by_accession(acc: str) -> dict | None:
    term = urllib.parse.quote(f"{acc}[Assembly Accession]")
    data = _get_json(f"{EUTILS}/esearch.fcgi?db=assembly&retmode=json&term={term}")
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    summary = _get_json(f"{EUTILS}/esummary.fcgi?db=assembly&retmode=json&id={ids[0]}")
    for uid in summary.get("result", {}).get("uids", []):
        r = summary["result"].get(uid) or {}
        if str(r.get("assemblyaccession", "")).startswith(("GCF_", "GCA_")):
            return r
    return None


def _assembly(species: str, strain: str, synonyms: list[str] | None = None) -> dict | None:
    tid = None
    for name in [species] + list(synonyms or []):
        tid = _taxid(name)
        if tid:
            break
    if not tid:
        return None
    ids = []
    for label in ("reference genome", "representative genome"):
        term = urllib.parse.quote(
            f'txid{tid}[Organism:exp] AND "{label}"[Filter]')
        data = _get_json(f"{EUTILS}/esearch.fcgi?db=assembly&retmode=json&term={term}&retmax=100")
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            break
    if not ids:
        # último recurso: buscar RefSeq en una página más amplia
        term = urllib.parse.quote(f"txid{tid}[Organism:exp]")
        data = _get_json(f"{EUTILS}/esearch.fcgi?db=assembly&retmode=json&term={term}&retmax=500")
        ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    # esummary acepta varios ids
    summary = _get_json(f"{EUTILS}/esummary.fcgi?db=assembly&retmode=json&id={','.join(ids[:100])}")
    candidates = []
    for uid in ids[:100]:
        r = summary.get("result", {}).get(uid) or {}
        acc = str(r.get("assemblyaccession") or "")
        status = str(r.get("assemblystatus") or "")
        name = str(r.get("assemblyname") or "")
        organism = str(r.get("speciesname") or "")
        ftp = str(r.get("ftppath_refseq") or r.get("ftppath_genbank") or "")
        if not ftp and r.get("ftppath_assembly_rpt"):
            ftp = str(r["ftppath_assembly_rpt"]).rsplit("/", 1)[0]
        if not acc.startswith(("GCF_", "GCA_")) or not ftp:
            continue
        score = 0
        if acc.startswith("GCF_"):
            score += 2
        if status.lower() == "complete genome":
            score += 4
        if "reference genome" in str(r.get("refseq_category") or "").lower():
            score += 6
        if strain and strain.lower() in (name + " " + organism).lower():
            score += 8
        candidates.append((score, acc, r))
    if not candidates:
        return None
    return max(candidates, key=lambda x: (x[0], x[1]))[2]


def _download(url: str, dest: Path, timeout: int = 120) -> None:
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=timeout) as src, tmp.open("wb") as out:
        shutil.copyfileobj(src, out)
    tmp.replace(dest)


def fetch_one(rec: dict, dry_run: bool = False) -> dict:
    rid = rec["id"]
    base = OUT / rid
    meta_path = base / "assembly.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if rec.get("accession") and meta.get("accession") != rec["accession"]:
            meta = {}
    else:
        meta = {}
    if not meta.get("ftp"):
        asm = (_assembly_by_accession(rec["accession"])
               if rec.get("accession")
               else _assembly(rec["species"], rec.get("strain", ""),
                              rec.get("synonyms") or []))
        if not asm:
            return {"id": rid, "status": "not_found", "reason": "sin ensamblado RefSeq"}
        ftp = str(asm.get("ftppath_refseq") or asm.get("ftppath_genbank") or "")
        meta = {"accession": asm.get("assemblyaccession"), "assembly_name": asm.get("assemblyname"),
                "species": asm.get("speciesname"), "status": asm.get("assemblystatus"),
                "ftp": ftp}
        if not dry_run:
            base.mkdir(parents=True, exist_ok=True)
            meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    ftp = str(rec.get("ftp_path") or meta["ftp"]).rstrip("/")
    acc = meta["accession"]
    fna = base / f"{rid}_genomic.fna.gz"
    gff = base / f"{rid}_genomic.gff.gz"
    if dry_run:
        return {"id": rid, "status": "planned", "accession": acc, "reason": "dry-run"}
    base.mkdir(parents=True, exist_ok=True)
    for suffix, dest in (("_genomic.fna.gz", fna), ("_genomic.gff.gz", gff)):
        if dest.exists() and dest.stat().st_size > 1000:
            continue
        _download(f"{ftp}/{acc}_{meta['assembly_name']}{suffix}", dest)
        time.sleep(0.4)
    rec.update({"accession": acc, "status": "downloaded",
                "genome_path": str(fna.relative_to(ROOT)),
                "annotation_path": str(gff.relative_to(ROOT))})
    return {"id": rid, "status": "downloaded", "accession": acc,
            "genome": str(fna), "annotation": str(gff)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", type=int, choices=(1, 2), default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    rows = json.loads(MANIFEST.read_text(encoding="utf-8"))
    todo = [r for r in rows if r.get("status") not in ("available", "annotated")
            and (a.level is None or r.get("level") == a.level)]
    if a.limit:
        todo = todo[:a.limit]
    print(f"reference genomes pendientes={len(todo)} nivel={a.level}")
    for i, rec in enumerate(todo, 1):
        try:
            res = fetch_one(rec, a.dry_run)
        except Exception as exc:  # noqa: BLE001
            res = {"id": rec.get("id"), "status": "error",
                   "reason": f"{type(exc).__name__}: {exc}"}
        print(f"  [{i}/{len(todo)}] {res}")
        if not a.dry_run:
            MANIFEST.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n",
                                encoding="utf-8")
        time.sleep(0.7)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
