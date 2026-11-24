"""Auditoría CPU de homología remota para proteínas de los caminos B/C.

"Remota" aquí significa homología distante (UniRef/TrEMBL/InterPro/eggNOG/Pfam),
no un servidor remoto. La auditoría usa solo bases locales con DIAMOND/HMMER:
no hace peticiones de red, no descarga bases, no reclasifica A/B/C y no cierra
función. Solo añade evidencia ``remote_*`` y escribe
``<genoma>.remote_homology.tsv``.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

SOURCES = ("uniref", "trembl", "interpro", "eggnog", "pfam")
SEQUENCE_PRIORITY = ("uniref", "trembl", "eggnog")
PROFILE_PRIORITY = ("pfam", "interpro")
REMOTE_FIELDS = ("remote_status", "remote_source", "remote_hit",
                 "remote_evalue", "remote_note")
TSV_COLUMNS = ("genome", "protein_id", "path", *REMOTE_FIELDS)

EVALUE = "1e-5"
DIAMOND_TIMEOUT = 7200
HMMER_TIMEOUT = 3600
_HMM_COMPANIONS = (".h3f", ".h3i", ".h3m", ".h3p")
_DB_DIR_HINTS = ("uniref", "trembl", "interpro", "eggnog", "pfam",
                 "db", "database")
_SEQ_SUFFIXES = (".fasta", ".fa", ".faa", ".fasta.gz", ".fa.gz")


def _candidate_files(data_root: Path) -> list[Path]:
    """Lista archivos candidatos sin recorrer todo ``data_root``.

    Mira ``data_root/db``, ``data_root`` y ``data_root/databases``; solo baja
    un nivel y únicamente en carpetas con nombre de base de datos. Así no se
    recorren los cientos de carpetas de resultados del proyecto.
    """
    root = Path(data_root)
    bases = [root / "db", root, root / "databases"]
    seen: set[str] = set()
    out: list[Path] = []
    for base in bases:
        if not base.is_dir():
            continue
        try:
            entries = sorted(base.iterdir())
        except OSError:
            continue
        for entry in entries:
            candidates: list[Path] = []
            if entry.is_file():
                candidates.append(entry)
            elif entry.is_dir() and any(h in entry.name.lower()
                                        for h in _DB_DIR_HINTS):
                try:
                    candidates.extend(p for p in sorted(entry.iterdir())
                                      if p.is_file())
                except OSError:
                    continue
            for path in candidates:
                try:
                    key = str(path.resolve())
                except OSError:
                    key = str(path)
                if key not in seen:
                    seen.add(key)
                    out.append(path)
    return out


def _pick(files: list[Path], predicate) -> Path | None:
    for path in files:
        if predicate(path):
            return path
    return None


def _is_seq_file(path: Path) -> bool:
    return path.name.lower().endswith(_SEQ_SUFFIXES)


def _is_hmm(path: Path) -> bool:
    return path.suffix.lower() == ".hmm"


def _pressed(hmm: Path) -> bool:
    return all(Path(str(hmm) + ext).is_file() for ext in _HMM_COMPANIONS)


def _info(path: Path | None, engine: str, role: str,
          note: str = "") -> dict[str, Any]:
    if path is None:
        return {"available": False, "path": "", "kind": "", "engine": "",
                "role": role, "searchable": False, "note": ""}
    kind = engine or "data"
    searchable = False
    if engine == "hmmscan" and not _pressed(path):
        note = "HMM sin hmmpress (.h3f/.h3i/.h3m/.h3p); no se puede usar"
    elif engine and shutil.which(engine) is None:
        note = f"{engine} no está en PATH"
    elif engine:
        searchable = True
    return {"available": True, "path": str(path), "kind": kind,
            "engine": engine, "role": role, "searchable": searchable,
            "note": note}


def _pfam_rank(path: Path) -> tuple[int, str]:
    name = path.name.lower()
    if name.startswith("pfam-a"):
        rank = 0
    elif "pfam-a" in name:
        rank = 1
    elif "pfam" in name:
        rank = 2
    elif "panel_focus" in name:
        rank = 3
    elif "panel" in name:
        rank = 4
    else:
        rank = 5
    return rank, str(path)


def _detect(name: str, files: list[Path]) -> dict[str, Any]:
    low = lambda path: path.name.lower()  # noqa: E731

    if name == "uniref":
        path = _pick(files, lambda p: "uniref" in low(p)
                     and low(p).endswith(".dmnd"))
        if path:
            return _info(path, "diamond", "sequence")
        path = _pick(files, lambda p: "uniref" in low(p) and _is_seq_file(p))
        if path:
            return _info(path, "", "sequence",
                         "UniRef sin índice DIAMOND (.dmnd); no consultada")
        return _info(None, "", "sequence")

    if name == "trembl":
        path = _pick(files, lambda p: "trembl" in low(p)
                     and low(p).endswith(".dmnd"))
        if path:
            return _info(path, "diamond", "sequence")
        path = _pick(files, lambda p: "trembl" in low(p) and _is_seq_file(p))
        if path:
            return _info(path, "", "sequence",
                         "TrEMBL sin índice DIAMOND (.dmnd); no consultada")
        return _info(None, "", "sequence")

    if name == "interpro":
        path = _pick(files, lambda p: "interpro" in low(p) and _is_hmm(p))
        if path:
            return _info(path, "hmmscan", "profile")
        path = _pick(files, lambda p: "interpro" in low(p)
                     and low(p).endswith(".dmnd"))
        if path:
            return _info(path, "diamond", "profile")
        path = _pick(files, lambda p: "interpro" in low(p) and _is_seq_file(p))
        if path:
            return _info(path, "", "profile",
                         "InterPro sin índice consultable; no consultada")
        return _info(None, "", "profile")

    if name == "eggnog":
        path = _pick(files, lambda p: "eggnog" in low(p)
                     and low(p).endswith(".dmnd"))
        if path:
            return _info(path, "diamond", "sequence")
        path = _pick(files, lambda p: "eggnog" in low(p)
                     and (low(p).endswith(".db") or low(p).endswith(".tsv")
                          or low(p).endswith(".tsv.gz")))
        if path:
            return _info(path, "", "sequence",
                         "eggNOG detectado sin índice DIAMOND; no consultado")
        return _info(None, "", "sequence")

    if name == "pfam":
        hmms = [p for p in files if _is_hmm(p)
                and ("pfam" in low(p) or "pfam" in p.parent.name.lower())]
        if hmms:
            hmms.sort(key=_pfam_rank)
            return _info(hmms[0], "hmmscan", "profile")
        return _info(None, "", "profile")

    raise KeyError(name)


def available(data_root: str | Path) -> dict[str, dict[str, Any]]:
    """Detecta las bases locales disponibles. No usa la red.

    Devuelve ``{fuente: info}`` para UniRef, TrEMBL, InterPro, eggNOG y Pfam.
    Cada ``info`` tiene ``available``, ``path``, ``kind``, ``engine``,
    ``role``, ``searchable`` y ``note``. ``searchable`` es ``False`` si falta
    el binario, si el HMM no está presionado o si la base no se puede
    consultar offline (p. ej. un ``eggnog.db`` sin DIAMOND).
    """
    files = _candidate_files(Path(data_root))
    return {name: _detect(name, files) for name in SOURCES}


def _tail(text: str, limit: int = 400) -> str:
    return (text or "").strip()[-limit:]


def _run_diamond(faa: Path, db: Path, out: Path,
                 threads: int) -> dict[str, dict[str, Any]]:
    exe = shutil.which("diamond")
    if exe is None:
        raise RuntimeError("diamond no está en PATH")
    cmd = [
        exe, "blastp", "-d", str(db), "-q", str(faa), "-o", str(out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue",
        "bitscore", "stitle", "--evalue", EVALUE, "--max-target-seqs", "1",
        "--threads", str(max(1, int(threads))), "--quiet",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=DIAMOND_TIMEOUT)
    if proc.returncode != 0:
        raise RuntimeError(_tail(proc.stderr or proc.stdout or "DIAMOND falló"))
    return _parse_diamond(out)


def _parse_diamond(path: Path) -> dict[str, dict[str, Any]]:
    hits: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return hits
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 6:
            continue
        try:
            rec = {"hit": cols[1], "pident": float(cols[2]),
                   "evalue": float(cols[4]), "bitscore": float(cols[5])}
        except ValueError:
            continue
        old = hits.get(cols[0])
        if old is None or (rec["evalue"], -rec["bitscore"]) < (
                old["evalue"], -old["bitscore"]):
            hits[cols[0]] = rec
    return hits


def _has_ga(hmm: Path) -> bool:
    """True si el HMM declara umbrales GA (p. ej. Pfam-A o panel curado)."""
    try:
        with hmm.open("r", encoding="utf-8", errors="replace") as fh:
            for _ in range(200):
                line = fh.readline()
                if not line:
                    break
                if line.startswith("GA "):
                    return True
    except OSError:
        return False
    return False


def _run_hmmscan(faa: Path, hmm: Path, out: Path,
                 threads: int) -> dict[str, dict[str, Any]]:
    exe = shutil.which("hmmscan")
    if exe is None:
        raise RuntimeError("hmmscan no está en PATH")
    cmd = [exe, "--domtblout", str(out), "--noali", "--cpu",
           str(max(1, int(threads)))]
    if _has_ga(hmm):
        cmd.append("--cut_ga")
    else:
        cmd.extend(["-E", EVALUE, "--domE", EVALUE])
    cmd.extend([str(hmm), str(faa)])
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=HMMER_TIMEOUT)
    if proc.returncode != 0:
        raise RuntimeError(_tail(proc.stderr or proc.stdout or "hmmscan falló"))
    return _parse_domtbl(out)


def _parse_domtbl(path: Path) -> dict[str, dict[str, Any]]:
    hits: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return hits
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 23:
            continue
        try:
            evalue = float(fields[12])
            score = float(fields[13])
            full_e = float(fields[6])
        except ValueError:
            continue
        rec = {"hit": fields[1] if fields[1] not in ("", "-") else fields[0],
               "evalue": evalue, "full_evalue": full_e, "score": score,
               "description": " ".join(fields[22:])}
        query = fields[3]
        old = hits.get(query)
        if old is None or (rec["evalue"], -rec["score"]) < (
                old["evalue"], -old["score"]):
            hits[query] = rec
    return hits


def _search(source: str, faa: Path, info: dict[str, Any], tmp: Path,
            threads: int) -> dict[str, dict[str, Any]]:
    db = Path(info["path"])
    engine = info["engine"]
    if engine == "diamond":
        return _run_diamond(faa, db, tmp / f"{source}.m8", threads)
    if engine == "hmmscan":
        return _run_hmmscan(faa, db, tmp / f"{source}.domtbl", threads)
    raise RuntimeError(f"motor no soportado para {source}: {engine!r}")


def _write_faa(seqs: dict[str, str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for pid, seq in seqs.items():
            fh.write(f">{pid}\n{seq}\n")


def _path_of(row: dict[str, Any]) -> str:
    return str(row.get("path", "")).strip().upper()


def _protein_id(row: dict[str, Any]) -> str:
    return str(row.get("protein_id") or row.get("id") or "")


def _clean(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _assign(row: dict[str, Any], status: str, source: str = "", hit: str = "",
            evalue: Any = "", note: str = "") -> None:
    row["remote_status"] = status
    row["remote_source"] = source
    row["remote_hit"] = hit
    row["remote_evalue"] = evalue
    row["remote_note"] = note


def _hit_note(source: str, engine: str, hit: dict[str, Any]) -> str:
    if engine == "diamond":
        detail = (f"pident={hit.get('pident', 0):.1f}%, "
                  f"bits={hit.get('bitscore', 0):.0f}")
        tool = "DIAMOND"
    else:
        detail = _clean(hit.get("description", ""))[:80] or "dominio"
        tool = "hmmscan"
    return (f"{tool} vs {source}: {detail}; homología local, no función; "
            "A/B/C sin cambios")


def _not_searchable_reason(info: dict[str, dict[str, Any]]) -> str:
    reasons: list[str] = []
    for name in SOURCES:
        item = info[name]
        if item["available"]:
            reasons.append(f"{name}: {item['note'] or 'sin motor utilizable'}")
    return "; ".join(reasons)[:300]


def _write_tsv(path: Path, rows: list[dict[str, Any]], genome: str) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t",
                                restval="", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            rec = {"genome": _clean(genome),
                   "protein_id": _clean(_protein_id(row)),
                   "path": _path_of(row)}
            for field in REMOTE_FIELDS:
                rec[field] = _clean(row.get(field, ""))
            writer.writerow(rec)
    return path


def run(rows: list[dict[str, Any]], proteins: list[dict[str, Any]],
        outdir: str | Path, genome: str, data_root: str | Path,
        threads: int = 4) -> dict[str, Any]:
    """Audita homología remota local solo para filas B/C.

    Añade ``remote_status``, ``remote_source``, ``remote_hit``,
    ``remote_evalue`` y ``remote_note`` a las filas B/C (en sitio), escribe
    ``<genoma>.remote_homology.tsv`` y devuelve un resumen. Las filas A no se
    tocan. Si no hay base local, deja ``remote_status="not_available"`` y no
    falla. Nunca cambia ``path``/``label`` ni cierra función.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bc_rows = [row for row in rows if _path_of(row) in ("B", "C")]

    seqs: dict[str, str] = {}
    for protein in proteins:
        pid = protein.get("id") or protein.get("protein_id")
        sequence = protein.get("aa") or protein.get("sequence") or ""
        if pid and sequence:
            seqs[str(pid)] = str(sequence)

    info = available(data_root)
    detected = [name for name in SOURCES if info[name]["available"]]
    searched: list[str] = []
    errors: list[str] = []
    seq_hits: dict[str, dict[str, Any]] = {}
    prof_hits: dict[str, dict[str, Any]] = {}
    seq_source = prof_source = ""

    seq_plan = [name for name in SEQUENCE_PRIORITY
                if info[name]["searchable"] and info[name]["role"] == "sequence"]
    prof_plan = [name for name in PROFILE_PRIORITY
                 if info[name]["searchable"] and info[name]["role"] == "profile"]

    if bc_rows and not detected:
        note = ("sin DB local (UniRef/TrEMBL/InterPro/eggNOG/Pfam); "
                "auditoría no ejecutada")
        for row in bc_rows:
            _assign(row, "not_available", note=note)
    elif bc_rows:
        query_ids: list[str] = []
        for row in bc_rows:
            pid = _protein_id(row)
            if pid in seqs and pid not in query_ids:
                query_ids.append(pid)
        with tempfile.TemporaryDirectory(prefix="remote_homology_") as td:
            tmp = Path(td)
            if query_ids and seq_plan:
                faa = tmp / "query.faa"
                _write_faa({pid: seqs[pid] for pid in query_ids}, faa)
                for source in seq_plan:
                    try:
                        seq_hits = _search(source, faa, info[source], tmp,
                                           threads)
                        seq_source = source
                        searched.append(source)
                        break
                    except Exception as exc:  # noqa: BLE001
                        errors.append(f"{source}: {type(exc).__name__}: {exc}")
            remaining = [pid for pid in query_ids if pid not in seq_hits]
            if remaining and prof_plan:
                faa = tmp / "remaining.faa"
                _write_faa({pid: seqs[pid] for pid in remaining}, faa)
                for source in prof_plan:
                    try:
                        prof_hits = _search(source, faa, info[source], tmp,
                                            threads)
                        prof_source = source
                        searched.append(source)
                        break
                    except Exception as exc:  # noqa: BLE001
                        errors.append(f"{source}: {type(exc).__name__}: {exc}")

        for row in bc_rows:
            pid = _protein_id(row)
            if pid not in seqs:
                _assign(row, "no_sequence",
                        note="sin secuencia AA para esta proteína; no se buscó")
            elif pid in seq_hits:
                hit = seq_hits[pid]
                _assign(row, "hit", seq_source, hit["hit"], hit["evalue"],
                        _hit_note(seq_source, info[seq_source]["engine"], hit))
            elif pid in prof_hits:
                hit = prof_hits[pid]
                _assign(row, "hit", prof_source, hit["hit"], hit["evalue"],
                        _hit_note(prof_source, info[prof_source]["engine"], hit))
            elif searched:
                note = (f"sin hit en {', '.join(searched)}; no prueba ausencia "
                        "de homología")
                if errors:
                    note += f" (fallos: {'; '.join(errors)[-200:]})"
                _assign(row, "no_hit", note=note)
            elif errors:
                _assign(row, "error",
                        note=("búsqueda local falló: "
                              + "; ".join(errors))[-300:])
            else:
                reason = _not_searchable_reason(info) or "motores no disponibles"
                _assign(row, "not_searched",
                        note=("DB local detectada pero sin búsqueda: "
                              + reason)[:300])

    tsv = _write_tsv(outdir / f"{genome}.remote_homology.tsv", bc_rows, genome)
    status_counts: dict[str, int] = {}
    for row in bc_rows:
        status = str(row.get("remote_status", ""))
        status_counts[status] = status_counts.get(status, 0) + 1
    n_searched = status_counts.get("hit", 0) + status_counts.get("no_hit", 0)
    return {
        "ok": True,
        "genome": genome,
        "tsv": str(tsv),
        "columns": list(TSV_COLUMNS),
        "n_bc": len(bc_rows),
        "n_hit": status_counts.get("hit", 0),
        "n_searched": n_searched,
        "n_status": status_counts,
        "searched": searched,
        "sources": info,
        "warnings": errors,
        "rows": bc_rows,
    }
