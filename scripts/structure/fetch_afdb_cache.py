#!/usr/bin/env python3
"""Caché AFDB/PAE por accesión única (reanudable, sin duplicados).

- Reúne accesiones UniProt de los genes de los 400.
- Excluye las que ya tienen PDB/PAE local.
- Descarga PDB y PAE JSON de AlphaFold DB.
- Catálogo TSV con hash, tamaño, estado, origen y duplicados.
- Compresión opcional (gzip) con descompresión transparente al leer.
- ``--limit`` para pilotos; ``--pae`` para bajar también PAE.
- ``--dry-run`` para ver el plan sin red ni escrituras.
- No descarga dos veces la misma accesión.

Ejemplos:
  # Piloto conservador (siempre <= 10 accesiones nuevas)
  python scripts/structure/fetch_afdb_cache.py --limit 10 --pae \
      --workers 1 --delay 1.5

  # Ver el plan sin descargar ni escribir
  python scripts/structure/fetch_afdb_cache.py --limit 10 --dry-run

  # Guardar las descargas nuevas comprimidas (PDB -> .pdb.gz,
  # PAE -> .pae.json.gz). Los archivos ya existentes no se tocan.
  python scripts/structure/fetch_afdb_cache.py --limit 10 --pae --compress
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import re
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# Columnas del catálogo. Las 11 primeras son compatibles con el catálogo
# anterior; el resto son aditivas y los lectores viejos pueden ignorarlas.
CAT_COLS = [
    "accession", "status", "pdb", "pae", "sha256", "pdb_bytes", "pae_bytes",
    "pdb_url", "pae_url", "error", "updated",
    "sequence_sha256", "duplicate_of", "source", "compressed",
]

SOURCE_AFDB = "alphafold_db"
SOURCE_LOCAL = "local_cache"
STATUS_PERMANENT = "not_found"
STATUS_RATE_LIMITED = "rate_limited"

PDB_NAME = "{acc}.pdb"
PAE_NAME = "{acc}.pae.json"

# Traducción mínima de residuos para deduplicar por SEQRES cuando la API no
# aporta la secuencia (por ejemplo, catálogos antiguos sin metadatos).
_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z", "UNK": "X",
}

__all__ = [
    "API",
    "CAT_COLS",
    "apply_dedup",
    "build_local_row",
    "cached_path",
    "fetch",
    "is_permanent",
    "load_catalog",
    "main",
    "metadata_sequence",
    "read_cached_bytes",
    "read_cached_text",
    "sequence_from_pdb",
    "sequence_sha256",
    "write_cached",
    "write_catalog",
]


# ---------------------------------------------------------------------------
# Rutas, compresión y hashes
# ---------------------------------------------------------------------------
def cached_path(targets: Path, acc: str, kind: str = "pdb") -> Path | None:
    """Devuelve el archivo local de ``acc``, comprimido o no.

    ``kind`` es ``"pdb"`` o ``"pae"``. Se prefiere el archivo sin comprimir
    para no cambiar el comportamiento de lectores antiguos. Devuelve ``None``
    si no existe ninguno.
    """
    if kind == "pdb":
        names = (PDB_NAME.format(acc=acc), PDB_NAME.format(acc=acc) + ".gz")
    elif kind == "pae":
        names = (PAE_NAME.format(acc=acc), PAE_NAME.format(acc=acc) + ".gz")
    else:
        raise ValueError(f"tipo de archivo desconocido: {kind!r}")
    for name in names:
        path = Path(targets) / name
        if path.is_file():
            return path
    return None


def read_cached_bytes(path: str | os.PathLike[str]) -> bytes:
    """Lee un archivo local descomprimiéndolo de forma transparente."""
    p = Path(path)
    if p.suffix == ".gz":
        with gzip.open(p, "rb") as fh:
            return fh.read()
    return p.read_bytes()


def read_cached_text(path: str | os.PathLike[str], encoding: str = "utf-8",
                     errors: str = "replace") -> str:
    """Lee un archivo de texto (``.gz`` o plano) y lo decodifica."""
    return read_cached_bytes(path).decode(encoding, errors=errors)


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def write_cached(path: str | os.PathLike[str], data: bytes,
                 compress: bool = False) -> Path:
    """Escribe ``data`` de forma atómica y opcionalmente con gzip.

    Si ``compress`` es verdadero y el nombre no termina en ``.gz``, se añade
    ``.gz`` (``X.pdb`` -> ``X.pdb.gz``). Devuelve la ruta final. Nunca borra ni
    reescribe el archivo sin comprimir de una llamada anterior distinta.
    """
    p = Path(path)
    if compress and p.suffix != ".gz":
        p = p.with_name(p.name + ".gz")
    payload = gzip.compress(data, compresslevel=6, mtime=0) if p.suffix == ".gz" else data
    _atomic_write(p, payload)
    return p


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | os.PathLike[str]) -> str:
    """SHA-256 del contenido útil (descomprimido) del archivo."""
    return sha256_bytes(read_cached_bytes(path))


def normalize_sequence(sequence: Any) -> str:
    """Normaliza una secuencia proteica para comparar (A-Z, mayúsculas)."""
    if not isinstance(sequence, str):
        return ""
    return re.sub(r"[^A-Z]", "", sequence.upper())


def sequence_sha256(sequence: Any) -> str:
    """SHA-256 de la secuencia normalizada; ``""`` si no hay secuencia."""
    normalized = normalize_sequence(sequence)
    if not normalized:
        return ""
    return sha256_bytes(normalized.encode("ascii"))


def sequence_from_pdb(text: str) -> str:
    """Reconstruye la secuencia desde registros ``SEQRES`` de un PDB.

    Es un respaldo para catálogos antiguos sin metadatos de la API. Si el PDB
    no tiene ``SEQRES`` devuelve ``""``. Los residuos desconocidos se escriben
    como ``X``.
    """
    chains: dict[str, list[str]] = {}
    order: list[str] = []
    for line in text.splitlines():
        if not line.startswith("SEQRES") or len(line) < 20:
            continue
        chain = line[11:12].strip() or "_"
        if chain not in chains:
            chains[chain] = []
            order.append(chain)
        block = line[19:70]
        for i in range(0, len(block), 4):
            resname = block[i:i + 4].strip().upper()
            if resname:
                chains[chain].append(_AA3.get(resname, "X"))
    return "".join("".join(chains[chain]) for chain in order)


def _compressed_flag(*paths: Path | None) -> str:
    return "true" if any(p is not None and p.suffix == ".gz" for p in paths) else "false"


def _is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "si", "sí"}


# ---------------------------------------------------------------------------
# Catálogo
# ---------------------------------------------------------------------------
def _blank_row(acc: str) -> dict[str, str]:
    row = {col: "" for col in CAT_COLS}
    row["accession"] = acc
    row["status"] = "pending"
    row["compressed"] = "false"
    return row


def normalize_row(row: Mapping[str, Any]) -> dict[str, str]:
    """Completa columnas faltantes de catálogos antiguos y normaliza flags."""
    out = {col: "" for col in CAT_COLS}
    for col in CAT_COLS:
        value = row.get(col, "")
        out[col] = "" if value is None else str(value)
    out["accession"] = out["accession"].strip()
    out["compressed"] = "true" if _is_true(out["compressed"]) else "false"
    return out


def load_catalog(path: Path) -> dict[str, dict[str, str]]:
    """Carga el catálogo TSV; tolera cabeceras antiguas sin columnas nuevas."""
    path = Path(path)
    out: dict[str, dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if acc:
                out[acc] = normalize_row(row)
    return out


def write_catalog(path: Path, rows: Mapping[str, Mapping[str, Any]]) -> None:
    """Escribe el catálogo de forma atómica, ordenado por accesión."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CAT_COLS, delimiter="\t",
                                extrasaction="ignore", restval="")
        writer.writeheader()
        for acc in sorted(rows):
            writer.writerow(normalize_row(rows[acc]))
    os.replace(tmp, path)


def is_permanent(row: Mapping[str, Any] | None) -> bool:
    """``True`` solo para 404/not_found, que no se reintenta en cada lote."""
    if not row:
        return False
    if str(row.get("status", "")) == STATUS_PERMANENT:
        return True
    err = str(row.get("error", ""))
    return "404" in err and ("HTTP 404" in err or "HTTP Error 404" in err)


def apply_dedup(rows: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Marca ``duplicate_of`` por hash de secuencia (mutación in situ).

    La accesión canónica es la que tiene modelo local (``ok`` o
    ``already_present``); si hay empate se elige la menor lexicográficamente,
    de modo que el resultado es determinista y estable entre ejecuciones. La
    deduplicación por accesión ya la garantiza el diccionario del catálogo.
    """
    groups: dict[str, list[str]] = {}
    for row in rows.values():
        row.setdefault("duplicate_of", "")
    for acc, row in rows.items():
        digest = str(row.get("sequence_sha256") or "").strip()
        if digest:
            groups.setdefault(digest, []).append(acc)

    desired: dict[str, str] = {}
    for accessions in groups.values():
        if len(accessions) < 2:
            continue

        def _rank(acc: str) -> tuple[int, str]:
            row = rows[acc]
            has_model = (row.get("status") in {"ok", "already_present"}
                         and bool(row.get("pdb")))
            return (0 if has_model else 1, acc)

        canonical = min(accessions, key=_rank)
        for acc in accessions:
            desired[acc] = "" if acc == canonical else canonical

    for acc, row in rows.items():
        if acc in desired:
            row["duplicate_of"] = desired[acc]
            continue
        previous = str(row.get("duplicate_of") or "").strip()
        if not previous:
            continue
        target = rows.get(previous)
        same_sequence = (
            target is not None
            and row.get("sequence_sha256")
            and target.get("sequence_sha256") == row.get("sequence_sha256")
        )
        if previous == acc or not same_sequence:
            row["duplicate_of"] = ""
    return rows


# ---------------------------------------------------------------------------
# Secuencias y filas locales
# ---------------------------------------------------------------------------
def metadata_sequence(item: Mapping[str, Any]) -> str:
    """Extrae la secuencia UniProt de los metadatos AFDB, si existe."""
    for key in ("uniprotSequence", "sequence", "seq"):
        value = item.get(key)
        normalized = normalize_sequence(value)
        if normalized:
            return normalized
    return ""


def build_local_row(acc: str, targets: Path,
                    existing: Mapping[str, Any] | None = None) -> dict[str, str] | None:
    """Fila de catálogo para un PDB ya presente, sin tocar la red.

    Reutiliza hash y secuencia del catálogo si están; si faltan, calcula el
    SHA-256 del PDB y la secuencia desde ``SEQRES``. Devuelve ``None`` si no
    hay PDB local.
    """
    pdb = cached_path(targets, acc, "pdb")
    if pdb is None:
        return None
    pae = cached_path(targets, acc, "pae")
    previous = normalize_row(existing or {})
    row = _blank_row(acc)
    row["status"] = "already_present"
    row["pdb"] = str(pdb)
    row["pae"] = str(pae) if pae is not None else ""
    row["pdb_bytes"] = str(pdb.stat().st_size)
    row["pae_bytes"] = str(pae.stat().st_size) if pae is not None else "0"
    row["pdb_url"] = previous.get("pdb_url", "")
    row["pae_url"] = previous.get("pae_url", "")
    row["sha256"] = previous.get("sha256", "") or sha256_file(pdb)
    digest = previous.get("sequence_sha256", "")
    if not digest:
        try:
            digest = sequence_sha256(sequence_from_pdb(read_cached_text(pdb)))
        except OSError:
            digest = ""
    row["sequence_sha256"] = digest
    row["duplicate_of"] = previous.get("duplicate_of", "")
    row["source"] = previous.get("source") or SOURCE_LOCAL
    row["compressed"] = _compressed_flag(pdb, pae)
    row["error"] = ""
    row["updated"] = _now()
    return row


def _merge_local(row: dict[str, str], local: Mapping[str, Any] | None) -> dict[str, str]:
    """Completa una fila de error/404 con los archivos locales que sí existan."""
    if not local:
        return row
    for key in ("pdb", "pae", "sha256", "pdb_bytes", "pae_bytes",
                "sequence_sha256", "source", "duplicate_of"):
        if local.get(key) and not row.get(key):
            row[key] = str(local[key])
    # ``compressed`` arranca en "false"; no debe tapar un local comprimido.
    if local.get("compressed") == "true":
        row["compressed"] = "true"
    return row


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


# ---------------------------------------------------------------------------
# Red, throttling y backoff
# ---------------------------------------------------------------------------
def _http_bytes(url: str, timeout: float) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310
        return response.read()


def _throttle(delay: float) -> None:
    """Pausa configurable antes de cada solicitud de red."""
    if delay and delay > 0:
        time.sleep(delay)


def _parse_retry_after(value: Any) -> float | None:
    """Acepta ``Retry-After`` como segundos o como fecha HTTP."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        seconds = float(text)
    except ValueError:
        try:
            moment = parsedate_to_datetime(text)
            seconds = moment.timestamp() - time.time()
        except (TypeError, ValueError, OverflowError):
            return None
    return max(0.0, seconds)


def backoff_delay(attempt: int, base: float, cap: float,
                  jitter: float = 0.0) -> float:
    """Backoff exponencial con tope y jitter opcional (``attempt`` >= 1)."""
    delay = min(max(0.0, cap), max(0.0, base) * (2 ** max(0, attempt - 1)))
    if jitter > 0 and delay > 0:
        delay = min(max(0.0, cap), delay + random.uniform(0.0, delay * jitter))
    return max(0.0, delay)


def _retry_wait(exc: urllib.error.HTTPError, attempt: int, base: float,
                cap: float, retry_after_max: float,
                jitter: float) -> float:
    headers = getattr(exc, "headers", None)
    retry_after = _parse_retry_after(headers.get("Retry-After")
                                     if headers is not None else None)
    if retry_after is not None:
        return min(retry_after, max(0.0, retry_after_max))
    return backoff_delay(attempt, base, cap, jitter)


# ---------------------------------------------------------------------------
# Descarga
# ---------------------------------------------------------------------------
def _row_ok(acc: str, pdb: Path, pae: Path | None, pdb_url: str,
            pae_url: str, sequence_digest: str,
            existing: Mapping[str, Any]) -> dict[str, str]:
    row = _blank_row(acc)
    row["status"] = "ok"
    row["pdb"] = str(pdb)
    row["pae"] = str(pae) if pae is not None else ""
    row["sha256"] = existing.get("sha256") or sha256_file(pdb)
    row["pdb_bytes"] = str(pdb.stat().st_size)
    row["pae_bytes"] = str(pae.stat().st_size) if pae is not None else "0"
    row["pdb_url"] = pdb_url
    row["pae_url"] = pae_url
    if not sequence_digest:
        # Respaldo si la API no expone la secuencia (o cambia de campo).
        try:
            sequence_digest = sequence_sha256(sequence_from_pdb(read_cached_text(pdb)))
        except OSError:
            sequence_digest = ""
    row["sequence_sha256"] = sequence_digest or existing.get("sequence_sha256", "")
    row["duplicate_of"] = existing.get("duplicate_of", "")
    row["source"] = SOURCE_AFDB
    row["compressed"] = _compressed_flag(pdb, pae)
    row["error"] = ""
    row["updated"] = _now()
    return row


def fetch(acc: str, targets: Path,
          catalog: Mapping[str, Mapping[str, Any]] | None = None,
          timeout: float = 25, retries: int = 3, want_pae: bool = False,
          delay: float = 0.0, compress: bool = False,
          backoff_base: float = 2.0, backoff_max: float = 120.0,
          retry_after_max: float = 120.0,
          backoff_jitter: float = 0.25) -> dict[str, str]:
    """Descarga (o reutiliza) PDB/PAE de una accesión.

    - Si ya existen los archivos requeridos, no hace ninguna petición.
    - ``404`` devuelve ``status="not_found"`` y no se reintenta.
    - ``429`` respeta ``Retry-After`` y aplica backoff exponencial; si se
      agotan los intentos devuelve ``status="rate_limited"`` (temporal).
    - ``compress=True`` solo afecta a archivos nuevos; lo ya descargado se
      reutiliza sin reescribirlo.
    """
    targets = Path(targets)
    existing = dict((catalog or {}).get(acc) or {})
    if not retries or retries < 1:
        retries = 1

    local = build_local_row(acc, targets, existing)
    if local is not None and (not want_pae or local.get("pae")):
        return local
    if local is None and is_permanent(existing):
        # 404 ya conocido: no se vuelve a golpear la API en cada lote.
        return normalize_row(existing)

    pdb_path = Path(local["pdb"]) if local and local.get("pdb") else None
    pae_path = Path(local["pae"]) if local and local.get("pae") else None
    sequence_digest = (local or {}).get("sequence_sha256", "") or existing.get(
        "sequence_sha256", "")
    pdb_url = existing.get("pdb_url", "")
    pae_url = existing.get("pae_url", "")

    # Si solo falta el PAE y el catálogo ya conoce su URL, no hace falta
    # volver a pedir metadatos: reanudación con el mínimo de red posible.
    need_metadata = not (pdb_path is not None and want_pae
                         and pae_path is None and pae_url)

    last_error = ""
    final_status = "error"
    for attempt in range(1, retries + 1):
        try:
            if need_metadata:
                _throttle(delay)
                payload = _http_bytes(API.format(acc=acc), timeout)
                meta = json.loads(payload.decode("utf-8"))
                if not isinstance(meta, list) or not meta:
                    raise ValueError("respuesta AFDB vacía o sin entradas")
                item = meta[0] if isinstance(meta[0], dict) else {}
                pdb_url = str(item.get("pdbUrl") or pdb_url)
                pae_url = str(item.get("paeDocUrl") or item.get("paeUrl")
                              or pae_url)
                sequence = metadata_sequence(item)
                if sequence:
                    sequence_digest = sequence_sha256(sequence)
                need_metadata = False

            if pdb_path is None:
                if not pdb_url:
                    raise ValueError("AFDB sin pdbUrl")
                _throttle(delay)
                data = _http_bytes(pdb_url, timeout * 3)
                pdb_path = write_cached(targets / PDB_NAME.format(acc=acc),
                                        data, compress)
            if want_pae and pae_path is None:
                if not pae_url:
                    raise ValueError("AFDB sin paeDocUrl")
                _throttle(delay)
                data = _http_bytes(pae_url, timeout * 3)
                pae_path = write_cached(targets / PAE_NAME.format(acc=acc),
                                        data, compress)
            hint = dict(existing)
            if local:
                for key, value in local.items():
                    if value:
                        hint[key] = value
            else:
                # El PDB se acaba de descargar: no reutilizar el hash de un
                # catálogo que apuntaba a un archivo ya desaparecido.
                hint.pop("sha256", None)
            return _row_ok(acc, pdb_path, pae_path, pdb_url, pae_url,
                           sequence_digest, hint)
        except urllib.error.HTTPError as exc:
            last_error = f"HTTPError: HTTP {exc.code}: {exc.reason}"
            if exc.code == 404:
                row = _blank_row(acc)
                row["status"] = STATUS_PERMANENT
                row["error"] = last_error
                row["pdb_url"] = pdb_url
                row["pae_url"] = pae_url
                row["source"] = SOURCE_AFDB
                row["updated"] = _now()
                # Si el PDB se escribió y el 404 fue del PAE, conservar la ruta.
                local = build_local_row(acc, targets, existing) or local
                return _merge_local(row, local)
            if exc.code == 429:
                final_status = STATUS_RATE_LIMITED
                if attempt >= retries:
                    break
                time.sleep(_retry_wait(exc, attempt, backoff_base, backoff_max,
                                       retry_after_max, backoff_jitter))
                continue
            if 500 <= exc.code < 600:
                final_status = "error"
                if attempt >= retries:
                    break
                time.sleep(_retry_wait(exc, attempt, backoff_base, backoff_max,
                                       retry_after_max, backoff_jitter))
                continue
            # Otros 4xx (400, 401, 403...) no se arreglan reintentando.
            final_status = "error"
            break
        except (urllib.error.URLError, TimeoutError, OSError, ValueError,
                json.JSONDecodeError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            final_status = "error"
            if attempt >= retries:
                break
            time.sleep(backoff_delay(attempt, backoff_base, backoff_max,
                                     backoff_jitter))

    row = _blank_row(acc)
    row["status"] = final_status
    row["error"] = last_error
    row["pdb_url"] = pdb_url
    row["pae_url"] = pae_url
    row["sequence_sha256"] = sequence_digest
    row["updated"] = _now()
    row["source"] = SOURCE_AFDB
    # Un fallo del PAE no debe ocultar el PDB que sí se alcanzó a guardar.
    local = build_local_row(acc, targets, existing) or local
    return _merge_local(row, local)


# ---------------------------------------------------------------------------
# Entrada y orquestación
# ---------------------------------------------------------------------------
def accessions(data_root: Path) -> list[str]:
    """Accesiones únicas, ordenadas, de ``data/regulation/*/genes.tsv``."""
    seen: set[str] = set()
    for path in sorted((Path(data_root) / "regulation").glob("*/genes.tsv")):
        try:
            with path.open(encoding="utf-8", errors="replace", newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = (row.get("acc") or "").strip()
                    if acc:
                        seen.add(acc)
        except OSError:
            continue
    return sorted(seen)


def _missing(acc: str, targets: Path, want_pae: bool) -> bool:
    if cached_path(targets, acc, "pdb") is None:
        return True
    return want_pae and cached_path(targets, acc, "pae") is None


def _print_summary(catalog_path: Path, rows: Mapping[str, Mapping[str, Any]],
                   pending_total: int) -> None:
    counts = Counter(str(row.get("status", "")) for row in rows.values())
    duplicates = sum(1 for row in rows.values() if row.get("duplicate_of"))
    print(
        f"catálogo={catalog_path} ok={counts.get('ok', 0)} "
        f"already_present={counts.get('already_present', 0)} "
        f"error={counts.get('error', 0)} "
        f"rate_limited={counts.get('rate_limited', 0)} "
        f"not_found={counts.get(STATUS_PERMANENT, 0)} "
        f"duplicados={duplicates} pendientes_lote={pending_total}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Caché AFDB/PAE reanudable y deduplicada por accesión.")
    parser.add_argument("--data-root", type=Path, default=ROOT / "data")
    parser.add_argument("--targets", type=Path,
                        default=ROOT / "data" / "metabolism" / "targets" / "afdb25k")
    parser.add_argument("--catalog", type=Path,
                        default=ROOT / "data" / "metabolism" / "targets"
                        / "afdb_catalog.tsv")
    parser.add_argument("--limit", type=int, default=0,
                        help="máximo de accesiones nuevas en este lote; 0 = todas")
    parser.add_argument("--pae", action="store_true", help="descargar PAE JSON")
    parser.add_argument("--pae-only", action="store_true",
                        help="solo PDB ya presentes a los que falta PAE")
    parser.add_argument("--compress", action="store_true",
                        help="guardar descargas nuevas como .pdb.gz/.json.gz")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--delay", type=float, default=0.5,
                        help="segundos antes de cada solicitud por worker")
    parser.add_argument("--timeout", type=float, default=25)
    parser.add_argument("--retries", type=int, default=3,
                        help="intentos totales por accesión")
    parser.add_argument("--backoff-base", type=float, default=2.0)
    parser.add_argument("--backoff-max", type=float, default=120.0)
    parser.add_argument("--retry-after-max", type=float, default=120.0,
                        help="tope al valor de Retry-After")
    parser.add_argument("--flush-every", type=int, default=500,
                        help="guardar el catálogo cada N accesiones (def. 500)")
    parser.add_argument("--dry-run", action="store_true",
                        help="mostrar el plan sin red ni escrituras")
    args = parser.parse_args(argv)
    if args.pae_only:
        args.pae = True

    if args.limit < 0:
        parser.error("--limit debe ser >= 0")
    if args.workers < 1:
        parser.error("--workers debe ser >= 1")
    if args.retries < 1:
        parser.error("--retries debe ser >= 1")
    if args.delay < 0:
        parser.error("--delay debe ser >= 0")
    if args.flush_every < 1:
        parser.error("--flush-every debe ser >= 1")

    catalog = load_catalog(args.catalog)
    accs = accessions(args.data_root)

    if not args.dry_run:
        # Reconciliación local sin red: completa hash/secuencia y detecta
        # archivos ya presentes, incluso si el catálogo estaba incompleto.
        for acc in accs:
            if cached_path(args.targets, acc, "pdb") is not None:
                local = build_local_row(acc, args.targets, catalog.get(acc))
                if local is not None:
                    catalog[acc] = local
        apply_dedup(catalog)
        write_catalog(args.catalog, catalog)

    if args.pae_only:
        pending = [acc for acc in accs
                   if cached_path(args.targets, acc, "pdb") is not None
                   and cached_path(args.targets, acc, "pae") is None
                   and not is_permanent(catalog.get(acc))]
    else:
        pending = [acc for acc in accs
                   if _missing(acc, args.targets, args.pae)
                   and not is_permanent(catalog.get(acc))]
    pending_total = len(pending)
    if args.limit:
        pending = pending[:args.limit]

    print(f"accesiones totales={len(accs)} pendientes={pending_total} "
          f"lote={len(pending)} pae={args.pae} pae_only={args.pae_only} "
          f"workers={args.workers} delay={args.delay}")

    if args.dry_run:
        for acc in pending:
            print(f"  pendiente {acc}")
        print("dry-run: no se descarga nada ni se escribe el catálogo")
        _print_summary(args.catalog, catalog, len(pending))
        return 0

    done = 0
    if pending:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    fetch, acc, args.targets, catalog, args.timeout,
                    args.retries, args.pae, args.delay, args.compress,
                    args.backoff_base, args.backoff_max, args.retry_after_max,
                ): acc
                for acc in pending
            }
            for future in as_completed(futures):
                acc = futures[future]
                try:
                    row = future.result()
                except Exception as exc:  # noqa: BLE001 - no abortar el lote
                    row = _blank_row(acc)
                    row["status"] = "error"
                    row["error"] = f"{type(exc).__name__}: {exc}"
                    row["updated"] = _now()
                catalog[acc] = row
                done += 1
                if done % args.flush_every == 0 or done == len(pending):
                    # Escritura periódica sin recalcular dedup en cada fila;
                    # los archivos ya están en disco y son reanudables.
                    write_catalog(args.catalog, catalog)
                    print(f"  {done}/{len(pending)}")

    apply_dedup(catalog)
    write_catalog(args.catalog, catalog)
    _print_summary(args.catalog, catalog, pending_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
