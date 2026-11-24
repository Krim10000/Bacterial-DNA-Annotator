"""Estado estructural del camino B.

No pliega ni filtra por defecto: busca modelos AFDB locales y deja un
manifiesto. Un modelo disponible no es evidencia de función ni de complejo;
solo habilita el siguiente paso (filtros/estructura).

``run`` conserva el contrato histórico (manifiesto + ``structural_status``).
``evaluate`` es el paso explícito de filtros: solo evalúa filas B con modelo
AFDB local y, sin una base Foldseek local, deja las filas como no evaluadas en
lugar de inventar resultados negativos.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

BOOL_FILTER_FIELDS = (
    "filter_f1_foldseek",
    "filter_f2_enzymm",
    "filter_f3_pocket",
    "filter_f4_alphafill",
    "filter_f5_active_site",
)
FILTER_FIELDS = (
    "filter_status",
    *BOOL_FILTER_FIELDS,
    "filter_verdict",
    "filter_limitation",
    "filter_note",
)
FILTER_COLUMNS = (
    "genome",
    "protein_id",
    "uniprot",
    "pdb",
    *FILTER_FIELDS,
)

NO_LOCAL_DB_NOTE = (
    "foldseek_db no existe; no se ejecutaron F1-F5 ni APIs externas. "
    "Las filas quedan sin evaluar: esto no es evidencia negativa."
)
NO_LOCAL_DB_LIMITATION = (
    "Sin base Foldseek local no hay evaluación estructural; no se puede "
    "afirmar presencia ni ausencia de señal."
)
NO_MODEL_LIMITATION = (
    "Sin modelo AFDB local no hay estructura para evaluar; no implica "
    "ausencia de función."
)
MAX_MODELS_LIMITATION = (
    "No evaluado por el límite max_models; no implica ausencia de señal."
)
FILTER_ERROR_LIMITATION = (
    "Fallo técnico durante la evaluación; no interpretar como resultado "
    "negativo."
)


def _model_dir(data_root: Path) -> Path:
    return Path(data_root) / "metabolism" / "targets" / "afdb25k"


def available(data_root: Path) -> tuple[bool, str]:
    d = _model_dir(data_root)
    if d.exists() and any(d.glob("*.pdb")):
        return True, ""
    return False, f"sin modelos AFDB en {d}"


def run(rows: list[dict], outdir: Path, genome: str,
        data_root: Path) -> dict[str, Any]:
    ok, why = available(data_root)
    models = _model_dir(data_root)
    manifest = []
    for r in rows:
        if r.get("path") != "B":
            continue
        acc = str(r.get("uniprot") or "").split("|")[-1]
        pdb = models / f"{acc}.pdb"
        if acc and pdb.exists():
            r["structural_status"] = "afdb_model_available"
            r["structural_evidence"] = "local_afdb_model"
            r["structural_note"] = ("Modelo AFDB disponible; no evaluado con "
                                    "filtros, no confirma función ni complejo.")
            manifest.append({"genome": genome, "protein_id": r.get("protein_id", ""),
                             "uniprot": acc, "pdb": str(pdb)})
        else:
            r["structural_status"] = "no_afdb_model"
            r["structural_evidence"] = "none"
            r["structural_note"] = ("Sin modelo AFDB local; se necesita descarga "
                                    "o plegado para evaluar estructura.")
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{genome}.structure_b_manifest.tsv"
    cols = ["genome", "protein_id", "uniprot", "pdb"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader(); w.writerows(manifest)
    return {"ok": ok, "reason": why, "n_available": len(manifest),
            "manifest": str(path)}


def _local_db(foldseek_db: str | Path | None) -> Path | None:
    """Devuelve la base local si existe (archivo o prefijo de base Foldseek)."""
    if foldseek_db is None or str(foldseek_db).strip() == "":
        return None
    p = Path(foldseek_db)
    if p.exists():
        return p
    # Foldseek usa un prefijo: ``db`` puede no existir pero sí ``db.dbtype``.
    # Se exigen sufijos propios de una base para no confundir un nombre
    # parecido (p. ej. ``db_backup``) con una base disponible.
    try:
        for suffix in (".dbtype", ".index", ".lookup", ".source"):
            if p.with_name(p.name + suffix).exists():
                return p
    except OSError:
        return None
    return None


def _uniprot(row: dict) -> str:
    return str(row.get("uniprot") or "").split("|")[-1]


def _clear_filter_fields(row: dict) -> None:
    for field in FILTER_FIELDS:
        row.pop(field, None)


def _set_not_evaluated(row: dict, status: str, note: str,
                       limitation: str) -> None:
    _clear_filter_fields(row)
    row["filter_status"] = status
    for field in BOOL_FILTER_FIELDS:
        row[field] = ""
    row["filter_verdict"] = "not_evaluated"
    row["filter_limitation"] = limitation
    row["filter_note"] = note


def _tsv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _record(genome: str, row: dict, uniprot: str, pdb: str) -> dict[str, str]:
    rec = {
        "genome": genome,
        "protein_id": str(row.get("protein_id", "")),
        "uniprot": uniprot,
        "pdb": pdb,
    }
    for field in FILTER_FIELDS:
        rec[field] = _tsv_value(row.get(field, ""))
    return rec


def _write_filters(outdir: Path, genome: str,
                   records: list[dict[str, str]]) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{genome}.structure_b_filters.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(FILTER_COLUMNS))
        w.writeheader(); w.writerows(records)
    return path


def _lookup_hits(hits_by: dict, pdb: Path, uniprot: str) -> list[dict]:
    for key in (pdb.stem, pdb.name, str(pdb), uniprot):
        if key in hits_by:
            return hits_by[key] or []
    # Foldseek puede devolver el query con extensión o con ruta.
    for key, value in hits_by.items():
        if Path(str(key)).stem == pdb.stem:
            return value or []
    return []


def _batch_search(pdbs: list[Path], foldseek_db: Path, workdir: Path,
                  cpu: int) -> dict:
    """Importa Foldseek local solo en la rama que realmente lo necesita."""
    try:
        from .foldseek_local import batch_search
    except ImportError:  # el pipeline importa ``structural_b`` como módulo plano
        from foldseek_local import batch_search
    return batch_search(pdbs, foldseek_db, workdir,
                        threads=max(1, int(cpu)), max_seqs=5)


def _run_filters(pdb: Path, genome: str, uniprot: str,
                 foldseek_hits: list[dict]) -> dict:
    """Importa los filtros de forma perezosa para no tocar APIs por defecto."""
    try:
        from .structural_filters import run_filters
    except ImportError:  # el pipeline importa ``structural_b`` como módulo plano
        from structural_filters import run_filters
    return run_filters(pdb, genome=genome, uniprot=uniprot,
                       foldseek_hits=foldseek_hits)


def _apply_filter_result(row: dict, result: dict) -> None:
    _clear_filter_fields(row)
    for field, key in zip(BOOL_FILTER_FIELDS,
                          ("f1_foldseek", "f2_enzymm", "f3_pocket",
                           "f4_alphafill", "f5_active_site")):
        value = result.get(key)
        row[field] = "" if value is None else bool(value)
    error = str(result.get("error") or "")
    row["filter_status"] = ("evaluated_local_foldseek_partial" if error
                            else "evaluated_local_foldseek")
    row["filter_verdict"] = str(result.get("priority") or "unknown")
    limitation = str(result.get("limitation")
                     or "Filtros computacionales; no confirman función.")
    row["filter_limitation"] = (
        limitation + " F1 se resolvió con Foldseek local; F3/F4 dependen de "
        "binarios o API de run_filters."
    )
    note_parts = [
        f"n_pass={result.get('n_pass', '')}/5",
        f"gating={result.get('gating', '')}",
        f"quality_gate={result.get('quality_gate', '')}",
        f"foldseek_hits={result.get('foldseek_top', '')}",
        f"source={result.get('source', '')}",
    ]
    if error:
        note_parts.append(f"error={error}")
    row["filter_note"] = "; ".join(note_parts)


def evaluate(rows: list[dict], outdir: Path, genome: str, data_root: Path,
             foldseek_db: str | Path | None = None, cpu: int = 4,
             max_models: int = 25) -> dict[str, Any]:
    """Evalúa filtros estructurales solo para filas B con modelo AFDB local.

    Contrato:
    - Sin ``foldseek_db`` local no se ejecuta ningún filtro ni API externa:
      las filas quedan con ``filter_status="not_evaluated_no_local_db"``.
    - Con base local se usa ``foldseek_local.batch_search`` y se pasa el
      resultado a ``structural_filters.run_filters(foldseek_hits=...)``.
    - ``max_models`` limita cuántos modelos se evalúan; el excedente queda
      explícitamente como ``not_evaluated_max_models``.

    Nunca convierte una omisión en un resultado negativo: los campos F1-F5
    quedan vacíos cuando no se evaluaron. El manifiesto y ``structural_status``
    de ``run`` no se modifican.
    """
    outdir = Path(outdir)
    models = _model_dir(data_root)
    db = _local_db(foldseek_db)
    limit = max(0, int(max_models))
    b_rows: list[tuple[dict, str, Path]] = []
    with_model: list[tuple[dict, str, Path]] = []
    n_no_model = 0

    for row in rows:
        if row.get("path") != "B":
            continue
        acc = _uniprot(row)
        pdb = models / f"{acc}.pdb"
        b_rows.append((row, acc, pdb))
        if acc and pdb.is_file():
            with_model.append((row, acc, pdb))
        else:
            n_no_model += 1
            _set_not_evaluated(
                row, "no_afdb_model",
                f"sin modelo AFDB local en {pdb}; filtros no ejecutados",
                NO_MODEL_LIMITATION,
            )

    n_evaluated = 0
    n_no_local_db = 0
    n_max_models = 0
    n_errors = 0

    if db is None:
        n_no_local_db = len(with_model)
        note = (NO_LOCAL_DB_NOTE
                + (f" (ruta solicitada: {foldseek_db})" if foldseek_db else ""))
        for row, acc, pdb in with_model:
            _set_not_evaluated(row, "not_evaluated_no_local_db", note,
                               NO_LOCAL_DB_LIMITATION)
    elif limit == 0:
        n_max_models = len(with_model)
        for row, acc, pdb in with_model:
            _set_not_evaluated(
                row, "not_evaluated_max_models",
                "no evaluado: max_models=0",
                MAX_MODELS_LIMITATION,
            )
    else:
        selected = with_model[:limit]
        overflow = with_model[limit:]
        n_max_models = len(overflow)
        workdir = outdir / f"{genome}.structure_b_foldseek_work"
        if selected:
            try:
                batch = _batch_search([pdb for _, _, pdb in selected], db,
                                      workdir, cpu)
            except Exception as exc:  # noqa: BLE001
                batch = {"hits": {}, "error": f"{type(exc).__name__}: {exc}"}
        else:
            batch = {"hits": {}, "error": ""}
        batch_error = str(batch.get("error") or "")
        hits_by = batch.get("hits") or {}
        if batch_error:
            n_errors = len(selected)
            for row, acc, pdb in selected:
                _set_not_evaluated(
                    row, "not_evaluated_foldseek_error",
                    f"Foldseek local falló: {batch_error}; filtros no ejecutados",
                    FILTER_ERROR_LIMITATION,
                )
        else:
            for row, acc, pdb in selected:
                hits = _lookup_hits(hits_by, pdb, acc)
                try:
                    result = _run_filters(pdb, genome, acc, hits)
                except Exception as exc:  # noqa: BLE001
                    n_errors += 1
                    _set_not_evaluated(
                        row, "filter_error",
                        f"run_filters falló: {type(exc).__name__}: {exc}",
                        FILTER_ERROR_LIMITATION,
                    )
                else:
                    _apply_filter_result(row, result)
                    n_evaluated += 1
        for row, acc, pdb in overflow:
            _set_not_evaluated(
                row, "not_evaluated_max_models",
                f"no evaluado: max_models={limit}",
                MAX_MODELS_LIMITATION,
            )

    records = [_record(genome, row, acc, str(pdb))
               for row, acc, pdb in b_rows]
    path = _write_filters(outdir, genome, records)
    return {
        "ok": True,
        "filters": str(path),
        "path": str(path),
        "foldseek_db": str(db) if db is not None else "",
        "network": "not_used" if db is None else "delegated_to_run_filters",
        "n_b_with_model": len(with_model),
        "n_evaluated": n_evaluated,
        "n_no_model": n_no_model,
        "n_no_local_db": n_no_local_db,
        "n_skipped_max_models": n_max_models,
        "n_errors": n_errors,
    }
