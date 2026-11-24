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
LOCAL_ONLY_LIMITATION = (
    "Modo local sin red: F1 Foldseek local evaluado; F2 EnzyMM y F5 sitios "
    "activos locales solo si están disponibles. F3 P2Rank y F4 AlphaFill "
    "quedan sin evaluar; no se consultaron APIs externas. No confirma función "
    "ni cierra brechas."
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


def _enzymm_local(pdb: Path, workdir: Path,
                  cpu: int) -> tuple[bool | None, str, int, str]:
    """F2 local (EnzyMM) sin APIs; ``None`` significa no evaluado."""
    try:
        try:
            from .structural_annotation import scan_pdbs, summarize
        except ImportError:  # pipeline en modo módulo plano
            from structural_annotation import scan_pdbs, summarize
        scan = scan_pdbs([pdb], workdir, jobs=max(1, int(cpu)))
    except Exception as exc:  # noqa: BLE001
        return None, f"EnzyMM no evaluado ({type(exc).__name__}: {exc})", 0, ""
    error = str(scan.get("error") or "")
    if error:
        return None, f"EnzyMM no evaluado: {error}", 0, ""
    hits = summarize(scan.get("hits") or [])
    templates = ";".join(
        str(h.get("template_pdb_id"))
        for h in hits[:5] if h.get("template_pdb_id")
    )
    return bool(hits), "", len(hits), templates


def _active_site_local(data_root: Path, genome: str,
                       uniprot: str) -> tuple[bool | None, str, list[int]]:
    """F5 local desde los mapas de sitios activos; no consulta APIs."""
    if not genome or not uniprot:
        return None, "sin genome/uniprot para consultar sitios activos", []
    base = Path(data_root) / "metabolism" / "targets"
    paths = [
        base / "ACTIVE_SITE_STRUCTURAL_MAP.tsv",
        base / "gap_glycerate" / "ACTIVE_SITE_STRUCTURAL_MAP.tsv",
    ]
    existing = [p for p in paths if p.is_file()]
    if not existing:
        return None, "mapa local de sitios activos ausente", []
    positions: list[int] = []
    for p in existing:
        try:
            lines = p.read_text(encoding="utf-8").splitlines()
        except OSError as exc:
            return None, f"no se pudo leer {p}: {exc}", []
        for line in lines[1:]:
            c = line.split("\t")
            if (len(c) > 10 and c[0] == genome and c[1] == uniprot
                    and c[2] == "Active site"
                    and str(c[10]).lower() == "true"):
                try:
                    positions.append(int(c[8]))
                except (TypeError, ValueError):
                    continue
    positions = sorted(set(positions))
    return bool(positions), "", positions


def _apply_local_only_result(row: dict, pdb: Path, genome: str, uniprot: str,
                             hits: list[dict], data_root: Path, cpu: int,
                             workdir: Path) -> None:
    """F1/F2/F5 locales; F3/F4 quedan vacíos y no se toca ninguna API."""
    _clear_filter_fields(row)
    top = (hits or [{}])[0]
    f2, f2_note, enz_hits, templates = _enzymm_local(pdb, workdir, cpu)
    f5, f5_note, positions = _active_site_local(data_root, genome, uniprot)
    row["filter_status"] = "evaluated_local_no_external"
    row["filter_f1_foldseek"] = bool(hits)
    row["filter_f2_enzymm"] = "" if f2 is None else f2
    row["filter_f3_pocket"] = ""
    row["filter_f4_alphafill"] = ""
    row["filter_f5_active_site"] = "" if f5 is None else f5
    row["filter_verdict"] = "local_only_partial"
    limitation = LOCAL_ONLY_LIMITATION
    if f2 is None or f5 is None:
        missing = []
        if f2 is None:
            missing.append("F2 no disponible")
        if f5 is None:
            missing.append("F5 no disponible")
        limitation += " " + " ".join(missing) + "."
    row["filter_limitation"] = limitation
    note_parts = [
        "modo=local_sin_red",
        f"foldseek_hits={top.get('pdb_id', '')}",
        f"foldseek_bits={top.get('bits', '')}",
        "P2Rank/AlphaFill/APIs=no_consultados",
    ]
    if f2 is None:
        note_parts.append(f"f2_no_evaluado={f2_note}")
    else:
        note_parts.append(f"enzymm_hits={enz_hits}")
        if templates:
            note_parts.append(f"enzymm_templates={templates}")
    if f5 is None:
        note_parts.append(f"f5_no_evaluado={f5_note}")
    else:
        note_parts.append("active_positions=" + ";".join(map(str, positions)))
    row["filter_note"] = "; ".join(note_parts)


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
             max_models: int = 25, external: bool = False) -> dict[str, Any]:
    """Evalúa filtros estructurales solo para filas B con modelo AFDB local.

    Contrato:
    - Sin ``foldseek_db`` local no se ejecuta ningún filtro ni API externa:
      las filas quedan con ``filter_status="not_evaluated_no_local_db"``.
    - Con base local y ``external=False`` (por defecto) solo se usan fuentes
      locales: F1 Foldseek local, F2 EnzyMM y F5 sitios activos locales si
      están disponibles. F3/F4 quedan vacíos y no se consultan P2Rank,
      AlphaFill ni APIs externas.
    - Con ``external=True`` se conserva el comportamiento delegado:
      ``structural_filters.run_filters(foldseek_hits=...)`` (puede usar
      P2Rank/AlphaFill).
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
    n_evaluated_local_only = 0
    n_evaluated_delegated = 0
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
                if external:
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
                        n_evaluated_delegated += 1
                else:
                    enz_workdir = (outdir / f"{genome}.structure_b_enzymm_work"
                                   / pdb.stem)
                    try:
                        _apply_local_only_result(
                            row, pdb, genome, acc, hits, data_root, cpu,
                            enz_workdir,
                        )
                    except Exception as exc:  # noqa: BLE001
                        n_errors += 1
                        _set_not_evaluated(
                            row, "filter_error",
                            f"filtros locales fallaron: "
                            f"{type(exc).__name__}: {exc}",
                            FILTER_ERROR_LIMITATION,
                        )
                    else:
                        n_evaluated += 1
                        n_evaluated_local_only += 1
        for row, acc, pdb in overflow:
            _set_not_evaluated(
                row, "not_evaluated_max_models",
                f"no evaluado: max_models={limit}",
                MAX_MODELS_LIMITATION,
            )

    records = [_record(genome, row, acc, str(pdb))
               for row, acc, pdb in b_rows]
    path = _write_filters(outdir, genome, records)
    if db is None:
        network = "not_used"
    elif external:
        network = "delegated_to_run_filters"
    else:
        network = "local_only"
    return {
        "ok": True,
        "filters": str(path),
        "path": str(path),
        "foldseek_db": str(db) if db is not None else "",
        "external": bool(external),
        "network": network,
        "n_b_with_model": len(with_model),
        "n_evaluated": n_evaluated,
        "n_evaluated_local_only": n_evaluated_local_only,
        "n_evaluated_delegated": n_evaluated_delegated,
        "n_no_model": n_no_model,
        "n_no_local_db": n_no_local_db,
        "n_skipped_max_models": n_max_models,
        "n_errors": n_errors,
    }
