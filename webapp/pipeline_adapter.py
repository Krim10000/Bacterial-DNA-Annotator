"""Adaptador en memoria entre la GUI y el pipeline CLI (`webapp/pipeline.py`).

Ejecuta la misma ruta de capas que el CLI sin pasar por Streamlit ni por un
FASTA temporal: recibe ``(filename, data)``, escribe los mismos artefactos en
``outdir`` y devuelve el resultado enriquecido con el estado de cada capa.

La secuencia es:

    annotate -> TSS/TTS precalculados -> ncRNA -> taxonomía -> functional
    -> dominios C -> contexto C -> B estructural -> complejos
    (firma/modelo/transferencia) -> regulación -> membrana -> oriC
    -> metabolismo opcional -> report

No reimplementa ningún algoritmo: delega en los módulos del proyecto. Una
capa que no puede ejecutarse (falta un binario, una base de datos o un panel)
devuelve ``ok=False`` con ``status`` ``not_available``/``error`` y ``reason``;
no se inventan filas ni archivos para esas capas.

Los artefactos se escriben con ``write_outputs`` del CLI y con las funciones
de escritura de cada capa. ``res["outputs"]`` mapea nombres cortos
(``gff``, ``protein_paths``, ``complexes``, ``ncrna``, ``membrane``, ...) a
las rutas realmente creadas.

TSS/TTS se cargan de tablas locales ya calculadas (``webapp/tss_layer.py``):
no se descargan, no se predicen y no se transfieren coordenadas entre genomas.
Si no hay tabla para el genoma, ``res["tss_tts"]["status"]`` es
``not_available`` y no se escribe ningún TSV de posiciones.
"""
from __future__ import annotations

import gzip
import importlib.util
import shutil
import time
from pathlib import Path
from typing import Any

if __package__:  # importado como paquete: ``webapp.pipeline_adapter``
    from . import oric, tss_layer
    from .annotator import annotate, build_gff, build_report, safe_name
    from .complex_model import run as run_complex_model
    from .complex_transfer import load as load_complex_transfer
    from .complexes import enrich_gaps, enrich_rows, infer as infer_complexes
    from .complexes import write as write_complexes
    from .context_layer import run as run_context
    from .domain_layer import run as run_domains
    from .functional import annotate as annotate_functional
    from .functional import write_paths
    from .membrane_layer import run as run_membrane
    from .metabolism import build_network
    from .ncrna_layer import run as run_ncrna
    from .pipeline import write_outputs
    from .regulation_layer import DEFAULT_MAX_GAP, run as run_regulation
    from .structural_b import evaluate as evaluate_structural_b
    from .structural_b import run as run_structural_b
    from .taxonomy_layer import run as run_taxonomy
else:  # ejecución con ``webapp/`` en sys.path, como el CLI
    import oric
    import tss_layer
    from annotator import annotate, build_gff, build_report, safe_name
    from complex_model import run as run_complex_model
    from complex_transfer import load as load_complex_transfer
    from complexes import enrich_gaps, enrich_rows, infer as infer_complexes
    from complexes import write as write_complexes
    from context_layer import run as run_context
    from domain_layer import run as run_domains
    from functional import annotate as annotate_functional
    from functional import write_paths
    from membrane_layer import run as run_membrane
    from metabolism import build_network
    from ncrna_layer import run as run_ncrna
    from pipeline import write_outputs
    from regulation_layer import DEFAULT_MAX_GAP, run as run_regulation
    from structural_b import evaluate as evaluate_structural_b
    from structural_b import run as run_structural_b
    from taxonomy_layer import run as run_taxonomy

#: Opciones equivalentes a los valores por defecto del CLI. El metabolismo se
#: deja apagado porque es opcional y puede lanzar DIAMOND contra Swiss-Prot.
DEFAULT_OPTIONS: dict[str, Any] = {
    "cds": True,
    "trna": True,
    "rrna": True,
    "ghost": True,
    "deserts": True,
    "esm": False,
    "esm_top": 0,
    "min_region_len": 150,
    "min_aa": 60,
    "ncrna": True,
    "taxonomy": True,
    "functional": True,
    "domains": True,
    "context": True,
    "structural_b": True,
    "structural_b_filters": True,
    "complexes": True,
    "regulation": True,
    "membrane": True,
    "oric": True,
    "metabolic": False,
    "operon_max_gap": DEFAULT_MAX_GAP,
    "reference": "",
    "accession": "",
    "foldseek_local_db": None,
    "structure_b_external": False,
    "panel_hmm": None,
    "transport_db": None,
    "producer_db": None,
}

#: Capas que se pueden desactivar con el nombre ``no_<capa>`` del CLI o de la
#: GUI. ``no_cds``/``no_trna``/... son los flags reales de ``pipeline.py``.
_NEGATIVE_FLAGS = {
    "no_cds": "cds",
    "no_trna": "trna",
    "no_rrna": "rrna",
    "no_ghost": "ghost",
    "no_deserts": "deserts",
    "no_ncrna": "ncrna",
    "no_taxonomy": "taxonomy",
    "no_functional": "functional",
    "no_domains": "domains",
    "no_context": "context",
    "no_structural_b": "structural_b",
    "no_structural_b_filters": "structural_b_filters",
    "no_complexes": "complexes",
    "no_regulation": "regulation",
    "no_membrane": "membrane",
    "no_oric": "oric",
    "no_metabolic": "metabolic",
}

#: Forma mínima de cada capa cuando se desactiva o no está disponible. Evita
#: que la GUI tenga que distinguir entre una lista vacía y una clave ausente.
_LAYER_SHAPES: dict[str, dict[str, Any]] = {
    "ncrna": {"rows": [], "n": 0},
    "taxonomy": {"rows": [], "n_rows": 0},
    "functional": {
        "rows": [], "hits": {}, "paths": {"A": 0, "B": 0, "C": 0},
        "warnings": [], "n_hits": 0, "match_rate": 0.0,
    },
    "domains": {"rows": [], "n": 0},
    "context_c": {"rows": [], "n": 0},
    "structural_b": {"manifest": "", "n_available": 0},
    "structural_b_filters": {
        "n_evaluated": 0, "n_no_model": 0, "n_no_local_db": 0,
    },
    "complexes": {
        "complexes": [], "by_protein": {}, "n_complete": 0,
        "n_incomplete": 0, "n_candidate": 0, "n_ambiguous": 0,
    },
    "regulation": {
        "rows": [], "promoter_rows": [], "n_operons": 0, "n_promoters": 0,
    },
    "tss_tts": {"n_tss": 0, "n_tts": 0, "n_total": 0, "sources": []},
    "membrane": {"rows": [], "n_rows": 0},
    "oric": {},
}

#: Nombres de archivo que produce el CLI y sus capas (con ``{g}`` = genome).
_OUTPUT_NAMES = {
    "gff": "{g}.gff3",
    "proteins": "{g}.proteins.faa",
    "candidates": "{g}.candidates.tsv",
    "report": "{g}.report.md",
    "summary": "{g}.summary.json",
    "protein_paths": "{g}.protein_paths.tsv",
    "functional": "{g}.functional.tsv",
    "functional_json": "{g}.functional_annotation.json",
    "complexes": "{g}.complexes.tsv",
    "gaps": "{g}.gaps.tsv",
    "ncrna": "{g}.ncrna.tsv",
    "membrane": "{g}.membrane.tsv",
    "membrane_gff": "{g}.membrane.gff3",
    "taxonomy": "{g}.taxonomy.tsv",
    "operons": "{g}.operons_auto.tsv",
    "promoters": "{g}.promoters_auto.tsv",
    "tss_tts": "{g}.tss_tts.tsv",
    "domains": "{g}.domains.tsv",
    "context_c": "{g}.context_c.tsv",
    "structural_b_manifest": "{g}.structure_b_manifest.tsv",
    "structural_b_filters": "{g}.structure_b_filters.tsv",
    "oric": "{g}.oric.json",
    "oric_gff": "{g}.oric.gff3",
    "esmfold_closure": "{g}.esmfold_closure.tsv",
    "structural_evidence": "{g}.structural_evidence.tsv",
    "structural_hypotheses": "{g}.structural_hypotheses.tsv",
    "structural_hypotheses_json": "{g}.structural_hypotheses.json",
}

def _as_fasta_bytes(data: bytes | bytearray | memoryview | str) -> bytes:
    """Normaliza la carga a bytes y descomprime gzip por magic number."""
    if isinstance(data, str):
        data = data.encode("utf-8")
    elif isinstance(data, (bytearray, memoryview)):
        data = bytes(data)
    if data[:2] == b"\x1f\x8b":
        return gzip.decompress(data)
    return data


def _failure(exc: BaseException) -> dict[str, Any]:
    """Registro de una capa que lanzó una excepción, sin inventar filas."""
    missing = isinstance(exc, (FileNotFoundError, ImportError, ModuleNotFoundError))
    return {
        "ok": False,
        "status": "not_available" if missing else "error",
        "reason": f"{type(exc).__name__}: {exc}",
    }


def _disabled(reason: str) -> dict[str, Any]:
    return {"ok": False, "status": "disabled", "reason": reason}


def _unavailable(reason: str) -> dict[str, Any]:
    return {"ok": False, "status": "not_available", "reason": reason}


def _normalize(result: Any, default_status: str = "not_available") -> dict[str, Any]:
    """Completa ``ok``/``status``/``reason`` sin reescribir lo que ya trae."""
    if not isinstance(result, dict):
        return {
            "ok": bool(result),
            "status": "ok" if result else default_status,
            "reason": "" if result else f"resultado no válido: {type(result).__name__}",
        }
    result.setdefault("reason", "")
    if "status" not in result:
        result["status"] = "ok" if result.get("ok") else default_status
    if "ok" not in result:
        result["ok"] = result["status"] not in (
            "error", "not_available", "not_evaluated", "unavailable", "disabled")
    return result


def _call(fn, *args, **kwargs) -> dict[str, Any]:
    """Ejecuta una capa real y convierte cualquier excepción en un registro."""
    try:
        return _normalize(fn(*args, **kwargs))
    except Exception as exc:  # noqa: BLE001 - la capa degrada, no tumba el FASTA
        return _failure(exc)


def _has_module(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def _tool_availability() -> dict[str, bool]:
    """Comprueba herramientas sin importar torch/esm (find_spec es barato)."""
    return {
        "pyrodigal": _has_module("pyrodigal"),
        "trna": shutil.which("tRNAscan-SE") is not None,
        # ``run_rrna`` ejecuta barrnap; cmsearch se usa en la capa ncRNA.
        "rrna": shutil.which("barrnap") is not None,
        "esm": _has_module("torch") and _has_module("esm"),
    }


def _normalize_options(options: dict[str, Any] | None) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Aplica defaults, alias y flags negativos; desactiva herramientas ausentes.

    Respeta lo que el usuario pidió, pero si una herramienta externa no está
    instalada apaga esa opción y deja constancia en ``degraded``. Así el
    adaptador no lanza ``ImportError``/``FileNotFoundError`` a la GUI.
    """
    opts = dict(DEFAULT_OPTIONS)
    if options:
        if "max_gap" in options and "operon_max_gap" not in options:
            opts["operon_max_gap"] = options["max_gap"]
        opts.update(options)

    for flag, layer in _NEGATIVE_FLAGS.items():
        if opts.pop(flag, False):
            opts[layer] = False
    # La GUI puede enviar ``esm_top`` aunque ESM esté apagado.
    try:
        esm_top = int(opts.get("esm_top") or 0)
    except (TypeError, ValueError):
        esm_top = 0
    opts["esm_top"] = max(0, esm_top) if opts.get("esm") else 0
    for key, default in (("min_region_len", 150), ("min_aa", 60),
                         ("operon_max_gap", DEFAULT_MAX_GAP)):
        try:
            opts[key] = int(opts.get(key, default))
        except (TypeError, ValueError):
            opts[key] = default
    opts["reference"] = str(opts.get("reference") or "")
    opts["accession"] = str(opts.get("accession") or "")
    for key in ("cds", "trna", "rrna", "ghost", "deserts", "esm", "ncrna",
                "taxonomy", "functional", "domains", "context", "structural_b",
                "structural_b_filters", "complexes", "regulation", "membrane",
                "oric", "metabolic", "structure_b_external"):
        opts[key] = bool(opts.get(key))

    degraded: list[dict[str, str]] = []
    tools = _tool_availability()
    if opts["cds"] and not tools["pyrodigal"]:
        opts["cds"] = False
        degraded.append({
            "layer": "annotate", "option": "cds", "status": "not_available",
            "reason": "pyrodigal no está instalado; CDS desactivado",
        })
    if opts["trna"] and not tools["trna"]:
        opts["trna"] = False
        degraded.append({
            "layer": "annotate", "option": "trna", "status": "not_available",
            "reason": "tRNAscan-SE no está en PATH; tRNA desactivado",
        })
    if opts["rrna"] and not tools["rrna"]:
        opts["rrna"] = False
        degraded.append({
            "layer": "annotate", "option": "rrna", "status": "not_available",
            "reason": "barrnap no está en PATH; rRNA desactivado",
        })
    if opts["esm"] and not tools["esm"]:
        opts["esm"] = False
        opts["esm_top"] = 0
        degraded.append({
            "layer": "annotate", "option": "esm", "status": "not_available",
            "reason": "torch/fair-esm no están instalados; ESM desactivado",
        })
    return opts, degraded


def _failed_annotation(
    filename: str,
    outdir: Path,
    data_root: Path,
    threads: int,
    options: dict[str, Any],
    degraded: list[dict[str, str]],
    exc: BaseException,
) -> dict[str, Any]:
    """Resultado mínimo cuando falla la capa base (sin GFF ni proteínas)."""
    reason = f"{type(exc).__name__}: {exc}"
    record = _failure(exc)
    result: dict[str, Any] = {
        "name": filename,
        "genome": safe_name(filename),
        "ok": False,
        "status": record["status"],
        "reason": reason,
        "error": reason,
        "features": [],
        "proteins": [],
        "sequences": {},
        "bp": 0,
        "gc": 0.0,
        "coverage": 0.0,
        "annotated_bp": 0,
        "candidates": [],
        "deserts": [],
        "stats": {},
        "gff": "",
        "proteins_fasta": "",
        "candidates_tsv": "",
        "report_md": "",
        "outdir": str(outdir),
        "data_root": str(data_root),
        "threads": threads,
        "options": dict(options),
        "degraded": degraded,
        "outputs": {},
    }
    for layer in _LAYER_SHAPES:
        result[layer] = dict(record)
    result["metabolic"] = None
    result["metabolic_status"] = dict(record)
    _ensure_shapes(result)
    result["layer_status"] = {
        layer: record["status"] for layer in (*_LAYER_SHAPES, "metabolic")
    }
    return result


def _ensure_shapes(res: dict[str, Any]) -> None:
    """Rellena listas/contadores vacíos; nunca añade filas inventadas."""
    for layer, defaults in _LAYER_SHAPES.items():
        value = res.get(layer)
        if not isinstance(value, dict):
            value = _unavailable(f"capa sin resultado: {layer}")
            res[layer] = value
        for key, default in defaults.items():
            value.setdefault(key, default)


def _collect_outputs(outdir: Path, genome: str) -> dict[str, str]:
    """Devuelve solo archivos que existen realmente en ``outdir``."""
    outputs: dict[str, str] = {}
    for label, pattern in _OUTPUT_NAMES.items():
        path = outdir / pattern.format(g=genome)
        if path.is_file():
            outputs[label] = str(path)
    if "operons" in outputs:
        outputs["regulation"] = outputs["operons"]
    known = {Path(value) for value in outputs.values()}
    extra = sorted(str(p) for p in outdir.glob(f"{genome}.*")
                   if p.is_file() and p not in known)
    if extra:
        outputs["extra"] = extra
    return outputs


#: Columnas del TSV TSS/TTS que escribe el adaptador. Son las del registro
#: canónico de ``tss_layer``; no se añade ninguna coordenada calculada.
_TSS_TTS_COLUMNS = (
    "type", "contig", "pos", "end", "strand", "source", "class",
    "location", "mechanism", "detail",
)


def _tsv_cell(value: Any) -> str:
    """Celda TSV en una línea (sin tabs ni saltos embebidos)."""
    if value is None:
        return ""
    return " ".join(str(value).split())


def _write_tss_tts(loaded: dict[str, Any], outdir: Path,
                   genome: str) -> str:
    """Escribe las marcas cargadas; sin posiciones no crea archivo.

    Devuelve la ruta escrita o ``""``. Las filas son exactamente las que
    devolvió ``tss_layer.load``: no se interpolan ni se aproximan posiciones.
    """
    rows = list(loaded.get("tss") or []) + list(loaded.get("tts") or [])
    if not rows:
        return ""
    path = outdir / f"{genome}.tss_tts.tsv"
    lines = ["\t".join(_TSS_TTS_COLUMNS)]
    for record in rows:
        lines.append("\t".join(
            _tsv_cell(record.get(column)) for column in _TSS_TTS_COLUMNS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _tss_attach_payload(features: list[dict[str, Any]],
                        error: str = "") -> dict[str, Any]:
    """Conteos del enlace TSS/TTS a features, sin fabricar marcas."""
    items = [feature for feature in features or []
             if isinstance(feature, dict)]
    return {
        "ok": not error,
        "max_distance": tss_layer.MAX_ATTACH_DISTANCE,
        "n_features": len(items),
        "n_with_tss_upstream": sum(
            1 for feature in items
            if isinstance(feature.get("tss_upstream"), dict)),
        "n_with_tts_downstream": sum(
            1 for feature in items
            if isinstance(feature.get("tts_downstream"), dict)),
        "n_with_tss_inside": sum(
            1 for feature in items
            if int(feature.get("n_tss_inside") or 0) > 0),
        "n_with_tts_inside": sum(
            1 for feature in items
            if int(feature.get("n_tts_inside") or 0) > 0),
        "error": error,
    }


def _tss_tts_payload(loaded: dict[str, Any], features: list[dict[str, Any]],
                     attach_error: str = "") -> dict[str, Any]:
    """Resume la carga TSS/TTS para ``res["tss_tts"]``.

    ``loaded`` es el dict de ``tss_layer.load`` (o el registro de fallo que
    devuelve ``_call``). Se conservan ``status``, conteos, ``sources``,
    ``reason`` y los límites; un genoma sin tabla queda en ``not_available``
    con listas vacías y sin coordenadas inventadas.
    """
    if isinstance(loaded.get("tss"), list) and isinstance(loaded.get("tts"), list):
        payload = tss_layer.summary(loaded)
        payload["available"] = (
            payload.get("status") != tss_layer.STATUS_NOT_AVAILABLE)
    else:
        # La carga lanzó una excepción: no hay filas ni posiciones.
        payload = {
            "genome": str(loaded.get("genome") or ""),
            "status": str(loaded.get("status") or "not_available"),
            "available": False,
            "reason": str(loaded.get("reason") or ""),
            "n_tss": 0,
            "n_tts": 0,
            "n_total": 0,
            "n_tss_unique_positions": 0,
            "n_tts_unique_positions": 0,
            "n_duplicates": 0,
            "n_rejected": 0,
            "n_shared_positions": 0,
            "by_strand": {"+": 0, "-": 0, ".": 0},
            "tss_classes": {},
            "tts_classes": {},
            "tss_locations": {},
            "tts_locations": {},
            "tss_sources": {},
            "tts_sources": {},
            "mechanisms": {},
            "contigs": [],
            "pos_min": None,
            "pos_max": None,
            "sources": list(loaded.get("sources") or []),
            "limitation": tss_layer.LIMITATION,
        }
    payload["ok"] = payload.get("status") in {
        tss_layer.STATUS_AVAILABLE, tss_layer.STATUS_PARTIAL,
    }
    payload["genome"] = str(payload.get("genome") or loaded.get("genome") or "")
    payload["data_root"] = str(loaded.get("data_root") or "")
    payload["searched"] = list(loaded.get("searched") or [])
    payload["limitation"] = (
        payload.get("limitation") or tss_layer.LIMITATION)
    # Las listas normalizadas son las mismas filas que escribe el TSV; en un
    # genoma sin tabla quedan vacías, nunca rellenadas con una aproximación.
    payload["tss"] = list(loaded.get("tss") or [])
    payload["tts"] = list(loaded.get("tts") or [])
    payload["limits"] = {
        "attach_max_distance_bp": tss_layer.MAX_ATTACH_DISTANCE,
        "predicts_positions": False,
        "downloads": False,
        "cross_genome_transfer": False,
        "detail": tss_layer.LIMITATION,
    }
    payload["attach"] = _tss_attach_payload(features, attach_error)
    if attach_error:
        reason = payload.get("reason") or "enlace TSS/TTS falló"
        payload.update(
            ok=False, status="error", available=False,
            reason=f"{reason}; {attach_error}",
            attach_error=attach_error,
        )
    return payload


def run_in_memory(
    filename: str,
    data: bytes,
    outdir: Path,
    data_root: Path,
    options: dict[str, Any] | None = None,
    threads: int = 4,
) -> dict[str, Any]:
    """Ejecuta el pipeline completo sobre un FASTA ya cargado en memoria.

    Parameters
    ----------
    filename:
        Nombre original (se usa para ``name`` y para ``safe_name``/genoma).
    data:
        Contenido FASTA, comprimido con gzip o no.
    outdir:
        Directorio donde se escriben GFF, proteínas, TSV/JSON y el informe.
    data_root:
        Raíz de datos del proyecto (``data/``): Swiss-Prot, Rfam, taxonomía,
        complejos, modelo ORCA, etc. Puede no existir; las capas degradan.
    options:
        Opciones del CLI/GUI. Acepta tanto claves positivas (``ncrna``,
        ``functional``, ...) como los flags ``no_*`` del CLI. Con ``None`` se
        usan los defaults de :data:`DEFAULT_OPTIONS` (metabolismo apagado).
    threads:
        CPUs para las capas que lo aceptan.

    Returns
    -------
    dict
        Resultado de ``annotate`` con ``tss_tts``, ``ncrna``, ``taxonomy``,
        ``functional``, ``domains``, ``context_c``, ``structural_b``,
        ``structural_b_filters``, ``complexes``, ``regulation``, ``membrane``,
        ``oric``, ``metabolic`` (``None`` si no se pidió) y ``report_md``,
        además de ``outputs`` y ``layer_status``. Cada capa expone ``ok``,
        ``status`` y ``reason``. ``tss_tts`` es un resumen de la carga local
        (``status``, conteos, ``sources``, ``reason``, ``limitation``,
        ``limits`` y ``attach``); las filas se escriben en
        ``<genoma>.tss_tts.tsv`` solo si existen.
    """
    started = time.perf_counter()
    outdir = Path(outdir)
    data_root = Path(data_root)
    threads = max(1, int(threads or 1))
    opts, degraded = _normalize_options(options)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = _as_fasta_bytes(data)
    annotate_options = {
        "cds": opts["cds"],
        "trna": opts["trna"],
        "rrna": opts["rrna"],
        "ghost": opts["ghost"],
        "deserts": opts["deserts"],
        "esm": opts["esm"],
        "esm_top": opts["esm_top"],
        "min_region_len": opts["min_region_len"],
        "min_aa": opts["min_aa"],
    }
    try:
        res = annotate(str(filename), raw, annotate_options)
    except ValueError:
        # FASTA vacío/inválido: es un error de entrada, no una dependencia.
        raise
    except Exception as exc:  # noqa: BLE001 - herramienta/base ausente o rota
        return _failed_annotation(filename, outdir, data_root, threads, opts,
                                  degraded, exc)
    res["threads"] = threads
    res["outdir"] = str(outdir)
    res["data_root"] = str(data_root)
    res["options"] = dict(opts)
    res["degraded"] = degraded

    # TSS/TTS precalculados: solo lectura local. No se descarga nada ni se
    # predice una posición; un genoma sin tabla queda en ``not_available``.
    tss_loaded = _call(tss_layer.load, res["genome"], data_root)
    attach_error = ""
    try:
        tss_layer.attach(
            res.get("features") or [],
            tss_loaded.get("tss") or [],
            tss_loaded.get("tts") or [],
        )
    except Exception as exc:  # noqa: BLE001 - el enlace no tumba el pipeline
        attach_error = f"{type(exc).__name__}: {exc}"
    res["tss_tts"] = _tss_tts_payload(
        tss_loaded, res.get("features") or [], attach_error)
    try:
        _write_tss_tts(tss_loaded, outdir, res["genome"])
    except OSError as exc:  # pragma: no cover - disco/permiso
        res["tss_tts"]["write_error"] = f"{type(exc).__name__}: {exc}"

    # 1) ncRNA. ``run_ncrna`` muta ``features``; si encontró algo se rehace el
    #    GFF para que el archivo final incluya esos features (igual que CLI).
    if opts["ncrna"]:
        nc = _call(run_ncrna, res.get("sequences") or {},
                   res.get("features") or [], outdir, res["genome"], data_root,
                   cpu=threads)
        if nc.get("ok") and nc.get("n"):
            try:
                res["gff"] = build_gff(res["sequences"], res["features"])
            except Exception as exc:  # noqa: BLE001
                nc["gff_error"] = f"{type(exc).__name__}: {exc}"
    else:
        nc = _disabled("ncRNA desactivado por opción")
    res["ncrna"] = nc

    # 2) Taxonomía local precalculada.
    if opts["taxonomy"]:
        tax = _call(run_taxonomy, res["genome"], data_root, outdir)
    else:
        tax = _disabled("taxonomía desactivada por opción")
    res["taxonomy"] = tax

    # 3) Anotación funcional (DIAMOND/Swiss-Prot). Si no hay proteínas no se
    #    ejecuta y las capas derivadas quedan explícitamente sin evaluar.
    functional_ready = False
    if not opts["functional"]:
        func = _disabled("anotación funcional desactivada por opción")
    elif not res.get("proteins"):
        func = _unavailable("sin proteínas: CDS desactivado o sin ORFs")
    else:
        func = _call(annotate_functional, res["proteins"], data_root, outdir,
                     res["genome"], threads=threads)
        functional_ready = True
    res["functional"] = func
    rows = list(func.get("rows") or [])
    hits = func.get("hits") or {}

    # 4) Dominios C.
    if not opts["domains"]:
        domains = _disabled("dominios desactivados por opción")
    elif not functional_ready:
        domains = _unavailable("requiere anotación funcional")
    else:
        domains = _call(run_domains, rows, res.get("proteins") or [], outdir,
                        res["genome"], data_root, cpu=threads)
    res["domains"] = domains

    # 5) Contexto C (hidropatía/señal/lipobox y TMHMM opcional).
    if not opts["context"]:
        context = _disabled("contexto C desactivado por opción")
    elif not functional_ready:
        context = _unavailable("requiere anotación funcional")
    else:
        context = _call(run_context, rows, res.get("proteins") or [], outdir,
                        res["genome"], cpu=threads)
    res["context_c"] = context

    # 6) Estado estructural B y filtros F1-F5.
    if not opts["structural_b"]:
        structural_b = _disabled("estado estructural B desactivado por opción")
    elif not functional_ready:
        structural_b = _unavailable("requiere anotación funcional")
    else:
        structural_b = _call(run_structural_b, rows, outdir, res["genome"],
                             data_root)
    res["structural_b"] = structural_b

    if not opts["structural_b_filters"]:
        filters = _disabled("filtros estructurales desactivados por opción")
    elif not functional_ready:
        filters = _unavailable("requiere anotación funcional")
    else:
        filters = _call(
            evaluate_structural_b, rows, outdir, res["genome"], data_root,
            foldseek_db=opts.get("foldseek_local_db"), cpu=threads,
            external=opts["structure_b_external"],
        )
    res["structural_b_filters"] = filters

    # 7) Complejos: firmas curadas + modelo K-12 + transferencia por ortología.
    if not opts["complexes"]:
        complexes = _disabled("complejos desactivados por opción")
    elif not functional_ready:
        complexes = _unavailable("requiere anotación funcional")
    else:
        try:
            comp = infer_complexes(res.get("features") or [], hits,
                                   res["genome"])
            model = _call(run_complex_model, res.get("features") or [], hits,
                          outdir, res["genome"], data_root)
            if model.get("ok"):
                comp["complexes"] += model.get("complexes") or []
                comp["by_protein"].update(model.get("by_protein") or {})
                comp["n_candidate"] = model.get("n_candidate", 0)
            transfer = _call(load_complex_transfer, data_root, res["genome"],
                             res.get("features") or [])
            if transfer.get("ok"):
                comp["complexes"] += transfer.get("complexes") or []
                comp["by_protein"].update(transfer.get("by_protein") or {})
                comp["n_candidate"] = (comp.get("n_candidate", 0)
                                       + transfer.get("n_candidate", 0))
            enrich_rows(rows, comp)
            write_paths(rows, outdir, res["genome"])
            comp["protein_rows"] = rows
            comp["model"] = {
                "ok": bool(model.get("ok")),
                "status": model.get("status", "ok" if model.get("ok")
                                    else "not_available"),
                "reason": model.get("reason", ""),
            }
            comp["transfer"] = {
                "ok": bool(transfer.get("ok")),
                "status": transfer.get("status", "ok" if transfer.get("ok")
                                       else "not_available"),
                "reason": transfer.get("reason", ""),
            }
            write_complexes(comp, outdir)
            comp.update(ok=True, status="ok", reason="")
            complexes = comp
        except Exception as exc:  # noqa: BLE001
            complexes = _failure(exc)
    res["complexes"] = complexes

    # 8) Regulación putativa (operones/promotores por reglas).
    if not opts["regulation"]:
        regulation = _disabled("regulación desactivada por opción")
    elif not res.get("features"):
        regulation = _unavailable("sin features CDS para construir operones")
    else:
        regulation = _call(run_regulation, res["features"], outdir,
                           res["genome"], max_gap=opts["operon_max_gap"],
                           sequences=res.get("sequences"))
    res["regulation"] = regulation

    # 9) Membrana/transporte (hmmscan + panel Pfam).
    if not opts["membrane"]:
        membrane = _disabled("membrana desactivada por opción")
    elif not res.get("gff"):
        membrane = _unavailable("sin GFF para evaluar membrana")
    else:
        membrane = _call(run_membrane, res, outdir, data_root, cpu=threads)
    res["membrane"] = membrane

    # 10) oriC (ORCA + GC skew; DoriC solo si se pasa accession).
    if not opts["oric"]:
        oric_result = _disabled("oriC desactivado por opción")
    elif not res.get("sequences"):
        oric_result = _unavailable("sin secuencias de contigs")
    else:
        try:
            oric_result = _normalize(
                oric.run(res, data_root, accession=opts.get("accession", "")))
            oric.write(res, oric_result, outdir)
        except Exception as exc:  # noqa: BLE001
            oric_result = _failure(exc)
    res["oric"] = oric_result

    # 11) Metabolismo opcional.
    metabolic = None
    if not opts["metabolic"]:
        res["metabolic"] = None
        res["metabolic_status"] = _disabled(
            "metabolismo desactivado (capa opcional)")
    elif not res.get("proteins"):
        res["metabolic"] = None
        res["metabolic_status"] = _unavailable(
            "sin proteínas: no se puede construir la red metabólica")
    else:
        try:
            panel = Path(opts["panel_hmm"]) if opts.get("panel_hmm") else (
                data_root / "db" / "pfam" / "panel_focus.hmm")
            tdb = Path(opts["transport_db"]) if opts.get("transport_db") else (
                data_root / "db" / "transport.dmnd")
            pdb = Path(opts["producer_db"]) if opts.get("producer_db") else (
                data_root / "db" / "producer_family.dmnd")
            metabolic = build_network(
                res["proteins"], data_root, opts["reference"] or None,
                threads=threads,
                panel_hmm=panel if panel.exists() else None,
                transport_db=tdb if tdb.exists() else None,
                producer_db=pdb if pdb.exists() else None,
            )
            comp = res.get("complexes")
            if isinstance(comp, dict) and comp.get("ok"):
                metabolic["gaps"] = enrich_gaps(
                    metabolic.get("gaps", []), comp, data_root)
            res["metabolic"] = metabolic
            res["metabolic_status"] = {"ok": True, "status": "ok", "reason": ""}
        except Exception as exc:  # noqa: BLE001
            metabolic = None
            res["metabolic"] = None
            res["metabolic_status"] = _failure(exc)

    # 12) Informe final y artefactos, con las funciones reales del CLI.
    try:
        res["report_md"] = build_report(res)
    except Exception as exc:  # noqa: BLE001
        res["report_error"] = f"{type(exc).__name__}: {exc}"
        res.setdefault("report_md", "")

    _ensure_shapes(res)
    try:
        write_outputs(res, outdir, metabolic)
    except Exception as exc:  # noqa: BLE001
        res["write_error"] = f"{type(exc).__name__}: {exc}"

    res["layer_status"] = {
        layer: (res.get(layer) or {}).get("status", "not_available")
        for layer in _LAYER_SHAPES
    }
    res["layer_status"]["metabolic"] = (
        res.get("metabolic_status") or {}).get("status", "not_available")
    res["outputs"] = _collect_outputs(outdir, res["genome"])
    res["elapsed"] = round(time.perf_counter() - started, 3)
    res["ok"] = not (res.get("write_error") or res.get("report_error"))
    res["status"] = "ok" if res["ok"] else "error"
    return res


__all__ = ["DEFAULT_OPTIONS", "run_in_memory"]
