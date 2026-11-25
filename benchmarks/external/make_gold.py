#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construye un gold pequeño e independiente desde un GFF3/GBFF local.

El gold no se descarga y no se deriva de Bakta/Prokka/DFAST: se extrae una
ventana de un ensamblado y su anotación de referencia local (por ejemplo NCBI
RefSeq). La ventana se reindexa a 1 y se escribe junto a un manifiesto con
SHA256, conteos y un comando ``--run`` reproducible por genoma.

Uso mínimo::

    .venv/bin/python benchmarks/external/make_gold.py \
      --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
      --annotation data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
      --out-dir benchmarks/external/gold/chlorobaculum_60kb \
      --name chlorobaculum_60kb

Salidas dentro de ``--out-dir``:

- ``genome.fna``: ventana FASTA (sin ``.gz``) con coordenadas reindexadas.
- ``reference.gff3``: CDS/tRNA/rRNA de la ventana, 1-based inclusivo.
- ``gold_manifest.json``: procedencia, SHA256, selección y comando ``--run``.
- ``run_benchmark.sh``: ejecuta el benchmark con ``--run`` para este genoma.

El gold hereda los criterios y errores de la anotación de origen; no es un gold
curado a mano ni una verdad biológica absoluta. No descarga imágenes ni bases.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import quote

MODULE_DIR = Path(__file__).resolve().parent
ROOT = MODULE_DIR.parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import external_benchmark as ext  # noqa: E402

GOLD_SCHEMA = ext.GOLD_SCHEMA
PROFILE_SCHEMA = ext.PROFILE_SCHEMA
GOLD_LIMITATIONS = [
    "El gold es una ventana extraída de una anotación de referencia local; hereda sus criterios y errores.",
    "No es un gold curado a mano ni una verdad biológica independiente de toda base de datos.",
    "Solo incluye CDS/tRNA/rRNA; no evalúa ncRNA, tmRNA, pseudogenes, dominios ni rutas.",
    "Las coordenadas del gold están reindexadas a 1 dentro de la ventana; no coinciden con el ensamblado completo.",
    "Una herramienta ausente no es un cero: el runner la deja not_available sin métricas.",
]


# ---------------------------------------------------------------------------
# Entrada de secuencias
# ---------------------------------------------------------------------------
def _utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def read_fasta_sequences(path: str | Path) -> dict[str, dict[str, str]]:
    """Lee un FASTA plano/gzip conservando la secuencia (no solo la longitud)."""
    path = Path(path)
    if not path.is_file():
        raise ext.InputError(f"no existe el FASTA: {path}")
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    chunks: list[str] = []
    with ext.open_text(path) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    current["sequence"] = "".join(chunks)
                    records[current["id"]] = current
                header = line[1:].strip()
                identifier = header.split()[0] if header else f"contig_{len(records) + 1}"
                if identifier in records:
                    raise ext.InputError(
                        f"ID FASTA duplicado en {path}: {identifier} (línea {line_number})"
                    )
                current = {"id": identifier, "description": header, "sequence": ""}
                chunks = []
                continue
            if not line.strip():
                continue
            if current is None:
                raise ext.InputError(
                    f"FASTA inválido: secuencia antes de la cabecera (línea {line_number})"
                )
            chunks.append(line.strip())
    if current is not None:
        current["sequence"] = "".join(chunks)
        records[current["id"]] = current
    if not records:
        raise ext.InputError(f"el FASTA no contiene contigs: {path}")
    return records


def read_genbank_sequences(path: str | Path) -> dict[str, dict[str, str]]:
    """Extrae las secuencias de un GBFF con Biopython (dependencia opcional)."""
    try:
        from Bio import SeqIO  # noqa: PLC0415 - dependencia opcional
    except ImportError as exc:  # pragma: no cover - depende del entorno
        raise ext.DependencyError(
            "Biopython es necesario para extraer secuencias de GBFF; "
            "instalar con: pip install biopython"
        ) from exc
    records: dict[str, dict[str, str]] = {}
    with ext.open_text(path) as fh:
        for record in SeqIO.parse(fh, "genbank"):
            identifier = str(record.id or record.name or "").strip()
            if not identifier:
                continue
            records[identifier] = {
                "id": identifier,
                "description": str(record.description or identifier),
                "sequence": str(record.seq),
            }
    if not records:
        raise ext.InputError(f"el GBFF no contiene secuencias: {path}")
    return records


def load_sequences(
    annotation_path: Path,
    annotation_format: str,
    genome_path: Path | None,
) -> tuple[dict[str, dict[str, str]], str]:
    if genome_path is not None:
        return read_fasta_sequences(genome_path), str(genome_path)
    if annotation_format == "genbank":
        return read_genbank_sequences(annotation_path), str(annotation_path)
    raise ext.InputError(
        "para una anotación GFF3 hay que indicar --genome con el FASTA del ensamblado"
    )


# ---------------------------------------------------------------------------
# Selección de ventana
# ---------------------------------------------------------------------------
def _safe_name(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "gold"


def _gff_value(value: str) -> str:
    """Escapa un valor de atributo GFF3 sin romper ``clave=valor;...``."""
    return quote(str(value), safe=" /:+-_.()[]'")


def choose_window(
    features: Sequence[ext.Feature],
    contig_length: int,
    *,
    max_bp: int,
    flank: int,
    min_features: int,
) -> dict[str, Any]:
    """Elige la ventana con más features completos que quepa en ``max_bp``.

    Se ancla en el inicio de un CDS (o feature si no hay CDS) y se recorre la
    lista ordenada. El resultado es determinista para el mismo GFF/GBFF.
    """
    ordered = sorted(features, key=lambda item: (item.start, item.end, item.type))
    anchors = [item for item in ordered if item.type == "CDS"] or ordered
    if not anchors:
        raise ext.InputError("no hay features CDS/tRNA/rRNA para construir el gold")

    best: dict[str, Any] | None = None
    for anchor in anchors:
        start = max(1, int(anchor.start) - flank)
        start = min(start, max(1, contig_length - max_bp + 1))
        end = min(contig_length, start + max_bp - 1)
        selected = [
            item for item in ordered if item.start >= start and item.end <= end
        ]
        if not selected:
            continue
        score = (len(selected), -start)
        if best is None or score > best["score"]:
            best = {
                "score": score,
                "start": start,
                "end": end,
                "features": selected,
            }
        if len(selected) >= min_features and anchor is anchors[0]:
            # El primer ancla ya cumple el mínimo: no hace falta seguir.
            break

    if best is None:
        raise ext.InputError("no se pudo seleccionar una ventana con features completos")
    if len(best["features"]) < min_features:
        best["warning"] = (
            f"la ventana solo contiene {len(best['features'])} features "
            f"(mínimo solicitado {min_features}); se conserva la mejor disponible"
        )
    return best


# ---------------------------------------------------------------------------
# Escritura del gold
# ---------------------------------------------------------------------------
def _write_fasta(path: Path, seqid: str, sequence: str, header: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{seqid} {header}".strip()]
    for index in range(0, len(sequence), 80):
        lines.append(sequence[index : index + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_gff(
    path: Path,
    *,
    seqid: str,
    contig_length: int,
    features: Sequence[ext.Feature],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "##gff-version 3",
        f"##sequence-region {seqid} 1 {contig_length}",
        f"# gold extraído de una anotación local; ventana reindexada a 1",
        f"# generado por benchmarks/external/make_gold.py",
    ]
    seen_ids: Counter[str] = Counter()
    for index, feature in enumerate(
        sorted(features, key=lambda item: (item.start, item.end, item.type)), 1
    ):
        attrs = dict(feature.attrs or {})
        identifier = attrs.get("ID", "").strip()
        if not identifier:
            identifier = f"gold_{index:05d}"
        seen_ids[identifier] += 1
        if seen_ids[identifier] > 1:
            identifier = f"{identifier}_{seen_ids[identifier]}"
        attrs["ID"] = identifier
        attr_text = ";".join(
            f"{_gff_value(key)}={_gff_value(value)}"
            for key, value in attrs.items()
            if str(value).strip()
        )
        lines.append(
            "\t".join(
                [
                    seqid,
                    "gold_extract",
                    feature.type,
                    str(feature.start),
                    str(feature.end),
                    ".",
                    feature.strand if feature.strand in {"+", "-"} else ".",
                    ".",
                    attr_text,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _relative(path: Path) -> str:
    # abspath normaliza sin seguir symlinks (.venv/bin/python -> python3).
    candidate = Path(os.path.abspath(str(path)))
    root = Path(os.path.abspath(str(ROOT)))
    try:
        return candidate.relative_to(root).as_posix()
    except ValueError:
        return str(candidate)


def _python_executable() -> str:
    venv_python = ROOT / ".venv" / "bin" / "python"
    if venv_python.is_file() and os.access(venv_python, os.X_OK):
        return _relative(venv_python)
    return sys.executable or "python3"


def _run_command(
    name: str,
    genome_path: Path,
    reference_path: Path,
    *,
    tools: str,
    timeout: float,
    threads: int,
) -> list[str]:
    report_dir = Path("reports") / "external_benchmark" / name
    return [
        _python_executable(),
        "benchmarks/external/external_benchmark.py",
        "--genome",
        _relative(genome_path),
        "--reference",
        _relative(reference_path),
        "--tools",
        tools,
        "--run",
        "--threads",
        str(threads),
        "--timeout",
        str(timeout),
        "--save-gffs",
        (report_dir / "evidence").as_posix(),
        "--json",
        (report_dir / "benchmark.json").as_posix(),
        "--md",
        (report_dir / "benchmark.md").as_posix(),
    ]


def _write_run_script(path: Path, command: Sequence[str]) -> None:
    if not command:
        raise ext.InputError("comando --run vacío")
    rest = " ".join(_shell_quote(item) for item in command[1:])
    root_fallback = _shell_quote(str(Path(os.path.abspath(str(ROOT)))))
    text = f"""#!/usr/bin/env bash
# Ejecución reproducible con --run para este gold.
# No descarga herramientas, imágenes ni bases. Bakta DB: exportar BAKTA_DB.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
if [[ -n "${{PROJECT_ROOT:-}}" ]]; then
    ROOT="$PROJECT_ROOT"
elif ROOT_GIT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    ROOT="$ROOT_GIT"
else
    ROOT="$(cd "$SCRIPT_DIR/../../../.." 2>/dev/null && pwd || true)"
fi
if [[ -z "${{ROOT:-}}" || ! -f "$ROOT/benchmarks/external/external_benchmark.py" ]]; then
    ROOT={root_fallback}
fi
cd "$ROOT"
PYTHON="${{PYTHON:-$ROOT/.venv/bin/python}}"
if [[ ! -x "$PYTHON" ]]; then
    PYTHON="${{PYTHON_FALLBACK:-python3}}"
fi
exec "$PYTHON" {rest} "$@"
"""
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)


def _shell_quote(value: str) -> str:
    text = str(value)
    if not text:
        return "''"
    if re.match(r"^[A-Za-z0-9_./:=+-]+$", text):
        return text
    return "'" + text.replace("'", "'\\''") + "'"


def build_gold(
    *,
    annotation_path: str | Path,
    out_dir: str | Path,
    genome_path: str | Path | None = None,
    name: str | None = None,
    max_bp: int = 60_000,
    flank: int = 200,
    min_features: int = 3,
    contig: str | None = None,
    exclude_pseudo: bool = True,
    force: bool = False,
    tools: str = "all",
    timeout: float = 3600.0,
    threads: int = 1,
) -> dict[str, Any]:
    """Construye el gold, escribe FASTA/GFF/manifiesto/script y devuelve el manifiesto."""
    if max_bp < 100:
        raise ext.InputError("--max-bp debe ser >= 100")
    if flank < 0:
        raise ext.InputError("--flank no puede ser negativo")
    if min_features < 1:
        raise ext.InputError("--min-features debe ser >= 1")
    annotation_path = Path(annotation_path)
    out_dir = Path(out_dir)
    if not annotation_path.is_file():
        raise ext.InputError(f"no existe la anotación: {annotation_path}")
    genome_path = Path(genome_path) if genome_path is not None else None
    if genome_path is not None and not genome_path.is_file():
        raise ext.InputError(f"no existe el genoma FASTA: {genome_path}")
    if out_dir.exists() and not out_dir.is_dir():
        raise ext.InputError(f"la ruta de salida existe y no es un directorio: {out_dir}")
    if out_dir.exists() and any(out_dir.iterdir()) and not force:
        raise ext.InputError(
            f"el directorio de salida no está vacío: {out_dir} (usar --force)"
        )

    annotation_format = ext.detect_format(annotation_path)
    annotation = ext.load_annotations(
        annotation_path, exclude_pseudo=exclude_pseudo
    )
    if not annotation.features:
        raise ext.InputError("la anotación no tiene CDS/tRNA/rRNA evaluables")

    sequences, sequence_source = load_sequences(
        annotation_path, annotation_format, genome_path
    )
    contigs = [
        {"id": identifier, "description": record["description"], "length": len(record["sequence"])}
        for identifier, record in sequences.items()
    ]
    genome_info = ext.GenomeInfo(
        path=str(genome_path or annotation_path),
        sha256=ext.sha256_file(genome_path or annotation_path),
        contigs=contigs,
    )
    resolver = ext.build_seqid_resolver(genome_info)

    mapped: list[ext.Feature] = []
    skipped: list[dict[str, Any]] = []
    for feature in annotation.features:
        seqid = resolver.resolve(feature.seqid)
        record = sequences.get(seqid)
        if record is None:
            skipped.append(
                {"seqid": feature.seqid, "reason": "seqid fuera del FASTA"}
            )
            continue
        length = len(record["sequence"])
        if feature.end > length or feature.start < 1:
            skipped.append(
                {
                    "seqid": feature.seqid,
                    "start": feature.start,
                    "end": feature.end,
                    "reason": f"coordenadas fuera del contig ({length} bp)",
                }
            )
            continue
        mapped.append(feature.with_seqid(seqid))

    if not mapped:
        raise ext.InputError("ningún feature de la anotación cae dentro del FASTA")

    counts_by_contig = Counter(feature.seqid for feature in mapped)
    if contig is not None:
        target = resolver.resolve(contig)
        if target not in sequences:
            raise ext.InputError(f"el contig solicitado no existe en el FASTA: {contig}")
        if counts_by_contig.get(target, 0) == 0:
            raise ext.InputError(f"el contig solicitado no tiene features: {contig}")
    else:
        target = sorted(
            counts_by_contig, key=lambda item: (-counts_by_contig[item], item)
        )[0]
    contig_features = [feature for feature in mapped if feature.seqid == target]
    contig_sequence = sequences[target]["sequence"]
    window = choose_window(
        contig_features,
        len(contig_sequence),
        max_bp=max_bp,
        flank=flank,
        min_features=min_features,
    )
    window_start = int(window["start"])
    window_end = int(window["end"])
    window_sequence = contig_sequence[window_start - 1 : window_end]
    shifted = [
        ext.Feature(
            seqid=target,
            type=feature.type,
            start=feature.start - window_start + 1,
            end=feature.end - window_start + 1,
            strand=feature.strand,
            attrs=feature.attrs,
            source=feature.source,
        )
        for feature in window["features"]
    ]
    shifted_counts = {
        category: sum(1 for feature in shifted if feature.type == category)
        for category in ext.CATEGORIES
    }
    if shifted_counts["CDS"] == 0:
        raise ext.InputError("la ventana seleccionada no contiene ningún CDS")

    out_dir.mkdir(parents=True, exist_ok=True)
    gold_fasta = out_dir / "genome.fna"
    gold_gff = out_dir / "reference.gff3"
    manifest_path = out_dir / "gold_manifest.json"
    run_script = out_dir / "run_benchmark.sh"
    gold_name = _safe_name(name or out_dir.name)
    header = (
        f"gold_name={gold_name} source={target} "
        f"window={window_start}-{window_end} of={len(contig_sequence)}bp"
    )
    _write_fasta(gold_fasta, target, window_sequence, header)
    _write_gff(
        gold_gff,
        seqid=target,
        contig_length=len(window_sequence),
        features=shifted,
    )

    command = _run_command(
        gold_name,
        gold_fasta,
        gold_gff,
        tools=tools,
        timeout=timeout,
        threads=threads,
    )
    _write_run_script(run_script, command)

    warnings: list[str] = []
    if window.get("warning"):
        warnings.append(str(window["warning"]))
    if skipped:
        warnings.append(f"{len(skipped)} feature(s) omitidos por seqid o coordenadas")
    if annotation.warnings:
        warnings.extend(annotation.warnings)
    if sequence_source == str(annotation_path) and annotation_format == "genbank":
        warnings.append(
            "la secuencia se extrajo del propio GBFF; no se reutilizó un FASTA externo"
        )

    manifest: dict[str, Any] = {
        "schema": GOLD_SCHEMA,
        "profile_schema": PROFILE_SCHEMA,
        "name": gold_name,
        "generated_at": _utc_now(),
        "source": {
            "annotation": {
                "path": ext.display_path(annotation_path),
                "format": annotation_format,
                "file": ext.file_fingerprint(annotation_path),
                "features_total": len(annotation.features),
                "excluded_pseudo": annotation.excluded_pseudo,
            },
            "genome": {
                "path": ext.display_path(genome_path) if genome_path else None,
                "sequence_source": ext.display_path(sequence_source),
                "file": ext.file_fingerprint(genome_path)
                if genome_path
                else ext.file_fingerprint(annotation_path),
                "sha256": genome_info.sha256,
                "contigs": len(contigs),
            },
        },
        "selection": {
            "seqid": target,
            "source_contig_length": len(contig_sequence),
            "window_start": window_start,
            "window_end": window_end,
            "window_bp": len(window_sequence),
            "flank": flank,
            "max_bp": max_bp,
            "contig_features_total": len(contig_features),
            "features_selected": len(shifted),
            "counts": shifted_counts,
            "skipped": skipped[:100],
            "skipped_total": len(skipped),
        },
        "outputs": {
            "genome": {
                **ext.file_fingerprint(gold_fasta),
                "sha256": ext.sha256_file(gold_fasta),
            },
            "reference": {
                **ext.file_fingerprint(gold_gff),
                "sha256": ext.sha256_file(gold_gff),
            },
            "run_script": ext.display_path(run_script),
            "manifest": ext.display_path(manifest_path),
        },
        "counts": shifted_counts,
        "warnings": warnings,
        "method": (
            "Ventana de una anotación local reindexada a 1; features CDS/tRNA/rRNA "
            "completamente contenidos en la ventana. No se deriva de ninguna "
            "herramienta evaluada."
        ),
        "limitations": list(GOLD_LIMITATIONS),
        "run": {
            "schema": PROFILE_SCHEMA,
            "tools": (
                list(ext.TOOLS)
                if tools.strip().lower() == "all"
                else [item.strip() for item in tools.split(",") if item.strip()]
            ),
            "run": True,
            "threads": threads,
            "timeout_seconds": timeout,
            "command": command,
            "shell": " ".join(_shell_quote(item) for item in command),
            "script": ext.display_path(run_script),
            "working_directory": "raíz del repositorio",
            "reports": {
                "json": f"reports/external_benchmark/{gold_name}/benchmark.json",
                "md": f"reports/external_benchmark/{gold_name}/benchmark.md",
                "evidence": f"reports/external_benchmark/{gold_name}/evidence",
            },
            "bakta_db": (
                "Bakta requiere una DB local; exportar BAKTA_DB o añadir "
                "--bakta-db. No se descarga."
            ),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Construye un gold pequeño desde un GFF3/GBFF local y un comando "
            "--run reproducible. No descarga datos ni ejecuta herramientas."
        )
    )
    parser.add_argument("--genome", type=Path, default=None, help="FASTA local (necesario para GFF3)")
    parser.add_argument(
        "--annotation",
        "--reference",
        dest="annotation",
        type=Path,
        required=True,
        help="GFF3 o GBFF local de referencia",
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="directorio del gold")
    parser.add_argument("--name", default=None, help="nombre lógico (por defecto el directorio)")
    parser.add_argument("--max-bp", type=int, default=60_000, help="tamaño máximo de la ventana")
    parser.add_argument("--flank", type=int, default=200, help="bp de flanco antes del primer CDS")
    parser.add_argument(
        "--min-features", type=int, default=3, help="mínimo deseado de features en la ventana"
    )
    parser.add_argument("--contig", default=None, help="seqid a extraer (por defecto el más denso)")
    parser.add_argument(
        "--keep-pseudo", action="store_true", help="no excluir CDS pseudo"
    )
    parser.add_argument("--force", action="store_true", help="sobrescribir un directorio no vacío")
    parser.add_argument("--tools", default="all", help="herramientas del comando --run")
    parser.add_argument("--threads", type=int, default=1, help="hilos del comando --run")
    parser.add_argument("--timeout", type=float, default=3600.0, help="timeout del comando --run")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads debe ser >= 1")
    if args.timeout <= 0:
        parser.error("--timeout debe ser > 0")
    try:
        manifest = build_gold(
            annotation_path=args.annotation,
            out_dir=args.out_dir,
            genome_path=args.genome,
            name=args.name,
            max_bp=args.max_bp,
            flank=args.flank,
            min_features=args.min_features,
            contig=args.contig,
            exclude_pseudo=not args.keep_pseudo,
            force=args.force,
            tools=args.tools,
            timeout=args.timeout,
            threads=args.threads,
        )
    except ext.BenchmarkError as exc:
        print(f"[make-gold] error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"[make-gold] error de E/S: {exc}", file=sys.stderr)
        return 2
    print(
        f"[make-gold] {manifest['name']}: "
        f"{manifest['selection']['seqid']}:"
        f"{manifest['selection']['window_start']}-{manifest['selection']['window_end']} "
        f"({manifest['selection']['window_bp']} bp, "
        f"CDS={manifest['counts']['CDS']} tRNA={manifest['counts']['tRNA']} "
        f"rRNA={manifest['counts']['rRNA']})"
    )
    print(f"[make-gold] gold: {ext.display_path(args.out_dir)}")
    print("[make-gold] comando --run:")
    print("  " + " ".join(_shell_quote(item) for item in manifest["run"]["command"]))
    print(f"[make-gold] o ejecutar: bash {manifest['run']['script']}")
    return 0


__all__ = [
    "GOLD_SCHEMA",
    "GOLD_LIMITATIONS",
    "build_gold",
    "build_parser",
    "choose_window",
    "load_sequences",
    "main",
    "read_fasta_sequences",
    "read_genbank_sequences",
]


if __name__ == "__main__":
    raise SystemExit(main())
