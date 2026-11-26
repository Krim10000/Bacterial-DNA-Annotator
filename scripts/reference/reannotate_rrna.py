#!/usr/bin/env python3
"""Reanota solo la capa rRNA de los genomas de referencia ya anotados.

Contexto: el GFF del sistema se generó con ``barrnap`` 0.9. En esta máquina su
base ``bac.hmm`` no incluye el modelo 23S y su E-value por defecto (1e-6)
descarta copias de 5S. La referencia RefSeq/NCBI anota 16S/23S/5S con
Infernal/Rfam. Este script vuelve a llamar **solo rRNA** sobre el FASTA del
genoma con ``webapp.rrna_calling`` (modelos Rfam locales) y reescribe las filas
``rRNA`` del GFF del sistema y el conteo ``features.rRNA`` del resumen.

No lee ni usa el GFF de referencia para llamar features: la referencia solo se
consulta después, para el informe de antes/después. Por defecto es un dry-run;
con ``--apply`` modifica ``results/`` (nunca ``data/reference_genomes/``).

Salidas: ``reports/rrna_reannotation.tsv`` y ``.md``.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "reference"))

from webapp.rrna_calling import (  # noqa: E402
    DEFAULT_RRNA_CM,
    call_rrna_fasta,
    format_gff_line,
    parse_barrnap_gff,
    rrna_class,
)

import coordinate_parity as cp  # noqa: E402
import reference_parity as rp  # noqa: E402

ELIGIBLE = ("annotated", "available")
CLASSES_BY_TYPE = ("16S", "23S", "5S")


# ---------------------------------------------------------------------------
# Entradas / GFF
# ---------------------------------------------------------------------------
def read_manifest(path: Path | str) -> list[dict[str, Any]]:
    return rp.load_manifest(path)


def _fasta_text(path: Path) -> str:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8", errors="replace")


_RFAM_CLASSES = {
    "RF00177": "16S",
    "RF02541": "23S",
    "RF00001": "5S",
}


def _raw_rrna_rows(gff_path: Path) -> list[dict[str, Any]]:
    """Lee las filas rRNA de un GFF conservando producto y acceso Rfam.

    ``coordinate_parity.parse_gff`` colapsa por ``ID`` y no conserva el
    producto; para diagnosticar 16S/23S/5S hace falta el atributo original.
    """
    opener = rp.open_text
    rows: list[dict[str, Any]] = []
    with opener(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "rRNA":
                continue
            attrs = fields[8]
            product = ""
            match = re.search(r"(?:^|;)product=([^;]+)", attrs)
            if match:
                product = match.group(1).strip()
            rfam = ""
            match = re.search(r"RF\d{5}", attrs)
            if match:
                rfam = match.group(0)
            try:
                start, end = int(fields[3]), int(fields[4])
            except (TypeError, ValueError):
                continue
            rows.append(
                {
                    "contig": fields[0],
                    "start": min(start, end),
                    "end": max(start, end),
                    "strand": fields[6] if fields[6] in ("+", "-") else ".",
                    "product": product,
                    "rfam": rfam,
                    "model": "",
                    "source": "gff",
                }
            )
    return rows


def _old_rrna_features(gff_path: Path) -> list[dict[str, Any]]:
    return _raw_rrna_rows(gff_path)


def _count_by_type(features: Sequence[dict[str, Any]]) -> dict[str, int]:
    counts = {kind: 0 for kind in CLASSES_BY_TYPE}
    for feature in features:
        kind = rrna_class(feature)
        if kind in counts:
            counts[kind] += 1
        else:
            counts.setdefault("other", 0)
            counts["other"] += 1
    return counts


def _reference_rrna_counts(record: dict[str, Any], root: Path) -> dict[str, int]:
    raw = str(record.get("annotation_path") or "").strip()
    path = rp._resolve(root, raw)
    if path is None or not path.is_file():
        return {kind: 0 for kind in CLASSES_BY_TYPE}
    return _count_by_type(_raw_rrna_rows(path))


def _rewrite_gff(
    gff_path: Path,
    new_features: Sequence[dict[str, Any]],
) -> None:
    """Sustituye solo las filas ``rRNA`` conservando el resto del GFF."""
    pre: list[str] = []
    body: list[str] = []
    post: list[str] = []
    in_fasta = False
    for line in gff_path.read_text(encoding="utf-8", errors="replace").splitlines(
        keepends=True
    ):
        if in_fasta:
            post.append(line)
            continue
        if line.startswith("##FASTA"):
            in_fasta = True
            post.append(line)
            continue
        if line.startswith("#") or not line.strip():
            pre.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3 and fields[2] == "rRNA":
            continue
        body.append(line if line.endswith("\n") else line + "\n")
    for feature in new_features:
        body.append(format_gff_line(feature))
    gff_path.write_text("".join(pre + body + post), encoding="utf-8")


def _update_summary(summary_path: Path, n_rrna: int) -> None:
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    features = data.get("features")
    if not isinstance(features, dict):
        raise ValueError(f"resumen sin bloque 'features': {summary_path}")
    features["rRNA"] = int(n_rrna)
    summary_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Ejecución por genoma
# ---------------------------------------------------------------------------
def _cache_key(genome: str, cm_path: Path, strategy: str) -> str:
    digest = hashlib.sha256(cm_path.read_bytes()).hexdigest()[:16]
    return f"{genome}.{digest}.{strategy}"


def _load_cache(cache_dir: Path, key: str) -> dict[str, Any] | None:
    path = cache_dir / f"{key}.json"
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or "features" not in data:
        return None
    return data


def _save_cache(cache_dir: Path, key: str, data: dict[str, Any]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / f"{key}.json").write_text(
        json.dumps(data, ensure_ascii=False), encoding="utf-8"
    )


def reannotate_one(
    record: dict[str, Any],
    *,
    root: Path,
    cm_path: Path,
    cache_dir: Path | None,
    apply: bool,
    cpu: int,
    timeout: int,
    use_cache: bool,
    strategy: str,
) -> dict[str, Any]:
    genome = str(record.get("id") or "")
    row: dict[str, Any] = {
        "genome": genome,
        "species": str(record.get("species") or ""),
        "strain": str(record.get("strain") or ""),
        "state": "error",
        "tool": "",
        "strategy": strategy,
        "before_total": None,
        "before_16S": None,
        "before_23S": None,
        "before_5S": None,
        "after_total": None,
        "after_16S": None,
        "after_23S": None,
        "after_5S": None,
        "ref_total": None,
        "ref_16S": None,
        "ref_23S": None,
        "ref_5S": None,
        "elapsed_s": 0.0,
        "error": "",
    }
    started = time.perf_counter()
    try:
        genome_path = rp._resolve(root, record.get("genome_path"))
        if genome_path is None or not genome_path.is_file():
            raise FileNotFoundError("falta genome_path del genoma")
        system_path, system_display = cp.system_gff_path(record, root)
        if system_path is None:
            raise FileNotFoundError(f"falta el GFF del sistema: {system_display}")
        summary_raw = str(record.get("system_summary") or "").strip()
        summary_path = rp._resolve(root, summary_raw) if summary_raw else None
        if summary_path is not None and not summary_path.is_file():
            raise FileNotFoundError(f"falta el resumen del sistema: {summary_raw}")
        if summary_path is None:
            row["error"] = "sin system_summary: solo se reescribe el GFF"

        old_features = _old_rrna_features(system_path)
        before = _count_by_type(old_features)
        reference = _reference_rrna_counts(record, root)
        row.update(
            {
                "before_total": len(old_features),
                "before_16S": before["16S"],
                "before_23S": before["23S"],
                "before_5S": before["5S"],
                "ref_total": sum(reference.values()),
                "ref_16S": reference["16S"],
                "ref_23S": reference["23S"],
                "ref_5S": reference["5S"],
            }
        )

        key = _cache_key(genome, cm_path, strategy)
        cached = _load_cache(cache_dir, key) if (cache_dir and use_cache) else None
        if cached is not None:
            result = cached
        else:
            with tempfile.TemporaryDirectory(prefix=f"rrna_{genome}_") as td:
                fasta = Path(td) / f"{genome}.fna"
                fasta.write_text(_fasta_text(genome_path), encoding="utf-8")
                result = call_rrna_fasta(
                    fasta,
                    cm_path=cm_path,
                    cpu=cpu,
                    timeout=timeout,
                    strategy=strategy,
                )
            if cache_dir:
                _save_cache(cache_dir, key, result)

        new_features = result["features"]
        after = _count_by_type(new_features)
        row.update(
            {
                "state": "ok" if result.get("tool") not in ("", "none") else "empty",
                "tool": str(result.get("tool") or ""),
                "strategy": str(result.get("strategy") or strategy),
                "after_total": len(new_features),
                "after_16S": after["16S"],
                "after_23S": after["23S"],
                "after_5S": after["5S"],
            }
        )
        if apply:
            _rewrite_gff(system_path, new_features)
            if summary_path is not None:
                _update_summary(summary_path, len(new_features))
    except Exception as exc:  # noqa: BLE001 - se reporta por genoma
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["state"] = "error"
    row["elapsed_s"] = round(time.perf_counter() - started, 1)
    return row


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------
TSV_COLUMNS = (
    "genome",
    "species",
    "strain",
    "state",
    "tool",
    "strategy",
    "before_total",
    "before_16S",
    "before_23S",
    "before_5S",
    "after_total",
    "after_16S",
    "after_23S",
    "after_5S",
    "ref_total",
    "ref_16S",
    "ref_23S",
    "ref_5S",
    "elapsed_s",
    "error",
)


def write_tsv(rows: Sequence[dict[str, Any]], path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in TSV_COLUMNS})


def _md(value: Any) -> str:
    return "—" if value in (None, "") else str(value).replace("|", "\\|")


def render_markdown(rows: Sequence[dict[str, Any]], *, applied: bool) -> str:
    before_total = sum(int(row["before_total"] or 0) for row in rows)
    after_total = sum(int(row["after_total"] or 0) for row in rows)
    ref_total = sum(int(row["ref_total"] or 0) for row in rows)
    affected = [
        row for row in rows
        if row["before_total"] != row["after_total"]
        or row["before_23S"] != row["after_23S"]
        or row["before_5S"] != row["after_5S"]
    ]
    lines = [
        "# Reanotación de la capa rRNA en el panel de referencia",
        "",
        f"**Generado:** `{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}`  ",
        f"**Modo:** {'aplicado a results/' if applied else 'dry-run'}  ",
        f"**Genomas procesados:** {len(rows)}  ",
        f"**Genomas con cambio:** {len(affected)}",
        "",
        "## Causa raíz",
        "",
        "El GFF del sistema se generó con `barrnap` 0.9. En esta máquina su base "
        "`bac.hmm` solo contiene los modelos `16S_rRNA` (RF00177) y `5S_rRNA` "
        "(RF00001): **falta el modelo 23S** (`LSU_rRNA_bacteria`, RF02541). "
        "Además, el umbral por defecto de barrnap (`-E 1e-6`) descarta copias de "
        "5S cuyo mejor hit queda por encima de ese E-value (p. ej. E=0.023 en "
        "H. pylori). La referencia RefSeq/NCBI anota los tres tipos con "
        "Infernal/Rfam.",
        "",
        "La corrección vuelve a llamar **solo rRNA** sobre el FASTA del genoma "
        "con los modelos Rfam locales (`data/rfam_models/bac_rRNA.cm`) y "
        "`cmsearch --cut_ga`; no usa el GFF de referencia como entrada. Si "
        "`cmsearch` no está disponible, usa `barrnap --evalue 10` como respaldo.",
        "",
        "## Totales",
        "",
        "| | antes | después | referencia |",
        "|---|---:|---:|---:|",
        f"| rRNA | {before_total} | {after_total} | {ref_total} |",
        "",
    ]
    by_ref = {kind: sum(int(row[f"ref_{kind}"] or 0) for row in rows)
              for kind in CLASSES_BY_TYPE}
    by_before = {kind: sum(int(row[f"before_{kind}"] or 0) for row in rows)
                 for kind in CLASSES_BY_TYPE}
    by_after = {kind: sum(int(row[f"after_{kind}"] or 0) for row in rows)
                for kind in CLASSES_BY_TYPE}
    lines += [
        "| tipo | antes | después | referencia |",
        "|---|---:|---:|---:|",
    ]
    for kind in CLASSES_BY_TYPE:
        lines.append(
            f"| {kind} | {by_before[kind]} | {by_after[kind]} | {by_ref[kind]} |"
        )
    lines += [
        "",
        "## Por genoma",
        "",
        "| genoma | estado | motor | estrategia | antes 16S/23S/5S | "
        "después 16S/23S/5S | referencia 16S/23S/5S |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['state'])} | "
            f"{_md(row['tool'])} | {_md(row.get('strategy'))} | "
            f"{_md(row['before_16S'])}/{_md(row['before_23S'])}/"
            f"{_md(row['before_5S'])} | "
            f"{_md(row['after_16S'])}/{_md(row['after_23S'])}/"
            f"{_md(row['after_5S'])} | "
            f"{_md(row['ref_16S'])}/{_md(row['ref_23S'])}/"
            f"{_md(row['ref_5S'])} |"
        )
    errors = [row for row in rows if row["state"] == "error"]
    if errors:
        lines += ["", "## Errores", "", "| genoma | error |", "|---|---|"]
        for row in errors:
            lines.append(f"| `{_md(row['genome'])}` | {_md(row['error'])} |")
    lines += [
        "",
        "## Límites",
        "",
        "- Solo se reescriben las filas `rRNA` del GFF del sistema y "
        "`features.rRNA` del resumen; CDS, tRNA y el resto de capas no se "
        "tocan.",
        "- La referencia no se usa para llamar features: solo para comparar "
        "conteos en este informe.",
        "- `cmsearch` es sensible al modelo y al corte GA de Rfam; un rRNA muy "
        "divergente o fragmentado puede quedar fuera. No se imputan copias.",
        "- La reanotación es la única modificación de `results/`; no se toca "
        "`data/reference_genomes/`.",
        "",
        "## Reproducir",
        "",
        "```bash",
        ".venv/bin/python scripts/reference/reannotate_rrna.py --apply --jobs 4",
        ".venv/bin/python scripts/reference/reference_parity.py",
        ".venv/bin/python scripts/reference/coordinate_parity.py",
        "```",
        "",
    ]
    return "\n".join(lines)


def write_markdown(rows: Sequence[dict[str, Any]], path: Path | str, *,
                   applied: bool) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_markdown(rows, applied=applied), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--cm", type=Path, default=None,
                        help=f"modelos Rfam (default: data/rfam_models/bac_rRNA.cm)")
    parser.add_argument("--apply", action="store_true",
                        help="reescribe GFF del sistema y resumen (default: dry-run)")
    parser.add_argument("--jobs", type=int, default=1,
                        help="genomas en paralelo (default: 1)")
    parser.add_argument("--cpu", type=int, default=1,
                        help="CPUs por cmsearch (default: 1)")
    parser.add_argument("--timeout", type=int, default=3600,
                        help="timeout por genoma en segundos (default: 3600)")
    parser.add_argument(
        "--strategy",
        choices=["exact", "balanced", "recall", "fast"],
        default="exact",
        help=(
            "exact=CM 16S/23S/5S; balanced=CM 16S/5S+HMM 23S; "
            "recall=CM 5S+HMM 16S/23S; fast=HMM en los tres (default: exact)"
        ),
    )
    parser.add_argument("--only", action="append", default=[],
                        help="restringe a un id (repetible)")
    parser.add_argument("--cache-dir", type=Path, default=None,
                        help="caché de llamadas (default: .cache/rrna_calls)")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--tsv", type=Path, default=None)
    parser.add_argument("--md", type=Path, default=None)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = (args.root or ROOT).resolve()
    manifest = (args.manifest or root / "data" / "reference_genomes.json").resolve()
    cm_path = (args.cm or root / DEFAULT_RRNA_CM).resolve()
    cache_dir = args.cache_dir or (root / ".cache" / "rrna_calls")
    tsv_path = (args.tsv or root / "reports" / "rrna_reannotation.tsv").resolve()
    md_path = (args.md or root / "reports" / "rrna_reannotation.md").resolve()

    if not cm_path.is_file():
        print(f"ERROR: no existe el panel de modelos: {cm_path}", file=sys.stderr)
        return 2
    try:
        records = read_manifest(manifest)
    except (OSError, EOFError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    only = {str(item) for item in args.only}
    selected = [
        record for record in records
        if str(record.get("status") or "").strip().lower() in ELIGIBLE
        and (not only or str(record.get("id")) in only)
    ]
    if args.jobs > 1 and args.cpu > 1:
        print(
            "AVISO: con --jobs > 1 conviene --cpu 1 para no saturar la máquina",
            file=sys.stderr,
        )
    print(f"genomas a reanotar: {len(selected)}", flush=True)

    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {
            pool.submit(
                reannotate_one,
                record,
                root=root,
                cm_path=cm_path,
                cache_dir=cache_dir,
                apply=args.apply,
                cpu=max(1, args.cpu),
                timeout=max(1, args.timeout),
                use_cache=not args.no_cache,
                strategy=args.strategy,
            ): record
            for record in selected
        }
        for done, future in enumerate(as_completed(futures), 1):
            row = future.result()
            rows.append(row)
            print(
                f"  [{done}/{len(selected)}] {row['genome']}: "
                f"{row['tool']} {row['before_total']} -> {row['after_total']} "
                f"(ref {row['ref_total']}) en {row['elapsed_s']}s",
                flush=True,
            )

    rows.sort(key=lambda row: str(row["genome"]))
    write_tsv(rows, tsv_path)
    write_markdown(rows, md_path, applied=args.apply)
    before = sum(int(row["before_total"] or 0) for row in rows)
    after = sum(int(row["after_total"] or 0) for row in rows)
    errors = sum(1 for row in rows if row["state"] == "error")
    print(
        f"rRNA total: {before} -> {after}; errores: {errors}; "
        f"modo: {'apply' if args.apply else 'dry-run'}"
    )
    print(f"TSV: {tsv_path}")
    print(f"MD:  {md_path}")
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
