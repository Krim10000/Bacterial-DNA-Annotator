#!/usr/bin/env python3
"""migrate_hairpin_schema.py — corrige el esquema antiguo de horquillas.

En los TSV antiguos de ``data/operon_finder`` el campo ``strand`` contiene
``"intergénica"``. Eso es una **localización**, no una hebra: una horquilla
intergénica puede estar en la hebra ``+`` o ``-`` según el gen que la flanquea.

La convención del pipeline (la misma que usa ``annotate_hairpins.py`` y que el
navegador aplica en memoria con ``_hairpin_strands``) es:

* ``location`` conserva la localización (``"intergénica"``);
* ``strand`` es la hebra del gen izquierdo (``left_gene``) en
  ``data/regulation/<genoma>/genes.tsv``; si no existe, la del derecho;
* si no hay gen flanqueante, ``strand=""`` y
  ``strand_source="sin gen flanqueante"``. No se inventa una hebra.

El script conserva el orden y los valores de las columnas existentes y añade
``location`` y ``strand_source`` al final. Escribe de forma atómica (temporal en
el mismo directorio + ``os.replace``), es idempotente, admite ``--dry-run`` e
incluye ``data/operon_finder/excluded/`` cuando existe. No usa red.

Uso
---
    python scripts/regulation/migrate_hairpin_schema.py --dry-run
    python scripts/regulation/migrate_hairpin_schema.py
    python scripts/regulation/migrate_hairpin_schema.py \
        --dir data/operon_finder --genes-root data/regulation
"""
from __future__ import annotations

import argparse
import csv
import os
import stat
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DIR = ROOT / "data" / "operon_finder"
DEFAULT_GENES_ROOT = ROOT / "data" / "regulation"

LOCATION_INTERGENIC = "intergénica"
SOURCE_LEFT = "gen izquierdo"
SOURCE_RIGHT = "gen derecho"
SOURCE_NONE = "sin gen flanqueante"
SOURCE_ANNOTATION = "anotación"
NEW_COLUMNS = ("location", "strand_source")
HAIRPIN_SUFFIX = "_hairpins.tsv"


def normalize_strand(value) -> str:
    """Normaliza una hebra a ``"+"``/``"-"``; lo desconocido queda en ``""``.

    Acepta las convenciones de la anotación (``1``/``-1``) y variantes de
    texto. Es deliberadamente estricta: no adivina hebras a partir de otros
    campos.
    """
    text = str(value or "").strip().lower().replace("−", "-")
    if text in {"1", "+", "+1", "plus", "forward", "fwd", "directa"}:
        return "+"
    if text in {"-1", "-", "minus", "reverse", "rev", "reversa", "inversa"}:
        return "-"
    return ""


@dataclass
class GeneIndex:
    """Índice de hebras por identificador de gen.

    Un identificador que aparece con hebras distintas (p. ej. un ``gene_name``
    repetido en parálogos de hebras opuestas) se marca ambiguo y no resuelve:
    preferimos dejar ``strand=""`` antes que inventar una hebra.
    """

    strands: dict[str, str] = field(default_factory=dict)
    ambiguous: set[str] = field(default_factory=set)

    def resolve(self, identifier) -> str:
        key = str(identifier or "").strip()
        if not key or key in self.ambiguous:
            return ""
        return self.strands.get(key, "")


def _safe_genes_path(genome: str, genes_root: Path) -> Path | None:
    """``genes_root/<genoma>/genes.tsv`` sin permitir rutas fuera de la raíz."""
    genome = str(genome or "").strip()
    if not genome or genome in {".", ".."} or "/" in genome or "\\" in genome:
        return None
    return genes_root / genome / "genes.tsv"


def load_gene_index(path: Path | None) -> GeneIndex:
    """Carga hebras por ``locus_tag``, ``gene_id`` y ``gene_name``.

    Un gen puede ser referenciado por cualquiera de las tres columnas; los
    TSV de horquillas usan ``locus_tag`` en K-12 y ``gene_id`` en los genomas
    anotados con pyrodigal.
    """
    index = GeneIndex()
    if path is None or not path.is_file():
        return index
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            strand = normalize_strand(row.get("strand"))
            if not strand:
                continue
            for key in ("locus_tag", "gene_id", "gene_name"):
                value = str(row.get(key) or "").strip()
                if not value or value in index.ambiguous:
                    continue
                previous = index.strands.get(value)
                if previous is None:
                    index.strands[value] = strand
                elif previous != strand:
                    index.strands.pop(value, None)
                    index.ambiguous.add(value)
    return index


def migrate_row(row: dict, genes: GeneIndex) -> dict:
    """Devuelve una copia de ``row`` con ``strand``/``location``/``strand_source``.

    Si ``strand`` ya es una hebra real se conserva (normalizada) y no se toca
    la localización; así una segunda pasada no cambia nada. El texto antiguo
    del campo ``strand`` se conserva como ``location``.
    """
    old = str(row.get("strand") or "").strip()
    old_strand = normalize_strand(old)

    location = str(row.get("location") or "").strip()
    if not location:
        # El texto antiguo era la localización. Si no hay texto, la tabla de
        # horquillas solo contiene regiones intergénicas.
        location = old if old and not old_strand else LOCATION_INTERGENIC

    if old_strand:
        strand = old_strand
        source = str(row.get("strand_source") or "").strip() or SOURCE_ANNOTATION
    else:
        left = genes.resolve(row.get("left_gene"))
        if left:
            strand, source = left, SOURCE_LEFT
        else:
            right = genes.resolve(row.get("right_gene"))
            if right:
                strand, source = right, SOURCE_RIGHT
            else:
                strand, source = "", SOURCE_NONE

    out = dict(row)
    out["strand"] = strand
    out["location"] = location
    out["strand_source"] = source
    return out


def _read_rows(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    return fields, rows


def _line_terminator(path: Path) -> str:
    """Conserva el fin de línea del fichero (los TSV actuales usan CRLF)."""
    with path.open("rb") as fh:
        head = fh.read(65536)
    return "\r\n" if b"\r\n" in head else "\n"


def _atomic_write(path: Path, fields: list[str], rows: list[dict],
                  terminator: str) -> None:
    """Escribe en un temporal del mismo directorio y hace ``os.replace``.

    Conserva el modo del fichero original. Si algo falla, el temporal se borra
    y el TSV original queda intacto.
    """
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=fields, delimiter="\t",
                lineterminator=terminator)
            writer.writeheader()
            writer.writerows(rows)
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, stat.S_IMODE(path.stat().st_mode))
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


@dataclass
class FileStats:
    path: Path
    rows: int
    changed_rows: int
    written: bool
    strands: Counter
    sources: Counter
    genomes_without_genes: set[str]

    @property
    def changed(self) -> bool:
        return self.written or self.changed_rows > 0


def migrate_file(path: Path, genes_root: Path, dry_run: bool = False) -> FileStats:
    """Migra un ``*_hairpins.tsv``. No escribe nada si ``dry_run`` es True."""
    fields, rows = _read_rows(path)
    if not fields:
        return FileStats(path, 0, 0, False, Counter(), Counter(), set())
    original_fields = list(fields)
    for column in ("strand", *NEW_COLUMNS):
        if column not in fields:
            fields.append(column)

    # El genoma sale del nombre del fichero; si una fila trae otro, se respeta.
    default_genome = path.name[:-len(HAIRPIN_SUFFIX)] if path.name.endswith(
        HAIRPIN_SUFFIX) else path.stem
    cache: dict[str, GeneIndex] = {}
    missing: set[str] = set()
    new_rows: list[dict] = []
    changed_rows = 0
    for row in rows:
        genome = str(row.get("genome") or "").strip() or default_genome
        if genome not in cache:
            genes_path = _safe_genes_path(genome, genes_root)
            cache[genome] = load_gene_index(genes_path)
            if genes_path is None or not genes_path.is_file():
                missing.add(genome)
        migrated = migrate_row(row, cache[genome])
        if any(row.get(k, "") != migrated.get(k, "") for k in fields):
            changed_rows += 1
        new_rows.append(migrated)

    added_columns = [c for c in fields if c not in original_fields]
    needs_write = bool(changed_rows) or bool(added_columns)
    written = False
    if needs_write:
        if not dry_run:
            _atomic_write(path, fields, new_rows, _line_terminator(path))
        written = True

    return FileStats(
        path=path,
        rows=len(rows),
        changed_rows=changed_rows,
        written=written,
        strands=Counter(r.get("strand", "") for r in new_rows),
        sources=Counter(r.get("strand_source", "") for r in new_rows),
        genomes_without_genes=missing,
    )


def collect_files(base: Path) -> list[Path]:
    """Ficheros ``*_hairpins.tsv`` de ``base`` y de ``base/excluded`` si existe."""
    if not base.is_dir():
        return []
    files = sorted(base.glob(f"*{HAIRPIN_SUFFIX}"))
    excluded = base / "excluded"
    if excluded.is_dir():
        files.extend(sorted(excluded.glob(f"*{HAIRPIN_SUFFIX}")))
    return files


def _display(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _dir_stats(results: list[FileStats], label: str, verb: str) -> None:
    print(f"    {label}: {len(results)} archivos, "
          f"{sum(r.rows for r in results)} filas, "
          f"{sum(1 for r in results if r.changed)} {verb}")


def _print_summary(results: list[FileStats], base: Path, genes_root: Path,
                   dry_run: bool) -> None:
    verb = "a cambiar" if dry_run else "reescritos"
    changed_files = sum(1 for r in results if r.changed)
    total_rows = sum(r.rows for r in results)
    changed_rows = sum(r.changed_rows for r in results)
    strands = Counter()
    sources = Counter()
    missing: set[str] = set()
    for r in results:
        strands.update(r.strands)
        sources.update(r.sources)
        missing |= r.genomes_without_genes

    excluded = base / "excluded"

    print("migrate_hairpin_schema: "
          + ("dry-run (no se escribe nada)" if dry_run
             else "escritura atómica"))
    print(f"  directorio: {_display(base)}"
          + (" (+ excluded)" if excluded.is_dir() else ""))
    print(f"  genes: {_display(genes_root)}/<genoma>/genes.tsv")
    print()
    print(f"  archivos: {len(results)} escaneados, "
          f"{changed_files} {verb}, {len(results) - changed_files} sin cambios")
    print(f"  filas: {total_rows} total, {changed_rows} "
          f"{'a cambiar' if dry_run else 'cambiadas'}")
    print(f"  hebra: + {strands.get('+', 0)}, - {strands.get('-', 0)}, "
          f"sin hebra {strands.get('', 0)}")
    detail = (f"gen izquierdo {sources.get(SOURCE_LEFT, 0)}, "
              f"gen derecho {sources.get(SOURCE_RIGHT, 0)}, "
              f"sin gen flanqueante {sources.get(SOURCE_NONE, 0)}")
    if sources.get(SOURCE_ANNOTATION):
        detail += f", anotación {sources[SOURCE_ANNOTATION]}"
    print(f"  fuente: {detail}")
    if missing:
        names = ", ".join(sorted(missing))
        print(f"  genomas sin genes.tsv ({len(missing)}): {names}")

    if excluded.is_dir():
        print("  por directorio:")
        _dir_stats([r for r in results if r.path.parent == base],
                   "principal", verb)
        _dir_stats([r for r in results if r.path.parent == excluded],
                   "excluded", verb)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Migra strand='intergénica' a location + strand real.")
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR,
                    help="carpeta con *_hairpins.tsv (default data/operon_finder)")
    ap.add_argument("--genes-root", type=Path, default=DEFAULT_GENES_ROOT,
                    help="raíz de data/regulation (default data/regulation)")
    ap.add_argument("--dry-run", action="store_true",
                    help="calcula el resumen sin escribir ningún fichero")
    args = ap.parse_args(argv)

    if not args.dir.is_dir():
        print(f"no existe el directorio: {args.dir}", file=sys.stderr)
        return 2
    files = collect_files(args.dir)
    if not files:
        print(f"no hay *{HAIRPIN_SUFFIX} en {args.dir}", file=sys.stderr)
        return 0

    results: list[FileStats] = []
    errors = 0
    for path in files:
        try:
            results.append(migrate_file(path, args.genes_root, args.dry_run))
        except (OSError, csv.Error, UnicodeDecodeError) as exc:
            errors += 1
            print(f"ERROR {path}: {exc}", file=sys.stderr)

    _print_summary(results, args.dir, args.genes_root, args.dry_run)
    if errors:
        print(f"{errors} ficheros con error; no se modificaron", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
