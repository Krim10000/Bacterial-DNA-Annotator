"""Tests de ``scripts/reference/coordinate_parity.py`` con datos sintéticos.

No usan red, GPU ni binarios externos. Cubren paridad exacta, tolerancia de
bordes, solape recíproco, TP/FP/FN y tasas, features multiparte, GFF con gzip,
hebras/contigs distintos, coordenadas invertidas en hebra negativa, ausencia de
GFF, GFF sin clases objetivo, secciones ``##FASTA``, filtrado del manifiesto,
escritura de TSV/Markdown y no modificación de las entradas.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REF_DIR = ROOT / "scripts" / "reference"
if str(REF_DIR) not in sys.path:
    sys.path.insert(0, str(REF_DIR))

import coordinate_parity as cp  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
Row = tuple[str, str, int, int, str, str]


def make_row(
    kind: str,
    start: int,
    end: int,
    fid: str,
    *,
    seqid: str = "ctg",
    strand: str = "+",
) -> Row:
    return (seqid, kind, start, end, strand, fid)


def write_gff(path: Path, rows: list[Row], *, gz: bool = False) -> Path:
    """Escribe un GFF3 mínimo con filas (seqid, tipo, inicio, fin, hebra, ID)."""
    lines = ["##gff-version 3"]
    for seqid, kind, start, end, strand, fid in rows:
        phase = "0" if kind == "CDS" else "."
        attrs = f"ID={fid}" if fid else "."
        lines.append(
            "\t".join(
                [
                    seqid,
                    "test",
                    kind,
                    str(start),
                    str(end),
                    ".",
                    strand,
                    phase,
                    attrs,
                ]
            )
        )
    text = "\n".join(lines) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _rel(root: Path, path: Path) -> str:
    return str(path.relative_to(root))


def setup_genome(
    root: Path,
    *,
    genome: str = "g1",
    status: str = "annotated",
    system_rows: list[Row] | None = None,
    reference_rows: list[Row] | None = None,
    system_gz: bool = False,
    annotation_path: bool = True,
    system_result: bool = True,
    system_dir_name: str | None = None,
) -> tuple[Path, Path, Path, Path]:
    """Crea manifiesto y GFF sintéticos; devuelve sus rutas."""
    sys_dir = root / "results" / (system_dir_name or genome)
    suffix = ".gff3.gz" if system_gz else ".gff3"
    system_path = sys_dir / f"{genome}{suffix}"
    reference_path = root / "data" / "ref" / f"{genome}_genomic.gff.gz"
    if system_rows is not None:
        write_gff(system_path, system_rows, gz=system_gz)
    if reference_rows is not None:
        write_gff(reference_path, reference_rows, gz=True)
    record = {
        "id": genome,
        "species": "Synthetic species",
        "strain": "T",
        "status": status,
        "annotation_path": _rel(root, reference_path) if annotation_path else None,
        "system_result": _rel(root, sys_dir) if system_result else None,
    }
    manifest = root / "data" / "reference_genomes.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps([record]), encoding="utf-8")
    return root, manifest, system_path, reference_path


def exact_rows() -> list[Row]:
    return [
        make_row("CDS", 100, 399, "cds1"),
        make_row("CDS", 500, 799, "cds2", strand="-"),
        make_row("tRNA", 900, 970, "trna1"),
        make_row("rRNA", 1000, 1150, "rrna1"),
    ]


# ---------------------------------------------------------------------------
# Paridad y tolerancia
# ---------------------------------------------------------------------------
def test_exact_parity_is_ok_and_metrics_are_perfect(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=rows, reference_rows=rows
    )

    result = cp.build_report(manifest, root, generated_at="2026-01-01T00:00:00Z")

    assert result["generated_at"] == "2026-01-01T00:00:00Z"
    assert result["summary"] == {
        "candidates": 1,
        "comparable": 1,
        "ok": 1,
        "partial": 0,
        "not_available": 0,
        "skipped": 0,
    }
    row = result["rows"][0]
    assert row["state"] == "ok"
    assert row["classes"]["CDS"]["state"] == "ok"
    assert row["classes"]["tRNA"]["state"] == "ok"
    assert row["classes"]["rRNA"]["state"] == "ok"

    pooled = row["all"]
    assert pooled["system_features"] == 4
    assert pooled["reference_features"] == 4
    assert pooled["tp"] == 4
    assert pooled["fp"] == 0
    assert pooled["fn"] == 0
    assert pooled["precision"] == 1.0
    assert pooled["recall"] == 1.0
    assert pooled["f1"] == 1.0
    assert pooled["exact_edges"] == 4
    assert pooled["exact_pct"] == 100.0
    assert pooled["mean_start_diff"] == 0
    assert pooled["median_start_diff"] == 0
    assert pooled["mean_end_diff"] == 0
    assert pooled["median_end_diff"] == 0
    assert pooled["max_abs_edge_diff"] == 0


def test_shift_within_tolerance_is_matched_and_ok(tmp_path: Path) -> None:
    system = [make_row("CDS", 103, 397, "cds1")]
    reference = [make_row("CDS", 100, 400, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root, tolerance=5)
    metrics = result["rows"][0]["all"]
    assert result["rows"][0]["state"] == "ok"
    assert metrics["tp"] == 1
    assert metrics["fp"] == 0
    assert metrics["fn"] == 0
    assert metrics["exact_edges"] == 0
    assert metrics["exact_pct"] == 0.0
    assert metrics["mean_start_diff"] == 3.0
    assert metrics["median_start_diff"] == 3.0
    assert metrics["mean_end_diff"] == -3.0
    assert metrics["median_end_diff"] == -3.0
    assert metrics["max_abs_edge_diff"] == 3

    stricter = cp.build_report(manifest, root, tolerance=2)
    strict_metrics = stricter["rows"][0]["all"]
    assert stricter["rows"][0]["state"] == "partial"
    assert strict_metrics["tp"] == 1
    assert strict_metrics["max_abs_edge_diff"] == 3


def test_edge_tolerance_rescues_small_feature(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 103, "cds1")]
    reference = [make_row("CDS", 98, 108, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    tolerant = cp.build_report(manifest, root, tolerance=5)
    metrics = tolerant["rows"][0]["all"]
    assert tolerant["rows"][0]["state"] == "ok"
    assert metrics["tp"] == 1
    assert metrics["max_abs_edge_diff"] == 5
    assert metrics["exact_edges"] == 0

    strict = cp.build_report(manifest, root, tolerance=4)
    strict_metrics = strict["rows"][0]["all"]
    assert strict["rows"][0]["state"] == "partial"
    assert strict_metrics["tp"] == 0
    assert strict_metrics["fp"] == 1
    assert strict_metrics["fn"] == 1
    assert strict_metrics["precision"] == 0.0
    assert strict_metrics["recall"] == 0.0
    assert strict_metrics["f1"] == 0.0


def test_reciprocal_overlap_threshold_is_configurable(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 199, "cds1")]
    reference = [make_row("CDS", 160, 259, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    default = cp.build_report(manifest, root)  # solape 40 %, umbral 50 %
    default_metrics = default["rows"][0]["all"]
    assert default_metrics["tp"] == 0
    assert default_metrics["fp"] == 1
    assert default_metrics["fn"] == 1

    loose = cp.build_report(manifest, root, min_reciprocal_overlap=0.3)
    loose_metrics = loose["rows"][0]["all"]
    assert loose_metrics["tp"] == 1
    assert loose_metrics["precision"] == 1.0
    assert loose_metrics["recall"] == 1.0
    assert loose_metrics["max_abs_edge_diff"] == 60
    assert loose["rows"][0]["state"] == "partial"


# ---------------------------------------------------------------------------
# TP/FP/FN y tasas
# ---------------------------------------------------------------------------
def test_fp_fn_precision_recall_and_f1(tmp_path: Path) -> None:
    system = [
        make_row("CDS", 100, 199, "cds1"),
        make_row("CDS", 300, 399, "cds2"),
        make_row("tRNA", 700, 770, "trna1"),
    ]
    reference = [
        make_row("CDS", 100, 199, "ref_cds1"),
        make_row("CDS", 500, 599, "ref_cds2"),
        make_row("rRNA", 800, 900, "ref_rrna1"),
    ]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    pooled = row["all"]
    assert row["state"] == "partial"
    assert pooled["tp"] == 1
    assert pooled["fp"] == 2
    assert pooled["fn"] == 2
    assert pooled["precision"] == 0.3333
    assert pooled["recall"] == 0.3333
    assert pooled["f1"] == 0.3333

    cds = row["classes"]["CDS"]
    assert cds["tp"] == 1
    assert cds["fp"] == 1
    assert cds["fn"] == 1
    assert cds["precision"] == 0.5
    assert cds["recall"] == 0.5

    trna = row["classes"]["tRNA"]
    assert trna["state"] == "partial"
    assert trna["system_features"] == 1
    assert trna["reference_features"] == 0
    assert trna["tp"] == 0
    assert trna["fp"] == 1
    assert trna["fn"] == 0
    assert trna["precision"] == 0.0
    assert trna["recall"] is None

    rrna = row["classes"]["rRNA"]
    assert rrna["system_features"] == 0
    assert rrna["reference_features"] == 1
    assert rrna["tp"] == 0
    assert rrna["fp"] == 0
    assert rrna["fn"] == 1
    assert rrna["precision"] is None
    assert rrna["recall"] == 0.0


def test_strand_mismatch_is_not_matched(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 199, "cds1", strand="+")]
    reference = [make_row("CDS", 100, 199, "ref1", strand="-")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    assert result["rows"][0]["state"] == "partial"
    metrics = result["rows"][0]["all"]
    assert metrics["tp"] == 0
    assert metrics["fp"] == 1
    assert metrics["fn"] == 1


def test_contig_mismatch_is_not_matched(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 199, "cds1", seqid="ctgA")]
    reference = [make_row("CDS", 100, 199, "ref1", seqid="ctgB")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    metrics = result["rows"][0]["all"]
    assert metrics["tp"] == 0
    assert metrics["fp"] == 1
    assert metrics["fn"] == 1


# ---------------------------------------------------------------------------
# Multiparte y coordenadas
# ---------------------------------------------------------------------------
def test_multipart_system_cds_is_collapsed(tmp_path: Path) -> None:
    system = [
        make_row("CDS", 100, 150, "cds1"),
        make_row("CDS", 200, 250, "cds1"),
    ]
    reference = [make_row("CDS", 100, 250, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    metrics = result["rows"][0]["all"]
    assert metrics["system_features"] == 1
    assert metrics["reference_features"] == 1
    assert metrics["tp"] == 1
    assert metrics["exact_edges"] == 1
    assert result["rows"][0]["state"] == "ok"


def test_multipart_reference_cds_is_collapsed(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 250, "cds1")]
    reference = [
        make_row("CDS", 100, 150, "ref1"),
        make_row("CDS", 200, 250, "ref1"),
    ]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    metrics = result["rows"][0]["all"]
    assert metrics["system_features"] == 1
    assert metrics["reference_features"] == 1
    assert metrics["tp"] == 1


def test_same_id_on_different_contigs_is_not_collapsed(tmp_path: Path) -> None:
    system = [
        make_row("CDS", 100, 150, "cds1", seqid="ctgA"),
        make_row("CDS", 100, 150, "cds1", seqid="ctgB"),
    ]
    reference = [
        make_row("CDS", 100, 150, "refA", seqid="ctgA"),
        make_row("CDS", 100, 150, "refB", seqid="ctgB"),
    ]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    metrics = result["rows"][0]["all"]
    assert metrics["system_features"] == 2
    assert metrics["reference_features"] == 2
    assert metrics["tp"] == 2
    assert result["rows"][0]["state"] == "ok"


def test_minus_strand_reversed_coordinates_are_normalized(
    tmp_path: Path,
) -> None:
    # El pipeline puede emitir la hebra negativa con inicio > fin.
    system = [make_row("CDS", 250, 200, "cds1", strand="-")]
    reference = [make_row("CDS", 200, 250, "ref1", strand="-")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    parsed = cp.parse_gff(_sys)
    assert parsed.features[0].start == 200
    assert parsed.features[0].end == 250

    result = cp.build_report(manifest, root)
    metrics = result["rows"][0]["all"]
    assert metrics["tp"] == 1
    assert metrics["exact_edges"] == 1
    assert result["rows"][0]["state"] == "ok"


# ---------------------------------------------------------------------------
# Ausencias y entradas incompletas
# ---------------------------------------------------------------------------
def test_missing_system_gff_is_not_available(tmp_path: Path) -> None:
    reference = [make_row("CDS", 100, 199, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=None, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "falta el GFF del sistema" in row["reason"]
    assert row["all"]["tp"] is None
    assert row["all"]["system_features"] is None
    assert result["summary"]["not_available"] == 1
    assert result["summary"]["comparable"] == 0


def test_missing_reference_gff_is_not_available(tmp_path: Path) -> None:
    system = [make_row("CDS", 100, 199, "cds1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=None
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "falta el GFF de referencia" in row["reason"]
    assert row["all"]["tp"] is None


def test_manifest_without_annotation_path_is_not_available(
    tmp_path: Path,
) -> None:
    system = [make_row("CDS", 100, 199, "cds1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path,
        system_rows=system,
        reference_rows=[make_row("CDS", 100, 199, "ref1")],
        annotation_path=False,
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "no declara annotation_path" in row["reason"]


def test_reference_without_target_classes_is_not_available(
    tmp_path: Path,
) -> None:
    system = [make_row("CDS", 100, 199, "cds1")]
    reference = [
        ("ctg", "gene", 100, 199, "+", "gene1"),
        ("ctg", "exon", 100, 150, "+", "exon1"),
        ("ctg", "ncRNA", 300, 400, "+", "ncrna1"),
    ]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "no contiene CDS/tRNA/rRNA" in row["reason"]


def test_empty_system_gff_is_comparable_with_all_fn(tmp_path: Path) -> None:
    root, manifest, _sys, _ref = setup_genome(
        tmp_path,
        system_rows=[],
        reference_rows=[make_row("CDS", 100, 199, "ref1")],
    )

    result = cp.build_report(manifest, root)
    row = result["rows"][0]
    assert row["state"] == "partial"
    metrics = row["all"]
    assert metrics["system_features"] == 0
    assert metrics["reference_features"] == 1
    assert metrics["tp"] == 0
    assert metrics["fp"] == 0
    assert metrics["fn"] == 1
    assert metrics["precision"] is None
    assert metrics["recall"] == 0.0


def test_gff_comments_fasta_and_invalid_lines_are_handled(
    tmp_path: Path,
) -> None:
    path = tmp_path / "weird.gff3"
    path.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "# comentario",
                "ctg\ttest\tCDS\t100\t199\t.\t+\t0\tID=cds1",
                "linea invalida",
                "ctg\ttest\tCDS\tnope\t200\t.\t+\t0\tID=bad",
                "ctg\ttest\tCDS\t300\tbad\t.\t+\t0\tID=bad2",
                "##FASTA",
                ">ctg",
                "ACGTACGT",
                "ctg\ttest\tCDS\t900\t999\t.\t+\t0\tID=fasta_cds",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    parsed = cp.parse_gff(path)
    assert len(parsed.features) == 1
    assert parsed.features[0].feature_id == "cds1"
    assert parsed.invalid_lines == 3
    assert parsed.raw_counts["CDS"] == 1


def test_gzip_system_and_reference_gffs(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path,
        system_rows=rows,
        reference_rows=rows,
        system_gz=True,
    )

    result = cp.build_report(manifest, root)
    assert result["rows"][0]["state"] == "ok"
    assert result["rows"][0]["all"]["tp"] == 4


def test_declared_system_result_directory(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, system_path, _ref = setup_genome(
        tmp_path,
        system_rows=rows,
        reference_rows=rows,
        system_dir_name="custom_system_out",
    )

    result = cp.build_report(manifest, root)
    assert result["rows"][0]["state"] == "ok"
    assert result["rows"][0]["system_gff"] == _rel(root, system_path)


# ---------------------------------------------------------------------------
# Manifiesto y agregados
# ---------------------------------------------------------------------------
def test_status_filtering_skips_downloaded(tmp_path: Path) -> None:
    rows = exact_rows()
    root = tmp_path
    records = []
    for genome, status in (
        ("annotated_one", "annotated"),
        ("downloaded_one", "downloaded"),
    ):
        sys_dir = root / "results" / genome
        ref_path = root / "data" / "ref" / f"{genome}_genomic.gff.gz"
        write_gff(sys_dir / f"{genome}.gff3", rows)
        write_gff(ref_path, rows, gz=True)
        records.append(
            {
                "id": genome,
                "species": "Synthetic species",
                "strain": "T",
                "status": status,
                "annotation_path": _rel(root, ref_path),
                "system_result": _rel(root, sys_dir),
            }
        )
    manifest = root / "data" / "reference_genomes.json"
    manifest.write_text(json.dumps(records), encoding="utf-8")

    result = cp.build_report(manifest, root)
    assert result["summary"]["candidates"] == 1
    assert result["summary"]["skipped"] == 1
    assert [row["genome"] for row in result["rows"]] == ["annotated_one"]


def test_available_status_is_included(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, status="available", system_rows=rows, reference_rows=rows
    )

    result = cp.build_report(manifest, root)
    assert result["summary"]["candidates"] == 1
    assert result["rows"][0]["manifest_status"] == "available"
    assert result["rows"][0]["state"] == "ok"


def test_aggregate_micro_metrics(tmp_path: Path) -> None:
    root = tmp_path
    rows = [make_row("CDS", 100, 199, "cds1"), make_row("CDS", 300, 399, "cds2")]
    _r, manifest, _s, _ref = setup_genome(
        root, genome="exact", system_rows=rows, reference_rows=rows
    )
    system2 = [make_row("CDS", 100, 199, "cds1")]
    reference2 = [make_row("CDS", 500, 599, "ref1")]
    rec2 = {
        "id": "extra",
        "species": "Synthetic species",
        "strain": "T",
        "status": "annotated",
        "annotation_path": "data/ref/extra_genomic.gff.gz",
        "system_result": "results/extra",
    }
    write_gff(root / "results" / "extra" / "extra.gff3", system2)
    write_gff(root / "data" / "ref" / "extra_genomic.gff.gz", reference2, gz=True)
    manifest.write_text(
        json.dumps([json.loads(manifest.read_text())[0], rec2]),
        encoding="utf-8",
    )

    result = cp.build_report(manifest, root)
    agg = result["aggregate"]["CDS"]
    assert agg["system_features"] == 3
    assert agg["reference_features"] == 3
    assert agg["tp"] == 2
    assert agg["fp"] == 1
    assert agg["fn"] == 1
    assert agg["precision"] == 0.6667
    assert agg["recall"] == 0.6667
    assert agg["f1"] == 0.6667
    assert agg["exact_edges"] == 2


# ---------------------------------------------------------------------------
# Salidas y CLI
# ---------------------------------------------------------------------------
def test_tsv_and_markdown_outputs(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=rows, reference_rows=rows
    )
    result = cp.build_report(manifest, root)
    tsv_path = root / "reports" / "coordinate.tsv"
    md_path = root / "reports" / "coordinate.md"

    written_tsv, written_md = cp.write_reports(result, tsv_path, md_path)
    assert written_tsv == tsv_path
    assert written_md == md_path
    assert tsv_path.is_file()
    assert md_path.is_file()

    with tsv_path.open(encoding="utf-8", newline="") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert tuple(table[0].keys()) == cp.TSV_COLUMNS
    assert [row["class"] for row in table] == ["CDS", "tRNA", "rRNA", "ALL"]
    assert table[-1]["state"] == "ok"
    assert table[-1]["tp"] == "4"
    assert table[-1]["precision"] == "1.0000"
    assert table[-1]["exact_pct"] == "100.00"
    assert table[-1]["system_gff"].endswith("g1.gff3")

    markdown = md_path.read_text(encoding="utf-8")
    assert "# Paridad por coordenadas" in markdown
    assert "## Resultados por genoma" in markdown
    assert "## Agregado (micro)" in markdown
    assert "## Límites" in markdown
    assert "`g1`" in markdown
    assert "100.00%" in markdown


def test_not_available_rows_are_written_with_blanks(tmp_path: Path) -> None:
    reference = [make_row("CDS", 100, 199, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=None, reference_rows=reference
    )
    result = cp.build_report(manifest, root)
    tsv_path = root / "reports" / "coordinate.tsv"
    cp.write_tsv(result["rows"], tsv_path)

    with tsv_path.open(encoding="utf-8", newline="") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert len(table) == 4
    assert all(row["state"] == "not_available" for row in table)
    assert all(row["tp"] == "" for row in table)
    assert all(row["system_features"] == "" for row in table)


def test_inputs_are_not_modified(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, system_path, reference_path = setup_genome(
        tmp_path, system_rows=rows, reference_rows=rows
    )
    before = {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in (manifest, system_path, reference_path)
    }

    result = cp.build_report(manifest, root)
    cp.write_reports(
        result, root / "reports" / "coordinate.tsv", root / "reports" / "coordinate.md"
    )

    for path, (data, mtime) in before.items():
        assert path.read_bytes() == data
        assert path.stat().st_mtime_ns == mtime


def test_cli_writes_reports(tmp_path: Path, capsys) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=rows, reference_rows=rows
    )
    tsv_path = root / "out" / "coord.tsv"
    md_path = root / "out" / "coord.md"

    exit_code = cp.main(
        [
            "--root",
            str(root),
            "--manifest",
            str(manifest),
            "--tsv",
            str(tsv_path),
            "--md",
            str(md_path),
        ]
    )

    assert exit_code == 0
    assert tsv_path.is_file()
    assert md_path.is_file()
    captured = capsys.readouterr()
    assert "1 comparables" in captured.out
    assert "TSV:" in captured.out


def test_cli_tolerance_option_changes_state(tmp_path: Path) -> None:
    system = [make_row("CDS", 103, 397, "cds1")]
    reference = [make_row("CDS", 100, 400, "ref1")]
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=system, reference_rows=reference
    )

    tolerant_tsv = root / "out" / "tolerant.tsv"
    strict_tsv = root / "out" / "strict.tsv"
    common = [
        "--root",
        str(root),
        "--manifest",
        str(manifest),
        "--md",
        str(root / "out" / "report.md"),
    ]
    assert cp.main(common + ["--tolerance", "5", "--tsv", str(tolerant_tsv)]) == 0
    assert cp.main(common + ["--tolerance", "2", "--tsv", str(strict_tsv)]) == 0

    with tolerant_tsv.open(encoding="utf-8", newline="") as fh:
        tolerant = [row for row in csv.DictReader(fh, delimiter="\t") if row["class"] == "ALL"]
    with strict_tsv.open(encoding="utf-8", newline="") as fh:
        strict = [row for row in csv.DictReader(fh, delimiter="\t") if row["class"] == "ALL"]
    assert tolerant[0]["state"] == "ok"
    assert strict[0]["state"] == "partial"


def test_cli_rejects_negative_tolerance() -> None:
    assert cp.main(["--tolerance", "-1"]) == 2


def test_cli_rejects_invalid_overlap_fraction() -> None:
    assert cp.main(["--min-reciprocal-overlap", "1.5"]) == 2


def test_invalid_parameters_in_build_report(tmp_path: Path) -> None:
    rows = exact_rows()
    root, manifest, _sys, _ref = setup_genome(
        tmp_path, system_rows=rows, reference_rows=rows
    )

    with pytest.raises(ValueError):
        cp.build_report(manifest, root, tolerance=-1)
    with pytest.raises(ValueError):
        cp.build_report(manifest, root, min_reciprocal_overlap=2.0)
