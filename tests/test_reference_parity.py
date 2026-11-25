"""Tests de ``scripts/reference/reference_parity.py`` con datos sintéticos.

No usan red, GPU ni binarios externos. Cubren paridad exacta, deltas,
porcentajes, CDS multiparte, gzip, ausencia de resumen/GFF, resumen ilegible,
GFF sin clases objetivo y el filtrado del manifiesto.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF_DIR = ROOT / "scripts" / "reference"
if str(REF_DIR) not in sys.path:
    sys.path.insert(0, str(REF_DIR))

import reference_parity as rp  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_gff(
    path: Path,
    rows: list[tuple[str, int, int, str]],
    *,
    gz: bool = False,
) -> Path:
    """Escribe un GFF3 mínimo con filas (tipo, inicio, fin, atributos)."""
    lines = ["##gff-version 3", "##sequence-region ctg 1 10000"]
    for ftype, start, end, attrs in rows:
        phase = "0" if ftype == "CDS" else "."
        lines.append(
            "\t".join(
                ["ctg", "test", ftype, str(start), str(end), ".", "+", phase, attrs]
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


def write_summary(path: Path, *, cds: int, trna: int, rrna: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "genome": path.parent.name,
                "features": {"CDS": cds, "tRNA": trna, "rRNA": rrna, "ncRNA": 0},
            }
        ),
        encoding="utf-8",
    )
    return path


def default_rows(*, cds: int = 3, trna: int = 1, rrna: int = 1) -> list[tuple[str, int, int, str]]:
    rows: list[tuple[str, int, int, str]] = []
    start = 100
    for i in range(cds):
        rows.append(("CDS", start, start + 299, f"ID=cds{i};product=protein {i}"))
        start += 400
    for i in range(trna):
        rows.append(("tRNA", start, start + 70, f"ID=tRNA{i};product=tRNA-Ala"))
        start += 100
    for i in range(rrna):
        rows.append(("rRNA", start, start + 100, f"ID=rRNA{i};product=16S ribosomal RNA"))
        start += 200
    return rows


def setup_root(
    tmp_path: Path,
    records: list[dict],
    *,
    summary_counts: tuple[int, int, int] | None = (3, 1, 1),
    gff_rows: list[tuple[str, int, int, str]] | None = None,
    gz: bool = True,
    write_summary_file: bool = True,
    write_gff_file: bool = True,
) -> Path:
    root = tmp_path / "proj"
    (root / "data").mkdir(parents=True, exist_ok=True)
    (root / "data" / "reference_genomes.json").write_text(
        json.dumps(records), encoding="utf-8"
    )
    if summary_counts is not None and write_summary_file:
        write_summary(
            root / "results" / "demo" / "demo.summary.json",
            cds=summary_counts[0],
            trna=summary_counts[1],
            rrna=summary_counts[2],
        )
    if write_gff_file:
        write_gff(
            root / "data" / "reference_genomes" / "demo" / "demo_genomic.gff.gz",
            gff_rows if gff_rows is not None else default_rows(),
            gz=gz,
        )
    return root


def demo_record(**overrides) -> dict:
    record = {
        "id": "demo",
        "species": "Demo bacterium",
        "strain": "D1",
        "status": "annotated",
        "annotation_path": "data/reference_genomes/demo/demo_genomic.gff.gz",
        "system_summary": "results/demo/demo.summary.json",
    }
    record.update(overrides)
    return record


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Paridad y métricas
# ---------------------------------------------------------------------------
def test_exact_counts_ok(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()])
    result = rp.build_report(root / "data" / "reference_genomes.json", root)

    assert result["summary"] == {
        "candidates": 1,
        "selected": 1,
        "ok": 1,
        "partial": 0,
        "not_available": 0,
        "skipped": 0,
    }
    row = result["rows"][0]
    assert row["state"] == "ok"
    assert row["reason"] == "conteos idénticos en CDS, tRNA y rRNA"
    for kind in rp.CLASSES:
        cell = row["classes"][kind]
        assert cell == {"system": 3 if kind == "CDS" else 1,
                        "reference": 3 if kind == "CDS" else 1,
                        "delta": 0, "pct": 100.0}


def test_partial_deltas_and_percentages(tmp_path: Path) -> None:
    ref_rows = default_rows(cds=8, trna=2, rrna=2)
    root = setup_root(
        tmp_path,
        [demo_record()],
        summary_counts=(10, 2, 1),
        gff_rows=ref_rows,
    )
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]

    assert row["state"] == "partial"
    assert row["classes"]["CDS"]["delta"] == 2
    assert row["classes"]["CDS"]["pct"] == 125.0
    assert row["classes"]["tRNA"]["delta"] == 0
    assert row["classes"]["tRNA"]["pct"] == 100.0
    assert row["classes"]["rRNA"]["delta"] == -1
    assert row["classes"]["rRNA"]["pct"] == 50.0
    assert "CDS=+2" in row["reason"]
    assert "rRNA=-1" in row["reason"]
    assert result["summary"]["partial"] == 1


def test_multipart_cds_is_counted_as_one_feature(tmp_path: Path) -> None:
    ref_rows = [
        ("CDS", 100, 200, "ID=cds1;product=DNA polymerase"),
        ("CDS", 300, 400, "ID=cds1;product=DNA polymerase"),
        ("CDS", 500, 799, "ID=cds2;product=hypothetical protein"),
        ("tRNA", 900, 970, "ID=trna1;product=tRNA-Ala"),
        ("rRNA", 1000, 1100, "ID=rrna1;product=16S ribosomal RNA"),
    ]
    root = setup_root(
        tmp_path,
        [demo_record()],
        summary_counts=(2, 1, 1),
        gff_rows=ref_rows,
    )
    parsed = rp.count_gff_features(
        root / "data" / "reference_genomes" / "demo" / "demo_genomic.gff.gz"
    )
    assert parsed["raw_counts"]["CDS"] == 3
    assert parsed["counts"] == {"CDS": 2, "tRNA": 1, "rRNA": 1}

    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    assert result["rows"][0]["state"] == "ok"


def test_plain_gff_is_supported(tmp_path: Path) -> None:
    root = setup_root(
        tmp_path,
        [demo_record(annotation_path="data/reference_genomes/demo/demo.gff3")],
    )
    # El helper escribió .gff.gz; se descomprime a .gff3 para el caso plano.
    gz_path = root / "data" / "reference_genomes" / "demo" / "demo_genomic.gff.gz"
    plain_path = root / "data" / "reference_genomes" / "demo" / "demo.gff3"
    plain_path.write_bytes(gzip.decompress(gz_path.read_bytes()))
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    assert result["rows"][0]["state"] == "ok"


def test_reference_without_a_class_is_partial_and_pct_blank(tmp_path: Path) -> None:
    ref_rows = [row for row in default_rows() if row[0] != "tRNA"]
    root = setup_root(
        tmp_path,
        [demo_record()],
        summary_counts=(3, 0, 1),
        gff_rows=ref_rows,
    )
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "partial"
    assert row["classes"]["tRNA"]["reference"] == 0
    assert row["classes"]["tRNA"]["pct"] is None
    assert "la referencia no tiene tRNA" in row["reason"]
    # La clase sin referencia no se agrega al micro-promedio.
    assert result["aggregate"]["tRNA"] == {
        "system": 0,
        "reference": 0,
        "delta": 0,
        "pct": None,
    }
    assert result["aggregate"]["CDS"]["pct"] == 100.0


# ---------------------------------------------------------------------------
# Ausencias y errores de entrada
# ---------------------------------------------------------------------------
def test_missing_system_summary_is_not_available(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()], write_summary_file=False)
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "falta el resumen del sistema" in row["reason"]
    assert row["classes"]["CDS"]["system"] is None
    assert row["classes"]["CDS"]["reference"] is None


def test_missing_reference_gff_is_not_available(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()], write_gff_file=False)
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "falta el GFF de referencia" in row["reason"]
    assert row["classes"]["CDS"]["system"] == 3
    assert row["classes"]["CDS"]["reference"] is None


def test_missing_annotation_path_is_not_available(tmp_path: Path) -> None:
    record = demo_record()
    record.pop("annotation_path")
    root = setup_root(tmp_path, [record])
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "no declara annotation_path" in row["reason"]


def test_invalid_summary_is_not_available(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()])
    summary = root / "results" / "demo" / "demo.summary.json"
    summary.write_text("{not json", encoding="utf-8")
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "resumen ilegible" in row["reason"]


def test_gff_without_target_classes_is_not_available(tmp_path: Path) -> None:
    root = setup_root(
        tmp_path,
        [demo_record()],
        gff_rows=[("gene", 100, 400, "ID=gene1")],
    )
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    row = result["rows"][0]
    assert row["state"] == "not_available"
    assert "no contiene CDS/tRNA/rRNA" in row["reason"]


# ---------------------------------------------------------------------------
# Selección del manifiesto y CLI
# ---------------------------------------------------------------------------
def test_selection_skips_downloaded_and_missing_summary(tmp_path: Path) -> None:
    records = [
        demo_record(),
        demo_record(id="available_no_summary", status="available", system_summary=""),
        demo_record(id="downloaded_with_summary", status="downloaded"),
        demo_record(id="annotated_no_summary", system_summary=""),
    ]
    root = setup_root(tmp_path, records)
    result = rp.build_report(root / "data" / "reference_genomes.json", root)

    assert [row["genome"] for row in result["rows"]] == ["demo"]
    assert result["summary"]["selected"] == 1
    assert result["summary"]["skipped"] == 3
    assert {item["genome"] for item in result["eligible_skipped"]} == {
        "available_no_summary",
        "annotated_no_summary",
    }
    assert any(
        item["genome"] == "downloaded_with_summary"
        and "fuera de annotated/available" in item["reason"]
        for item in result["skipped"]
    )


def test_cli_writes_tsv_and_markdown(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()])
    code = rp.main(["--root", str(root)])
    assert code == 0

    tsv_path = root / "reports" / "reference_parity.tsv"
    md_path = root / "reports" / "reference_parity.md"
    assert tsv_path.is_file()
    assert md_path.is_file()

    rows = read_tsv(tsv_path)
    assert len(rows) == 1
    assert rows[0]["genome"] == "demo"
    assert rows[0]["state"] == "ok"
    assert rows[0]["cds_system"] == "3"
    assert rows[0]["cds_reference"] == "3"
    assert rows[0]["cds_delta"] == "0"
    assert rows[0]["cds_pct"] == "100.00"
    assert rows[0]["trna_pct"] == "100.00"

    text = md_path.read_text(encoding="utf-8")
    assert "# Paridad de anotación: sistema vs GFF de referencia" in text
    assert "| `demo` | ok |" in text
    assert "conteos idénticos en CDS, tRNA y rRNA" in text


def test_cli_missing_manifest_returns_2(tmp_path: Path, capsys) -> None:
    root = tmp_path / "empty"
    root.mkdir()
    code = rp.main(["--root", str(root)])
    assert code == 2
    assert "ERROR:" in capsys.readouterr().err
    assert not (root / "reports" / "reference_parity.tsv").exists()


def test_tsv_leaves_blanks_for_not_available(tmp_path: Path) -> None:
    root = setup_root(tmp_path, [demo_record()], write_summary_file=False)
    result = rp.build_report(root / "data" / "reference_genomes.json", root)
    tsv_path = root / "reports" / "reference_parity.tsv"
    md_path = root / "reports" / "reference_parity.md"
    rp.write_reports(result, tsv_path, md_path)

    rows = read_tsv(tsv_path)
    assert len(rows) == 1
    assert rows[0]["state"] == "not_available"
    for kind in ("cds", "trna", "rrna"):
        assert rows[0][f"{kind}_system"] == ""
        assert rows[0][f"{kind}_reference"] == ""
        assert rows[0][f"{kind}_delta"] == ""
        assert rows[0][f"{kind}_pct"] == ""
    assert "falta el resumen del sistema" in rows[0]["reason"]
    assert "not_available" in md_path.read_text(encoding="utf-8")
