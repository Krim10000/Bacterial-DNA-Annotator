"""Contrato de la auditoría de cobertura de transferencia K-12.

Los TSV son sintéticos: no se leen los 400 genomas ni se descarga nada.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.complexes.transfer_coverage_audit import (
    CATALOG_AVAILABLE,
    NOT_AVAILABLE,
    TRANSFER_ROLE,
    TSV_COLUMNS,
    build_audit,
    discover_genomes,
    main,
    render_markdown,
    write_tsv,
)

TRANSFER_COLUMNS = [
    "complex_id", "complex_name", "n_subunits", "n_mapped", "frac_mapped",
    "identity_mean", "score", "genes_k12", "genes_genoma", "observed_n",
    "expected_n", "stoich_status", "mass_kda", "phylo_jacc", "ecs",
    "evidence_class", "eco", "limitation",
]
MODEL_COLUMNS = [
    "complex_id", "name", "contig", "start", "end", "subunits", "expected",
    "missing", "status", "evidence", "confidence", "source", "ecs", "flag",
    "limitation",
]
GOLD_COLUMNS = [
    "complex_id", "complex_name", "subunits", "n_genes", "locus_tags",
    "gene_names", "evidence_code", "evidence_class", "properties", "source",
]


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _transfer_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in TRANSFER_COLUMNS}
    row.update({
        "complex_id": "CPX-1",
        "complex_name": "Test complex",
        "n_subunits": "2",
        "n_mapped": "2",
        "frac_mapped": "1.0",
        "identity_mean": "80.0",
        "genes_k12": "b0001,b0002",
        "genes_genoma": "g1,g2",
        "observed_n": "2",
        "stoich_status": "unknown",
        "evidence_class": "documentado_inferido",
        "eco": "ECO:0005544",
        "limitation": "RBH sin cobertura de alineamiento en los TSV actuales.",
    })
    row.update(overrides)
    return row


def _model_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in MODEL_COLUMNS}
    row.update({
        "complex_id": "GCF_AAA:cx1",
        "name": "Candidato",
        "subunits": "g1:g1;g2:g2",
        "status": "complex_candidate",
        "evidence": "modelo K-12",
        "confidence": "alta",
        "source": "modelo K-12 (build_candidates.py)",
        "flag": "yellow",
        "limitation": (
            "Sospecha por modelo genómico K-12; no demuestra interacción "
            "física ni estequiometría. No cierra brechas."
        ),
    })
    row.update(overrides)
    return row


def make_tree(tmp_path: Path) -> Path:
    """Dos genomas: uno sin catálogo y otro con catálogo independiente."""
    root = tmp_path
    cx = root / "data" / "complexes"

    # Catálogo K-12: solo mapea b-numbers a nombres de gen para las firmas.
    _tsv(cx / "gold_k12.tsv", GOLD_COLUMNS, [
        {"complex_id": "CPX-1", "complex_name": "ATP synthase",
         "locus_tags": "b0001,b0002,b0003,b0004",
         "gene_names": "atpA,atpB,atpC,atpD"},
        {"complex_id": "CPX-2", "complex_name": "ATP synthase partial",
         "locus_tags": "b0010,b0011,b0012",
         "gene_names": "atpA,atpB,atpC"},
        {"complex_id": "CPX-3", "complex_name": "Other",
         "locus_tags": "b0020,b0021", "gene_names": "fooA,fooB"},
    ])

    # GCF_AAA: transferencia + modelo, sin catálogo independiente.
    _tsv(cx / "transfer_GCF_AAA.tsv", TRANSFER_COLUMNS, [
        _transfer_row(
            complex_id="CPX-1", complex_name="ATP synthase",
            n_subunits="4", n_mapped="4", frac_mapped="1.0",
            identity_mean="90.0", genes_k12="b0001,b0002,b0003,b0004",
            genes_genoma="g1,g2,g3,g4", observed_n="4", expected_n="4",
            stoich_status="complete", ecs="7.1.2.2",
        ),
        _transfer_row(
            complex_id="CPX-2", complex_name="ATP synthase partial",
            n_subunits="4", n_mapped="3", frac_mapped="0.75",
            identity_mean="80.0", genes_k12="b0010,b0011,b0012",
            genes_genoma="g5,g6,g7", observed_n="3", expected_n="4",
            stoich_status="incomplete", ecs="",
        ),
        _transfer_row(
            complex_id="CPX-3", complex_name="Other",
            n_subunits="2", n_mapped="2", frac_mapped="1.0",
            identity_mean="70.0", genes_k12="b0020,b0021",
            genes_genoma="g8,g9", observed_n="2", expected_n="",
            stoich_status="unknown", ecs="",
        ),
    ])
    _tsv(cx / "model_GCF_AAA.complexes.tsv", MODEL_COLUMNS, [
        _model_row(complex_id="GCF_AAA:cx1"),
        _model_row(
            complex_id="GCF_AAA:cx2",
            subunits="g3:g3;g4:g4;g5:g5",
            status="complex_complete",
        ),
        _model_row(
            complex_id="GCF_AAA:cx3", subunits="g6:g6",
            status="complex_incomplete",
        ),
    ])

    # GCF_BBB: transferencia + modelo + catálogo independiente.
    _tsv(cx / "transfer_GCF_BBB.tsv", TRANSFER_COLUMNS, [
        _transfer_row(
            complex_id="CPX-9", complex_name="Other",
            n_subunits="3", n_mapped="2", frac_mapped="0.6667",
            identity_mean="60.0", genes_k12="b0020,b0021",
            genes_genoma="g1,g2", observed_n="2", expected_n="3",
            stoich_status="incomplete",
        ),
    ])
    _tsv(cx / "model_GCF_BBB.complexes.tsv", MODEL_COLUMNS, [
        _model_row(complex_id="GCF_BBB:cx1"),
    ])
    _write(
        cx / "catalog_GCF_BBB.tsv",
        "complex_id\tcomplex_name\tlocus_tags\n"
        "CPX-B1\tCatálogo B\tg1,g2\n"
        "CPX-B2\tCatálogo B2\tg3,g4\n",
    )
    return root


def _by_genome(rows: list[dict]) -> dict[str, dict]:
    return {row["genome"]: row for row in rows}


def test_audit_per_genome_counts_and_states(tmp_path):
    root = make_tree(tmp_path)
    rows = build_audit(root)
    by = _by_genome(rows)

    assert discover_genomes(root / "data" / "complexes") == [
        "GCF_AAA", "GCF_BBB",
    ]
    assert [row["genome"] for row in rows] == ["GCF_AAA", "GCF_BBB"]

    a = by["GCF_AAA"]
    assert a["transfer_role"] == TRANSFER_ROLE
    assert a["transfer_complexes"] == 3
    assert a["transfer_suspected"] == 2
    assert a["transfer_candidate"] == 1
    assert a["transfer_stoich_complete"] == 1
    assert a["transfer_stoich_incomplete"] == 1
    assert a["transfer_stoich_unknown"] == 1
    assert a["transfer_full_subunits"] == 2
    assert a["transfer_with_ec"] == 1
    assert a["subunits_expected"] == 10
    assert a["subunits_mapped"] == 9
    assert a["subunit_coverage"] == pytest.approx(0.9)
    assert a["mean_frac_mapped"] == pytest.approx((1.0 + 0.75 + 1.0) / 3)
    assert a["mean_identity"] == pytest.approx(80.0)
    assert a["signature_status"] == "evaluated"
    assert a["signature_matches"] == 2
    assert a["signature_rows"] == 2
    assert a["signature_complete"] == 1
    assert a["signature_incomplete"] == 1
    assert a["signature_ids"] == "ATP_SYNTHASE"
    assert a["model_complexes"] == 3
    assert a["model_candidates"] == 1
    assert a["model_complete"] == 1
    assert a["model_incomplete"] == 1
    assert a["model_subunits"] == 2 + 3 + 1
    assert a["independent_catalog"] == NOT_AVAILABLE
    assert a["catalog_status"] == NOT_AVAILABLE
    assert a["catalog_complexes"] == 0
    assert "sin catálogo ComplexPortal independiente" in a["limitations"]
    assert "no es validación" in a["limitations"]

    b = by["GCF_BBB"]
    assert b["catalog_status"] == CATALOG_AVAILABLE
    assert b["independent_catalog"] == "data/complexes/catalog_GCF_BBB.tsv"
    assert b["catalog_complexes"] == 2
    assert b["transfer_complexes"] == 1
    assert b["transfer_candidate"] == 1
    assert b["subunit_coverage"] == pytest.approx(2 / 3)
    assert b["model_complexes"] == 1
    assert b["signature_matches"] == 0
    assert "catálogo independiente disponible" in b["limitations"]
    assert "no se usa como validación" in b["limitations"]


def test_k12_gold_is_not_independent_catalog(tmp_path):
    root = make_tree(tmp_path)
    a = _by_genome(build_audit(root))["GCF_AAA"]
    assert (root / "data" / "complexes" / "gold_k12.tsv").is_file()
    assert a["independent_catalog"] == NOT_AVAILABLE
    assert a["catalog_status"] == NOT_AVAILABLE
    assert "gold_k12" not in a["independent_catalog"]


def test_empty_catalog_is_not_available(tmp_path):
    root = make_tree(tmp_path)
    _write(
        root / "data" / "complexes" / "catalog_GCF_BBB.tsv",
        "complex_id\tcomplex_name\tlocus_tags\n",
    )
    b = _by_genome(build_audit(root))["GCF_BBB"]
    assert b["catalog_status"] == NOT_AVAILABLE
    assert b["independent_catalog"] == NOT_AVAILABLE
    assert b["catalog_complexes"] == 0


def test_taxid_named_catalog_is_found(tmp_path):
    root = tmp_path
    cx = root / "data" / "complexes"
    _tsv(cx / "transfer_GCF_CCC.tsv", TRANSFER_COLUMNS, [
        _transfer_row(complex_id="CPX-C1", genes_k12="b0020,b0021"),
    ])
    _write(
        root / "data" / "taxonomy" / "full" / "GCF_CCC.json",
        '{"genome": "GCF_CCC", "ncbi_taxid": 555}',
    )
    _write(
        cx / "gold_555.tsv",
        "complex_id\tcomplex_name\tlocus_tags\nCPX-T1\tCatálogo taxid\tg1,g2\n",
    )

    row = build_audit(root)[0]
    assert row["catalog_status"] == CATALOG_AVAILABLE
    assert row["catalog_complexes"] == 1
    assert row["independent_catalog"] == "data/complexes/gold_555.tsv"


def test_tsv_and_markdown_contract(tmp_path):
    root = make_tree(tmp_path)
    rows = build_audit(root)
    out = tmp_path / "out" / "coverage.tsv"
    write_tsv(rows, out)

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(TSV_COLUMNS)
    assert len(lines) == len(rows) + 1
    assert all(len(line.split("\t")) == len(TSV_COLUMNS) for line in lines)

    md = render_markdown(rows, "2026-01-02T03:04:05Z")
    assert "# Cobertura de la transferencia de complejos K-12" in md
    assert "Generado: `2026-01-02T03:04:05Z`" in md
    assert "| genomas auditados | 2 |" in md
    assert "| con catálogo ComplexPortal independiente | 1 |" in md
    assert "| sin catálogo independiente (`not_available`) | 1 |" in md
    assert "| complejos transferidos | 4 |" in md
    assert "| transferidos con `frac_mapped >= 0.8` (sospechosos) | 2 |" in md
    assert "| candidatos del modelo K-12 | 4 |" in md
    assert "| coincidencias con firmas curadas | 2 |" in md
    assert "| firmas completas / incompletas | 1 / 1 |" in md
    assert "## Detalle por genoma" in md
    assert "| `GCF_AAA` | `not_available` | 3 | 2 | 0.9000 | 1 / 1 | 3 | 1 / 1 / 1 |" in md
    assert "data/complexes/catalog_GCF_BBB.tsv" in md


def test_main_writes_reports_and_does_not_modify_inputs(tmp_path):
    root = make_tree(tmp_path)
    inputs = sorted((root / "data").rglob("*.tsv"))
    before = {path: path.read_bytes() for path in inputs}

    out_tsv = tmp_path / "out" / "coverage.tsv"
    out_md = tmp_path / "out" / "coverage.md"
    assert main([
        "--root", str(root),
        "--out-tsv", str(out_tsv),
        "--out-md", str(out_md),
        "--generated-at", "2026-01-02T03:04:05Z",
    ]) == 0

    assert out_tsv.is_file()
    assert out_md.is_file()
    assert out_tsv.read_text(encoding="utf-8").count("\n") == 3
    assert "Generado: `2026-01-02T03:04:05Z`" in out_md.read_text(
        encoding="utf-8"
    )
    assert {path: path.read_bytes() for path in inputs} == before


def test_model_only_genome_and_empty_tree(tmp_path):
    root = make_tree(tmp_path)
    (root / "data" / "complexes" / "transfer_GCF_BBB.tsv").unlink()
    b = _by_genome(build_audit(root))["GCF_BBB"]
    assert b["transfer_complexes"] == 0
    assert b["model_complexes"] == 1
    assert "sin complejos transferidos" in b["limitations"]

    empty = tmp_path / "empty"
    assert build_audit(empty) == []
    md = render_markdown([], "2026-01-02T03:04:05Z")
    assert "| genomas auditados | 0 |" in md
    assert "| sin catálogo independiente (`not_available`) | 0 |" in md


def test_blank_values_do_not_crash(tmp_path):
    root = tmp_path
    cx = root / "data" / "complexes"
    _tsv(cx / "transfer_GCF_DDD.tsv", TRANSFER_COLUMNS, [
        _transfer_row(
            complex_id="CPX-D1", n_subunits="", n_mapped="",
            frac_mapped="", identity_mean="", genes_k12="",
            observed_n="", expected_n="", stoich_status="",
        ),
    ])
    _tsv(cx / "model_GCF_DDD.complexes.tsv", MODEL_COLUMNS, [
        _model_row(complex_id="GCF_DDD:cx1", subunits="", status=""),
    ])

    row = build_audit(root)[0]
    assert row["transfer_complexes"] == 1
    assert row["transfer_candidate"] == 1
    assert row["transfer_stoich_unknown"] == 1
    assert row["subunits_expected"] == 0
    assert row["subunit_coverage"] == 0.0
    assert row["model_other"] == 1
    assert row["model_subunits"] == 0
    # No se inventa un catálogo ni un estado de validación.
    assert row["catalog_status"] == NOT_AVAILABLE
    assert row["transfer_role"] == TRANSFER_ROLE
