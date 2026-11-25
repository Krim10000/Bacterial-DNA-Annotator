"""Contrato del inventario ComplexPortal: local, determinista y sin inventar taxid."""
from __future__ import annotations

import json
from pathlib import Path

from scripts.complexes.inventory_complexportal_400 import (
    NOT_AVAILABLE,
    STATUS_CATALOG_MISSING,
    STATUS_RAW_ONLY,
    STATUS_READY,
    STATUS_TAXID_MISSING,
    STATUS_TAXONOMY_MISSING,
    TSV_COLUMNS,
    build_inventory,
    discover_genomes,
    main,
    render_markdown,
    write_tsv,
)

GENES_HEADER = (
    "genome\tcontig\tstart\tend\tstrand\tgene_id\tlocus_tag\tgene_name\tacc\n"
)


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _genes(root: Path, genome: str, acc: str | None = "P1") -> None:
    acc_value = acc or ""
    _write(
        root / "data" / "regulation" / genome / "genes.tsv",
        GENES_HEADER
        + f"{genome}\tctg\t1\t100\t1\t{genome}_1\t\tgeneA\t{acc_value}\n"
        + f"{genome}\tctg\t200\t300\t1\t{genome}_2\t\tgeneB\t{acc_value}\n",
    )


def _taxonomy(
    root: Path,
    genome: str,
    *,
    species: str,
    taxid: int | None = None,
    layer: str = "full",
) -> None:
    data = {"genome": genome, "species": species}
    if taxid is not None:
        data["ncbi_taxid"] = taxid
        data["evidence"] = {"ncbi": {"taxid": taxid}}
    _write(
        root / "data" / "taxonomy" / layer / f"{genome}.json",
        json.dumps(data),
    )


def _functional(root: Path, genome: str) -> None:
    _write(
        root / "data" / "annotations_ab_initio" / genome / "functional.tsv",
        f"{genome}_1\tsp|P1|TEST\t99\t100\t0\t0\t1\t100\t1\t100\t1e-9\t200\n",
    )


def _catalog(path: Path, taxid: int) -> None:
    _write(
        path,
        "#Complex ac\tRecommended name\tTaxonomy identifier\tMembers\n"
        f"CPX-1\tTest complex\t{taxid}\tP1(1)\n",
    )


def make_tree(tmp_path: Path, *, k12_raw: bool = True) -> Path:
    root = tmp_path
    # GCF_AAA: taxid + functional.tsv, sin catálogo.
    _genes(root, "GCF_AAA")
    _taxonomy(root, "GCF_AAA", species="Alpha beta", taxid=111)
    _functional(root, "GCF_AAA")

    # GCF_BBB: especie pero sin taxid.
    _genes(root, "GCF_BBB")
    _taxonomy(root, "GCF_BBB", species="Beta gamma")

    # GCF_CCC: taxid en evidence.ncbi, JSON en data/taxonomy/ (fallback).
    _genes(root, "GCF_CCC")
    _write(
        root / "data" / "taxonomy" / "GCF_CCC.json",
        json.dumps({
            "genome": "GCF_CCC",
            "species": "Gamma delta",
            "evidence": {"ncbi": {"taxid": 333}},
        }),
    )

    # GCF_DDD: sin JSON, sin ac, sin catálogo.
    _genes(root, "GCF_DDD", acc=None)

    # GCF_RAW: solo TSV crudo.
    _genes(root, "GCF_RAW")
    _taxonomy(root, "GCF_RAW", species="Raw species", taxid=444)
    _catalog(root / "data" / "complexes" / "raw" / "complexportal_444.tsv", 444)

    # GCF_GOLD: catálogo gold ya mapeado.
    _genes(root, "GCF_GOLD")
    _taxonomy(root, "GCF_GOLD", species="Gold species", taxid=555)
    _write(root / "data" / "complexes" / "gold_555.tsv",
           "complex_id\tcomplex_name\tlocus_tags\nCPX-2\tGold\tgeneA,geneB\n")

    # K-12: sin JSON, especie del GBFF, taxid del TSV crudo.
    _genes(root, "ecoli_k12")
    _write(
        root / "data" / "regulation" / "ecoli_k12" / "genes_meta.json",
        json.dumps({"source": {"gbff": "data/test_ecoli.gbff"}}),
    )
    _write(
        root / "data" / "test_ecoli.gbff",
        "LOCUS       TEST\n"
        "DEFINITION  Escherichia coli str. K-12 substr. MG1655.\n"
        "SOURCE      Escherichia coli str. K-12 substr. MG1655\n"
        "  ORGANISM  Escherichia coli str. K-12 substr. MG1655\n",
    )
    if k12_raw:
        _catalog(
            root / "data" / "complexes" / "raw" / "complexportal_83333.tsv",
            83333,
        )
    _write(root / "data" / "complexes" / "gold_k12.tsv",
           "complex_id\tcomplex_name\tlocus_tags\nCPX-3\tK12\tb0001,b0002\n")

    # Directorio oculto sin genes.tsv válido: no debe entrar.
    _write(root / "data" / "regulation" / ".tmp" / "genes.tsv", GENES_HEADER)
    return root


def _by_genome(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["genome"]: row for row in rows}


def test_build_inventory_contract_and_statuses(tmp_path):
    root = make_tree(tmp_path)
    rows = build_inventory(root)
    by = _by_genome(rows)

    assert discover_genomes(root) == sorted(by)
    assert ".tmp" not in by
    assert [row["genome"] for row in rows] == sorted(by)
    assert all(tuple(row) == TSV_COLUMNS for row in rows)

    assert by["GCF_AAA"]["status"] == STATUS_CATALOG_MISSING
    assert by["GCF_AAA"]["taxid"] == "111"
    assert by["GCF_AAA"]["local_catalog"] == NOT_AVAILABLE
    assert "functional.tsv" in by["GCF_AAA"]["reason"]

    assert by["GCF_BBB"]["status"] == STATUS_TAXID_MISSING
    assert by["GCF_BBB"]["taxid"] == NOT_AVAILABLE
    assert by["GCF_BBB"]["species"] == "Beta gamma"

    assert by["GCF_CCC"]["taxid"] == "333"
    assert by["GCF_CCC"]["status"] == STATUS_CATALOG_MISSING

    assert by["GCF_DDD"]["status"] == STATUS_TAXONOMY_MISSING
    assert by["GCF_DDD"]["species"] == NOT_AVAILABLE
    assert by["GCF_DDD"]["taxid"] == NOT_AVAILABLE
    assert "no disponible" in by["GCF_DDD"]["reason"]

    assert by["GCF_RAW"]["status"] == STATUS_RAW_ONLY
    assert by["GCF_RAW"]["local_catalog"] == (
        "data/complexes/raw/complexportal_444.tsv"
    )

    assert by["GCF_GOLD"]["status"] == STATUS_READY
    assert by["GCF_GOLD"]["local_catalog"] == "data/complexes/gold_555.tsv"

    k12 = by["ecoli_k12"]
    assert k12["status"] == STATUS_READY
    assert k12["taxid"] == "83333"
    assert k12["species"] == "Escherichia coli str. K-12 substr. MG1655"
    assert k12["local_catalog"] == "data/complexes/gold_k12.tsv"
    assert "taxid leído del catálogo" in k12["reason"]
    assert "genes.tsv acc" in k12["reason"]

    counts = {
        status: sum(1 for row in rows if row["status"] == status)
        for status in {
            STATUS_READY,
            STATUS_RAW_ONLY,
            STATUS_CATALOG_MISSING,
            STATUS_TAXID_MISSING,
            STATUS_TAXONOMY_MISSING,
        }
    }
    assert counts == {
        STATUS_READY: 2,
        STATUS_RAW_ONLY: 1,
        STATUS_CATALOG_MISSING: 2,
        STATUS_TAXID_MISSING: 1,
        STATUS_TAXONOMY_MISSING: 1,
    }


def test_k12_without_raw_catalog_does_not_invent_taxid(tmp_path):
    root = make_tree(tmp_path, k12_raw=False)
    by = _by_genome(build_inventory(root))
    k12 = by["ecoli_k12"]

    assert k12["taxid"] == NOT_AVAILABLE
    assert k12["status"] == STATUS_TAXID_MISSING
    assert k12["local_catalog"] == NOT_AVAILABLE
    assert "83333" not in " ".join(k12.values())


def test_tsv_and_markdown_contract(tmp_path):
    root = make_tree(tmp_path)
    rows = build_inventory(root)
    out = tmp_path / "reports" / "inventory.tsv"
    write_tsv(rows, out)

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(TSV_COLUMNS)
    assert len(lines) == len(rows) + 1
    assert all(len(line.split("\t")) == len(TSV_COLUMNS) for line in lines)

    md = render_markdown(rows, "2026-01-02T03:04:05Z")
    assert "# Inventario de catálogos ComplexPortal" in md
    assert "2026-01-02T03:04:05Z" in md
    assert "| genomas | 7 |" in md
    assert "| con catálogo local | 3 |" in md
    assert "| panel con catálogo local | 2 |" in md
    assert "## Inventario completo" in md
    assert "| genome | species | taxid | local_catalog | status | reason |" in md
    assert "`ecoli_k12`" in md
    assert "`catalog_missing`" in md


def test_main_writes_reports_with_fixed_timestamp(tmp_path):
    root = make_tree(tmp_path)
    out_tsv = tmp_path / "out" / "inventory.tsv"
    out_md = tmp_path / "out" / "inventory.md"

    assert main([
        "--root", str(root),
        "--out-tsv", str(out_tsv),
        "--out-md", str(out_md),
        "--generated-at", "2026-01-02T03:04:05Z",
    ]) == 0

    assert out_tsv.is_file()
    assert out_md.is_file()
    assert out_tsv.read_text(encoding="utf-8").count("\n") == 8
    assert "Generado: `2026-01-02T03:04:05Z`" in out_md.read_text(
        encoding="utf-8"
    )


def test_empty_tree_is_valid(tmp_path):
    assert build_inventory(tmp_path) == []
    md = render_markdown([], "2026-01-02T03:04:05Z")
    assert "No hay ningún catálogo ComplexPortal local." in md
    assert "| genomas | 0 |" in md
