"""Tests de la capa de taxonomía local.

No usan red, GPU ni datos reales: construyen JSON sintéticos en ``tmp_path``
y comprueban el contrato de columnas, la ausencia de invención de especie y
el estado ``not_available`` cuando falta el genoma o el JSON es inválido.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from webapp.taxonomy_layer import (
    FIELDS,
    METHOD_GTDB,
    METHOD_NCBI,
    METHOD_SPECIES,
    RANK_ORDER,
    available,
    inventory,
    normalize_taxonomy,
    run,
)

GTDB_FULL = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
    "s__Escherichia coli"
)


def _record(**overrides) -> dict:
    data = {
        "genome": "GCF_TEST.1",
        "status": "assigned",
        "confidence": "high",
        "species": "Escherichia coli",
        "species_ncbi": "Escherichia coli",
        "species_gtdb": "Escherichia coli",
        "closest_ref": "GCF_REF.1",
        "ani": 98.7,
        "gtdb_lineage": GTDB_FULL,
        "ncbi_lineage": GTDB_FULL,
        "evidence": {"conflicts": []},
    }
    data.update(overrides)
    return data


def _data_root(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    (data / "taxonomy" / "full").mkdir(parents=True)
    return data


def _write_record(data_root: Path, genome: str, data: dict) -> Path:
    path = data_root / "taxonomy" / "full" / f"{genome}.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_available_without_taxonomy_dir(tmp_path: Path):
    empty = tmp_path / "empty"
    empty.mkdir()
    ok, why = available(empty)
    assert ok is False
    assert "taxonom" in why.lower()


def test_available_detects_local_maps(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_TEST.1", _record())
    (data / "taxonomy" / "ani_hits.tsv").write_text(
        "genome\tref\tani\n", encoding="utf-8")
    (data / "selected.tsv").write_text(
        "assembly_accession\nGCF_TEST.1\n", encoding="utf-8")

    ok, why = available(data)
    assert ok is True
    assert why == ""

    info = inventory(data)
    assert info["available"] is True
    assert info["n_full"] == 1
    assert info["n_taxonomy"] == 1
    assert "taxonomy/ani_hits.tsv" in [m["path"] for m in info["ani_maps"]]
    assert "selected.tsv" in [m["path"] for m in info["genome_maps"]]


def test_available_accepts_minimal_lineage_record(tmp_path: Path):
    """`available` no exige el campo ``genome`` si el JSON tiene linaje."""
    data = _data_root(tmp_path)
    path = data / "taxonomy" / "full" / "GCF_MIN.1.json"
    path.write_text(json.dumps({"gtdb_lineage": GTDB_FULL}), encoding="utf-8")

    ok, why = available(data)
    assert ok is True
    assert why == ""

    result = run("GCF_MIN.1", data, tmp_path / "out")
    assert result["status"] == "ok"
    assert result["n_rows"] == len(RANK_ORDER) - 1


def test_run_valid_json_writes_expected_fields(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_TEST.1", _record())
    outdir = tmp_path / "out"

    result = run("GCF_TEST.1", data, outdir)
    assert result["status"] == "ok"
    assert result["ok"] is True
    expected_ranks = [rank for rank in RANK_ORDER if rank != "strain"]
    assert result["n_rows"] == len(expected_ranks)
    assert result["taxonomy_status"] == "assigned"
    assert result["species"] == "Escherichia coli"
    assert result["ani"] == 98.7

    tsv = Path(result["tsv"])
    assert tsv == outdir / "GCF_TEST.1.taxonomy.tsv"
    rows = _read_tsv(tsv)
    assert list(rows[0]) == list(FIELDS)
    assert [row["rank"] for row in rows] == expected_ranks
    assert rows[0]["taxon"] == "Bacteria"
    assert rows[-1]["taxon"] == "Escherichia coli"
    assert {row["method"] for row in rows} == {METHOD_GTDB}
    assert {row["confidence"] for row in rows} == {"high"}
    assert {row["source"] for row in rows} == {
        "taxonomy/full/GCF_TEST.1.json"
    }


def test_missing_genome_is_not_available(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_OTHER.1", _record(genome="GCF_OTHER.1"))
    outdir = tmp_path / "out"

    result = run("GCF_MISSING.1", data, outdir)
    assert result["status"] == "not_available"
    assert result["ok"] is False
    assert result["rows"] == []
    assert result["n_rows"] == 0
    assert "GCF_MISSING.1" in result["reason"]
    assert not (outdir / "GCF_MISSING.1.taxonomy.tsv").exists()


def test_species_fallback_uses_declared_field_only(tmp_path: Path):
    data = _data_root(tmp_path)
    record = _record(
        confidence="medium",
        species_ncbi=None,
        gtdb_lineage="unknown",
        ncbi_lineage=(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__"
        ),
        species="Escherichia coli",
    )
    _write_record(data, "GCF_TEST.1", record)

    result = run("GCF_TEST.1", data, tmp_path / "out")
    assert result["status"] == "ok"
    rows = _read_tsv(Path(result["tsv"]))
    assert [row["method"] for row in rows[:-1]] == [METHOD_NCBI] * 6
    assert rows[-1]["rank"] == "species"
    assert rows[-1]["taxon"] == "Escherichia coli"
    assert rows[-1]["method"] == METHOD_SPECIES


def test_no_species_is_invented_when_json_lacks_it(tmp_path: Path):
    data = _data_root(tmp_path)
    record = _record(
        species=None,
        species_ncbi=None,
        species_gtdb=None,
        ncbi_lineage=None,
        gtdb_lineage=(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia"
        ),
    )
    _write_record(data, "GCF_TEST.1", record)

    result = run("GCF_TEST.1", data, tmp_path / "out")
    assert result["status"] == "ok"
    assert result["species"] is None
    rows = _read_tsv(Path(result["tsv"]))
    assert "species" not in [row["rank"] for row in rows]
    assert all(row["taxon"] for row in rows)


def test_invalid_json_is_not_available(tmp_path: Path):
    data = _data_root(tmp_path)
    path = data / "taxonomy" / "full" / "GCF_BAD.1.json"
    path.write_text("{not valid json", encoding="utf-8")
    outdir = tmp_path / "out"

    result = run("GCF_BAD.1", data, outdir)
    assert result["status"] == "not_available"
    assert "JSON" in result["reason"]
    assert result["rows"] == []
    assert not (outdir / "GCF_BAD.1.taxonomy.tsv").exists()


def test_normalize_rows_have_exactly_the_contract_fields():
    rows = normalize_taxonomy(_record(), source="taxonomy/full/GCF_TEST.1.json")
    assert rows
    assert all(set(row) == set(FIELDS) for row in rows)
    assert [row["rank"] for row in rows] == [rank for rank in RANK_ORDER
                                              if rank != "strain"]
