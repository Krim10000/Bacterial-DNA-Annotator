"""Tests del golden set TSS-sigma (``data/sigma/gold/``).

Cubren:

- contrato del esquema y cabeceras;
- validador limpio y CLI;
- conteos congelados del snapshot;
- procedencia K-12 contra ``PromoterSet.tsv`` y ``gt_promoters.tsv``;
- procedencia *B. subtilis* contra ``gt_tss.tsv``;
- deteccion de filas invalidas (coordenada, hebra, sigma, PMID, duplicados).

Los tests de procedencia se saltan si no estan las fuentes primarias.
"""
from __future__ import annotations

import csv
import importlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

V = importlib.import_module("scripts.sigma.40_validate_gold")

GOLD = ROOT / "data" / "sigma" / "gold"
SCHEMA_PATH = GOLD / "schema.json"
K12_PATH = GOLD / "ecoli_k12" / "tss_sigma.tsv"
BSU_PATH = GOLD / "subtilis" / "tss_sigma.tsv"
K12_SOURCE = ROOT / "data" / "regulondb" / "PromoterSet.tsv"
K12_GT = ROOT / "data" / "regulondb" / "gt" / "gt_promoters.tsv"
BSU_SOURCE = ROOT / "data" / "dbtbs" / "gt_tss.tsv"

COLUMNS = [
    "species", "genome_id", "accession", "contig", "pos", "strand",
    "sigma", "gene", "evidence", "pmid", "source", "confidence", "notes",
]
SIGMAS = {"sigma24", "sigma28", "sigma32", "sigma38", "sigma54", "sigma70"}
BSU_MAP = {"A": "sigma70", "B": "sigma38", "D": "sigma28",
           "I": "sigma32", "L": "sigma54"}
PMID_SPLIT = re.compile(r"[;,]\s*")

pytestmark = pytest.mark.skipif(
    not GOLD.exists(),
    reason="golden set de sigma no incluido en esta instalación",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def clean(value: str) -> str:
    return " ".join(value.replace("\t", " ").replace("\r", " ").replace("\n", " ").split())


def norm_pmids(value: str) -> str:
    result: list[str] = []
    for part in PMID_SPLIT.split(value or ""):
        part = part.strip()
        if part and part not in result:
            result.append(part)
    return ";".join(result)


def extract_pubmed(evidence: str) -> str:
    result: list[str] = []
    for match in re.finditer(r"PubMed:\s*([0-9][0-9,\s;]*)", evidence):
        for token in re.findall(r"\d+", match.group(1)):
            if token not in result:
                result.append(token)
    return ";".join(result)


@pytest.fixture(scope="module")
def schema() -> dict:
    return V.load_schema(SCHEMA_PATH)


@pytest.fixture(scope="module")
def report(schema: dict) -> dict:
    return V.validate_gold(schema, GOLD)


# ---------------------------------------------------------------------------
# Esquema y archivos
# ---------------------------------------------------------------------------
def test_schema_contract(schema: dict) -> None:
    assert schema["name"] == "sigma-gold"
    assert schema["version"] == "1.0.0"
    assert [col["name"] for col in schema["columns"]] == COLUMNS
    assert schema["primary_key"] == [
        "genome_id", "contig", "pos", "strand", "sigma",
    ]
    assert set(schema["sigma_allowlist"]) == SIGMAS
    assert schema["genomes"]["ecoli_k12"]["contig"] == "NC_000913.3"
    assert schema["genomes"]["subtilis"]["contig"] == "basu168"
    assert schema["contains"]["mechanism_labels"] is False
    assert schema["contains"]["negative_set"] is False
    # No hay columna de caso/mecanismo.
    assert "case" not in COLUMNS and "mechanism" not in COLUMNS


def test_gold_files_have_exact_header() -> None:
    for path in (K12_PATH, BSU_PATH):
        assert path.is_file(), path
        header, rows = V.read_tsv(path)
        assert header == COLUMNS
        assert rows
        assert all(len(row) == len(COLUMNS) for _, row in rows)
        assert all(any(cell.strip() for cell in row) for _, row in rows)


def test_validator_accepts_gold(report: dict) -> None:
    assert report["ok"] is True
    assert report["errors"] == []
    assert report["warnings"] == []
    assert report["summary"]["rows"] == 5540
    assert report["summary"]["files"] == 2
    assert report["summary"]["genomes"] == {"ecoli_k12": 3314, "subtilis": 2226}
    assert report["summary"]["duplicates"] == 0


def test_frozen_counts(report: dict) -> None:
    assert report["summary"]["sigma"] == {
        "sigma24": 522,
        "sigma28": 214,
        "sigma32": 325,
        "sigma38": 423,
        "sigma54": 110,
        "sigma70": 3946,
    }
    assert report["summary"]["confidence"] == {
        "C": 73,
        "S": 928,
        "W": 2307,
        "curated": 2232,
    }
    assert report["summary"]["source"] == {
        "DBTBS": 2226,
        "RegulonDB": 3314,
    }
    assert report["summary"]["pmid"] == {"with": 3743, "total": 5540}
    assert report["summary"]["gene"] == {"with": 1545, "total": 5540}


# ---------------------------------------------------------------------------
# Invariantes de datos
# ---------------------------------------------------------------------------
def test_primary_key_is_unique() -> None:
    rows = load_rows(K12_PATH) + load_rows(BSU_PATH)
    keys = [
        (row["genome_id"], row["contig"], row["pos"], row["strand"], row["sigma"])
        for row in rows
    ]
    assert len(keys) == 5540
    assert len(keys) == len(set(keys))


def test_coordinates_strands_sigmas_and_sources() -> None:
    for path in (K12_PATH, BSU_PATH):
        for row in load_rows(path):
            assert int(row["pos"]) > 0
            assert row["strand"] in {"+", "-"}
            assert row["sigma"] in SIGMAS
            assert row["sigma"] not in {"sigma19", ""}
            assert row["confidence"] in {"S", "W", "C", "curated"}
            assert row["source"] in {"RegulonDB", "DBTBS"}
            assert row["evidence"].strip()
            assert row["genome_id"] in {"ecoli_k12", "subtilis"}


def test_pmid_format_and_coverage() -> None:
    k12 = load_rows(K12_PATH)
    bsu = load_rows(BSU_PATH)
    for row in k12 + bsu:
        if not row["pmid"]:
            continue
        parts = [part.strip() for part in PMID_SPLIT.split(row["pmid"]) if part.strip()]
        assert parts
        assert all(part.isdigit() and int(part) > 0 for part in parts)
        assert len(parts) == len(set(parts))
    assert sum(1 for row in k12 if row["pmid"]) == 3314
    assert sum(1 for row in bsu if row["pmid"]) == 429


def test_gene_is_honest_about_missing_values() -> None:
    k12 = load_rows(K12_PATH)
    bsu = load_rows(BSU_PATH)
    assert sum(1 for row in k12 if row["gene"]) == 1545
    assert all(row["gene"] == "" for row in bsu)
    # No se usan placeholders como NA/None para rellenar.
    assert not any(row["gene"] in {"NA", "None", "-"} for row in k12 + bsu)


def test_k12_multi_sigma_positions_are_documented() -> None:
    positions = Counter(
        (row["pos"], row["strand"]) for row in load_rows(K12_PATH)
    )
    assert sum(1 for count in positions.values() if count > 1) == 124
    assert sum(count - 1 for count in positions.values() if count > 1) == 125


# ---------------------------------------------------------------------------
# Procedencia K-12
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (K12_SOURCE.is_file() and K12_GT.is_file()),
    reason="sin fuentes de RegulonDB en este checkout",
)
def test_k12_provenance_matches_sources() -> None:
    source_rows = load_rows(K12_SOURCE)
    gt_rows = {row["promoter_id"]: row for row in load_rows(K12_GT)}
    expected: dict[tuple[int, str, str], dict] = {}

    for row in source_rows:
        sigma = row["5)sigmaFactor"].strip()
        if sigma not in SIGMAS:
            continue
        try:
            pos = int(row["4)posTSS"])
        except (TypeError, ValueError):
            continue
        if pos <= 0:
            continue
        strand = {"forward": "+", "reverse": "-"}.get(
            row["3)strand"].strip().lower()
        )
        if strand is None:
            continue
        key = (pos, strand, sigma)
        pid = row["1)id"].strip()
        if key in expected:
            expected[key]["also"].append(pid)
            continue
        confidence = row["15)confidenceLevel"].strip()
        expected[key] = {
            "promoter_id": pid,
            "gene": clean(row["7)firstGeneName"]),
            "evidence": clean(row["13)promoterEvidence"]),
            "pmid": norm_pmids(row["16)pmids"]),
            "confidence": confidence if confidence in {"S", "W", "C"} else "curated",
            "also": [],
        }

    gold = {
        (int(row["pos"]), row["strand"], row["sigma"]): row
        for row in load_rows(K12_PATH)
    }
    assert len(expected) == 3314
    assert set(gold) == set(expected)

    for key, exp in expected.items():
        row = gold[key]
        assert row["species"] == "Escherichia coli K-12 MG1655"
        assert row["genome_id"] == "ecoli_k12"
        assert row["accession"] == "GCF_000005845.2"
        assert row["contig"] == "NC_000913.3"
        assert row["source"] == "RegulonDB"
        assert row["gene"] == exp["gene"]
        assert row["evidence"] == exp["evidence"]
        assert row["pmid"] == exp["pmid"]
        assert row["confidence"] == exp["confidence"]
        assert f"regulondb_id={exp['promoter_id']}" in row["notes"]

        gt = gt_rows[exp["promoter_id"]]
        assert int(gt["tss"]) == key[0]
        assert gt["strand"].strip() == key[1]
        assert gt["sigma_factor"].strip() == key[2]
        assert gt["first_gene"].strip() == exp["gene"]

    duplicate_key = (3946780, "+", "sigma28")
    assert duplicate_key in gold
    notes = gold[duplicate_key]["notes"]
    assert "also_regulondb_id=RDBECOLIPMC02191" in notes
    assert "also_name=trpTp2" in notes


# ---------------------------------------------------------------------------
# Procedencia B. subtilis
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not BSU_SOURCE.is_file(),
    reason="sin gt_tss.tsv de DBTBS en este checkout",
)
def test_subtilis_provenance_matches_source() -> None:
    source_rows = load_rows(BSU_SOURCE)
    expected: dict[tuple[int, str, str], dict] = {}

    for row in source_rows:
        label = row["sigma"].strip()
        if label not in BSU_MAP:
            continue
        pos = int(row["pos"])
        strand = row["strand"].strip()
        if pos <= 0 or strand not in {"+", "-"}:
            continue
        key = (pos, strand, BSU_MAP[label])
        assert key not in expected, f"duplicado inesperado en DBTBS: {key}"
        expected[key] = {
            "label": label,
            "evidence": clean(row["evidence"]),
            "pmid": extract_pubmed(row["evidence"]),
        }

    gold = {
        (int(row["pos"]), row["strand"], row["sigma"]): row
        for row in load_rows(BSU_PATH)
    }
    assert len(expected) == 2226
    assert set(gold) == set(expected)

    for key, exp in expected.items():
        row = gold[key]
        assert row["species"] == "Bacillus subtilis 168"
        assert row["genome_id"] == "subtilis"
        assert row["accession"] == "GCF_000009045.1"
        assert row["contig"] == "basu168"
        assert row["source"] == "DBTBS"
        assert row["gene"] == ""
        assert row["evidence"] == exp["evidence"]
        assert row["pmid"] == exp["pmid"]
        assert row["confidence"] == "curated"
        assert f"dbtbs_label={exp['label']}" in row["notes"]

    # Etiquetas excluidas: no deben aparecer como sigma normalizado.
    excluded = {"?", "E", "F", "G", "H", "K", "W", "X", "M", "YlaC"}
    source_labels = {row["sigma"].strip() for row in source_rows}
    assert excluded <= source_labels
    assert not (set(gold_key[2] for gold_key in gold) & excluded)


# ---------------------------------------------------------------------------
# Validador: casos invalidos y CLI
# ---------------------------------------------------------------------------
def test_validator_catches_bad_rows(tmp_path: Path, schema: dict) -> None:
    bad_path = tmp_path / "bad.tsv"
    good = [
        "Escherichia coli K-12 MG1655", "ecoli_k12", "GCF_000005845.2",
        "NC_000913.3", "100", "+", "sigma70", "thrL", "COMP-AINF:W",
        "14529615", "RegulonDB", "S", "ok",
    ]
    broken = [
        "Escherichia coli K-12 MG1655", "ecoli_k12", "GCF_BAD",
        "NC_999999.1", "0", "*", "sigma19", "", "", "abc", "Otra", "high",
        "rota",
    ]
    with bad_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerow(broken)
        writer.writerow(good)
        writer.writerow(good)  # duplicado de clave
        writer.writerow(good[:-1])  # 12 columnas
        writer.writerow([])  # linea vacia

    result = V.validate_file(bad_path, schema, "ecoli_k12")
    text = "\n".join(result["errors"])

    assert result["duplicates"] == 1
    for token in (
        "pos", "strand", "sigma", "contig", "source", "confidence",
        "pmid", "evidence", "clave duplicada", "columnas", "linea vacia",
        "accession",
    ):
        assert token in text, token


def test_cli_json_mode() -> None:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "sigma" / "40_validate_gold.py"),
         "--json"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert payload["ok"] is True
    assert payload["errors"] == []
    assert payload["summary"]["rows"] == 5540


def test_cli_text_mode() -> None:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "sigma" / "40_validate_gold.py")],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "OK: gold valido" in proc.stdout
    assert "filas: 5540" in proc.stdout


# ---------------------------------------------------------------------------
# Documentacion
# ---------------------------------------------------------------------------
def test_docs_explain_gold_limits_and_sigma_cases() -> None:
    docs_path = ROOT / "docs" / "SIGMA_GOLD.md"
    assert docs_path.is_file()
    text = docs_path.read_text(encoding="utf-8")
    for token in (
        "SIGMA_CASES", "SIGMA_BENCHMARK", "5.540", "3.314", "2.226",
        "schema.json", "40_validate_gold.py", "test_sigma_gold.py",
        "Qué es el gold", "Qué **no** es el gold", "Límites",
    ):
        assert token in text, token
