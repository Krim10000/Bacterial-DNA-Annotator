"""Tests de métricas y secciones sintéticas de ``tools/benchmark_layers.py``.

No usan red ni los datos reales del proyecto: crean un proyecto mínimo en
``tmp_path`` para comprobar que las métricas y los ``not_evaluable`` son
reproducibles.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.benchmark_layers import (  # noqa: E402
    best_match,
    build_data,
    circular_distance,
    distance_summary,
    greedy_one_to_one,
    jaccard,
    render_markdown,
    set_coverage,
    summarize_abc,
    summarize_complexes,
    summarize_membrane,
    summarize_oric,
)


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Métricas puras
# ---------------------------------------------------------------------------
def test_jaccard_and_coverage():
    assert jaccard(set(), set()) == 0.0
    assert jaccard({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)
    out = set_coverage({"a", "b"}, {"b", "c"})
    assert out["n_shared"] == 1
    assert out["coverage_of_gold"] == 0.5
    assert out["fraction_of_pred_in_gold"] == 0.5
    assert out["jaccard"] == pytest.approx(1 / 3, abs=1e-3)


def test_circular_distance_wraps_and_uses_same_convention_as_orca():
    # Caso real: 229 y 1.072.919 en un contig de 1.072.950 pb.
    assert circular_distance(229, 1_072_919, 1_072_950) == 260
    assert circular_distance(1, 1_000, 1_000) == 1
    assert circular_distance(500, 500, 1_000) == 0
    with pytest.raises(ValueError):
        circular_distance(1, 2, 0)


def test_distance_summary_ignores_none_and_counts_thresholds():
    out = distance_summary([0, 1_000, 50_000, None, 2_000_000])
    assert out["n"] == 4
    assert out["min_bp"] == 0
    assert out["max_bp"] == 2_000_000
    assert out["within_1000bp"]["n"] == 2
    assert out["within_10000bp"]["n"] == 2
    assert out["within_50000bp"]["n"] == 3
    assert out["beyond_max_threshold"]["n"] == 1
    assert distance_summary([])["n"] == 0


def test_best_match_and_greedy_one_to_one_are_deterministic():
    refs = {
        "g1": {"a", "b", "c"},
        "g2": {"a", "b"},
    }
    match = best_match({"a", "b"}, refs)
    # g2 es subconjunto exacto de la query; gana con Jaccard 1.0.
    assert match["id"] == "g2"
    assert match["jaccard"] == 1.0

    # Empate exacto: el desempate es lexicográfico por id.
    tie_refs = {"g1": {"a", "b", "c"}, "g2": {"a", "b", "d"}}
    assert best_match({"a", "b"}, tie_refs)["id"] == "g1"

    preds = {"c1": {"a", "b", "c"}, "c2": {"a", "b"}}
    matches = greedy_one_to_one(preds, refs, threshold=0.5)
    assert {m["pred_id"] for m in matches} == {"c1", "c2"}
    assert {m["gold_id"] for m in matches} == {"g1", "g2"}
    # Ningún cluster se reutiliza.
    assert len({m["pred_id"] for m in matches}) == len(matches)


# ---------------------------------------------------------------------------
# A/B/C sintético
# ---------------------------------------------------------------------------
def _abc_fixture(tmp_path: Path) -> Path:
    header = [
        "protein_id", "path", "label", "evidence_strength", "evidence_toward",
        "subunit_hint", "subunit_evidence", "interpretation", "uniprot", "gene",
        "product", "ec", "organism", "pident", "coverage", "evalue", "next_step",
    ]
    write_tsv(
        tmp_path / "data/annotations_ab_initio/g1/protein_paths.tsv",
        header,
        [
            ["p1", "A", "A_named_function", "strong", "named_function", "none", "", "", "U1", "g1", "prod", "", "org", "100", "1", "0", ""],
            ["p2", "B", "B_uncharacterized_homolog", "moderate", "uncharacterized_family", "none", "", "", "U2", "g2", "DUF", "", "org", "50", "0.8", "1e-10", ""],
            ["p3", "C", "C_no_swissprot_hit", "none", "", "none", "", "", "", "", "", "", "org", "", "", "", ""],
        ],
    )
    write_tsv(
        tmp_path / "data/annotations_ab_initio/g2/protein_paths.tsv",
        header,
        [
            ["p4", "A", "A_named_function", "strong", "named_function", "named_subunit", "", "", "U4", "g4", "prod", "", "org", "100", "1", "0", ""],
        ],
    )
    write_text(
        tmp_path / "data/annotations_ab_initio/batch_summary.json",
        json.dumps({"totals": {"cds": 4}}),
    )
    return tmp_path


def test_summarize_abc_counts_and_not_evaluable(tmp_path: Path):
    root = _abc_fixture(tmp_path)
    section = summarize_abc(root)
    assert section["status"] == "partial"
    metrics = section["metrics"]
    assert metrics["genomes"] == 2
    assert metrics["proteins"] == 4
    assert metrics["labels"]["A_named_function"]["n"] == 2
    assert metrics["labels"]["B_uncharacterized_homolog"]["n"] == 1
    assert metrics["labels"]["C_no_swissprot_hit"]["n"] == 1
    assert metrics["batch_summary"]["matches_protein_rows"] is True
    assert metrics["with_uniprot"]["n"] == 3
    assert metrics["subunit_hint_named"]["n"] == 1
    assert any("función" in item["metric"].lower() for item in section["not_evaluable"])
    assert metrics["gold_domains_or_function"]["found"] is False


# ---------------------------------------------------------------------------
# Membrana sintética
# ---------------------------------------------------------------------------
def test_summarize_membrane_composition_and_structural_validation(tmp_path: Path):
    header = ["genome", "contig", "start", "end", "strand", "protein_id", "type",
              "family", "substrate", "evidence", "confidence"]
    tsv = write_tsv(
        tmp_path / "data/membrane/full/GCF_test.tsv",
        header,
        [
            ["GCF_test", "ctg", "1", "30", "+", "p1", "tm", "PF1", "", "kd", "high"],
            ["GCF_test", "ctg", "40", "90", "-", "p2", "transporter", "PF2", "sugar", "hmmscan", "medium"],
            ["GCF_test", "ctg", "100", "150", "+", "p3", "signal_peptide", "PF3", "", "sp", "low"],
        ],
    )
    assert tsv.exists()
    write_text(
        tmp_path / "data/membrane/validation_full.json",
        json.dumps({
            "ok": True,
            "processed_ok": True,
            "checked_at": "2026-01-01T00:00:00Z",
            "totals": {
                "genomes": 1, "rows": 3, "gff_features": 3,
                "errors": 0, "pending": 0, "genomes_ok": 1,
            },
        }),
    )
    section = summarize_membrane(tmp_path)
    metrics = section["metrics"]
    assert metrics["genomes"] == 1
    assert metrics["rows"] == 3
    assert metrics["types"]["tm"]["n"] == 1
    assert metrics["transporters"]["with_substrate"] == 1
    assert metrics["validation"]["consistent_with_tables"] is True
    assert metrics["validation"]["ok"] is True
    assert any("localización" in item["metric"] for item in section["not_evaluable"])


# ---------------------------------------------------------------------------
# Complejos sintéticos
# ---------------------------------------------------------------------------
def _complex_fixture(tmp_path: Path) -> Path:
    gold_header = ["complex_id", "complex_name", "subunits", "n_genes",
                   "locus_tags", "gene_names", "evidence_code", "evidence_class",
                   "properties", "source"]
    write_tsv(
        tmp_path / "data/complexes/gold_k12.tsv",
        gold_header,
        [
            ["G1", "Complejo 1", "3", "3", "a,b,c", "a,b,c", "", "experimental", "", "test"],
            ["G2", "Complejo 2", "2", "2", "d,e", "d,e", "", "inferida", "", "test"],
        ],
    )
    cluster_header = ["cluster_id", "genome", "n_genes", "genes", "gene_names",
                      "mean_score", "min_score", "confidence", "evidence", "evidencia",
                      "curado_id", "curado_nombre", "jaccard", "cobertura_cluster",
                      "curado_evidencia"]
    write_tsv(
        tmp_path / "data/complexes/clusters_ecoli_k12_evidencia.tsv",
        cluster_header,
        [
            ["C1", "ecoli_k12", "3", "a,b,c", "a,b,c", "1", "1", "alta", "", "documentado", "G1", "Complejo 1", "1", "1", "experimental"],
            ["C2", "ecoli_k12", "2", "d,e", "d,e", "1", "1", "media", "", "documentado", "G2", "Complejo 2", "1", "1", "inferida"],
            ["C3", "ecoli_k12", "2", "x,y", "x,y", "1", "1", "alta", "", "nuevo", "", "", "0", "0", ""],
        ],
    )
    return tmp_path


def test_summarize_complexes_internal_metrics(tmp_path: Path):
    root = _complex_fixture(tmp_path)
    section = summarize_complexes(root)
    metrics = section["metrics"]
    assert section["status"] == "internal_no_blind"
    assert metrics["benchmark_type"] == "internal_reconstruction_no_blind"
    assert metrics["blind"] is False
    assert metrics["gold"]["n"] == 2
    assert metrics["clusters"]["n"] == 3
    block = metrics["thresholds"]["0.8"]
    assert block["any_cluster_covers_gold"]["n"] == 2
    assert block["one_to_one"]["pairs"] == 2
    assert block["one_to_one"]["gold_recall"] == 1.0
    assert block["one_to_one"]["cluster_concordance"] == pytest.approx(0.6667, abs=1e-3)
    assert block["one_to_one"]["gene_recall"] == 1.0
    assert block["high_confidence_with_best_jaccard_ge"]["n"] == 1
    assert any("ciega" in item["metric"].lower() for item in section["not_evaluable"])


# ---------------------------------------------------------------------------
# oriC sintético
# ---------------------------------------------------------------------------
def _oric_fixture(tmp_path: Path) -> Path:
    header = ["genome", "status", "position", "confidence", "validated_by", "reason",
              "attempted", "lineage", "organism", "phylum", "dnaA", "orca_pos",
              "orca_conf", "doric", "skew", "boxes"]
    write_tsv(
        tmp_path / "data/orca_predictions_fix/oriC_status.tsv",
        header,
        [
            ["g1", "validated", "100", "0.9", "ORCA + DoriC + dnaA", "", "", "", "", "", "", "100", "0.9", "100", "", ""],
            ["g2", "consensus", "200", "0.8", "ORCA + proximidad dnaA", "", "", "", "", "", "", "200", "0.8", "", "", ""],
            ["g3", "anomalous", "900", "0.7", "ORCA", "", "", "", "", "", "", "900", "0.7", "100", "", ""],
        ],
    )
    write_text(
        tmp_path / "data/annotations_ab_initio/g1/annotation.gff3",
        "##gff-version 3\n##sequence-region chr 1 1000\n",
    )
    write_text(
        tmp_path / "data/annotations_ab_initio/g2/annotation.gff3",
        "##gff-version 3\n##sequence-region chr 1 1000\n",
    )
    write_text(
        tmp_path / "data/annotations_ab_initio/g3/annotation.gff3",
        "##gff-version 3\n##sequence-region chr 1 1000\n",
    )
    return tmp_path


def test_summarize_oric_separates_doric_from_consensus(tmp_path: Path):
    root = _oric_fixture(tmp_path)
    section = summarize_oric(root)
    metrics = section["metrics"]
    assert metrics["status_counts"]["validated"]["n"] == 1
    assert metrics["status_counts"]["consensus"]["n"] == 1
    assert metrics["status_counts"]["anomalous"]["n"] == 1
    assert metrics["doric_available"]["n"] == 2
    assert metrics["validated_anchored"]["n"] == 1
    assert metrics["validated_anchored"]["median_bp"] == 0
    assert metrics["anomalous_with_doric"]["median_bp"] == 200
    assert metrics["consensus"]["with_doric"] == 0
    assert metrics["consensus"]["distance"] == "not_evaluable"
    assert metrics["validated_anchored"]["n_using_doric"] == 1
    assert any("independiente" in item["metric"].lower() for item in section["not_evaluable"])


# ---------------------------------------------------------------------------
# Informe y build_data
# ---------------------------------------------------------------------------
def test_build_data_and_markdown_include_not_evaluable(tmp_path: Path):
    _abc_fixture(tmp_path)
    data = build_data(tmp_path, ("abc",))
    assert data["schema"] == "layer-benchmark/v1"
    assert "abc" in data["sections"]
    md = render_markdown(data)
    assert "# Benchmark reproducible por capas" in md
    assert "A/B/C" in md
    assert "No evaluable" in md
    assert "no evaluable" in md.lower()
    # El JSON sí conserva el contrato not_evaluable.
    assert data["sections"]["abc"]["not_evaluable"]
