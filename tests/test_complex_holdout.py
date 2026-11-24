"""Pruebas sintéticas del holdout de complejos contra ComplexPortal.

El fixture reproduce el problema real: solo hay catálogo K-12 y el modelo se
entrenó en K-12. Por eso el holdout independiente debe ser ``not_evaluable`` y
los cálculos K-12 deben quedar marcados como ``in_sample``.
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

from tools import complex_holdout as ch  # noqa: E402


# --------------------------------------------------------------------------- #
# Fixture
# --------------------------------------------------------------------------- #
def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


GOLD_FIELDS = [
    "complex_id", "complex_name", "subunits", "n_genes", "locus_tags",
    "gene_names", "evidence_code", "evidence_class", "properties", "source",
]
RAW_FIELDS = [
    "#Complex ac", "Recommended name", "Taxonomy identifier",
    "Identifiers (and stoichiometry) of molecules in complex",
    "Evidence Code", "Expanded participant list",
]
CLUSTER_FIELDS = [
    "cluster_id", "genome", "n_genes", "genes", "gene_names", "mean_score",
    "min_score", "confidence", "evidence", "evidencia", "curado_id",
]
TRANSFER_FIELDS = [
    "complex_id", "complex_name", "genes_k12", "genes_genoma",
]


def build_fixture(tmp_path: Path, with_independent: bool = False) -> dict[str, Path]:
    complexes = tmp_path / "complexes"
    raw = complexes / "raw"
    complexes.mkdir()
    raw.mkdir(parents=True)

    # 8 genes K-12 con uniprot P1..P8.
    k12 = complexes / "k12_genes.tsv"
    _write_tsv(k12, ["contig", "start", "end", "strand", "locus_tag",
                     "gene_name", "ec", "uniprot"], [
        {"contig": "NC_000913", "start": 100 * i, "end": 100 * i + 90,
         "strand": "+", "locus_tag": f"b{i}", "gene_name": f"g{i}",
         "ec": "", "uniprot": f"P{i}"}
        for i in range(1, 9)
    ])

    gold = complexes / "gold_k12.tsv"
    _write_tsv(gold, GOLD_FIELDS, [
        {"complex_id": "CPX-G1", "complex_name": "G1", "subunits": "2",
         "n_genes": 2, "locus_tags": "b1,b2", "gene_names": "g1,g2",
         "evidence_code": "ECO:0005547", "evidence_class": "inferida",
         "properties": "Heterodimer", "source": "ComplexPortal"},
        {"complex_id": "CPX-G2", "complex_name": "G2", "subunits": "2",
         "n_genes": 2, "locus_tags": "b3,b4", "gene_names": "g3,g4",
         "evidence_code": "ECO:0000353", "evidence_class": "experimental",
         "properties": "Heterodimer", "source": "ComplexPortal"},
    ])

    # Raw: G1/G2 en gold; H1 raw-only con 2 genes; H2 raw-only con 1 gen.
    _write_tsv(raw / "complexportal_83333.tsv", RAW_FIELDS, [
        {"#Complex ac": "CPX-G1", "Recommended name": "G1",
         "Taxonomy identifier": "83333",
         "Identifiers (and stoichiometry) of molecules in complex":
             "P1(1)|P2(1)",
         "Evidence Code": "ECO:0005547(inference)",
         "Expanded participant list": "P1(1)|P2(1)"},
        {"#Complex ac": "CPX-G2", "Recommended name": "G2",
         "Taxonomy identifier": "83333",
         "Identifiers (and stoichiometry) of molecules in complex":
             "P3(1)|P4(1)",
         "Evidence Code": "ECO:0000353(physical)",
         "Expanded participant list": "P3(1)|P4(1)"},
        {"#Complex ac": "CPX-H1", "Recommended name": "Holdout H1",
         "Taxonomy identifier": "83333",
         "Identifiers (and stoichiometry) of molecules in complex":
             "P5(1)|P6(1)",
         "Evidence Code": "ECO:0000353(physical)",
         "Expanded participant list": "P5(1)|P6(1)"},
        {"#Complex ac": "CPX-H2", "Recommended name": "Holdout H2 single",
         "Taxonomy identifier": "83333",
         "Identifiers (and stoichiometry) of molecules in complex":
             "P7(1)",
         "Evidence Code": "ECO:0000353(physical)",
         "Expanded participant list": "P7(1)"},
    ])

    model = complexes / "model_k12.genomic.json"
    model.write_text(json.dumps({
        "trained_genome": "ecoli_k12",
        "gold": str(gold),
        "gold_complexes": 2,
        "n_pos": 4,
        "n_neg": 10,
        "features": ["same_operon"],
        "version": 1,
    }, ensure_ascii=False), encoding="utf-8")

    clusters = complexes / "clusters_ecoli_k12_evidencia.tsv"
    _write_tsv(clusters, CLUSTER_FIELDS, [
        {"cluster_id": "c1", "genome": "ecoli_k12", "n_genes": 2,
         "genes": "b1,b2", "gene_names": "g1,g2", "mean_score": 1.0,
         "min_score": 1.0, "confidence": "alta", "evidence": "operón",
         "evidencia": "documentado", "curado_id": "CPX-G1"},
        {"cluster_id": "c2", "genome": "ecoli_k12", "n_genes": 2,
         "genes": "b3,b4", "gene_names": "g3,g4", "mean_score": 1.0,
         "min_score": 1.0, "confidence": "alta", "evidence": "operón",
         "evidencia": "documentado", "curado_id": "CPX-G2"},
        {"cluster_id": "c3", "genome": "ecoli_k12", "n_genes": 2,
         "genes": "b7,b8", "gene_names": "g7,g8", "mean_score": 0.2,
         "min_score": 0.2, "confidence": "media", "evidence": "adyacente",
         "evidencia": "nuevo", "curado_id": ""},
        {"cluster_id": "c4", "genome": "ecoli_k12", "n_genes": 2,
         "genes": "b5,b6", "gene_names": "g5,g6", "mean_score": 0.3,
         "min_score": 0.3, "confidence": "baja", "evidence": "adyacente",
         "evidencia": "nuevo", "curado_id": ""},
    ])

    _write_tsv(complexes / "transfer_g1.tsv", TRANSFER_FIELDS, [
        {"complex_id": "CPX-G1", "complex_name": "G1",
         "genes_k12": "b1,b2", "genes_genoma": "g1,g2"},
        {"complex_id": "CPX-G2", "complex_name": "G2",
         "genes_k12": "b3,b4", "genes_genoma": "g3,g4"},
        {"complex_id": "CPX-G2", "complex_name": "G2",
         "genes_k12": "b3,b4", "genes_genoma": "g3,g4"},
    ])

    if with_independent:
        # Catálogo de un genoma distinto de K-12 + predicciones con el mismo
        # namespace de genes. Debe habilitar el holdout por genoma.
        gold_g2 = complexes / "gold_g2.tsv"
        _write_tsv(gold_g2, GOLD_FIELDS, [
            {"complex_id": "CPX-X1", "complex_name": "X1", "subunits": "2",
             "n_genes": 2, "locus_tags": "x1,x2", "gene_names": "x1,x2",
             "evidence_code": "ECO:0000353", "evidence_class": "experimental",
             "properties": "Heterodimer", "source": "ComplexPortal"},
        ])
        _write_tsv(complexes / "clusters_g2.tsv", CLUSTER_FIELDS, [
            {"cluster_id": "g2c1", "genome": "g2", "n_genes": 2,
             "genes": "x1,x2", "gene_names": "x1,x2", "mean_score": 1.0,
             "min_score": 1.0, "confidence": "alta", "evidence": "operón",
             "evidencia": "documentado", "curado_id": ""},
        ])

    return {
        "complexes": complexes,
        "raw": raw,
        "k12": k12,
        "model": model,
        "clusters": clusters,
        "out_json": tmp_path / "reports" / "holdout.json",
        "out_md": tmp_path / "reports" / "holdout.md",
    }


def run_analysis(tmp_path: Path, with_independent: bool = False) -> dict:
    paths = build_fixture(tmp_path, with_independent=with_independent)
    return ch.analyze(
        complexes_dir=paths["complexes"],
        raw_dir=paths["raw"],
        k12_genes_path=paths["k12"],
        model_json=paths["model"],
        clusters_path=paths["clusters"],
        out_json=paths["out_json"],
        out_md=paths["out_md"],
    )


# --------------------------------------------------------------------------- #
# Parseo y métricas
# --------------------------------------------------------------------------- #
def test_parse_participants_and_proteoforms():
    assert ch.parse_participants("P76052(1)|P77357(1)") == [
        ("P76052", 1), ("P77357", 1)]
    assert ch.parse_participants("P0AAJ8-PRO_0000042273(0)") == [
        ("P0AAJ8", 0)]
    assert ch.parse_participants("-") == []
    assert ch.parse_participants("") == []


def test_jaccard_and_best_match():
    a = frozenset({"b1", "b2"})
    b = frozenset({"b1", "b2"})
    c = frozenset({"b2", "b3"})
    assert ch.jaccard(a, b) == 1.0
    assert ch.jaccard(a, c) == 1 / 3
    assert ch.jaccard(a, frozenset()) == 0.0

    catalog = [
        ch.CatalogComplex("gold_k12", "CPX-1", "one", "ecoli_k12", a, 2, 2,
                          "", "ComplexPortal"),
        ch.CatalogComplex("gold_k12", "CPX-2", "two", "ecoli_k12", c, 2, 2,
                          "", "ComplexPortal"),
    ]
    jac, cid = ch.best_match(frozenset({"b1", "b2"}), catalog)
    assert jac == 1.0
    assert cid == "CPX-1"


def test_catalog_inventory_and_genome_aliases(tmp_path: Path):
    paths = build_fixture(tmp_path)
    uniprot = ch.load_uniprot_to_locus(paths["k12"])
    catalogs, inventory = ch.discover_catalogs(
        paths["complexes"], paths["raw"], uniprot)
    assert len(catalogs) == 6  # 2 gold + 4 raw
    by_catalog = {item["catalogo"]: item for item in inventory}
    assert by_catalog["gold_k12"]["complejos"] == 2
    assert by_catalog["gold_k12"]["en_gold"] == 2
    assert by_catalog["complexportal_83333"]["complejos"] == 4
    assert by_catalog["complexportal_83333"]["con_2_o_mas_genes_mapeados"] == 3
    assert by_catalog["complexportal_83333"]["con_menos_de_2_genes_mapeados"] == 1
    assert ch.genome_from_gold_name(Path("gold_k12.tsv")) == "ecoli_k12"


# --------------------------------------------------------------------------- #
# Holdout y contexto
# --------------------------------------------------------------------------- #
def test_independent_holdout_is_not_evaluable_without_other_catalog(tmp_path: Path):
    report = run_analysis(tmp_path)
    holdout = report["holdout"]
    assert holdout["estado"] == "not_evaluable"
    assert holdout["evaluable_independiente"] is False
    assert holdout["catalogos_independientes"] == []
    assert holdout["estrategias"]["genoma"]["estado"] == "not_evaluable"
    assert "K-12" in holdout["motivo"]

    assert report["procedencia_modelo"]["trained_genome"] == "ecoli_k12"
    assert report["procedencia_modelo"]["gold_complexes_encontrados"] == 2
    assert report["evaluaciones"]["in_sample_k12"]["estado"] == "in_sample_k12_no_holdout"


def test_in_sample_k12_metrics(tmp_path: Path):
    report = run_analysis(tmp_path)
    metrics = report["evaluaciones"]["in_sample_k12"]["metricas"]
    principal = metrics["principal"]
    assert metrics["predicciones"] == 4
    assert metrics["catalogo_elegible"] == 2
    assert principal["predicciones_tp"] == 2
    assert principal["catalogo_recuperado"] == 2
    assert principal["precision"] == 0.5
    assert principal["recall"] == 1.0
    assert principal["f1"] == pytest.approx(2 / 3, abs=1e-4)

    by_conf = report["evaluaciones"]["in_sample_k12"]["por_confianza"]
    assert by_conf["alta"]["precision"] == 1.0
    assert by_conf["media"]["precision"] == 0.0
    assert by_conf["baja"]["precision"] == 0.0

    by_ev = report["evaluaciones"]["in_sample_k12"]["por_evidencia"]
    assert by_ev["documentado"]["predicciones"] == 2
    assert by_ev["documentado"]["precision"] == 1.0


def test_complex_holdout_excludes_single_gene_and_reports_recall(tmp_path: Path):
    report = run_analysis(tmp_path)
    block = report["evaluaciones"]["holdout_complejo_k12"]
    assert block["estado"] == "exploratorio_no_independiente"
    assert block["catalogo_holdout"] == 2
    assert block["catalogo_holdout_evaluable"] == 1
    assert block["excluidos_menos_de_2_genes"] == 1
    assert "CPX-H2" in block["ids_excluidos"]
    assert block["recall"]["0.5"]["recall"] == 1.0
    assert block["detalle"][0]["complex_id"] == "CPX-H1"
    assert block["detalle"][0]["mejor_prediccion"] == "c4"
    assert block["precision"]["estado"] == "not_evaluable"


def test_transfer_coverage_is_not_precision(tmp_path: Path):
    report = run_analysis(tmp_path)
    block = report["evaluaciones"]["transferencias_k12_cobertura"]
    assert block["estado"] == "cobertura_derivada_no_holdout"
    assert block["filas_transferencia"] == 3
    assert block["complejos_transferidos"] == 2
    assert block["cobertura_complejos_gold"] == 1.0
    assert block["cobertura_genes_gold"] == 1.0
    assert block["complejos_gold_completos"] == 2
    assert "no precisión" in block["advertencia"]


def test_genome_holdout_evaluable_with_independent_catalog(tmp_path: Path):
    report = run_analysis(tmp_path, with_independent=True)
    holdout = report["holdout"]
    assert holdout["estado"] == "evaluable"
    assert holdout["evaluable_independiente"] is True
    assert holdout["catalogos_independientes"] == ["g2"]

    evaluation = holdout["estrategias"]["genoma"]["evaluaciones"][0]
    assert evaluation["estado"] == "evaluable"
    assert evaluation["metricas"]["principal"]["precision"] == 1.0
    assert evaluation["metricas"]["principal"]["recall"] == 1.0

    # El catálogo independiente no contamina el in-sample K-12.
    assert report["evaluaciones"]["in_sample_k12"]["metricas"]["catalogo_elegible"] == 2


# --------------------------------------------------------------------------- #
# Salidas
# --------------------------------------------------------------------------- #
def test_markdown_and_json_outputs(tmp_path: Path):
    paths = build_fixture(tmp_path)
    code = ch.main([
        "--complexes-dir", str(paths["complexes"]),
        "--raw-dir", str(paths["raw"]),
        "--k12-genes", str(paths["k12"]),
        "--model-json", str(paths["model"]),
        "--clusters", str(paths["clusters"]),
        "--out-json", str(paths["out_json"]),
        "--out-md", str(paths["out_md"]),
    ])
    assert code == 0
    assert paths["out_json"].is_file()
    assert paths["out_md"].is_file()

    report = json.loads(paths["out_json"].read_text(encoding="utf-8"))
    assert report["holdout"]["estado"] == "not_evaluable"
    assert report["evaluaciones"]["in_sample_k12"]["metricas"]["principal"]["precision"] == 0.5

    markdown = paths["out_md"].read_text(encoding="utf-8")
    assert "Veredicto: `not_evaluable`" in markdown
    assert "no es holdout" in markdown.lower()
    assert "Holdout por genoma" in markdown
    assert "Holdout por complejo" in markdown
    assert "no usa red" in markdown.lower()


def test_no_catalog_is_not_evaluable(tmp_path: Path):
    # Sin gold ni raw el holdout sigue siendo not_evaluable, no un cero.
    complexes = tmp_path / "empty" / "complexes"
    raw = complexes / "raw"
    complexes.mkdir(parents=True)
    raw.mkdir()
    report = ch.analyze(
        complexes_dir=complexes,
        raw_dir=raw,
        k12_genes_path=complexes / "missing.tsv",
        model_json=complexes / "missing.json",
        clusters_path=complexes / "missing_clusters.tsv",
        out_json=tmp_path / "j.json",
        out_md=tmp_path / "m.md",
    )
    assert report["holdout"]["estado"] == "not_evaluable"
    assert report["evaluaciones"]["in_sample_k12"]["estado"] == "not_evaluable"
