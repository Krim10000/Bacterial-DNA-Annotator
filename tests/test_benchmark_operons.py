"""Tests de ``tools/benchmark_operons.py`` con predicciones/gold sintéticos.

No usan red ni los datos reales del proyecto: cada caso crea un mini proyecto
en ``tmp_path``. Comprueban las dos vistas (gen e intervalo), la separación por
hebra, el emparejamiento GFF y los ``not_evaluable`` honestos.
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

from tools.benchmark_operons import (  # noqa: E402
    Operon,
    evaluate,
    load_gff,
    main,
    normalize_strand,
    reciprocal_overlap,
    render_markdown,
    split_genes,
)

HEADER = ["operon_id", "contig", "strand", "start", "end", "genes", "n_genes"]


def write_tsv(path: Path, rows: list[list[object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)
    return path


def write_pred(path: Path, rows: list[list[object]]) -> Path:
    """Predicción con hebras 1/-1 para probar la normalización."""
    norm = []
    for oid, contig, strand, start, end, genes, n in rows:
        strand = {1: "1", -1: "-1"}.get(strand, strand)
        norm.append([oid, contig, strand, start, end, genes, n])
    return write_tsv(path, norm)


def make_case(tmp_path: Path) -> tuple[Path, Path]:
    # g7/g8 son genes válidos del universo gold que solo aparecen en
    # operones de un gen: unirlos en la predicción es un falso positivo que
    # debe verse en la vista por gen y en la de intervalo.
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "g1,g2", 2],
        ["G2", "chr1", "+", 1000, 1400, "g3,g4", 2],
        ["G3", "chr1", "-", 2000, 2400, "g5,g6", 2],
        ["G4", "chr1", "+", 3000, 3100, "g7", 1],
        ["G5", "chr1", "+", 3200, 3300, "g8", 1],
    ])
    pred = write_pred(tmp_path / "pred.tsv", [
        ["P1", "chr1", 1, 100, 300, "g1;g2", 2],
        ["P2", "chr1", 1, 1000, 1400, "g3;g4", 2],
        ["P3", "chr1", 1, 3000, 3400, "g7;g8", 2],   # falso positivo
        ["P4", "chr1", -1, 2000, 2400, "g5;g6", 2],
    ])
    return gold, pred


def test_metrics_by_gene_and_interval_split_by_strand(tmp_path: Path):
    gold, pred = make_case(tmp_path)
    res = evaluate(gold, pred, genome="synthetic", min_genes=2, root=tmp_path)

    assert res["status"] == "evaluated"
    assert res["evaluation_level"] == "full"
    assert res["counts"]["gt_operons_filtered"] == 3
    assert res["counts"]["pred_operons_filtered"] == 4
    assert res["counts"]["shared_genes"] == 8

    # Global: 3 pares verdaderos, 4 predichos y 3 gold.
    all_b = res["by_strand"]["all"]
    assert all_b["status"] == "evaluated"
    assert all_b["gene_pair"]["tp"] == 3
    assert all_b["gene_pair"]["precision"] == pytest.approx(0.75)
    assert all_b["gene_pair"]["recall"] == pytest.approx(1.0)
    assert all_b["gene_pair"]["f1"] == pytest.approx(0.8571, abs=1e-4)
    assert all_b["gene_exact"]["f1"] == pytest.approx(0.8571, abs=1e-4)
    assert all_b["interval_overlap"]["f1"] == pytest.approx(0.8571, abs=1e-4)
    assert all_b["interval_exact"]["n_exact"] == 3

    # Hebra +: el falso positivo baja precisión pero no recall.
    plus = res["by_strand"]["+"]
    assert plus["status"] == "evaluated"
    assert plus["gene_pair"]["precision"] == pytest.approx(2 / 3, abs=1e-4)
    assert plus["gene_pair"]["recall"] == pytest.approx(1.0)
    assert plus["gene_pair"]["f1"] == pytest.approx(0.8, abs=1e-4)
    assert plus["gene_exact"]["precision"] == pytest.approx(2 / 3, abs=1e-4)
    assert plus["interval_overlap"]["f1"] == pytest.approx(0.8, abs=1e-4)

    # Hebra -: perfecta y aislada del falso positivo de +.
    minus = res["by_strand"]["-"]
    assert minus["status"] == "evaluated"
    for key in ("gene_pair", "gene_exact", "interval_overlap",
                "interval_exact"):
        assert minus[key]["f1"] == pytest.approx(1.0)
    assert minus["interval_exact"]["n_exact"] == 1


def test_interval_overlap_accepts_shift_but_exact_interval_does_not(
        tmp_path: Path):
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "g1,g2", 2],
        ["G2", "chr1", "+", 1000, 1400, "g3,g4", 2],
    ])
    pred = write_tsv(tmp_path / "pred.tsv", [
        ["P1", "chr1", "+", 120, 310, "g1,g2", 2],
        ["P2", "chr1", "+", 1050, 1450, "g3,g4", 2],
    ])
    res = evaluate(gold, pred, genome="synthetic", min_genes=2,
                   min_overlap=0.5, root=tmp_path)
    plus = res["by_strand"]["+"]
    assert plus["gene_pair"]["f1"] == pytest.approx(1.0)
    assert plus["gene_exact"]["f1"] == pytest.approx(1.0)
    assert plus["interval_overlap"]["f1"] == pytest.approx(1.0)
    assert plus["interval_overlap"]["n_matched"] == 2
    assert plus["interval_exact"]["f1"] == pytest.approx(0.0)
    assert plus["interval_exact"]["n_exact"] == 0


def test_unmapped_gene_ids_keep_interval_metrics_evaluable(tmp_path: Path):
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "b0001,b0002", 2],
    ])
    pred = write_tsv(tmp_path / "pred.tsv", [
        ["P1", "chr1", "+", 100, 300, "cds_1;cds_2", 2],
    ])
    res = evaluate(gold, pred, genome="synthetic", min_genes=2, root=tmp_path)
    assert res["status"] == "evaluated"
    assert res["evaluation_level"] == "interval_only"
    assert res["counts"]["shared_genes"] == 0
    assert res["by_strand"]["all"]["gene_pair"]["status"] == "not_evaluable"
    assert (res["by_strand"]["all"]["gene_pair"]["reason"]
            == "sin_operones_predichos_con_genes_mapeados")
    assert res["by_strand"]["all"]["interval_overlap"]["status"] == "evaluated"
    assert res["by_strand"]["all"]["interval_overlap"]["f1"] == pytest.approx(1.0)
    assert any("Sin IDs de gen compartidos" in w for w in res["warnings"])


def test_no_gold_is_not_evaluable_with_reason(tmp_path: Path):
    json_path = tmp_path / "out.json"
    md_path = tmp_path / "out.md"
    code = main(["--root", str(tmp_path), "--genome", "synthetic",
                 "--json", str(json_path), "--md", str(md_path), "--quiet"])
    assert code == 0
    res = json.loads(json_path.read_text())
    assert res["status"] == "not_evaluable"
    assert res["reason_code"] == "gold_not_found"
    assert "No se encontró un gold local" in res["reason"]
    assert "No se usa red" in res["reason"]
    assert res["gold"]["searched"]
    assert md_path.exists()
    md = md_path.read_text()
    assert "no evaluable" in md.lower()
    assert "no se inventa" in md.lower()


def test_no_prediction_is_not_evaluable(tmp_path: Path):
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "g1,g2", 2],
    ])
    json_path = tmp_path / "out.json"
    code = main(["--root", str(tmp_path), "--gold", str(gold),
                 "--genome", "synthetic", "--json", str(json_path),
                 "--md", str(tmp_path / "out.md"), "--quiet"])
    assert code == 0
    res = json.loads(json_path.read_text())
    assert res["status"] == "not_evaluable"
    assert res["reason_code"] == "prediction_not_found"


def test_local_autodetection_prefers_genome_gold_and_operons_auto(
        tmp_path: Path):
    gold = tmp_path / "data" / "regulation" / "synthetic" / "gold_operons.tsv"
    write_tsv(gold, [["G1", "chr1", "+", 100, 300, "g1,g2", 2]])
    pred = tmp_path / "results" / "synthetic.operons_auto.tsv"
    write_tsv(pred, [["P1", "chr1", "+", 100, 300, "g1;g2", 2]])
    json_path = tmp_path / "out.json"
    code = main(["--root", str(tmp_path), "--genome", "synthetic",
                 "--json", str(json_path), "--md", str(tmp_path / "out.md"),
                 "--quiet"])
    assert code == 0
    res = json.loads(json_path.read_text())
    assert res["status"] == "evaluated"
    assert res["evaluation_level"] == "full"
    assert res["gold"]["path"] == "data/regulation/synthetic/gold_operons.tsv"
    assert res["pred"]["path"] == "results/synthetic.operons_auto.tsv"


def test_gff_prediction_with_genes_attribute_and_children(tmp_path: Path):
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "g1,g2", 2],
        ["G2", "chr1", "-", 1000, 1400, "g3,g4", 2],
    ])
    gff = tmp_path / "pred.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\toperon\t100\t300\t.\t+\t.\tID=op1;genes=g1,g2\n"
        "chr1\ttest\toperon\t1000\t1400\t.\t-\t.\tID=op2\n"
        "chr1\ttest\tCDS\t1000\t1100\t.\t-\t.\tID=cds3;Parent=op2;gene=g3\n"
        "chr1\ttest\tCDS\t1200\t1400\t.\t-\t.\tID=cds4;Parent=op2;gene=g4\n",
        encoding="utf-8")
    parsed = load_gff(gff)
    assert parsed.fmt == "gff"
    assert len(parsed.operons) == 2
    op2 = next(o for o in parsed.operons if o.operon_id == "op2")
    assert op2.genes == frozenset({"g3", "g4"})
    assert (op2.start, op2.end, op2.strand) == (1000, 1400, "-")

    res = evaluate(gold, gff, genome="synthetic", min_genes=2, root=tmp_path)
    assert res["status"] == "evaluated"
    assert res["pred"]["fmt"] == "gff"
    assert res["by_strand"]["+"]["gene_exact"]["f1"] == pytest.approx(1.0)
    assert res["by_strand"]["-"]["gene_exact"]["f1"] == pytest.approx(1.0)


def test_alias_file_maps_gene_names_to_bnumbers(tmp_path: Path):
    gold_dir = tmp_path / "data" / "regulation" / "synthetic"
    gold = write_tsv(gold_dir / "gold_operons.tsv", [
        ["G1", "chr1", "+", 100, 300, "b0001,b0002", 2],
    ])
    alias = gold_dir / "gt_gene_aliases.tsv"
    alias.write_text("name\tbnumber\nname1\tb0001\nname2\tb0002\n",
                     encoding="utf-8")
    pred = write_tsv(tmp_path / "pred.tsv", [
        ["P1", "chr1", "+", 100, 300, "name1;name2", 2],
    ])
    res = evaluate(gold, pred, genome="synthetic", min_genes=2, root=tmp_path)
    assert res["evaluation_level"] == "full"
    assert res["counts"]["aliased_pred_genes"] == 2
    assert res["by_strand"]["all"]["gene_pair"]["f1"] == pytest.approx(1.0)
    assert res["parameters"]["alias"].endswith("gt_gene_aliases.tsv")


def test_cli_explicit_paths_write_json_and_markdown(tmp_path: Path):
    gold, pred = make_case(tmp_path)
    json_path, md_path = tmp_path / "report.json", tmp_path / "report.md"
    code = main(["--gold", str(gold), "--pred", str(pred),
                 "--genome", "synthetic", "--json", str(json_path),
                 "--md", str(md_path), "--quiet"])
    assert code == 0
    res = json.loads(json_path.read_text())
    assert res["status"] == "evaluated"
    assert res["command"].startswith("python tools/benchmark_operons.py")
    md = md_path.read_text()
    assert "Benchmark de operones" in md
    assert "Métricas por hebra" in md
    assert "No es un holdout ciego" in md
    assert render_markdown(res) == md


def test_min_genes_filter_and_pure_helpers(tmp_path: Path):
    gold = write_tsv(tmp_path / "gold.tsv", [
        ["G1", "chr1", "+", 100, 300, "g1,g2", 2],
        ["G2", "chr1", "+", 500, 600, "g3", 1],
    ])
    pred = write_tsv(tmp_path / "pred.tsv", [
        ["P1", "chr1", "+", 100, 300, "g1,g2", 2],
        ["P2", "chr1", "+", 500, 600, "g3", 1],
    ])
    strict = evaluate(gold, pred, genome="s", min_genes=2, root=tmp_path)
    permissive = evaluate(gold, pred, genome="s", min_genes=1, root=tmp_path)
    assert strict["counts"]["gt_operons_filtered"] == 1
    assert permissive["counts"]["gt_operons_filtered"] == 2

    assert split_genes("a,b;c|d e") == ["a", "b", "c", "d", "e"]
    assert normalize_strand("1") == "+"
    assert normalize_strand("-1") == "-"
    assert normalize_strand("?") == "unknown"
    a = Operon("a", "chr1", "+", 100, 200, frozenset({"g1"}))
    b = Operon("b", "chr1", "+", 150, 250, frozenset({"g2"}))
    c = Operon("c", "chr1", "+", 500, 600, frozenset({"g3"}))
    assert reciprocal_overlap(a, b) == pytest.approx(51 / 101)
    assert reciprocal_overlap(a, c) == 0.0
