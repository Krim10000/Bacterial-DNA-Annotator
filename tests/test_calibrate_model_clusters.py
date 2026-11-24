"""Pruebas sintéticas de agregación/prioridad de clusters del modelo K-12.

No usan los datos reales ni red. Construyen un mini-dataset con dos contigs,
seis clusters, cuatro transferencias, cinco firmas y tres complejos gold.
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

from tools import calibrate_model_clusters as cal  # noqa: E402


# --------------------------------------------------------------------------- #
# Fixture sintético
# --------------------------------------------------------------------------- #
def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


MODEL_FIELDS = [
    "complex_id", "name", "contig", "start", "end", "subunits", "expected",
    "missing", "status", "evidence", "confidence", "source", "ecs", "flag",
    "limitation",
]
TRANSFER_FIELDS = [
    "complex_id", "complex_name", "n_subunits", "n_mapped", "frac_mapped",
    "identity_mean", "score", "genes_k12", "genes_genoma", "observed_n",
    "expected_n", "stoich_status", "mass_kda", "phylo_jacc", "ecs",
    "evidence_class", "eco", "limitation",
]
SIGNATURE_FIELDS = [
    "cluster_id", "genome", "n_genes", "genes", "gene_names", "mean_score",
    "min_score", "confidence", "evidence", "evidencia", "curado_id",
    "curado_nombre", "jaccard", "cobertura_cluster", "curado_evidencia",
]
GOLD_FIELDS = [
    "complex_id", "complex_name", "subunits", "n_genes", "locus_tags",
    "gene_names", "evidence_code", "evidence_class", "properties", "source",
]
K12_FIELDS = ["contig", "start", "end", "strand", "locus_tag", "gene_name",
              "ec", "uniprot"]


def build_fixture(tmp_path: Path) -> dict[str, Path]:
    complexes = tmp_path / "complexes"
    reg = tmp_path / "regulation" / "g1"
    complexes.mkdir()
    reg.mkdir(parents=True)

    genes = [
        ("g1", "", "c1", 100, 200, "a"),
        ("g2", "", "c1", 300, 400, "b"),
        ("g3", "", "c1", 500, 600, "c"),
        ("g4", "", "c1", 700, 800, "d"),
        ("g5", "", "c1", 900, 1000, "e"),
        ("g6", "", "c1", 2000, 2100, "f"),
        ("g7", "", "c1", 3000, 3100, "g"),
        ("g8", "tag8", "c2", 100, 200, "h"),
        ("g9", "tag9", "c2", 300, 400, "i"),
    ]
    _write_tsv(
        reg / "genes.tsv",
        ["gene_id", "locus_tag", "contig", "start", "end", "gene_name"],
        [{"gene_id": gid, "locus_tag": tag, "contig": contig,
          "start": start, "end": end, "gene_name": name}
         for gid, tag, contig, start, end, name in genes])

    def model_row(cid, contig, start, end, subunits, confidence, evidence):
        return {
            "complex_id": cid, "name": f"cand {cid}", "contig": contig,
            "start": start, "end": end, "subunits": subunits,
            "expected": "", "missing": "", "status": "complex_candidate",
            "evidence": evidence, "confidence": confidence,
            "source": "modelo K-12 (build_candidates.py)", "ecs": "",
            "flag": "yellow", "limitation": "sospecha",
        }

    model_rows = [
        model_row("m1", "c1", 100, 400, "a:g1;b:g2", "alta",
                  "modelo K-12 (model_k12.genomic.json); operón;adyacente; "
                  "score medio=0.9; score mínimo=0.8"),
        model_row("m2", "c1", 100, 600, "a:g1;b:g2;c:g3", "alta",
                  "modelo K-12 (model_k12.genomic.json); operón;adyacente; "
                  "score medio=0.9"),
        model_row("m3", "c1", 2000, 3100, "a:g6;b:g7", "media",
                  "modelo K-12 (model_k12.genomic.json); operón"),
        model_row("m4", "c1", 700, 1000, "a:g4;b:g5", "baja",
                  "modelo K-12 (model_k12.genomic.json); adyacente"),
        model_row("m5", "c2", 100, 400, "a:tag8;b:tag9", "alta",
                  "modelo K-12 (model_k12.genomic.json); operón;adyacente"),
        model_row("m6", "c1", 5000, 5100, "a:u1;b:u2", "alta",
                  "modelo K-12 (model_k12.genomic.json); operón"),
    ]
    _write_tsv(complexes / "model_g1.complexes.tsv", MODEL_FIELDS, model_rows)

    def transfer_row(cid, k12, target, n=2, frac=1.0):
        return {
            "complex_id": cid, "complex_name": f"complex {cid}",
            "n_subunits": n, "n_mapped": n, "frac_mapped": frac,
            "identity_mean": 90.0, "score": 0.9, "genes_k12": k12,
            "genes_genoma": target, "observed_n": n, "expected_n": n,
            "stoich_status": "complete", "mass_kda": 10.0, "phylo_jacc": 0.5,
            "ecs": "", "evidence_class": "documentado_inferido",
            "eco": "ECO:0005544", "limitation": "RBH",
        }

    transfer_rows = [
        transfer_row("CPX-1", "b1,b2", "g1,g2"),
        transfer_row("CPX-2", "b4,b5", "g4,g5"),
        transfer_row("CPX-3", "b8,b9", "tag8,tag9"),
        transfer_row("CPX-4", "b7", "g7", n=1, frac=0.5),
    ]
    _write_tsv(complexes / "transfer_g1.tsv", TRANSFER_FIELDS, transfer_rows)

    def signature_row(cid, genes, confidence, curado=""):
        return {
            "cluster_id": cid, "genome": "ecoli_k12",
            "n_genes": len(genes.split(",")), "genes": genes,
            "gene_names": "", "mean_score": 1.0, "min_score": 1.0,
            "confidence": confidence, "evidence": "operón",
            "evidencia": "documentado", "curado_id": curado,
            "curado_nombre": "", "jaccard": 1.0, "cobertura_cluster": 1.0,
            "curado_evidencia": "inferida",
        }

    signature_rows = [
        signature_row("sig1", "b1,b2", "alta", "CPX-1"),
        signature_row("sig2", "b4,b5", "alta", "CPX-2"),
        signature_row("sig3", "b8,b9", "alta", "CPX-3"),
        signature_row("sig4", "b6,b7", "media"),
        signature_row("sig5", "b1,b2,b3", "baja"),
    ]
    signatures = tmp_path / "signatures.tsv"
    _write_tsv(signatures, SIGNATURE_FIELDS, signature_rows)

    gold = tmp_path / "gold.tsv"
    _write_tsv(gold, GOLD_FIELDS, [
        {"complex_id": cid, "complex_name": cid, "subunits": "2",
         "n_genes": 2, "locus_tags": tags, "gene_names": "",
         "evidence_code": "ECO:0005547", "evidence_class": "inferida",
         "properties": "-", "source": "ComplexPortal"}
        for cid, tags in [("CPX-1", "b1,b2"), ("CPX-2", "b4,b5"),
                          ("CPX-3", "b8,b9")]
    ])

    k12 = tmp_path / "k12_genes.tsv"
    k12_rows = []
    coords = {
        "b1": (100, 200), "b2": (300, 400), "b3": (500, 600),
        "b4": (700, 800), "b5": (900, 1000), "b6": (2000, 2100),
        "b7": (3000, 3100), "b8": (5000, 5100), "b9": (5200, 5300),
    }
    for tag, (start, end) in coords.items():
        k12_rows.append({
            "contig": "NC_000913", "start": start, "end": end,
            "strand": "+", "locus_tag": tag, "gene_name": tag,
            "ec": "", "uniprot": "",
        })
    _write_tsv(k12, K12_FIELDS, k12_rows)

    return {
        "complexes": complexes,
        "regulation": tmp_path / "regulation",
        "signatures": signatures,
        "gold": gold,
        "k12": k12,
    }


def run_analysis(tmp_path: Path, **kwargs) -> dict:
    paths = build_fixture(tmp_path)
    return cal.analyze(
        complexes_dir=paths["complexes"],
        regulation_dir=paths["regulation"],
        signatures_path=paths["signatures"],
        gold_path=paths["gold"],
        k12_genes_path=paths["k12"],
        **kwargs,
    )


def _by_id(report: dict) -> dict[str, dict]:
    return {row["complex_id"]: row for row in report["priorizacion"]["top"]}


# --------------------------------------------------------------------------- #
# Utilidades puras
# --------------------------------------------------------------------------- #
def test_parse_gene_ids_and_aliases():
    assert cal.parse_gene_ids("atpA:g1;atpB:g2") == ["g1", "g2"]
    assert cal.parse_gene_ids("rna-TC_RS1,TC_RS2") == ["rna-TC_RS1", "TC_RS2"]
    assert cal.parse_gene_ids("g1;g1;g2") == ["g1", "g2"]
    assert cal.parse_gene_ids("") == []


def test_size_bin_label():
    assert cal.size_bin_label(2) == "2"
    assert cal.size_bin_label(4) == "4"
    assert cal.size_bin_label(5) == "5-6"
    assert cal.size_bin_label(6) == "5-6"
    assert cal.size_bin_label(7) == "7-9"
    assert cal.size_bin_label(10) == "10+"
    assert cal.size_bin_label(95) == "10+"


def test_interval_overlap_fraction():
    assert cal.interval_overlap_fraction(100, 200, 150, 250) == pytest.approx(0.5)
    assert cal.interval_overlap_fraction(100, 200, 100, 200) == 1.0
    assert cal.interval_overlap_fraction(100, 200, 201, 300) == 0.0
    assert cal.interval_overlap_fraction(0, 100, 10, 20) == 0.0
    assert cal.interval_overlap_fraction(100, 200, 120, 130) == pytest.approx(0.1)


def test_evidence_tokens_collapse_scores():
    tokens = cal.evidence_tokens(
        "modelo K-12 (model_k12.genomic.json); operón;adyacente; "
        "score medio=0.9; score mínimo=0.8")
    assert tokens == ["modelo K-12", "operón", "adyacente",
                      "score medio", "score mínimo"]


def test_priority_score_prefers_confidence_and_overlap():
    high = cal.priority_score("alta", 2, 1.0, 1.0, 1.0, 1.0)
    low = cal.priority_score("baja", 2, 0.0, 0.0, 0.0, 0.0)
    assert high > low
    assert cal.priority_tier(high) == "P1"
    assert cal.priority_tier(low) == "P3"


# --------------------------------------------------------------------------- #
# Agregación
# --------------------------------------------------------------------------- #
def test_aggregation_by_confidence_size_and_genome(tmp_path: Path):
    report = run_analysis(tmp_path)
    totals = report["totales"]
    assert totals["genomas"] == 1
    assert totals["clusters"] == 6
    assert totals["confianza"] == {"alta": 4, "media": 1, "baja": 1}
    assert totals["subunidades_listadas"] == 13
    assert totals["clusters_sin_genes_resueltos"] == 1

    assert report["por_confianza"]["alta"]["clusters"] == 4
    assert report["por_confianza"]["media"]["clusters"] == 1
    assert report["por_confianza"]["baja"]["clusters"] == 1
    assert report["por_confianza"]["alta"]["tamano"]["mean"] == pytest.approx(2.25)

    assert report["por_tamano"]["2"]["clusters"] == 5
    assert report["por_tamano"]["3"]["clusters"] == 1

    genome = report["por_genoma"]["g1"]
    assert genome["clusters"] == 6
    assert genome["confianza"]["alta"] == 4
    assert genome["solapamiento"]["has_transfer_gene"] == 5


def test_evidence_aggregation(tmp_path: Path):
    report = run_analysis(tmp_path)
    top = {row["token"]: row["clusters"]
           for row in report["evidencia"]["top_global"]}
    assert top["modelo K-12"] == 6
    assert top["operón"] == 5
    assert top["adyacente"] == 4
    assert top["score medio"] == 2
    assert report["agregado_global"]["evidencia_top"][0]["clusters"] >= 5


# --------------------------------------------------------------------------- #
# Solapamiento
# --------------------------------------------------------------------------- #
def test_transfer_overlap_gene_and_interval(tmp_path: Path):
    report = run_analysis(tmp_path)
    by_id = _by_id(report)
    m1 = by_id["m1"]
    assert m1["t_best_id"] == "CPX-1"
    assert m1["metrics"]["t_gene_cov"] == 1.0
    assert m1["metrics"]["t_gene_jac"] == 1.0
    assert m1["metrics"]["t_iv_frac"] == 1.0
    assert m1["flags"]["strong_transfer_gene"] is True

    m2 = by_id["m2"]
    assert m2["t_best_id"] == "CPX-1"
    assert m2["metrics"]["t_gene_cov"] == pytest.approx(2 / 3, abs=1e-4)
    assert m2["metrics"]["t_gene_jac"] == pytest.approx(2 / 3, abs=1e-4)
    assert m2["metrics"]["t_iv_frac"] == pytest.approx(0.6, abs=1e-4)

    m3 = by_id["m3"]
    assert m3["t_best_id"] == "CPX-4"
    assert m3["metrics"]["t_gene_cov"] == 0.5
    assert m3["metrics"]["t_iv_frac"] == pytest.approx(100 / 1100, abs=1e-4)
    assert m3["flags"]["strong_transfer_iv"] is False

    m6 = by_id["m6"]
    assert m6["metrics"]["t_gene_cov"] == 0.0
    assert m6["flags"]["has_transfer_gene"] is False

    agg = report["agregado_global"]["solapamiento"]
    assert agg["has_transfer_gene"] == 5
    assert agg["strong_transfer_gene"] == 5
    assert agg["has_transfer_iv"] == 5
    assert agg["strong_transfer_iv"] == 4
    assert agg["has_bridge"] == 5
    assert agg["bridge_signature_backed"] == 4


def test_signature_overlap_through_bridge(tmp_path: Path):
    report = run_analysis(tmp_path)
    by_id = _by_id(report)

    m1 = by_id["m1"]
    assert m1["s_best_id"] == "sig1"
    assert m1["metrics"]["s_gene_jac"] == 1.0
    assert m1["metrics"]["s_iv_frac"] == 1.0

    m2 = by_id["m2"]
    assert m2["metrics"]["s_gene_jac"] == 1.0
    assert m2["s_best_id"] in {"sig1", "sig5"}

    m3 = by_id["m3"]
    assert m3["s_best_id"] == "sig4"
    assert m3["metrics"]["s_gene_jac"] == 0.5
    assert m3["metrics"]["s_iv_frac"] == 1.0

    m6 = by_id["m6"]
    assert m6["flags"]["has_bridge"] is False
    assert m6["flags"]["has_signature_gene"] is False
    assert m6["metrics"]["s_gene_jac"] == 0.0

    firmas = report["solapamiento"]["firmas"]
    assert firmas["firmas"] == 5
    assert firmas["clusters_con_puente_k12"] == 5
    assert firmas["clusters_sin_puente_k12"] == 1
    assert firmas["clusters_con_match_firma_por_gen"] == 5
    assert firmas["firmas_con_match_por_gen"] == 4


def test_colocalized_interval_filter(tmp_path: Path):
    # Con umbral 250 pb, ningún intervalo del fixture (300 pb o más) es
    # colocalizado: el intervalo amplio debe seguir existiendo como diagnóstico.
    report = run_analysis(tmp_path, max_colocal_span=250)
    by_id = _by_id(report)
    m1 = by_id["m1"]
    assert m1["metrics"]["t_iv_frac"] == 0.0
    assert m1["metrics"]["t_iv_frac_any"] == 1.0
    assert m1["flags"]["has_transfer_iv"] is False
    assert m1["flags"]["has_transfer_iv_any"] is True

    transfer = report["solapamiento"]["transferencias"]
    # CPX-4 tiene un solo gen (span 100 pb), por eso sí es colocalizado.
    assert transfer["filas_con_intervalo_colocalizado"] == 1
    assert transfer["filas_sin_intervalo_colocalizado"] == 3
    firmas = report["solapamiento"]["firmas"]
    assert firmas["firmas_con_intervalo_colocalizado"] == 0
    assert firmas["clusters_con_puente_k12_colocalizado"] == 1
    assert report["inputs"]["max_span_colocalizado"] == 250


def test_gold_coverage(tmp_path: Path):
    report = run_analysis(tmp_path)
    gold = report["solapamiento"]["gold"]
    assert gold["complejos_gold"] == 3
    assert gold["complejos_con_transferencia"] == 3
    assert gold["complejos_transferidos"] == 4
    assert gold["complejos_transferidos_en_gold"] == 3
    assert gold["firmas_con_curado_id"] == 3
    assert gold["firmas_sin_curado_id"] == 2
    assert gold["firmas_curado_id_en_gold"] == 3


# --------------------------------------------------------------------------- #
# Prioridad
# --------------------------------------------------------------------------- #
def test_priority_ranks_high_confidence_with_transfer_first(tmp_path: Path):
    report = run_analysis(tmp_path)
    top = report["priorizacion"]["top"]
    assert [row["complex_id"] for row in top[:3]] == ["m1", "m5", "m2"]
    assert top[0]["priority_tier"] == "P1"
    assert top[0]["metrics"]["priority_score"] > top[2]["metrics"]["priority_score"]

    by_id = _by_id(report)
    # m4 (baja, cobertura total) supera a m3 (media, cobertura parcial).
    assert by_id["m4"]["rank"] < by_id["m3"]["rank"]
    # m6 (alta sin anclaje) queda debajo de m4 y m3.
    assert by_id["m6"]["rank"] == 6
    assert report["priorizacion"]["conteo_tiers"] == {"P1": 3, "P2": 1, "P3": 2}


def test_priority_top_n_is_respected(tmp_path: Path):
    report = run_analysis(tmp_path, top_n=3)
    assert report["priorizacion"]["top_n"] == 3
    assert len(report["priorizacion"]["top"]) == 3
    assert report["priorizacion"]["top"][0]["rank"] == 1


# --------------------------------------------------------------------------- #
# Salidas
# --------------------------------------------------------------------------- #
def test_markdown_and_json_outputs(tmp_path: Path):
    paths = build_fixture(tmp_path)
    out_json = tmp_path / "reports" / "calibration.json"
    out_md = tmp_path / "reports" / "calibration.md"
    code = cal.main([
        "--complexes-dir", str(paths["complexes"]),
        "--regulation-dir", str(paths["regulation"]),
        "--signatures", str(paths["signatures"]),
        "--gold", str(paths["gold"]),
        "--k12-genes", str(paths["k12"]),
        "--out-json", str(out_json),
        "--out-md", str(out_md),
    ])
    assert code == 0
    assert out_json.is_file()
    assert out_md.is_file()

    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["totales"]["clusters"] == 6
    assert report["priorizacion"]["top"][0]["complex_id"] == "m1"

    markdown = out_md.read_text(encoding="utf-8")
    assert "# Calibración de clusters del modelo K-12" in markdown
    assert "## 5. Limitaciones" in markdown
    assert "m1" in markdown
    assert "priorización" in markdown.lower()
    assert "no es una probabilidad" in markdown.lower()


def test_limit_genomes(tmp_path: Path):
    report = run_analysis(tmp_path, limit_genomes=1)
    assert report["inputs"]["genomas_modelo"] == 1
    assert report["totales"]["clusters"] == 6
