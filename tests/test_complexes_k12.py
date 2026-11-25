"""Pruebas del buscador de complejos de `scripts/complexes`.

Nota de nombre: este fichero NO es `tests/test_complexes.py`, que cubre
`webapp/complexes.py` (complejos por firmas de nombres). Aquí se prueba el
módulo de `scripts/complexes/`, que puntúa pares y valida contra ComplexPortal.

Cubren lo que puede romperse en silencio: el parser de GenBank (una
continuación de `/product` que empieza con dígito no debe contaminar la
location), el mapeo de ComplexPortal, las features y el filtro top-k mutuo.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "complexes"))

import complexes_lib as L  # noqa: E402

GBFF = """\
LOCUS       TEST0001                  900 bp    DNA     linear   BCT 01-JAN-2020
ACCESSION   TEST0001
FEATURES             Location/Qualifiers
     CDS             1..300
                     /gene="thrA"
                     /locus_tag="b0001"
                     /product="fused protein with a wrapped line
                     1"
                     /db_xref="UniProtKB/Swiss-Prot:P00561"
     CDS             400..699
                     /gene="thrB"
                     /locus_tag="b0002"
                     /EC_number="2.7.1.39"
                     /db_xref="UniProtKB/Swiss-Prot:P00547"
     CDS             complement(750..899)
                     /locus_tag="b0003"
                     /db_xref="UniProtKB/Swiss-Prot:P00547"
ORIGIN
//
"""


def _genes(tmp_path: Path) -> dict[str, L.Gene]:
    p = tmp_path / "test.gbff"
    p.write_text(GBFF)
    return L.read_genes_gbff(p)


def test_gbff_location_no_se_contamina(tmp_path):
    genes = _genes(tmp_path)
    assert genes["b0001"].start == 1
    assert genes["b0001"].end == 300          # no 3001 por la línea "1\""
    assert genes["b0001"].uniprot == "P00561"
    assert genes["b0002"].ec == "2.7.1.39"
    assert genes["b0003"].strand == "-"
    assert genes["b0003"].start == 750


def test_features_genomicas(tmp_path):
    genes = _genes(tmp_path)
    ctx = L.build_context(genes, {"b0001": "op1", "b0002": "op1"})
    f = L.features_row("b0001", "b0002", ctx)
    assert f["same_operon"] == 1.0
    assert f["adjacent"] == 1.0
    assert f["log_gap"] > 0                   # 400-300-1 = 99 pb
    g = L.features_row("b0002", "b0003", ctx)
    assert g["same_operon"] == 0.0
    assert g["adjacent"] == 0.0               # hebras distintas
    assert L.pair_geometry("b0003", "b0001", ctx)["upstream"] == "b0001"


def test_pair_key_es_simetrica():
    assert L.pair_key("b0002", "b0001") == L.pair_key("b0001", "b0002") == ("b0001", "b0002")


def test_complexportal_ignora_ids_numericos(tmp_path):
    import csv
    rows = [
        ["Complex ac", "name", "aliases", "tax", "participants", "ev", "exp",
         "go", "x", "desc", "props", "assembly", "lig", "dis", "ago", "ant",
         "comment", "source", "expanded"],
        ["CPX-1", "complejo A", "-", "83333", "P00561(1)|P00547(1)", "-", "-",
         "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"],
        ["CPX-2", "homómero", "-", "83333", "P00561(2)", "-", "-", "-", "-",
         "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"],
        ["CPX-3", "con ligando", "-", "83333", "12345(1)|P00561(1)", "-", "-",
         "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"],
    ]
    p = tmp_path / "cp.tsv"
    with open(p, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)
    cx, stats = L.read_complexportal(p, {"P00561": "b0001", "P00547": "b0002"})
    assert stats["total"] == 3
    assert stats["homomers"] == 2             # CPX-2 y CPX-3 (ligando + 1 proteína)
    assert [c["complex_id"] for c in cx] == ["CPX-1"]
    assert cx[0]["evidence_class"] == "otra"


def test_clasifica_evidencia_complexportal():
    assert L.classify_evidence("ECO:0000353(physical interaction)") == "experimental"
    assert L.classify_evidence("ECO:0005547(background knowledge)") == "inferida"
    assert L.classify_evidence("ECO:0005544(orthology)") == "inferida"
    assert L.classify_evidence("") == "otra"


def test_cluster_top2_mutuo_no_encadena():
    # estrella: a es el hub; d solo conecta con a, que ya tiene 2 mejores
    edges = {L.pair_key("a", "b"): 0.95, L.pair_key("a", "c"): 0.90,
             L.pair_key("a", "d"): 0.85}
    cl = L.cluster_edges(edges, top_k=2, min_p=0.1)
    assert cl == [["a", "b", "c"]]


def test_cluster_dos_pares_disjuntos():
    edges = {L.pair_key("a", "b"): 0.8, L.pair_key("c", "d"): 0.7}
    cl = L.cluster_edges(edges, top_k=2, min_p=0.1)
    assert sorted(cl) == [["a", "b"], ["c", "d"]]


def test_modelo_aprende_y_predice():
    import pytest
    pytest.importorskip("sklearn")
    names = ["x"]
    X = [[0.0], [0.1], [0.2], [0.9], [1.0], [1.1]]
    y = [1, 1, 1, 0, 0, 0]
    m = L.fit_model(X, y, names)
    assert L.predict(m, [0.0]) > 0.5
    assert L.predict(m, [1.1]) < 0.5
    assert m["features"] == ["x"]


def test_candidate_pairs_prioriza_operon(tmp_path):
    genes = _genes(tmp_path)
    ctx = L.build_context(genes, {"b0001": "op1", "b0002": "op1", "b0003": "op1"})
    cand = set(L.candidate_pairs(ctx))
    assert L.pair_key("b0001", "b0002") in cand
    assert L.pair_key("b0001", "b0003") in cand
