"""Pruebas de la capa TSS/TTS: normalización, descubrimiento y ausencia.

No usan red, GPU ni bases externas: todos los casos sintéticos se escriben en
``tmp_path``. La prueba final, opcional, comprueba la tabla real de *B.
subtilis* si el repositorio la tiene.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from webapp import tss_layer as T


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Normalización de valores
# --------------------------------------------------------------------------
def test_normalize_strand_pos_source_class():
    assert T.normalize_strand("+") == "+"
    assert T.normalize_strand(1) == "+"
    assert T.normalize_strand("forward") == "+"
    assert T.normalize_strand("-1") == "-"
    assert T.normalize_strand("reverse") == "-"
    assert T.normalize_strand("?") == "."

    assert T.normalize_pos("1,234") == 1234
    assert T.normalize_pos("250.0") == 250
    assert T.normalize_pos(0) is None
    assert T.normalize_pos("-5") is None
    assert T.normalize_pos("12.5") is None
    assert T.normalize_pos("abc") is None

    assert T.normalize_source("transferido_K12") == T.SOURCE_TRANSFER
    assert T.normalize_source("bsgatlas") == T.SOURCE_BSGATLAS
    assert T.normalize_source("alto rendimiento") == T.SOURCE_HT

    assert T.normalize_class("A", T.TSS) == "sigmaA"
    assert T.normalize_class("SigA", T.TSS) == "sigmaA"
    assert T.normalize_class("sigma70", T.TSS) == "sigma70"
    assert T.normalize_class("intergenic", T.TSS) == "unknown"
    assert T.normalize_class("RpoD", T.TSS) == "unknown"
    assert T.normalize_class("A_fuerte", T.TTS) == "strong"
    assert T.normalize_class("B_media", T.TTS) == "medium"
    assert T.normalize_class("C_debil", T.TTS) == "weak"

    assert T.normalize_location("intergénica") == "intergenic"
    assert T.normalize_location("antisentido") == "antisense"
    assert T.normalize_location("borde_de_gen") == "gene_edge"
    assert T.normalize_mechanism("rho-independent") == "rho_independent"


# --------------------------------------------------------------------------
# TSV combinado: normalización, duplicados y filas inválidas
# --------------------------------------------------------------------------
def test_load_combined_tsv_normalizes_and_dedupes(tmp_path: Path):
    _write(
        tmp_path / "dbtbs" / "mygenome_tss_tts.tsv",
        "type\tcontig\tpos\tstrand\tclass\tsource\tlocation\tend\n"
        "TSS\tc1\t100\tforward\tA\tdbtbs\tintergenic\t\n"
        "TSS\tc1\t100\t+\tSigA\tDBTBS\tintergenic\t\n"
        "TTS\tc1\t500\treverse\tA_fuerte\tbsgatlas\tgene_edge\t520\n"
        "TTS\tc1\t0\t+\tC_debil\tBSGatlas\t\t\n"
        "TSS\tc1\tabc\t+\tA\tDBTBS\t\t\n",
    )
    loaded = T.load("mygenome", tmp_path)

    assert loaded["status"] == T.STATUS_AVAILABLE
    assert loaded["available"] is True
    assert loaded["n_tss"] == 1
    assert loaded["n_tts"] == 1
    assert loaded["n_duplicates"] == 1
    assert loaded["n_rejected"] == 2
    assert loaded["reason"] == ""

    tss = loaded["tss"][0]
    assert tss["type"] == T.TSS
    assert tss["contig"] == "c1"
    assert tss["pos"] == 100
    assert tss["strand"] == "+"
    assert tss["class"] == "sigmaA"
    assert tss["source"] == T.SOURCE_DBTBS
    assert tss["location"] == "intergenic"

    tts = loaded["tts"][0]
    assert tts["type"] == T.TTS
    assert tts["pos"] == 500
    assert tts["end"] == 520
    assert tts["strand"] == "-"
    assert tts["class"] == "strong"
    assert tts["source"] == T.SOURCE_BSGATLAS
    assert tts["location"] == "gene_edge"
    # El archivo aportó 2 filas TSS; una se descartó como duplicada exacta.
    assert loaded["sources"][0]["n_tss"] == 2
    assert loaded["sources"][0]["n_tts"] == 1
    assert loaded["sources"][0]["n_rejected"] == 2


# --------------------------------------------------------------------------
# JSON con registros (no un resumen)
# --------------------------------------------------------------------------
def test_load_json_records(tmp_path: Path):
    _write(
        tmp_path / "regulation" / "g1" / "tss_tts.json",
        json.dumps({
            "tss": [
                {"contig": "c1", "pos": "250", "strand": "plus",
                 "sigma": "sigma70"},
            ],
            "tts": [
                {"contig": "c1", "pos": 900, "strand": "minus",
                 "class": "B_media", "source": "precalculado"},
            ],
        }),
    )
    loaded = T.load("g1", tmp_path)

    assert loaded["status"] == T.STATUS_AVAILABLE
    assert loaded["n_tss"] == 1
    assert loaded["n_tts"] == 1
    assert loaded["tss"][0]["pos"] == 250
    assert loaded["tss"][0]["strand"] == "+"
    assert loaded["tss"][0]["class"] == "sigma70"
    assert loaded["tss"][0]["source"] == T.SOURCE_PRECOMPUTED
    assert loaded["tts"][0]["strand"] == "-"
    assert loaded["tts"][0]["class"] == "medium"
    assert loaded["tts"][0]["source"] == T.SOURCE_PRECOMPUTED


def test_json_summary_without_positions_is_not_available(tmp_path: Path):
    _write(
        tmp_path / "dbtbs" / "mygenome_tss_tts.json",
        json.dumps({"regiones": 10, "total": 20,
                    "por_tipo": {"TSS": 10, "TTS": 10}}),
    )
    loaded = T.load("mygenome", tmp_path)

    assert loaded["status"] == T.STATUS_NOT_AVAILABLE
    assert loaded["available"] is False
    assert loaded["tss"] == []
    assert loaded["tts"] == []
    assert loaded["n_tss"] == 0 and loaded["n_tts"] == 0
    assert "no se predicen" in loaded["reason"]
    assert "resumen" in loaded["reason"]


# --------------------------------------------------------------------------
# Archivos separados en data/regulation/<genoma>/ y promotores ignorados
# --------------------------------------------------------------------------
def test_load_regulation_separate_files_ignores_promoters(tmp_path: Path):
    regulation = tmp_path / "regulation" / "g2"
    _write(
        regulation / "tss.tsv",
        "pos\tstrand\tsigma\tsource\n"
        "120\t+\tA\tdbtbs\n",
    )
    _write(
        regulation / "tts.tsv",
        "start\tend\tstrand\tclass\tsource\n"
        "400\t430\treverse\tC_debil\tbsgatlas\n",
    )
    # Un archivo de promotores putativos no es un TSS y no debe cargarse.
    _write(
        regulation / "promoters.tsv",
        "pos\tstrand\tsigma\tsource\n"
        "999\t+\tA\trule_based\n",
    )
    # Un resumen JSON no debe tapar los archivos separados con posiciones.
    _write(
        regulation / "tss_tts_meta.json",
        json.dumps({"total": 2, "por_tipo": {"TSS": 1, "TTS": 1}}),
    )
    loaded = T.load("g2", tmp_path)

    assert loaded["status"] == T.STATUS_AVAILABLE
    assert loaded["n_tss"] == 1
    assert loaded["n_tts"] == 1
    assert loaded["tss"][0]["pos"] == 120
    assert loaded["tss"][0]["class"] == "sigmaA"
    assert loaded["tss"][0]["source"] == T.SOURCE_DBTBS
    assert loaded["tts"][0]["pos"] == 400
    assert loaded["tts"][0]["end"] == 430
    assert loaded["tts"][0]["class"] == "weak"
    assert loaded["tts"][0]["source"] == T.SOURCE_BSGATLAS
    assert all("promoters" not in source["path"]
               for source in loaded["sources"])


# --------------------------------------------------------------------------
# Un genoma nuevo no hereda los datos de otra especie
# --------------------------------------------------------------------------
def test_new_genome_does_not_borrow_another_species(tmp_path: Path):
    _write(
        tmp_path / "dbtbs" / "subtilis_tss_tts.tsv",
        "type\tcontig\tpos\tstrand\tclass\tlocation\n"
        "TSS\tc1\t100\t+\tA\tintergenic\n"
        "TTS\tc1\t500\t+\tA_fuerte\tintergenic\n",
    )
    loaded = T.load("genoma_nuevo", tmp_path)

    assert loaded["status"] == T.STATUS_NOT_AVAILABLE
    assert loaded["n_tss"] == 0 and loaded["n_tts"] == 0
    assert "no se predicen" in loaded["reason"]

    # El alias sí encuentra la tabla de la misma especie.
    same = T.load("bsu", tmp_path)
    assert same["status"] == T.STATUS_AVAILABLE
    assert same["n_tss"] == 1 and same["n_tts"] == 1


def test_empty_data_root_is_not_available(tmp_path: Path):
    loaded = T.load("cualquiera", tmp_path)
    assert loaded["status"] == T.STATUS_NOT_AVAILABLE
    assert loaded["tss"] == [] and loaded["tts"] == []
    assert loaded["sources"] == []
    assert "no se predicen" in loaded["reason"]


def test_regulation_meta_only_does_not_block_dbtbs(tmp_path: Path):
    _write(
        tmp_path / "regulation" / "g4" / "tss_tts_meta.json",
        json.dumps({"total": 2}),
    )
    _write(
        tmp_path / "dbtbs" / "g4_tss_tts.tsv",
        "type\tcontig\tpos\tstrand\tclass\n"
        "TSS\tc1\t100\t+\tA\n"
        "TTS\tc1\t500\t+\tA_fuerte\n",
    )
    loaded = T.load("g4", tmp_path)
    assert loaded["status"] == T.STATUS_AVAILABLE
    assert loaded["n_tss"] == 1 and loaded["n_tts"] == 1
    assert loaded["sources"][0]["path"].endswith("g4_tss_tts.tsv")


def test_partial_when_only_one_type_is_present(tmp_path: Path):
    _write(
        tmp_path / "regulation" / "g5" / "tss.tsv",
        "pos\tstrand\tsigma\n100\t+\tA\n",
    )
    loaded = T.load("g5", tmp_path)
    assert loaded["status"] == T.STATUS_PARTIAL
    assert loaded["n_tss"] == 1 and loaded["n_tts"] == 0
    assert loaded["available"] is True
    assert "no hay TTS" in loaded["reason"]


# --------------------------------------------------------------------------
# attach: vecino más cercano por contig/hebra y ventana
# --------------------------------------------------------------------------
def test_attach_nearest_upstream_downstream_and_window():
    tss = [
        {"type": T.TSS, "contig": "c1", "pos": 100, "strand": "+",
         "source": "DBTBS", "class": "sigmaA"},
        {"type": T.TSS, "contig": "c1", "pos": 600, "strand": "+",
         "source": "DBTBS", "class": "sigmaA"},
        {"type": T.TSS, "contig": "c1", "pos": 1500, "strand": "+",
         "source": "DBTBS", "class": "sigmaA"},
        {"type": T.TSS, "contig": "c1", "pos": 900, "strand": "-",
         "source": "DBTBS", "class": "sigmaA"},
    ]
    tts = [
        {"type": T.TTS, "contig": "c1", "pos": 900, "strand": "+",
         "source": "BSGatlas", "class": "strong"},
        {"type": T.TTS, "contig": "c1", "pos": 100, "strand": "-",
         "source": "BSGatlas", "class": "medium"},
    ]
    features = [
        {"id": "plus", "contig": "c1", "start": 500, "end": 700,
         "strand": "+"},
        {"id": "minus", "contig": "c1", "start": 500, "end": 700,
         "strand": "-"},
        {"id": "lejos", "contig": "c1", "start": 5000, "end": 5100,
         "strand": "+"},
        {"id": "otro_contig", "contig": "c3", "start": 500, "end": 700,
         "strand": "+"},
        {"id": "sin_hebra", "contig": "c1", "start": 500, "end": 700,
         "strand": "?"},
    ]
    returned = T.attach(features, tss, tts)
    assert returned is features
    by_id = {feature["id"]: feature for feature in features}

    plus = by_id["plus"]
    assert plus["tss_upstream"]["pos"] == 100
    assert plus["tss_upstream_distance"] == 400
    assert plus["tts_downstream"]["pos"] == 900
    assert plus["tts_downstream_distance"] == 200
    assert plus["n_tss_inside"] == 1  # el TSS de 600 cae dentro
    assert plus["n_tts_inside"] == 0

    minus = by_id["minus"]
    assert minus["tss_upstream"]["pos"] == 900
    assert minus["tss_upstream_distance"] == 200
    assert minus["tts_downstream"]["pos"] == 100
    assert minus["tts_downstream_distance"] == 400

    lejos = by_id["lejos"]
    assert lejos["tss_upstream"] is None
    assert lejos["tts_downstream"] is None

    assert by_id["otro_contig"]["tss_upstream"] is None
    assert by_id["otro_contig"]["tts_downstream"] is None
    # Una hebra desconocida no se enlaza a ciegas.
    assert "tss_upstream" not in by_id["sin_hebra"]


def test_attach_respects_max_distance():
    features = [{"id": "x", "contig": "c1", "start": 1000, "end": 1100,
                 "strand": "+"}]
    tss = [{"type": T.TSS, "contig": "c1", "pos": 100, "strand": "+"}]
    tts = [{"type": T.TTS, "contig": "c1", "pos": 1500, "strand": "+"}]
    T.attach(features, tss, tts, max_distance=100)
    assert features[0]["tss_upstream"] is None
    assert features[0]["tts_downstream"] is None
    T.attach(features, tss, tts, max_distance=1000)
    assert features[0]["tss_upstream_distance"] == 900
    assert features[0]["tts_downstream_distance"] == 400


# --------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------
def test_summary_from_load_and_from_lists(tmp_path: Path):
    _write(
        tmp_path / "dbtbs" / "mygenome_tss_tts.tsv",
        "type\tcontig\tpos\tstrand\tclass\tsource\tlocation\n"
        "TSS\tc1\t100\t+\tA\tdbtbs\tintergenic\n"
        "TTS\tc1\t500\t-\tA_fuerte\tbsgatlas\tgene_edge\n",
    )
    loaded = T.load("mygenome", tmp_path)
    result = T.summary(loaded)

    assert result["status"] == T.STATUS_AVAILABLE
    assert result["n_tss"] == 1
    assert result["n_tts"] == 1
    assert result["n_total"] == 2
    assert result["tss_classes"] == {"sigmaA": 1}
    assert result["tts_classes"] == {"strong": 1}
    assert result["tss_locations"] == {"intergenic": 1}
    assert result["tts_locations"] == {"gene_edge": 1}
    assert result["by_strand"]["+"] == 1
    assert result["by_strand"]["-"] == 1
    assert result["pos_min"] == 100
    assert result["pos_max"] == 500

    direct = T.summary(loaded["tss"], loaded["tts"])
    assert direct["n_tss"] == 1
    assert direct["n_tts"] == 1
    assert direct["genome"] == ""

    empty = T.summary([], [])
    assert empty["status"] == T.STATUS_NOT_AVAILABLE
    assert empty["n_total"] == 0


# --------------------------------------------------------------------------
# Tabla real de B. subtilis, solo si el repositorio la tiene
# --------------------------------------------------------------------------
def test_real_subtilis_table_if_present():
    root = Path(__file__).resolve().parents[1] / "data"
    table = root / "dbtbs" / "subtilis_tss_tts.tsv"
    if not table.is_file():
        pytest.skip("no está data/dbtbs/subtilis_tss_tts.tsv")
    loaded = T.load("subtilis", root)
    assert loaded["status"] == T.STATUS_AVAILABLE
    assert loaded["n_tss"] == 3390
    assert loaded["n_tts"] == 3009
    assert loaded["tss"][0]["source"] == T.SOURCE_BSGATLAS
    assert loaded["tts"][0]["source"] == T.SOURCE_BSGATLAS
    assert loaded["n_duplicates"] == 0
    assert loaded["sources"][0]["path"].endswith("subtilis_tss_tts.tsv")
