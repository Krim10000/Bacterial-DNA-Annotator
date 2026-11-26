"""Pruebas puras de regiones oriC y extracción de secuencia.

No usan Streamlit, pandas, red ni los datos de 400 genomas: solo secuencias y
regiones sintéticas. El contrato del índice se comprueba aparte y se salta si
no se ha generado ``data/oric_regions/oriC_regions.json``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from webapp import oric_regions as ori


def test_normalize_region_keeps_published_endpoints():
    """No se reescribe el endpoint publicado, aunque la longitud difiera 1 pb."""
    span = ori.normalize_region(100, 199, length=100, genome_length=1000)
    assert span is not None
    assert span.segments == ((100, 199),)
    assert span.length == 100
    assert span.length_delta == 0
    assert not span.wrap

    # Convención exclusiva: end=200 con length=100. Se conserva 200 y se
    # declara la diferencia; no se recorta a 199.
    span = ori.normalize_region(100, 200, length=100, genome_length=1000)
    assert span is not None
    assert span.segments == ((100, 200),)
    assert span.length == 101
    assert span.length_delta == 1


def test_normalize_region_wraps_origin():
    span = ori.normalize_region(95, 10, length=16, genome_length=100)
    assert span is not None
    assert span.wrap
    assert span.segments == ((95, 100), (1, 10))
    assert span.length == 16
    assert span.length_delta == 0


def test_normalize_region_end_zero_uses_reported_length():
    """``end == 0`` es un centinela; la longitud publicada manda."""
    span = ori.normalize_region(4144817, 0, length=634,
                                genome_length=4145451)
    assert span is not None
    assert not span.wrap
    assert span.segments == ((4144817, 4145450),)
    assert span.length == 634
    assert span.length_delta == 0

    # Si la longitud cruza el origen, se parte en dos segmentos.
    span = ori.normalize_region(95, 0, length=10, genome_length=100)
    assert span is not None
    assert span.wrap
    assert span.segments == ((95, 100), (1, 4))
    assert span.length == 10


def test_normalize_region_rejects_invalid_or_ambiguous():
    assert ori.normalize_region(0, 100, genome_length=1000) is None
    assert ori.normalize_region(100, 1001, genome_length=1000) is None
    assert ori.normalize_region(100, 200, length=0, genome_length=1000) is None
    # Sin longitud del genoma no se puede resolver un cruce de origen.
    assert ori.normalize_region(95, 10, length=16, genome_length=None) is None


def test_extract_region_sequence_non_wrap():
    seq = "AAAACCCCGGTT"
    assert ori.extract_region_sequence({"c": seq}, "c", 2, 5) == "AAAC"
    assert ori.extract_region_sequence({"c": seq}, "c", 1, 1) == "A"
    assert ori.extract_region_sequence({"c": seq}, "c", 12, 12) == "T"


def test_extract_region_sequence_wraps_origin():
    seq = "AAAACCCCGG"
    # 8..10 + 1..2 = CGG + AA
    assert ori.extract_region_sequence({"c": seq}, "c", 8, 2) == "CGGAA"


def test_extract_region_sequence_end_zero_by_length():
    seq = "AAAACCCCGG"
    # end=0 con length=2: el intervalo es 8..9, no 8..10.
    assert ori.extract_region_sequence(
        {"c": seq}, "c", 8, 0, length=2) == "CG"


def test_extract_region_sequence_invalid_returns_empty():
    seq = "AAAACCCCGG"
    assert ori.extract_region_sequence({"c": seq}, "c", 0, 5) == ""
    assert ori.extract_region_sequence({"c": seq}, "c", 5, 20) == ""
    assert ori.extract_region_sequence({"a": seq, "b": seq}, "z", 1, 3) == ""
    assert ori.extract_region_sequence({}, "c", 1, 3) == ""


def test_extract_region_sequence_single_contig_fallback():
    seq = "AAAACCCCGG"
    assert ori.extract_region_sequence({"unico": seq}, "", 1, 4) == "AAAA"


def test_region_contains_and_distance_wrap():
    span = ori.normalize_region(95, 10, length=16, genome_length=100)
    assert span is not None
    assert span.contains(95) and span.contains(100)
    assert span.contains(1) and span.contains(10)
    assert not span.contains(50)
    assert span.distance_to(50) == 40
    assert span.distance_to(94) == 1
    assert span.distance_to(11) == 1


def test_select_doric_region_prefers_agreement_with_point():
    records = [
        {"start": 100, "end": 200, "length": 101,
         "comparative_support": "Singleton"},
        {"start": 900, "end": 950, "length": 51,
         "comparative_support": "Recurrent"},
    ]
    chosen = ori.select_doric_region(records, genome_length=1000, point=950)
    assert chosen is not None
    assert chosen["selection"] == "concordante_con_el_punto"
    assert chosen["selected"]["start"] == 900
    assert chosen["selected"]["point_inside"] is True
    assert len(chosen["alternatives"]) == 1


def test_select_doric_region_uses_closest_when_point_is_outside():
    records = [
        {"start": 100, "end": 200, "length": 101},
        {"start": 800, "end": 900, "length": 101},
    ]
    chosen = ori.select_doric_region(records, genome_length=1000, point=850)
    assert chosen is not None
    assert chosen["selection"] == "concordante_con_el_punto"
    assert chosen["selected"]["start"] == 800

    chosen = ori.select_doric_region(records, genome_length=1000, point=750)
    assert chosen is not None
    assert chosen["selection"] == "mas_cercana_al_punto"
    assert chosen["selected"]["start"] == 800
    assert chosen["selected"]["distance_bp"] == 50


def test_select_doric_region_without_point_prefers_recurrence():
    records = [
        {"start": 100, "end": 200, "length": 101,
         "comparative_support": "Singleton"},
        {"start": 800, "end": 900, "length": 101,
         "comparative_support": "Recurrent"},
    ]
    chosen = ori.select_doric_region(records, genome_length=1000, point=None)
    assert chosen is not None
    assert chosen["selection"] == "mayor_recurrencia"
    assert chosen["selected"]["start"] == 800


def test_resolve_index_entry_exposes_compatible_keys():
    entry = {
        "genome": "G", "pos": 12, "genome_length": 100,
        "confidence": 0.9, "source": "DoriC-anchored",
        "span_source": "DoriC", "replicon": "NC_1",
        "region": {"start": 95, "end": 10, "length_reported": 16},
        "point": {"kind": "DoriC (inicio publicado)"},
    }
    out = ori.resolve_index_entry(entry)
    assert out is not None
    assert out["pos"] == 12
    assert out["start"] == 95 and out["end"] == 10
    assert out["span_source"] == "DoriC"
    assert out["region"]["wrap"] is True
    assert out["region"]["segments"] == [[95, 100], [1, 10]]
    assert out["length"] == 100  # compatibilidad: longitud del genoma


def test_resolve_index_entry_point_only_has_no_interval():
    entry = {
        "genome": "G", "pos": 42, "genome_length": 1000,
        "confidence": 0.8, "source": "ORCA",
        "span_source": "ORCA (punto)", "region": None,
    }
    out = ori.resolve_index_entry(entry)
    assert out is not None
    assert out["start"] == out["end"] == out["pos"] == 42
    assert out["region"] is None
    assert out["span_source"] == "ORCA (punto)"


INDEX = Path(__file__).resolve().parents[1] / "data" / "oric_regions" / "oriC_regions.json"


@pytest.mark.skipif(not INDEX.exists(), reason="índice de regiones no generado")
def test_index_contract():
    data = json.loads(INDEX.read_text(encoding="utf-8"))
    genomes = data["genomes"]
    assert len(genomes) == 400
    assert data["counts"]["genomes"] == 400
    assert data["counts"]["with_region"] == 321
    assert data["counts"]["point_only"] == 79
    seen = 0
    for genome, entry in genomes.items():
        assert entry["span_source"] in ("DoriC", "ORCA (punto)")
        assert 1 <= entry["pos"] <= entry["genome_length"]
        if entry["span_source"] == "ORCA (punto)":
            assert entry["start"] == entry["end"] == entry["pos"]
            assert entry["region"] is None
            assert entry["span_length"] == 1
            seen += 1
            continue
        region = entry["region"]
        span = ori.normalize_region(
            region["start"], region["end"],
            length=region["length_reported"],
            genome_length=entry["genome_length"])
        assert span is not None, genome
        assert span.segments == tuple(tuple(s) for s in region["segments"])
        assert span.length == region["length"]
        assert region["length_delta"] in (None, -1, 0, 1)
        if "posible plásmido" in " ".join(entry.get("notes") or []):
            # La GT de plásmido se descarta: o queda el punto, o una fila
            # cromosómica local elegida por recurrencia (nunca doric_ground_truth).
            assert region["selection"] != "doric_ground_truth"
    assert seen == data["counts"]["point_only"]
    assert data["counts"]["gt_discarded_not_chromosome"] == 5


@pytest.mark.skipif(not INDEX.exists(), reason="índice de regiones no generado")
def test_load_region_from_index_wrap_and_point():
    data_root = INDEX.parents[2]  # .../data
    entry = ori.load_region(data_root, "GCF_000006685.1_ASM668v1")
    assert entry is not None
    assert entry["span_source"] == "DoriC"
    assert entry["region"]["wrap"] is True
    assert entry["region"]["segments"] == [[1072919, 1072950], [1, 265]]

    point = ori.load_region(data_root, "GCF_000508245.1_ASM50824v1")
    assert point is not None
    assert point["span_source"] == "ORCA (punto)"
    assert point["start"] == point["end"] == point["pos"] == 600000


def test_fallback_reads_prediction_schema(tmp_path):
    """El esquema de ``orca_predictions`` usa ``orca_position``, no ``chosen``."""
    pred_dir = tmp_path / "data" / "orca_predictions"
    pred_dir.mkdir(parents=True)
    (pred_dir / "G.json").write_text(json.dumps({
        "genome": "G",
        "replicon": "NC_1",
        "orca_position": 150,
        "orca_confidence": 0.75,
        "doric_truth": {"orics": [{"start": 100, "end": 200,
                                   "length": 101}]},
    }), encoding="utf-8")
    entry = ori.fallback_from_prediction(tmp_path / "data", "G")
    assert entry is not None
    assert entry["pos"] == 150
    assert entry["start"] == 100 and entry["end"] == 200
    assert entry["span_source"] == "DoriC"
    assert entry["confidence"] == 0.75


def test_fallback_point_only_without_region(tmp_path):
    pred_dir = tmp_path / "data" / "orca_predictions"
    pred_dir.mkdir(parents=True)
    (pred_dir / "G.json").write_text(json.dumps({
        "genome": "G", "replicon": "NC_1",
        "orca_position": 42, "orca_confidence": 0.5,
        "doric_truth": None,
    }), encoding="utf-8")
    entry = ori.fallback_from_prediction(tmp_path / "data", "G")
    assert entry is not None
    assert entry["start"] == entry["end"] == entry["pos"] == 42
    assert entry["span_source"] == "ORCA (punto)"
    assert entry["region"] is None
