"""Estructura secundaria de ARN: predicción determinista, no 3D.

Casos conocidos con ViennaRNA (si está instalado) y casos sin herramienta:
el respaldo Nussinov no inventa pares y ``allow_fallback=False`` devuelve
``not_available``.
"""
from __future__ import annotations

import json

import pytest

from webapp import rna_structure as rs

# tRNA-Phe de S. cerevisiae (76 nt), cloverleaf canónico.
TRNA_PHE = (
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
)

HAIRPIN = "GGGGAAAACCCC"


def _require_vienna():
    if rs._import_vienna() is None:
        pytest.skip("ViennaRNA/RNAfold no está instalado")


def _assert_valid(result: rs.RnaStructure) -> None:
    assert len(result.dot_bracket) == len(result.sequence)
    assert set(result.dot_bracket) <= {".", "(", ")"}
    pairs = rs.pairs_from_dot_bracket(result.dot_bracket, result.sequence)
    assert pairs == result.pairs
    assert result.score == len(pairs)
    assert (result.source == "not_available") or result.available
    for i, j in result.pairs:
        assert j - i - 1 >= rs.MIN_HAIRPIN


def _hairpin_loops(dot_bracket: str) -> list[tuple[int, int]]:
    """Pares cuyo interior está completamente sin emparejar."""
    return [
        (i, j)
        for i, j in rs.pairs_from_dot_bracket(dot_bracket)
        if all(dot_bracket[k] == "." for k in range(i + 1, j))
    ]


def test_normalize_accepts_dna_and_whitespace():
    assert rs.normalize_sequence("  ggg gaa\naac ccc\n") == "GGGGAAAACCCC"
    # ADN -> ARN, sin espacios internos.
    assert rs.normalize_sequence("GGGGAAAACCCC") == "GGGGAAAACCCC"


def test_vienna_known_hairpin():
    _require_vienna()
    result = rs.predict_structure(HAIRPIN, backend="vienna")
    assert result.source == "vienna"
    assert result.dot_bracket == "((((....))))"
    assert result.pairs == ((0, 11), (1, 10), (2, 9), (3, 8))
    assert result.score == 4
    assert result.energy_kcal_mol is not None
    assert result.energy_kcal_mol < 0
    assert result.reason == ""
    _assert_valid(result)


def test_vienna_trna_known_is_deterministic_and_has_three_loops():
    _require_vienna()
    first = rs.predict_structure(TRNA_PHE, backend="vienna")
    second = rs.predict_structure(TRNA_PHE, backend="vienna")
    assert first == second
    assert first.source == "vienna"
    assert first.sequence == TRNA_PHE
    assert first.score == len(first.pairs)
    assert first.score >= 15
    assert len(_hairpin_loops(first.dot_bracket)) >= 3
    assert first.energy_kcal_mol is not None and first.energy_kcal_mol < 0
    _assert_valid(first)


def test_nussinov_known_hairpin():
    result = rs.predict_structure(HAIRPIN, backend="nussinov")
    assert result.source == "nussinov"
    assert result.dot_bracket == "((((....))))"
    assert result.pairs == ((0, 11), (1, 10), (2, 9), (3, 8))
    assert result.score == 4
    assert result.energy_kcal_mol is None
    assert "no una energía" in result.reason
    _assert_valid(result)


def test_nussinov_is_deterministic():
    a = rs.predict_structure(TRNA_PHE, backend="nussinov")
    b = rs.predict_structure(TRNA_PHE, backend="nussinov")
    assert a == b
    assert a.source == "nussinov"
    _assert_valid(a)


def test_explicit_nussinov_does_not_depend_on_fallback_flag():
    result = rs.predict_structure(HAIRPIN, backend="nussinov", allow_fallback=False)
    assert result.source == "nussinov"
    assert result.dot_bracket == "((((....))))"


def test_nussinov_no_pairs_for_poly_a():
    result = rs.predict_structure("AAAAAAAA", backend="nussinov")
    assert result.dot_bracket == "........"
    assert result.pairs == ()
    assert result.score == 0
    assert result.energy_kcal_mol is None


def test_nussinov_short_sequence_cannot_form_loop():
    # Con tallo-bucle mínimo 3, 4 nt no pueden formar ningún par.
    result = rs.predict_structure("GCGC", backend="nussinov")
    assert result.dot_bracket == "...."
    assert result.pairs == ()
    assert result.score == 0


def test_auto_without_vienna_uses_nussinov(monkeypatch):
    monkeypatch.setattr(rs, "_import_vienna", lambda: None)
    result = rs.predict_structure(HAIRPIN)
    assert result.source == "nussinov"
    assert result.dot_bracket == "((((....))))"
    assert result.energy_kcal_mol is None
    assert "no está instalado" in result.reason
    _assert_valid(result)


def test_without_tool_and_no_fallback_is_not_available(monkeypatch):
    monkeypatch.setattr(rs, "_import_vienna", lambda: None)
    result = rs.predict_structure(HAIRPIN, backend="vienna", allow_fallback=False)
    assert result.source == "not_available"
    assert result.dot_bracket == ""
    assert result.pairs == ()
    assert result.score is None
    assert result.energy_kcal_mol is None
    assert "no está instalado" in result.reason
    assert not result.available


def test_auto_without_tool_and_no_fallback_is_not_available(monkeypatch):
    monkeypatch.setattr(rs, "_import_vienna", lambda: None)
    result = rs.predict_structure(HAIRPIN, allow_fallback=False)
    assert result.source == "not_available"
    assert result.dot_bracket == ""
    assert result.pairs == ()
    assert "desactivado" in result.reason


def test_vienna_fold_failure_falls_back_or_declares_not_available(monkeypatch):
    class Broken:
        @staticmethod
        def fold(_sequence):
            raise RuntimeError("modelo roto")

    monkeypatch.setattr(rs, "_import_vienna", lambda: Broken)
    fallback = rs.predict_structure(HAIRPIN, backend="vienna")
    assert fallback.source == "nussinov"
    assert "falló" in fallback.reason

    unavailable = rs.predict_structure(HAIRPIN, backend="vienna", allow_fallback=False)
    assert unavailable.source == "not_available"
    assert unavailable.dot_bracket == ""
    assert "falló" in unavailable.reason


def test_empty_sequence_is_not_available():
    for value in ("", "   \n\t"):
        result = rs.predict_structure(value)
        assert result.source == "not_available"
        assert result.dot_bracket == ""
        assert result.reason == "secuencia vacía"


@pytest.mark.parametrize("sequence", ["ACGN", "ACGU R", "ACGU-", "ACGU*", "ACGU1"])
def test_ambiguous_or_invalid_bases_are_rejected(sequence):
    with pytest.raises(ValueError, match="base no soportada"):
        rs.predict_structure(sequence)


def test_non_string_sequence_is_rejected():
    with pytest.raises(TypeError):
        rs.predict_structure(None)  # type: ignore[arg-type]


def test_invalid_backend_is_rejected():
    with pytest.raises(ValueError, match="backend desconocido"):
        rs.predict_structure(HAIRPIN, backend="nupack")


def test_hard_length_limit():
    with pytest.raises(ValueError, match="límite de la API"):
        rs.predict_structure("A" * (rs.MAX_SEQUENCE_LENGTH + 1), backend="nussinov")


def test_fallback_length_limit_is_not_available(monkeypatch):
    monkeypatch.setattr(rs, "_import_vienna", lambda: None)
    result = rs.predict_structure("A" * (rs.FALLBACK_MAX_LENGTH + 1))
    assert result.source == "not_available"
    assert result.dot_bracket == ""
    assert "solo cubre" in result.reason


def test_pairs_from_dot_bracket_roundtrip():
    dot = "((..)).."
    assert rs.pairs_from_dot_bracket(dot) == ((0, 5), (1, 4))
    assert rs.dot_bracket_from_pairs(len(dot), [(0, 5), (1, 4)]) == dot


def test_pairs_from_dot_bracket_rejects_malformed():
    with pytest.raises(ValueError, match="sin '\\('"):
        rs.pairs_from_dot_bracket("..)")
    with pytest.raises(ValueError, match="sin cerrar"):
        rs.pairs_from_dot_bracket("((..")
    with pytest.raises(ValueError, match="no soportado"):
        rs.pairs_from_dot_bracket("[[..]]")
    with pytest.raises(ValueError, match="longitud distinta"):
        rs.pairs_from_dot_bracket("()", "ACG")
    with pytest.raises(ValueError, match="no canónico"):
        rs.pairs_from_dot_bracket("()", "AA")


def test_dot_bracket_from_pairs_rejects_crossing_and_duplicates():
    with pytest.raises(ValueError, match="pseudoknot"):
        rs.dot_bracket_from_pairs(8, [(0, 5), (1, 6)])
    with pytest.raises(ValueError, match="repetida"):
        rs.dot_bracket_from_pairs(8, [(0, 5), (0, 6)])
    with pytest.raises(ValueError, match="fuera de rango"):
        rs.dot_bracket_from_pairs(4, [(0, 4)])


def test_as_dict_is_json_serializable_and_has_expected_keys():
    result = rs.predict_structure(HAIRPIN, backend="nussinov")
    payload = result.as_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert set(payload) == {
        "sequence",
        "dot_bracket",
        "energy_kcal_mol",
        "score",
        "source",
        "pairs",
        "reason",
    }
    assert payload["pairs"] == [[0, 11], [1, 10], [2, 9], [3, 8]]
