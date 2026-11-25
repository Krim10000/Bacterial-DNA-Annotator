"""Tests de ``webapp/pae_complex.py`` con PAE y mapeos sintéticos.

No usan red, GPU ni datos descargados. Cubren:

- resumen por par de subunidades y ambas direcciones (A->B y B->A);
- varias cadenas por subunidad (homómeros) y agregación ponderada;
- estados ``ok``, ``not_available`` e ``invalid``;
- ausencia de PAE, cadenas ausentes y mapeos explícitos;
- JSON inválido, matriz no NxN, cadenas desconocidas o solapadas;
- rechazo de la resolución por índice ordinal para no inventar cadenas.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from webapp.pae_complex import (
    LIMITATION,
    NO_CHAIN_REASON,
    NO_MAPPING_REASON,
    STATUSES,
    STATUS_INVALID,
    STATUS_NOT_AVAILABLE,
    STATUS_OK,
    analyze_complex,
    interface_pair,
    subunits_from_chain_subunits,
)

# Matriz 4x4 asimétrica: A(0,1) -> B(2,3) con PAE 9; B -> A con PAE 1.
INTERFACE = [
    [0, 2, 9, 9],
    [2, 0, 9, 9],
    [1, 1, 0, 3],
    [1, 1, 3, 0],
]

CHAINS_3 = {"A": [0, 1], "B": [2, 3], "C": [4, 5]}


def three_chain_matrix() -> list[list[float]]:
    """Matriz 6x6 con bloques A/B/C de PAE distinto y asimétrico."""
    values = {
        ("A", "B"): 2.0,
        ("B", "A"): 8.0,
        ("A", "C"): 4.0,
        ("C", "A"): 6.0,
        ("B", "C"): 3.0,
        ("C", "B"): 7.0,
    }
    matrix = [[0.0] * 6 for _ in range(6)]
    for (chain_a, chain_b), value in values.items():
        for i in CHAINS_3[chain_a]:
            for j in CHAINS_3[chain_b]:
                matrix[i][j] = value
    return matrix


def write_json(tmp_path: Path, name: str, payload) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def pair_for(result: dict, subunit_a: str, subunit_b: str) -> dict:
    for pair in result["pairs"]:
        if {pair["subunit_a"], pair["subunit_b"]} == {subunit_a, subunit_b}:
            return pair
    raise AssertionError(f"no está el par {subunit_a}/{subunit_b}")


def assert_without_values(pair: dict) -> None:
    assert pair["mean_pae"] is None
    assert pair["max_pae"] is None
    assert pair["n_pairs"] == 0
    for key in ("a_to_b", "b_to_a"):
        direction = pair[key]
        assert direction["mean_pae"] is None
        assert direction["max_pae"] is None
        assert direction["n_pairs"] == 0
        for chain_pair in direction["chain_pairs"]:
            if chain_pair["status"] != STATUS_OK:
                assert chain_pair["mean_pae"] is None
                assert chain_pair["max_pae"] is None


# ---------------------------------------------------------------------------
# Cálculo direccional y resumen por par
# ---------------------------------------------------------------------------
def test_three_subunits_directional_means_and_pair_summaries():
    result = analyze_complex(
        three_chain_matrix(),
        {"Alpha": ["A"], "Beta": ["B"], "Gamma": ["C"]},
        chains=CHAINS_3,
    )

    assert result["status"] == STATUS_OK
    assert result["ok"] is True
    assert result["available"] is True
    assert result["reason"] == "ok"
    assert result["n_subunits"] == 3
    assert len(result["pairs"]) == 3
    assert result["subunits"] == {
        "Alpha": ["A"],
        "Beta": ["B"],
        "Gamma": ["C"],
    }

    alpha_beta = pair_for(result, "Alpha", "Beta")
    assert alpha_beta["status"] == STATUS_OK
    assert alpha_beta["a_to_b"]["mean_pae"] == pytest.approx(2.0)
    assert alpha_beta["a_to_b"]["max_pae"] == pytest.approx(2.0)
    assert alpha_beta["b_to_a"]["mean_pae"] == pytest.approx(8.0)
    assert alpha_beta["b_to_a"]["max_pae"] == pytest.approx(8.0)
    # Media combinada ponderada: 4 celdas a 2 y 4 celdas a 8.
    assert alpha_beta["mean_pae"] == pytest.approx(5.0)
    assert alpha_beta["max_pae"] == pytest.approx(8.0)
    assert alpha_beta["n_pairs"] == 8

    beta_gamma = pair_for(result, "Beta", "Gamma")
    assert beta_gamma["a_to_b"]["mean_pae"] == pytest.approx(3.0)
    assert beta_gamma["b_to_a"]["mean_pae"] == pytest.approx(7.0)
    assert beta_gamma["mean_pae"] == pytest.approx(5.0)
    assert beta_gamma["max_pae"] == pytest.approx(7.0)

    summary = result["summary"]
    assert summary["n_pairs"] == 3
    assert summary["n_ok"] == 3
    assert summary["n_not_available"] == 0
    assert summary["n_invalid"] == 0
    assert summary["n_residue_pairs"] == 24
    # Seis bloques de 4 celdas con medias 2, 8, 4, 6, 3 y 7.
    assert summary["mean_pae"] == pytest.approx(5.0)
    assert summary["max_pae"] == pytest.approx(8.0)


def test_direction_a_to_b_differs_from_b_to_a():
    result = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    pair = result["pairs"][0]
    assert pair["a_to_b"]["mean_pae"] == pytest.approx(9.0)
    assert pair["b_to_a"]["mean_pae"] == pytest.approx(1.0)
    assert pair["a_to_b"]["mean_pae"] != pair["b_to_a"]["mean_pae"]


def test_limitation_mentions_interaction_and_stoichiometry():
    assert "no demuestra interacción" in LIMITATION
    assert "estequiometría" in LIMITATION
    assert "no identifica qué residuos" in LIMITATION
    result = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["limitation"] == LIMITATION


def test_multi_chain_subunit_pools_cross_chain_blocks():
    result = analyze_complex(
        three_chain_matrix(),
        {"Homomer": ["A", "B"], "Other": ["C"]},
        chains=CHAINS_3,
    )
    assert result["status"] == STATUS_OK
    pair = result["pairs"][0]

    a_to_b = pair["a_to_b"]
    assert a_to_b["status"] == STATUS_OK
    assert a_to_b["n_chain_pairs"] == 2
    assert a_to_b["expected_chain_pairs"] == 2
    assert a_to_b["n_pairs"] == 8
    # A->C = 4 (4 celdas) y B->C = 3 (4 celdas).
    assert a_to_b["mean_pae"] == pytest.approx(3.5)
    assert a_to_b["max_pae"] == pytest.approx(4.0)
    assert [cp["resolved_chain_a"] for cp in a_to_b["chain_pairs"]] == ["A", "B"]

    b_to_a = pair["b_to_a"]
    assert b_to_a["mean_pae"] == pytest.approx(6.5)
    assert b_to_a["max_pae"] == pytest.approx(7.0)
    assert pair["mean_pae"] == pytest.approx(5.0)
    assert pair["max_pae"] == pytest.approx(7.0)
    assert pair["n_pairs"] == 16


def test_partial_chain_failure_does_not_publish_incomplete_mean():
    # Z no existe: A->C sí se puede calcular, Z->C no. La dirección debe
    # quedar invalid y sin media agregada, aunque el detalle conserve A->C.
    result = analyze_complex(
        three_chain_matrix(),
        {"Homomer": ["A", "Z"], "Other": ["C"]},
        chains=CHAINS_3,
    )
    assert result["status"] == STATUS_INVALID
    pair = result["pairs"][0]
    assert pair["status"] == STATUS_INVALID
    assert pair["mean_pae"] is None
    a_to_b = pair["a_to_b"]
    assert a_to_b["status"] == STATUS_INVALID
    assert a_to_b["mean_pae"] is None
    assert a_to_b["n_pairs"] == 0
    assert a_to_b["n_chain_pairs"] == 2
    assert a_to_b["n_chain_pairs_ok"] == 1
    assert [cp["status"] for cp in a_to_b["chain_pairs"]] == [
        STATUS_OK,
        STATUS_INVALID,
    ]
    assert a_to_b["chain_pairs"][0]["mean_pae"] == pytest.approx(4.0)
    assert "no encontrada" in a_to_b["chain_pairs"][1]["reason"]


def test_mixed_pairs_summary_only_uses_ok_pairs():
    result = analyze_complex(
        three_chain_matrix(),
        {"Alpha": ["A"], "Beta": ["B"], "Gamma": ["Z"]},
        chains=CHAINS_3,
    )
    assert result["status"] == STATUS_INVALID
    assert result["summary"]["n_pairs"] == 3
    assert result["summary"]["n_ok"] == 1
    assert result["summary"]["n_invalid"] == 2
    # Solo el par Alpha/Beta aporta 8 pares de residuos.
    assert result["summary"]["n_residue_pairs"] == 8
    assert result["summary"]["mean_pae"] == pytest.approx(5.0)
    assert result["summary"]["max_pae"] == pytest.approx(8.0)


# ---------------------------------------------------------------------------
# Formatos de PAE y cadenas
# ---------------------------------------------------------------------------
def test_path_with_token_chain_ids_matches_matrix(tmp_path: Path):
    path = write_json(
        tmp_path,
        "complejo.pae.json",
        {
            "pae": three_chain_matrix(),
            "token_chain_ids": ["A", "A", "B", "B", "C", "C"],
            "ptm": 0.8,
            "iptm": 0.7,
        },
    )
    mapping = {"Alpha": ["A"], "Beta": ["B"], "Gamma": ["C"]}
    from_path = analyze_complex(path, mapping)
    from_matrix = analyze_complex(three_chain_matrix(), mapping, chains=CHAINS_3)

    assert from_path["status"] == STATUS_OK
    assert from_path["source_kind"] == "file"
    assert from_path["summary"] == from_matrix["summary"]
    assert from_path["pairs"][0]["mean_pae"] == pytest.approx(5.0)


def test_afdb_monomer_without_chains_is_not_available(tmp_path: Path):
    path = write_json(
        tmp_path,
        "AF-P1.pae.json",
        {"predicted_aligned_error": INTERFACE, "max_predicted_aligned_error": 31},
    )
    result = analyze_complex(path, {"A": ["A"], "B": ["B"]})

    assert result["status"] == STATUS_NOT_AVAILABLE
    assert result["ok"] is False
    assert result["available"] is True
    assert "no declara cadenas" in result["reason"]
    assert len(result["pairs"]) == 1
    pair = result["pairs"][0]
    assert pair["status"] == STATUS_NOT_AVAILABLE
    assert pair["a_to_b"]["reason"] == NO_CHAIN_REASON
    assert pair["a_to_b"]["expected_chain_pairs"] == 1
    assert_without_values(pair)
    assert result["summary"]["n_ok"] == 0
    assert result["summary"]["mean_pae"] is None


def test_explicit_chains_unlock_afdb_format(tmp_path: Path):
    path = write_json(
        tmp_path,
        "AF-P1.pae.json",
        {"predicted_aligned_error": INTERFACE},
    )
    result = analyze_complex(
        path,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_OK
    assert result["pairs"][0]["a_to_b"]["mean_pae"] == pytest.approx(9.0)
    assert result["pairs"][0]["b_to_a"]["mean_pae"] == pytest.approx(1.0)


def test_explicit_chains_series_is_accepted():
    result = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["B"]},
        chains=["A", "A", "B", "B"],
    )
    assert result["status"] == STATUS_OK
    assert result["pairs"][0]["mean_pae"] == pytest.approx(5.0)


def test_explicit_chains_override_file_metadata(tmp_path: Path):
    path = write_json(
        tmp_path,
        "metadatos.pae.json",
        {"pae": INTERFACE, "token_chain_ids": ["X", "X", "Y", "Y"]},
    )
    result = analyze_complex(
        path,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_OK
    assert result["pairs"][0]["a_to_b"]["chain_pairs"][0]["resolved_chain_a"] == "A"


# ---------------------------------------------------------------------------
# Ausencias y errores
# ---------------------------------------------------------------------------
def test_missing_file_not_available_and_no_values(tmp_path: Path):
    result = analyze_complex(
        tmp_path / "no_existe.pae.json",
        {"A": ["A"], "B": ["B"]},
    )
    assert result["status"] == STATUS_NOT_AVAILABLE
    assert result["ok"] is False
    assert result["available"] is False
    assert "ausente" in result["reason"]
    assert len(result["pairs"]) == 1
    assert_without_values(result["pairs"][0])
    assert result["summary"]["n_not_available"] == 1
    assert result["summary"]["mean_pae"] is None


def test_missing_mapping_not_available():
    result = analyze_complex(INTERFACE, None)
    assert result["status"] == STATUS_NOT_AVAILABLE
    assert result["reason"] == NO_MAPPING_REASON
    assert result["pairs"] == []
    assert result["n_subunits"] == 0
    assert result["summary"]["mean_pae"] is None
    assert "estequiometría" in result["limitation"]


def test_invalid_json_file_is_invalid_not_not_available(tmp_path: Path):
    path = tmp_path / "roto.pae.json"
    path.write_text("{esto no es json", encoding="utf-8")
    result = analyze_complex(path, {"A": ["A"], "B": ["B"]})

    assert result["status"] == STATUS_INVALID
    assert result["available"] is True
    assert "PAE no utilizable" in result["reason"]
    assert "JSON inválido" in result["reason"]
    assert result["pairs"][0]["status"] == STATUS_INVALID
    assert_without_values(result["pairs"][0])


def test_non_square_matrix_file_is_invalid(tmp_path: Path):
    path = write_json(
        tmp_path,
        "cuadrada.pae.json",
        {"pae": [[0, 1, 2], [1, 0]]},
    )
    result = analyze_complex(path, {"A": ["A"], "B": ["B"]})
    assert result["status"] == STATUS_INVALID
    assert "columnas" in result["reason"]
    assert_without_values(result["pairs"][0])


def test_non_square_matrix_source_is_invalid():
    result = analyze_complex(
        [[0, 1], [1]],
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0], "B": [1]},
    )
    assert result["status"] == STATUS_INVALID
    assert "PAE no utilizable" in result["reason"]
    assert result["summary"]["mean_pae"] is None


@pytest.mark.parametrize(
    "mapping, needle",
    [
        ({"A": ["A"]}, "al menos 2 subunidades"),
        ({"A": [], "B": ["B"]}, "sin cadenas"),
        ({"A": ["A"], "B": ["A"]}, "dos subunidades"),
        ({"A": ["A", "A"], "B": ["B"]}, "cadena duplicada"),
        ({"A": ["A"], "B": [""]}, "vacío"),
        ({"A": ["A"], "B": [True]}, "booleano"),
        ({"A": ["A"], "B": {"chains": ["B"]}}, "se esperaba una cadena"),
        ({"A": ["A"], "B": ["B"], " A ": ["C"]}, "subunidad duplicada"),
        ("A", "debe ser un objeto"),
        ({}, "al menos 2 subunidades"),
    ],
)
def test_malformed_mapping_is_invalid(mapping, needle: str):
    result = analyze_complex(
        INTERFACE, mapping, chains={"A": [0, 1], "B": [2, 3]}
    )
    assert result["status"] == STATUS_INVALID
    assert "mapeo de subunidades inválido" in result["reason"]
    assert needle in result["reason"]
    assert result["pairs"] == []
    assert result["summary"]["mean_pae"] is None


def test_unknown_chain_is_invalid():
    result = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["Z"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_INVALID
    assert "no encontrada" in result["reason"]
    pair = result["pairs"][0]
    assert pair["status"] == STATUS_INVALID
    assert_without_values(pair)


def test_numeric_chain_id_order_fallback_is_rejected():
    # "0"/"1" no son cadenas declaradas: pae_utils los resolvería por orden a
    # A/B. La capa de complejos debe rechazarlo para no inventar cadenas.
    result = analyze_complex(
        INTERFACE,
        {"X": ["0"], "Y": ["1"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_INVALID
    assert "se rechaza la resolución por orden" in result["reason"]
    pair = result["pairs"][0]
    assert pair["status"] == STATUS_INVALID
    assert_without_values(pair)


def test_numeric_chain_ids_that_exist_are_accepted():
    result = analyze_complex(
        INTERFACE,
        {"X": ["0"], "Y": ["1"]},
        chains={"0": [0, 1], "1": [2, 3]},
    )
    assert result["status"] == STATUS_OK
    assert result["pairs"][0]["a_to_b"]["mean_pae"] == pytest.approx(9.0)
    assert result["pairs"][0]["a_to_b"]["chain_pairs"][0]["resolved_chain_a"] == "0"


def test_chain_ids_resolve_case_insensitively_when_unique():
    result = analyze_complex(
        INTERFACE,
        {"X": ["a"], "Y": ["b"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_OK
    pair = result["pairs"][0]
    assert pair["a_to_b"]["chain_pairs"][0]["resolved_chain_a"] == "A"
    assert pair["b_to_a"]["chain_pairs"][0]["resolved_chain_a"] == "B"


def test_overlapping_chain_indices_are_invalid():
    result = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [1, 2]},
    )
    assert result["status"] == STATUS_INVALID
    assert "ya asignados" in result["reason"]
    assert_without_values(result["pairs"][0])


def test_all_statuses_belong_to_declared_set(tmp_path: Path):
    assert set(STATUSES) == {STATUS_OK, STATUS_NOT_AVAILABLE, STATUS_INVALID}
    ok = analyze_complex(
        INTERFACE,
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    missing = analyze_complex(
        tmp_path / "inexistente.pae.json", {"A": ["A"], "B": ["B"]}
    )
    invalid = analyze_complex(INTERFACE, {"A": ["A"]})
    for result in (ok, missing, invalid):
        assert result["status"] in STATUSES


# ---------------------------------------------------------------------------
# interface_pair
# ---------------------------------------------------------------------------
def test_interface_pair_returns_both_directions():
    result = interface_pair(
        INTERFACE,
        "A",
        "B",
        {"A": ["A"], "B": ["B"]},
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["status"] == STATUS_OK
    assert result["a_to_b"]["mean_pae"] == pytest.approx(9.0)
    assert result["b_to_a"]["mean_pae"] == pytest.approx(1.0)
    assert result["mean_pae"] == pytest.approx(5.0)
    assert result["max_pae"] == pytest.approx(9.0)
    assert result["source_kind"] == "matrix"
    assert result["subunits"] == {"A": ["A"], "B": ["B"]}


def test_interface_pair_unknown_or_same_subunit_is_invalid():
    mapping = {"A": ["A"], "B": ["B"]}
    chains = {"A": [0, 1], "B": [2, 3]}

    unknown = interface_pair(INTERFACE, "A", "Z", mapping, chains=chains)
    assert unknown["status"] == STATUS_INVALID
    assert "no encontrada" in unknown["reason"]
    assert unknown["a_to_b"]["status"] == STATUS_INVALID
    assert unknown["mean_pae"] is None

    same = interface_pair(INTERFACE, "A", "A", mapping, chains=chains)
    assert same["status"] == STATUS_INVALID
    assert "misma subunidad" in same["reason"]
    assert same["mean_pae"] is None


def test_interface_pair_without_mapping_is_not_available():
    result = interface_pair(INTERFACE, "A", "B", None)
    assert result["status"] == STATUS_NOT_AVAILABLE
    assert result["reason"] == NO_MAPPING_REASON
    assert result["a_to_b"]["status"] == STATUS_NOT_AVAILABLE
    assert result["mean_pae"] is None


# ---------------------------------------------------------------------------
# Inversión del mapeo cadena -> subunidad
# ---------------------------------------------------------------------------
def test_subunits_from_chain_subunits_groups_homomer():
    grouped = subunits_from_chain_subunits(
        {"A": "RpoA", "B": "RpoB", "C": "RpoA"}
    )
    assert grouped == {"RpoA": ["A", "C"], "RpoB": ["B"]}

    result = analyze_complex(
        three_chain_matrix(),
        grouped,
        chains=CHAINS_3,
    )
    assert result["status"] == STATUS_OK
    pair = pair_for(result, "RpoA", "RpoB")
    assert pair["a_to_b"]["n_chain_pairs"] == 2
    # RpoA = A + C; A->B = 2 y C->B = 7, cuatro celdas cada uno.
    assert pair["a_to_b"]["mean_pae"] == pytest.approx(4.5)
    assert pair["a_to_b"]["max_pae"] == pytest.approx(7.0)


@pytest.mark.parametrize(
    "mapping, needle",
    [
        ({"A": "RpoA", " A": "RpoA"}, "duplicada"),
        ({"A": ""}, "vacía"),
        ({"A": 3}, "texto"),
        ({}, "vacío"),
        ("A", "debe ser un objeto"),
    ],
)
def test_subunits_from_chain_subunits_rejects_invalid(mapping, needle: str):
    with pytest.raises(ValueError) as excinfo:
        subunits_from_chain_subunits(mapping)
    assert needle in str(excinfo.value)
