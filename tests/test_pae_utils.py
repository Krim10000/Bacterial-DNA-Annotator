"""Tests de ``webapp/pae_utils.py`` con JSON PAE sintéticos.

No usan red, GPU ni datos descargados. Cubren:
- formato AFDB (``predicted_aligned_error``);
- AlphaFold 3 / Server / ColabFold (``pae`` + ``token_chain_ids``);
- matriz directa y lista de modelos;
- resúmenes por residuo/par e interfaz direccional;
- JSON inválido, matriz no NxN, archivo ausente y metadatos de cadena.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from webapp.pae_utils import (
    LIMITATION,
    interface_pae,
    load_matrix,
    load_pae,
    max_pae,
    mean_pae,
    pair_summary,
    residue_summary,
    summarize,
)

# Matriz 3x3 simétrica y pequeña para cuentas exactas.
SMALL = [
    [0, 2, 4],
    [2, 0, 6],
    [4, 6, 0],
]

# Matriz 4x4 asimétrica: bloque A (0,1) -> B (2,3) con PAE alto,
# B -> A con PAE bajo. Sirve para probar la dirección de la interfaz.
INTERFACE = [
    [0, 2, 9, 9],
    [2, 0, 9, 9],
    [1, 1, 0, 3],
    [1, 1, 3, 0],
]


def write_json(tmp_path: Path, name: str, payload) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# load_matrix
# ---------------------------------------------------------------------------
def test_load_matrix_afdb_format(tmp_path: Path):
    path = write_json(
        tmp_path,
        "AF-P1.pae.json",
        {
            "predicted_aligned_error": SMALL,
            "max_predicted_aligned_error": 31,
        },
    )
    matrix = load_matrix(path)
    assert matrix == [[0.0, 2.0, 4.0], [2.0, 0.0, 6.0], [4.0, 6.0, 0.0]]
    assert all(isinstance(value, float) for row in matrix for value in row)


@pytest.mark.parametrize(
    "payload",
    [
        {"pae": SMALL},
        {"pae_matrix": SMALL},
        {"aligned_error": SMALL},
        {"matrix": SMALL},
        SMALL,
        [{"pae": SMALL}],
    ],
)
def test_load_matrix_accepts_variants(tmp_path: Path, payload):
    path = write_json(tmp_path, "variant.pae.json", payload)
    assert load_matrix(path) == [[0.0, 2.0, 4.0], [2.0, 0.0, 6.0], [4.0, 6.0, 0.0]]


def test_load_matrix_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="no encontrado"):
        load_matrix(tmp_path / "no_existe.pae.json")


def test_load_matrix_invalid_json_raises(tmp_path: Path):
    path = tmp_path / "roto.pae.json"
    path.write_text("{esto no es json", encoding="utf-8")
    with pytest.raises(ValueError) as excinfo:
        load_matrix(path)
    assert "JSON inválido" in str(excinfo.value)


@pytest.mark.parametrize(
    "payload, needle",
    [
        ({"pae": []}, "vacía"),
        ({"pae": [[0, 1], [1]]}, "columnas"),
        ({"pae": [[0, 1, 2], [1, 0, 2]]}, "columnas"),
        ({"pae": [[0, 1], 2]}, "fila 2"),
        ([1, 2, 3], "fila 1"),
        ({"pae": [[0, "x"], [1, 0]]}, "número"),
        ({"pae": [[0, -1], [1, 0]]}, "negativo"),
        ({"pae": [[0, float("nan")], [1, 0]]}, "no finito"),
        ({"no_pae": []}, "sin matriz"),
        ("texto", "objeto o una lista"),
    ],
)
def test_load_matrix_rejects_invalid_payloads(tmp_path: Path, payload, needle: str):
    path = write_json(tmp_path, "invalido.pae.json", payload)
    with pytest.raises(ValueError) as excinfo:
        load_matrix(path)
    assert needle in str(excinfo.value)


def test_load_matrix_does_not_silently_truncate(tmp_path: Path):
    path = write_json(
        tmp_path,
        "truncada.pae.json",
        {"predicted_aligned_error": [[0, 1, 2], [1, 0], [2, 1, 0]]},
    )
    with pytest.raises(ValueError):
        load_matrix(path)


# ---------------------------------------------------------------------------
# summarize
# ---------------------------------------------------------------------------
def test_summarize_valid(tmp_path: Path):
    path = write_json(
        tmp_path,
        "valido.pae.json",
        {"predicted_aligned_error": SMALL, "max_predicted_aligned_error": 31},
    )
    summary = summarize(path)
    assert set(summary) == {
        "available", "ok", "n", "mean", "max", "reason", "limitation"
    }
    assert summary["available"] is True
    assert summary["ok"] is True
    assert summary["n"] == 3
    assert summary["mean"] == pytest.approx(24 / 9)
    assert summary["max"] == 6
    assert summary["reason"] == "ok"
    assert summary["limitation"] == LIMITATION
    assert "no demuestra interacción" in summary["limitation"]


def test_summarize_missing_does_not_invent(tmp_path: Path):
    summary = summarize(tmp_path / "ausente.pae.json")
    assert summary["available"] is False
    assert summary["ok"] is False
    assert summary["n"] == 0
    assert summary["mean"] is None
    assert summary["max"] is None
    assert "ausente" in summary["reason"]
    assert summary["limitation"] == LIMITATION


def test_summarize_invalid_json(tmp_path: Path):
    path = tmp_path / "roto.pae.json"
    path.write_text("{", encoding="utf-8")
    summary = summarize(path)
    assert summary["available"] is True  # el archivo existe
    assert summary["ok"] is False
    assert summary["n"] == 0
    assert summary["mean"] is None
    assert summary["max"] is None
    assert "JSON inválido" in summary["reason"]


def test_summarize_non_square_matrix(tmp_path: Path):
    path = write_json(
        tmp_path, "cuadrada.pae.json", {"pae": [[0, 1, 2], [3, 4, 5]]}
    )
    summary = summarize(path)
    assert summary["available"] is True
    assert summary["ok"] is False
    assert summary["n"] == 0
    assert summary["mean"] is None
    assert "columnas" in summary["reason"]


def test_summarize_without_pae_key(tmp_path: Path):
    path = write_json(tmp_path, "otro.json", {"plddt": [90, 80, 70]})
    summary = summarize(path)
    assert summary["available"] is True
    assert summary["ok"] is False
    assert "sin matriz PAE" in summary["reason"]


# ---------------------------------------------------------------------------
# load_pae / metadatos
# ---------------------------------------------------------------------------
def test_load_pae_afdb_has_no_chain_metadata(tmp_path: Path):
    path = write_json(
        tmp_path,
        "AF-P1.pae.json",
        {"predicted_aligned_error": SMALL, "max_predicted_aligned_error": 31},
    )
    info = load_pae(path)
    assert info["n"] == 3
    assert info["source_key"] == "predicted_aligned_error"
    assert info["model_index"] is None
    assert info["chains"] == {}
    assert info["chain_order"] == []
    assert "no declara cadenas" in info["chain_reason"]
    assert info["max_predicted_aligned_error"] == 31
    assert info["matrix"] == load_matrix(path)


def test_load_pae_af3_metadata_and_chains(tmp_path: Path):
    path = write_json(
        tmp_path,
        "af3.pae.json",
        {
            "pae": INTERFACE,
            "token_chain_ids": ["A", "A", "B", "B"],
            "token_plddts": [91, 92, 93, 94],
            "ptm": 0.81,
            "iptm": 0.72,
        },
    )
    info = load_pae(path)
    assert info["source_key"] == "pae"
    assert info["chains"] == {"A": [0, 1], "B": [2, 3]}
    assert info["chain_order"] == ["A", "B"]
    assert info["chain_reason"] == ""
    assert info["plddt"] == [91.0, 92.0, 93.0, 94.0]
    assert info["ptm"] == 0.81
    assert info["iptm"] == 0.72


def test_load_pae_uses_first_model_with_pae(tmp_path: Path):
    path = write_json(
        tmp_path,
        "colabfold.pae.json",
        [
            {"plddt": [99, 99, 99]},
            {"pae": SMALL, "plddt": [80, 81, 82]},
        ],
    )
    info = load_pae(path)
    assert info["model_index"] == 1
    assert info["plddt"] == [80.0, 81.0, 82.0]
    assert info["matrix"] == [[0.0, 2.0, 4.0], [2.0, 0.0, 6.0], [4.0, 6.0, 0.0]]


def test_load_pae_ignores_invalid_optional_metadata(tmp_path: Path):
    path = write_json(
        tmp_path,
        "meta.pae.json",
        {"pae": SMALL, "plddt": [1, 2], "ptm": "alto", "iptm": None},
    )
    info = load_pae(path)
    assert info["plddt"] is None
    assert info["ptm"] is None
    assert info["iptm"] is None


def test_load_pae_bad_chain_series_does_not_break_matrix(tmp_path: Path):
    path = write_json(
        tmp_path,
        "cadenas_mal.pae.json",
        {"pae": [[0, 1], [1, 0]], "token_chain_ids": ["A"]},
    )
    info = load_pae(path)
    assert info["n"] == 2
    assert info["chains"] == {}
    assert "inválidos" in info["chain_reason"]
    # El resumen del PAE sigue siendo válido; lo que falta son las cadenas.
    assert summarize(path)["ok"] is True


# ---------------------------------------------------------------------------
# Estadísticas
# ---------------------------------------------------------------------------
def test_mean_max_and_reject_invalid():
    assert mean_pae(SMALL) == pytest.approx(24 / 9)
    assert max_pae(SMALL) == 6
    with pytest.raises(ValueError):
        mean_pae([[0, 1], [1]])
    with pytest.raises(ValueError):
        max_pae([[0, -1], [1, 0]])


def test_residue_summary_excludes_self():
    rows = residue_summary(SMALL)
    assert [row["residue"] for row in rows] == [1, 2, 3]
    assert rows[0]["mean_pae"] == pytest.approx(3.0)
    assert rows[0]["max_pae"] == 4
    assert rows[0]["min_pae"] == 2
    assert rows[1]["mean_pae"] == pytest.approx(4.0)
    assert rows[2]["mean_pae"] == pytest.approx(5.0)


def test_pair_summary_ordered_off_diagonal():
    summary = pair_summary(SMALL)
    assert summary["n"] == 3
    assert summary["n_pairs"] == 6
    assert summary["mean_pae"] == pytest.approx(4.0)
    assert summary["max_pae"] == 6
    assert summary["min_pae"] == 2
    assert summary["median_pae"] == pytest.approx(4.0)
    assert summary["highest"] == {"i": 2, "j": 3, "pae": 6.0}
    assert summary["lowest"] == {"i": 1, "j": 2, "pae": 2.0}


def test_single_residue_has_no_pairs():
    rows = residue_summary([[0]])
    assert rows == [
        {"residue": 1, "mean_pae": None, "max_pae": None, "min_pae": None}
    ]
    summary = pair_summary([[0]])
    assert summary["n_pairs"] == 0
    assert summary["mean_pae"] is None
    assert summary["highest"] is None
    assert summary["lowest"] is None


# ---------------------------------------------------------------------------
# interface_pae
# ---------------------------------------------------------------------------
def test_interface_pae_is_directional(tmp_path: Path):
    path = write_json(
        tmp_path,
        "interfaz.pae.json",
        {"pae": INTERFACE, "token_chain_ids": ["A", "A", "B", "B"]},
    )
    a_to_b = interface_pae(path, "A", "B")
    assert a_to_b["available"] is True
    assert a_to_b["ok"] is True
    assert a_to_b["chain_a"] == "A"
    assert a_to_b["chain_b"] == "B"
    assert a_to_b["n_a"] == 2
    assert a_to_b["n_b"] == 2
    assert a_to_b["n_pairs"] == 4
    assert a_to_b["mean"] == pytest.approx(9.0)
    assert a_to_b["max"] == 9
    assert a_to_b["median"] == pytest.approx(9.0)
    assert a_to_b["mean_pae"] == pytest.approx(9.0)
    assert a_to_b["max_pae"] == 9
    assert a_to_b["reason"] == "ok"
    assert a_to_b["limitation"] == LIMITATION

    b_to_a = interface_pae(path, "B", "A")
    assert b_to_a["ok"] is True
    assert b_to_a["mean"] == pytest.approx(1.0)
    assert b_to_a["max"] == 1


def test_interface_pae_accepts_lowercase_and_order_index(tmp_path: Path):
    path = write_json(
        tmp_path,
        "interfaz.pae.json",
        {"pae": INTERFACE, "token_chain_ids": ["A", "A", "B", "B"]},
    )
    assert interface_pae(path, "a", "b")["mean"] == pytest.approx(9.0)
    assert interface_pae(path, 0, 1)["mean"] == pytest.approx(9.0)


def test_interface_pae_explicit_chains_on_matrix():
    result = interface_pae(
        INTERFACE,
        "A",
        "B",
        chains={"A": [0, 1], "B": [2, 3]},
    )
    assert result["ok"] is True
    assert result["mean"] == pytest.approx(9.0)

    ranged = interface_pae(
        INTERFACE,
        "A",
        "B",
        chains={"A": {"start": 0, "end": 1}, "B": {"start": 2, "end": 3}},
    )
    assert ranged["ok"] is True
    assert ranged["mean"] == pytest.approx(9.0)

    listed = interface_pae(
        INTERFACE,
        "A",
        "B",
        chains=[
            {"chain_id": "A", "indices": [0, 1]},
            {"chain_id": "B", "indices": [2, 3]},
        ],
    )
    assert listed["ok"] is True
    assert listed["mean"] == pytest.approx(9.0)

    series = interface_pae(
        INTERFACE,
        "A",
        "B",
        chains=["A", "A", "B", "B"],
    )
    assert series["ok"] is True
    assert series["mean"] == pytest.approx(9.0)


def test_interface_pae_without_chain_metadata(tmp_path: Path):
    path = write_json(
        tmp_path,
        "AF-P1.pae.json",
        {"predicted_aligned_error": SMALL},
    )
    result = interface_pae(path, "A", "B")
    assert result["available"] is True
    assert result["ok"] is False
    assert result["mean"] is None
    assert result["max"] is None
    assert "no declara cadenas" in result["reason"]


def test_interface_pae_missing_file(tmp_path: Path):
    result = interface_pae(tmp_path / "no_existe.pae.json", "A", "B")
    assert result["available"] is False
    assert result["ok"] is False
    assert result["mean"] is None
    assert "ausente" in result["reason"]


def test_interface_pae_invalid_json(tmp_path: Path):
    path = tmp_path / "roto.pae.json"
    path.write_text("no json", encoding="utf-8")
    result = interface_pae(path, "A", "B")
    assert result["available"] is True
    assert result["ok"] is False
    assert "no utilizable" in result["reason"]


def test_interface_pae_unknown_and_same_chain(tmp_path: Path):
    path = write_json(
        tmp_path,
        "interfaz.pae.json",
        {"pae": INTERFACE, "token_chain_ids": ["A", "A", "B", "B"]},
    )
    unknown = interface_pae(path, "A", "Z")
    assert unknown["ok"] is False
    assert "no encontrada" in unknown["reason"]

    same = interface_pae(path, "A", "A")
    assert same["ok"] is False
    assert "misma cadena" in same["reason"]


def test_interface_pae_rejects_overlapping_chains():
    result = interface_pae(
        [[0, 1, 2], [1, 0, 2], [2, 2, 0]],
        "A",
        "B",
        chains={"A": [0, 1], "B": [1, 2]},
    )
    assert result["ok"] is False
    assert result["mean"] is None
    assert "ya asignados" in result["reason"]


def test_interface_pae_empty_matrix_source():
    result = interface_pae([], "A", "B", chains={"A": [0], "B": [1]})
    assert result["ok"] is False
    assert "vacía" in result["reason"]
