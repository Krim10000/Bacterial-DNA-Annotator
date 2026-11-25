"""Pruebas de la ficha de PAE de complejos (``webapp/complex_pae_ui.py``).

Todo es sintético: matrices PAE mínimas, JSON locales y un Streamlit falso.
No hay red, GPU, ni archivos descargados.

Cubren:

- búsqueda local por ID de complejo y por el conjunto de subunidades;
- ruta declarada (incluida relativa a ``data/`` y a la raíz del proyecto);
- ausencia de PAE, ausencia de cadenas y ausencia de mapeo;
- mapeo explícito desde la fila, un sidecar o el propio JSON del PAE;
- rechazo de mapeos inventados por igualdad de ID, orden o nombre;
- tabla con media/máximo ``A -> B`` y ``B -> A``;
- el aviso de que el PAE no demuestra interacción ni estequiometría;
- render con Streamlit falso para ``ok``, ``not_available`` e ``invalid``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import complex_pae_ui as ui


# Matriz 5x5: cadenas A (filas/columnas 0-2) y B (3-4).
# A -> B: valores 4 y 6  -> media 5.0, máximo 6.0
# B -> A: valores 8, 9, 10 -> media 9.0, máximo 10.0
MATRIX = [
    [0, 0, 0, 4, 6],
    [0, 0, 0, 4, 6],
    [0, 0, 0, 4, 6],
    [8, 9, 10, 0, 0],
    [8, 9, 10, 0, 0],
]
CHAIN_IDS = ["A", "A", "A", "B", "B"]


def _root(tmp_path: Path) -> Path:
    """Crea un ``data/`` mínimo con ``metabolism/targets``."""
    root = tmp_path / "data"
    (root / "metabolism" / "targets").mkdir(parents=True, exist_ok=True)
    return root


def _write_pae(
    path: Path,
    *,
    matrix: list[list[float]] | None = None,
    chain_ids: list[str] | None = None,
    with_chains: bool = True,
    extra: dict | None = None,
) -> Path:
    """Escribe un PAE local sintético (formato AF3/ColabFold: ``pae``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data: dict = {"pae": matrix if matrix is not None else MATRIX}
    if with_chains:
        data["token_chain_ids"] = CHAIN_IDS if chain_ids is None else chain_ids
    if extra:
        data.update(extra)
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _row(**overrides) -> dict:
    row = {
        "complex_id": "CX1",
        "name": "Complejo sintético",
        "subunits": "cyda:NC_1;cydb:NC_2",
    }
    row.update(overrides)
    return row


def _ok_mapping() -> dict:
    return {"chain_subunits": {"A": "SubA", "B": "SubB"}}


# ---------------------------------------------------------------------------
# Parseo de subunidades
# ---------------------------------------------------------------------------
def test_complex_subunits_accepts_project_string_format():
    assert ui.complex_subunits(
        {"subunits": "cyda:NC_1; cydb:NC_2"}
    ) == ["cyda", "cydb"]


def test_complex_subunits_accepts_structured_list():
    labels = ui.complex_subunits({
        "subunits": [
            {"gene": "atpA", "locus_tag": "b1", "uniprot": "P1"},
            {"uniprot": "P2"},
        ]
    })
    assert labels == ["atpA", "P2"]


def test_complex_subunits_uses_dict_keys_only():
    # Los valores de un dict de subunidades no se interpretan como cadenas.
    assert ui.complex_subunits(
        {"subunits": {"RpoA": ["A"], "RpoB": ["B"]}}
    ) == ["RpoA", "RpoB"]


def test_complex_subunits_empty_and_none():
    assert ui.complex_subunits(None) == []
    assert ui.complex_subunits({}) == []
    assert ui.complex_subunits({"subunits": "  "}) == []


def test_complex_subunits_does_not_mutate_row():
    row = {"subunits": "A;B"}
    before = dict(row)
    ui.complex_subunits(row)
    assert row == before


# ---------------------------------------------------------------------------
# Búsqueda local
# ---------------------------------------------------------------------------
def test_find_pae_by_complex_id(tmp_path):
    root = _root(tmp_path)
    path = _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")

    found = ui.find_complex_pae(_row(), root)

    assert found["status"] == ui.STATUS_OK
    assert found["path"] == path
    assert found["declared"] is False


def test_find_pae_by_sanitized_complex_id(tmp_path):
    root = _root(tmp_path)
    path = _write_pae(
        root / "metabolism" / "targets" /
        "ATP_SYNTHASE_NC_002620.2_691287.pae.json"
    )
    row = _row(
        complex_id="ATP_SYNTHASE:NC_002620.2:691287",
        subunits="atpa:NC_1;atpb:NC_2",
    )

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_OK
    assert found["path"] == path


def test_find_pae_by_all_subunits(tmp_path):
    root = _root(tmp_path)
    path = _write_pae(
        root / "metabolism" / "targets" /
        "P37494__P59653.pae.json"
    )
    row = _row(complex_id="OTRO", subunits="P37494;P59653")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_OK
    assert found["path"] == path


def test_single_subunit_filename_is_not_a_complex_pae(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "P37494.pae.json")
    row = _row(complex_id="OTRO", subunits="P37494;P59653")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert found["path"] is None


def test_ambiguous_pae_is_invalid_and_never_chosen_at_random(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "P1__P2.pae.json")
    _write_pae(root / "metabolism" / "targets" / "P1_P2.pae.json")
    row = _row(complex_id="OTRO", subunits="P1;P2")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_INVALID
    assert found["path"] is None
    assert "no se elige uno al azar" in found["reason"]


def test_declared_path_relative_to_data_root(tmp_path):
    root = _root(tmp_path)
    path = _write_pae(
        root / "metabolism" / "targets" / "custom.pae.json"
    )
    row = _row(pae="metabolism/targets/custom.pae.json")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_OK
    assert found["path"] == path
    assert found["declared"] is True


def test_declared_path_relative_to_project_root(tmp_path):
    root = _root(tmp_path)
    path = _write_pae(root / "metabolism" / "targets" / "custom.pae.json")
    row = _row(pae="data/metabolism/targets/custom.pae.json")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_OK
    assert found["path"] == path


def test_declared_path_missing_reports_not_available(tmp_path):
    root = _root(tmp_path)
    row = _row(pae="metabolism/targets/missing.pae.json")

    found = ui.find_complex_pae(row, root)
    view = ui.complex_pae_view(row, root)

    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert "no encontrada" in found["reason"]
    assert view["status"] == ui.STATUS_NOT_AVAILABLE
    assert view["pae_declared"] is True
    assert view["pae_found"] is False


def test_remote_declared_path_is_not_downloaded(tmp_path):
    root = _root(tmp_path)
    row = _row(pae="https://example.org/CX1.pae.json")

    found = ui.find_complex_pae(row, root)

    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert found["path"] is None
    assert "no descarga" in found["reason"]


def test_no_pae_in_targets_is_not_available(tmp_path):
    root = _root(tmp_path)

    found = ui.find_complex_pae(_row(), root)

    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert found["path"] is None
    assert "no hay PAE local" in found["reason"]


def test_no_targets_directory_is_not_available(tmp_path):
    root = tmp_path / "data"
    found = ui.find_complex_pae(_row(), root)
    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert "no existe el directorio" in found["reason"]


def test_no_data_root_is_not_available():
    found = ui.find_complex_pae(_row(), None)
    assert found["status"] == ui.STATUS_NOT_AVAILABLE
    assert "data_root" in found["reason"]


# ---------------------------------------------------------------------------
# Mapeo explícito y estados
# ---------------------------------------------------------------------------
def test_view_ok_with_row_chain_subunits(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")

    view = ui.complex_pae_view(_row(**_ok_mapping()), root)

    assert view["status"] == ui.STATUS_OK
    assert view["ok"] is True
    assert view["pae_found"] is True
    assert view["mapping"] == {"SubA": ["A"], "SubB": ["B"]}
    assert view["mapping_source"] == "complex_row"
    assert view["summary"]["n_ok"] == 1
    assert view["summary"]["mean_pae"] == pytest.approx(7.0)
    assert view["summary"]["max_pae"] == pytest.approx(10.0)
    assert view["summary"]["n_residue_pairs"] == 12

    pair = view["rows"][0]
    assert pair["subunit_a"] == "SubA"
    assert pair["subunit_b"] == "SubB"
    assert pair["chains_a"] == ["A"]
    assert pair["chains_b"] == ["B"]
    assert pair["mean_a_to_b"] == pytest.approx(5.0)
    assert pair["max_a_to_b"] == pytest.approx(6.0)
    assert pair["mean_b_to_a"] == pytest.approx(9.0)
    assert pair["max_b_to_a"] == pytest.approx(10.0)


def test_complex_pae_view_accepts_pandas_series(tmp_path):
    import pandas as pd

    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    row = pd.Series(_row(**_ok_mapping()))

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_OK
    assert view["complex_id"] == "CX1"
    assert view["rows"][0]["mean_a_to_b"] == pytest.approx(5.0)


def test_view_ok_with_row_subunit_chains(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    row = _row(subunit_chains={"SubA": ["A"], "SubB": ["B"]})

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_OK
    assert view["mapping_source"] == "complex_row"
    assert view["rows"][0]["mean_a_to_b"] == pytest.approx(5.0)


def test_view_ok_with_subunits_dict_as_explicit_mapping(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    row = _row(subunits={"SubA": ["A"], "SubB": ["B"]})

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_OK
    assert view["mapping"] == {"SubA": ["A"], "SubB": ["B"]}
    assert view["mapping_source"] == "complex_row"
    assert view["rows"][0]["mean_b_to_a"] == pytest.approx(9.0)


def test_view_ok_with_sidecar_mapping(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    (root / "metabolism" / "targets" / "CX1.chains.json").write_text(
        json.dumps({"chain_subunits": {"A": "SubA", "B": "SubB"}}),
        encoding="utf-8",
    )

    view = ui.complex_pae_view(_row(), root)

    assert view["status"] == ui.STATUS_OK
    assert view["mapping_source"] == "sidecar"
    assert view["mapping"] == {"SubA": ["A"], "SubB": ["B"]}


def test_view_ok_with_mapping_embedded_in_pae_json(tmp_path):
    root = _root(tmp_path)
    _write_pae(
        root / "metabolism" / "targets" / "CX1.pae.json",
        extra={"chain_subunits": {"A": "SubA", "B": "SubB"}},
    )

    view = ui.complex_pae_view(_row(), root)

    assert view["status"] == ui.STATUS_OK
    assert view["mapping_source"] == "pae_json"
    assert view["analysis"]["status"] == ui.STATUS_OK


def test_declared_row_dict_carries_path_and_mapping(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "custom.pae.json")
    row = _row(
        pae={
            "path": "metabolism/targets/custom.pae.json",
            "chain_subunits": {"A": "SubA", "B": "SubB"},
        }
    )

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_OK
    assert view["pae_path"].endswith("custom.pae.json")
    assert view["mapping_source"] == "complex_row"


def test_missing_mapping_is_not_available_and_never_inferred(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    # Las cadenas del PAE se llaman A/B y las subunidades también A/B.
    # Aun así, la igualdad de ID no es un mapeo declarado.
    row = _row(complex_id="CX1", subunits="A;B")

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_NOT_AVAILABLE
    assert view["ok"] is False
    assert view["analysis"] is None
    assert view["mapping"] is None
    assert view["declared_chains"] == ["A", "B"]
    assert "falta el mapeo cadena-subunidad" in view["reason"]
    assert "no se infiere" in view["reason"]
    assert view["rows"] == []


def test_mapping_with_empty_subunit_is_invalid(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    row = _row(chain_subunits={"A": "", "B": "SubB"})

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_INVALID
    assert "mapeo cadena-subunidad inválido" in view["reason"]


def test_unknown_chain_is_invalid(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    row = _row(chain_subunits={"A": "SubA", "C": "SubB"})

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_INVALID
    assert view["analysis"] is not None
    assert view["analysis"]["status"] == ui.STATUS_INVALID


def test_pae_without_chains_is_not_available(tmp_path):
    root = _root(tmp_path)
    _write_pae(
        root / "metabolism" / "targets" / "CX1.pae.json",
        with_chains=False,
    )

    view = ui.complex_pae_view(_row(**_ok_mapping()), root)

    assert view["status"] == ui.STATUS_NOT_AVAILABLE
    assert view["analysis"]["status"] == ui.STATUS_NOT_AVAILABLE
    assert "cadena" in view["reason"].casefold()


def test_corrupt_pae_is_invalid(tmp_path):
    root = _root(tmp_path)
    path = root / "metabolism" / "targets" / "CX1.pae.json"
    path.write_text("{no es json", encoding="utf-8")

    view = ui.complex_pae_view(_row(**_ok_mapping()), root)

    assert view["status"] == ui.STATUS_INVALID
    assert "PAE no utilizable" in view["reason"]


def test_homomer_mapping_groups_two_chains(tmp_path):
    root = _root(tmp_path)
    # Cadenas: A (SubA, 2 residuos), C (SubA, 2), B (SubB, 2).
    matrix = [[0.0] * 6 for _ in range(6)]
    for i in (0, 1):          # A -> B
        for j in (4, 5):
            matrix[i][j] = 4.0
    for i in (2, 3):          # C -> B
        for j in (4, 5):
            matrix[i][j] = 6.0
    for i in (4, 5):          # B -> A
        for j in (0, 1):
            matrix[i][j] = 8.0
    for i in (4, 5):          # B -> C
        for j in (2, 3):
            matrix[i][j] = 10.0
    _write_pae(
        root / "metabolism" / "targets" / "CX1.pae.json",
        matrix=matrix,
        chain_ids=["A", "A", "C", "C", "B", "B"],
    )
    row = _row(chain_subunits={"A": "SubA", "C": "SubA", "B": "SubB"})

    view = ui.complex_pae_view(row, root)

    assert view["status"] == ui.STATUS_OK
    assert view["mapping"] == {"SubA": ["A", "C"], "SubB": ["B"]}
    pair = view["rows"][0]
    assert pair["chains_a"] == ["A", "C"]
    assert pair["mean_a_to_b"] == pytest.approx(5.0)
    assert pair["max_a_to_b"] == pytest.approx(6.0)
    assert pair["mean_b_to_a"] == pytest.approx(9.0)
    assert pair["max_b_to_a"] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# Tabla
# ---------------------------------------------------------------------------
def test_table_rows_have_both_directions_and_arrows(tmp_path):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")
    view = ui.complex_pae_view(_row(**_ok_mapping()), root)

    rows = ui.table_rows(view)

    assert len(rows) == 1
    row = rows[0]
    assert row["Par"] == "SubA / SubB"
    assert row["Cadenas A"] == "A"
    assert row["Cadenas B"] == "B"
    assert row["A→B media (Å)"] == "5.0"
    assert row["A→B máx (Å)"] == "6.0"
    assert row["B→A media (Å)"] == "9.0"
    assert row["B→A máx (Å)"] == "10.0"
    assert row["Pares (residuos)"] == 12
    assert row["Estado"] == "ok"
    assert row["Motivo"] == ""


def test_table_rows_empty_when_no_analysis():
    view = ui.complex_pae_view(_row(), None)
    assert ui.table_rows(view) == []


def test_caveat_mentions_interaction_and_stoichiometry():
    text = ui.CAVEAT.casefold()
    assert "no demuestra interacción" in text or "no demuestra interaccion" in text
    assert "estequiometría" in text or "estequiometria" in text
    assert "no se infiere" in text


# ---------------------------------------------------------------------------
# Render con Streamlit falso
# ---------------------------------------------------------------------------
class _FakeStreamlit:
    """Registra llamadas; no necesita un runtime real de Streamlit."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple, dict]] = []

    def _record(self, name: str, *args, **kwargs) -> None:
        self.calls.append((name, args, kwargs))

    def markdown(self, *args, **kwargs) -> None:
        self._record("markdown", *args, **kwargs)

    def caption(self, *args, **kwargs) -> None:
        self._record("caption", *args, **kwargs)

    def info(self, *args, **kwargs) -> None:
        self._record("info", *args, **kwargs)

    def warning(self, *args, **kwargs) -> None:
        self._record("warning", *args, **kwargs)

    def dataframe(self, data, **kwargs) -> None:
        self._record("dataframe", data, **kwargs)

    def names(self) -> list[str]:
        return [name for name, _args, _kwargs in self.calls]

    def texts(self) -> list[str]:
        return [
            str(args[0]) for _name, args, _kwargs in self.calls if args
        ]

    def dataframes(self) -> list[tuple]:
        return [args for name, args, _kwargs in self.calls
                if name == "dataframe"]

    def dataframe_keys(self) -> list[str]:
        return [kwargs.get("key", "") for name, _args, kwargs in self.calls
                if name == "dataframe"]


@pytest.fixture
def fake_st(monkeypatch):
    fake = _FakeStreamlit()
    monkeypatch.setattr(ui, "st", fake)
    return fake


def test_render_ok_draws_table_and_caveat(tmp_path, fake_st):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")

    view = ui.render_complex_pae(
        _row(**_ok_mapping()), root, key_prefix="cx_"
    )

    assert view["status"] == ui.STATUS_OK
    assert "dataframe" in fake_st.names()
    assert fake_st.dataframe_keys()[0].startswith("cx_complex_pae_")
    table = fake_st.dataframes()[0][0]
    assert table[0]["A→B media (Å)"] == "5.0"
    assert table[0]["B→A máx (Å)"] == "10.0"
    joined = " ".join(fake_st.texts()).casefold()
    assert "no demuestra interacción" in joined
    assert "estequiometría" in joined
    assert "info" not in fake_st.names()
    assert "warning" not in fake_st.names()


def test_render_not_available_uses_status_message(tmp_path, fake_st):
    root = _root(tmp_path)

    view = ui.render_complex_pae(_row(), root)

    assert view["status"] == ui.STATUS_NOT_AVAILABLE
    infos = [str(args[0]) for name, args, _kwargs in fake_st.calls
             if name == "info" and args]
    assert infos
    assert "not_available" in infos[0]
    assert "dataframe" not in fake_st.names()
    assert "estequiometría" in " ".join(fake_st.texts()).casefold()


def test_render_invalid_uses_warning(tmp_path, fake_st):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "P1__P2.pae.json")
    _write_pae(root / "metabolism" / "targets" / "P1_P2.pae.json")
    row = _row(complex_id="OTRO", subunits="P1;P2")

    view = ui.render_complex_pae(row, root)

    assert view["status"] == ui.STATUS_INVALID
    warnings = [str(args[0]) for name, args, _kwargs in fake_st.calls
                if name == "warning" and args]
    assert warnings
    assert "invalid" in warnings[0]
    assert "dataframe" not in fake_st.names()


def test_render_key_prefix_namespaces_dataframe(tmp_path, fake_st):
    root = _root(tmp_path)
    _write_pae(root / "metabolism" / "targets" / "CX1.pae.json")

    ui.render_complex_pae(_row(**_ok_mapping()), root, key_prefix="a_")
    ui.render_complex_pae(_row(**_ok_mapping()), root, key_prefix="b_")

    keys = fake_st.dataframe_keys()
    assert len(keys) == 2
    assert keys[0].startswith("a_complex_pae_")
    assert keys[1].startswith("b_complex_pae_")
    assert keys[0] != keys[1]
