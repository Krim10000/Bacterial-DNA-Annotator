"""Pruebas de la vista de hipótesis estructurales (``webapp/hypotheses_ui.py``).

Cubren la función pura ``safe_foldseek_detail`` (evidencia en dict, JSON,
texto plano tipo ``q9plu2`` y basura), la clasificación de confianza y su
explicación, la priorización de filas y un humo de ``render_hypotheses`` con un
Streamlit falso que falla si se le pasa un string a ``st.json``.

No hay datos reales ni PDB descargados: todo es sintético.
"""
from __future__ import annotations

import copy
import itertools
import json
from pathlib import Path

import pytest

import hypotheses_ui as ui
import structural_hypotheses


def _row(**overrides) -> dict:
    row = {
        "genome": "GCF_TEST",
        "protein_id": "P1",
        "uniprot": "Q9PLU2",
        "pdb": "",
        "hypothesis": "posible EC 3.1.1.1",
        "ec": "3.1.1.1",
        "cofactor": "",
        "pocket": "",
        "template_mcsa": "1abc",
        "foldseek_top": "q9plu2",
        "f1_foldseek": True,
        "f2_enzymm": True,
        "f3_pocket": False,
        "f4_alphafill": False,
        "f5_active_site": False,
        "n_sources": 2,
        "confidence": "orientative",
        "review_decision": "pending",
        "limitation": "Hipótesis computacional; no confirma función.",
        "evidence_json": "q9plu2",
    }
    row.update(overrides)
    return row


# ---------------------------------------------------------------------------
# safe_foldseek_detail
# ---------------------------------------------------------------------------
def test_detail_from_dict_with_row():
    row = _row(foldseek_top="", evidence_json={
        "row": {"foldseek_top": "1abc", "foldseek_bits": 120},
        "sources": ["Foldseek"],
    })
    out = ui.safe_foldseek_detail(row)
    assert out["kind"] == "structured"
    assert out["has_detail"] is True
    assert out["detail"] == {"foldseek_top": "1abc", "foldseek_bits": 120}
    assert out["top"] == "1abc"
    assert out["text"] == ""


def test_detail_from_json_string():
    evidence = json.dumps({"row": {"foldseek_bits": 88,
                                   "active_positions": "10;20"}})
    out = ui.safe_foldseek_detail(_row(evidence_json=evidence))
    assert out["kind"] == "structured"
    assert out["detail"]["active_positions"] == "10;20"
    assert out["top"] == "q9plu2"


def test_top_falls_back_to_evidence_row():
    out = ui.safe_foldseek_detail(_row(
        foldseek_top="",
        evidence_json={"row": {"foldseek_top": "9zzz"}},
    ))
    assert out["top"] == "9zzz"
    assert out["detail"] == {"foldseek_top": "9zzz"}


@pytest.mark.parametrize("plain", ["q9plu2", "  q9plu2  ", "no hay json"])
def test_plain_text_is_shown_as_text_never_as_json(plain):
    out = ui.safe_foldseek_detail(_row(evidence_json=plain, foldseek_top=""))
    assert isinstance(out["detail"], dict)
    assert out["detail"] == {}
    assert out["has_detail"] is False
    assert out["kind"] == "text"
    assert out["text"] == plain.strip()
    assert "texto plano" in out["message"]


def test_json_encoded_string_is_text_not_object():
    out = ui.safe_foldseek_detail(_row(evidence_json='"q9plu2"'))
    assert out["detail"] == {}
    assert out["kind"] == "text"
    assert out["text"] == "q9plu2"


def test_json_list_is_scalar_text():
    out = ui.safe_foldseek_detail(_row(evidence_json=json.dumps(["a", "b"])))
    assert out["detail"] == {}
    assert out["kind"] == "scalar"
    assert "a" in out["text"]


def test_missing_evidence_keeps_top():
    out = ui.safe_foldseek_detail(_row(evidence_json=None,
                                       foldseek_top="7xyz"))
    assert out["kind"] == "top_only"
    assert out["top"] == "7xyz"
    assert out["detail"] == {}
    assert out["has_detail"] is False


def test_missing_everything_is_empty():
    out = ui.safe_foldseek_detail(_row(evidence_json="", foldseek_top=""))
    assert out["kind"] == "empty"
    assert out["detail"] == {}
    assert out["has_detail"] is False
    assert out["message"]


def test_evidence_dict_without_row_is_kept():
    out = ui.safe_foldseek_detail(_row(evidence_json={
        "sources": ["Foldseek"], "meta": {"genome": "G"},
    }))
    assert out["detail"] == {"sources": ["Foldseek"],
                             "meta": {"genome": "G"}}
    assert out["kind"] == "structured"


def test_evidence_row_as_plain_string_is_text():
    out = ui.safe_foldseek_detail(_row(evidence_json={"row": "q9plu2"}))
    assert out["detail"] == {}
    assert out["kind"] == "text"
    assert out["text"] == "q9plu2"


@pytest.mark.parametrize("bad_row", [None, "cadena", 42, 3.5, ["a"], ["x", 1]])
def test_invalid_row_never_raises(bad_row):
    out = ui.safe_foldseek_detail(bad_row)
    assert out["kind"] == "invalid"
    assert out["detail"] == {}
    assert out["has_detail"] is False
    assert out["message"]


@pytest.mark.parametrize("weird", [b"bytes", 12, 1.5, True, {"row": ""}])
def test_weird_evidence_json_never_raises(weird):
    out = ui.safe_foldseek_detail(_row(evidence_json=weird))
    assert isinstance(out["detail"], dict)
    assert isinstance(out["text"], str)
    assert out["kind"] in {"scalar", "text", "empty", "top_only", "structured"}


def test_safe_detail_is_pure_and_returns_fresh_dicts():
    row = _row(evidence_json={"row": {"a": 1}})
    before = copy.deepcopy(row)
    first = ui.safe_foldseek_detail(row)
    second = ui.safe_foldseek_detail(row)
    assert row == before
    assert first == second
    assert first["detail"] is not second["detail"]


# ---------------------------------------------------------------------------
# Clasificación y explicación de confianza
# ---------------------------------------------------------------------------
def _flags(**on) -> dict:
    row = {field: False for field in ui._FLAG_FIELDS}
    row.update(on)
    return row


@pytest.mark.parametrize("flags, expected", [
    ({"f2_enzymm": True, "f3_pocket": True, "f4_alphafill": True},
     "supported_hypothesis"),
    ({"f2_enzymm": True, "f5_active_site": True}, "high_priority"),
    ({"f2_enzymm": True}, "orientative"),
    ({"f1_foldseek": True, "f3_pocket": True}, "structural_only"),
    ({"f1_foldseek": True, "f4_alphafill": True}, "structural_only"),
    ({"f1_foldseek": True}, "uncharacterized"),
    ({}, "uncharacterized"),
])
def test_classify_confidence_matches_documented_rule(flags, expected):
    assert ui.classify_confidence(_flags(**flags)) == expected


def test_classification_matches_pipeline_build_for_all_flag_combinations():
    """La GUI no debe divergir de ``structural_hypotheses.build``."""
    for combo in itertools.product([False, True], repeat=5):
        row = dict(zip(ui._FLAG_FIELDS, combo))
        built = structural_hypotheses.build(dict(row), {"protein_id": "P"})
        assert ui.classify_confidence(row) == built["confidence"], row


def test_conflict_has_priority_over_support():
    row = _flags(f2_enzymm=True, f3_pocket=True, f4_alphafill=True)
    row["confidence"] = "conflict"
    assert ui.classify_confidence(row) == "conflict"
    row = _flags(f2_enzymm=True, f3_pocket=True, f4_alphafill=True)
    row["review_decision"] = "conflict"
    assert ui.classify_confidence(row) == "conflict"
    row = _flags(f2_enzymm=True, f3_pocket=True, f4_alphafill=True)
    row["conflict"] = "true"
    assert ui.classify_confidence(row) == "conflict"


def test_conflict_string_false_is_not_a_conflict():
    row = _flags(f2_enzymm=True)
    row["conflict"] = "False"
    assert ui.classify_confidence(row) == "orientative"


@pytest.mark.parametrize("value", ["True", "true", "1", "sí", 1, True])
def test_truthy_flags_from_pipeline_and_tsv(value):
    row = _flags(f2_enzymm=value)
    assert ui.classify_confidence(row) == "orientative"


@pytest.mark.parametrize("value", ["False", "false", "0", "", 0, False, None])
def test_falsey_flags_from_pipeline_and_tsv(value):
    assert ui.classify_confidence(_flags(f2_enzymm=value)) == "uncharacterized"


def test_classify_confidence_invalid_row():
    assert ui.classify_confidence(None) == "uncharacterized"
    assert ui.classify_confidence("cadena") == "uncharacterized"


@pytest.mark.parametrize("stored", [
    "supported_hypothesis", "high_priority", "conflict", "orientative",
    "structural_only", "uncharacterized",
])
def test_confidence_label_keeps_known_label(stored):
    assert ui.confidence_label({"confidence": stored}) == stored


def test_confidence_label_computes_when_missing_or_unknown():
    assert ui.confidence_label({"f2_enzymm": True}) == "orientative"
    assert ui.confidence_label(
        {"confidence": "raro", "f2_enzymm": True}
    ) == "orientative"
    assert ui.confidence_label(
        {"confidence": "  SUPPORTED_HYPOTHESIS  "}
    ) == "supported_hypothesis"


def test_confidence_label_invalid_row():
    assert ui.confidence_label(None) == "uncharacterized"
    assert ui.confidence_label(42) == "uncharacterized"


def test_confidence_info_explains_every_label():
    for key in ("supported_hypothesis", "high_priority", "conflict",
                "orientative", "structural_only", "uncharacterized"):
        info = ui.confidence_info(key)
        assert info["key"] == key
        assert info["emoji"]
        assert info["label"]
        assert info["rule"]
        assert info["meaning"]
        assert info["caveat"]


def test_confidence_info_normalizes_case_and_spaces():
    info = ui.confidence_info("  SUPPORTED_HYPOTHESIS ")
    assert info["key"] == "supported_hypothesis"
    assert info["raw"] == "SUPPORTED_HYPOTHESIS"


def test_confidence_info_unknown_falls_back_without_inventing():
    info = ui.confidence_info("otra_cosa")
    assert info["key"] == "unknown"
    assert info["raw"] == "otra_cosa"
    assert info["meaning"]
    assert ui.confidence_info(None)["key"] == "unknown"


def test_conflict_explanation_is_honest():
    info = ui.confidence_info("conflict")
    text = (info["meaning"] + " " + info["caveat"]).lower()
    assert "no" in text
    assert "brecha" in text
    assert "anotación" in text


def test_supported_explanation_does_not_claim_confirmation():
    info = ui.confidence_info("supported_hypothesis")
    text = (info["meaning"] + " " + info["caveat"]).lower()
    assert "no confirma" in text or "no son un ensayo" in text


def test_confidence_catalog_is_complete():
    required = {"emoji", "label", "rule", "meaning", "caveat"}
    for key, info in ui.CONFIDENCE_INFO.items():
        assert required <= set(info), key
        assert info["meaning"].strip()
        assert info["rule"].strip()
        assert info["caveat"].strip()


def test_fields_guide_explains_f1_to_f5():
    guide = {item["Filtro"]: item for item in ui.FIELDS_GUIDE}
    assert set(guide) == {"F1", "F2", "F3", "F4", "F5"}
    for item in guide.values():
        assert item["Fuente"]
        assert item["Qué mide"]
        assert item["Qué no prueba"]


# ---------------------------------------------------------------------------
# Priorización
# ---------------------------------------------------------------------------
def test_prioritize_rows_orders_by_confidence_then_sources():
    rows = [
        {"protein_id": "low", "confidence": "uncharacterized", "n_sources": 0},
        {"protein_id": "sup", "confidence": "supported_hypothesis",
         "n_sources": 3},
        {"protein_id": "ori2", "confidence": "orientative", "n_sources": 2},
        {"protein_id": "ori1", "confidence": "orientative", "n_sources": 1},
        {"protein_id": "conf", "confidence": "conflict", "n_sources": 1},
    ]
    order = [row["protein_id"] for row in ui.prioritize_rows(rows)]
    assert order == ["sup", "conf", "ori2", "ori1", "low"]


def test_prioritize_rows_compares_sources_numerically():
    rows = [
        {"protein_id": "a", "confidence": "orientative", "n_sources": "3"},
        {"protein_id": "b", "confidence": "orientative", "n_sources": "12"},
        {"protein_id": "c", "confidence": "orientative", "n_sources": ""},
    ]
    assert [row["protein_id"] for row in ui.prioritize_rows(rows)] == [
        "b", "a", "c",
    ]


def test_prioritize_rows_ignores_non_dicts_and_does_not_mutate():
    rows = ["basura", None, {"protein_id": "x", "confidence": "orientative"}]
    before = copy.deepcopy(rows)
    out = ui.prioritize_rows(rows)
    assert rows == before
    assert [row["protein_id"] for row in out] == ["x"]


# ---------------------------------------------------------------------------
# Humo de render_hypotheses con Streamlit falso
# ---------------------------------------------------------------------------
class _FakeColumn:
    def __init__(self, parent: "_FakeStreamlit") -> None:
        self._parent = parent

    def metric(self, label, value, **kwargs):
        self._parent.metric(label, value, **kwargs)

    def download_button(self, *args, **kwargs):
        self._parent.download_button(*args, **kwargs)


class _FakeExpander:
    def __init__(self, parent: "_FakeStreamlit") -> None:
        self._parent = parent

    def __enter__(self):
        return self._parent

    def __exit__(self, *exc_info):
        return False


class _FakeStreamlit:
    """Registra llamadas y reproduce los errores reales que importan."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple, dict]] = []
        self.session_state: dict = {}

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))

    def caption(self, *args, **kwargs):
        self._record("caption", *args, **kwargs)

    def markdown(self, *args, **kwargs):
        self._record("markdown", *args, **kwargs)

    def info(self, *args, **kwargs):
        self._record("info", *args, **kwargs)

    def warning(self, *args, **kwargs):
        self._record("warning", *args, **kwargs)

    def success(self, *args, **kwargs):
        self._record("success", *args, **kwargs)

    def divider(self, *args, **kwargs):
        self._record("divider", *args, **kwargs)

    def write(self, *args, **kwargs):
        self._record("write", *args, **kwargs)

    def text(self, body, **kwargs):
        if not isinstance(body, str):
            raise TypeError(f"st.text recibió {type(body).__name__}")
        self._record("text", body, **kwargs)

    def json(self, data, **kwargs):
        # Regresión del bug real: st.json("q9plu2") no debe ocurrir.
        if not isinstance(data, dict):
            raise TypeError(f"st.json recibió {type(data).__name__}")
        self._record("json", data, **kwargs)

    def dataframe(self, data, **kwargs):
        self._record("dataframe", data, **kwargs)

    def download_button(self, *args, **kwargs):
        self._record("download_button", *args, **kwargs)
        return True

    def metric(self, label, value, **kwargs):
        self._record("metric", label, value, **kwargs)

    def columns(self, n, **kwargs):
        self._record("columns", n, **kwargs)
        return [_FakeColumn(self) for _ in range(n)]

    def expander(self, label, **kwargs):
        self._record("expander", label, **kwargs)
        return _FakeExpander(self)

    def selectbox(self, label, options, index=0, **kwargs):
        self._record("selectbox", label, options, index=index, **kwargs)
        opts = list(options)
        return opts[index] if 0 <= index < len(opts) else opts[0]

    def downloaded_files(self) -> list[str]:
        return [kwargs.get("file_name", "")
                for name, _args, kwargs in self.calls
                if name == "download_button"]

    def texts(self) -> list[str]:
        return [args[0] for name, args, _kwargs in self.calls if name == "text"]

    def json_payloads(self) -> list[dict]:
        return [args[0] for name, args, _kwargs in self.calls if name == "json"]

    def metric_values(self) -> dict[str, object]:
        return {args[0]: args[1] for name, args, _kwargs in self.calls
                if name == "metric"}


class _FakeViewer:
    def __init__(self) -> None:
        self.calls: list[dict] = []

    def render(self, pdb, key="", caption=""):
        self.calls.append({"pdb": str(pdb), "key": key, "caption": caption})


@pytest.fixture
def fake_ui(monkeypatch):
    fake = _FakeStreamlit()
    viewer = _FakeViewer()
    monkeypatch.setattr(ui, "st", fake)
    monkeypatch.setattr(ui, "structure_viewer", viewer)
    return fake, viewer


def _pdb_file(tmp_path: Path) -> Path:
    path = tmp_path / "candidate.pdb"
    path.write_text(
        "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504  1.00  0.00\n"
        "HETATM    2  ZN  ZN  A 101       1.000   2.000   3.000  1.00  0.00\n"
        "END\n",
        encoding="utf-8",
    )
    return path


def test_render_plain_text_evidence_does_not_pass_string_to_json(fake_ui):
    fake, viewer = fake_ui
    ui.render_hypotheses({
        "genome": "G1",
        "structural_hypotheses": [_row(evidence_json="q9plu2", pdb="")],
    })
    assert "q9plu2" in fake.texts()
    assert fake.json_payloads() == []
    assert viewer.calls == []
    assert f"G1_structural_hypotheses.tsv" in fake.downloaded_files()


def test_render_structured_evidence_uses_json_dict(fake_ui):
    fake, _viewer = fake_ui
    evidence = json.dumps({"row": {"foldseek_top": "1abc",
                                   "foldseek_bits": 99}})
    ui.render_hypotheses({
        "genome": "G2",
        "structural_hypotheses": [_row(evidence_json=evidence)],
    })
    payloads = fake.json_payloads()
    assert payloads
    assert payloads[0] == {"foldseek_top": "1abc", "foldseek_bits": 99}


def test_render_missing_evidence_shows_no_detail(fake_ui):
    fake, _viewer = fake_ui
    ui.render_hypotheses({
        "genome": "G3",
        "structural_hypotheses": [_row(evidence_json=None, foldseek_top="")],
    })
    captions = [args[0] for name, args, _kwargs in fake.calls
                if name == "caption"]
    assert any("Sin detalle de Foldseek" in str(c) for c in captions)
    assert fake.json_payloads() == []


def test_render_with_mixed_rows_metrics_and_pdb(fake_ui, tmp_path):
    fake, viewer = fake_ui
    pdb = _pdb_file(tmp_path)
    rows = [
        _row(protein_id="P1", confidence="supported_hypothesis",
             n_sources=3, pdb=str(pdb)),
        _row(protein_id="P2", confidence="conflict", n_sources=1,
             evidence_json="q9plu2"),
        _row(protein_id="P3", confidence=None, n_sources="2",
             f1_foldseek=False, f2_enzymm=False, f3_pocket=False,
             f4_alphafill=False, f5_active_site=False),
    ]
    ui.render_hypotheses({"genome": "G4", "structural_hypotheses": rows})

    metrics = fake.metric_values()
    assert metrics["Hipótesis"] == 3
    assert metrics["Con ≥2 fuentes"] == 2
    assert metrics["Apoyadas"] == 1
    assert metrics["Conflictos"] == 1
    assert metrics["Alta prioridad"] == 0
    assert viewer.calls and viewer.calls[0]["pdb"] == str(pdb)
    assert "G4_structural_hypotheses.tsv" in fake.downloaded_files()
    assert pdb.name in fake.downloaded_files()
    assert fake.json_payloads() == []


def test_render_weird_rows_never_raise(fake_ui):
    fake, viewer = fake_ui
    ui.render_hypotheses({
        "genome": "G5",
        "structural_hypotheses": [
            "basura",
            {"protein_id": "P1", "evidence_json": ["no", "es", "dict"],
             "n_sources": "no-numero"},
            {"protein_id": "P2", "evidence_json": {"row": "q9plu2"}},
        ],
    })
    assert viewer.calls == []
    assert any(name == "download_button" for name, _args, _kwargs in fake.calls)


def test_render_without_rows_returns_info(fake_ui):
    fake, viewer = fake_ui
    ui.render_hypotheses({"genome": "G6", "structural_hypotheses": []})
    infos = [args[0] for name, args, _kwargs in fake.calls if name == "info"]
    assert any("Sin hipótesis estructurales" in str(x) for x in infos)
    assert viewer.calls == []


def test_render_non_list_hypotheses_does_not_raise(fake_ui):
    fake, viewer = fake_ui
    ui.render_hypotheses({"genome": "G7", "structural_hypotheses": "texto"})
    assert viewer.calls == []
    assert any(name == "info" for name, _args, _kwargs in fake.calls)
