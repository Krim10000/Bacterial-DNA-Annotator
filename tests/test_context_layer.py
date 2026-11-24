"""Contexto del camino C: heurísticas declaradas y contrato de filas.

Secuencias sintéticas: hidrofóbica (TM), señal péptido Sec/SPI, lipoproteína
(lipobox) y soluble. No se usan herramientas externas ni red.
"""
from __future__ import annotations

import csv
from pathlib import Path

from webapp.context_layer import (
    CONTEXT_FIELDS,
    CONTEXT_NOTE,
    context_call,
    hydrophobic_segments,
    run,
    signal_peptide_call,
)

# Segmento hidrofóbico interno (no N-terminal): debe dar membrane_like.
HYDROPHOBIC = ("MKTAYIAKQRQISFVKSHFSRQL"
               "LLLVVILLLVVILLLVVILLL"
               "MKTAYIAKQRQISFVKSHFSRQL")
# Región h N-terminal + corte Sec/SPI (A-S-A) sin C posterior: secreted_like.
SIGNAL = ("MKKLLLLLLLLLLMLVL" + "ASA" + "QASAT"
          + "MKTAYIAKQRQISFVKSHFSRQL")
# Señal de lipoproteína (lipobox LAGC) dentro de los primeros 40 aa:
# ancla de membrana, no secreción soluble.
LIPOBOX = ("MKKLLLLLLLLLVL" + "LAGC" + "QASAT"
           + "MKTAYIAKQRQISFVKSHFSRQL")
SOLUBLE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN"
# Señal Sec realista sin ninguna ventana KD19 sobre umbral: no debe convertirse
# en membrane_like solo por tener señal péptido.
SHORT_SIGNAL = "MKQSTIALALLPLLFTPVTKA" + "RTPEMPVLENR"


def _write_check(tmp_path: Path, rows: list[dict], proteins: list[dict],
                 genome: str = "g") -> tuple[dict, Path]:
    res = run(rows, proteins, tmp_path, genome)
    return res, tmp_path / f"{genome}.context_c.tsv"


def test_hydrophobic_is_membrane_like_without_function(tmp_path: Path):
    rows = [{"protein_id": "p1", "path": "C"}]
    res, tsv = _write_check(tmp_path, rows, [{"id": "p1", "aa": HYDROPHOBIC}])
    assert res["ok"] is True
    assert res["n"] == 1
    row = rows[0]
    assert row["tm_segments"] >= 1
    assert row["context_hypothesis"] == "membrane_like"
    assert row["lipobox"] == "no"
    assert row["signal_peptide"] in ("no", "uncertain")
    assert row["context_note"] == CONTEXT_NOTE
    assert "no asigna función" in row["context_note"]
    assert "no equivale a TMHMM" in row["context_note"]
    assert tsv.exists()

    with tsv.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert reader.fieldnames == ["genome", "protein_id", *CONTEXT_FIELDS]
        data = list(reader)
    assert len(data) == 1
    assert data[0]["genome"] == "g"
    assert data[0]["protein_id"] == "p1"
    assert data[0]["context_hypothesis"] == "membrane_like"


def test_signal_peptide_is_secreted_like(tmp_path: Path):
    rows = [{"protein_id": "p2", "path": "C"}]
    res, tsv = _write_check(tmp_path, rows, [{"id": "p2", "aa": SIGNAL}])
    assert res["ok"] is True
    row = rows[0]
    assert row["signal_peptide"] == "yes"
    assert row["lipobox"] == "no"
    assert row["context_hypothesis"] == "secreted_like"
    assert "cleavage=si" in row["context_evidence"]
    assert tsv.exists()


def test_lipobox_is_membrane_anchored(tmp_path: Path):
    rows = [{"protein_id": "p3", "path": "C"}]
    res, _ = _write_check(tmp_path, rows, [{"id": "p3", "aa": LIPOBOX}])
    assert res["ok"] is True
    row = rows[0]
    assert row["lipobox"] == "yes"
    assert row["signal_peptide"] == "yes"
    assert row["context_hypothesis"] == "membrane_like"
    assert "lipobox" in row["context_evidence"]


def test_short_signal_without_tm_is_still_secreted_like(tmp_path: Path):
    rows = [{"protein_id": "p5", "path": "C"}]
    res, _ = _write_check(tmp_path, rows, [{"id": "p5", "aa": SHORT_SIGNAL}])
    assert res["ok"] is True
    row = rows[0]
    assert row["signal_peptide"] == "yes"
    assert row["tm_segments"] == 0
    assert row["context_hypothesis"] == "secreted_like"


def test_soluble_sequence(tmp_path: Path):
    rows = [{"protein_id": "p4", "path": "C"}]
    res, _ = _write_check(tmp_path, rows, [{"id": "p4", "aa": SOLUBLE}])
    assert res["ok"] is True
    row = rows[0]
    assert row["tm_segments"] == 0
    assert row["signal_peptide"] == "no"
    assert row["lipobox"] == "no"
    assert row["context_hypothesis"] == "soluble_like"


def test_only_c_rows_are_touched_and_path_never_changes(tmp_path: Path):
    rows = [
        {"protein_id": "a", "path": "A", "label": "A_named_function"},
        {"protein_id": "b", "path": "B", "label": "B_uncharacterized_homolog"},
        {"protein_id": "c", "path": "C", "label": "C_no_swissprot_hit"},
    ]
    proteins = [
        {"id": "a", "aa": HYDROPHOBIC},
        {"id": "b", "aa": SIGNAL},
        {"id": "c", "aa": SOLUBLE},
    ]
    res, tsv = _write_check(tmp_path, rows, proteins, genome="solo_c")
    assert res["ok"] is True
    assert res["n"] == 1
    assert [r["path"] for r in rows] == ["A", "B", "C"]
    assert rows[0]["label"] == "A_named_function"
    assert rows[1]["label"] == "B_uncharacterized_homolog"
    for key in CONTEXT_FIELDS:
        assert key not in rows[0]
        assert key not in rows[1]
        assert key in rows[2]
    with tsv.open(encoding="utf-8") as fh:
        data = list(csv.DictReader(fh, delimiter="\t"))
    assert [d["protein_id"] for d in data] == ["c"]


def test_no_c_rows_is_ok_and_writes_header(tmp_path: Path):
    rows = [{"protein_id": "a", "path": "A"}]
    res, tsv = _write_check(tmp_path, rows, [{"id": "a", "aa": SOLUBLE}])
    assert res["ok"] is True
    assert res["n"] == 0
    assert tsv.read_text(encoding="utf-8").strip().split("\t") == [
        "genome", "protein_id", *CONTEXT_FIELDS]


def test_c_row_without_sequence_is_unknown(tmp_path: Path):
    rows = [{"protein_id": "c", "path": "C"}]
    res, _ = _write_check(tmp_path, rows, [{"id": "c", "aa": ""}])
    assert res["ok"] is True
    assert rows[0]["context_hypothesis"] == "unknown"
    assert rows[0]["context_evidence"] == "sin_secuencia_para_la_proteina_C"


def test_unusable_outdir_returns_reason(tmp_path: Path):
    blocker = tmp_path / "file"
    blocker.write_text("x", encoding="utf-8")
    res = run([{"protein_id": "c", "path": "C"}],
              [{"id": "c", "aa": SOLUBLE}], blocker, "g")
    assert res["ok"] is False
    assert "no se pudo crear" in res["reason"]


def test_helper_segments_and_signal_states():
    assert hydrophobic_segments(SOLUBLE) == []
    segs = hydrophobic_segments(HYDROPHOBIC)
    assert len(segs) == 1
    start, end, mean_max = segs[0]
    assert start < end
    assert mean_max >= 1.6
    assert signal_peptide_call(SIGNAL)["status"] == "yes"
    assert signal_peptide_call(SOLUBLE)["status"] == "no"
    assert context_call(SOLUBLE)["context_hypothesis"] == "soluble_like"
    assert context_call("")["context_hypothesis"] == "unknown"
