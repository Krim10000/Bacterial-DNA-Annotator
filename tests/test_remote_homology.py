"""Auditoría CPU de homología remota B/C.

Cubre el contrato mínimo: sin DB local no falla, fila A intacta, columnas del
TSV y ausencia de red incluso cuando se ejecuta una búsqueda local.
"""
from __future__ import annotations

import csv
import socket
from pathlib import Path

import pytest

from webapp import remote_homology as rh


def _rows() -> list[dict[str, str]]:
    return [
        {"protein_id": "pA", "path": "A", "label": "A_named_function"},
        {"protein_id": "pB", "path": "B",
         "label": "B_uncharacterized_homolog"},
        {"protein_id": "pC", "path": "C", "label": "C_no_swissprot_hit"},
    ]


def _proteins() -> list[dict[str, str]]:
    return [
        {"id": "pA", "aa": "M" * 120},
        {"id": "pB", "aa": "M" * 100},
        {"id": "pC", "aa": "M" * 80},
    ]


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


def _press(hmm: Path) -> None:
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        _touch(Path(str(hmm) + ext))


def _completed() -> object:
    return type("Completed", (), {"returncode": 0, "stdout": "",
                                  "stderr": ""})()


def test_no_local_db_marks_not_available_without_failing(tmp_path: Path):
    rows = _rows()
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path,
                    threads=2)

    assert result["ok"] is True
    assert result["n_bc"] == 2
    assert result["n_hit"] == 0
    assert result["n_status"] == {"not_available": 2}
    assert result["searched"] == []
    for row in rows[1:]:
        assert row["remote_status"] == "not_available"
        assert row["remote_source"] == ""
        assert row["remote_hit"] == ""
        assert row["remote_evalue"] == ""
        assert "sin DB local" in row["remote_note"]

    info = rh.available(tmp_path)
    assert set(info) == set(rh.SOURCES)
    assert all(not item["available"] for item in info.values())

    tsv = tmp_path / "out" / "g1.remote_homology.tsv"
    assert tsv.is_file()
    data, header = _read_tsv(tsv)
    assert header == list(rh.TSV_COLUMNS)
    assert [row["protein_id"] for row in data] == ["pB", "pC"]
    assert {row["remote_status"] for row in data} == {"not_available"}


def test_a_row_is_not_modified(tmp_path: Path):
    rows = _rows()
    before = dict(rows[0])
    rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert rows[0] == before
    assert not any(key.startswith("remote_") for key in rows[0])


def test_tsv_column_contract(tmp_path: Path):
    rows = _rows()
    rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)
    data, header = _read_tsv(tmp_path / "out" / "g1.remote_homology.tsv")

    assert header == ["genome", "protein_id", "path", "remote_status",
                      "remote_source", "remote_hit", "remote_evalue",
                      "remote_note"]
    assert header == list(rh.TSV_COLUMNS)
    assert len(data) == 2
    for row in data:
        assert set(row) == set(header)
        assert row["genome"] == "g1"
        assert row["path"] in {"B", "C"}
        assert row["remote_status"] == "not_available"


def test_available_detects_the_five_local_sources(tmp_path: Path,
                                                  monkeypatch):
    monkeypatch.setattr(
        rh.shutil, "which",
        lambda name: f"/opt/{name}" if name in ("diamond", "hmmscan") else None)
    db = tmp_path / "db"
    _touch(db / "uniref90.dmnd")
    _touch(db / "uniprot_trembl.dmnd")
    interpro = db / "interpro" / "interpro.hmm"
    _touch(interpro)
    _press(interpro)
    _touch(db / "eggnog" / "eggnog_proteins.dmnd")
    pfam = db / "pfam" / "Pfam-A.hmm"
    _touch(pfam)
    _press(pfam)

    info = rh.available(tmp_path)
    for name in rh.SOURCES:
        assert info[name]["available"] is True, name
        assert info[name]["searchable"] is True, name
    assert info["uniref"]["engine"] == "diamond"
    assert info["trembl"]["engine"] == "diamond"
    assert info["interpro"]["engine"] == "hmmscan"
    assert info["eggnog"]["engine"] == "diamond"
    assert info["pfam"]["engine"] == "hmmscan"
    assert info["uniref"]["role"] == "sequence"
    assert info["pfam"]["role"] == "profile"


def test_local_search_never_uses_network(tmp_path: Path, monkeypatch):
    _touch(tmp_path / "db" / "uniref90.dmnd")
    monkeypatch.setattr(
        rh.shutil, "which",
        lambda name: "/opt/diamond" if name == "diamond" else None)
    calls: list[list[str]] = []

    def fake_run(cmd, **kwargs):
        calls.append(list(cmd))
        out = Path(cmd[cmd.index("-o") + 1])
        out.write_text("", encoding="utf-8")
        return _completed()

    def deny(*args, **kwargs):
        raise AssertionError("la auditoría intentó usar la red")

    monkeypatch.setattr(rh.subprocess, "run", fake_run)
    monkeypatch.setattr(socket, "socket", deny)
    monkeypatch.setattr(socket, "create_connection", deny)

    rows = _rows()
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert result["ok"] is True
    assert result["searched"] == ["uniref"]
    assert result["n_status"] == {"no_hit": 2}
    assert calls and calls[0][0] == "/opt/diamond"
    for row in rows[1:]:
        assert row["remote_status"] == "no_hit"
        assert "sin hit en uniref" in row["remote_note"]


def test_hit_is_evidence_only_and_path_stays(tmp_path: Path, monkeypatch):
    _touch(tmp_path / "db" / "uniref90.dmnd")
    monkeypatch.setattr(
        rh.shutil, "which",
        lambda name: "/opt/diamond" if name == "diamond" else None)

    def fake_run(cmd, **kwargs):
        out = Path(cmd[cmd.index("-o") + 1])
        out.write_text(
            "pC\tUniRef90_P12345\t78.5\t80\t1e-30\t210\tremota\n",
            encoding="utf-8")
        return _completed()

    monkeypatch.setattr(rh.subprocess, "run", fake_run)
    rows = _rows()
    before = {row["protein_id"]: dict(row) for row in rows}
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert result["n_hit"] == 1
    c = rows[2]
    assert c["remote_status"] == "hit"
    assert c["remote_source"] == "uniref"
    assert c["remote_hit"] == "UniRef90_P12345"
    assert float(c["remote_evalue"]) == pytest.approx(1e-30)
    assert c["path"] == "C" and c["label"] == "C_no_swissprot_hit"
    assert "no función" in c["remote_note"]
    assert rows[1]["remote_status"] == "no_hit"
    assert rows[1]["path"] == "B"
    for row in rows[1:]:
        added = set(row) - set(before[row["protein_id"]])
        assert added == set(rh.REMOTE_FIELDS)
        assert row["path"] == before[row["protein_id"]]["path"]


def test_pfam_hmmscan_audits_remaining_rows(tmp_path: Path, monkeypatch):
    pfam = tmp_path / "db" / "pfam" / "Pfam-A.hmm"
    _touch(pfam)
    _press(pfam)
    monkeypatch.setattr(
        rh.shutil, "which",
        lambda name: "/opt/hmmscan" if name == "hmmscan" else None)
    calls: list[list[str]] = []

    def fake_run(cmd, **kwargs):
        calls.append(list(cmd))
        out = Path(cmd[cmd.index("--domtblout") + 1])
        out.write_text(
            "PF00001.1 PF00001 100 pC - 80 1e-25 120 0 1 1 1e-25 1e-25 "
            "120 0 1 60 1 60 1 60 0.95 dominio_de_prueba\n",
            encoding="utf-8")
        return _completed()

    monkeypatch.setattr(rh.subprocess, "run", fake_run)
    rows = _rows()
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert result["searched"] == ["pfam"]
    assert result["n_hit"] == 1
    c = rows[2]
    assert c["remote_status"] == "hit"
    assert c["remote_source"] == "pfam"
    assert c["remote_hit"] == "PF00001"
    assert float(c["remote_evalue"]) == pytest.approx(1e-25)
    assert rows[1]["remote_status"] == "no_hit"
    assert calls and "--domtblout" in calls[0]


def test_eggnog_without_diamond_is_not_searched(tmp_path: Path):
    _touch(tmp_path / "db" / "eggnog.db")
    rows = _rows()
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert result["ok"] is True
    assert result["n_status"] == {"not_searched": 2}
    for row in rows[1:]:
        assert row["remote_status"] == "not_searched"
        assert "eggnog" in row["remote_note"]


def test_no_bc_rows_writes_header_only(tmp_path: Path):
    rows = [{"protein_id": "pA", "path": "A"}]
    result = rh.run(rows, _proteins(), tmp_path / "out", "g1", tmp_path)

    assert result["ok"] is True
    assert result["n_bc"] == 0
    data, header = _read_tsv(tmp_path / "out" / "g1.remote_homology.tsv")
    assert data == []
    assert header == list(rh.TSV_COLUMNS)
