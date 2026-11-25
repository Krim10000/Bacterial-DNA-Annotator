"""Estructuras comprimidas y resolución por catálogo en la ficha del gen.

Las pruebas usan PDB/mmCIF y catálogos sintéticos en ``tmp_path``. No leen ni
escriben ``data/metabolism/targets`` y no levantan Streamlit: reemplazan la
interfaz y el visor por grabadores mínimos.
"""
from __future__ import annotations

import csv
import gzip
from contextlib import nullcontext
from pathlib import Path

import pytest

from webapp import gene_card, structure_viewer

PDB_TEXT = (
    "HEADER    SYNTHETIC TEST\n"
    "ATOM      1  CA  MET A   1      11.104  13.207   2.000  1.00 20.00           C\n"
    "ATOM      2  CA  ALA A   2      12.104  14.207   3.000  1.00 30.00           C\n"
    "END\n"
)
CIF_TEXT = (
    "data_TEST\n"
    "#\n"
    "loop_\n"
    "_atom_site.group_PDB\n"
    "_atom_site.id\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "ATOM 1 11.104 13.207 2.000\n"
    "#\n"
)


def _gzip(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode("utf-8")))
    return path


def _write_catalog(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t",
            fieldnames=["accession", "source", "path", "status"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# Lectura transparente en el visor
# ---------------------------------------------------------------------------
def test_read_pdb_keeps_plain_file_and_raw_text(tmp_path: Path):
    plain = tmp_path / "plain.pdb"
    plain.write_text(PDB_TEXT, encoding="utf-8")

    assert structure_viewer._read_pdb(plain) == PDB_TEXT
    # Si no es un archivo, se interpreta como texto de coordenadas.
    assert structure_viewer._read_pdb(PDB_TEXT) == PDB_TEXT


def test_read_pdb_decompresses_pdb_and_cif_gzip(tmp_path: Path):
    pdb_gz = _gzip(tmp_path / "model.pdb.gz", PDB_TEXT)
    cif_gz = _gzip(tmp_path / "model.cif.gz", CIF_TEXT)

    assert structure_viewer._read_pdb(pdb_gz) == PDB_TEXT
    assert structure_viewer._read_pdb(cif_gz) == CIF_TEXT
    assert structure_viewer._structure_format(pdb_gz, PDB_TEXT) == "pdb"
    assert structure_viewer._structure_format(cif_gz, CIF_TEXT) == "cif"
    # También se detecta mmCIF por contenido, sin depender del sufijo.
    assert structure_viewer._structure_format("model", CIF_TEXT) == "cif"


class _FakeView:
    def __init__(self) -> None:
        self.models: list[tuple[str, str]] = []

    def addModel(self, text, fmt) -> None:
        self.models.append((str(text), str(fmt)))

    def setStyle(self, *args, **kwargs) -> None:
        return None

    def addStyle(self, *args, **kwargs) -> None:
        return None

    def setBackgroundColor(self, *args, **kwargs) -> None:
        return None

    def zoomTo(self, *args, **kwargs) -> None:
        return None


class _FakePy3Dmol:
    def __init__(self) -> None:
        self.view_obj = _FakeView()

    def view(self, **kwargs) -> _FakeView:
        return self.view_obj


class _FakeStmol:
    def __init__(self) -> None:
        self.shown: list[_FakeView] = []

    def showmol(self, view, **kwargs) -> None:
        self.shown.append(view)


class _FakeSt:
    def __init__(self) -> None:
        self.infos: list[str] = []
        self.captions: list[str] = []
        self.codes: list[str] = []
        self.downloads: list[dict] = []

    def info(self, body, **kwargs) -> None:
        self.infos.append(str(body))

    def caption(self, body, **kwargs) -> None:
        self.captions.append(str(body))

    def code(self, body, **kwargs) -> None:
        self.codes.append(str(body))

    def download_button(self, label, data=None, file_name=None, mime=None,
                        key=None, **kwargs) -> None:
        self.downloads.append({
            "label": str(label), "data": data, "file_name": str(file_name),
            "mime": mime, "key": key,
        })


def test_render_adds_decompressed_gzip_model_to_py3dmol(
        monkeypatch, tmp_path: Path):
    pdb_gz = _gzip(tmp_path / "model.pdb.gz", PDB_TEXT)
    fake_st = _FakeSt()
    fake_py = _FakePy3Dmol()
    fake_stmol = _FakeStmol()
    monkeypatch.setattr(structure_viewer, "st", fake_st)
    monkeypatch.setattr(structure_viewer, "py3Dmol", fake_py)
    monkeypatch.setattr(structure_viewer, "stmol", fake_stmol)
    monkeypatch.setattr(structure_viewer, "available", lambda: (True, ""))

    structure_viewer.render(pdb_gz)

    assert fake_py.view_obj.models == [(PDB_TEXT, "pdb")]
    assert fake_st.downloads[0]["data"] == PDB_TEXT
    assert fake_st.downloads[0]["file_name"] == "pdb_viewer.pdb"
    assert fake_st.infos == []


def test_render_detects_cif_gzip_for_py3dmol(monkeypatch, tmp_path: Path):
    cif_gz = _gzip(tmp_path / "model.cif.gz", CIF_TEXT)
    fake_st = _FakeSt()
    fake_py = _FakePy3Dmol()
    monkeypatch.setattr(structure_viewer, "st", fake_st)
    monkeypatch.setattr(structure_viewer, "py3Dmol", fake_py)
    monkeypatch.setattr(structure_viewer, "stmol", _FakeStmol())
    monkeypatch.setattr(structure_viewer, "available", lambda: (True, ""))

    structure_viewer.render(cif_gz)

    assert fake_py.view_obj.models == [(CIF_TEXT, "cif")]
    assert fake_st.downloads[0]["file_name"] == "pdb_viewer.cif"
    assert fake_st.downloads[0]["mime"] == "chemical/x-cif"


# ---------------------------------------------------------------------------
# Ficha del gen: grabador mínimo
# ---------------------------------------------------------------------------
class _CardSt(_FakeSt):
    def __init__(self) -> None:
        super().__init__()
        self.markdowns: list[str] = []
        self.buttons: list[str] = []
        self.click = False

    def markdown(self, body, **kwargs) -> None:
        self.markdowns.append(str(body))

    def button(self, label, key=None, **kwargs) -> bool:
        self.buttons.append(str(label))
        return self.click

    def expander(self, label, **kwargs):
        return nullcontext()


def _gene(**overrides) -> dict:
    gene = {
        "genome": "demo",
        "contig": "c1",
        "start": 1,
        "end": 9,
        "strand": "+",
        "gene_id": "g1",
        "locus_tag": "TAG_1",
        "gene_name": "dnaA",
        "product": "Synthetic test protein",
        "status": "consensus",
        "confidence": "medium",
        "uniprot": "TESTACC",
    }
    gene.update(overrides)
    return gene


@pytest.fixture
def run_card(monkeypatch, tmp_path: Path):
    def _run(gene: dict, data_root: Path, click: bool = True):
        fake = _CardSt()
        fake.click = click
        viewer_calls: list[tuple[Path, dict]] = []
        monkeypatch.setattr(gene_card, "st", fake)
        monkeypatch.setattr(gene_card, "_language", lambda: "es")
        monkeypatch.setattr(
            gene_card.structure_viewer, "render",
            lambda pdb, **kwargs: viewer_calls.append((Path(pdb), kwargs)),
        )
        gene_card.render_gene_card(
            gene, "demo", data_root, sequences={"c1": "ATGAAATAA"},
            key_prefix="t",
        )
        return fake, viewer_calls

    return _run


def _pdb_download(fake: _CardSt) -> dict:
    matches = [d for d in fake.downloads if str(d["file_name"]).endswith(
        (".pdb", ".cif", ".mmcif"))]
    assert len(matches) == 1
    return matches[0]


# ---------------------------------------------------------------------------
# Orden de resolución: pdb declarado > catálogo > AFDB/ESMFold
# ---------------------------------------------------------------------------
def test_declared_plain_pdb_is_unchanged(run_card, tmp_path: Path):
    pdb = tmp_path / "declared.pdb"
    pdb.write_text(PDB_TEXT, encoding="utf-8")

    fake, viewer = run_card(_gene(pdb=str(pdb)), tmp_path)

    assert viewer and viewer[0][0] == pdb
    download = _pdb_download(fake)
    assert download["data"] == PDB_TEXT.encode("utf-8")
    assert download["file_name"] == "declared.pdb"
    assert download["mime"] == "chemical/x-pdb"


def test_declared_pdb_gzip_precedes_catalog(run_card, tmp_path: Path):
    declared = _gzip(tmp_path / "declared.pdb.gz", PDB_TEXT)
    targets = tmp_path / "data" / "metabolism" / "targets"
    catalog_model = targets / "models" / "TESTACC.cif.gz"
    _gzip(catalog_model, CIF_TEXT)
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "AFDB",
          "path": "data/metabolism/targets/models/TESTACC.cif.gz",
          "status": "ok"}],
    )

    fake, viewer = run_card(
        _gene(pdb=str(declared), uniprot="sp|TESTACC|TEST"), tmp_path / "data",
    )

    assert viewer and viewer[0][0] == declared
    download = _pdb_download(fake)
    assert download["data"] == PDB_TEXT.encode("utf-8")
    # La descarga se ofrece descomprimida y con el nombre sin .gz.
    assert download["file_name"] == "declared.pdb"
    assert download["mime"] == "chemical/x-pdb"


def test_catalog_resolves_accession_to_cif_gzip(run_card, tmp_path: Path):
    data_root = tmp_path / "data"
    targets = data_root / "metabolism" / "targets"
    model = _gzip(targets / "models" / "TESTACC.cif.gz", CIF_TEXT)
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "AFDB",
          "path": "data/metabolism/targets/models/TESTACC.cif.gz",
          "status": "ok"}],
    )

    fake, viewer = run_card(
        _gene(uniprot="sp|TESTACC|TEST"), data_root,
    )

    assert viewer and viewer[0][0] == model
    download = _pdb_download(fake)
    assert download["data"] == CIF_TEXT.encode("utf-8")
    assert download["file_name"] == "TESTACC.cif"
    assert download["mime"] == "chemical/x-cif"


def test_catalog_path_relative_to_catalog_directory(run_card, tmp_path: Path):
    data_root = tmp_path / "data"
    targets = data_root / "metabolism" / "targets"
    model = targets / "TESTACC.pdb"
    model.parent.mkdir(parents=True, exist_ok=True)
    model.write_text(PDB_TEXT, encoding="utf-8")
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "ESMFold",
          "path": "TESTACC.pdb", "status": "ok"}],
    )

    _fake, viewer = run_card(
        _gene(uniprot="TESTACC"), data_root,
    )

    assert viewer and viewer[0][0] == model


def test_catalog_precedes_afdb_fallback(run_card, tmp_path: Path):
    data_root = tmp_path / "data"
    targets = data_root / "metabolism" / "targets"
    catalog_model = _gzip(targets / "models" / "TESTACC.cif.gz", CIF_TEXT)
    (targets / "afdb25k").mkdir(parents=True, exist_ok=True)
    (targets / "afdb25k" / "TESTACC.pdb").write_text(
        PDB_TEXT, encoding="utf-8")
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "AFDB",
          "path": "models/TESTACC.cif.gz", "status": "ok"}],
    )

    _fake, viewer = run_card(_gene(uniprot="TESTACC"), data_root)

    assert viewer and viewer[0][0] == catalog_model


def test_missing_catalog_file_falls_back_to_compressed_afdb(run_card,
                                                            tmp_path: Path):
    data_root = tmp_path / "data"
    targets = data_root / "metabolism" / "targets"
    fallback = _gzip(targets / "afdb25k" / "TESTACC.pdb.gz", PDB_TEXT)
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "AFDB",
          "path": "models/does_not_exist.pdb", "status": "missing"}],
    )

    fake, viewer = run_card(_gene(uniprot="TESTACC"), data_root)

    assert viewer and viewer[0][0] == fallback
    assert _pdb_download(fake)["data"] == PDB_TEXT.encode("utf-8")


def test_missing_declared_pdb_falls_through_to_catalog(run_card,
                                                        tmp_path: Path):
    data_root = tmp_path / "data"
    targets = data_root / "metabolism" / "targets"
    model = targets / "TESTACC.pdb"
    model.parent.mkdir(parents=True, exist_ok=True)
    model.write_text(PDB_TEXT, encoding="utf-8")
    _write_catalog(
        targets / "structure_catalog.tsv",
        [{"accession": "TESTACC", "source": "ESMFold",
          "path": "TESTACC.pdb", "status": "ok"}],
    )

    _fake, viewer = run_card(
        _gene(pdb="missing.pdb", uniprot="TESTACC"), data_root,
    )

    assert viewer and viewer[0][0] == model


def test_no_structure_available_reports_it(run_card, tmp_path: Path):
    fake, viewer = run_card(
        _gene(uniprot="NO_CATALOG_ENTRY"), tmp_path / "empty",
    )

    assert viewer == []
    assert fake.buttons == []
    assert not [
        d for d in fake.downloads
        if str(d["file_name"]).endswith((".pdb", ".cif", ".mmcif"))
    ]
    assert any("No hay PDB local" in note for note in fake.infos)
