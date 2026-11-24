"""Ficha de gen: extracción, traducción y descargas condicionadas.

Las pruebas no levantan Streamlit: reemplazan ``gene_card.st`` por un grabador
mínimo y comprueban qué se muestra y qué descargas se ofrecen. Así se verifica
que la ficha no invente secuencia, traducción ni estructura cuando faltan.
"""
from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path

import pytest

from webapp import gene_card


class FakeSt:
    """Grabador de llamadas de Streamlit suficiente para la ficha."""

    def __init__(self) -> None:
        self.markdowns: list[str] = []
        self.captions: list[str] = []
        self.infos: list[str] = []
        self.warnings: list[str] = []
        self.codes: list[str] = []
        self.downloads: list[dict] = []
        self.buttons: list[str] = []
        self.click = False

    def markdown(self, body, **kwargs) -> None:
        self.markdowns.append(str(body))

    def caption(self, body, **kwargs) -> None:
        self.captions.append(str(body))

    def info(self, body, **kwargs) -> None:
        self.infos.append(str(body))

    def warning(self, body, **kwargs) -> None:
        self.warnings.append(str(body))

    def code(self, body, **kwargs) -> None:
        self.codes.append(str(body))

    def download_button(self, label, data=None, file_name=None, mime=None,
                        key=None, **kwargs) -> None:
        self.downloads.append({
            "label": str(label), "data": data, "file_name": str(file_name),
            "mime": mime, "key": key,
        })

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
        "product": "Chromosomal replication initiator protein DnaA",
        "gene_type": "CDS",
        "source": "genbank",
        "status": "consensus",
        "confidence": "medium",
    }
    gene.update(overrides)
    return gene


def _all_text(fake: FakeSt) -> str:
    return "\n".join(
        fake.markdowns + fake.captions + fake.infos + fake.warnings + fake.codes
    )


@pytest.fixture
def run_card(monkeypatch, tmp_path: Path):
    """Ejecuta la ficha con un Streamlit falso y devuelve lo registrado."""

    def _run(gene: dict, sequences=None, click: bool = False,
             prefix: str = "t", data_root: Path | None = None):
        fake = FakeSt()
        fake.click = click
        viewer_calls: list[tuple[str, dict]] = []
        monkeypatch.setattr(gene_card, "st", fake)
        monkeypatch.setattr(gene_card, "_language", lambda: "es")
        assert gene_card.structure_viewer is not None
        monkeypatch.setattr(
            gene_card.structure_viewer, "render",
            lambda pdb, **kwargs: viewer_calls.append((str(pdb), kwargs)),
        )
        gene_card.render_gene_card(
            gene, "demo", data_root or tmp_path,
            sequences=sequences, key_prefix=prefix,
        )
        return fake, viewer_calls

    return _run


# --------------------------------------------------------------- extracción
def test_extract_gene_dna_plus_strand():
    seq = "ATGAAACCCGGGTTT"          # 15 pb
    assert gene_card.extract_gene_dna(seq, 1, 15, "+") == seq
    assert gene_card.extract_gene_dna(seq, 4, 6, "+") == "AAA"
    assert gene_card.extract_gene_dna(seq, 12, 15, 1) == "GTTT"


def test_extract_gene_dna_minus_strand_is_reverse_complement():
    assert gene_card.extract_gene_dna("ATGAAACCC", 1, 9, "-") == "GGGTTTCAT"
    assert gene_card.extract_gene_dna("AAAACCCCGGGG", 2, 6, -1) == "GGTTT"


def test_extract_gene_dna_rejects_invalid_coordinates_or_empty_sequence():
    seq = "ATGC"
    assert gene_card.extract_gene_dna("", 1, 3, "+") == ""
    assert gene_card.extract_gene_dna(seq, 0, 3, "+") == ""
    assert gene_card.extract_gene_dna(seq, 3, 2, "+") == ""
    assert gene_card.extract_gene_dna(seq, 1, 99, "+") == ""


# --------------------------------------------------------------- traducción
def test_translate_dna_standard_code_and_edge_cases():
    assert gene_card.translate_dna("ATGAAATAA") == "MK*"
    assert gene_card.translate_dna("atgaaa") == "MK"
    assert gene_card.translate_dna("AUGAAA") == "MK"       # ARN -> ADN
    assert gene_card.translate_dna("ATGNNN") == "MX"       # ambiguo -> X
    assert gene_card.translate_dna("ATGA") == "M"          # codón incompleto
    assert gene_card.translate_dna("ATGTAAGCC", to_stop=True) == "M"
    assert gene_card.translate_dna("") == ""


# --------------------------------------------------- A/B/C sin inventar
def test_abc_path_only_uses_declared_values():
    assert gene_card._abc_path({}) == ""
    assert gene_card._abc_path({"path": "b"}) == "B"
    assert gene_card._abc_path({"label": "B_uncharacterized_homolog"}) == "B"
    assert gene_card._abc_path({"path": "not_evaluated"}) == "not_evaluated"
    assert gene_card._abc_path({"path": "weird"}) == "weird"


# ------------------------------------------------------- sin secuencia/PDB
def test_render_without_sequence_reports_it(run_card):
    fake, viewer = run_card(_gene(), sequences=None)
    text = "\n".join(fake.infos)
    assert "Secuencia de ADN no disponible" in text
    assert "No se genera traducción" in text
    assert "Estructura PDB no disponible" in text
    assert fake.downloads == []
    assert fake.buttons == []
    assert viewer == []


def test_render_without_pdb_reports_it(run_card):
    fake, viewer = run_card(_gene(), sequences={"c1": "ATGAAATAA"})
    files = [item["file_name"] for item in fake.downloads]
    assert any(name.endswith(".fna") for name in files)
    assert any(name.endswith(".faa") for name in files)
    assert not any(name.endswith(".pdb") for name in files)
    assert any("Estructura PDB no disponible" in note for note in fake.infos)
    assert fake.buttons == []
    assert viewer == []
    # La traducción es automática; debe declararse como no validada.
    assert "no está validada experimentalmente" in _all_text(fake)


def test_render_with_sequence_reports_invalid_coordinates(run_card):
    fake, _ = run_card(_gene(end=99), sequences={"c1": "ATGC"})
    assert any("exceden el contig" in note for note in fake.infos)
    assert fake.downloads == []


# ---------------------------------------------------------- PDB disponible
def test_render_with_pdb_offers_download_and_viewer(run_card, tmp_path: Path):
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n",
        encoding="utf-8",
    )
    fake, viewer = run_card(
        _gene(pdb=str(pdb)), sequences={"c1": "ATGAAATAA"}, click=True,
    )
    pdb_downloads = [d for d in fake.downloads
                     if str(d["file_name"]).endswith(".pdb")]
    assert len(pdb_downloads) == 1
    assert pdb_downloads[0]["data"] == pdb.read_bytes()
    assert viewer and viewer[0][0] == str(pdb)


# ----------------------------------------- metadatos declarados vs inventados
def test_render_shows_declared_metadata(run_card):
    gene = _gene(
        path="A", label="A_named_function", confidence="medium",
        status="consensus", evidence_strength="moderate", ec="2.7.2.4",
        complex_name="ATP synthase (F1Fo)",
        complex_status="complex_suspected", complex_role="subunit",
    )
    fake, _ = run_card(gene, sequences={"c1": "ATGAAATAA"})
    text = _all_text(fake)
    for expected in (
        "A_named_function", "medium", "consensus", "moderate", "2.7.2.4",
        "ATP synthase (F1Fo)", "complex_suspected", "subunit",
    ):
        assert expected in text, expected
    assert "Homología Swiss-Prot" in text


def test_render_does_not_invent_missing_confidence_or_abcs(run_card):
    gene = _gene(confidence=None, status=None)
    fake, _ = run_card(gene, sequences={"c1": "ATGAAATAA"})
    text = _all_text(fake)
    assert "A/B/C no declarado" in text
    assert "Confianza/estado no declarados" in text
    assert "high" not in text.lower()
