"""Pruebas del árbol taxonómico local (``webapp/taxonomy_tree_ui.py``).

Usan JSON sintéticos en ``tmp_path``: no hay red, GPU ni taxonomía real. Se
comprueban los linajes parciales, el respaldo NCBI, la deduplicación
full/piloto, el caso sin taxonomía y el resaltado textual.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

import taxonomy_tree_ui


GTDB_ECOLI = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
    "s__Escherichia coli"
)
GTDB_KLEB = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;"
    "s__Klebsiella pneumoniae"
)
GTDB_BACILLUS = (
    "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;"
    "f__Bacillaceae;g__Bacillus;s__Bacillus subtilis"
)
GTDB_SALMONELLA = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;"
    "s__Salmonella enterica"
)


def _data_root(tmp_path: Path) -> Path:
    root = tmp_path / "data"
    (root / "taxonomy" / "full").mkdir(parents=True)
    return root


def _write_record(
    data_root: Path,
    genome: str,
    data: dict,
    *,
    layer: str = "full",
) -> Path:
    directory = data_root / "taxonomy" / ("full" if layer == "full" else "")
    path = directory / f"{genome}.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _record(genome: str, lineage: str, **overrides) -> dict:
    data = {
        "genome": genome,
        "gtdb_lineage": lineage,
        "species": "",
        "confidence": "high",
        "status": "assigned",
        "ani": None,
        "closest_ref": None,
    }
    data.update(overrides)
    return data


@pytest.fixture()
def three_genomes(tmp_path: Path) -> Path:
    data = _data_root(tmp_path)
    _write_record(data, "GCF_A.1", _record(
        "GCF_A.1", GTDB_ECOLI, species="Escherichia coli", ani=98.7,
        closest_ref="GCF_REF.1",
    ))
    _write_record(data, "GCF_B.1", _record(
        "GCF_B.1", GTDB_KLEB, species="Klebsiella pneumoniae", ani=96.1,
    ))
    _write_record(data, "GCF_C.1", _record(
        "GCF_C.1", GTDB_BACILLUS, species="Bacillus subtilis",
    ))
    return data


# ------------------------------------------------------------------- parser
def test_parse_lineage_ignores_placeholders_and_unknown_tokens():
    parsed = taxonomy_tree_ui.parse_lineage(
        "d__Bacteria;p__;x__Ignorado;c__Gammaproteobacteria;"
        "g__Escherichia;s__Escherichia coli"
    )
    assert parsed == [
        ("domain", "Bacteria"),
        ("class", "Gammaproteobacteria"),
        ("genus", "Escherichia"),
        ("species", "Escherichia coli"),
    ]
    assert taxonomy_tree_ui.parse_lineage("unknown") == []
    assert taxonomy_tree_ui.parse_lineage(None) == []


def test_parse_lineage_accepts_mapping_and_plain_names():
    assert taxonomy_tree_ui.parse_lineage(
        {"domain": "Bacteria", "genus": "Escherichia", "species": ""}
    ) == [("domain", "Bacteria"), ("genus", "Escherichia")]
    assert taxonomy_tree_ui.parse_lineage(
        "Bacteria;Pseudomonadota;Gammaproteobacteria"
    ) == [
        ("domain", "Bacteria"),
        ("phylum", "Pseudomonadota"),
        ("class", "Gammaproteobacteria"),
    ]


# --------------------------------------------------------------- load_lineages
def test_load_lineages_parses_full_records(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)

    assert list(lineages) == ["GCF_A.1", "GCF_B.1", "GCF_C.1"]
    first = lineages["GCF_A.1"]
    assert list(first["ranks"]) == list(taxonomy_tree_ui.RANKS)
    assert first["ranks"]["domain"] == "Bacteria"
    assert first["ranks"]["species"] == "Escherichia coli"
    assert first["lineage"] == GTDB_ECOLI
    assert first["species"] == "Escherichia coli"
    assert first["method"] == "gtdb"
    assert first["species_from_json"] is False
    assert first["confidence"] == "high"
    assert first["status"] == "assigned"
    assert first["ani"] == 98.7
    assert first["closest_ref"] == "GCF_REF.1"
    assert first["source"] == "taxonomy/full/GCF_A.1.json"
    assert first["layer"] == "full"
    assert first["n_ranks"] == len(taxonomy_tree_ui.RANKS)
    assert first["missing_ranks"] == ()
    assert first["external"] is False


def test_load_lineages_accepts_project_root(three_genomes: Path):
    """``data_root`` puede ser el proyecto que contiene ``data/``."""
    lineages = taxonomy_tree_ui.load_lineages(three_genomes.parent)
    assert len(lineages) == 3
    assert lineages["GCF_A.1"]["source"] == "taxonomy/full/GCF_A.1.json"


def test_load_lineages_falls_back_to_ncbi_and_declared_species(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_NCBI.1", {
        "genome": "GCF_NCBI.1",
        "gtdb_lineage": "unknown",
        "ncbi_lineage": (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__"
        ),
        "species": "Escherichia coli",
        "confidence": "medium",
    })
    _write_record(data, "GCF_SPECIES.1", {
        "genome": "GCF_SPECIES.1",
        "gtdb_lineage": "unknown",
        "species": "Candidatus Solo especie",
    })

    lineages = taxonomy_tree_ui.load_lineages(data)

    ncbi = lineages["GCF_NCBI.1"]
    assert ncbi["method"] == "ncbi"
    assert ncbi["species"] == "Escherichia coli"
    assert ncbi["species_from_json"] is True
    assert ncbi["ranks"]["species"] == "Escherichia coli"
    assert ncbi["n_ranks"] == 7

    species_only = lineages["GCF_SPECIES.1"]
    assert species_only["method"] == "species_json"
    assert species_only["ranks"] == {"species": "Candidatus Solo especie"}
    assert species_only["n_ranks"] == 1
    assert species_only["missing_ranks"] == tuple(
        rank for rank in taxonomy_tree_ui.RANKS if rank != "species"
    )


def test_load_lineages_skips_aggregates_invalid_json_and_duplicates(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_FULL.1", _record(
        "GCF_FULL.1", GTDB_ECOLI, species="Escherichia coli",
    ))
    # El piloto duplica el genoma: debe ganar el registro de full/.
    _write_record(data, "GCF_FULL.1", _record(
        "GCF_FULL.1", GTDB_BACILLUS, species="Bacillus subtilis",
    ), layer="pilot")
    # Genoma solo piloto: sí se agrega.
    _write_record(data, "GCF_PILOT.1", _record(
        "GCF_PILOT.1", GTDB_BACILLUS, species="Bacillus subtilis",
    ), layer="pilot")
    # Resumen agregado: no es un registro taxonómico.
    (data / "taxonomy" / "full" / "summary.json").write_text(
        json.dumps({"n_genomes": 400, "coverage": {}}), encoding="utf-8"
    )
    (data / "taxonomy" / "full" / "broken.json").write_text(
        "{not valid json", encoding="utf-8"
    )

    lineages = taxonomy_tree_ui.load_lineages(data)

    assert list(lineages) == ["GCF_FULL.1", "GCF_PILOT.1"]
    assert lineages["GCF_FULL.1"]["lineage"] == GTDB_ECOLI
    assert lineages["GCF_FULL.1"]["layer"] == "full"
    assert lineages["GCF_PILOT.1"]["layer"] == "pilot"

    info = taxonomy_tree_ui.inventory(data)
    assert info["available"] is True
    assert info["n_full"] == 1
    assert info["n_pilot"] == 1
    assert info["n_lineages"] == 2
    assert info["reason"] == ""


def test_load_lineages_without_taxonomy(tmp_path: Path):
    empty = tmp_path / "empty"
    empty.mkdir()

    assert taxonomy_tree_ui.load_lineages(empty) == {}

    info = taxonomy_tree_ui.inventory(empty)
    assert info["available"] is False
    assert "taxonom" in info["reason"].lower()
    assert info["n_lineages"] == 0

    tree = taxonomy_tree_ui.build_tree({})
    assert tree["n_genomes"] == 0
    assert tree["children"] == []
    assert taxonomy_tree_ui.find_path(tree, genome="GCF_X.1") == []


# ------------------------------------------------------------------ build_tree
def test_build_tree_creates_full_hierarchy(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(lineages)

    assert tree["name"] == taxonomy_tree_ui.ROOT_NAME
    assert tree["rank"] == "root"
    assert tree["n_genomes"] == 3
    assert tree["n_external"] == 0
    assert [node["name"] for node in tree["children"]] == ["Bacteria"]

    bacteria = tree["children"][0]
    assert bacteria["rank"] == "domain"
    assert bacteria["n_genomes"] == 3
    assert [node["name"] for node in bacteria["children"]] == [
        "Bacillota", "Pseudomonadota",
    ]

    pseudo = next(node for node in bacteria["children"]
                  if node["name"] == "Pseudomonadota")
    assert pseudo["n_genomes"] == 2
    assert [node["name"] for node in pseudo["children"]] == [
        "Gammaproteobacteria",
    ]
    enterobacterales = pseudo["children"][0]["children"][0]
    assert enterobacterales["name"] == "Enterobacterales"
    assert {node["name"] for node in enterobacterales["children"]} == {
        "Enterobacteriaceae",
    }
    family = enterobacterales["children"][0]
    assert {node["name"] for node in family["children"]} == {
        "Escherichia", "Klebsiella",
    }

    ecoli = taxonomy_tree_ui.find_path(tree, genome="GCF_A.1")
    assert ecoli == [
        "Bacteria", "Pseudomonadota", "Gammaproteobacteria",
        "Enterobacterales", "Enterobacteriaceae", "Escherichia",
        "Escherichia coli",
    ]


def test_build_tree_does_not_invent_missing_ranks(tmp_path: Path):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_PART.1", {
        "genome": "GCF_PART.1",
        "ncbi_lineage": (
            "d__Bacteria;p__Nitrospirota;c__Nitrospiria;"
            "o__;f__;g__Candidatus Magnetomonas;s__"
        ),
        "species": "Candidatus Magnetomonas plexicatena",
    })

    lineages = taxonomy_tree_ui.load_lineages(data)
    record = lineages["GCF_PART.1"]
    assert "order" not in record["ranks"]
    assert "family" not in record["ranks"]
    assert "order" in record["missing_ranks"]

    tree = taxonomy_tree_ui.build_tree(lineages)
    path = taxonomy_tree_ui.find_path(tree, genome="GCF_PART.1")
    assert path == [
        "Bacteria", "Nitrospirota", "Nitrospiria",
        "Candidatus Magnetomonas",
        "Candidatus Magnetomonas plexicatena",
    ]
    all_names = []

    def collect(node):
        all_names.append(node["name"])
        for child in node["children"]:
            collect(child)

    collect(tree)
    assert "unknown" not in {name.casefold() for name in all_names}
    assert "" not in all_names


def test_build_tree_dedupes_shared_ancestors(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(lineages)
    # Solo un nodo Bacteria (dominio) aunque haya tres genomas.
    assert [node["name"] for node in tree["children"]] == ["Bacteria"]
    domain = tree["children"][0]
    assert [node["name"] for node in domain["children"]] == [
        "Bacillota", "Pseudomonadota",
    ]
    all_domains = []

    def collect(node):
        if node["rank"] == "domain":
            all_domains.append(node["name"])
        for child in node["children"]:
            collect(child)

    collect(tree)
    assert all_domains == ["Bacteria"]
    pseudo = next(node for node in domain["children"]
                  if node["name"] == "Pseudomonadota")
    assert len(pseudo["genomes"]) == 0
    assert len(pseudo["children"]) == 1


# ------------------------------------------------------------------- find_path
def test_find_path_by_genome_and_lineage(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(lineages)

    assert taxonomy_tree_ui.find_path(tree, genome="GCF_B.1")[-1] == (
        "Klebsiella pneumoniae"
    )
    assert taxonomy_tree_ui.find_path(tree, lineage=GTDB_ECOLI) == (
        taxonomy_tree_ui.find_path(tree, genome="GCF_A.1")
    )
    # Coincidencia sin distinguir mayúsculas.
    assert taxonomy_tree_ui.find_path(
        tree, lineage=GTDB_ECOLI.lower()
    )[-1] == "Escherichia coli"
    # Diccionario por rango sin dominio: se ubica por el primer nombre.
    assert taxonomy_tree_ui.find_path(
        tree, lineage={"species": "Escherichia coli"}
    )[-1] == "Escherichia coli"
    # Prefijo parcial: devuelve lo que sí coincide.
    assert taxonomy_tree_ui.find_path(
        tree, lineage="d__Bacteria;p__Pseudomonadota;c__NoExiste"
    ) == ["Bacteria", "Pseudomonadota"]
    assert taxonomy_tree_ui.find_path(tree, genome="GCF_MISSING.1") == []
    assert taxonomy_tree_ui.find_path(tree, lineage="d__NoExiste") == []


# --------------------------------------------------------------- rama externa
def test_build_tree_adds_external_new_lineage(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(
        lineages,
        new_genome="GCF_NEW.1",
        new_lineage=GTDB_SALMONELLA,
    )

    assert tree["n_genomes"] == 3
    assert tree["n_external"] == 1
    assert tree["has_external"] is True
    assert tree["new_genome"] == "GCF_NEW.1"
    assert tree["new_lineage"] == GTDB_SALMONELLA
    path = taxonomy_tree_ui.find_path(tree, genome="GCF_NEW.1")
    assert path == [
        "Bacteria", "Pseudomonadota", "Gammaproteobacteria",
        "Enterobacterales", "Enterobacteriaceae", "Salmonella",
        "Salmonella enterica",
    ]
    assert taxonomy_tree_ui.find_path(tree, lineage=GTDB_SALMONELLA) == path


def test_build_tree_does_not_duplicate_existing_new_genome(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(
        lineages,
        new_genome="GCF_A.1",
        new_lineage=GTDB_BACILLUS,
    )
    assert tree["n_genomes"] == 3
    assert tree["n_external"] == 0
    assert tree["has_external"] is False
    assert taxonomy_tree_ui.find_path(tree, genome="GCF_A.1")[-1] == (
        "Escherichia coli"
    )


# --------------------------------------------------------------- texto/streamlit
def test_tree_lines_highlights_selected_and_new(three_genomes: Path):
    lineages = taxonomy_tree_ui.load_lineages(three_genomes)
    tree = taxonomy_tree_ui.build_tree(
        lineages,
        new_genome="GCF_NEW.1",
        new_lineage=GTDB_SALMONELLA,
    )

    text = taxonomy_tree_ui.tree_text(
        tree,
        selected="GCF_A.1",
        new_genome="GCF_NEW.1",
        new_lineage=GTDB_SALMONELLA,
    )
    assert "◀ seleccionado" in text
    assert "★ nuevo" in text
    assert "Escherichia coli" in text
    assert text.count("Taxonomía") == 1

    pruned = taxonomy_tree_ui.tree_text(
        tree,
        selected="GCF_A.1",
        new_genome="GCF_NEW.1",
        new_lineage=GTDB_SALMONELLA,
        only_paths=True,
    )
    assert "Bacillota" not in pruned
    assert "Escherichia coli" in pruned
    assert "Salmonella enterica" in pruned
    assert "rama" in pruned

    truncated = taxonomy_tree_ui.tree_lines(tree, max_lines=3)
    assert len(truncated) == 4
    assert truncated[-1] == "… (árbol truncado)"


class _FakeCtx:
    def __init__(self, calls):
        self.calls = calls

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            return _FakeCtx(self.calls)

        return recorder


class _FakeStreamlit:
    """Doble mínimo de Streamlit; permite ejecutar ``render`` sin runtime."""

    def __init__(self, selected=None, view="Ruta seleccionada"):
        self.calls: list = []
        self.selected = selected
        self.view = view

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "columns":
                count = int(args[0]) if args else 1
                return [_FakeCtx(self.calls) for _ in range(count)]
            if name == "selectbox":
                if self.selected is not None:
                    return self.selected
                options = list(args[1]) if len(args) > 1 else []
                return options[0] if options else None
            if name == "radio":
                return self.view
            return _FakeCtx(self.calls)

        return recorder


def test_render_without_taxonomy_shows_empty_state(tmp_path: Path, monkeypatch):
    empty = tmp_path / "empty"
    empty.mkdir()
    fake = _FakeStreamlit()
    monkeypatch.setitem(sys.modules, "streamlit", fake)

    result = taxonomy_tree_ui.render(empty)

    assert result["available"] is False
    assert result["n_genomes"] == 0
    assert result["path"] == []
    names = [call[0] for call in fake.calls]
    assert "info" in names
    texts = " ".join(
        str(call[1][0])
        for call in fake.calls
        if call[0] in {"warning", "caption", "info"} and call[1]
    )
    assert "precalculada" in texts.casefold()
    assert "No se consulta la red" in texts or "no se consulta la red" in texts


def test_render_with_data_highlights_selected_and_new(
    three_genomes: Path, monkeypatch
):
    fake = _FakeStreamlit(selected="GCF_A.1")
    monkeypatch.setitem(sys.modules, "streamlit", fake)

    result = taxonomy_tree_ui.render(
        three_genomes,
        new_genome="GCF_NEW.1",
        new_lineage=GTDB_SALMONELLA,
    )

    assert result["available"] is True
    assert result["n_genomes"] == 3
    assert result["selected"] == "GCF_A.1"
    assert result["path"][-1] == "Escherichia coli"
    assert result["new_path"][-1] == "Salmonella enterica"
    assert result["external"] is False
    assert result["tree"]["n_external"] == 1

    selectboxes = [call for call in fake.calls if call[0] == "selectbox"]
    assert len(selectboxes) == 1
    assert len(selectboxes[0][1][1]) == 4  # 3 locales + el genoma nuevo
    assert "GCF_NEW.1" in selectboxes[0][1][1]

    metrics = [call for call in fake.calls if call[0] == "metric"]
    assert any(call[1][0] == "Genoma" and call[1][1] == "GCF_A.1"
               for call in metrics)
    assert any(call[1][0] == "Especie" and call[1][1] == "Escherichia coli"
               for call in metrics)

    codes = [str(call[1][0]) for call in fake.calls if call[0] == "code"]
    assert any("◀ seleccionado" in code for code in codes)
    assert any("★ nuevo" in code for code in codes)
    assert any("Escherichia coli" in code for code in codes)
    assert any("Salmonella enterica" in code for code in codes)

    frames = [call for call in fake.calls if call[0] == "dataframe"]
    assert frames
    ranks = [row["Rango"] for row in frames[0][1][0]]
    assert ranks[0] == "Dominio"
    assert ranks[-1] == "Especie"


def test_render_new_genome_without_lineage_warns(tmp_path: Path, monkeypatch):
    data = _data_root(tmp_path)
    _write_record(data, "GCF_A.1", _record(
        "GCF_A.1", GTDB_ECOLI, species="Escherichia coli",
    ))
    fake = _FakeStreamlit(selected="GCF_A.1")
    monkeypatch.setitem(sys.modules, "streamlit", fake)

    result = taxonomy_tree_ui.render(data, new_genome="GCF_NEW.1")

    assert result["new_path"] == []
    warnings = " ".join(
        str(call[1][0])
        for call in fake.calls
        if call[0] == "warning" and call[1]
    )
    assert "new_lineage" in warnings
