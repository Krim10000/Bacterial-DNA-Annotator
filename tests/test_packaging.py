"""Validación offline del empaquetado Bioconda + Galaxy.

No usa red ni ejecuta conda/planemo: solo parsea los artefactos y comprueba
que sean coherentes con ``pyproject.toml`` y con el CLI real del pipeline.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

try:  # Python >= 3.11
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
META = ROOT / "packaging" / "bioconda" / "meta.yaml"
BUILD = ROOT / "packaging" / "bioconda" / "build.sh"
GALAXY = ROOT / "packaging" / "galaxy" / "bacterial_dna_annotator.xml"
DOC = ROOT / "docs" / "PACKAGING.md"

GALAXY_FLAGS = (
    "--outdir",
    "--threads",
    "--no-trna",
    "--no-rrna",
    "--no-functional",
    "--no-membrane",
    "--no-oric",
    "--no-metabolic",
)


def _pyproject() -> dict:
    with PYPROJECT.open("rb") as fh:
        return tomllib.load(fh)


def _meta() -> dict:
    with META.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _dist_name(spec: str) -> str:
    """Nombre de distribución de una especificación PEP 508/conda."""
    return re.split(r"[\s<>=!~;\[\]]", spec.strip(), maxsplit=1)[0].lower().replace(
        "_", "-"
    ).replace(".", "-")


def test_meta_yaml_parses_and_matches_pyproject():
    project = _pyproject()["project"]
    meta = _meta()

    assert meta["package"]["name"] == project["name"]
    assert str(meta["package"]["version"]) == project["version"]

    source = meta["source"]
    source_root = (META.parent / source["path"]).resolve()
    assert source_root == ROOT
    assert (source_root / "pyproject.toml").is_file()
    # Sin red no se puede calcular un SHA256: la receta local no debe llevarlo.
    assert "sha256" not in source
    assert not re.search(r"\b[0-9a-f]{64}\b", META.read_text(encoding="utf-8"), re.I)


def test_meta_yaml_entry_point_and_license():
    project = _pyproject()["project"]
    meta = _meta()

    assert meta["build"]["noarch"] == "python"
    assert meta["build"]["script"] == "build.sh"
    expected_entry_point = (
        f"dna-annotate = {project['scripts']['dna-annotate']}"
    )
    assert expected_entry_point in meta["build"]["entry_points"]

    about = meta["about"]
    assert about["license"] == project["license"]["text"] == "MIT"
    assert about["license_family"] == "MIT"
    assert (ROOT / about["license_file"]).is_file()
    assert (ROOT / "LICENSE").read_text(encoding="utf-8").startswith("MIT License")

    assert "dna-annotate --help" in meta["test"]["commands"]


def test_meta_yaml_declares_real_runtime_dependencies():
    project = _pyproject()["project"]
    meta = _meta()

    declared = {_dist_name(spec) for spec in meta["requirements"]["run"]}
    for dependency in project["dependencies"]:
        assert _dist_name(dependency) in declared, dependency

    assert "python" in declared
    # Importado por webapp.esmfold_tools y webapp.structure_evidence, aunque
    # no esté declarado en pyproject.toml.
    assert "requests" in declared


def test_build_script_is_valid_and_handles_ghost():
    assert BUILD.is_file()
    result = subprocess.run(
        ["bash", "-n", str(BUILD)], capture_output=True, text=True, timeout=30
    )
    assert result.returncode == 0, result.stderr

    text = BUILD.read_text(encoding="utf-8")
    assert "pip install . --no-deps --no-build-isolation" in text
    assert "scripts/ghost/ghost" in text
    assert "SITE_DIR" in text
    assert os.access(BUILD, os.X_OK)


def test_galaxy_xml_parses_and_has_required_interface():
    project = _pyproject()["project"]
    root = ET.parse(GALAXY).getroot()

    assert root.tag == "tool"
    assert root.get("id") == "bacterial_dna_annotator"
    assert root.get("version") == project["version"]

    inputs = root.findall("./inputs/param")
    fasta_inputs = [
        p for p in inputs
        if p.get("name") == "input" and p.get("type") == "data"
    ]
    assert len(fasta_inputs) == 1
    assert fasta_inputs[0].get("format") == "fasta"

    outputs = {d.get("name"): d for d in root.findall("./outputs/data")}
    assert {"gff", "tsv", "json"} <= set(outputs)
    assert outputs["gff"].get("format") == "gff3"
    assert outputs["tsv"].get("format") == "tsv"
    assert outputs["json"].get("format") == "json"
    assert outputs["gff"].get("from_work_dir") == "out/input.gff3"
    assert outputs["tsv"].get("from_work_dir") == "out/input.candidates.tsv"
    assert outputs["json"].get("from_work_dir") == "out/input.summary.json"

    requirements = root.findall("./requirements/requirement")
    assert any(
        (r.text or "").strip() == project["name"]
        and r.get("version") == project["version"]
        for r in requirements
    )

    command = root.findtext("command") or ""
    assert "python -m webapp.pipeline" in command
    assert "input.fna" in command
    assert "--outdir out" in command
    for flag in GALAXY_FLAGS:
        assert flag in command, flag


def test_galaxy_flags_exist_in_the_real_cli():
    result = subprocess.run(
        [sys.executable, "-m", "webapp.pipeline", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    for flag in GALAXY_FLAGS:
        assert flag in result.stdout, flag


def test_galaxy_output_names_match_pipeline():
    pipeline = (ROOT / "webapp" / "pipeline.py").read_text(encoding="utf-8")
    for name in ('f"{g}.gff3"', 'f"{g}.candidates.tsv"', 'f"{g}.summary.json"'):
        assert name in pipeline, name


def test_packaging_doc_records_offline_status():
    text = DOC.read_text(encoding="utf-8").lower()
    assert "bioconda" in text
    assert "galaxy" in text
    assert "sin publicar" in text
    assert "no se ha publicado" in text
