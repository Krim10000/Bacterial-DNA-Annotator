"""Validación estática del entorno Docker para el benchmark externo.

No construye imágenes, no usa red y no ejecuta Bakta/Prokka/DFAST. Solo
comprueba que la receta y la documentación existen y declaran los comandos,
herramientas y límites que el benchmark necesita.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "benchmark.Dockerfile"
DOC = ROOT / "docs" / "BENCHMARK_ENV.md"
BENCHMARK = ROOT / "tools" / "benchmark_external.py"


def read(path: Path) -> str:
    assert path.is_file(), f"falta {path}"
    return path.read_text(encoding="utf-8")


def test_files_exist():
    assert DOCKERFILE.is_file()
    assert DOC.is_file()
    assert BENCHMARK.is_file()


def test_dockerfile_is_cpu_and_uses_micromamba():
    text = read(DOCKERFILE)
    assert text.startswith("# syntax=docker/dockerfile:1")
    assert "FROM mambaorg/micromamba:latest" in text
    assert "python=" in text
    assert "USER root" in text
    assert "WORKDIR /work" in text
    assert 'VOLUME ["/data"]' in text
    assert "CMD [" in text
    lower = text.lower()
    for forbidden in ("nvidia", "cuda", "gpu"):
        assert forbidden not in lower, f"la imagen debe ser CPU: {forbidden}"


def test_dockerfile_attempts_all_three_tools_in_separate_envs():
    text = read(DOCKERFILE)
    for tool in ("bakta", "prokka", "dfast"):
        assert f"-n {tool}" in text
        assert f"micromamba create" in text
        assert f"/usr/local/bin/{tool}" in text
        assert f"'exec /usr/local/bin/micromamba run -n {tool} {tool} \"$@\"'" in text
    assert text.count("micromamba create -y -n") == 3
    assert "conda-forge" in text and "bioconda" in text


def test_dockerfile_versions_are_explicitly_unpinned_by_default():
    text = read(DOCKERFILE)
    assert 'ARG BAKTA_VERSION=""' in text
    assert 'ARG PROKKA_VERSION=""' in text
    assert 'ARG DFAST_VERSION=""' in text
    assert "BAKTA_SPEC=" in text
    assert "PROKKA_SPEC=" in text
    assert "DFAST_SPEC=" in text
    # No se inventan versiones fijas en los defaults.
    assert "BAKTA_VERSION=1." not in text
    assert "PROKKA_VERSION=1." not in text
    assert "DFAST_VERSION=1." not in text


def test_dockerfile_does_not_download_databases():
    text = read(DOCKERFILE)
    lower = text.lower()
    for forbidden in ("wget", "curl", "git clone", "bakta_db download",
                      "download the database", "descarga la base"):
        assert forbidden not in lower, f"no debe descargar bases: {forbidden}"
    assert "no se descarga" in lower or "no descarga" in lower


def test_dockerfile_copies_benchmark_and_runs_help():
    text = read(DOCKERFILE)
    assert "COPY tools/benchmark_external.py" in text
    assert "/opt/benchmark/tools/benchmark_external.py" in text
    assert 'CMD ["python", "/opt/benchmark/tools/benchmark_external.py", "--help"]' in text


def test_docs_explain_build_run_and_no_execution_claim():
    text = read(DOC)
    assert "docker build" in text
    assert "docker/benchmark.Dockerfile" in text
    assert "docker run" in text
    assert "python /work/tools/benchmark_external.py" in text
    assert "--genome" in text
    assert "--reference-gff" in text
    assert "--tool all --run" in text
    assert "--bakta-db" in text
    assert "no se construyó ni se ejecutó" in text
    assert "El build requiere red" in text
    assert "No hay build ejecutado aquí" in text


def test_docs_cover_versions_licenses_and_limits():
    text = read(DOC)
    assert "versiones quedan sin pinear" in text
    assert "última versión" in text
    assert "GPL-3.0" in text
    assert "sin licencia declarada" in text
    assert "Bakta" in text and "Prokka" in text and "DFAST" in text
    assert "CPU" in text
    assert "bases de datos" in text
    assert "no evalúa función" in text or "no evalúa función" in text.lower()


def test_docs_commands_match_benchmark_cli():
    benchmark = read(BENCHMARK)
    for flag in ("--genome", "--reference-gff", "--tool", "--run",
                 "--bakta-db", "--md", "--json"):
        assert flag in benchmark, f"benchmark_external.py no declara {flag}"
        assert flag in read(DOC), f"la documentación no usa {flag}"


def test_docs_do_not_claim_tools_were_executed():
    text = read(DOC).lower()
    assert "no se construyó ni se ejecutó" in text
    # No debe aparecer una afirmación positiva de ejecución.
    for phrase in ("se ejecutó bakta", "se ejecutó prokka", "se ejecutó dfast",
                   "construcción exitosa", "build exitoso"):
        assert phrase not in text
