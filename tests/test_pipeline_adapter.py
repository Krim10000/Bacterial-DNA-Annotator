"""Pruebas del adaptador en memoria GUI -> pipeline CLI.

Usan un FASTA mínimo y un ``data_root`` vacío: no deben lanzar, deben devolver
las claves de todas las capas y marcar como ``not_available``/``disabled`` lo
que no se pudo ejecutar. No se invoca DIAMOND ni ninguna base real.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

from webapp import pipeline_adapter
from webapp.pipeline_adapter import run_in_memory

SEQUENCE = "ACGTAGCTAGCTAGGCTA" * 20
MINIMAL_FASTA = (">contig1 test\n" + SEQUENCE + "\n").encode("utf-8")

#: Todas las capas de la ruta unificada que el resultado debe exponer.
EXPECTED_LAYERS = (
    "ncrna",
    "taxonomy",
    "functional",
    "domains",
    "context_c",
    "structural_b",
    "structural_b_filters",
    "complexes",
    "regulation",
    "membrane",
    "oric",
)

#: Apaga las herramientas externas de anotación para que la prueba no dependa
#: de pyrodigal/tRNAscan-SE/barrnap. Las capas opcionales siguen activas.
NO_EXTERNAL_TOOLS = {
    "cds": False,
    "trna": False,
    "rrna": False,
    "ghost": False,
    "deserts": False,
    "esm": False,
    "metabolic": False,
}


def _run(tmp_path: Path, options: dict | None = None, name: str = "mini.fna"):
    outdir = tmp_path / "out"
    data_root = tmp_path / "data"
    data_root.mkdir(exist_ok=True)
    res = run_in_memory(name, MINIMAL_FASTA, outdir, data_root,
                        options=options, threads=1)
    return res, outdir, data_root


def test_minimal_fasta_without_data_root_degrades_without_raising(tmp_path):
    res, outdir, _ = _run(tmp_path, options=NO_EXTERNAL_TOOLS)

    assert res["genome"] == "mini"
    assert res["bp"] == len(SEQUENCE)
    assert res["gff"].startswith("##gff-version 3")
    assert res["report_md"].startswith("# Informe de anotación")

    for layer in EXPECTED_LAYERS:
        assert layer in res
        assert isinstance(res[layer], dict)
        assert "ok" in res[layer]
        assert "status" in res[layer]
        assert "reason" in res[layer]

    # Con data_root vacío ningún resultado se inventa: cada capa declara por
    # qué no pudo ejecutarse (o se desactivó explícitamente).
    assert res["ncrna"]["status"] == "not_available"
    assert res["taxonomy"]["status"] == "not_available"
    assert res["functional"]["status"] == "not_available"
    assert "proteínas" in res["functional"]["reason"]
    assert res["domains"]["status"] == "not_available"
    assert res["context_c"]["status"] == "not_available"
    assert res["structural_b"]["status"] == "not_available"
    assert res["complexes"]["status"] == "not_available"
    assert res["membrane"]["status"] == "not_available"
    assert res["metabolic"] is None
    assert res["metabolic_status"]["status"] == "disabled"

    # Capas con datos disponibles localmente sí se pueden ejecutar sin red.
    assert res["regulation"]["status"] in {"ok", "not_available"}
    assert res["oric"]["status"] in {"unavailable", "not_available", "error"}

    # Salidas base del CLI, escritas por write_outputs.
    assert (outdir / "mini.gff3").is_file()
    assert (outdir / "mini.proteins.faa").is_file()
    assert (outdir / "mini.candidates.tsv").is_file()
    assert (outdir / "mini.report.md").is_file()
    assert (outdir / "mini.summary.json").is_file()
    assert {"gff", "proteins", "candidates", "report", "summary"} <= set(res["outputs"])

    summary = json.loads((outdir / "mini.summary.json").read_text())
    assert summary["genome"] == "mini"
    assert summary["bp"] == len(SEQUENCE)

    # No debe haber artefactos de capas no disponibles (no se inventan).
    assert not (outdir / "mini.ncrna.tsv").exists()
    assert not (outdir / "mini.taxonomy.tsv").exists()
    assert not (outdir / "mini.membrane.tsv").exists()
    assert "ncrna" not in res["outputs"]
    assert "membrane" not in res["outputs"]


def test_flags_negativos_y_gzip_respetan_opciones(tmp_path):
    outdir = tmp_path / "out"
    data_root = tmp_path / "data"
    data_root.mkdir()
    options = {
        "no_cds": True,
        "no_trna": True,
        "no_rrna": True,
        "no_ghost": True,
        "no_deserts": True,
        "no_ncrna": True,
        "no_taxonomy": True,
        "no_metabolic": True,
    }
    res = run_in_memory("mini.fna", gzip.compress(MINIMAL_FASTA), outdir,
                        data_root, options=options, threads=1)

    assert res["genome"] == "mini"
    assert res["options"]["cds"] is False
    assert res["ncrna"]["status"] == "disabled"
    assert res["taxonomy"]["status"] == "disabled"
    assert res["metabolic"] is None
    assert not (outdir / "mini.ncrna.tsv").exists()
    assert not (outdir / "mini.taxonomy.tsv").exists()
    assert "gff" in res["outputs"]


def test_capa_con_dependencia_ausente_no_lanza(monkeypatch, tmp_path):
    def missing_cmsearch(*_args, **_kwargs):
        raise FileNotFoundError("cmsearch no está en PATH")

    monkeypatch.setattr(pipeline_adapter, "run_ncrna", missing_cmsearch)
    res, _, _ = _run(tmp_path, options=NO_EXTERNAL_TOOLS)

    assert res["ncrna"]["ok"] is False
    assert res["ncrna"]["status"] == "not_available"
    assert "cmsearch" in res["ncrna"]["reason"]
    # La capa que falla no impide el resto del pipeline ni el informe.
    assert res["report_md"]
    assert "summary" in res["outputs"]


def test_metabolismo_opcional_sin_datos_no_lanza(tmp_path):
    options = {**NO_EXTERNAL_TOOLS, "metabolic": True}
    res, _, _ = _run(tmp_path, options=options)

    assert res["metabolic"] is None
    assert res["metabolic_status"]["ok"] is False
    assert res["metabolic_status"]["status"] == "not_available"
    assert "proteínas" in res["metabolic_status"]["reason"]
    assert res["layer_status"]["metabolic"] == "not_available"


def test_fallo_de_la_capa_base_devuelve_error_sin_lanzar(monkeypatch, tmp_path):
    def broken_annotate(*_args, **_kwargs):
        raise RuntimeError("pyrodigal falló")

    monkeypatch.setattr(pipeline_adapter, "annotate", broken_annotate)
    res, _, _ = _run(tmp_path, options=NO_EXTERNAL_TOOLS)

    assert res["ok"] is False
    assert res["status"] == "error"
    assert "pyrodigal" in res["reason"]
    assert res["outputs"] == {}
    for layer in EXPECTED_LAYERS:
        assert res[layer]["status"] == "error"


def test_anotacion_base_con_pyrodigal_si_esta_disponible(tmp_path):
    import pytest

    pytest.importorskip("pyrodigal")

    # ORF sintético con codón de inicio/parada; CDS activo, herramientas
    # externas y ORCA apagados para que la prueba sea rápida y local.
    orf = "ATG" + "GCT" * 200 + "TAA"
    spacer = "A" * 60
    fasta = (">contig1 test\n" + (spacer + orf + spacer) * 3 + "\n").encode()
    outdir = tmp_path / "out"
    data_root = tmp_path / "data"
    data_root.mkdir()
    options = {
        **NO_EXTERNAL_TOOLS,
        "cds": True,
        "oric": False,
        "regulation": True,
    }
    res = run_in_memory("cds.fna", fasta, outdir, data_root, options=options,
                        threads=1)

    assert res["bp"] > 0
    assert res["gff"].startswith("##gff-version 3")
    assert res["oric"]["status"] == "disabled"
    # Sin Swiss-Prot local, la capa funcional degrada, pero no lanza y deja
    # constancia del motivo en lugar de fabricar hits.
    assert res["functional"]["status"] in {"not_available", "not_evaluated"}
    assert (outdir / "cds.summary.json").is_file()
