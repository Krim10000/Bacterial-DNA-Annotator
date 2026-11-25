"""Pruebas de ejecución real-posible del benchmark externo.

Todo es local y sintético:

- no se descargan imágenes, bases ni datos de red;
- no se ejecutan Bakta/Prokka/DFAST reales;
- los binarios son scripts Python falsos en ``tmp_path``.

Se cubren detección/ausencia, validación de Bakta DB, gold pequeño desde
GFF3/GBFF local, manifiesto ``--run`` reproducible, SHA256/tiempos y salida de
informes en ``reports/external_benchmark/``.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_DIR = ROOT / "benchmarks" / "external"
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import external_benchmark as ext  # noqa: E402
import make_gold  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers sintéticos
# ---------------------------------------------------------------------------
def write_fasta(path: Path, *, contig: str = "ctg", length: int = 4000) -> Path:
    sequence = ("ACGT" * ((length // 4) + 1))[:length]
    path.write_text(
        f">{contig} synthetic contig\n{sequence}\n", encoding="utf-8"
    )
    return path


def write_gff(path: Path, features: list[tuple]) -> Path:
    lines = ["##gff-version 3", "##sequence-region ctg 1 4000"]
    for seqid, ftype, start, end, strand, attrs in features:
        lines.append(
            f"{seqid}\ttest\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_genbank(path: Path, *, length: int = 600) -> Path:
    sequence = ("ACGT" * ((length // 4) + 1))[:length]
    text = f"""LOCUS       {"ctg":<16}{length:>12} bp    DNA     linear   BCT 01-JAN-2026
DEFINITION  synthetic record.
ACCESSION   ctg
VERSION     ctg.1
FEATURES             Location/Qualifiers
     source          1..{length}
                     /organism="synthetic bacterium"
     CDS             1..90
                     /product="DNA polymerase"
                     /gene="polA"
                     /locus_tag="CTG_0001"
     CDS             complement(120..300)
                     /product="hypothetical protein"
                     /locus_tag="CTG_0002"
     tRNA            400..470
                     /product="tRNA-Ala"
ORIGIN
        1 {sequence}
//
"""
    path.write_text(text, encoding="utf-8")
    return path


def synthetic_inputs(tmp_path: Path) -> tuple[Path, Path]:
    genome = write_fasta(tmp_path / "genome.fna")
    reference = write_gff(
        tmp_path / "reference.gff3",
        [
            ("ctg", "CDS", 101, 400, "+", "ID=c1;product=DNA polymerase"),
            ("ctg", "CDS", 501, 800, "-", "ID=c2;product=hypothetical protein"),
            ("ctg", "tRNA", 900, 970, "+", "ID=t1;product=tRNA-Ala"),
            ("ctg", "rRNA", 1100, 1400, "+", "ID=r1;product=16S ribosomal RNA"),
        ],
    )
    return genome, reference


def write_fake_tool(
    path: Path,
    *,
    binary_name: str = "bakta",
    version: str = "bakta 1.12.2",
    write_output: bool = True,
    fail: bool = False,
) -> Path:
    """Binario falso; escribe un GFF3 compatible con el genoma sintético."""
    path.mkdir(parents=True, exist_ok=True)
    target = path / binary_name
    body = (
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "if '--version' in args:\n"
        f"    print({version!r})\n"
        "    sys.exit(0)\n"
    )
    if fail:
        body += "sys.exit(3)\n"
    else:
        body += (
            "out = None\n"
            "for flag in ('--output', '--outdir', '--out'):\n"
            "    if flag in args:\n"
            "        out = args[args.index(flag) + 1]\n"
            "if out is None:\n"
            "    sys.exit(2)\n"
            "outdir = pathlib.Path(out)\n"
            "outdir.mkdir(parents=True, exist_ok=True)\n"
        )
        if write_output:
            body += (
                "(outdir / 'result.gff3').write_text(\n"
                "    '##gff-version 3\\n'\n"
                "    'ctg\\ttest\\tCDS\\t101\\t400\\t.\\t+\\t.\\tID=q1;product=DNA polymerase\\n'\n"
                "    'ctg\\ttest\\ttRNA\\t900\\t970\\t.\\t+\\t.\\tID=q2;product=tRNA-Ala\\n'\n"
                ")\n"
            )
        body += "sys.exit(0)\n"
    target.write_text(body, encoding="utf-8")
    target.chmod(0o755)
    return target


def make_bakta_db(path: Path, *, version: str = "5.1", nested: bool = False) -> Path:
    db = path / "db" if nested else path
    db.mkdir(parents=True, exist_ok=True)
    (db / "version.json").write_text(
        json.dumps({"tool": "bakta", "version": version}) + "\n", encoding="utf-8"
    )
    return path


# ---------------------------------------------------------------------------
# Detección de binarios y ausencia
# ---------------------------------------------------------------------------
def test_all_tools_missing_are_not_available_without_execution(
    tmp_path: Path, monkeypatch
):
    genome, reference = synthetic_inputs(tmp_path)
    called: list[object] = []

    def forbidden(*args, **kwargs):
        called.append((args, kwargs))
        raise AssertionError("no se debe ejecutar un binario ausente")

    monkeypatch.setattr(ext.subprocess, "run", forbidden)
    data = ext.run_benchmark(
        genome,
        reference,
        tools="all",
        run=True,
        which=lambda _name: None,
        environ={},
    )

    assert data["not_available"] == ["bakta", "prokka", "dfast"]
    assert data["environment"]["available_tools"] == []
    assert data["environment"]["bakta_db"]["status"] == "not_found"
    for name in ext.TOOLS:
        result = data["results"][name]
        assert result["status"] == "not_available"
        assert result["comparison"] is None
        assert result["execution"] is None
        assert "shutil.which" in result["reason"]
        assert result["query_gff"] is None
    assert called == []
    # El JSON no inventa métricas para herramientas ausentes.
    assert all(value["comparison"] is None for value in data["results"].values())
    assert "F1" not in json.dumps(data["results"])


def test_detect_tool_requires_executable_file(tmp_path: Path):
    plain = tmp_path / "prokka"
    plain.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    plain.chmod(0o644)
    missing_exec = ext.detect_tool("prokka", which=lambda _name: str(plain))
    assert missing_exec["available"] is False
    assert missing_exec["executable"] is False
    assert "ejecutable" in missing_exec["reason"]

    plain.chmod(0o755)
    ready = ext.detect_tool("prokka", which=lambda _name: str(plain))
    assert ready["available"] is True
    assert ready["executable"] is True
    assert ready["reason"] == ""


def test_detect_tool_missing_is_honest():
    missing = ext.detect_tool("dfast", which=lambda _name: None)
    assert missing["available"] is False
    assert missing["path"] is None
    assert "shutil.which" in missing["reason"]


# ---------------------------------------------------------------------------
# Bakta DB
# ---------------------------------------------------------------------------
def test_detect_bakta_db_from_env_with_version_and_sha(tmp_path: Path):
    db = make_bakta_db(tmp_path / "bakta_db", version="5.1.2")
    detected = ext.detect_bakta_db(environ={"BAKTA_DB": str(db)})

    assert detected["status"] == "ok"
    assert detected["valid"] is True
    assert detected["available"] is True
    assert detected["source"] == "BAKTA_DB"
    assert detected["version"] == "5.1.2"
    assert detected["version_source"] == "version.json"
    assert "version.json" in detected["markers_found"]
    marker = detected["marker_files"]["version.json"]
    assert marker["exists"] is True
    assert len(marker["sha256"]) == 64


def test_detect_bakta_db_explicit_argument_precedes_env(tmp_path: Path):
    explicit = make_bakta_db(tmp_path / "explicit", version="5.0")
    env = make_bakta_db(tmp_path / "env", version="4.0")
    detected = ext.detect_bakta_db(
        explicit, environ={"BAKTA_DB": str(env)}
    )
    assert detected["source"] == "argument"
    assert detected["version"] == "5.0"
    assert str(explicit) in (detected["path"] or "")


def test_detect_bakta_db_nested_layout(tmp_path: Path):
    db = make_bakta_db(tmp_path / "bakta_db", version="5.2", nested=True)
    detected = ext.detect_bakta_db(db, environ={})
    assert detected["status"] == "ok"
    assert detected["valid"] is True
    assert detected["version"] == "5.2"


def test_detect_bakta_db_missing_and_incomplete_are_not_valid(tmp_path: Path):
    missing = ext.detect_bakta_db(tmp_path / "does_not_exist", environ={})
    assert missing["status"] == "not_found"
    assert missing["valid"] is False
    assert missing["path"] is None
    assert "no la descarga" in missing["reason"] or "no se encontró" in missing["reason"]

    empty = tmp_path / "empty_db"
    empty.mkdir()
    incomplete = ext.detect_bakta_db(empty, environ={})
    assert incomplete["status"] == "incomplete"
    assert incomplete["valid"] is False
    assert "sin marcadores" in incomplete["reason"]

    file_path = tmp_path / "not_a_dir"
    file_path.write_text("x", encoding="utf-8")
    invalid = ext.detect_bakta_db(file_path, environ={})
    assert invalid["status"] == "invalid"
    assert invalid["valid"] is False


def test_bakta_run_without_database_does_not_execute_or_invent_metrics(
    tmp_path: Path, monkeypatch
):
    genome, reference = synthetic_inputs(tmp_path)
    fake = write_fake_tool(tmp_path)
    called: list[object] = []

    def forbidden(*args, **kwargs):
        called.append((args, kwargs))
        raise AssertionError("Bakta no debe ejecutarse sin DB válida")

    monkeypatch.setattr(ext.subprocess, "run", forbidden)
    data = ext.run_benchmark(
        genome,
        reference,
        tools="bakta",
        run=True,
        which=lambda _name: str(fake),
        bakta_db=tmp_path / "missing_db",
        environ={},
        probe_versions=False,
    )
    result = data["results"]["bakta"]
    assert result["status"] == "not_available"
    assert result["comparison"] is None
    assert result["execution"] is None
    assert result["bakta_db"]["valid"] is False
    assert "Bakta DB" in result["reason"]
    assert called == []
    assert data["not_available"] == ["bakta"]


def test_bakta_run_with_valid_db_executes_and_records_evidence(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    fake = write_fake_tool(tmp_path)
    db = make_bakta_db(tmp_path / "bakta_db")
    evidence = tmp_path / "evidence"
    report_json = tmp_path / "reports" / "benchmark.json"
    report_md = tmp_path / "reports" / "benchmark.md"

    data = ext.run_benchmark(
        genome,
        reference,
        tools="bakta",
        run=True,
        which=lambda _name: str(fake),
        bakta_db=db,
        environ={},
        probe_versions=True,
        save_gffs=evidence,
        command=["python", "external_benchmark.py", "--run"],
        report_json=report_json,
        report_md=report_md,
    )
    result = data["results"]["bakta"]

    assert result["status"] == "ok", result
    assert result["source"] == "executed"
    assert result["comparison"]["aggregate"]["tp"] == 2
    assert result["execution"]["status"] == "ok"
    assert result["execution"]["runtime_seconds"] >= 0
    assert result["execution"]["started_at"]
    assert result["execution"]["finished_at"]
    assert "--db" in result["execution"]["command"]
    command_db = result["execution"]["command"][
        result["execution"]["command"].index("--db") + 1
    ]
    assert command_db == str(db.resolve())
    assert len(result["execution"]["output"]["sha256"]) == 64
    assert result["execution"]["output"]["size_bytes"] > 0
    assert result["version"]["status"] == "ok"
    assert result["version"]["version"] == "1.12.2"
    assert result["timings"]["total_seconds"] >= 0
    assert result["query_gff"]["file"]["sha256"] == result["execution"]["output"]["sha256"]
    assert (evidence / "bakta" / "command.txt").is_file()
    assert (evidence / "bakta" / "result.gff3").is_file()

    assert data["environment"]["bakta_db"]["valid"] is True
    assert data["environment"]["bakta_db"]["version"] == "5.1"
    assert data["validation"]["timing"]["total_seconds"] >= 0
    assert data["validation"]["reports"]["json"] == str(report_json)
    assert data["validation"]["reports"]["md"] == str(report_md)
    assert data["validation"]["versions"]["bakta"] == "1.12.2"
    assert data["not_available"] == []


def test_bakta_failed_execution_has_no_metrics(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    fake = write_fake_tool(tmp_path, fail=True)
    db = make_bakta_db(tmp_path / "bakta_db")
    data = ext.run_benchmark(
        genome,
        reference,
        tools="bakta",
        run=True,
        which=lambda _name: str(fake),
        bakta_db=db,
        environ={},
        probe_versions=False,
    )
    result = data["results"]["bakta"]
    assert result["status"] == "failed"
    assert result["comparison"] is None
    assert result["execution"]["returncode"] == 3


# ---------------------------------------------------------------------------
# Gold pequeño y reproducible
# ---------------------------------------------------------------------------
def test_gold_from_local_gff_is_small_reindexed_and_has_sha256(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    out = tmp_path / "gold"
    manifest = make_gold.build_gold(
        genome_path=genome,
        annotation_path=reference,
        out_dir=out,
        name="synthetic_gold",
        max_bp=1000,
        flank=10,
        min_features=2,
    )

    gold_fasta = out / "genome.fna"
    gold_gff = out / "reference.gff3"
    manifest_path = out / "gold_manifest.json"
    run_script = out / "run_benchmark.sh"
    assert gold_fasta.is_file() and gold_gff.is_file()
    assert manifest_path.is_file() and run_script.is_file()
    assert os.access(run_script, os.X_OK)
    assert manifest["schema"] == ext.GOLD_SCHEMA
    assert manifest["name"] == "synthetic_gold"
    assert manifest["counts"]["CDS"] >= 2
    assert manifest["selection"]["window_bp"] <= 1000
    assert manifest["source"]["annotation"]["file"]["sha256"] == ext.sha256_file(reference)
    assert manifest["source"]["genome"]["file"]["sha256"] == ext.sha256_file(genome)
    assert len(manifest["outputs"]["genome"]["sha256"]) == 64
    assert len(manifest["outputs"]["reference"]["sha256"]) == 64

    gold_set = ext.load_annotations(gold_gff)
    assert gold_set.counts == manifest["counts"]
    assert all(feature.start >= 1 for feature in gold_set.features)
    assert all(
        feature.end <= manifest["selection"]["window_bp"]
        for feature in gold_set.features
    )
    # La ventana está reindexada: no conserva las coordenadas del contig original.
    original = ext.load_annotations(reference)
    first_original = min(feature.start for feature in original.features)
    assert min(feature.start for feature in gold_set.features) < first_original

    # El comando --run es reproducible y guarda informes bajo reports/external_benchmark/.
    command = manifest["run"]["command"]
    assert command[0]
    assert "--run" in command
    assert "--genome" in command and "--reference" in command
    assert manifest["run"]["reports"]["json"].startswith(
        "reports/external_benchmark/"
    )
    script_text = run_script.read_text(encoding="utf-8")
    assert "--run" in script_text
    assert "BAKTA_DB" in script_text
    check = subprocess.run(
        ["bash", "-n", str(run_script)], capture_output=True, text=True
    )
    assert check.returncode == 0, check.stderr


def test_gold_window_sequence_matches_source_slice(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    source = genome.read_text(encoding="utf-8").splitlines()[1]
    out = tmp_path / "gold"
    manifest = make_gold.build_gold(
        genome_path=genome,
        annotation_path=reference,
        out_dir=out,
        name="slice_gold",
        max_bp=1200,
        flank=50,
        min_features=2,
    )
    gold_sequence = "".join(
        out.joinpath("genome.fna").read_text(encoding="utf-8").splitlines()[1:]
    )
    start = manifest["selection"]["window_start"]
    end = manifest["selection"]["window_end"]
    assert gold_sequence == source[start - 1 : end]
    assert len(gold_sequence) == manifest["selection"]["window_bp"]


def test_gold_from_local_genbank_without_external_fasta(tmp_path: Path):
    gbff = write_genbank(tmp_path / "reference.gbff")
    out = tmp_path / "gold_gbff"
    manifest = make_gold.build_gold(
        annotation_path=gbff,
        out_dir=out,
        name="gbff_gold",
        max_bp=600,
        flank=0,
        min_features=2,
    )
    assert manifest["source"]["genome"]["path"] is None
    assert manifest["source"]["genome"]["sequence_source"]
    assert manifest["counts"]["CDS"] >= 1
    gold = ext.load_annotations(out / "reference.gff3")
    assert gold.counts == manifest["counts"]
    assert gold.counts["tRNA"] == 1


def test_gold_rejects_nonempty_output_without_force(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    out = tmp_path / "gold"
    out.mkdir()
    (out / "keep.txt").write_text("x", encoding="utf-8")
    with pytest.raises(ext.InputError):
        make_gold.build_gold(
            genome_path=genome,
            annotation_path=reference,
            out_dir=out,
            name="x",
        )


# ---------------------------------------------------------------------------
# Informes, validación y ausencia de red
# ---------------------------------------------------------------------------
def test_default_reports_live_under_reports_external_benchmark(tmp_path: Path, monkeypatch):
    genome, reference = synthetic_inputs(tmp_path)
    query = write_gff(
        tmp_path / "query.gff3",
        [
            ("ctg", "CDS", 101, 400, "+", "ID=q1;product=DNA polymerase"),
            ("ctg", "tRNA", 900, 970, "+", "ID=q2;product=tRNA-Ala"),
        ],
    )
    monkeypatch.chdir(tmp_path)
    code = ext.main(
        [
            "--genome",
            str(genome),
            "--reference",
            str(reference),
            "--tools",
            "prokka",
            "--query-gff",
            f"prokka={query}",
            "--quiet",
        ]
    )
    assert code == 0
    default_json = tmp_path / "reports" / "external_benchmark" / "benchmark.json"
    default_md = tmp_path / "reports" / "external_benchmark" / "benchmark.md"
    assert default_json.is_file() and default_md.is_file()
    payload = json.loads(default_json.read_text(encoding="utf-8"))
    assert payload["results"]["prokka"]["source"] == "precomputed"
    assert payload["validation"]["reports"]["json"].endswith(
        "reports/external_benchmark/benchmark.json"
    )
    assert "## Validación" in default_md.read_text(encoding="utf-8")


def test_missing_tools_report_exact_command_and_blocking_reason(tmp_path: Path):
    genome, reference = synthetic_inputs(tmp_path)
    data = ext.run_benchmark(
        genome,
        reference,
        tools="all",
        run=True,
        which=lambda _name: None,
        environ={},
        report_json=tmp_path / "benchmark.json",
        report_md=tmp_path / "benchmark.md",
    )
    command = data["reproduction"]["run_command"]
    assert "--run" in command
    assert "--genome" in command and "--reference" in command
    assert any("Bakta DB" in warning for warning in data["warnings"])
    checks = {item["name"]: item for item in data["validation"]["checks"]}
    assert checks["bakta_db"]["status"] == "missing"
    assert checks["tool_bakta"]["status"] == "missing"
    assert checks["result_bakta"]["status"] == "missing"
    markdown = ext.render_markdown(data)
    assert "not_available" in markdown
    assert "no disponible" in markdown
    assert "## Validación" in markdown
    assert "--run" in markdown


def test_runner_source_has_no_network_downloads():
    source = (MODULE_DIR / "external_benchmark.py").read_text(encoding="utf-8").lower()
    for token in (
        "requests.",
        "urlopen(",
        "urlretrieve(",
        "curl ",
        "wget ",
        "docker pull",
        "bakta_db download",
    ):
        assert token not in source, f"el runner no debe contener {token!r}"


def test_parse_version_text():
    assert ext.parse_version_text("bakta 1.12.2") == "1.12.2"
    assert ext.parse_version_text("Prokka version 1.15.6") == "1.15.6"
    assert ext.parse_version_text("dfast v1.6.2") == "v1.6.2"
    assert ext.parse_version_text("sin versión") is None


def test_committed_gold_sample_is_consistent_and_runnable():
    """El gold real versionado es local, pequeño y su manifiesto cuadra."""
    gold_dir = MODULE_DIR / "gold" / "chlorobaculum_tepidum_60kb"
    assert gold_dir.is_dir()
    manifest = json.loads(
        (gold_dir / "gold_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["schema"] == ext.GOLD_SCHEMA
    assert manifest["name"] == "chlorobaculum_tepidum_60kb"
    assert manifest["selection"]["window_bp"] == 60_000
    assert manifest["counts"]["CDS"] == 55
    assert manifest["counts"]["tRNA"] == 1

    genome = gold_dir / "genome.fna"
    reference = gold_dir / "reference.gff3"
    assert manifest["outputs"]["genome"]["sha256"] == ext.sha256_file(genome)
    assert manifest["outputs"]["reference"]["sha256"] == ext.sha256_file(reference)
    parsed = ext.load_annotations(reference)
    assert parsed.counts == manifest["counts"]
    assert parsed.format == "gff3"

    command = manifest["run"]["command"]
    assert "--run" in command and "--genome" in command and "--reference" in command
    assert manifest["run"]["script"].endswith("run_benchmark.sh")
    assert (gold_dir / "run_benchmark.sh").is_file()
    assert os.access(gold_dir / "run_benchmark.sh", os.X_OK)

    # El gold se puede comparar consigo mismo sin invocar ninguna herramienta.
    data = ext.run_benchmark(
        genome,
        reference,
        tools="prokka",
        run=False,
        which=lambda _name: None,
        environ={},
        query_gffs={"prokka": reference},
    )
    assert data["results"]["prokka"]["status"] == "ok"
    assert data["results"]["prokka"]["comparison"]["aggregate"]["f1"] == 1.0
    assert data["not_available"] == []


def test_cli_run_with_fake_bakta_on_path_and_default_reports(
    tmp_path: Path, monkeypatch
):
    """El camino real del CLI detecta PATH, valida BAKTA_DB y ejecuta ``--run``."""
    genome, reference = synthetic_inputs(tmp_path)
    bin_dir = tmp_path / "bin"
    write_fake_tool(bin_dir, binary_name="bakta", version="bakta 1.12.2")
    db = make_bakta_db(tmp_path / "bakta_db", version="5.1.2")
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("BAKTA_DB", str(db))
    monkeypatch.chdir(tmp_path)

    code = ext.main(
        [
            "--genome",
            str(genome),
            "--reference",
            str(reference),
            "--tools",
            "bakta",
            "--run",
            "--threads",
            "1",
            "--quiet",
        ]
    )
    assert code == 0

    report_dir = tmp_path / "reports" / "external_benchmark"
    payload = json.loads((report_dir / "benchmark.json").read_text(encoding="utf-8"))
    result = payload["results"]["bakta"]
    assert result["status"] == "ok", result
    assert result["source"] == "executed"
    assert result["comparison"]["aggregate"]["tp"] == 2
    assert payload["environment"]["bakta_db"]["source"] == "BAKTA_DB"
    assert payload["environment"]["bakta_db"]["version"] == "5.1.2"
    assert payload["validation"]["bakta_db"]["valid"] is True
    assert result["version"]["version"] == "1.12.2"
    command = result["execution"]["command"]
    assert "--db" in command
    assert command[command.index("--db") + 1] == str(db.resolve())
    assert (report_dir / "benchmark.md").is_file()
