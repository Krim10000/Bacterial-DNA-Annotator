"""Tests de ``tools/benchmark_external.py`` con GFF sintéticos y binarios ausentes.

No usan red, no descargan bases y no ejecutan bakta/prokka/dfast. Las pruebas de
ejecución usan el propio intérprete de Python en local.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import tools.benchmark_external as ext  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_gff(
    path: Path,
    features: list[tuple[str, str, int, int, str, str]],
    gz: bool = False,
    regions: dict[str, int] | None = None,
) -> Path:
    """Escribe un GFF3 mínimo: (seqid, type, start, end, strand, attrs)."""
    lines = ["##gff-version 3"]
    for contig, length in (regions or {}).items():
        lines.append(f"##sequence-region {contig} 1 {length}")
    for seqid, ftype, start, end, strand, attrs in features:
        lines.append(
            f"{seqid}\ttest\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
        )
    text = "\n".join(lines) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def synthetic_gffs(tmp_path: Path) -> tuple[Path, Path]:
    ref = write_gff(
        tmp_path / "reference.gff3",
        [
            ("ctg", "CDS", 100, 300, "+", "ID=cds1;product=DNA polymerase"),
            ("ctg", "CDS", 400, 600, "-", "ID=cds2;product=hypothetical protein"),
            ("ctg", "CDS", 1000, 1100, "+", "ID=cds3;pseudo=true;product=fragment"),
            ("ctg", "tRNA", 700, 770, "+", "ID=trna1;product=tRNA-Ala"),
            ("ctg", "rRNA", 800, 900, "+", "ID=rrna1;product=16S ribosomal RNA"),
        ],
        regions={"ctg": 2000},
    )
    query = write_gff(
        tmp_path / "query.gff3",
        [
            ("ctg", "CDS", 100, 300, "+", "ID=q1;product=DNA polymerase"),
            ("ctg", "CDS", 405, 600, "-", "ID=q2;product=hypothetical protein"),
            ("ctg", "CDS", 1200, 1300, "+", "ID=q3;product=unknown protein"),
            ("ctg", "tRNA", 702, 770, "+", "ID=q4;product=tRNA-Ala"),
        ],
        regions={"ctg": 2000},
    )
    return ref, query


# ---------------------------------------------------------------------------
# Parseo
# ---------------------------------------------------------------------------
def test_parse_gff_counts_excludes_pseudo_and_reads_gzip(tmp_path: Path):
    path = write_gff(
        tmp_path / "x.gff.gz",
        [
            ("ctg", "CDS", 10, 100, "+", "ID=c1;product=DNA polymerase"),
            ("ctg", "CDS", 200, 300, "+", "ID=c2;pseudo=true;product=fragment"),
            ("ctg", "tRNA", 400, 470, "-", "ID=t1;product=tRNA-Ala"),
        ],
        gz=True,
        regions={"ctg": 1000},
    )
    parsed = ext.parse_gff(path)
    assert parsed["counts"] == {"CDS": 1, "tRNA": 1, "rRNA": 0}
    assert parsed["excluded_pseudo"] == 1
    assert parsed["contigs"] == ["ctg"]
    assert parsed["sequence_regions"]["ctg"] == 1000
    assert parsed["features"][0]["attrs"]["product"] == "DNA polymerase"


def test_hypothetical_classifier_is_conservative():
    assert ext.is_hypothetical_product("hypothetical protein")
    assert ext.is_hypothetical_product("Uncharacterized protein")
    assert ext.is_hypothetical_product("protein of unknown function")
    assert ext.is_hypothetical_product("DUF1234 domain protein")
    assert not ext.is_hypothetical_product("DNA polymerase")
    assert not ext.is_hypothetical_product("putative transporter")
    assert not ext.is_hypothetical_product("")


# ---------------------------------------------------------------------------
# Comparación
# ---------------------------------------------------------------------------
def test_compare_gffs_synthetic_metrics(tmp_path: Path):
    ref, query = synthetic_gffs(tmp_path)
    comp = ext.compare_gffs(ref, query, threshold=0.5)
    assert comp["status"] == "ok"
    assert comp["reference_counts"] == {"CDS": 2, "tRNA": 1, "rRNA": 1}
    assert comp["query_counts"] == {"CDS": 3, "tRNA": 1, "rRNA": 0}

    cds = comp["categories"]["CDS"]
    assert (cds["tp"], cds["fp"], cds["fn"]) == (2, 1, 0)
    assert cds["precision"] == pytest.approx(0.6667, abs=1e-3)
    assert cds["recall"] == 1.0
    assert cds["f1"] == pytest.approx(0.8, abs=1e-3)
    assert cds["exact_matches"] == 1

    trna = comp["categories"]["tRNA"]
    assert (trna["tp"], trna["fp"], trna["fn"]) == (1, 0, 0)
    assert trna["exact_matches"] == 0

    rrna = comp["categories"]["rRNA"]
    assert (rrna["tp"], rrna["fp"], rrna["fn"]) == (0, 0, 1)
    assert rrna["precision"] is None
    assert rrna["recall"] == 0.0

    agg = comp["aggregate"]
    assert agg["reference_n"] == 4
    assert agg["query_n"] == 4
    assert (agg["tp"], agg["fp"], agg["fn"]) == (3, 1, 1)
    assert agg["precision"] == 0.75
    assert agg["recall"] == 0.75
    assert agg["f1"] == 0.75

    attrs = cds["attributes"]
    assert attrs["matched_pairs"] == 2
    assert attrs["product"]["literal_agreement"] == 1.0
    assert attrs["query_attribute_coverage"]["product"]["n"] == 3

    assert cds["hypothetical"]["reference"]["pct_hypothetical"] == 50.0
    assert cds["hypothetical"]["query"]["pct_hypothetical"] == pytest.approx(66.67, abs=0.01)
    assert cds["hypothetical"]["query"]["hypothetical"] == 2


def test_compare_missing_reference_is_not_evaluable(tmp_path: Path):
    _, query = synthetic_gffs(tmp_path)
    comp = ext.compare_gffs(tmp_path / "missing.gff3", query)
    assert comp["status"] == "not_evaluable"
    assert comp["categories"] == {}
    assert comp["not_evaluable"]


def test_compare_without_product_is_partial(tmp_path: Path):
    ref = write_gff(
        tmp_path / "ref.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=r1;product=DNA polymerase")],
    )
    query = write_gff(
        tmp_path / "q.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=q1")],
    )
    comp = ext.compare_gffs(ref, query)
    assert comp["status"] == "partial"
    assert any("product" in item["metric"] for item in comp["not_evaluable"])


# ---------------------------------------------------------------------------
# Detección y ejecución
# ---------------------------------------------------------------------------
def test_detect_tools_missing():
    detected = ext.detect_tools(which=lambda _name: None)
    assert set(detected) == {"bakta", "prokka", "dfast"}
    assert all(not item["available"] for item in detected.values())
    assert all("shutil.which" in item["reason"] for item in detected.values())


def test_benchmark_tool_missing_binary_does_not_execute(tmp_path: Path, monkeypatch):
    ref, _ = synthetic_gffs(tmp_path)
    called: list[object] = []

    def fake_run(*args, **kwargs):
        called.append((args, kwargs))
        raise AssertionError("no debe ejecutarse un binario ausente")

    monkeypatch.setattr(ext.subprocess, "run", fake_run)
    result = ext.benchmark_tool(
        "bakta",
        genome=tmp_path / "genome.fna",
        reference_gff=ref,
        run=True,
        which=lambda _name: None,
    )
    assert result["status"] == "not_evaluable"
    assert result["available"] is False
    assert result["run"] is None
    assert result["not_evaluable"]
    assert called == []


def test_benchmark_tool_available_without_run_does_not_execute(tmp_path: Path, monkeypatch):
    ref, _ = synthetic_gffs(tmp_path)
    called: list[object] = []
    monkeypatch.setattr(ext.subprocess, "run", lambda *a, **k: called.append((a, k)))
    result = ext.benchmark_tool(
        "prokka",
        genome=tmp_path / "genome.fna",
        reference_gff=ref,
        run=False,
        which=lambda _name: "/bin/true",
    )
    assert result["status"] == "not_run"
    assert result["run"]["status"] == "not_run"
    assert called == []


def test_execute_command_measures_time_and_rss():
    result = ext.execute_command(
        [sys.executable, "-c", "print('ok')"],
        timeout=30,
    )
    assert result["status"] == "ok"
    assert result["returncode"] == 0
    assert result["runtime_seconds"] >= 0
    assert "ok" in result["stdout_tail"]
    # ru_maxrss existe en Linux; en otras plataformas puede ser None.
    assert result["max_rss_kb"] is None or result["max_rss_kb"] >= 0


# ---------------------------------------------------------------------------
# Informe
# ---------------------------------------------------------------------------
def test_build_data_and_markdown_with_missing_tools(tmp_path: Path):
    ref, _ = synthetic_gffs(tmp_path)
    data = ext.build_data(
        genome=tmp_path / "genome.fna",
        reference_gff=ref,
        tool="all",
        run=True,
        which=lambda _name: None,
    )
    assert data["schema"] == "external-benchmark/v1"
    assert data["reference"]["status"] == "ok"
    assert data["environment"]["available_tools"] == []
    assert all(r["status"] == "not_evaluable" for r in data["results"].values())
    assert data["not_evaluable"]
    md = ext.render_markdown(data)
    assert "# Benchmark externo contra anotadores" in md
    assert "no evaluable" in md.lower()
    assert "mejor" not in md.lower()
    assert "bakta" in md.lower()




def test_benchmark_tool_run_with_fake_binary_and_comparison(tmp_path: Path):
    """Ejecuta un binario falso local para probar el camino --run y la comparación."""
    fake = tmp_path / "fake_bakta.py"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "out = None\n"
        "for i, arg in enumerate(args):\n"
        "    if arg in {'--output', '--outdir', '--out'} and i + 1 < len(args):\n"
        "        out = args[i + 1]\n"
        "if out is None:\n"
        "    sys.exit(2)\n"
        "outdir = pathlib.Path(out); outdir.mkdir(parents=True, exist_ok=True)\n"
        "(outdir / 'result.gff3').write_text(\n"
        "    '##gff-version 3\\n'\n"
        "    'ctg\\ttest\\tCDS\\t1\\t100\\t.\\t+\\t.\\tID=q1;product=DNA polymerase\\n'\n"
        ")\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    ref = write_gff(
        tmp_path / "ref.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=r1;product=DNA polymerase")],
    )
    genome = tmp_path / "genome.fna"
    genome.write_text(">ctg\nACGTACGT\n", encoding="utf-8")

    result = ext.benchmark_tool(
        "bakta",
        genome=genome,
        reference_gff=ref,
        run=True,
        which=lambda _name: str(fake),
        timeout=30,
    )
    assert result["status"] == "ok"
    assert result["run"]["status"] == "ok"
    assert result["run"]["runtime_seconds"] >= 0
    assert result["run"]["gff"].endswith("result.gff3")
    comp = result["comparison"]
    assert comp["aggregate"]["tp"] == 1

    data = {
        "schema": "external-benchmark/v1",
        "generated_at": "2026-01-01T00:00:00Z",
        "project_root": str(tmp_path),
        "scope": "prueba",
        "parameters": {"genome": "genome.fna", "reference_gff": "ref.gff3",
                       "tool": "bakta", "run": True, "overlap_threshold": 0.5},
        "environment": {"tools": {}, "available_tools": [], "missing_tools": []},
        "reference": ext.summarize_reference(ref),
        "results": {"bakta": result},
        "not_evaluable": [],
        "limitations": [],
    }
    md = ext.render_markdown(data)
    assert "#### Features" in md
    assert "#### Hipotéticas" in md
    assert "Bakta" in md


def test_main_writes_json_and_markdown(tmp_path: Path, monkeypatch):
    ref, _ = synthetic_gffs(tmp_path)
    monkeypatch.setattr(ext.shutil, "which", lambda _name: None)
    md = tmp_path / "out.md"
    js = tmp_path / "out.json"
    rc = ext.main([
        "--genome", str(tmp_path / "genome.fna"),
        "--reference-gff", str(ref),
        "--tool", "all",
        "--md", str(md),
        "--json", str(js),
        "--quiet",
    ])
    assert rc == 0
    assert md.is_file() and js.is_file()
    payload = json.loads(js.read_text(encoding="utf-8"))
    assert payload["schema"] == "external-benchmark/v1"
    assert payload["results"]["bakta"]["status"] == "not_evaluable"
    assert payload["command"][1].endswith("benchmark_external.py")
