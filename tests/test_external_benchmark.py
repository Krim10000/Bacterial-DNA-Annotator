"""Tests del benchmark externo CPU (Bakta, Prokka, DFAST).

Todo es local y sintético: no se descargan imágenes, bases ni datos, y no se
ejecutan Bakta/Prokka/DFAST. Las pruebas de ejecución usan un binario falso de
Python escrito en ``tmp_path``.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_DIR = ROOT / "benchmarks" / "external"
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import external_benchmark as ext  # noqa: E402

DOCKERFILE = MODULE_DIR / "Dockerfile"
COMPOSE = MODULE_DIR / "docker-compose.yml"
PREPARE = MODULE_DIR / "prepare_tools.sh"
DOC = ROOT / "docs" / "EXTERNAL_BENCHMARK_RUN.md"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_gff(
    path: Path,
    features: list[tuple[str, str, int, int, str, str]],
    *,
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


def write_genome(path: Path, contig: str = "ctg", length: int = 2000) -> Path:
    path.write_text(f">{contig} synthetic contig\n" + "ACGT" * (length // 4) + "\n",
                    encoding="utf-8")
    return path


def write_genbank(path: Path, *, locus: str = "ctg", length: int = 60) -> Path:
    """GenBank mínimo con 1 CDS, 1 CDS pseudo y 1 tRNA."""
    sequence = ("ACGT" * ((length // 4) + 1))[:length]
    text = f"""LOCUS       {locus:<16}{length:>12} bp    DNA     linear   BCT 01-JAN-2026
DEFINITION  synthetic record.
ACCESSION   {locus}
VERSION     {locus}.1
KEYWORDS    .
SOURCE      synthetic
  ORGANISM  synthetic bacterium
            Bacteria.
FEATURES             Location/Qualifiers
     source          1..{length}
                     /organism="synthetic bacterium"
     CDS             1..30
                     /product="DNA polymerase"
                     /gene="polA"
                     /locus_tag="CTAG_0001"
     CDS             complement(31..60)
                     /pseudo
                     /product="hypothetical protein"
     tRNA            40..55
                     /product="tRNA-Ala"
ORIGIN
        1 {sequence}
//
"""
    path.write_text(text, encoding="utf-8")
    return path


def synthetic_pair(tmp_path: Path) -> tuple[Path, Path]:
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


def fake_tool(tmp_path: Path, *, fail: bool = False, version: str = "fake 1.2.3") -> Path:
    """Binario falso que escribe un GFF3 y responde a --version."""
    body = (
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        f"if '--version' in args:\n"
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
            "(outdir / 'result.gff3').write_text(\n"
            "    '##gff-version 3\\n'\n"
            "    'ctg\\ttest\\tCDS\\t100\\t300\\t.\\t+\\t.\\tID=q1;product=DNA polymerase\\n'\n"
            "    'ctg\\ttest\\ttRNA\\t700\\t770\\t.\\t+\\t.\\tID=q2;product=tRNA-Ala\\n'\n"
            ")\n"
            "sys.exit(0)\n"
        )
    path = tmp_path / "fake_tool.py"
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


# ---------------------------------------------------------------------------
# Parseo GFF3 / GenBank / FASTA
# ---------------------------------------------------------------------------
def test_parse_gff_excludes_pseudo_and_reads_gzip(tmp_path: Path):
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
    assert parsed.counts == {"CDS": 1, "tRNA": 1, "rRNA": 0}
    assert parsed.excluded_pseudo == 1
    assert parsed.contigs == ["ctg"]
    assert parsed.sequence_regions["ctg"] == 1000
    assert parsed.features[0].attr("product") == "DNA polymerase"


def test_parse_gff_merges_multisegment_cds(tmp_path: Path):
    path = write_gff(
        tmp_path / "frameshift.gff3",
        [
            ("ctg", "CDS", 1984037, 1984105, "+", "ID=cds1;product=release factor"),
            ("ctg", "CDS", 1984107, 1985132, "+", "ID=cds1;product=release factor"),
            ("ctg", "CDS", 1, 100, "+", "ID=cds2;product=other"),
        ],
    )
    parsed = ext.parse_gff(path)
    assert parsed.counts["CDS"] == 2
    merged = next(feature for feature in parsed.features if feature.attr("ID") == "cds1")
    assert (merged.start, merged.end) == (1984037, 1985132)
    assert any("multi-segmento" in warning for warning in parsed.warnings)


def test_parse_gff_keeps_pseudo_when_requested(tmp_path: Path):
    path = write_gff(
        tmp_path / "x.gff3",
        [("ctg", "CDS", 1, 90, "+", "ID=c1;pseudogene=true;product=fragment")],
    )
    parsed = ext.parse_gff(path, exclude_pseudo=False)
    assert parsed.counts["CDS"] == 1
    assert parsed.excluded_pseudo == 0


def test_parse_genbank_and_hypothetical_classifier(tmp_path: Path):
    path = write_genbank(tmp_path / "x.gbff")
    parsed = ext.parse_genbank(path)
    assert parsed.format == "genbank"
    assert parsed.counts == {"CDS": 1, "tRNA": 1, "rRNA": 0}
    assert parsed.excluded_pseudo == 1
    assert list(parsed.sequence_regions.values()) == [60]
    cds = next(feature for feature in parsed.features if feature.type == "CDS")
    assert (cds.start, cds.end, cds.strand) == (1, 30, "+")
    assert cds.attr("product") == "DNA polymerase"
    assert cds.attr("gene") == "polA"
    assert cds.attr("locus_tag") == "CTAG_0001"

    assert ext.is_hypothetical("hypothetical protein")
    assert ext.is_hypothetical("Uncharacterized protein")
    assert ext.is_hypothetical("protein of unknown function")
    assert ext.is_hypothetical("DUF1234 domain protein")
    assert not ext.is_hypothetical("putative transporter")
    assert not ext.is_hypothetical("DNA polymerase")
    assert not ext.is_hypothetical("")


def test_detect_format_auto(tmp_path: Path):
    gff = write_gff(tmp_path / "a.gff3", [("ctg", "CDS", 1, 10, "+", "ID=c")])
    gbk = write_genbank(tmp_path / "b.gbff")
    assert ext.detect_format(gff) == "gff3"
    assert ext.detect_format(gbk) == "genbank"


def test_read_fasta_and_seqid_resolver(tmp_path: Path):
    fasta = tmp_path / "g.fna"
    fasta.write_text(
        ">NC_000913.3 Escherichia coli str. K-12\nACGTACGT\n"
        ">plasmid_1\nACGT\n",
        encoding="utf-8",
    )
    genome = ext.read_fasta(fasta)
    assert [item["id"] for item in genome.contigs] == ["NC_000913.3", "plasmid_1"]
    assert genome.total_length == 12
    assert len(genome.sha256 or "") == 64

    resolver = ext.build_seqid_resolver(genome)
    assert resolver.resolve("NC_000913.3") == "NC_000913.3"
    assert resolver.resolve("NC_000913") == "NC_000913.3"
    assert resolver.resolve("ref|NC_000913.3|") == "NC_000913.3"
    assert resolver.resolve("NC_000913.3 Escherichia coli") == "NC_000913.3"
    assert resolver.resolve("otro_contig") == "otro_contig"


# ---------------------------------------------------------------------------
# Comparación y métricas
# ---------------------------------------------------------------------------
def test_compare_synthetic_metrics_and_differences(tmp_path: Path):
    ref, query = synthetic_pair(tmp_path)
    comp = ext.compare_gffs(ref, query)

    assert comp["status"] == "ok"
    assert comp["reference_counts"] == {"CDS": 2, "tRNA": 1, "rRNA": 1}
    assert comp["query_counts"] == {"CDS": 3, "tRNA": 1, "rRNA": 0}

    cds = comp["categories"]["CDS"]
    assert (cds["reference_n"], cds["query_n"]) == (2, 3)
    assert (cds["tp"], cds["fp"], cds["fn"]) == (2, 1, 0)
    assert cds["precision"] == pytest.approx(0.6667, abs=1e-3)
    assert cds["recall"] == 1.0
    assert cds["f1"] == pytest.approx(0.8, abs=1e-3)
    assert cds["exact_matches"] == 1
    assert cds["boundary_shifted"] == 1

    trna = comp["categories"]["tRNA"]
    assert (trna["tp"], trna["fp"], trna["fn"]) == (1, 0, 0)
    assert trna["exact_matches"] == 0

    rrna = comp["categories"]["rRNA"]
    assert (rrna["tp"], rrna["fp"], rrna["fn"]) == (0, 0, 1)
    assert rrna["precision"] is None
    assert rrna["recall"] == 0.0
    assert rrna["f1"] == 0.0

    aggregate = comp["aggregate"]
    assert aggregate["reference_n"] == 4
    assert aggregate["query_n"] == 4
    assert (aggregate["tp"], aggregate["fp"], aggregate["fn"]) == (3, 1, 1)
    assert aggregate["precision"] == 0.75
    assert aggregate["recall"] == 0.75
    assert aggregate["f1"] == 0.75

    differences = comp["differences"]
    assert differences["counts"]["false_negatives"] == 1
    assert differences["counts"]["false_positives"] == 1
    assert differences["counts"]["boundary_shifts"] == 2
    assert differences["false_negatives"][0]["type"] == "rRNA"
    assert differences["false_positives"][0]["product"] == "unknown protein"
    assert {item["start_delta"] for item in differences["boundary_shifts"]} == {5, 2}


def test_compare_attributes_and_hypothetical(tmp_path: Path):
    ref, query = synthetic_pair(tmp_path)
    comp = ext.compare_gffs(ref, query)
    attributes = comp["categories"]["CDS"]["attributes"]
    assert attributes["matched_pairs"] == 2
    assert attributes["product"]["both_present"] == 2
    assert attributes["product"]["literal_matches"] == 2
    assert attributes["product"]["literal_agreement"] == 1.0
    assert attributes["product"]["token_jaccard_mean"] == 1.0

    hypothetical = comp["categories"]["CDS"]["hypothetical"]
    assert hypothetical["reference"]["pct_hypothetical"] == 50.0
    assert hypothetical["query"]["pct_hypothetical"] == pytest.approx(66.67, abs=0.01)
    assert hypothetical["query"]["hypothetical"] == 2


def test_compare_identical_is_perfect(tmp_path: Path):
    ref, _ = synthetic_pair(tmp_path)
    comp = ext.compare_gffs(ref, ref)
    assert comp["aggregate"]["precision"] == 1.0
    assert comp["aggregate"]["recall"] == 1.0
    assert comp["aggregate"]["f1"] == 1.0
    assert comp["differences"]["counts"]["false_negatives"] == 0
    assert comp["differences"]["counts"]["false_positives"] == 0


def test_compare_strand_mismatch_is_diagnosed(tmp_path: Path):
    ref = write_gff(
        tmp_path / "ref.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=r1;product=DNA polymerase")],
    )
    query = write_gff(
        tmp_path / "q.gff3",
        [("ctg", "CDS", 1, 100, "-", "ID=q1;product=DNA polymerase")],
    )
    comp = ext.compare_gffs(ref, query)
    assert comp["aggregate"]["tp"] == 0
    assert comp["differences"]["counts"]["strand_mismatches"] == 1
    assert comp["differences"]["strand_mismatches"][0]["reference"]["strand"] == "+"
    assert comp["differences"]["strand_mismatches"][0]["query"]["strand"] == "-"


def test_compare_no_common_contigs_is_not_comparable(tmp_path: Path):
    ref = write_gff(
        tmp_path / "ref.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=r1;product=DNA polymerase")],
    )
    query = write_gff(
        tmp_path / "q.gff3",
        [("otro", "CDS", 1, 100, "+", "ID=q1;product=DNA polymerase")],
    )
    comp = ext.compare_gffs(ref, query)
    assert comp["status"] == "not_comparable"
    assert comp["categories"] == {}
    assert comp["aggregate"] == {}
    assert "contigs en común" in comp["reason"]


def test_compare_empty_query_is_not_comparable(tmp_path: Path):
    ref = write_gff(
        tmp_path / "ref.gff3",
        [("ctg", "CDS", 1, 100, "+", "ID=r1;product=DNA polymerase")],
    )
    query = write_gff(tmp_path / "q.gff3", [])
    comp = ext.compare_gffs(ref, query)
    assert comp["status"] == "not_comparable"
    assert "no tiene features" in comp["reason"]


def test_seqid_aliasing_with_genbank_reference(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    gbk = write_genbank(tmp_path / "ref.gbff")
    query = write_gff(
        tmp_path / "q.gff3",
        [
            ("ctg", "CDS", 1, 30, "+", "ID=q1;product=DNA polymerase"),
            ("ctg", "tRNA", 40, 55, "+", "ID=q2;product=tRNA-Ala"),
        ],
    )
    genome = ext.read_fasta(fasta)
    resolver = ext.build_seqid_resolver(genome)
    reference = ext.load_annotations(gbk)
    query_set = ext.load_annotations(query)
    comp = ext.compare_annotations(reference, query_set, resolver=resolver)
    assert comp["status"] == "ok"
    assert comp["aggregate"]["tp"] == 2
    assert comp["aggregate"]["fn"] == 0


def test_reciprocal_overlap_threshold():
    a = ext.Feature("ctg", "CDS", 1, 100, "+")
    b = ext.Feature("ctg", "CDS", 51, 150, "+")
    assert ext.interval_overlap(1, 100, 51, 150) == 50
    assert ext.reciprocal_overlap(a, b) == pytest.approx(0.5)
    matches, unmatched_ref, unmatched_query = ext.match_features([a], [b], threshold=0.5)
    assert len(matches) == 1
    assert unmatched_ref == [] and unmatched_query == []


# ---------------------------------------------------------------------------
# Detección, estados y ejecución
# ---------------------------------------------------------------------------
def test_missing_tools_are_not_available_and_do_not_execute(
    tmp_path: Path, monkeypatch
):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    called: list[object] = []

    def forbidden(*args, **kwargs):
        called.append((args, kwargs))
        raise AssertionError("no se debe ejecutar un binario ausente")

    monkeypatch.setattr(ext.subprocess, "run", forbidden)
    data = ext.run_benchmark(
        fasta, ref, tools="all", run=True, which=lambda _name: None
    )

    assert data["schema"] == ext.SCHEMA
    assert data["not_available"] == ["bakta", "prokka", "dfast"]
    for name in ext.TOOLS:
        result = data["results"][name]
        assert result["status"] == "not_available"
        assert result["comparison"] is None
        assert result["execution"] is None
        assert "shutil.which" in result["reason"]
    assert called == []


def test_available_but_not_run_is_not_run(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    data = ext.run_benchmark(
        fasta, ref, tools="prokka", run=False, which=lambda _name: "/bin/true"
    )
    result = data["results"]["prokka"]
    assert result["status"] == "not_run"
    assert result["comparison"] is None
    assert data["not_available"] == []


def test_precomputed_query_gff_compares_without_binary(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, query = synthetic_pair(tmp_path)
    data = ext.run_benchmark(
        fasta,
        ref,
        tools="prokka",
        run=False,
        which=lambda _name: None,
        query_gffs={"prokka": query},
    )
    result = data["results"]["prokka"]
    assert result["status"] == "ok"
    assert result["source"] == "precomputed"
    assert result["available"] is False
    assert result["comparison"]["aggregate"]["tp"] == 3
    assert result["query_gff"]["sha256"] == ext.sha256_file(query)
    assert data["not_available"] == []


def test_execute_fake_tool_and_compare(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    fake = fake_tool(tmp_path)
    data = ext.run_benchmark(
        fasta,
        ref,
        tools="prokka",
        run=True,
        which=lambda _name: str(fake),
        probe_versions=False,
    )
    result = data["results"]["prokka"]
    assert result["status"] == "ok", result
    assert result["source"] == "executed"
    assert result["execution"]["status"] == "ok"
    assert result["execution"]["runtime_seconds"] >= 0
    assert result["comparison"]["aggregate"]["tp"] == 2
    assert result["query_gff"]["counts"] == {"CDS": 1, "tRNA": 1, "rRNA": 0}


def test_version_probe_is_recorded(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    fake = fake_tool(tmp_path)
    data = ext.run_benchmark(
        fasta,
        ref,
        tools="prokka",
        run=True,
        which=lambda _name: str(fake),
        probe_versions=True,
    )
    result = data["results"]["prokka"]
    assert result["version"]["status"] == "ok"
    assert "fake 1.2.3" in result["version"]["text"]
    assert result["status"] == "ok"


def test_failed_tool_has_no_metrics(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    fake = fake_tool(tmp_path, fail=True)
    data = ext.run_benchmark(
        fasta,
        ref,
        tools="bakta",
        run=True,
        which=lambda _name: str(fake),
        probe_versions=False,
    )
    result = data["results"]["bakta"]
    assert result["status"] == "failed"
    assert result["comparison"] is None
    assert result["execution"]["returncode"] == 3


def test_build_command_is_cpu_single_thread(tmp_path: Path):
    genome = write_genome(tmp_path / "g.fna")
    out = tmp_path / "out"
    bakta = ext.build_command("bakta", "/opt/bakta", genome, out, bakta_db=tmp_path / "db")
    assert "--threads" in bakta and bakta[bakta.index("--threads") + 1] == "1"
    assert "--db" in bakta and "--skip-plot" in bakta
    assert bakta[-1] == str(genome.resolve())

    prokka = ext.build_command("prokka", "/opt/prokka", genome, out)
    assert prokka[prokka.index("--cpus") + 1] == "1"
    assert "--force" in prokka

    dfast = ext.build_command("dfast", "/opt/dfast", genome, out)
    assert dfast[dfast.index("--cpu") + 1] == "1"
    assert "--genome" in dfast


def test_find_output_annotation_prefers_gff(tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "annotation.gbk").write_text("LOCUS x\n", encoding="utf-8")
    (outdir / "annotation.gff").write_text("##gff-version 3\n", encoding="utf-8")
    chosen, candidates = ext.find_output_annotation(outdir)
    assert chosen is not None and chosen.name == "annotation.gff"
    assert set(candidates) == {"annotation.gbk", "annotation.gff"}


# ---------------------------------------------------------------------------
# CLI, JSON, Markdown
# ---------------------------------------------------------------------------
def test_cli_writes_json_and_markdown(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, query = synthetic_pair(tmp_path)
    json_path = tmp_path / "out" / "benchmark.json"
    md_path = tmp_path / "out" / "benchmark.md"
    code = ext.main(
        [
            "--genome",
            str(fasta),
            "--reference",
            str(ref),
            "--tools",
            "prokka",
            "--query-gff",
            f"prokka={query}",
            "--json",
            str(json_path),
            "--md",
            str(md_path),
            "--quiet",
        ]
    )
    assert code == 0
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["schema"] == ext.SCHEMA
    assert payload["results"]["prokka"]["status"] == "ok"
    assert payload["parameters"]["tools"] == ["prokka"]
    markdown = md_path.read_text(encoding="utf-8")
    assert "Benchmark externo CPU" in markdown
    assert "F1" in markdown
    assert "Limitaciones" in markdown


def test_cli_missing_reference_returns_2(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    code = ext.main(
        [
            "--genome",
            str(fasta),
            "--reference",
            str(tmp_path / "missing.gff3"),
            "--tools",
            "bakta",
            "--quiet",
        ]
    )
    assert code == 2


def test_markdown_marks_not_available_without_metrics(tmp_path: Path):
    fasta = write_genome(tmp_path / "g.fna")
    ref, _ = synthetic_pair(tmp_path)
    data = ext.run_benchmark(fasta, ref, tools="all", run=False, which=lambda _n: None)
    markdown = ext.render_markdown(data)
    assert "not_available" in markdown
    assert "No se ejecutó" in markdown or "no se ejecutó" in markdown
    assert "F1 global" in markdown
    assert "no declara superioridad" in markdown


# ---------------------------------------------------------------------------
# Receta CPU y documentación
# ---------------------------------------------------------------------------
def test_dockerfile_exists_and_is_cpu_only():
    assert DOCKERFILE.is_file()
    text = DOCKERFILE.read_text(encoding="utf-8")
    lower = text.lower()
    assert text.startswith("# syntax=docker/dockerfile")
    assert "FROM mambaorg/micromamba:" in text
    for forbidden in ("nvidia", "cuda", "gpu"):
        assert forbidden not in lower
    assert "micromamba create -y -n bakta" in text
    assert "micromamba create -y -n prokka" in text
    assert "micromamba create -y -n dfast" in text
    assert "COPY benchmarks/external/external_benchmark.py" in text
    assert "ENTRYPOINT" not in text  # permite usar python /opt/... y --help
    for forbidden in ("wget ", "curl ", "git clone", "bakta_db download"):
        assert forbidden not in lower


def test_dockerfile_pins_versions_and_databases_are_not_included():
    text = DOCKERFILE.read_text(encoding="utf-8")
    assert "ARG BAKTA_VERSION=1.12.2" in text
    assert "ARG PROKKA_VERSION=1.15.6" in text
    assert "ARG DFAST_VERSION=1.6.2" in text
    assert "VOLUME [\"/data\"]" in text
    assert "BAKTA_DB" in text
    assert "biopython" in text


def test_compose_and_prepare_script():
    assert COMPOSE.is_file()
    compose = COMPOSE.read_text(encoding="utf-8")
    assert "benchmarks/external/Dockerfile" in compose
    assert "context: ../.." in compose
    assert "BAKTA_DB" in compose
    assert "nvidia" not in compose.lower()

    assert PREPARE.is_file()
    prepare = PREPARE.read_text(encoding="utf-8")
    for tool in ("bakta", "prokka", "dfast"):
        assert tool in prepare
    assert "micromamba" in prepare
    assert "BAKTA_DB" not in prepare or "no" in prepare.lower()
    result = __import__("subprocess").run(
        ["bash", "-n", str(PREPARE)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_docs_exist_and_explain_usage_and_limits():
    assert DOC.is_file()
    text = DOC.read_text(encoding="utf-8")
    for token in (
        "not_available",
        "docker build",
        "docker compose",
        "--genome",
        "--reference",
        "--tools",
        "--run",
        "--query-gff",
        "--bakta-db",
        "SHA256",
        "F1",
        "precisión",
        "recall",
        "Bakta DB",
        "CPU",
        "Limitaciones",
    ):
        assert token in text, f"falta {token!r} en {DOC}"
    lower = text.lower()
    assert "no se construyó ni se ejecutó" in lower
    assert "no se descargan" in lower or "no descarga" in lower
    assert "esta receta" in lower or "receta" in lower
    # No debe afirmar resultados de herramientas que no se ejecutaron.
    for phrase in ("se ejecutó bakta", "se ejecutó prokka", "se ejecutó dfast"):
        assert phrase not in lower
