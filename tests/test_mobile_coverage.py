"""Tests de ``tools/mobile_coverage.py`` con árboles sintéticos.

No leen los datos reales, no ejecutan HMMER/DefenseFinder/bases externas y no
escriben fuera de ``tmp_path``. Cubren el caso vacío, la ausencia de un genoma
del inventario, la diferencia entre "0 arrays evaluados" y "CRISPR no evaluado",
y un árbol sin inventario.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import tools.mobile_coverage as mc  # noqa: E402

TYPED_HEADER = (
    "genome\tcontig\tstart\tend\tstrand\telement_type\tsubtype\tgene\tevidence\tconfidence"
)
INVENTORY_HEADER = (
    "genome\tphylum\torganism\tcontigs\tgenome_bp\thas_fasta\tmissing_inputs"
)
PRIOR_HEADER = "genome\ttype\tstart\tend\tname\tsource"
SPACER_HEADER = "genome\tcontig\tstart\tend\tspacer"


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _typed(root: Path, genome: str, types: list[str]) -> Path:
    rows = []
    for i, etype in enumerate(types, start=1):
        rows.append(
            f"{genome}\tc1\t{i * 100}\t{i * 100 + 50}\t+\t{etype}\t"
            f"sub{i}\tg{i}\tsynthetic\tmedium"
        )
    return _write(
        root / "mobile_typing" / "full" / f"{genome}.tsv",
        TYPED_HEADER + "\n" + "\n".join(rows) + "\n",
    )


def _gff(root: Path, genome: str) -> Path:
    text = (
        "##gff-version 3\n"
        f"##sequence-region c1 1 100000\n"
        f"c1\ttest\tgene\t1\t50\t.\t+\t.\tID={genome}_g1\n"
    )
    return _write(root / "mobile_typing" / "gff" / "full" / f"{genome}.gff3", text)


def _spacers(root: Path, genome: str, rows: int = 1) -> Path:
    body = "".join(
        f"{genome}\tc1\t{100 + i}\t{130 + i}\tACGT\n" for i in range(rows)
    )
    return _write(
        root / "mobile_typing" / "crispr_spacers" / "full" / f"{genome}.tsv",
        SPACER_HEADER + "\n" + body,
    )


def _fasta(root: Path, genome: str) -> Path:
    return _write(root / "raw" / f"{genome}_genomic.fna.gz", "synthetic\n")


def _inventory(root: Path, genomes: list[str], with_fasta: set[str] | None = None) -> Path:
    with_fasta = with_fasta or set()
    lines = [INVENTORY_HEADER]
    for genome in genomes:
        has_fasta = "yes" if genome in with_fasta else "no"
        missing = "" if has_fasta == "yes" else "fasta"
        lines.append(f"{genome}\tFiloSintetico\tOrganismo {genome}\t1\t100000\t{has_fasta}\t{missing}")
    return _write(
        root / "mobile_typing" / "full_genomes.tsv", "\n".join(lines) + "\n"
    )


def _prior(root: Path, rows: list[tuple[str, str]]) -> Path:
    lines = [PRIOR_HEADER]
    for i, (genome, etype) in enumerate(rows, start=1):
        lines.append(f"{genome}\t{etype}\t{i * 10}\t{i * 10 + 5}\tgen{i}\tsynthetic")
    return _write(
        root / "mobile_defense" / "features.tsv", "\n".join(lines) + "\n"
    )


def _files(root: Path) -> set[str]:
    return {path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()}


def test_empty_data_root_is_not_an_error(tmp_path: Path):
    root = tmp_path / "data"
    root.mkdir()

    report = mc.build_report(root)

    assert report["status"] == "no_data"
    assert report["mode"] == "none"
    assert report["inventory"]["source"] is None
    assert report["inventory"]["expected_genomes"] is None
    assert report["coverage"]["genomes_discovered"] == 0
    assert report["coverage"]["genomes_with_typed_mobile"] == 0
    assert report["coverage"]["rows_total"] == 0
    assert report["per_genome"] == []
    assert report["gaps"]["typed_missing"] == []
    assert report["gaps"]["to_reach_inventory"] is None
    assert any("No hay full_genomes.tsv" in w for w in report["warnings"])

    json_path = tmp_path / "out" / "coverage.json"
    md_path = tmp_path / "out" / "coverage.md"
    mc.write_reports(report, json_path, md_path)

    assert json.loads(json_path.read_text(encoding="utf-8"))["status"] == "no_data"
    markdown = md_path.read_text(encoding="utf-8")
    assert "no_data" in markdown
    assert "No se encontraron datos" in markdown


def test_inventory_coverage_and_not_available_per_genome(tmp_path: Path):
    root = tmp_path / "data"
    _inventory(root, ["GCF_A.1", "GCF_B.1"], with_fasta={"GCF_A.1"})
    _typed(
        root,
        "GCF_A.1",
        [
            "transposon",
            "transposon",
            "prophage",
            "defense",
            "crispr",
            "cas",
            "amr",
            "virulence",
        ],
    )
    _gff(root, "GCF_A.1")
    _spacers(root, "GCF_A.1", rows=1)
    _fasta(root, "GCF_A.1")
    _prior(root, [("GCF_A.1", "amr"), ("GCF_A.1", "mobile_element")])
    before = _files(root)

    report = mc.build_report(root)

    # El informe no escribe nada en el árbol de datos.
    assert _files(root) == before

    assert report["status"] == "gaps"
    assert report["inventory"]["expected_genomes"] == 2
    coverage = report["coverage"]
    assert coverage["genomes_with_typed_mobile"] == 1
    assert coverage["genomes_without_typed_mobile"] == 1
    assert coverage["typed_pct"] == 50.0
    assert coverage["genomes_coverage_full"] == 1
    assert coverage["genomes_coverage_none"] == 1
    assert coverage["rows_total"] == 8
    assert coverage["typed_by_category"] == {
        "mobile_element": 2,
        "prophage": 1,
        "defense": 1,
        "crispr": 2,
        "amr": 1,
        "virulence": 1,
    }
    assert coverage["categories"]["crispr"]["typed_rows"] == 2
    assert coverage["categories"]["crispr"]["typed_genomes"] == 1
    assert coverage["categories"]["amr"]["prior_rows"] == 1

    by_id = {item["genome"]: item for item in report["per_genome"]}
    a = by_id["GCF_A.1"]
    assert a["status"] == "available"
    assert a["coverage"] == "full"
    assert a["typed_rows"] == 8
    assert a["not_available"] == {}
    assert a["layers"]["crispr_scan"]["status"] == "available"
    assert a["layers"]["crispr_scan"]["rows"] == 1

    b = by_id["GCF_B.1"]
    assert b["status"] == "not_available"
    assert b["coverage"] == "none"
    assert b["reason"] == "genoma del inventario sin salida de mobile_typing"
    assert set(b["not_available"]) == {
        "typed_mobile",
        "gff",
        "crispr_scan",
        "prior_features",
        "fasta",
    }
    assert "sin TSV tipado" in b["not_available"]["typed_mobile"]

    gaps = report["gaps"]
    assert gaps["scope"] == "inventory"
    assert gaps["typed_missing"] == ["GCF_B.1"]
    assert gaps["to_reach_inventory"] == 1
    assert gaps["fasta_missing"] == ["GCF_B.1"]
    assert gaps["crispr_not_evaluated"] == ["GCF_B.1"]
    assert any("sin TSV tipado" in w for w in report["warnings"])

    json_path = tmp_path / "reports" / "mobile_coverage.json"
    md_path = tmp_path / "reports" / "mobile_coverage.md"
    mc.write_reports(report, json_path, md_path)
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["gaps"]["typed_missing"] == ["GCF_B.1"]
    markdown = md_path.read_text(encoding="utf-8")
    assert "GCF_B.1" in markdown
    assert "`not_available`" in markdown


def test_zero_arrays_evaluated_is_not_the_same_as_not_evaluated(tmp_path: Path):
    root = tmp_path / "data"
    _inventory(root, ["GCF_C.1", "GCF_D.1"], with_fasta={"GCF_C.1", "GCF_D.1"})
    for genome in ("GCF_C.1", "GCF_D.1"):
        _typed(root, genome, ["transposon", "defense"])
        _gff(root, genome)
        _fasta(root, genome)
    _spacers(root, "GCF_C.1", rows=0)  # evaluado, sin arrays
    # GCF_D.1: sin archivo de espaciadores -> CRISPR no disponible
    _prior(root, [("GCF_C.1", "defense"), ("GCF_D.1", "defense")])

    report = mc.build_report(root)
    by_id = {item["genome"]: item for item in report["per_genome"]}

    c = by_id["GCF_C.1"]
    assert c["status"] == "available"
    assert c["coverage"] == "full"
    assert c["layers"]["crispr_scan"]["status"] == "available"
    assert c["layers"]["crispr_scan"]["rows"] == 0
    assert c["layers"]["crispr_scan"]["evaluated"] is True
    assert "crispr_scan" not in c["not_available"]

    d = by_id["GCF_D.1"]
    assert d["status"] == "available"
    assert d["coverage"] == "partial"
    assert d["layers"]["crispr_scan"]["status"] == "not_available"
    assert d["layers"]["crispr_scan"]["evaluated"] is False
    assert "falta el TSV de espaciadores" in d["not_available"]["crispr_scan"]

    # D tiene FASTA, así que no aparece como "FASTA ausente"; sí como espaciador.
    assert report["gaps"]["fasta_missing"] == []
    assert report["gaps"]["spacer_missing"] == ["GCF_D.1"]


def test_without_inventory_discovers_genomes_from_artifacts(tmp_path: Path):
    root = tmp_path / "data"
    _typed(root, "GCF_X.1", ["transposon"])
    _prior(root, [("GCF_Y.1", "amr")])
    annotation = root / "annotations_ab_initio" / "GCF_Z.1"
    _write(annotation / "annotation.gff3", "##gff-version 3\n")
    _write(annotation / "proteins.faa", ">p1\nMKV\n")

    report = mc.build_report(root)

    assert report["status"] == "no_inventory"
    assert report["inventory"]["expected_genomes"] is None
    assert report["coverage"]["genomes_discovered"] == 3
    assert report["coverage"]["genomes_with_typed_mobile"] == 1
    assert report["gaps"]["scope"] == "universe"
    assert report["gaps"]["to_reach_inventory"] is None
    assert set(report["gaps"]["typed_missing"]) == {"GCF_Y.1", "GCF_Z.1"}
    assert any("No hay full_genomes.tsv" in w for w in report["warnings"])

    by_id = {item["genome"]: item for item in report["per_genome"]}
    assert by_id["GCF_X.1"]["in_inventory"] is False
    assert by_id["GCF_X.1"]["inputs"]["annotation_gff"] is False
    assert by_id["GCF_X.1"]["coverage"] == "partial"
    assert by_id["GCF_Y.1"]["status"] == "not_available"
    assert by_id["GCF_Y.1"]["inputs"]["annotation_gff"] is False
    assert "solo capa previa" in by_id["GCF_Y.1"]["reason"]
    assert by_id["GCF_Z.1"]["inputs"]["annotation_gff"] is True
    assert by_id["GCF_Z.1"]["inputs"]["proteins"] is True


def test_cli_writes_both_reports(tmp_path: Path):
    root = tmp_path / "data"
    root.mkdir()
    json_path = tmp_path / "reports" / "mobile_coverage.json"
    md_path = tmp_path / "reports" / "mobile_coverage.md"

    rc = mc.main(
        [
            "--data",
            str(root),
            "--json",
            str(json_path),
            "--md",
            str(md_path),
            "--quiet",
        ]
    )

    assert rc == 0
    assert json_path.is_file()
    assert md_path.is_file()
    assert json.loads(json_path.read_text(encoding="utf-8"))["status"] == "no_data"
    assert "Cobertura móvil" in md_path.read_text(encoding="utf-8")


def test_prior_layer_is_auxiliary_not_required_for_full_coverage(tmp_path: Path):
    root = tmp_path / "data"
    _inventory(root, ["GCF_E.1"], with_fasta={"GCF_E.1"})
    _typed(root, "GCF_E.1", ["transposon"])
    _gff(root, "GCF_E.1")
    _spacers(root, "GCF_E.1", rows=0)
    _fasta(root, "GCF_E.1")
    # Sin data/mobile_defense/features.tsv a propósito.

    report = mc.build_report(root)
    genome = report["per_genome"][0]

    assert genome["status"] == "available"
    assert genome["coverage"] == "full"
    assert genome["not_available"] == {
        "prior_features": "sin filas para este genoma en data/mobile_defense/features.tsv"
    }
    assert report["gaps"]["prior_missing"] == ["GCF_E.1"]
