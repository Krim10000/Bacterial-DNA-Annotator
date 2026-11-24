"""Tests de ``tools/annotation_parity.py`` con GFF3 sintéticos.

No usan red, GPU ni binarios externos. Cubren paridad exacta, genes
ausentes/extra, fragmentación, fusión, hebra, productos, ncRNA, ORFs pequeños,
descubrimiento de entradas, gzip y el caso sin referencia (``not_evaluable``).
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import tools.annotation_parity as ap  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_gff(
    path: Path,
    rows: list[tuple[str, str, int, int, str, str]],
    *,
    gz: bool = False,
    regions: dict[str, int] | None = None,
) -> Path:
    """Escribe un GFF3 mínimo: (seqid, type, start, end, strand, attrs)."""
    lines = ["##gff-version 3"]
    for contig, length in (regions or {}).items():
        lines.append(f"##sequence-region {contig} 1 {length}")
    for seqid, ftype, start, end, strand, attrs in rows:
        phase = "0" if ftype == "CDS" else "."
        lines.append(
            "\t".join(
                [seqid, "test", ftype, str(start), str(end), ".", strand, phase, attrs]
            )
        )
    text = "\n".join(lines) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def make_args(
    tmp_path: Path,
    ours: Path,
    reference: Path | None = None,
    *extra: str,
):
    argv = [
        "--ours",
        str(ours),
        "--root",
        str(tmp_path),
        "--md",
        str(tmp_path / "reports" / "parity.md"),
        "--json",
        str(tmp_path / "reports" / "parity.json"),
        *extra,
    ]
    if reference is not None:
        argv.extend(["--reference", str(reference)])
    args = ap.build_parser().parse_args(argv)
    args._argv = argv  # type: ignore[attr-defined]
    return args


def run_pair(
    tmp_path: Path,
    ref_rows: list[tuple[str, str, int, int, str, str]],
    ours_rows: list[tuple[str, str, int, int, str, str]],
    *extra: str,
    gz_reference: bool = False,
) -> dict:
    ref = write_gff(tmp_path / "reference.gff3", ref_rows, gz=gz_reference)
    ours = write_gff(tmp_path / "ours.gff3", ours_rows)
    return ap.run_parity(make_args(tmp_path, ours, ref, *extra))


def categories(data: dict) -> dict[str, int]:
    return data["differences"]["actionable"]["by_category"]


# ---------------------------------------------------------------------------
# Paridad exacta y métricas básicas
# ---------------------------------------------------------------------------
def test_exact_parity(tmp_path: Path) -> None:
    rows = [
        ("ctg", "CDS", 100, 399, "+", "ID=cds1;product=DNA polymerase"),
        ("ctg", "CDS", 500, 799, "-", "ID=cds2;product=hypothetical protein"),
        ("ctg", "tRNA", 900, 970, "+", "ID=trna1;product=tRNA-Ala"),
        ("ctg", "rRNA", 1000, 1100, "+", "ID=rrna1;product=16S ribosomal RNA"),
    ]
    data = run_pair(tmp_path, rows, rows)

    assert data["status"] == "evaluated"
    assert data["summary"]["CDS"]["matched"] == 2
    assert data["summary"]["CDS"]["exact"] == 2
    assert data["summary"]["CDS"]["recall"] == 1.0
    assert data["summary"]["CDS"]["precision"] == 1.0
    assert data["summary"]["CDS"]["exact_rate"] == 1.0
    assert data["summary"]["tRNA"]["matched"] == 1
    assert data["summary"]["rRNA"]["matched"] == 1
    assert data["differences"]["actionable"]["total"] == 0
    assert data["differences"]["review"]["total"] == 0
    assert data["products"]["ours"]["pct_hypothetical"] == 50.0
    assert data["products"]["paired"]["literal_agreement"] == 1.0
    assert data["products"]["paired"]["hypothetical_disagreements"] == 0


def test_coordinate_tolerance(tmp_path: Path) -> None:
    ref = [("ctg", "CDS", 100, 400, "+", "ID=r1;product=DNA polymerase")]
    ours = [("ctg", "CDS", 102, 400, "+", "ID=o1;product=DNA polymerase")]
    data = run_pair(tmp_path, ref, ours, "--coord-tolerance", "2")

    assert data["summary"]["CDS"]["exact"] == 0
    assert data["summary"]["CDS"]["exact_rate"] == 0.0
    assert data["summary"]["CDS"]["within_tolerance"] == 1
    assert data["summary"]["CDS"]["within_tolerance_rate"] == 1.0
    # 2 nt no supera la alerta de borde por defecto (30 nt).
    assert "boundary_shift" not in categories(data)


# ---------------------------------------------------------------------------
# Ausentes/extra y ORFs pequeños
# ---------------------------------------------------------------------------
def test_missing_extra_and_small_orfs(tmp_path: Path) -> None:
    ref = [
        ("ctg", "CDS", 100, 399, "+", "ID=r1;product=hypothetical protein"),
        ("ctg", "CDS", 1000, 1899, "+", "ID=r2;product=DNA polymerase"),
    ]
    ours = [
        ("ctg", "CDS", 500, 799, "+", "ID=o1;product=hypothetical protein"),
        ("ctg", "CDS", 1000, 1899, "+", "ID=o2;product=DNA polymerase"),
    ]
    data = run_pair(tmp_path, ref, ours)
    cats = categories(data)

    assert cats["missing_cds"] == 1
    assert cats["extra_cds"] == 1
    assert data["summary"]["CDS"]["matched"] == 1
    assert data["summary"]["CDS"]["recall"] == 0.5
    assert data["summary"]["CDS"]["precision"] == 0.5
    assert data["small_orfs"]["ours"]["small"] == 1
    assert data["small_orfs"]["reference"]["small"] == 1
    assert data["small_orfs"]["missing_small"] == 1
    assert data["small_orfs"]["extra_small"] == 1

    missing = next(
        i for i in data["differences"]["actionable"]["items"] if i["category"] == "missing_cds"
    )
    extra = next(
        i for i in data["differences"]["actionable"]["items"] if i["category"] == "extra_cds"
    )
    assert missing["severity"] == "medium"
    assert extra["severity"] == "low"
    assert "ORF pequeño" in extra["message"]


# ---------------------------------------------------------------------------
# Bordes, largo y estructura
# ---------------------------------------------------------------------------
def test_boundary_shift_and_length_discordance(tmp_path: Path) -> None:
    ref = [
        ("ctg", "CDS", 100, 999, "+", "ID=r1;product=DNA polymerase"),
        ("ctg", "CDS", 2000, 2899, "+", "ID=r2;product=RNA polymerase"),
    ]
    ours = [
        ("ctg", "CDS", 160, 999, "+", "ID=o1;product=DNA polymerase"),
        ("ctg", "CDS", 2000, 2500, "+", "ID=o2;product=RNA polymerase"),
    ]
    data = run_pair(tmp_path, ref, ours)
    cats = categories(data)

    assert cats["boundary_shift"] == 2
    assert cats["length_discordance"] == 1
    assert data["summary"]["CDS"]["exact"] == 0
    boundary_items = [
        i for i in data["differences"]["actionable"]["items"] if i["category"] == "boundary_shift"
    ]
    assert any(i["metrics"]["five_prime_delta_nt"] == 60 for i in boundary_items)
    assert any(i["severity"] == "medium" for i in boundary_items)


def test_multipart_cds_is_grouped_and_structure_reported(tmp_path: Path) -> None:
    ref_rows = [
        ("ctg", "CDS", 100, 200, "+", "ID=cds1;product=DNA polymerase"),
        ("ctg", "CDS", 300, 400, "+", "ID=cds1;product=DNA polymerase"),
    ]
    ours_rows = [("ctg", "CDS", 100, 400, "+", "ID=o1;product=DNA polymerase")]
    ref = write_gff(tmp_path / "reference.gff3", ref_rows)
    parsed_ref = ap.parse_gff(ref)
    assert parsed_ref["feature_counts"]["CDS"] == 1
    assert parsed_ref["multipart_cds"] == 1

    data = run_pair(tmp_path, ref_rows, ours_rows)
    assert data["summary"]["CDS"]["matched"] == 1
    assert data["summary"]["CDS"]["exact"] == 1
    assert data["summary"]["CDS"]["multipart_structure_mismatch"] == 1
    assert categories(data)["multipart_structure_mismatch"] == 1


# ---------------------------------------------------------------------------
# Fragmentación, fusión y hebra
# ---------------------------------------------------------------------------
def test_fragmented_gene_in_ours(tmp_path: Path) -> None:
    ref = [("ctg", "CDS", 1000, 1500, "+", "ID=r1;product=DNA polymerase")]
    ours = [
        ("ctg", "CDS", 1000, 1166, "+", "ID=o1;product=DNA polymerase"),
        ("ctg", "CDS", 1167, 1333, "+", "ID=o2;product=DNA polymerase"),
        ("ctg", "CDS", 1334, 1500, "+", "ID=o3;product=DNA polymerase"),
    ]
    data = run_pair(tmp_path, ref, ours)
    cats = categories(data)

    assert cats["fragmented_in_ours"] == 1
    assert cats.get("missing_cds", 0) == 0
    assert cats.get("extra_cds", 0) == 0
    item = next(
        i
        for i in data["differences"]["actionable"]["items"]
        if i["category"] == "fragmented_in_ours"
    )
    assert len(item["ours_features"]) == 3
    assert item["metrics"]["pieces"] == 3
    assert item["severity"] == "high"


def test_merged_genes_in_ours(tmp_path: Path) -> None:
    ref = [
        ("ctg", "CDS", 1000, 1100, "+", "ID=r1;product=protein A"),
        ("ctg", "CDS", 1200, 1300, "+", "ID=r2;product=protein B"),
    ]
    ours = [("ctg", "CDS", 1000, 1400, "+", "ID=o1;product=protein A")]
    data = run_pair(tmp_path, ref, ours)
    cats = categories(data)

    assert cats["merged_in_ours"] == 1
    assert cats.get("missing_cds", 0) == 0
    item = next(
        i
        for i in data["differences"]["actionable"]["items"]
        if i["category"] == "merged_in_ours"
    )
    assert len(item["reference_features"]) == 2


def test_strand_conflict(tmp_path: Path) -> None:
    ref = [("ctg", "CDS", 100, 400, "+", "ID=r1;product=DNA polymerase")]
    ours = [("ctg", "CDS", 100, 400, "-", "ID=o1;product=DNA polymerase")]
    data = run_pair(tmp_path, ref, ours)
    cats = categories(data)

    assert cats["strand_conflict"] == 1
    assert data["summary"]["CDS"]["matched"] == 0
    item = next(
        i for i in data["differences"]["actionable"]["items"] if i["category"] == "strand_conflict"
    )
    assert item["reference"]["strand"] == "+"
    assert item["ours"]["strand"] == "-"


# ---------------------------------------------------------------------------
# Productos y ncRNA
# ---------------------------------------------------------------------------
def test_products_hypothetical_and_ncrna(tmp_path: Path) -> None:
    ref = [
        ("ctg", "CDS", 100, 399, "+", "ID=r1;product=hypothetical protein"),
        ("ctg", "CDS", 500, 899, "+", "ID=r2;product=DNA polymerase"),
        ("ctg", "tmRNA", 1000, 1200, "+", "ID=nc1;product=tmRNA"),
    ]
    ours = [
        ("ctg", "CDS", 100, 399, "+", "ID=o1;product=ABC transporter"),
        ("ctg", "CDS", 500, 899, "+", "ID=o2;product=DNA-directed DNA polymerase"),
    ]
    data = run_pair(tmp_path, ref, ours)

    assert categories(data)["hypothetical_disagreement"] == 1
    assert categories(data)["missing_ncrna"] == 1
    assert data["ncrna"]["status"] == "ours_absent"
    assert data["ncrna"]["reference_total"] == 1
    assert data["ncrna"]["missing"] == 1
    assert data["ncrna"]["reference_types"] == {"tmRNA": 1}
    assert data["products"]["ours"]["pct_hypothetical"] == 0.0
    assert data["products"]["reference"]["pct_hypothetical"] == 50.0
    assert data["products"]["paired"]["hypothetical_disagreements"] == 1
    assert data["products"]["paired"]["product_mismatches"] == 0


def test_product_mismatch_is_review_not_actionable(tmp_path: Path) -> None:
    ref = [("ctg", "CDS", 100, 399, "+", "ID=r1;product=DNA polymerase")]
    ours = [("ctg", "CDS", 100, 399, "+", "ID=o1;product=ribosomal protein S1")]
    data = run_pair(tmp_path, ref, ours)

    assert data["differences"]["review"]["by_category"]["product_mismatch"] == 1
    assert data["products"]["paired"]["product_mismatches"] == 1
    assert categories(data).get("hypothetical_disagreement", 0) == 0


def test_no_products_does_not_invent_hypothetical(tmp_path: Path) -> None:
    rows = [("ctg", "CDS", 100, 399, "+", "ID=cds1")]
    data = run_pair(tmp_path, rows, rows)

    assert data["products"]["ours"]["status"] == "not_available"
    assert data["products"]["ours"]["pct_hypothetical"] is None
    assert data["products"]["reference"]["pct_hypothetical"] is None
    assert data["products"]["paired"]["status"] == "not_evaluable"


# ---------------------------------------------------------------------------
# Sin referencia: not_evaluable, sin números inventados
# ---------------------------------------------------------------------------
def test_without_reference_is_not_evaluable(tmp_path: Path) -> None:
    ours = write_gff(
        tmp_path / "ours.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=cds1;product=DNA polymerase")],
    )
    data = ap.run_parity(make_args(tmp_path, ours))

    assert data["status"] == "not_evaluable"
    assert "referencia" in data["reason"].lower()
    assert data["inputs"]["reference"] is None
    assert data["summary"]["CDS"]["reference_total"] is None
    assert data["summary"]["CDS"]["matched"] is None
    assert data["summary"]["CDS"]["recall"] is None
    assert data["summary"]["CDS"]["precision"] is None
    assert data["summary"]["CDS"]["f1"] is None
    assert data["coordinates"]["status"] == "not_evaluable"
    assert data["coordinates"]["micro_precision"] is None
    assert data["differences"]["actionable"]["total"] is None
    assert data["products"]["reference"] is None
    assert data["products"]["paired"]["status"] == "not_evaluable"
    assert data["small_orfs"]["reference"] is None
    assert data["small_orfs"]["missing_small"] is None
    assert data["ncrna"]["reference_total"] is None
    assert data["summary"]["CDS"]["ours_total"] == 1


def test_missing_reference_path_is_not_evaluable(tmp_path: Path) -> None:
    ours = write_gff(
        tmp_path / "ours.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=cds1;product=DNA polymerase")],
    )
    missing = tmp_path / "nope.gff3"
    data = ap.run_parity(make_args(tmp_path, ours, missing))

    assert data["status"] == "not_evaluable"
    assert "no existe" in data["reason"].lower()


def test_cli_writes_not_evaluable_reports(tmp_path: Path) -> None:
    ours = write_gff(
        tmp_path / "ours.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=cds1;product=DNA polymerase")],
    )
    md = tmp_path / "reports" / "parity.md"
    json_path = tmp_path / "reports" / "parity.json"
    rc = ap.main(
        [
            "--ours",
            str(ours),
            "--root",
            str(tmp_path),
            "--md",
            str(md),
            "--json",
            str(json_path),
        ]
    )

    assert rc == 0
    assert md.is_file()
    assert json_path.is_file()
    report = json.loads(json_path.read_text(encoding="utf-8"))
    assert report["status"] == "not_evaluable"
    text = md.read_text(encoding="utf-8")
    assert "No evaluable" in text
    assert "No se calculan precision/recall" in text


# ---------------------------------------------------------------------------
# Descubrimiento, gzip, functional.tsv y contigs
# ---------------------------------------------------------------------------
def test_autodiscovery_gzip_and_functional_tsv(tmp_path: Path) -> None:
    genome = "g1"
    ours_dir = tmp_path / "data" / "annotations_ab_initio" / genome
    write_gff(ours_dir / "annotation.gff3", [("ctg", "CDS", 100, 399, "+", "ID=o1")])
    (ours_dir / "functional.tsv").write_text(
        "o1\tsp|Q1|TEST DNA polymerase OS=E. coli OX=1\n", encoding="utf-8"
    )
    write_gff(
        tmp_path / "data" / "annotations" / f"{genome}.gff.gz",
        [("ctg", "CDS", 100, 399, "+", "ID=r1;product=DNA polymerase")],
        gz=True,
    )
    args = ap.build_parser().parse_args(
        [
            "--genome",
            genome,
            "--root",
            str(tmp_path),
            "--md",
            str(tmp_path / "reports" / "parity.md"),
            "--json",
            str(tmp_path / "reports" / "parity.json"),
        ]
    )
    args._argv = []  # type: ignore[attr-defined]
    data = ap.run_parity(args)

    assert data["status"] == "evaluated"
    assert data["inputs"]["ours"]["path"].endswith("annotation.gff3")
    assert data["inputs"]["reference"]["path"].endswith(".gff.gz")
    assert data["products"]["ours"]["source"] == "TSV vecino: functional.tsv"
    assert data["products"]["ours"]["with_product"] == 1
    assert data["summary"]["CDS"]["matched"] == 1
    assert data["products"]["paired"]["literal_agreement"] == 1.0


def test_out_fl_discovery(tmp_path: Path) -> None:
    genome = "outg"
    write_gff(
        tmp_path / "out_fl" / "annotation" / f"{genome}.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=o1;product=DNA polymerase")],
    )
    write_gff(
        tmp_path / "data" / "annotations" / f"{genome}.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=r1;product=DNA polymerase")],
    )
    args = ap.build_parser().parse_args(
        ["--genome", genome, "--root", str(tmp_path), "--md", str(tmp_path / "m.md"), "--json", str(tmp_path / "j.json")]
    )
    args._argv = []  # type: ignore[attr-defined]
    data = ap.run_parity(args)

    assert data["inputs"]["ours"]["path"].endswith(f"out_fl/annotation/{genome}.gff3")
    assert data["status"] == "evaluated"
    assert data["summary"]["CDS"]["matched"] == 1


def test_contig_mismatch_makes_coordinates_not_evaluable(tmp_path: Path) -> None:
    ref = [("refctg", "CDS", 100, 399, "+", "ID=r1;product=DNA polymerase")]
    ours = [("ourctg", "CDS", 100, 399, "+", "ID=o1;product=DNA polymerase")]
    data = run_pair(tmp_path, ref, ours)

    assert data["status"] == "evaluated"
    assert data["coordinates"]["status"] == "not_evaluable"
    assert data["summary"]["CDS"]["status"] == "not_evaluable_contigs"
    assert data["summary"]["CDS"]["matched"] is None
    assert data["summary"]["CDS"]["recall"] is None
    assert data["differences"]["actionable"]["by_category"]["contig_mismatch"] == 1
    assert data["small_orfs"]["missing_small"] == 0


def test_pseudo_cds_excluded_by_default(tmp_path: Path) -> None:
    path = write_gff(
        tmp_path / "pseudo.gff3",
        [
            ("ctg", "CDS", 100, 399, "+", "ID=ok;product=DNA polymerase"),
            ("ctg", "CDS", 500, 799, "+", "ID=pseudo;pseudo=true;product=fragment"),
        ],
    )
    parsed = ap.parse_gff(path)
    assert parsed["feature_counts"]["CDS"] == 1
    assert parsed["excluded_pseudo"] == 1
    kept = ap.parse_gff(path, exclude_pseudo=False)
    assert kept["feature_counts"]["CDS"] == 2


def test_ours_cds_over_pseudo_reference_is_review_not_extra(tmp_path: Path) -> None:
    ref = [
        ("ctg", "CDS", 100, 399, "+", "ID=r1;product=DNA polymerase"),
        ("ctg", "CDS", 500, 799, "+", "ID=pseudo;pseudo=true;product=fragment"),
    ]
    ours = [
        ("ctg", "CDS", 100, 399, "+", "ID=o1;product=DNA polymerase"),
        ("ctg", "CDS", 500, 799, "+", "ID=o2;product=hypothetical protein"),
    ]
    data = run_pair(tmp_path, ref, ours)

    assert data["summary"]["CDS"]["reference_total"] == 1
    assert data["summary"]["CDS"]["ours_total"] == 2
    assert data["differences"]["review"]["by_category"]["pseudo_in_reference"] == 1
    assert data["differences"]["actionable"]["by_category"].get("extra_cds", 0) == 0


def test_max_details_truncates_but_keeps_total(tmp_path: Path) -> None:
    ref = [("ctg", "CDS", 100, 399, "+", "ID=r0;product=protein 0")]
    ours_rows = [
        ("ctg", "CDS", 100 + i * 500, 399 + i * 500, "+", f"ID=o{i};product=protein {i}")
        for i in range(6)
    ]
    data = run_pair(tmp_path, ref, ours_rows, "--max-details", "2")

    assert data["differences"]["actionable"]["total"] == 5
    assert len(data["differences"]["actionable"]["items"]) == 2
    assert data["differences"]["actionable"]["truncated"] is True


def test_reference_without_comparable_features_is_not_evaluable(tmp_path: Path) -> None:
    ref = write_gff(
        tmp_path / "reference.gff3",
        [("ctg", "gene", 100, 399, "+", "ID=gene1;product=DNA polymerase")],
    )
    ours = write_gff(
        tmp_path / "ours.gff3",
        [("ctg", "CDS", 100, 399, "+", "ID=o1;product=DNA polymerase")],
    )
    data = ap.run_parity(make_args(tmp_path, ours, ref))

    assert data["status"] == "not_evaluable"
    assert "no contiene" in data["reason"]
