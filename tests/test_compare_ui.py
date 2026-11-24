"""Pruebas del comparador GFF3 (``webapp/compare_ui.py``).

Los GFF son sintéticos: no se descarga nada ni se depende de anotaciones
reales. Se comprueban los conteos, la cascada locus_tag → ID → coordenadas,
los genomas sin solapamiento y la tabla/TSV.
"""
from __future__ import annotations

import csv
import gzip
import io
import sys
from pathlib import Path

import pytest

import compare_ui


def _write_gff(tmp_path, name: str, rows: list[str]) -> Path:
    path = tmp_path / name
    path.write_text(
        "##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8"
    )
    return path


def _line(
    seqid: str,
    start: int,
    end: int,
    strand: str,
    attrs: str,
    ftype: str = "CDS",
) -> str:
    phase = "." if ftype != "CDS" else "0"
    return "\t".join(
        [seqid, "test", ftype, str(start), str(end), ".", strand, phase, attrs]
    )


def _gff_a() -> list[str]:
    return [
        _line("chr1", 100, 400, "+",
              "ID=gene-1;locus_tag=LT_0001;Name=gyrA;product=DNA gyrase",
              "gene"),
        _line("chr1", 100, 400, "+", "ID=cds-1;Parent=gene-1"),
        _line("chr1", 500, 800, "+",
              "ID=gene-2;locus_tag=LT_0002;Name=recA;product=RecA",
              "gene"),
        _line("chr1", 500, 800, "+", "ID=cds-2;Parent=gene-2"),
        _line("chr1", 900, 1100, "-",
              "ID=gene-3;locus_tag=LT_0003;Name=exclusive_a",
              "gene"),
        _line("chr1", 900, 1100, "-", "ID=cds-3;Parent=gene-3"),
    ]


def _gff_b() -> list[str]:
    return [
        _line("contigB", 200, 500, "+",
              "ID=gene-10;locus_tag=LT_0001;Name=gyrA;product=DNA gyrase",
              "gene"),
        _line("contigB", 200, 500, "+", "ID=cds-10;Parent=gene-10"),
        _line("contigB", 600, 900, "+",
              "ID=gene-20;locus_tag=LT_0002;Name=recA;product=RecA",
              "gene"),
        _line("contigB", 600, 900, "+", "ID=cds-20;Parent=gene-20"),
        _line("contigB", 1200, 1400, "-",
              "ID=gene-30;locus_tag=LT_0004;Name=exclusive_b",
              "gene"),
        _line("contigB", 1200, 1400, "-", "ID=cds-30;Parent=gene-30"),
    ]


# ------------------------------------------------------------------ parsing
def test_counts_and_shared_by_locus_tag(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", _gff_a())
    b = _write_gff(tmp_path, "b.gff3", _gff_b())

    result = compare_ui.compare_gff(a, b, label_a="K12", label_b="EC")

    counts = result["counts"]
    assert counts["genes_a"] == 3
    assert counts["genes_b"] == 3
    assert counts["cds_a"] == 3
    assert counts["compartidos"] == 2
    assert counts["compartidos_locus_tag"] == 2
    assert counts["compartidos_id"] == 0
    assert counts["compartidos_coordenadas"] == 0
    assert counts["exclusivos_a"] == 1
    assert counts["exclusivos_b"] == 1
    assert counts["genes_por_genoma"] == {"K12": 3, "EC": 3}

    assert result["labels"] == {"a": "K12", "b": "EC"}
    assert len(result["shared"]) == 2
    assert len(result["exclusive_a"]) == 1
    assert len(result["exclusive_b"]) == 1
    assert len(result["table"]) == 4
    assert list(result["table"].columns) == compare_ui.COLUMNS
    assert set(result["table"]["estado"]) == {"compartido", "solo K12", "solo EC"}

    row = result["table"][
        (result["table"]["estado"] == "compartido")
        & (result["table"]["locus_tag_a"] == "LT_0001")
    ].iloc[0]
    assert row["coincidencia"] == "locus_tag"
    assert row["gen_a"] == "LT_0001"
    assert row["gen_b"] == "LT_0001"
    assert row["inicio_a"] == 100
    assert row["inicio_b"] == 200
    assert row["delta_inicio"] == 100
    assert row["delta_fin"] == 100
    assert bool(row["misma_hebra"]) is True
    assert row["producto_a"] == "DNA gyrase"
    assert "no se verificó la secuencia" in row["nota"]


def test_multipart_cds_are_one_gene_and_span_aggregates(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 100, 600, "+",
              "ID=gene-2;locus_tag=LT_0002;Name=split", "gene"),
        _line("chr1", 100, 250, "+", "ID=cds-2a;Parent=gene-2"),
        _line("chr1", 400, 600, "+", "ID=cds-2b;Parent=gene-2"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr2", 100, 600, "+",
              "ID=gene-20;locus_tag=LT_0002;Name=whole", "gene"),
        _line("chr2", 100, 600, "+", "ID=cds-20;Parent=gene-20"),
    ])

    result = compare_ui.compare_gff(a, b)
    genes_a = result["genes_a"]
    assert len(genes_a) == 1
    gene = genes_a[0]
    assert gene["start"] == 100
    assert gene["end"] == 600
    assert gene["cds_count"] == 2
    assert gene["cds_bp"] == 151 + 201
    assert gene["id"] == "gene-2"
    assert result["counts"]["cds_a"] == 2
    assert result["counts"]["compartidos"] == 1


def test_nested_gene_mrna_cds_inherits_locus_tag_and_gene_id():
    text = (
        "##gff-version 3\n"
        + _line("chr", 100, 500, "+",
                "ID=gene-1;locus_tag=LT_1;Name=abc", "gene") + "\n"
        + _line("chr", 100, 500, "+",
                "ID=rna-1;Parent=gene-1;product=ABC protein", "mRNA") + "\n"
        + _line("chr", 100, 200, "+", "ID=cds-1;Parent=rna-1") + "\n"
        + _line("chr", 300, 500, "+", "ID=cds-2;Parent=rna-1") + "\n"
    )
    genes, warnings = compare_ui.parse_cds(io.StringIO(text))

    assert warnings == []
    assert len(genes) == 1
    gene = genes[0]
    assert gene["locus_tag"] == "LT_1"
    assert gene["gene_id"] == "gene-1"
    assert gene["id"] == "gene-1"
    assert gene["product"] == "ABC protein"
    assert gene["start"] == 100
    assert gene["end"] == 500
    assert gene["cds_count"] == 2


def test_shared_by_id_when_there_is_no_locus_tag(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=shared_id;product=alpha"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("other", 50, 70, "+", "ID=shared_id;product=beta"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["compartidos"] == 1
    assert result["counts"]["compartidos_id"] == 1
    row = result["table"].iloc[0]
    assert row["coincidencia"] == "ID"
    assert row["gen_a"] == "shared_id"
    assert row["producto_a"] == "alpha"
    assert row["producto_b"] == "beta"
    assert any("por ID de anotación" in w for w in result["warnings"])


def test_conflicting_locus_tags_block_id_matching(tmp_path):
    """Un ID igual no debe vencer dos locus_tag explícitamente distintos."""
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=same_id;locus_tag=A_0001"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr1", 200, 230, "+", "ID=same_id;locus_tag=B_0001"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["compartidos"] == 0
    assert result["counts"]["compartidos_id"] == 0
    assert result["counts"]["exclusivos_a"] == 1
    assert result["counts"]["exclusivos_b"] == 1


def test_shared_by_exact_coordinates_is_flagged_as_positional(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=only_a;product=alpha"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr1", 10, 30, "+", "ID=only_b;product=beta"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["compartidos"] == 1
    assert result["counts"]["compartidos_coordenadas"] == 1
    row = result["table"].iloc[0]
    assert row["coincidencia"] == "coordenadas"
    assert "coincidencia posicional" in row["nota"]
    assert "no implica ortología" in row["nota"]
    assert any("solo por coordenadas" in w for w in result["warnings"])
    assert any("no demuestra ortología" in w for w in result["warnings"])


def test_coordinates_do_not_match_across_contigs(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=a1"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr2", 10, 30, "+", "ID=b1"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["compartidos"] == 0
    assert result["counts"]["exclusivos_a"] == 1
    assert result["counts"]["exclusivos_b"] == 1


def test_genomes_without_overlap_have_no_shared_genes(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=a1;product=alpha"),
        _line("chr1", 100, 130, "+", "ID=a2;product=beta"),
    ])
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr1", 500, 530, "+", "ID=b1;product=gamma"),
        _line("chr1", 1000, 1030, "+", "ID=b2;product=delta"),
    ])

    result = compare_ui.compare_gff(a, b)
    counts = result["counts"]
    assert counts["genes_a"] == 2
    assert counts["genes_b"] == 2
    assert counts["compartidos"] == 0
    assert counts["exclusivos_a"] == 2
    assert counts["exclusivos_b"] == 2
    assert len(result["table"]) == 4
    assert list(result["table"]["coincidencia"]) == [None, None, None, None]
    assert set(result["table"]["estado"]) == {"solo A", "solo B"}
    assert any("No hay genes compartidos" in w for w in result["warnings"])


def test_fasta_section_is_ignored_and_percent_encoding_is_decoded(tmp_path):
    a = tmp_path / "a.gff3"
    a.write_text(
        "##gff-version 3\n"
        + _line("chr1", 10, 30, "+",
                "ID=a1;product=DNA%20gyrase%20subunit%20A")
        + "\n##FASTA\n>chr1\n"
        + _line("chr1", 999, 1100, "+", "ID=should_not_count")
        + "\n",
        encoding="utf-8",
    )
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr1", 10, 30, "+",
              "ID=b1;product=DNA%20gyrase%20subunit%20A"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["genes_a"] == 1
    assert result["genes_a"][0]["product"] == "DNA gyrase subunit A"
    assert result["counts"]["compartidos_coordenadas"] == 1


def test_malformed_lines_warn_but_do_not_abort(tmp_path):
    a = tmp_path / "a.gff3"
    a.write_text(
        "##gff-version 3\n"
        "linea con espacios\n"
        + _line("chr1", 10, 30, "+", "ID=ok")
        + "\n"
        + "\t".join(["chr1", "test", "CDS", "no_num", "30", ".", "+", "0",
                     "ID=bad"])
        + "\n",
        encoding="utf-8",
    )
    b = _write_gff(tmp_path, "b.gff3", [
        _line("chr1", 10, 30, "+", "ID=ok"),
    ])

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["genes_a"] == 1
    assert any("9 columnas" in w for w in result["warnings"])
    assert any("coordenadas no numéricas" in w for w in result["warnings"])


def test_reads_gzip_and_file_like_objects(tmp_path):
    a = tmp_path / "a.gff3.gz"
    with gzip.open(a, "wt", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write(_line("chr1", 10, 30, "+", "ID=z1") + "\n")
    b = io.BytesIO(
        ("##gff-version 3\n" + _line("chr1", 10, 30, "+", "ID=z2") + "\n")
        .encode("utf-8")
    )

    result = compare_ui.compare_gff(a, b)
    assert result["counts"]["genes_a"] == 1
    assert result["counts"]["genes_b"] == 1
    assert result["counts"]["compartidos_coordenadas"] == 1


def test_tsv_has_stable_columns_and_rows(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", _gff_a())
    b = _write_gff(tmp_path, "b.gff3", _gff_b())
    result = compare_ui.compare_gff(a, b)

    tsv = compare_ui.to_tsv(result)
    lines = tsv.splitlines()
    assert lines[0].split("\t") == compare_ui.COLUMNS
    assert len(lines) == 5
    rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
    assert len(rows) == 4
    assert rows[0]["estado"] == "compartido"
    assert "coincidencia" in rows[0]


def test_errors_and_validation(tmp_path):
    a = _write_gff(tmp_path, "a.gff3", [
        _line("chr1", 10, 30, "+", "ID=a1"),
    ])
    with pytest.raises(compare_ui.CompareError):
        compare_ui.compare_gff(tmp_path / "missing.gff3", a)
    with pytest.raises(compare_ui.CompareError):
        compare_ui.compare_gff(a, a, label_a="X", label_b="X")
    with pytest.raises(compare_ui.CompareError):
        compare_ui.to_tsv({})


def test_limits_state_no_sequence_no_network():
    assert "no es ortología ni función" in compare_ui.CAPTION_COMPARE
    text = " ".join(compare_ui.LIMITS).lower()
    assert "no lee ni alinea adn" in text
    assert "no consulta la red" in text
    assert "identidad de secuencia" in text


# --------------------------------------------------------------------- UI
class _FakeCtx:
    def __init__(self, calls: list) -> None:
        self.calls = calls

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "columns":
                count = int(args[0]) if args else 1
                return [_FakeCtx(self.calls) for _ in range(count)]
            return _FakeCtx(self.calls)

        return recorder


class _FakeStreamlit:
    def __init__(self, click=False, uploads=None) -> None:
        self.calls: list = []
        self.session_state: dict = {}
        self.click = click
        self.uploads = uploads or {}

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "columns":
                count = int(args[0]) if args else 1
                return [_FakeCtx(self.calls) for _ in range(count)]
            if name == "file_uploader":
                return self.uploads.get(kwargs.get("key"))
            if name == "text_input":
                return kwargs.get("value", "")
            if name == "button":
                return self.click
            return _FakeCtx(self.calls)

        return recorder


def test_render_compare_has_two_sources_and_explicit_caption(monkeypatch):
    fake = _FakeStreamlit(click=False)
    monkeypatch.setitem(sys.modules, "streamlit", fake)

    compare_ui.render_compare()

    uploaders = [c for c in fake.calls if c[0] == "file_uploader"]
    assert len(uploaders) == 2
    captions = " ".join(str(c[1][0]) for c in fake.calls if c[0] == "caption")
    assert "no es ortología ni función" in captions
    assert any(c[0] == "info" for c in fake.calls)


def test_render_compare_shows_table_and_tsv_download(monkeypatch):
    upload_a = io.BytesIO(
        ("##gff-version 3\n" + "\n".join(_gff_a()) + "\n").encode("utf-8")
    )
    upload_a.name = "a.gff3"
    upload_b = io.BytesIO(
        ("##gff-version 3\n" + "\n".join(_gff_b()) + "\n").encode("utf-8")
    )
    upload_b.name = "b.gff3"
    fake = _FakeStreamlit(
        click=True,
        uploads={"compare_file_a": upload_a, "compare_file_b": upload_b},
    )
    monkeypatch.setitem(sys.modules, "streamlit", fake)

    compare_ui.render_compare()

    metrics = [c for c in fake.calls if c[0] == "metric"]
    assert any(c[1][0] == "Genes compartidos" and c[1][1] == 2
               for c in metrics)
    downloads = [c for c in fake.calls if c[0] == "download_button"]
    assert len(downloads) == 1
    data = downloads[0][2]["data"]
    assert data.splitlines()[0].split("\t") == compare_ui.COLUMNS
    assert len(data.splitlines()) == 5
