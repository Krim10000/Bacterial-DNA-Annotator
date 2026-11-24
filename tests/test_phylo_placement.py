"""Tests de la colocación tentativa por vecinos de k-mers.

No usan red, GPU, FASTA reales ni las referencias de 400 genomas: construyen
secuencias, FASTA y metadatos sintéticos en ``tmp_path``. Comprueban el
contrato de ``sketch``/``nearest_genomes``/``placement``, la equivalencia de
los backends y que un vecino k-mer no se convierta en una especie inventada.
"""
from __future__ import annotations

import gzip
import json
import random
from pathlib import Path

import pytest

from webapp import phylo_placement as phylo


# ---------------------------------------------------------------------------
# Utilidades sintéticas
# ---------------------------------------------------------------------------

def _seq(length: int, seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _mutate(sequence: str, rate: float, seed: int) -> str:
    rng = random.Random(seed)
    bases = "ACGT"
    out = []
    for base in sequence:
        if rng.random() < rate:
            out.append(rng.choice([b for b in bases if b != base]))
        else:
            out.append(base)
    return "".join(out)


def _fasta_text(records: list[tuple[str, str]]) -> str:
    parts = []
    for name, sequence in records:
        parts.append(f">{name}\n")
        parts.extend(
            sequence[start:start + 70] + "\n"
            for start in range(0, len(sequence), 70)
        )
    return "".join(parts)


def _write_fasta(
    path: Path, records: list[tuple[str, str]], compress: bool = False
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = _fasta_text(records)
    if compress or path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="ascii") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="ascii")
    return path


def _data_root(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    (data / "raw").mkdir(parents=True, exist_ok=True)
    (data / "raw_excluded").mkdir(parents=True, exist_ok=True)
    return data


def _taxonomy_record(
    genome: str,
    lineage: str | None,
    species: str | None = None,
    confidence: str = "high",
    **extra: object,
) -> dict:
    record: dict = {"genome": genome, "confidence": confidence}
    if lineage is not None:
        record["gtdb_lineage"] = lineage
    if species is not None:
        record["species_gtdb"] = species
    record.update(extra)
    return record


def _write_taxonomy(root: Path, genome: str, record: dict) -> Path:
    folder = root / "taxonomy" / "full"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{genome}.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# sketch / sketch_fasta
# ---------------------------------------------------------------------------

def test_sketch_deterministic_bounded_and_strand_invariant():
    sequence = _seq(3000, seed=1)
    sketch = phylo.sketch(sequence)
    assert isinstance(sketch, set)
    assert 0 < len(sketch) <= phylo.DEFAULT_SKETCH_SIZE
    assert all(isinstance(value, int) for value in sketch)
    assert all(0 <= value < 2**64 for value in sketch)
    assert sketch == phylo.sketch(sequence)
    assert sketch == phylo.sketch(sequence.lower())
    assert sketch == phylo.sketch(_revcomp(sequence))
    assert sketch != phylo.sketch(_seq(3000, seed=2))


def test_sketch_short_and_ambiguous_sequences():
    assert phylo.sketch("ACGT") == set()
    assert phylo.sketch("N" * 500) == set()
    assert phylo.sketch("ACGTNACGT", k=6) == set()

    left = _seq(400, seed=3)
    right = _seq(400, seed=4)
    with_gap = left + "N" * 25 + right
    expected = (
        phylo.sketch(left, size=5000) | phylo.sketch(right, size=5000)
    )
    assert phylo.sketch(with_gap, size=5000) == expected


def test_sketch_treats_u_as_t_and_ignores_whitespace():
    assert phylo.sketch("acgtu", k=5) == phylo.sketch("ACGTT", k=5)
    assert phylo.sketch("AC GT\nACGT", k=3) == phylo.sketch("ACGTACGT", k=3)


def test_sketch_supports_k_above_32():
    sequence = _seq(600, seed=5)
    sketch = phylo.sketch(sequence, k=40, size=128)
    assert sketch
    assert len(sketch) <= 128
    assert sketch == phylo.sketch(sequence, k=40, size=128)
    assert sketch == phylo.sketch(_revcomp(sequence), k=40, size=128)


def test_sketch_parametres_are_validated():
    with pytest.raises(ValueError):
        phylo.sketch("ACGT", k=0)
    with pytest.raises(ValueError):
        phylo.sketch("ACGT", size=0)
    with pytest.raises(TypeError):
        phylo.sketch(12345)  # type: ignore[arg-type]


def test_sketch_fasta_plain_gzip_and_multirecord(tmp_path: Path):
    left = _seq(900, seed=11)
    right = _seq(800, seed=12)
    records = [("contig_1", left), ("contig_2", right)]
    plain = _write_fasta(tmp_path / "genome.fna", records)
    compressed = _write_fasta(tmp_path / "genome.fna.gz", records)

    sketch_plain = phylo.sketch_fasta(plain, size=5000)
    sketch_gzip = phylo.sketch_fasta(compressed, size=5000)
    expected = (
        phylo.sketch(left, size=5000) | phylo.sketch(right, size=5000)
    )
    assert sketch_plain == sketch_gzip == expected
    # Sin cruzar la frontera entre registros: el sketch es la unión.
    assert phylo.sketch(left + right, size=5000) != expected


def test_sketch_fasta_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        phylo.sketch_fasta(tmp_path / "nope.fna")


@pytest.mark.skipif(phylo._np is None, reason="numpy no está disponible")
@pytest.mark.parametrize(("length", "k"), [(700, 21), (3000, 21), (5000, 15)])
def test_numpy_and_python_backends_are_equivalent(length: int, k: int):
    sequence = phylo._normalize_sequence(_seq(length, seed=length))
    python_sketch = phylo._sketch_sequence_python(sequence, k, 512)
    numpy_sketch = phylo._sketch_sequence_numpy(sequence, k, 512)
    assert numpy_sketch is not None
    assert numpy_sketch == python_sketch

    with_gap = sequence[: length // 2] + b"N" * 7 + sequence[length // 2 + 7:]
    assert phylo._sketch_sequence_numpy(with_gap, k, 512) == (
        phylo._sketch_sequence_python(with_gap, k, 512)
    )


# ---------------------------------------------------------------------------
# nearest_genomes
# ---------------------------------------------------------------------------

def _panel(tmp_path: Path):
    """Tres referencias: idéntica, mutada y no relacionada."""
    root = _data_root(tmp_path)
    query = _seq(3200, seed=21)
    close = _mutate(query, rate=0.05, seed=22)
    far = _seq(3200, seed=23)

    genome_exact = "GCF_900000001.1_ASM1v1"
    genome_close = "GCF_900000002.1_ASM2v1"
    genome_far = "GCF_900000003.1_ASM3v1"

    _write_fasta(
        root / "raw" / f"{genome_exact}_genomic.fna",
        [("chromosome", query)],
    )
    _write_fasta(
        root / "raw" / f"{genome_close}_genomic.fna.gz",
        [("chromosome", close)],
    )
    _write_fasta(
        root / "raw_excluded" / f"{genome_far}_genomic.fna",
        [("chromosome", far)],
    )
    return root, query, close, far, genome_exact, genome_close, genome_far


def test_nearest_genomes_ranks_and_returns_local_lineage(tmp_path: Path):
    (root, query, _close, _far, genome_exact, genome_close,
     genome_far) = _panel(tmp_path)
    lineage = (
        "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
        "s__Escherichia coli"
    )
    _write_taxonomy(
        root,
        genome_exact,
        _taxonomy_record(
            genome_exact, lineage, species="Escherichia coli"
        ),
    )

    nearest = phylo.nearest_genomes(query, root, top_n=3)

    assert [item["genome"] for item in nearest] == [
        genome_exact,
        genome_close,
        genome_far,
    ]
    assert nearest[0]["jaccard"] == 1.0
    assert nearest[0]["accession"] == "GCF_900000001.1"
    assert nearest[0]["source"] == "raw"
    assert nearest[0]["shared_hashes"] > 0
    assert 0.0 < nearest[1]["jaccard"] < 1.0
    assert nearest[2]["jaccard"] == 0.0
    assert nearest[2]["source"] == "raw_excluded"

    top_lineage = nearest[0]["lineage"]
    assert top_lineage["species"] == "Escherichia coli"
    assert top_lineage["ranks"]["phylum"] == "Pseudomonadota"
    assert top_lineage["source_kind"] == "taxonomy_json"
    assert nearest[0]["lineage_text"] == lineage
    # Las referencias sin metadatos no heredan el linaje de la mejor.
    assert nearest[1]["lineage"] is None
    assert nearest[2]["lineage"] is None


def test_nearest_genomes_top_n_and_invalid_top_n(tmp_path: Path):
    root, query, _close, _far, genome_exact, _gc, _gf = _panel(tmp_path)
    nearest = phylo.nearest_genomes(query, root, top_n=1)
    assert len(nearest) == 1
    assert nearest[0]["genome"] == genome_exact
    with pytest.raises(ValueError):
        phylo.nearest_genomes(query, root, top_n=0)


def test_nearest_genomes_accepts_path_sequence_and_list(tmp_path: Path):
    root, query, _close, _far, genome_exact, _gc, _gf = _panel(tmp_path)
    query_path = _write_fasta(tmp_path / "query.fna", [("query", query)])

    from_path = phylo.nearest_genomes(query_path, root, top_n=1)
    from_text = phylo.nearest_genomes(str(query_path), root, top_n=1)
    from_sequence = phylo.nearest_genomes(query, root, top_n=1)
    from_list = phylo.nearest_genomes([query], root, top_n=1)
    from_fasta_text = phylo.nearest_genomes(
        ">query\n" + query + "\n", root, top_n=1
    )

    for result in (from_path, from_text, from_sequence, from_list,
                   from_fasta_text):
        assert len(result) == 1
        assert result[0]["genome"] == genome_exact


def test_nearest_genomes_empty_without_references(tmp_path: Path):
    root = _data_root(tmp_path)
    assert phylo.nearest_genomes(_seq(1000), root) == []


def test_nearest_genomes_raises_for_missing_query_path(tmp_path: Path):
    root = _data_root(tmp_path)
    with pytest.raises(FileNotFoundError):
        phylo.nearest_genomes(str(tmp_path / "missing.fna"), root)


# ---------------------------------------------------------------------------
# placement
# ---------------------------------------------------------------------------

def test_placement_ok_contract(tmp_path: Path):
    root, query, _close, _far, genome_exact, _gc, _gf = _panel(tmp_path)
    lineage = (
        "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
        "s__Escherichia coli"
    )
    _write_taxonomy(
        root,
        genome_exact,
        _taxonomy_record(
            genome_exact, lineage, species="Escherichia coli"
        ),
    )

    result = phylo.placement(query, root)

    assert result["status"] == phylo.STATUS_OK
    assert result["confidence"] == "low"
    assert result["method"] == "kmer_nearest_not_ani"
    assert "ANI" in result["limitation"]
    assert "identificación taxonómica" in result["limitation"]
    assert result["confidence_reason"]
    assert result["reason"] == ""
    assert result["nearest"][0]["genome"] == genome_exact
    assert result["nearest"][0]["jaccard"] == 1.0
    assert result["lineage"]["species"] == "Escherichia coli"
    assert result["lineage_text"] == lineage
    assert result["lineage_from"] == genome_exact
    assert result["best_jaccard"] == 1.0
    assert result["query"]["n_hashes"] > 0
    assert result["query"]["k"] == phylo.DEFAULT_K
    assert result["references"]["n_fastas"] == 3
    assert result["references"]["n_sketched"] == 3
    assert result["references"]["directories"] == ["raw", "raw_excluded"]
    # El contrato no expone ANI real: el método lo declara explícitamente.
    assert "ani" not in {key.lower() for key in result}


def test_placement_not_available_without_references(tmp_path: Path):
    root = _data_root(tmp_path)
    result = phylo.placement(_seq(1000), root)
    assert result["status"] == phylo.STATUS_NOT_AVAILABLE
    assert result["nearest"] == []
    assert result["lineage"] is None
    assert result["lineage_text"] is None
    assert result["lineage_from"] is None
    assert result["confidence"] == "low"
    assert result["method"] == "kmer_nearest_not_ani"
    assert result["limitation"]
    assert result["reason"]
    assert "FASTA" in result["reason"]
    assert result["references"]["n_fastas"] == 0


def test_placement_not_available_without_shared_kmers(tmp_path: Path):
    root = _data_root(tmp_path)
    genome = "GCF_900000010.1_ASM10v1"
    _write_fasta(
        root / "raw" / f"{genome}_genomic.fna",
        [("chromosome", "A" * 2000)],
    )
    result = phylo.placement("C" * 2000, root)
    assert result["status"] == phylo.STATUS_NOT_AVAILABLE
    assert result["nearest"] == []
    assert "compartidos" in result["reason"]
    # nearest_genomes sí puede mostrar el vecino con Jaccard 0.
    nearest = phylo.nearest_genomes("C" * 2000, root)
    assert nearest[0]["jaccard"] == 0.0


def test_placement_not_available_without_valid_query(tmp_path: Path):
    root, _query, _close, _far, _ge, _gc, _gf = _panel(tmp_path)
    for query in ("", "ACGT", "N" * 100):
        result = phylo.placement(query, root)
        assert result["status"] == phylo.STATUS_NOT_AVAILABLE
        assert result["nearest"] == []
        assert "k-mers" in result["reason"]


def test_placement_missing_path_is_not_available(tmp_path: Path):
    root = _data_root(tmp_path)
    result = phylo.placement(str(tmp_path / "missing.fna"), root)
    assert result["status"] == phylo.STATUS_NOT_AVAILABLE
    assert "No se pudo procesar" in result["reason"]


def test_placement_does_not_invent_species_when_json_lacks_it(tmp_path: Path):
    root = _data_root(tmp_path)
    query = _seq(2500, seed=31)
    genome = "GCF_900000020.1_ASM20v1"
    _write_fasta(
        root / "raw" / f"{genome}_genomic.fna",
        [("chromosome", query)],
    )
    lineage = (
        "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__"
    )
    _write_taxonomy(
        root,
        genome,
        _taxonomy_record(genome, lineage, confidence="high"),
    )

    result = phylo.placement(query, root)
    assert result["status"] == phylo.STATUS_OK
    assert result["lineage"]["ranks"]["genus"] == "Escherichia"
    assert result["lineage"]["species"] is None
    assert "species" not in result["lineage"]["ranks"]
    assert result["lineage"]["partial"] is False
    assert result["lineage"]["declared_confidence"] == "high"


def test_placement_table_fallback_does_not_derive_species(tmp_path: Path):
    root = _data_root(tmp_path)
    query = _seq(2500, seed=32)
    accession = "GCF_900000030.1"
    genome = f"{accession}_ASM30v1"
    _write_fasta(
        root / "raw" / f"{genome}_genomic.fna",
        [("chromosome", query)],
    )
    (root / "selected.tsv").write_text(
        "assembly_accession\tspecies_taxid\torganism_name\tftp_path\tphylum\n"
        f"{accession}\t562\tEscherichia coli str. K-12\t\tPseudomonadota\n",
        encoding="utf-8",
    )

    result = phylo.placement(query, root)
    assert result["status"] == phylo.STATUS_OK
    lineage = result["lineage"]
    assert lineage["source_kind"] == "metadata_table"
    assert lineage["partial"] is True
    assert lineage["ranks"] == {"phylum": "Pseudomonadota"}
    # El nombre de organismo no se convierte en especie por inferencia.
    assert lineage["species"] is None
    assert "species" not in lineage["ranks"]
    assert lineage["organism"] == "Escherichia coli str. K-12"


def test_placement_does_not_borrow_lineage_from_another_genome(tmp_path: Path):
    root = _data_root(tmp_path)
    query = _seq(2600, seed=41)
    close = _mutate(query, rate=0.02, seed=42)
    genome_exact = "GCF_900000041.1_ASM41v1"
    genome_close = "GCF_900000042.1_ASM42v1"
    _write_fasta(
        root / "raw" / f"{genome_exact}_genomic.fna",
        [("chromosome", query)],
    )
    _write_fasta(
        root / "raw" / f"{genome_close}_genomic.fna",
        [("chromosome", close)],
    )
    lineage = "d__Bacteria;p__Pseudomonadota;g__Other;s__Other species"
    _write_taxonomy(
        root,
        genome_close,
        _taxonomy_record(genome_close, lineage, species="Other species"),
    )

    result = phylo.placement(query, root)
    assert result["status"] == phylo.STATUS_OK
    assert result["nearest"][0]["genome"] == genome_exact
    assert result["nearest"][0]["lineage"] is None
    assert result["nearest"][1]["lineage"]["species"] == "Other species"
    # El linaje del segundo vecino no se atribuye al primero.
    assert result["lineage"] is None
    assert result["lineage_text"] is None


def test_placement_data_root_accepts_project_data_and_raw(tmp_path: Path):
    root, query, _close, _far, genome_exact, _gc, _gf = _panel(tmp_path)
    for data_root in (root.parent, root, root / "raw"):
        result = phylo.placement(query, data_root, top_n=1)
        assert result["status"] == phylo.STATUS_OK
        assert result["nearest"][0]["genome"] == genome_exact


def test_placement_skips_corrupt_reference(tmp_path: Path):
    root = _data_root(tmp_path)
    query = _seq(2000, seed=51)
    good = "GCF_900000051.1_ASM51v1"
    bad = "GCF_900000052.1_ASM52v1"
    _write_fasta(
        root / "raw" / f"{good}_genomic.fna",
        [("chromosome", query)],
    )
    corrupt = root / "raw" / f"{bad}_genomic.fna.gz"
    corrupt.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")

    result = phylo.placement(query, root)
    assert result["status"] == phylo.STATUS_OK
    assert result["nearest"][0]["genome"] == good
    assert len(result["references"]["skipped"]) == 1
    assert bad in result["references"]["skipped"][0]["file"]
