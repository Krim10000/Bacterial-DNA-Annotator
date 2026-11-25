import gzip
from pathlib import Path

from webapp.browser import _read_fasta_any, load_genome_sequences


def test_read_fasta_gz(tmp_path):
    p = tmp_path / "x.fna.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">c1 desc\nATGC\n>c2\nTTAA\n")
    assert _read_fasta_any(p) == {"c1": "ATGC", "c2": "TTAA"}


def test_load_reference_genome():
    import pytest
    root = Path("data")
    if not (root / "reference_genomes" / "paeruginosa").exists():
        pytest.skip("sin referencias descargadas")
    seqs = load_genome_sequences(root, "paeruginosa")
    assert seqs and sum(len(s) for s in seqs.values()) > 1_000_000
