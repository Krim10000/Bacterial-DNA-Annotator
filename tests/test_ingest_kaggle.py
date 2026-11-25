import csv
import zipfile
from pathlib import Path

from scripts.structure.ingest_kaggle_zips import ingest


def _zip(path: Path, name: str, acc: str, plddt: float):
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(f"{name}/metrics.tsv",
                   f"uniprot\tstatus\tlength\tplddt\tptm\tseconds\terrorr\n"
                   f"{acc}\tdone\t100\t{plddt}\t0.8\t1\t\n")
        z.writestr(f"{name}/pdb/{acc}.pdb",
                   "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")


def test_ingest_dedup_keeps_best(tmp_path):
    z1 = tmp_path / "a.zip"
    z2 = tmp_path / "b.zip"
    _zip(z1, "a", "P1", 0.5)
    _zip(z2, "b", "P1", 0.9)
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"
    r = ingest([z1], out, cat)
    assert r["added"] == 1
    r2 = ingest([z2], out, cat)
    assert r2["improved"] == 1
    rows = list(csv.DictReader(open(cat), delimiter="\t"))
    assert rows[0]["plddt"] == "0.900"
    assert rows[0]["source_zip"] == "b.zip"


def test_ingest_dry_run_writes_nothing(tmp_path):
    z = tmp_path / "a.zip"
    _zip(z, "a", "P2", 0.7)
    out = tmp_path / "out"
    cat = tmp_path / "catalog.tsv"
    r = ingest([z], out, cat, dry_run=True)
    assert r["added"] == 1
    assert not cat.exists()
    assert not (out / "P2.pdb").exists()
