"""Carga de transferencias K-12 como sospecha fuerte."""
from __future__ import annotations

import csv
from pathlib import Path

from webapp.complex_transfer import load


def test_transfer_strong_suspicion(tmp_path: Path):
    reg = tmp_path / "regulation" / "g"
    reg.mkdir(parents=True)
    with (reg / "genes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["gene_id", "contig", "start", "end", "strand"])
        w.writeheader()
        w.writerow({"gene_id": "g1", "contig": "c", "start": 1, "end": 100, "strand": "+"})
        w.writerow({"gene_id": "g2", "contig": "c", "start": 200, "end": 300, "strand": "+"})
    cx = tmp_path / "complexes"
    cx.mkdir()
    (cx / "transfer_g.tsv").write_text(
        "complex_id\tcomplex_name\tn_mapped\tfrac_mapped\tidentity_mean\t"
        "genes_genoma\texpected_n\tecs\n"
        "CPX-1\tTest complex\t2\t1.0\t80\tg1,g2\t2\t1.2.3.4\n", encoding="utf-8")
    out = load(tmp_path, "g", [])
    assert out["ok"]
    assert out["complexes"][0]["status"] == "complex_suspected"
    assert out["complexes"][0]["flag"] == "yellow"
    assert out["by_protein"]["g1"]["complex_id"] == "K12:CPX-1"
