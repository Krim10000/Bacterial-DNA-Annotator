"""Tests sintéticos de la llamada de rRNA (``webapp/rrna_calling.py``).

No usan red, GPU ni binarios externos: solo parsers, deduplicación y la
orquestación con ``run_cmsearch``/``run_barrnap`` sustituidos. Los casos
cubren 16S/23S/5S, hebra negativa, coordenadas invertidas, duplicados entre
barrnap y cmsearch, tipos/hebras distintos y el GFF de reanotación.
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
REF_DIR = ROOT / "scripts" / "reference"
if str(REF_DIR) not in sys.path:
    sys.path.insert(0, str(REF_DIR))

from webapp import rrna_calling as rc  # noqa: E402
import reannotate_rrna as rr  # noqa: E402

# tblout de cmsearch: columnas 0, 2, 3, 7, 8, 9, 13, 15.
CM_TBLOUT = """\
# cmsearch :: INFERNAL 1.1.4
# target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc score bias E-value inc description
NC_1 - SSU_rRNA_bacteria RF00177 - - - 100 1600 + no 1 0.55 21.9 0.0 0.0 ! -
NC_1 - LSU_rRNA_bacteria RF02541 - - - 1700 4600 + no 1 0.53 70.2 0.0 0.0 ! -
NC_1 - 5S_rRNA RF00001 - - - 4700 4819 - no 1 0.61 12.0 0.0 1e-10 ! -
NC_1 - Unknown_model RF99999 - - - bad 5000 + no 1 0.5 1.0 0.0 2.0 ! -
"""

BARRNAP_GFF = """\
##gff-version 3
NC_1\tbarrnap:0.9\trRNA\t101\t1598\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
NC_1\tbarrnap:0.9\trRNA\t4700\t4810\t4e-12\t-\t.\tName=5S_rRNA;product=5S ribosomal RNA
NC_1\tbarrnap:0.9\tgene\t10\t20\t.\t+\t.\tID=gene1
"""


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def test_parse_cmsearch_maps_models_and_types() -> None:
    features = rc.parse_cmsearch_tblout(CM_TBLOUT)
    assert len(features) == 3  # la fila con coordenada inválida se ignora
    by_type = {rc.rrna_class(f): f for f in features}
    assert set(by_type) == {"16S", "23S", "5S"}
    assert by_type["16S"]["product"] == "16S ribosomal RNA"
    assert by_type["16S"]["rfam"] == "RF00177"
    assert by_type["23S"]["rfam"] == "RF02541"
    assert by_type["5S"]["source"] == "cmsearch"
    # Hebra negativa con seq from > seq to: se normaliza y conserva la hebra.
    assert by_type["5S"]["start"] == 4700
    assert by_type["5S"]["end"] == 4819
    assert by_type["5S"]["strand"] == "-"
    assert by_type["5S"]["score"] == 12.0
    assert by_type["5S"]["evalue"] == 1e-10


def test_parse_barrnap_ignores_non_rrna_and_reads_product() -> None:
    features = rc.parse_barrnap_gff(BARRNAP_GFF)
    assert len(features) == 2
    assert features[0]["start"] == 101
    assert features[0]["end"] == 1598
    assert features[0]["product"] == "16S ribosomal RNA"
    assert features[1]["strand"] == "-"
    assert features[1]["model"] == "5S_rRNA"
    assert rc.rrna_class(features[1]) == "5S"


def test_rrna_class_uses_model_rfam_or_product() -> None:
    assert rc.rrna_class({"model": "LSU_rRNA_bacteria"}) == "23S"
    assert rc.rrna_class({"rfam": "RF02541"}) == "23S"
    assert rc.rrna_class({"product": "5S ribosomal RNA"}) == "5S"
    assert rc.rrna_class({"product": "rnase P"}) == "other"


# ---------------------------------------------------------------------------
# Deduplicación
# ---------------------------------------------------------------------------
def test_dedupe_prefers_cmsearch_over_barrnap() -> None:
    cm = rc.parse_cmsearch_tblout(CM_TBLOUT)
    barrnap = rc.parse_barrnap_gff(BARRNAP_GFF)
    merged = rc.dedupe_rrna(cm + barrnap)
    # 16S y 5S de barrnap se solapan con los de cmsearch: se quedan los CM.
    assert len(merged) == 3
    assert all(feature["source"] == "cmsearch" for feature in merged)
    assert {rc.rrna_class(f) for f in merged} == {"16S", "23S", "5S"}


def test_dedupe_does_not_collapse_different_types_or_strands() -> None:
    features = [
        {"contig": "c", "start": 100, "end": 200, "strand": "+",
         "product": "16S ribosomal RNA", "source": "barrnap"},
        {"contig": "c", "start": 100, "end": 200, "strand": "+",
         "product": "23S ribosomal RNA", "source": "barrnap"},
        {"contig": "c", "start": 100, "end": 200, "strand": "-",
         "product": "16S ribosomal RNA", "source": "barrnap"},
    ]
    assert len(rc.dedupe_rrna(features)) == 3


def test_dedupe_removes_overlapping_duplicate_of_same_type() -> None:
    features = [
        {"contig": "c", "start": 100, "end": 200, "strand": "+",
         "product": "16S ribosomal RNA", "source": "cmsearch", "score": 10.0},
        {"contig": "c", "start": 105, "end": 198, "strand": "+",
         "product": "16S ribosomal RNA", "source": "cmsearch", "score": 5.0},
    ]
    kept = rc.dedupe_rrna(features)
    assert len(kept) == 1
    assert kept[0]["start"] == 100


def test_assign_ids_orders_by_genome_position() -> None:
    features = [
        {"contig": "c", "start": 500, "end": 600, "strand": "+",
         "product": "16S ribosomal RNA"},
        {"contig": "c", "start": 100, "end": 200, "strand": "+",
         "product": "5S ribosomal RNA"},
    ]
    ordered = rc.assign_ids(features)
    assert [f["id"] for f in ordered] == ["rrna_1", "rrna_2"]
    assert ordered[0]["start"] == 100


def test_format_gff_line_is_one_based_inclusive() -> None:
    feature = {
        "contig": "c", "start": 10, "end": 20, "strand": "-",
        "product": "23S ribosomal RNA", "rfam": "RF02541",
        "source": "cmsearch", "id": "rrna_1",
    }
    line = rc.format_gff_line(feature)
    cols = line.rstrip("\n").split("\t")
    assert cols[0] == "c"
    assert cols[2] == "rRNA"
    assert cols[3] == "10"
    assert cols[4] == "20"
    assert cols[6] == "-"
    assert "Dbxref=Rfam:RF02541" in line


# ---------------------------------------------------------------------------
# Orquestación con herramientas simuladas
# ---------------------------------------------------------------------------
def test_call_prefers_cmsearch(monkeypatch, tmp_path: Path) -> None:
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    cm = tmp_path / "models.cm"
    cm.write_text("# fake\n", encoding="utf-8")
    calls = {"cm": 0, "barrnap": 0}

    def fake_cm(*_args, **_kwargs):
        calls["cm"] += 1
        return [{"contig": "c", "start": 1, "end": 100, "strand": "+",
                 "product": "16S ribosomal RNA", "source": "cmsearch"}]

    def fake_barrnap(*_args, **_kwargs):
        calls["barrnap"] += 1
        return []

    monkeypatch.setattr(rc.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rc, "run_cmsearch", fake_cm)
    monkeypatch.setattr(rc, "run_barrnap", fake_barrnap)
    result = rc.call_rrna_fasta(fasta, cm_path=cm, cpu=1)
    assert result["tool"] == "cmsearch"
    assert len(result["features"]) == 1
    assert calls == {"cm": 1, "barrnap": 0}


def test_call_falls_back_to_barrnap(monkeypatch, tmp_path: Path) -> None:
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    cm = tmp_path / "models.cm"
    cm.write_text("# fake\n", encoding="utf-8")

    def failing_cm(*_args, **_kwargs):
        raise OSError("cmsearch ausente")

    def fake_barrnap(*_args, **_kwargs):
        return [{"contig": "c", "start": 5, "end": 10, "strand": "-",
                 "product": "5S ribosomal RNA", "source": "barrnap"}]

    monkeypatch.setattr(rc.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rc, "run_cmsearch", failing_cm)
    monkeypatch.setattr(rc, "run_barrnap", fake_barrnap)
    result = rc.call_rrna_fasta(fasta, cm_path=cm, cpu=1)
    assert result["tool"] == "barrnap"
    assert len(result["features"]) == 1
    assert "cmsearch" in result["reason"]


def test_extract_cm_models_keeps_only_requested_cm_records(tmp_path: Path) -> None:
    cm = tmp_path / "panel.cm"
    cm.write_text(
        "INFERNAL1/a\nNAME  SSU_rRNA_bacteria\nACC RF00177\n//\n"
        "HMMER3/f\nNAME SSU_rRNA_bacteria\n//\n"
        "INFERNAL1/a\nNAME  LSU_rRNA_bacteria\nACC RF02541\n//\n"
        "INFERNAL1/a\nNAME  5S_rRNA\nACC RF00001\n//\n",
        encoding="utf-8",
    )
    text = rc.extract_cm_models(cm.read_text(), ["SSU_rRNA_bacteria", "5S_rRNA"])
    assert text.count("INFERNAL1/a") == 2
    # El bloque HMMER3 asociado se conserva como filtro p7 del CM.
    assert text.count("HMMER3/f") == 1
    assert "LSU_rRNA_bacteria" not in text
    assert "5S_rRNA" in text


def test_balanced_strategy_uses_cm_for_5s_and_hmmonly_for_lsu(
    monkeypatch, tmp_path: Path
) -> None:
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    cm = tmp_path / "panel.cm"
    cm.write_text(
        "INFERNAL1/a\nNAME  SSU_rRNA_bacteria\n//\n"
        "INFERNAL1/a\nNAME  LSU_rRNA_bacteria\n//\n"
        "INFERNAL1/a\nNAME  5S_rRNA\n//\n",
        encoding="utf-8",
    )
    seen = []

    def fake_models(_fasta, _cm, models, **kwargs):
        seen.append((tuple(models), kwargs.get("hmmonly", False)))
        return []

    monkeypatch.setattr(rc.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rc, "run_cmsearch_models", fake_models)
    rc.call_rrna_fasta(fasta, cm_path=cm, strategy="balanced")
    assert seen == [
        (("SSU_rRNA_bacteria", "5S_rRNA"), False),
        (("LSU_rRNA_bacteria",), True),
    ]


def test_partial_cmsearch_failure_keeps_successful_hits(
    monkeypatch, tmp_path: Path
) -> None:
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    cm = tmp_path / "panel.cm"
    cm.write_text("# fake\n", encoding="utf-8")
    calls = {"n": 0}

    def fake_models(*_args, **_kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            return [{"contig": "c", "start": 1, "end": 10, "strand": "+",
                     "product": "5S ribosomal RNA", "source": "cmsearch"}]
        raise OSError("falló el segundo modelo")

    monkeypatch.setattr(rc.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rc, "run_cmsearch_models", fake_models)
    result = rc.call_rrna_fasta(fasta, cm_path=cm, strategy="balanced")
    assert result["tool"] == "cmsearch"
    assert len(result["features"]) == 1
    assert "cmsearch falló" in result["reason"]


def test_annotator_run_rrna_uses_recall_strategy(monkeypatch) -> None:
    import webapp.annotator as ann

    captured = {}

    def fake_call(_fasta, **kwargs):
        captured.update(kwargs)
        return {"features": [{"contig": "c", "start": 1, "end": 10,
                              "strand": "+", "product": "5S ribosomal RNA"}]}

    monkeypatch.setattr(rc, "call_rrna_fasta", fake_call)
    features = ann.run_rrna({"c": "ACGT"})
    assert len(features) == 1
    assert captured["strategy"] == "recall"


# ---------------------------------------------------------------------------
# Reescritura del GFF del sistema
# ---------------------------------------------------------------------------
def _write_gff(path: Path, rows: list[str]) -> None:
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")


def test_rewrite_gff_replaces_only_rrna(tmp_path: Path) -> None:
    gff = tmp_path / "g.gff3"
    _write_gff(
        gff,
        [
            "c\tsystem\tCDS\t10\t99\t.\t+\t0\tID=cds_1",
            "c\tbarrnap\trRNA\t200\t300\t.\t+\t.\tID=rrna_1;product=16S ribosomal RNA",
            "c\tsystem\ttRNA\t400\t470\t.\t-\t.\tID=trna_1",
        ],
    )
    new = [
        {"contig": "c", "start": 195, "end": 310, "strand": "+",
         "product": "23S ribosomal RNA", "rfam": "RF02541",
         "source": "cmsearch", "id": "rrna_1"},
    ]
    rr._rewrite_gff(gff, new)
    text = gff.read_text(encoding="utf-8")
    assert "ID=cds_1" in text
    assert "ID=trna_1" in text
    assert "16S ribosomal RNA" not in text
    assert "23S ribosomal RNA" in text
    assert "ID=rrna_1" in text


def test_raw_rrna_rows_classifies_rfam(tmp_path: Path) -> None:
    gff = tmp_path / "g.gff3"
    _write_gff(
        gff,
        [
            "c\tx\trRNA\t1\t100\t.\t+\t.\tID=a;Dbxref=Rfam:RF00177",
            "c\tx\trRNA\t200\t300\t.\t-\t.\tID=b;Dbxref=Rfam:RF02541",
            "c\tx\trRNA\t400\t500\t.\t+\t.\tID=c;Dbxref=Rfam:RF00001",
        ],
    )
    rows = rr._raw_rrna_rows(gff)
    counts = rr._count_by_type(rows)
    assert counts == {"16S": 1, "23S": 1, "5S": 1}


def test_update_summary_only_changes_rrna(tmp_path: Path) -> None:
    summary = tmp_path / "g.summary.json"
    summary.write_text(
        json.dumps({"features": {"CDS": 10, "tRNA": 2, "rRNA": 3, "ncRNA": 1}}),
        encoding="utf-8",
    )
    rr._update_summary(summary, 7)
    data = json.loads(summary.read_text(encoding="utf-8"))
    assert data["features"] == {"CDS": 10, "tRNA": 2, "rRNA": 7, "ncRNA": 1}
