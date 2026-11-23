"""Complejos por vecindad: no unir genes lejanos, de otra hebra o con RNA."""
from __future__ import annotations

from webapp.complexes import enrich_gaps, enrich_rows, infer, neighborhoods


def _feat(pid, gene, start, end=500, strand="+", contig="c", typ="CDS"):
    return {"id": pid, "type": typ, "contig": contig, "strand": strand,
            "start": start, "end": start + end, "gene": gene}


def _hits(feats):
    return {f["id"]: {"gene": f["gene"], "ec": ""} for f in feats}


def test_complete_and_incomplete_atp():
    feats = [_feat("a", "atpA", 100), _feat("b", "atpB", 700),
             _feat("c", "atpC", 1300), _feat("d", "atpD", 1900),
             _feat("e", "atpE", 2500)]
    out = infer(feats, _hits(feats), "x")
    assert out["n_complete"] == 1
    assert out["complexes"][0]["status"] == "complex_complete"
    assert out["complexes"][0]["name"].startswith("ATP synthase")
    assert "7.1.2.2" in out["complexes"][0]["ecs"]
    feats2 = [f for f in feats if f["id"] != "d"]
    out2 = infer(feats2, _hits(feats2), "x")
    assert out2["complexes"][0]["status"] == "complex_incomplete"
    assert "atpd" in out2["complexes"][0]["missing"]


def test_rna_breaks_neighborhood():
    feats = [_feat("a", "atpA", 100), _feat("r", "trna", 700, typ="tRNA"),
             _feat("b", "atpB", 900), _feat("c", "atpC", 1500),
             _feat("d", "atpD", 2100)]
    groups = neighborhoods(feats)
    assert len(groups) == 2
    out = infer(feats, _hits(feats), "x")
    assert out["complexes"][0]["status"] == "complex_incomplete"
    assert "atpa" in out["complexes"][0]["missing"]


def test_different_strand_or_contig_not_mixed():
    feats = [_feat("a", "atpA", 100), _feat("b", "atpB", 700, strand="-"),
             _feat("c", "atpC", 1300, contig="other")]
    assert infer(feats, _hits(feats), "x")["complexes"] == []


def test_enrich_rows_and_gaps():
    feats = [_feat("a", "atpA", 100), _feat("b", "atpB", 700),
             _feat("c", "atpC", 1300), _feat("d", "atpD", 1900)]
    hits = {f["id"]: {"gene": f["gene"], "ec": "7.1.2.2"} for f in feats}
    out = infer(feats, hits, "x")
    rows = [{"protein_id": "a"}]
    enrich_rows(rows, out)
    assert rows[0]["complex_id"].startswith("ATP_SYNTHASE:")
    assert rows[0]["complex_role"] == "subunit"
    assert rows[0]["complex_flag"] in ("green", "yellow")
    gaps = [{"producers": "RHEA:1"}]
    enrich_gaps(gaps, {**out, "protein_rows": [{"complex_id": rows[0]["complex_id"],
                                                "ec": "7.1.2.2"}]}, None)
    assert gaps[0]["complex_status"] == "not_evaluated"  # sin rhea2ec no se inventa
