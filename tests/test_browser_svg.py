"""El SVG del navegador debe ser XML válido, escapar etiquetas y marcar los
elementos que se pueden pulsar."""
from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd

from browser import svg


def _tables() -> dict:
    genes = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 400, "strand": "1",
         "gene_id": "g1", "locus_tag": "b0001", "gene_name": "thrA",
         "product": "aspartokinase", "ec": "2.7.2.4"},
        {"contig": "c1", "start": 500, "end": 800, "strand": "-1",
         "gene_id": "g2", "locus_tag": "b0002", "gene_name": "a<b>c",
         "product": "producto con <script>alert(1)</script>", "ec": ""},
    ])
    operons = pd.DataFrame([
        {"contig": "c1", "start": 100, "end": 800, "operon_id": "op1",
         "label": "thr", "n_genes": 2},
    ])
    promoters = pd.DataFrame([
        {"operon_id": "op1", "label": "thr", "box35": "TTGACA",
         "box10": "TATAAT", "sigma": "70", "score": "0.9"},
    ])
    hairpins = pd.DataFrame([
        {"contig": "c1", "start": 820, "end": 860, "class": "A_fuerte",
         "dG": "-12.3", "stem": "8", "gc": "0.6"},
    ])
    return {"genes": genes, "operons": operons, "promoters": promoters,
            "hairpins": hairpins}


def test_svg_is_valid_xml_and_escapes_labels():
    out = svg(_tables(), "c1", 0, 1000)
    ET.fromstring(out)          # lanza si el XML está roto
    assert "<script>" not in out
    assert "&lt;b&gt;" in out


def test_clickable_elements_have_kind_and_id():
    out = svg(_tables(), "c1", 0, 1000)
    for kind in ("gene", "operon", "promoter", "hairpin"):
        assert f'data-kind="{kind}"' in out, kind
    assert 'data-id="b0001"' in out
    assert 'data-id="op1"' in out


def test_empty_tables_do_not_crash():
    empty = {k: pd.DataFrame() for k in ("genes", "operons", "promoters", "hairpins")}
    out = svg(empty, "c1", 0, 1000)
    ET.fromstring(out)
