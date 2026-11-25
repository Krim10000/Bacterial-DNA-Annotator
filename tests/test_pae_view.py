from webapp.pae_view import heatmap_svg


def test_heatmap_svg_valid():
    svg = heatmap_svg([[1, 2], [3, 4]])
    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert "PAE" in svg and "2&#215;2" not in svg
    import xml.etree.ElementTree as ET
    ET.fromstring(svg)


def test_heatmap_empty():
    svg = heatmap_svg([])
    assert "no disponible" in svg
    import xml.etree.ElementTree as ET
    ET.fromstring(svg)
