"""Pruebas de i18n de la comparación de genomas y del visor 3D.

Cubren las cadenas visibles nuevas de ``webapp/compare_ui.py`` y
``webapp/structure_viewer.py``: presencia en ``_source.json`` y en los
catálogos EN/ZH, placeholders y markdown intactos, nombres técnicos sin
traducir y traducción efectiva en tiempo de ejecución.

Los nombres de formatos, estados y campos técnicos (``PDB``, ``GFF3``,
``TSV``, ``ESMFold``/``AFDB``/``AlphaFold``, ``locus_tag``, ``ID``,
``Parent``) se quedan literales a propósito.
"""
from __future__ import annotations

import io
import json
import re
from pathlib import Path
from typing import Any

import pytest

import compare_ui
import i18n
import structure_viewer

ROOT = Path(__file__).resolve().parents[1]
LOCALES = ROOT / "webapp" / "locales"

#: Cadenas nuevas de esta ronda; el español es la clave del catálogo.
NEW_KEYS = {
    "Comparar GFF3",
    "Genes {genoma}",
    "⬇️ Descargar PDB",
    "Genes compartidos",
    "GFF3 del genoma B",
    "GFF3 del genoma A",
    "Advertencias ({n})",
    "Exclusivos {genoma}",
    "Comparar dos genomas",
    "Etiqueta del genoma A",
    "Etiqueta del genoma B",
    "Comparando anotaciones…",
    "No se pudo comparar: {error}",
    "Qué no prueba esta comparación",
    "Visor 3D no disponible: {error}",
    "… o ruta local al GFF3 del genoma B",
    "… o ruta local al GFF3 del genoma A",
    "⬇️ Descargar tabla de diferencias (TSV)",
    "El archivo no parece un PDB con coordenadas.",
    "Las etiquetas de los genomas deben ser distintas.",
    "Selecciona un GFF3 para A y otro para B, como archivo o como ruta.",
    "Coincidencias por método: locus_tag={locus_tag}, ID={id}, "
    "coordenadas={coordenadas}.",
    "No se calcula identidad de secuencia ni similitud. Un gen compartido "
    "por anotación no es necesariamente un ortólogo.",
    "Si los contigs tienen nombres distintos o los genomas no están "
    "alineados, las coordenadas no son comparables entre sí.",
    "Un gen exclusivo puede ser un gen no anotado, un fragmento, un cambio "
    "de coordenadas o una diferencia real entre los genomas.",
    "Vista 3D de coordenadas. Una estructura predicha (ESMFold/AFDB/"
    "AlphaFold) no confirma función, sustrato, interfaz ni cierre de brecha.",
    "La comparación usa solo la anotación GFF3 (locus_tag, ID, Parent y "
    "coordenadas); no lee ni alinea ADN o proteína y no consulta la red.",
    "Cuando compares, verás cuántos genes tiene cada genoma, cuáles "
    "comparten locus_tag/ID/coordenadas y cuáles son exclusivos. La tabla "
    "se puede descargar como TSV.",
    "Sube un GFF3 por genoma o escribe su ruta local. Si hay archivo "
    "subido, se usa ese; la ruta se usa solo cuando no hay archivo. "
    "Formatos: .gff3, .gff, .txt y .gz.",
    "Dos genes con las mismas coordenadas pueden ser distintos; dos "
    "ortólogos pueden tener coordenadas distintas. La coincidencia "
    "posicional no prueba ortología ni función.",
    "Comparar coordenadas **no es ortología ni función**: dos genes en la "
    "misma posición pueden ser genes distintos, y dos ortólogos pueden "
    "estar en coordenadas distintas. Sin secuencia no se calcula identidad; "
    "esta tabla compara anotaciones, no biología.",
}

#: Tokens que deben sobrevivir a la traducción tal cual.
TECHNICAL_TOKENS = (
    "PDB", "GFF3", "TSV", "ESMFold", "AFDB", "AlphaFold",
    "locus_tag", "ID", "Parent",
)


def _json(name: str) -> Any:
    return json.loads((LOCALES / name).read_text(encoding="utf-8"))


def _markdown_links(text: str) -> list[str]:
    return re.findall(r"\[[^\]]+\]\([^)]+\)", text)


@pytest.fixture(autouse=True)
def _reset_language():
    i18n.set_language("es")
    yield
    i18n.set_language("es")


def test_new_keys_are_in_source_and_both_catalogs():
    source = set(_json("_source.json"))
    assert NEW_KEYS <= source, sorted(NEW_KEYS - source)
    for lang in ("en", "zh"):
        catalog = _json(f"{lang}.json")
        missing = NEW_KEYS - set(catalog)
        assert not missing, f"{lang} sin traducir: {sorted(missing)}"
        assert not [k for k in NEW_KEYS if not str(catalog[k]).strip()]


def test_catalogs_cover_every_source_string():
    source = set(_json("_source.json"))
    for lang in ("en", "zh"):
        missing = source - set(_json(f"{lang}.json"))
        assert not missing, f"{lang} sin traducir: {sorted(missing)[:10]}"


def test_new_translations_keep_placeholders_and_markdown():
    for lang in ("en", "zh"):
        catalog = _json(f"{lang}.json")
        for key in NEW_KEYS:
            value = catalog[key]
            assert sorted(re.findall(r"\{\w+\}", key)) == sorted(
                re.findall(r"\{\w+\}", value)
            ), f"{lang}: placeholders en {key!r}"
            assert key.count("`") == value.count("`"), f"{lang}: backticks en {key!r}"
            assert key.count("**") == value.count("**"), f"{lang}: negritas en {key!r}"
            assert _markdown_links(key) == _markdown_links(value), (
                f"{lang}: enlaces markdown en {key!r}"
            )


def test_technical_names_stay_literal_in_new_translations():
    for lang in ("en", "zh"):
        catalog = _json(f"{lang}.json")
        for key in NEW_KEYS:
            value = catalog[key]
            for token in TECHNICAL_TOKENS:
                if token in key:
                    assert token in value, (
                        f"{lang}: {token!r} se tradujo en {key!r} -> {value!r}"
                    )


def test_tr_translates_new_keys_and_falls_back_to_spanish():
    key = "Genes compartidos"
    assert i18n.tr(key) == key  # español = clave
    i18n.set_language("en")
    assert i18n.tr(key) == "Shared genes"
    assert i18n.tr("Genes {genoma}").format(genoma="A") == "Genes A"
    i18n.set_language("zh")
    assert i18n.tr(key) == "共享基因"
    assert i18n.tr("Genes {genoma}").format(genoma="B") == "基因 B"
    assert i18n.tr("cadena que no existe") == "cadena que no existe"


# ------------------------------------------------------------------ runtime
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
    def __init__(self) -> None:
        self.calls: list = []

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "columns":
                count = int(args[0]) if args else 1
                return [_FakeCtx(self.calls) for _ in range(count)]
            return _FakeCtx(self.calls)

        return recorder


def _gff(start: int, attrs: str, seqid: str = "chr1") -> str:
    return "\t".join(
        [seqid, "test", "CDS", str(start), str(start + 99), ".", "+", "0", attrs]
    )


def _compare_result():
    a = io.StringIO(
        "##gff-version 3\n" + _gff(10, "ID=cds_a;locus_tag=LT_0001") + "\n"
    )
    # Otro contig: añade una advertencia visible sin cambiar el emparejamiento
    # por locus_tag, así se cubre la etiqueta del expander de advertencias.
    b = io.StringIO(
        "##gff-version 3\n"
        + _gff(10, "ID=cds_b;locus_tag=LT_0001", seqid="chr2")
        + "\n"
    )
    return compare_ui.compare_gff(a, b, label_a="K12", label_b="EC")


def test_compare_result_renders_catalog_text_without_touching_data():
    result = _compare_result()
    fake = _FakeStreamlit()
    i18n.set_language("en")

    compare_ui._render_result(fake, result)

    metrics = [call[1][0] for call in fake.calls if call[0] == "metric"]
    assert "Genes K12" in metrics
    assert "Shared genes" in metrics
    assert "Exclusive to EC" in metrics
    captions = [call[1][0] for call in fake.calls if call[0] == "caption"]
    assert any("Matches by method: locus_tag=1, ID=0" in c for c in captions)
    downloads = [call for call in fake.calls if call[0] == "download_button"]
    assert downloads[0][1][0] == i18n.tr(
        "⬇️ Descargar tabla de diferencias (TSV)"
    )
    expanders = [call[1][0] for call in fake.calls if call[0] == "expander"]
    assert any(label.startswith("Warnings (") for label in expanders)
    assert "What this comparison does not prove" in expanders
    writes = [call[1][0] for call in fake.calls if call[0] == "write"]
    assert writes[0] == "• " + i18n.tr(compare_ui.LIMITS[0])
    # Los datos de la tabla (estados, notas) no se traducen al renderizar.
    frames = [call[1][0] for call in fake.calls if call[0] == "dataframe"]
    assert set(frames[0]["estado"]) == {"compartido"}


def test_structure_viewer_messages_use_catalog(monkeypatch):
    calls: list = []

    class _St:
        def __getattr__(self, name):
            def recorder(*args, **kwargs):
                calls.append((name, args, kwargs))
                return None

            return recorder

    monkeypatch.setattr(structure_viewer, "st", _St())
    monkeypatch.setattr(structure_viewer, "available", lambda: (True, ""))

    i18n.set_language("en")
    structure_viewer.render("no coordinates here")

    infos = [call[1][0] for call in calls if call[0] == "info"]
    assert infos == [i18n.tr("El archivo no parece un PDB con coordenadas.")]
    assert "PDB" in infos[0]

    calls.clear()
    monkeypatch.setattr(
        structure_viewer, "available", lambda: (False, "py3Dmol missing")
    )
    structure_viewer.render("whatever")
    infos = [call[1][0] for call in calls if call[0] == "info"]
    assert infos == [
        i18n.tr("Visor 3D no disponible: {error}").format(
            error="py3Dmol missing"
        )
    ]
    assert "3D viewer unavailable: py3Dmol missing" == infos[0]
