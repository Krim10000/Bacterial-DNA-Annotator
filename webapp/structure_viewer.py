"""Visor 3D de PDB para Streamlit (py3Dmol/stmol).

Solo dibuja coordenadas. Una estructura predicha no confirma función ni cierra
una brecha; el visor muestra PDB ya generados o aportados por el usuario.
"""
from __future__ import annotations

import gzip
import html as _html
import json
import re
from pathlib import Path
from typing import Any

import streamlit as st

try:  # importación normal dentro de la app (webapp/ en sys.path)
    from i18n import tr  # type: ignore
except ImportError:  # importación como paquete: webapp.structure_viewer
    try:
        from .i18n import tr  # type: ignore
    except ImportError:  # pragma: no cover - permite usar el módulo aislado
        def tr(text: Any) -> Any:  # type: ignore
            return text

try:
    import py3Dmol
    import stmol
    _AVAILABLE = True
    _WHY = ""
except Exception as exc:  # pragma: no cover - depende del entorno
    _AVAILABLE = False
    _WHY = f"{type(exc).__name__}: {exc}"


def available() -> tuple[bool, str]:
    return _AVAILABLE, _WHY


def _read_pdb(pdb: str | Path) -> str:
    """Lee coordenadas de un archivo local, descomprimiendo ``.gz``.

    Acepta PDB/mmCIF planos o comprimidos con gzip. Si ``pdb`` no es un
    archivo existente se devuelve tal cual: el llamador puede pasar el texto
    de las coordenadas directamente.
    """
    p = Path(pdb)
    if not p.is_file():
        return str(pdb)
    if p.suffix.lower() == ".gz":
        try:
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                return fh.read()
        except (OSError, EOFError):
            # Un .gz corrupto se lee como texto: el visor avisará que no hay
            # coordenadas en vez de romper la ficha.
            return p.read_text(encoding="utf-8", errors="replace")
    return p.read_text(encoding="utf-8", errors="replace")


def _structure_format(pdb: str | Path, text: str = "") -> str:
    """Devuelve ``"cif"`` o ``"pdb"`` según el nombre o el contenido."""
    name = str(pdb).lower()
    if name.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        return "cif"
    if text.lstrip().startswith("data_"):
        return "cif"
    return "pdb"


def _interactive_html(view: Any, height: int, sites: list[int] | None = None) -> str:
    """Envuelve el HTML de py3Dmol con controles de zoom y pantalla completa."""
    raw = view._make_html()
    vid = f"viewer_{view.uniqueid}"
    buttons = (
        '<div style="display:flex;gap:8px;align-items:center;'
        'font-family:system-ui,sans-serif;margin:0 0 6px 0">'
        f'<button onclick="{vid}.zoom(1.2);{vid}.render();">🔍 +</button>'
        f'<button onclick="{vid}.zoom(1/1.2);{vid}.render();">🔍 −</button>'
        f'<button onclick="{vid}.zoomTo();{vid}.render();">Ajustar</button>'
        f'<button onclick="document.documentElement.requestFullscreen()">⛶ Pantalla completa</button>'
        '<span style="font-size:11px;color:#64748b">rueda del ratón / pellizco</span>'
        '</div>'
    )
    if sites:
        raw += (
            f'<script>$3Dmolpromise.then(function(){{{vid}.addStyle('
            f'{{resi:{json.dumps(sites)}}},'
            f'{{stick:{{colorscheme:"magenta",radius:0.28}}}});'
            f'{vid}.render();}});</script>'
        )
    # El contenedor ocupa el ancho disponible; el visor mantiene su tamaño
    # interno para que la cámara no quede diminuta.
    raw = re.sub(
        r'style="position: relative; width: \d+px; height: \d+px;"',
        f'style="position: relative; width:100%; height:{height}px;"',
        raw, count=1,
    )
    return buttons + raw


def render(pdb: str | Path, key: str = "pdb_viewer", height: int = 440,
           caption: str = "", bfactor: bool = False,
           interactive: bool = False, width: int = 1100,
           sites: list[int] | None = None) -> None:
    """Dibuja un PDB con cartoon y heteroátomos en stick."""
    ok, why = available()
    if not ok:
        st.info(tr("Visor 3D no disponible: {error}").format(error=why))
        return
    text = _read_pdb(pdb)
    if "ATOM" not in text and "HETATM" not in text:
        st.info(tr("El archivo no parece un PDB con coordenadas."))
        return
    fmt = _structure_format(pdb, text)
    view = py3Dmol.view(width=width, height=height)
    view.addModel(text, fmt)
    if bfactor:
        view.setStyle({"cartoon": {"colorscheme": {
            "prop": "b", "gradient": "roygb"}}})
    else:
        view.setStyle({"cartoon": {"color": "spectrum"}})
    view.addStyle({"hetflag": True}, {"stick": {"radius": 0.18}})
    view.setBackgroundColor("white")
    view.zoomTo()
    if interactive:
        try:
            import streamlit.components.v1 as components
            components.html(_interactive_html(view, height, sites),
                            height=height + 48, scrolling=False)
        except Exception:
            stmol.showmol(view, height=height, width=width)
    else:
        stmol.showmol(view, height=height, width=width)
    if caption:
        st.caption(caption)
    st.caption(tr(
        "Vista 3D de coordenadas. Una estructura predicha (ESMFold/AFDB/"
        "AlphaFold) no confirma función, sustrato, interfaz ni cierre de "
        "brecha."
    ))
    st.download_button(
        tr("⬇️ Descargar PDB"), text, file_name=f"{key}.{fmt}",
        mime="chemical/x-cif" if fmt == "cif" else "chemical/x-pdb",
        key=f"{key}_pdb_dl",
    )
