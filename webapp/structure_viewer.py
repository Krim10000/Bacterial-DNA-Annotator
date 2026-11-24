"""Visor 3D de PDB para Streamlit (py3Dmol/stmol).

Solo dibuja coordenadas. Una estructura predicha no confirma función ni cierra
una brecha; el visor muestra PDB ya generados o aportados por el usuario.
"""
from __future__ import annotations

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
    p = Path(pdb)
    if p.is_file():
        return p.read_text(encoding="utf-8", errors="replace")
    return str(pdb)


def render(pdb: str | Path, key: str = "pdb_viewer", height: int = 440,
           caption: str = "", bfactor: bool = False) -> None:
    """Dibuja un PDB con cartoon y heteroátomos en stick."""
    ok, why = available()
    if not ok:
        st.info(tr("Visor 3D no disponible: {error}").format(error=why))
        return
    text = _read_pdb(pdb)
    if "ATOM" not in text and "HETATM" not in text:
        st.info(tr("El archivo no parece un PDB con coordenadas."))
        return
    fmt = "cif" if (str(pdb).lower().endswith((".cif", ".mmcif"))
                     or text.lstrip().startswith("data_")) else "pdb"
    view = py3Dmol.view(width=760, height=height)
    view.addModel(text, fmt)
    if bfactor:
        view.setStyle({"cartoon": {"colorscheme": {
            "prop": "b", "gradient": "roygb"}}})
    else:
        view.setStyle({"cartoon": {"color": "spectrum"}})
    view.addStyle({"hetflag": True}, {"stick": {"radius": 0.18}})
    view.setBackgroundColor("white")
    view.zoomTo()
    stmol.showmol(view, height=height)
    if caption:
        st.caption(caption)
    st.caption(tr(
        "Vista 3D de coordenadas. Una estructura predicha (ESMFold/AFDB/"
        "AlphaFold) no confirma función, sustrato, interfaz ni cierre de "
        "brecha."
    ))
    st.download_button(
        tr("⬇️ Descargar PDB"), text, file_name=f"{key}.pdb",
        mime="chemical/x-pdb", key=f"{key}_pdb_dl",
    )
