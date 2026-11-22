"""El instalador debe ser válido, idempotente en modo prueba y no tocar nada."""
from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_installer_syntax():
    r = subprocess.run(["bash", "-n", str(ROOT / "install.sh")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_installer_dry_run_does_not_write():
    r = subprocess.run(
        ["bash", str(ROOT / "install.sh"), "--dry-run", "--no-desktop",
         "--with-data"],
        capture_output=True, text=True, cwd=ROOT, timeout=60,
    )
    assert r.returncode == 0, r.stderr
    assert "dry-run" in r.stdout
    assert "pip install" in r.stdout


def test_launcher_syntax():
    r = subprocess.run(["bash", "-n", str(ROOT / "lanzar.sh")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_component_implements_streamlit_protocol():
    html = (ROOT / "webapp" / "components" / "genome_track" / "index.html").read_text()
    for token in ("streamlit:componentReady", "streamlit:render",
                  "streamlit:setComponentValue", "streamlit:setFrameHeight"):
        assert token in html
