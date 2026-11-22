#!/usr/bin/env python3
"""Orquestador del módulo 01_annotation."""
from __future__ import annotations
import argparse, logging, os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS_ROOT = HERE.parent
LOG = logging.getLogger(HERE.name)

# (script, soporta --force)
STEPS: list[tuple[str, bool]] = [
    ("00_annotate_single.py", False),
    ("01_annotate_batch.py", True),
]

def run_step(script: str, force: bool, passes_force: bool) -> int:
    path = HERE / script
    if not path.exists():
        LOG.error("No encontrado: %s", path); return 127
    cmd = [sys.executable, str(path)]
    if force and passes_force:
        cmd.append("--force")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SCRIPTS_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    LOG.info("▶ %s", " ".join(cmd))
    return subprocess.run(cmd, env=env).returncode

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%H:%M:%S")
    for script, passes_force in STEPS:
        rc = run_step(script, args.force, passes_force)
        if rc != 0:
            LOG.error("✗ Falló %s (rc=%d)", script, rc); return rc
    LOG.info("✅ Módulo %s completado.", HERE.name); return 0

if __name__ == "__main__":
    sys.exit(main())
