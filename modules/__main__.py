"""CLI de los módulos.

    python -m modules --list
    python -m modules --describe annotation
    python -m modules --selftest annotation
    python -m modules --selftest all
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .registry import discover, get

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--describe", metavar="ID")
    ap.add_argument("--selftest", metavar="ID|all")
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    args = ap.parse_args()

    mods = discover()
    if args.list or not any((args.describe, args.selftest)):
        for m in mods:
            st = m.status(args.data_root)
            print(f"{st.badge()} {m.id:15s} {m.name:38s} "
                  f"opciones={len(m.options):2d}  {'GPU' if m.gpu else 'CPU'}")
        return 0

    if args.describe:
        m = get(args.describe)
        if not m:
            print("No existe:", args.describe)
            return 2
        print(json.dumps({
            "id": m.id, "name": m.name, "summary": m.summary,
            "version": m.version, "gpu": m.gpu, "experimental": m.experimental,
            "requires": list(m.requires), "docs": m.docs,
            "options": [o.__dict__ for o in m.options],
            "status": m.status(args.data_root).__dict__,
            "estimate": m.estimate(m.defaults()),
        }, ensure_ascii=False, indent=2))
        return 0

    targets = mods if args.selftest == "all" else [get(args.selftest)]
    if any(t is None for t in targets):
        print("No existe:", args.selftest)
        return 2
    failed = 0
    for m in targets:
        try:
            msg = m.selftest()
            print(f"✅ {m.id}: {msg}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"❌ {m.id}: {exc}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
