#!/usr/bin/env python3
"""04_build_domain_panel.py — descarga los HMM Pfam del lote y arma el panel.

Toma la unión de accesiones Pfam esperadas de `work/gap_domains.tsv` (65 en
este lote), descarga el HMM de cada familia desde la API de InterPro/Pfam
(`?annotation=hmm`, gzip), y concatena un único `work/panel_evidence.hmm`
presionado con `hmmpress`.

La descarga es reanudable: un `.hmm` ya presente y no vacío no se vuelve a
pedir. Si una familia no se puede descargar, se registra en
`work/panel_meta.json` y las brechas que dependan de ella quedan
`not_checked` en el eje de dominio (nunca `no_domain`).
"""
from __future__ import annotations

import gzip
import io
import json
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import WORK, read_tsv, write_json, write_tsv

HMM_DIR = WORK.parent / "work" / "pfam_hmm"
PANEL = WORK / "panel_evidence.hmm"
URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=hmm"
TRIES = 3


def fetch_one(acc: str) -> bool:
    path = HMM_DIR / f"{acc}.hmm"
    if path.exists() and path.stat().st_size > 0:
        return True
    url = URL.format(acc=acc)
    for attempt in range(TRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "dna_poc/evidence"})
            with urllib.request.urlopen(req, timeout=90) as fh:
                raw = fh.read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", "replace")
            if "HMMER3" not in text:
                raise ValueError("respuesta sin HMMER3")
            path.write_text(text)
            return True
        except Exception as exc:
            if attempt == TRIES - 1:
                print(f"  fallo {acc}: {exc}", file=sys.stderr)
            time.sleep(1.0 * (attempt + 1))
    return False


def parse_meta(text: str) -> dict:
    name = re.search(r"^NAME\s+(\S+)", text, re.M)
    acc = re.search(r"^ACC\s+(\S+)", text, re.M)
    desc = re.search(r"^DESC\s+(.+)$", text, re.M)
    leng = re.search(r"^LENG\s+(\d+)", text, re.M)
    return {"acc": acc.group(1) if acc else "", "name": name.group(1) if name else "",
            "desc": desc.group(1).strip() if desc else "",
            "leng": int(leng.group(1)) if leng else 0}


def main() -> int:
    HMM_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(WORK / "gap_domains.tsv")
    wanted = sorted({d for r in rows for d in r["domains"].split(";") if d})
    print(f"familias Pfam necesarias: {len(wanted)}")
    ok, failed = [], []
    for i, acc in enumerate(wanted, 1):
        (ok if fetch_one(acc) else failed).append(acc)
        if i % 10 == 0 or i == len(wanted):
            print(f"  {i}/{len(wanted)} ok={len(ok)} fallos={len(failed)}",
                  file=sys.stderr)

    meta = []
    blocks = []
    for acc in ok:
        text = (HMM_DIR / f"{acc}.hmm").read_text()
        m = parse_meta(text)
        m["requested_acc"] = acc
        meta.append(m)
        blocks.append(text if text.endswith("\n") else text + "\n")
    PANEL.write_text("".join(blocks))
    subprocess.run(["hmmpress", "-f", str(PANEL)], check=True,
                   stdout=subprocess.DEVNULL)
    write_tsv(WORK / "panel_meta.tsv", meta,
              ["requested_acc", "acc", "name", "desc", "leng"])
    stats = {
        "families_requested": len(wanted),
        "families_ok": len(ok),
        "families_failed": failed,
        "panel_models": len(blocks),
        "panel": str(PANEL),
    }
    write_json(WORK / "panel_meta.json", stats)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
