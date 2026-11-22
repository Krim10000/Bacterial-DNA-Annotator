#!/usr/bin/env python3
"""02_build_pfam_panel.py — extrae un panel Pfam curado del Pfam-A completo.

`hmmfetch -f` relee el fichero completo una vez por clave (medido: 4 claves →
3 min 24 s), inviable para 339 modelos. Aquí se hace **una sola pasada** en
binario sobre `.cache/pfam/Pfam-A.hmm` y se escriben solo los registros pedidos.

Entrada:  scripts/mobile/pfam_panel_spec.json  (curado a mano, por NAME)
          .cache/pfam/Pfam-A.hmm  (3,6 GB, HMMER3)
Salida:   data/mobile_typing/_work/pfam_panel/panel.hmm (+ índice hmmpress)
          data/mobile_typing/_work/pfam_panel/panel_manifest.tsv
          data/mobile_typing/_work/pfam_panel/panel_info.json
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SPEC = Path(__file__).resolve().parent / "pfam_panel_spec.json"
PFAM_HMM = ROOT / ".cache" / "pfam" / "Pfam-A.hmm"
PFAM_DAT = ROOT / ".cache" / "Pfam-A.hmm.dat.gz"
DEFAULT_OUT = ROOT / "data" / "mobile_typing" / "_work" / "pfam_panel"


def load_dat():
    """accession -> (name, description, ga_seq, ga_dom) desde el .dat de Pfam."""
    meta = {}
    name = acc = desc = ga = None
    with gzip.open(PFAM_DAT, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GF ID"):
                name = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC"):
                acc = line.split(None, 2)[2].strip().split(".")[0]
            elif line.startswith("#=GF DE"):
                desc = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF GA"):
                ga = line.split(None, 2)[2].strip()
            elif line.startswith("//"):
                if acc:
                    meta[acc] = (name or "", desc or "", ga or "")
                name = acc = desc = ga = None
    return meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    SPEC = args.spec.resolve()
    OUTDIR = args.outdir.resolve()
    spec = json.loads(SPEC.read_text())
    wanted = {}
    for category, accessions in spec["categories"].items():
        for acc in accessions:
            wanted.setdefault(acc, []).append(category)
    meta = load_dat()
    missing = sorted(a for a in wanted if a not in meta)
    if missing:
        print(f"ERROR: accesiones inexistentes en Pfam: {missing}", file=sys.stderr)
        return 1
    if not PFAM_HMM.exists():
        print(f"ERROR: falta {PFAM_HMM}", file=sys.stderr)
        return 1

    OUTDIR.mkdir(parents=True, exist_ok=True)
    panel = OUTDIR / "panel.hmm"
    found = set()
    acc_re = re.compile(rb"\nACC\s+(\S+)")
    buf = b""
    chunk_size = 8 * 1024 * 1024
    with PFAM_HMM.open("rb") as src, panel.open("wb") as dst:
        while True:
            chunk = src.read(chunk_size)
            if not chunk:
                break
            buf += chunk
            # los registros acaban en "\n//\n"; se conserva el resto
            parts = buf.split(b"\n//\n")
            buf = parts.pop()
            for record in parts:
                match = acc_re.search(record[:4096])
                if not match:
                    continue
                acc = match.group(1).decode().split(".")[0]
                if acc in wanted and acc not in found:
                    dst.write(record + b"\n//\n")
                    found.add(acc)
            # corta el buffer si ya se encontraron todos (no hace falta seguir)
            if len(found) == len(wanted):
                break
        # último registro sin "//" final
        if buf.strip() and len(found) < len(wanted):
            match = acc_re.search(buf[:4096])
            if match:
                acc = match.group(1).decode().split(".")[0]
                if acc in wanted and acc not in found:
                    dst.write(buf + b"\n//\n")
                    found.add(acc)

    missing = sorted(set(wanted) - found)
    if missing:
        print(f"ERROR: no se extrajeron {missing}", file=sys.stderr)
        return 1

    with (OUTDIR / "panel_manifest.tsv").open("w") as fh:
        fh.write("accession\tname\tdescription\tcategories\tga_seq\tga_dom\n")
        for acc in sorted(found):
            name, desc, ga = meta[acc]
            ga_seq, ga_dom = "", ""
            if ga:
                nums = [x.strip(";") for x in ga.split()]
                if len(nums) >= 2:
                    ga_seq, ga_dom = nums[0], nums[1]
            fh.write(f"{acc}\t{name}\t{desc}\t{','.join(wanted[acc])}\t{ga_seq}\t{ga_dom}\n")

    # índice HMMER (necesario para hmmscan)
    subprocess.run(["hmmpress", "-f", str(panel)], check=True, stdout=subprocess.DEVNULL)
    digest = hashlib.sha256(panel.read_bytes()).hexdigest()
    info = {
        "panel_spec": str(SPEC.relative_to(ROOT)),
        "pfam_source": str(PFAM_HMM.relative_to(ROOT)),
        "models": len(found),
        "categories": {c: len(v) for c, v in spec["categories"].items()},
        "panel_sha256": digest,
    }
    (OUTDIR / "panel_info.json").write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
    print(f"OK: {len(found)} modelos -> {panel.relative_to(ROOT)} (sha256 {digest[:16]}…)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
