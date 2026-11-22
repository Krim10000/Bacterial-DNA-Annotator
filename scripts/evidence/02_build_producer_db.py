#!/usr/bin/env python3
"""02_build_producer_db.py — extrae secuencias de los productores y arma la DB.

Para las 5.472 accesiones UniProt del lote:
  * lee su secuencia y su organismo (OX=) de `data/db/uniprot_sprot.fasta`;
  * clasifica eucariota/procariota con `data/taxdump` (ancestro 2759 =
    Eukaryota, resolviendo taxones obsoletos con `merged.dmp`);
  * escribe `work/producer_accs.faa` (todas), `work/acc_info.tsv` y
    `work/no_euk.faa`;
  * construye la DB DIAMOND `work/producer_db.dmnd` solo con productores no
    eucariotas: una enzima de pollo no es evidencia para una bacteria.

Salidas en `data/metabolism/evidence/work/`. No toca la tabla base.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (SPROT, TAXDUMP, WORK, read_tsv, sprot_acc, write_json,
                    write_tsv)

EUKARYOTA = 2759


def load_parents() -> dict[int, int]:
    parents: dict[int, int] = {}
    with open(TAXDUMP / "nodes.dmp") as fh:
        for line in fh:
            p = line.split("\t|")
            if len(p) >= 2:
                parents[int(p[0])] = int(p[1].strip())
    return parents


def load_merged() -> dict[int, int]:
    merged: dict[int, int] = {}
    path = TAXDUMP / "merged.dmp"
    if not path.exists():
        return merged
    with open(path) as fh:
        for line in fh:
            p = line.split("\t|")
            if len(p) >= 2:
                merged[int(p[0])] = int(p[1].strip())
    return merged


def make_is_euk(parents: dict[int, int]):
    memo: dict[int, bool] = {}

    def is_euk(taxid: int) -> bool:
        if taxid in memo:
            return memo[taxid]
        path, seen, cur = [], set(), taxid
        while cur not in memo and cur not in seen and cur != 1:
            seen.add(cur)
            path.append(cur)
            cur = parents.get(cur, 1)
        base = memo.get(cur, cur == EUKARYOTA)
        for node in reversed(path):
            base = base or node == EUKARYOTA
            memo[node] = base
        return memo.get(taxid, base)

    return is_euk


def read_sprot(accs: set[str]) -> dict[str, tuple[str, str]]:
    """acc -> (cabecera completa, secuencia) solo para las accesiones pedidas."""
    out: dict[str, tuple[str, str]] = {}
    acc = None
    header = ""
    chunks: list[str] = []
    with open(SPROT) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = (header, "".join(chunks))
                header = line[1:].strip()
                acc = sprot_acc(header.split()[0])
                acc = acc if acc in accs else None
                chunks = []
            elif acc is not None:
                chunks.append(line)
        if acc is not None:
            out[acc] = (header, "".join(chunks))
    return out


def main() -> int:
    rows = read_tsv(WORK / "selected_gaps.tsv")
    accs = sorted({a for r in rows for a in r["producer_accs"].split(";") if a})
    print(f"accesiones del lote: {len(accs)}")

    parents = load_parents()
    merged = load_merged()
    is_euk = make_is_euk(parents)
    seqs = read_sprot(set(accs))

    info: dict[str, dict] = {}
    faa_all = WORK / "producer_accs.faa"
    faa_prok = WORK / "no_euk.faa"
    with open(faa_all, "w") as all_out, open(faa_prok, "w") as prok_out:
        for acc in sorted(seqs):
            header, seq = seqs[acc]
            m = re.search(r"OX=(\d+)", header)
            taxid = merged.get(int(m.group(1)), int(m.group(1))) if m else 0
            euk = bool(taxid and is_euk(taxid))
            info[acc] = {"acc": acc, "organism": header, "taxid": taxid,
                         "euk": euk}
            block = f">{acc} {header}\n{seq}"
            all_out.write(block)
            if not euk:
                prok_out.write(block)

    missing = sorted(set(accs) - set(seqs))
    write_tsv(
        WORK / "acc_info.tsv",
        [{"acc": a, "taxid": info[a]["taxid"], "euk": int(info[a]["euk"]),
          "organism": info[a]["organism"]} for a in sorted(info)],
        ["acc", "taxid", "euk", "organism"],
    )

    db = WORK / "producer_db"
    subprocess.run(["diamond", "makedb", "--in", str(faa_prok), "-d", str(db),
                    "--quiet"], check=True)
    stats = {
        "accs_requested": len(accs),
        "accs_found_in_sprot": len(seqs),
        "accs_missing": len(missing),
        "missing_sample": missing[:20],
        "euk_accs_excluded": sum(1 for a in info if info[a]["euk"]),
        "prok_accs_in_db": sum(1 for a in info if not info[a]["euk"]),
        "db": str(db) + ".dmnd",
    }
    write_json(WORK / "producer_db.json", stats)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
