#!/usr/bin/env python3
"""10 — Estadísticas agregadas de la capa completa (para el informe/GUI).

Lee los 400 JSON de `data/taxonomy/full/` y `ani_hits_full.tsv` y escribe
`data/taxonomy/full/stats_full.json` con la cobertura, el desglose de acuerdos y
desacuerdos NCBI/GTDB y la distribución de ANI. No modifica los JSON.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"


def strip_suffix(x: str | None) -> str | None:
    return re.sub(r"_[A-Z]$", "", x) if x else None


def epithet_key(x: str | None) -> str | None:
    if not x:
        return None
    e = x.lower()
    for s in ("ensis", "ense", "um", "us", "a", "e", "i", "o", "u"):
        if e.endswith(s) and len(e) > len(s) + 2:
            return e[: -len(s)]
    return e


def binomial(name: str | None):
    if not name:
        return None, None
    parts = name.split()
    if parts and parts[0].lower() == "candidatus" and len(parts) >= 3:
        parts = parts[1:]
    return (parts[0], parts[1]) if len(parts) >= 2 else (parts[0], None)


def main() -> int:
    docs = []
    for p in sorted(FULL.glob("GCF_*.json")):
        docs.append(json.loads(p.read_text(encoding="utf-8")))

    ani_values = [d["ani"] for d in docs if d["ani"] is not None]
    bins = Counter()
    for a in ani_values:
        if a >= 99:
            bins[">=99"] += 1
        elif a >= 95:
            bins["95-99"] += 1
        elif a >= 90:
            bins["90-95"] += 1
        elif a >= 80:
            bins["80-90"] += 1
        elif a >= 70:
            bins["70-80"] += 1
        else:
            bins["<70"] += 1

    species_agree = species_renamed = species_conflict = 0
    genus_same = genus_renamed = 0
    conflicts_list = []
    for d in docs:
        g_ncbi, s_ncbi = binomial(d.get("species_ncbi"))
        g_gtdb, s_gtdb = binomial(d.get("species_gtdb"))
        if d.get("species_gtdb") is None:
            continue
        if strip_suffix(s_ncbi) == strip_suffix(s_gtdb):
            species_agree += 1
        elif epithet_key(strip_suffix(s_ncbi)) == epithet_key(strip_suffix(s_gtdb)):
            species_agree += 1
        else:
            species_conflict += 1
            conflicts_list.append({"genome": d["genome"],
                                   "species_ncbi": d.get("species_ncbi"),
                                   "species_gtdb": d.get("species_gtdb")})
        if strip_suffix(g_ncbi) == strip_suffix(g_gtdb):
            genus_same += 1
        else:
            genus_renamed += 1
        if d["evidence"]["conflicts"]:
            pass

    # pares con ANI >= 95 entre genomas distintos
    over95 = []
    with open(FULL / "ani_hits_full.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if float(r["ani"]) >= 95.0:
                over95.append(r)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_genomes": len(docs),
        "coverage": {
            "ncbi_metadata": sum(1 for d in docs if d.get("ncbi_taxid")),
            "gtdb_lineage": sum(1 for d in docs if d["gtdb_lineage"] != "unknown"),
            "gtdb_missing": sum(1 for d in docs if d["gtdb_lineage"] == "unknown"),
            "species_named": sum(1 for d in docs if d["species"] != "unknown"),
            "ani": len(ani_values),
            "ani_null": len(docs) - len(ani_values),
        },
        "ani_bins": dict(bins),
        "confidence": dict(Counter(d["confidence"] for d in docs)),
        "status": dict(Counter(d["status"] for d in docs)),
        "scope": dict(Counter(d["evidence"]["ani"]["scope_rank"] for d in docs)),
        "species_consensus": {
            "agree_with_gtdb": species_agree,
            "conflict_with_gtdb": species_conflict,
            "genus_same": genus_same,
            "genus_renamed": genus_renamed,
        },
        "ani_species_level_pairs": len(over95),
        "species_conflicts": conflicts_list,
    }
    L.write_json(FULL / "stats_full.json", payload)
    print(json.dumps(payload, indent=2, ensure_ascii=False)[:2500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
