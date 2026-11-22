"""Avisos de excepciones documentadas para gaps de un FASTA nuevo."""
import csv
from pathlib import Path


def load_index(data_root) -> dict:
    p = Path(data_root) / "metabolism" / "essential" / "essential_exceptions.tsv"
    out = {}
    if not p.exists():
        return out
    for r in csv.DictReader(open(p), delimiter="\t"):
        out.setdefault(r["map_id"], []).append(r)
    return out


def annotate(gaps: list[dict], data_root) -> list[dict]:
    idx = load_index(data_root)
    for g in gaps:
        paths = g.get("pathways", "") or ""
        hits = []
        for mid, rows in idx.items():
            if mid in paths:
                hits.extend(rows)
        if hits:
            g["known_exception"] = "; ".join(sorted({h["exception"] for h in hits}))
            g["exception_note"] = " | ".join(
                f"{h['exception']} (ej. {h['example_taxon']}): {h['implication']}" for h in hits[:3]
            )
        else:
            g["known_exception"] = ""
            g["exception_note"] = ""
    return gaps
