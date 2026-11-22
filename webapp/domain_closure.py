"""Cierre de brechas por dominio (A1).

Para una brecha: toma los productores, obtiene sus Pfam vía UniProt y escanea
las proteínas del genoma con esos dominios. Requiere el dominio específico
(no uno genérico) para proponer cierre.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
import urllib.request
from pathlib import Path

# Dominios demasiado genéricos para concluir un cierre por sí solos.
GENERIC = {"FAD_binding_4", "NAD_binding_1", "NAD_binding_2", "Rossmann-like",
           "Aldedh", "Aminotran_1_2", "Aminotran_3", "Aminotran_4", "Aminotran_5"}


def pfam_domains(acc: str, timeout: int = 20) -> set[str]:
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        d = json.load(urllib.request.urlopen(url, timeout=timeout))
    except Exception:
        return set()
    out = set()
    for x in d.get("uniProtKBCrossReferences", []):
        if x.get("database") == "Pfam":
            for p in x.get("properties", []):
                if p.get("key") == "EntryName":
                    out.add(p["value"])
    return out


def scan(faa: Path, domains: set[str], pfam_hmm: Path) -> list[dict]:
    if not domains:
        return []
    with tempfile.TemporaryDirectory() as td:
        nr = Path(td) / "names.txt"
        nr.write_text("\n".join(sorted(domains)))
        panel = Path(td) / "panel.hmm"
        with open(panel, "w") as fh:
            subprocess.run(["hmmfetch", "-f", str(pfam_hmm), str(nr)], stdout=fh, check=True)
        subprocess.run(["hmmpress", "-f", str(panel)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = Path(td) / "hits.domtbl"
        subprocess.run(["hmmscan", "--cut_ga", "--domtblout", str(out), "--noali",
                        str(panel), str(faa)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        hits = []
        for line in open(out):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 14:
                hits.append({"domain": f[0], "query": f[3],
                             "evalue": float(f[6]), "score": float(f[7]),
                             "specific": f[0] not in GENERIC})
        return hits


def close(accs: set[str], faa: Path, pfam_hmm: Path) -> dict:
    domains = set()
    for acc in list(accs)[:3]:
        domains |= pfam_domains(acc)
    hits = scan(faa, domains, pfam_hmm)
    specific = [h for h in hits if h["specific"]]
    return {"domains": sorted(domains), "hits": hits,
            "specific_hits": specific,
            "status": "candidate_closure" if specific else "no_conclusive"}
