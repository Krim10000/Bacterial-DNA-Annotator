#!/usr/bin/env python3
"""fetch_dbtbs.py — DBTBS (*B. subtilis*): descarga manual, no automatizable.

Estado: BLOQUEADO POR DISEÑO (verificado el 2026-11-21)
-------------------------------------------------------
DBTBS (https://dbtbs.hgc.jp/) no publica sus datos en una ruta estática. La
secuencia real es:

    GET https://dbtbs.hgc.jp/           -> 200, 326 bytes (solo un <frame>)
    GET https://dbtbs.hgc.jp/ver4/      -> 200, 291 bytes (solo un <frame>)
    GET https://dbtbs.hgc.jp/download.html
        -> 200, texto: formulario con Name*, E-mail*, Institution*, Purpose*
           y JS que hace:  frm.action = "/cgi-bin/download.cgi"; submit()
    GET https://dbtbs.hgc.jp/robots.txt -> 200:  "Disallow: /cgi-bin/"
    GET /data/ /files/ sitemap.xml      -> 404

Es decir: los datos existen **detrás de un formulario de registro** para el que
hay que dar nombre, correo, institución y propósito. Un script no puede
rellenarlo sin inventar una identidad, y `robots.txt` pide explicitamente no
automatizar `/cgi-bin/`. Por eso este script NO descarga.

Alternativa integrada
---------------------
Se integro **BSGatlas** (ver `fetch_bsgatlas.py`), que es una anotacion abierta
de *B. subtilis* 168 que **incluye DBTBS entre sus fuentes** y se descarga en
GFF3 sin registro. Para el objetivo del proyecto (segunda especie curada) cubre
lo mismo y ademas trae terminadores con energia y TSS con sigma.

Uso de este script
------------------
1. Si el dueño descarga DBTBS a mano (formulario), deja el fichero en
   `data/dbtbs/dbtbs_manual.tsv`.
2. Este script lo normaliza al esquema de `data/dbtbs/gt_operons.tsv`.

    python scripts/validation/fetch_dbtbs.py            # comprueba e informa
    python scripts/validation/fetch_dbtbs.py --convert  # normaliza el manual

Columnas reconocidas del fichero manual (se buscan por nombre, sin importar el
orden ni mayusculas): operon/operon_name, start/leftEndPos, end/rightEndPos,
strand, genes/gene_name, pmid/pmids, evidence/evidences.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "data" / "dbtbs"
MANUAL = DIR / "dbtbs_manual.tsv"


def norm_key(k: str) -> str:
    return k.strip().lower().replace(")", "").replace(" ", "_")


def convert() -> int:
    if not MANUAL.exists():
        print(f"falta {MANUAL.relative_to(ROOT)}")
        print("descarga manualmente desde https://dbtbs.hgc.jp/download.html "
              "(formulario) y guardalo con ese nombre")
        return 1
    rows = []
    with open(MANUAL) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        mapk = {norm_key(k): k for k in (rd.fieldnames or [])}

        def get(r, *names):
            for n in names:
                if n in mapk:
                    return r.get(mapk[n], "") or ""
            return ""

        for r in rd:
            try:
                a = int(get(r, "start", "leftendpos"))
                b = int(get(r, "end", "rightendpos"))
            except ValueError:
                continue
            genes = [g for g in get(r, "genes", "gene_name").replace(",", ";").split(";") if g.strip()]
            rows.append({
                "operon_id": get(r, "operon", "operon_name", "id"),
                "operon_name": get(r, "operon_name", "operon"),
                "strand": get(r, "strand"),
                "start": min(a, b), "end": max(a, b),
                "n_genes": len(genes), "genes": ",".join(genes),
                "gene_names": get(r, "gene_name", "genes"),
                "n_genes_tus": len(genes), "mapping_ok": "yes" if genes else "no",
                "confidence": get(r, "confidence"),
                "evidence": get(r, "evidence", "evidences", "pmids"),
            })
    if not rows:
        print("no se pudo interpretar ninguna fila; revisa las columnas")
        return 1
    out = DIR / "gt_operons_dbtds.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} operones -> {out.relative_to(ROOT)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--convert", action="store_true")
    a = ap.parse_args()
    print("DBTBS: acceso automatizado bloqueado (formulario + robots.txt).")
    print("  ruta alternativa ya integrada: BSGatlas (ver fetch_bsgatlas.py)")
    if a.convert:
        return convert()
    print(f"  si tienes el fichero manual: {MANUAL.relative_to(ROOT)} "
          "y ejecuta --convert")
    return 0


if __name__ == "__main__":
    sys.exit(main())
