#!/usr/bin/env python3
"""Efectores por minado de texto de UniProt, restringido a nuestros 400 organismos.

Para cada organismo (taxid) se consultan las entradas Swiss-Prot revisadas que son TFs
y se extrae el efector mencionado en el comentario de función (p. ej. "repressed by
L-arginine", "in the presence of xylose").

Provenance: `uniprot_text` (precisión baja, se reporta con el fragmento de texto).

Salida: data/allostery/work/uniprot_effectors.tsv
Caché por taxid en data/allostery/work/uniprot_cache/ (permite reanudar).

Uso: python scripts/allostery/20_uniprot_effector_mining.py [--limit N]
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "allostery" / "work"
CACHE = WORK / "uniprot_cache"
OUT = WORK / "uniprot_effectors.tsv"

# Vocabulario de efectores/metabolitos (ampliable). Se busca con límites de palabra.
VOCAB = """
xylose arabinose lactose allolactose galactose glucose gluconate glucuronate trehalose
rhamnose fucose mannose mannitol sorbitol glycerol glycerate pyruvate acetate lactate
citrate succinate fumarate malate oxaloacetate 2-oxoglutarate ketoglutarate isocitrate
aconitate glyoxylate formate propionate butyrate ethanolamine choline betaine
alanine arginine asparagine aspartate cysteine glutamine glutamate glycine histidine
isoleucine leucine lysine methionine phenylalanine proline serine threonine tryptophan
tyrosine valine ornithine citrulline putrescine spermidine cadaverine agmatine
adenine guanine hypoxanthine xanthine uracil cytosine cytidine uridine adenosine
inosine guanosine thymine nicotinamide nicotinate NAD FAD FMN
S-adenosylmethionine adenosylmethionine SAM S-adenosylhomocysteine
coenzyme A biotin thiamine riboflavin cobalamin heme
iron zinc copper nickel cobalt manganese cadmium mercury lead silver arsenite arsenate
formaldehyde methylglyoxal glyoxal salicylate benzoate phenylacetate anthranilate
tetracycline IPTG isopropyl-beta-D-thiogalactopyranoside
N-acetylglucosamine glucosamine muramate ascorbate
trehalose-6-phosphate glucose-6-phosphate fructose-6-phosphate
ribose deoxyribose sucrose maltose cellobiose melibiose
""".split()

PATTERNS = [
    r"(?:repressed|induced|activated|inhibited)\s+by\s+(?:the\s+)?([\w\-\u2019']+)",
    r"in\s+(?:the\s+)?presence\s+of\s+([\w\-\u2019']+)",
    r"in\s+response\s+to\s+([\w\-\u2019']+)",
    r"binds\s+([\w\-\u2019']+)",
    r"([\w\-\u2019']+)[\-\s]responsive",
]
PAT = re.compile("|".join(PATTERNS), re.I)
VOCAB_SET = {v.lower() for v in VOCAB}


def query_uniprot(taxid: str, size: int = 400) -> str:
    # Ojo: Swiss-Prot anota a nivel de ESPECIE, no de cepa. Con el taxid de cepa
    # (lo que trae assembly_summary) la consulta devuelve vacío.
    q = (f"organism_id:{taxid} AND reviewed:true AND "
         f"(repressor OR \"transcription factor\" OR \"transcriptional regulator\")")
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({"query": q, "fields": "accession,protein_name,"
                                     "cc_function,gene_names,organism_name",
                                     "format": "tsv", "size": size}))
    req = urllib.request.Request(url, headers={"User-Agent": "alloreg"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read().decode()


CUE = re.compile(r"(presence|response|repress|induc|activat|inhibit|bind|bound|"
                 r"effector|ligand|senses|sensing|regulat|control|mediat)", re.I)


def extract(function_text: str, protein_name: str) -> list[dict]:
    """Extrae el efector solo si la frase tiene contexto de regulación.

    Confianza alta: el patrón explícito captura el término.
    Confianza media: el término aparece en una frase con verbo de regulación.
    """
    hits = []
    for sentence in re.split(r"(?<=[.;])\s+", function_text or ""):
        if not CUE.search(sentence):
            continue
        low = sentence.lower()
        strong = []
        for m in PAT.finditer(sentence):
            term = (m.group(1) or "").strip().strip(".,;:").lower()
            if len(term) >= 3 and term in VOCAB_SET:
                strong.append(term)
        weak = [v for v in VOCAB_SET
                if len(v) >= 3 and re.search(rf"\b{re.escape(v)}\b", low)]
        if strong:
            for v in dict.fromkeys(strong):
                hits.append({"efector": v, "snippet": sentence[:180],
                             "confianza": "alta"})
        elif weak:
            for v in dict.fromkeys(weak[:2]):
                hits.append({"efector": v, "snippet": sentence[:180],
                             "confianza": "media"})
    return hits[:2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    orgs = json.loads((WORK / "orgs400.json").read_text())
    if args.limit:
        orgs = orgs[:args.limit]

    rows = []
    n_err = 0
    for i, o in enumerate(orgs, 1):
        tax = o.get("species_taxid") or o.get("taxid")
        if not tax:
            continue
        cf = CACHE / f"{tax}.tsv"
        if cf.exists():
            text = cf.read_text()
        else:
            try:
                text = query_uniprot(tax)
                cf.write_text(text)
            except Exception as e:
                n_err += 1
                if n_err <= 3:
                    print(f"  error {tax}: {str(e)[:80]}")
                continue
            time.sleep(0.2)
        lines = text.splitlines()
        if len(lines) < 2:
            continue
        header = lines[0].split("\t")
        idx = {h: k for k, h in enumerate(header)}
        for line in lines[1:]:
            f = line.split("\t")
            if len(f) < len(header):
                continue
            acc = f[idx.get("Entry", 0)]
            name = f[idx.get("Protein names", 1)] if "Protein names" in idx else ""
            genes = f[idx.get("Gene Names", 2)] if "Gene Names" in idx else ""
            fcol = idx.get("Function [CC]", idx.get("Function", 3))
            fn = f[fcol] if fcol < len(f) else ""
            for h in extract(fn, name):
                rows.append({"taxid": tax, "organismo": o["organismo"],
                             "tf": (genes.split()[0] if genes else name[:40]),
                             "efector": h["efector"], "confianza": h["confianza"],
                             "accesion": acc,
                             "snippet": h["snippet"], "fuente": "uniprot_text"})
        if i % 50 == 0:
            print(f"  {i}/{len(orgs)} organismos · filas {len(rows)}", flush=True)

    with open(OUT, "w", newline="") as fh:
        cols = ["taxid", "organismo", "tf", "efector", "confianza", "accesion",
                "snippet", "fuente"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    tfs = {r["tf"].lower() for r in rows}
    summary = {"organismos": len(orgs), "errores": n_err, "filas": len(rows),
               "tfs_distintos": len(tfs),
               "organismos_con_efector": len({r["organismo"] for r in rows}),
               "efectores_distintos": len({r["efector"] for r in rows})}
    (WORK / "uniprot_effectors.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
