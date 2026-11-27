#!/usr/bin/env python3
"""Scraper de ENTRAF (668 fichas) con caché y mapeo a nuestras especies."""
from __future__ import annotations
import argparse, csv, json, re, time, urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[2]
BASE="http://web.pcyt.unam.mx/EntrafDB"
CITATION=("ENTRAF: An Encyclopedia of well-annotated DNA-binding TRanscription factors "
          "in Bacteria and Archaea. Computational Biology Lab, IIMAS-UNAM Sede Mérida.")
ACK=("Computational Biology Lab | IIMAS-UNAM Sede Mérida; Km. 5.5 Carr. Sierra Papacal - "
     "Chuburna Pto. Tablaje 31257, Sierra Papacal, 97302 Mérida, Yuc; "
     "+52 (999) 406 00 03 ext 7610; ernesto.perez@iimas.unam.mx; twitter.com/eperez09")


def fetch(url: str) -> str:
    req=urllib.request.Request(url,headers={"User-Agent":"dna-token-bac/1.0 (academic)"})
    with urllib.request.urlopen(req,timeout=60) as r: return r.read().decode("utf-8","replace")


def norm(s):
    s=(s or "").lower(); s=re.sub(r"\b(candidatus|subsp\.?|str\.?|strain)\b"," ",s)
    s=re.sub(r"\([^)]*\)"," ",s); return " ".join(re.sub(r"[^a-z0-9 ]+"," ",s).split()[:2])


def parse_entry(html: str) -> dict:
    soup=BeautifulSoup(html,"html.parser")
    out={}
    for table in soup.find_all("table"):
        for tr in table.find_all("tr"):
            tds=tr.find_all("td")
            if len(tds)<2: continue
            key=" ".join(tds[0].get_text(" ",strip=True).split())
            val=" ".join(tds[1].get_text(" ",strip=True).split())
            if key and val: out[key]=val
    # evidence table rows without key/value
    for table in soup.find_all("table"):
        for tr in table.find_all("tr"):
            tds=tr.find_all("td")
            if len(tds)==1:
                val=" ".join(tds[0].get_text(" ",strip=True).split())
                if val and val.lower() not in {"pdb","dnase footprinting","emsa"}:
                    out.setdefault("Experimental evidences", "")
                    out["Experimental evidences"]=(out["Experimental evidences"]+"; "+val).strip("; ")
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,default=ROOT/"data/regulatory/external/entraf")
    ap.add_argument("--limit",type=int,default=0)
    ap.add_argument("--delay",type=float,default=0.05)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    gp=a.out/"general_table.html"
    if not gp.exists(): gp.write_text(fetch(BASE+"/general_table.html"),encoding="utf-8")
    links=sorted(set(re.findall(r"individual_tables/([^\"']+\.html)",gp.read_text(encoding="utf-8",errors="replace"))))
    if a.limit: links=links[:a.limit]
    print("links",len(links),flush=True)
    # species map
    by=defaultdict(list)
    for d in (ROOT/"data/taxonomy/full").glob("*.json"):
        try: sp=json.loads(d.read_text()).get("species","")
        except Exception: continue
        if sp: by[norm(sp)].append(d.stem)
    manual={"bacillus subtilis":["subtilis"],"escherichia coli":["ecoli_k12"],
            "pseudomonas aeruginosa":["paeruginosa"]}
    raw=a.out/"raw"; raw.mkdir(exist_ok=True); rows=[]
    for i,name in enumerate(links,1):
        p=raw/name
        if p.exists(): html=p.read_text(encoding="utf-8",errors="replace")
        else:
            html=fetch(f"{BASE}/individual_tables/{quote(name)}"); p.write_text(html,encoding="utf-8"); time.sleep(a.delay)
        d=parse_entry(html)
        org=d.get("Organism",""); matches=by.get(norm(org), manual.get(norm(org), []))
        rows.append({"source_db":"ENTRAF","source_version":"web-2026-11-27",
                     "entraf_id":d.get("ENTRAF ID","") or name.split("_")[0],
                     "gene_name":d.get("Gene name",""),"uniprot_id":d.get("Uniprot ID",""),
                     "organism":org,"sequence":d.get("Aminoacid sequence",""),"length":d.get("Length",""),
                     "function":d.get("Function Uniprot",""),"role":d.get("Role",""),
                     "family":d.get("Protein family",""),"pdb":d.get("Structure (PDB ID)",""),
                     "pubmed":d.get("Pubmed ID",""),"evidence":d.get("Experimental evidences",""),
                     "genome_match":";".join(matches),"match_type":"direct_species" if matches else "no_match",
                     "source_url":f"{BASE}/individual_tables/{name}",
                     "citation":CITATION,"acknowledgement":ACK,
                     "retrieved_at":datetime.now(timezone.utc).isoformat()})
        if i%50==0: print(i,"/",len(links),flush=True)
    out=a.out/"entraf_entries.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["source_db"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"entries":len(rows),"species_matched":len({r['genome_match'] for r in rows if r['genome_match']}),
             "rows_matched":sum(1 for r in rows if r['genome_match']),"with_pdb":sum(1 for r in rows if r['pdb']),
             "with_pubmed":sum(1 for r in rows if r['pubmed']),"out":str(out),"citation":CITATION,"acknowledgement":ACK}
    (a.out/"entraf_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8"); print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
