#!/usr/bin/env python3
"""Transfiere sitios activos del homólogo al ORF por alineación de secuencia.

No sustituye alineación estructural ni experimento.
"""
import csv, json
from collections import Counter, defaultdict
from pathlib import Path
D=Path(__file__).resolve().parents[2]; E=D/"data"/"metabolism"/"essential"; T=D/"data"/"metabolism"/"targets"

def fasta(path):
    out={}; name=None; buf=[]
    for line in open(path,errors="ignore"):
        if line.startswith(">"):
            if name: out[name]="".join(buf)
            h=line[1:].split()[0]; name=h.split("|")[1] if "|" in h else h; buf=[]
        else: buf.append(line.strip())
    if name: out[name]="".join(buf)
    return out

def nw(a,b):
    n,m=len(a),len(b); NEG=-2; MATCH=2; MIS=-1
    dp=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=dp[i-1][0]+NEG
    for j in range(1,m+1): dp[0][j]=dp[0][j-1]+NEG
    for i in range(1,n+1):
        ai=a[i-1]
        for j in range(1,m+1):
            dp[i][j]=max(dp[i-1][j-1]+(MATCH if ai==b[j-1] else MIS), dp[i-1][j]+NEG, dp[i][j-1]+NEG)
    i,j=n,m; mapr={}; pairs=[]
    while i>0 and j>0:
        if dp[i][j]==dp[i-1][j-1]+(MATCH if a[i-1]==b[j-1] else MIS):
            mapr[i]=j; pairs.append((a[i-1],b[j-1])); i-=1; j-=1
        elif dp[i][j]==dp[i-1][j]+NEG: i-=1
        else: j-=1
    return mapr,sum(1 for x,y in pairs if x==y),len(pairs)

refs=fasta(D/"data/db/uniprot_sprot.fasta")
def fasta_full(path):
    out={}; name=None; buf=[]
    for line in open(path,errors="ignore"):
        if line.startswith(">"):
            if name: out[name]="".join(buf)
            name=line[1:].strip(); buf=[]
        else: buf.append(line.strip())
    if name: out[name]="".join(buf)
    return out
cands=fasta_full(T/"structural_candidates.faa")
struct={r["target"]:r for r in csv.DictReader(open(T/"STRUCTURAL_CANDIDATES.tsv"),delimiter="\t")}
asites=list(csv.DictReader(open(E/"ACTIVE_SITES_ESSENTIAL.tsv"),delimiter="\t"))
rows=[]
for r in asites:
    if r["feature"]=="none": continue
    s=struct.get(r["target"])
    if not s: continue
    acc=r["acc"]; ref=refs.get(acc,""); qry=cands.get(f"{s['genome']}|{s['protein_id']}","")
    if not ref or not qry: continue
    m,matches,paired=nw(ref,qry); ident=round(matches/max(1,paired),3)
    pos=int(r["start"]) if r["start"] else 0
    q=m.get(pos); qaa=qry[q-1] if q else ""; raa=ref[pos-1] if pos and pos<=len(ref) else ""
    rows.append({**r,"candidate":s["protein_id"],"best_acc":s["best_acc"],"ref_pos":pos,"ref_aa":raa,
                 "qry_pos":q or "","qry_aa":qaa,"conserved":("yes" if q and raa==qaa else "no" if q else "unmapped"),
                 "seq_identity":ident,"ref_len":len(ref),"qry_len":len(qry)})
with open(E/"ACTIVE_SITE_TRANSFER.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
by=defaultdict(list)
for r in rows: by[(r["target"],r["candidate"])].append(r)
with open(E/"ACTIVE_SITE_TRANSFER.md","w") as fh:
    fh.write("# Transferencia de sitios activos por secuencia\n\n")
    fh.write("Nivel secuencia; no estructura ni experimento. `conserved` = mismo residuo en la alineación global.\n\n")
    fh.write("| target | candidato | features | conservadas | identidad |\n|---|---|---:|---:|---:|\n")
    for (t,c),rs in by.items():
        fh.write(f"| {t} | {c} | {len(rs)} | {sum(1 for x in rs if x['conserved']=='yes')} | {rs[0]['seq_identity']} |\n")
    fh.write("\n## Límites\n\n- Alineación global de secuencia, no estructural.\n- `Binding site` no es necesariamente catálisis.\n- No cierra la brecha.\n")
s={"targets":len(by),"features":len(rows),"conserved":sum(1 for r in rows if r["conserved"]=="yes"),"unmapped":sum(1 for r in rows if r["conserved"]=="unmapped")}
json.dump(s,open(E/"ACTIVE_SITE_TRANSFER.json","w"),indent=2)
print(s)
