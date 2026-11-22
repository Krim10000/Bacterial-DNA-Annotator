#!/usr/bin/env python3
"""ESMFold en Neurosnap para un FASTA pequeño.

Requiere ESM_API_KEY en el entorno o en .env. No imprime la clave.
"""
import argparse, hashlib, json, os, time
from pathlib import Path
import requests

BASE="https://neurosnap.ai/api"
def key():
    k=os.environ.get("ESM_API_KEY")
    if k: return k
    p=Path(".env")
    if p.exists():
        for line in p.read_text().splitlines():
            if line.startswith("ESM_API_KEY="):
                return line.split("=",1)[1].strip().strip('"').strip("'")
    raise SystemExit("Falta ESM_API_KEY")

def fasta(path):
    out=[]; name=None; buf=[]
    for line in open(path):
        if line.startswith(">"):
            if name: out.append((name,"".join(buf)))
            name=line[1:].strip(); buf=[]
        else: buf.append(line.strip())
    if name: out.append((name,"".join(buf)))
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fasta",required=True); ap.add_argument("--out",required=True)
    ap.add_argument("--limit",type=int,default=1); ap.add_argument("--offset",type=int,default=0)
    ap.add_argument("--note",default="dna_poc test")
    a=ap.parse_args(); K=key(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    for name,seq in fasta(a.fasta)[a.offset:a.offset+a.limit]:
        parts = name.split("|")
        g = "".join(c if c.isalnum() or c in "_-" else "_" for c in parts[0])[:28]
        p = "".join(c if c.isalnum() or c in "_-" else "_" for c in (parts[-1] if len(parts) > 1 else "protein"))[-18:]
        safe = f"{g}_{p}_{hashlib.sha1(name.encode()).hexdigest()[:6]}"[:55] or "protein"
        payload=json.dumps({"aa":{safe:seq},"dna":{},"rna":{}})
        files={"Input Sequences":(None,payload),"Number Recycles":(None,"6")}
        r=requests.post(f"{BASE}/job/submit/ESMFold",params={"note":a.note},
                        headers={"X-API-KEY":K},files=files,timeout=60)
        if not r.ok:
            print("submit error",r.status_code,r.text[:500],flush=True); r.raise_for_status()
        jid=r.json()
        if isinstance(jid,dict): jid=jid.get("job_id") or jid.get("jobId")
        print("submitted",jid,flush=True)
        for _ in range(120):
            s=requests.get(f"{BASE}/job/status/{jid}",headers={"X-API-KEY":K},timeout=30).json()
            st=s if isinstance(s,str) else (s.get("Status") or s.get("status"))
            if st in ("completed","Complete","stopped","failed"): break
            time.sleep(5)
        print("status",st,flush=True)
        data=requests.get(f"{BASE}/job/data/{jid}",headers={"X-API-KEY":K},timeout=30).json()
        (out/f"{safe}_job.json").write_text(json.dumps({"job_id":jid,"status":st,"data":data},indent=2))
        for fn,_ in data.get("out",[]):
            u=f"{BASE}/job/file/{jid}/out/{requests.utils.quote(fn)}"
            rr=requests.get(u,headers={"X-API-KEY":K},timeout=120,allow_redirects=True)
            rr.raise_for_status(); (out/fn).write_bytes(rr.content); print("saved",fn,flush=True)

if __name__=="__main__": main()
