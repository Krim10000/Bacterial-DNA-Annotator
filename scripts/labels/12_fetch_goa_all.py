#!/usr/bin/env python3
"""Descarga GAF/GOA disponibles por especie y los mapea por taxid."""
from __future__ import annotations
import gzip, json, re, subprocess, time
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def fetch(url: str, dest: Path):
    if dest.exists() and dest.stat().st_size > 0: return "cached"
    dest.parent.mkdir(parents=True,exist_ok=True)
    p=subprocess.run(["curl","-sL","--fail","--retry","3","--max-time","300",
                      "-A","dna-token-bac/1.0",url,"-o",str(dest)],
                     capture_output=True,text=True)
    if p.returncode: raise RuntimeError(p.stderr.strip() or f"curl {p.returncode}")
    return "ok"


def main() -> int:
    ext=ROOT/"data/labels/external"
    idx_path=ROOT/".scratch/gaf_index2.html"
    if not idx_path.exists():
        fetch("https://current.geneontology.org/annotations/gaf/",idx_path)
    idx=idx_path.read_text(encoding="utf-8",errors="replace")
    urls=sorted(set(re.findall(r'https?://[^"\']+?\.gaf\.gz',idx)))
    urls=[u for u in urls if u.endswith("-uniprot.gaf.gz") or u.endswith("-mod.gaf.gz")]
    print("files",len(urls),flush=True)
    errors=[]
    for i,u in enumerate(urls,1):
        name=u.rsplit("/",1)[-1]; dest=ext/"goa"/"by_code"/name
        try:
            fetch(u,dest)
        except Exception as exc:
            errors.append({"file":name,"error":str(exc)})
        if i%20==0: print(i,"/",len(urls),flush=True)
        time.sleep(0.05)
    # taxid from first data lines
    tax={}
    for p in sorted((ext/"goa"/"by_code").glob("*.gaf.gz")):
        try:
            with gzip.open(p,"rt",errors="replace") as fh:
                for line in fh:
                    if line.startswith("!"): continue
                    f=line.rstrip("\n").split("\t")
                    m=next((x for x in f if x.startswith("taxon:")),None)
                    if m: tax[p.name]=m.split(":",1)[1].split("|")[0]; break
        except Exception: pass
    # map local genomes by ncbi_taxid
    rows=[]
    for d in sorted((ROOT/"data/regulation").iterdir()):
        if not d.is_dir() or not (d/"genes.tsv").exists(): continue
        g=d.name; t=ROOT/"data/taxonomy/full"/f"{g}.json"; tid=""
        if t.exists():
            try: tid=str(json.loads(t.read_text()).get("ncbi_taxid") or "")
            except Exception: pass
        code=next((k for k,v in tax.items() if tid and v==tid), "")
        rows.append({"genome":g,"ncbi_taxid":tid,"goa_file":code,
                     "status":"ok" if code else "no_data"})
    (ext/"goa"/"goa_map.tsv").write_text("\t".join(["genome","ncbi_taxid","goa_file","status"])+"\n"+
        "\n".join("\t".join([r["genome"],r["ncbi_taxid"],r["goa_file"],r["status"]]) for r in rows)+"\n",encoding="utf-8")
    summary={"gaf_files":len(urls),"downloaded_taxids":len(tax),"genomes":len(rows),
             "mapped":sum(1 for r in rows if r["status"]=="ok"),"errors":errors}
    (ext/"goa"/"goa_fetch_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k!="errors"},indent=1)); print("errors",len(errors))
    return 0

if __name__=="__main__": raise SystemExit(main())
