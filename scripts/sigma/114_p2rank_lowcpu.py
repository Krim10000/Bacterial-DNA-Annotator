#!/usr/bin/env python3
"""P2Rank en segundo plano, solo con CPU < 90 %, nice 19, reanudable."""
from __future__ import annotations
import argparse, csv, json, os, subprocess, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def cpu_percent(prev=None):
    vals={}
    with open('/proc/stat') as fh:
        for line in fh:
            if line.startswith('cpu '):
                f=list(map(int,line.split()[1:])); total=sum(f); idle=f[3]+f[4]
                break
    if prev is None: return 0.0,(total,idle)
    dt=total-prev[0]; di=idle-prev[1]
    return (100.0*(dt-di)/dt if dt else 0.0),(total,idle)


def has_prediction(path: Path):
    return any(path.glob('*_predictions.csv')) if path.exists() else False


def build_queue(data_root: Path, outroot: Path):
    jobs=[]; seen=set()
    sources=[
        (data_root/"metabolism/targets/afdb25k","AFDB"),
        (data_root/"sigma/structures/pdb","PDB_experimental"),
    ]
    for base,origin in sources:
        if not base.exists(): continue
        for p in sorted(base.glob('*.pdb')):
            acc=p.stem
            if acc in seen: continue
            seen.add(acc)
            old1=data_root/"sigma/structures/p2rank"/acc
            old2=data_root/"sigma/structures/p2rank_sigma"
            out=outroot/acc
            if has_prediction(out) or has_prediction(old1): continue
            # no duplicar sigma ya cubiertos por p2rank_sigma
            if any((old2).glob(f"*__{acc}*/*_predictions.csv")): continue
            jobs.append({"accession":acc,"pdb":str(p),"origin":origin,"out":str(out)})
    q=outroot/"queue.tsv"
    with q.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(jobs[0]) if jobs else ["accession"],delimiter='\t'); w.writeheader(); w.writerows(jobs)
    return jobs,q


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--outroot",type=Path,default=ROOT/"data/sigma/structures/p2rank_lowcpu")
    ap.add_argument("--build",action="store_true")
    ap.add_argument("--daemon",action="store_true")
    ap.add_argument("--status",action="store_true")
    ap.add_argument("--threads",type=int,default=1)
    a=ap.parse_args(); a.outroot.mkdir(parents=True,exist_ok=True)
    status_path=a.outroot/"status.json"; qpath=a.outroot/"queue.tsv"
    if a.status:
        print(status_path.read_text() if status_path.exists() else "{}"); return
    if a.build:
        jobs,q=build_queue(a.data_root,a.outroot)
        status={"pid":0,"state":"queued","total":len(jobs),"done":0,"last":"","cpu":None,"updated":time.time()}
        status_path.write_text(json.dumps(status,indent=1)); print(json.dumps({"queue":len(jobs),"path":str(q)},indent=1)); return
    if not a.daemon: return
    jobs=list(csv.DictReader(open(qpath),delimiter='\t')) if qpath.exists() else []
    if not jobs:
        jobs,q=build_queue(a.data_root,a.outroot); jobs=list(jobs)
    done=0; prev=None
    status={"pid":os.getpid(),"state":"running","total":len(jobs),"done":0,"last":"","cpu":None,"updated":time.time()}
    status_path.write_text(json.dumps(status,indent=1))
    log=(a.outroot/"daemon.log").open("a")
    for job in jobs:
        out=Path(job["out"])
        if has_prediction(out):
            done+=1; continue
        while True:
            pct,prev=cpu_percent(prev)
            status.update({"cpu":round(pct,1),"updated":time.time()})
            status_path.write_text(json.dumps(status,indent=1))
            if pct < 90.0: break
            time.sleep(5)
        out.mkdir(parents=True,exist_ok=True)
        cmd=["nice","-n","19","p2rank","predict","-f",job["pdb"],"-o",str(out),"-threads",str(a.threads)]
        log.write(f"{time.time()} start {job['accession']} {job['origin']} cpu={pct:.1f}\n"); log.flush()
        try:
            subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        except Exception as exc:
            log.write(f"error {job['accession']} {exc}\n"); log.flush()
        done+=1
        status.update({"done":done,"last":job["accession"],"updated":time.time()})
        status_path.write_text(json.dumps(status,indent=1))
        log.write(f"{time.time()} done {job['accession']} {done}/{len(jobs)}\n"); log.flush()
    status.update({"state":"finished","done":done,"updated":time.time()})
    status_path.write_text(json.dumps(status,indent=1)); log.close()
    print(json.dumps(status,indent=1))
if __name__=="__main__": main()
