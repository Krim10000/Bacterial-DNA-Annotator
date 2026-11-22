"""Foldseek local por lotes contra una base preconstruida."""
from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path


def batch_search(pdbs: list[str | Path], target_db: str | Path,
                 workdir: str | Path, threads: int = 8,
                 max_seqs: int = 5) -> dict:
    """Devuelve {stem: [hits]} usando Foldseek local; no usa API."""
    t0 = time.perf_counter()
    pdbs = [Path(p) for p in pdbs if Path(p).exists()]
    if not pdbs:
        return {"hits": {}, "error": "sin PDB", "seconds": 0.0}
    work = Path(workdir)
    inp = work / "input"
    qdb = work / "query_db"
    res = work / "result_db"
    tmp = work / "tmp"
    m8 = work / "hits.m8"
    for p in (inp, tmp):
        p.mkdir(parents=True, exist_ok=True)
    for old in (qdb, res, m8):
        for f in work.glob(old.name + "*"):
            f.unlink(missing_ok=True)
    for i, p in enumerate(pdbs):
        link = inp / f"{p.stem}.pdb"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(p.resolve())
    exe = shutil.which("foldseek") or "foldseek"
    try:
        subprocess.run([exe, "createdb", str(inp), str(qdb), "--threads", str(threads)],
                       check=True, capture_output=True, text=True)
        subprocess.run([exe, "search", str(qdb), str(target_db), str(res), str(tmp),
                        "--threads", str(threads), "--max-seqs", str(max_seqs)],
                       check=True, capture_output=True, text=True)
        subprocess.run([exe, "convertalis", str(qdb), str(target_db), str(res), str(m8),
                        "--format-output", "query,target,fident,alnlen,evalue,bits"],
                       check=True, capture_output=True, text=True)
        hits: dict[str, list[dict]] = {}
        for line in m8.read_text().splitlines():
            c = line.split("\t")
            if len(c) < 6:
                continue
            q, target, seqid, aln, ev, bits = c[:6]
            acc = target.split("-")[1] if target.startswith("AF-") and "-" in target else target
            hits.setdefault(q, []).append({
                "pdb_id": acc.lower(), "target": target, "seq_id": seqid,
                "bits": bits, "evalue": ev, "aln_len": aln,
            })
        return {"hits": {k: v[:max_seqs] for k, v in hits.items()},
                "error": "", "seconds": round(time.perf_counter() - t0, 1)}
    except Exception as exc:
        return {"hits": {}, "error": f"{type(exc).__name__}: {exc}",
                "seconds": round(time.perf_counter() - t0, 1)}
