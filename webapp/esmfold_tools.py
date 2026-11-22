"""ESMFold como herramienta de cierre de brechas.

Dos rutas, siempre declaradas al usuario:
- local: requiere GPU/RAM/disco y dependencias;
- Neurosnap API: X-API-KEY, pago por crédito.

No ejecuta nada automáticamente: el llamador decide y muestra el aviso.
"""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import requests

BASE = "https://neurosnap.ai/api"
PROVIDER_LABEL = "API Neurosnap"
LOCAL_LABEL = "ESMFold local"


def _gpu() -> tuple[str, float | None]:
    """Devuelve (nombre, VRAM GB) leyendo nvidia-smi si existe."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return "", None
    try:
        out = subprocess.check_output(
            [exe, "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
            text=True, timeout=3, stderr=subprocess.DEVNULL,
        ).strip().splitlines()
        if not out:
            return "", None
        name, mb = out[0].rsplit(",", 1)
        return name.strip(), round(float(mb) / 1024, 1)
    except Exception:
        return "", None


def _ram_gb() -> float | None:
    try:
        import psutil  # type: ignore
        return round(psutil.virtual_memory().total / 1024**3, 1)
    except Exception:
        try:
            pages = os.sysconf("SC_PHYS_PAGES")
            size = os.sysconf("SC_PAGE_SIZE")
            return round(pages * size / 1024**3, 1)
        except Exception:
            return None


def detect_capability() -> dict[str, Any]:
    """Evalúa la máquina donde corre la GUI/pipeline (no el navegador)."""
    gpu, vram = _gpu()
    ram = _ram_gb()
    try:
        disk = round(shutil.disk_usage(Path.cwd()).free / 1024**3, 1)
    except Exception:
        disk = None
    try:
        import torch  # type: ignore
        torch_version = getattr(torch, "__version__", "?")
        cuda = bool(torch.cuda.is_available())
        if cuda and not gpu:
            gpu = torch.cuda.get_device_name(0)
            vram = round(torch.cuda.get_device_properties(0).total_memory / 1024**3, 1)
    except Exception:
        torch_version, cuda = "", False
    deps_ok = all(_has(p) for p in ("esm", "openfold", "torch"))
    misses = [p for p in ("esm", "openfold", "torch") if not _has(p)]
    reasons: list[str] = []
    if not gpu:
        reasons.append("sin GPU NVIDIA visible; ESMFold 3B en CPU es horas y no se recomienda")
    elif vram is not None and vram < 10:
        reasons.append(f"VRAM {vram:.1f} GB: insuficiente para ESMFold 3B")
    if ram is not None and ram < 16:
        reasons.append(f"RAM {ram:.1f} GB: insuficiente")
    if disk is not None and disk < 15:
        reasons.append(f"disco libre {disk:.1f} GB: insuficiente para pesos/caché")
    if not cuda:
        reasons.append("PyTorch sin CUDA o no instalado")
    if misses:
        reasons.append("faltan dependencias: " + ", ".join(misses))

    if gpu and (vram or 0) >= 16 and (ram or 0) >= 32 and (disk or 0) >= 20 and cuda:
        level = "recommended"
    elif gpu and (vram or 0) >= 10 and (ram or 0) >= 16 and (disk or 0) >= 15:
        level = "possible_slow"
    else:
        level = "not_recommended"
    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "ram_gb": ram,
        "disk_free_gb": disk,
        "gpu": gpu,
        "vram_gb": vram,
        "cuda": cuda,
        "torch": torch_version,
        "level": level,
        "deps_ok": deps_ok,
        "can_run": level in ("recommended", "possible_slow") and deps_ok,
        "reasons": reasons,
    }


def _has(module: str) -> bool:
    import importlib.util
    return importlib.util.find_spec(module) is not None


def capability_rows(cap: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {"item": "Sistema", "valor": cap["platform"]},
        {"item": "Python", "valor": cap["python"]},
        {"item": "CPU", "valor": cap["cpu_count"]},
        {"item": "RAM (GB)", "valor": cap["ram_gb"]},
        {"item": "Disco libre (GB)", "valor": cap["disk_free_gb"]},
        {"item": "GPU", "valor": cap["gpu"] or "no detectada"},
        {"item": "VRAM (GB)", "valor": cap["vram_gb"] if cap["vram_gb"] is not None else "n/d"},
        {"item": "CUDA", "valor": cap["cuda"]},
        {"item": "PyTorch", "valor": cap["torch"] or "no instalado"},
        {"item": "esm/openfold", "valor": "sí" if cap.get("deps_ok") else "faltan"},
        {"item": "Nivel", "valor": cap["level"]},
    ]
    return [{"item": r["item"], "valor": "n/d" if r["valor"] is None else str(r["valor"])} for r in rows]


def capability_markdown(cap: dict[str, Any]) -> str:
    if cap["level"] == "recommended":
        head = "✅ **La máquina puede correr ESMFold local** (GPU y recursos suficientes)."
    elif cap["level"] == "possible_slow":
        head = "⚠️ **ESMFold local es posible pero lento/ajustado**; no lo uses en lote."
    else:
        head = "❌ **Esta máquina no es apta para ESMFold local.**"
    reasons = "\n".join(f"- {r}" for r in cap["reasons"]) or "- sin bloqueos detectados"
    return f"{head}\n\n{reasons}"


def local_command(fasta: str, outdir: str = "esmfold_local") -> str:
    return (
        "pip install 'fair-esm[esmfold]' 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'\n"
        f"python scripts/structure/04_local_esmfold.py --fasta {fasta} --out {outdir}"
    )


def sequence_for_gap(gap: dict[str, Any], proteins: list[dict[str, Any]]) -> tuple[str, str]:
    """Devuelve (nombre, secuencia) del ORF de la brecha, si existe."""
    q = str(gap.get("query") or gap.get("homolog_query") or "").strip()
    for p in proteins or []:
        if str(p.get("id", "")) == q:
            return q, str(p.get("aa", ""))
    # fallback: query puede venir con prefijo/sufijo
    for p in proteins or []:
        pid = str(p.get("id", ""))
        if q and (q in pid or pid in q):
            return pid, str(p.get("aa", ""))
    return q, ""


def _key() -> str:
    k = os.environ.get("ESM_API_KEY", "")
    if k:
        return k
    for p in (Path.cwd() / ".env", Path(__file__).resolve().parents[1] / ".env"):
        if p.exists():
            for line in p.read_text().splitlines():
                if line.startswith("ESM_API_KEY="):
                    return line.split("=", 1)[1].strip().strip('"').strip("'")
    return ""


def has_key() -> bool:
    return bool(_key())


_LOCAL_MODEL = None


def local_fold(sequence: str, name: str, outdir: Path | str) -> dict[str, Any]:
    """Pliega con fair-esm local. Solo llamar si capability.can_run es True."""
    global _LOCAL_MODEL
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    try:
        import esm
        import torch
        if _LOCAL_MODEL is None:
            _LOCAL_MODEL = esm.pretrained.esmfold_v1().eval()
            if torch.cuda.is_available():
                _LOCAL_MODEL = _LOCAL_MODEL.cuda()
        seq = sequence.replace("*", "").upper()
        pdb = _LOCAL_MODEL.infer_pdb(seq)
        safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)[:45] or "protein"
        path = out / f"{safe}.pdb"
        path.write_text(pdb)
        vals = []
        for line in pdb.splitlines():
            if line.startswith("ATOM"):
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
        plddt = round(sum(vals) / len(vals), 2) if vals else None
        return {"job_id": "local", "status": "completed", "pdb": str(path),
                "plddt": plddt, "ptm": None, "provider": "local"}
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}", "provider": "local"}


def submit_to_neuros_nap(sequence: str, name: str, note: str = "dna_poc gap closure",
                       outdir: Path | str = ".", timeout: int = 300) -> dict[str, Any]:
    """Envía un ESMFold a Neurosnap, espera y descarga PDB/métricas.

    Devuelve dict con job_id/status/pdb/plddt/ptm/error. Nunca imprime la key.
    """
    k = _key()
    if not k:
        return {"error": "Falta ESM_API_KEY"}
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)[:45] or "protein"
    payload = json.dumps({"aa": {safe: sequence}, "dna": {}, "rna": {}})
    files = {"Input Sequences": (None, payload), "Number Recycles": (None, "6")}
    try:
        r = requests.post(f"{BASE}/job/submit/ESMFold", params={"note": note},
                          headers={"X-API-KEY": k}, files=files, timeout=60)
        if not r.ok:
            return {"error": f"submit {r.status_code}: {r.text[:200]}"}
        jid = r.json()
        if isinstance(jid, dict):
            jid = jid.get("job_id") or jid.get("jobId")
        status = "pending"
        deadline = time.time() + timeout
        while time.time() < deadline:
            s = requests.get(f"{BASE}/job/status/{jid}", headers={"X-API-KEY": k}, timeout=30).json()
            status = s if isinstance(s, str) else (s.get("Status") or s.get("status"))
            if status in ("completed", "Complete", "stopped", "failed"):
                break
            time.sleep(5)
        data = requests.get(f"{BASE}/job/data/{jid}", headers={"X-API-KEY": k}, timeout=30).json()
        (out / f"{safe}_job.json").write_text(json.dumps({"job_id": jid, "status": status, "data": data}, indent=2))
        pdb_path = ""
        plddt = ptm = None
        for fn, _size in data.get("out", []):
            r2 = requests.get(f"{BASE}/job/file/{jid}/out/{requests.utils.quote(fn)}",
                              headers={"X-API-KEY": k}, timeout=120, allow_redirects=True)
            r2.raise_for_status()
            target = out / fn
            target.write_bytes(r2.content)
            if fn.endswith(".pdb"):
                pdb_path = str(target)
            if fn.endswith("metrics.json"):
                m = json.loads(r2.content)
                p = m.get("plddt") or []
                plddt = round(sum(p) / len(p), 2) if p else None
                ptm = m.get("ptm")
        return {"job_id": jid, "status": status, "pdb": pdb_path,
                "plddt": plddt, "ptm": ptm, "provider": "neuros nap"}
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def close_gaps(gaps: list[dict[str, Any]], proteins: list[dict[str, Any]],
               provider: str, outdir: Path | str, limit: int = 3,
               note: str = "dna_poc gap closure") -> list[dict[str, Any]]:
    """Cierre por estructura para las primeras brechas con ORF.

    provider: 'neuros nap' | 'local'. No ejecuta local si no es apto.
    """
    out = Path(outdir)
    rows: list[dict[str, Any]] = []
    if provider == "local":
        cap = detect_capability()
        if not cap["can_run"]:
            raise RuntimeError("ESMFold local no viable: " + "; ".join(cap["reasons"]))
        if not cap.get("deps_ok"):
            raise RuntimeError("ESMFold local no tiene esm/openfold/torch instalados")
    done = 0
    for g in gaps:
        if done >= limit:
            break
        name, seq = sequence_for_gap(g, proteins)
        if not seq or len(seq) < 10:
            continue
        if provider == "neuros nap":
            res = submit_to_neuros_nap(seq, name, note=note, outdir=out)
        else:
            res = local_fold(seq, name, out)
        plddt = res.get("plddt")
        evidence = "structure_predicted" if (plddt or 0) >= 80 and not res.get("error") else "candidate"
        g["closure_status"] = evidence if not res.get("error") else "structure_error"
        g["structure_provider"] = provider
        g["structure_pdb"] = res.get("pdb", "")
        g["structure_plddt"] = plddt
        g["structure_ptm"] = res.get("ptm")
        g["structure_limitation"] = "Predicción; no confirmado experimentalmente."
        rows.append({
            "genome": g.get("genome", ""), "gap": g.get("gap", ""), "query": name,
            "provider": provider, "job_id": res.get("job_id", ""), "status": res.get("status", ""),
            "plddt": plddt, "ptm": res.get("ptm"), "pdb": res.get("pdb", ""),
            "uniprot": str(g.get("best_acc") or g.get("homolog_acc") or ""),
            "error": res.get("error", ""), "evidence": evidence,
            "limitation": "Predicción computacional; no cierra la brecha sin validación/experimento.",
        })
        done += 1
    return rows
