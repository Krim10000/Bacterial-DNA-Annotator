"""Sitios activos: precalculados para esenciales y bajo demanda para el resto."""
import csv, glob, json, urllib.request
from collections import defaultdict
from pathlib import Path

_CACHE = {}

def _r2a(data_root):
    if "r2a" in _CACHE: return _CACHE["r2a"]
    out = defaultdict(set)
    p = Path(data_root)/"metabolism"/"rhea2uniprot_sprot.tsv"
    if p.exists():
        for line in open(p):
            f = line.rstrip().split("\t")
            if len(f) >= 4 and f[0] != "RHEA_ID":
                out["RHEA:"+f[0]].add(f[3])
    _CACHE["r2a"] = out
    return out

def _info(data_root):
    if "info" in _CACHE: return _CACHE["info"]
    out = {}
    for f in glob.glob(str(Path(data_root)/"metabolism"/"evidence"/"work*"/"acc_info.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            out[r["acc"]] = r.get("euk", "")
    _CACHE["info"] = out
    return out

def _refg(data_root):
    if "refg" in _CACHE: return _CACHE["refg"]
    out = {}
    p = Path(data_root)/"out_fl"/"refined"/"refined_gaps.tsv"
    if p.exists():
        for r in csv.DictReader(open(p), delimiter="\t"):
            out[(r["genome"], r["gap"])] = r
    _CACHE["refg"] = out
    return out

def _real(data_root):
    if "real" in _CACHE: return _CACHE["real"]
    out = {}
    p = Path(data_root)/"metabolism"/"real_gaps.tsv"
    if p.exists():
        for r in csv.DictReader(open(p), delimiter="\t"):
            out[(r["genome"], r["gap"])] = r
    _CACHE["real"] = out
    return out

def essential_rows(data_root, map_id, gap):
    p = Path(data_root)/"metabolism"/"essential"/"ACTIVE_SITES_ESSENTIAL.tsv"
    if not p.exists(): return []
    return [r for r in csv.DictReader(open(p), delimiter="\t")
            if r["map_id"] == map_id and r["target"] == gap and r["feature"] != "none"]

def fetch_features(acc, cache_dir):
    cache_dir = Path(cache_dir); cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir/f"{acc}.json"
    if cache.exists():
        d = json.load(open(cache))
    else:
        req = urllib.request.Request(
            f"https://rest.uniprot.org/uniprotkb/{acc}.json",
            headers={"User-Agent": "dna-poc/1.0"})
        d = json.load(urllib.request.urlopen(req, timeout=20))
        json.dump(d, open(cache, "w"))
    name = ((d.get("proteinDescription", {}) or {}).get("recommendedName", {}) or {}).get("fullName", {}).get("value", "")
    ecs = ((d.get("proteinDescription", {}) or {}).get("recommendedName", {}) or {}).get("ecNumbers", []) or []
    ec = ";".join(x.get("value", "") for x in ecs)
    rows = []
    for f in d.get("features", []) or []:
        if f.get("type") not in ("Active site", "Binding site", "Metal binding", "Site"):
            continue
        loc = f.get("location", {}) or {}
        rows.append({"acc": acc, "protein": name, "ec": ec,
                     "feature": f.get("type", ""), "description": f.get("description", ""),
                     "start": (loc.get("start", {}) or {}).get("value", ""),
                     "end": (loc.get("end", {}) or {}).get("value", "")})
    return rows

def for_gap(data_root, genome, gap, fetch=True, max_accs=3):
    """Precalculado si es esencial; bajo demanda si no."""
    real = _real(data_root).get((genome, gap), {})
    paths = real.get("pathways", "")
    map_id = paths.split(" ")[0] if paths else ""
    if map_id:
        er = essential_rows(data_root, map_id, gap)
        if er:
            return {"source": "essential_precomputed", "rows": er, "map_id": map_id}
    if not fetch:
        return {"source": "not_fetched", "rows": [], "map_id": map_id}
    rr = _refg(data_root).get((genome, gap), {})
    ps = [x for x in rr.get("producers", "").split(";") if x]
    accs = []
    for rid in ps:
        accs += [a for a in _r2a(data_root).get(rid, set()) if _info(data_root).get(a) == "0"]
    accs = sorted(set(accs))[:max_accs]
    rows = []
    for acc in accs:
        try: rows += fetch_features(acc, Path(data_root)/"metabolism"/"active_sites_cache")
        except Exception: continue
    return {"source": "on_demand_uniprot", "rows": rows, "map_id": map_id}
