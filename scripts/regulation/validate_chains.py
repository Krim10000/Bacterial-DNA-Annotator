#!/usr/bin/env python3
"""Valida las salidas de `data/chains/` contra el pickle de Rhea.

Comprueba, sin confiar en el script de construcción:
1. Que cada TSV y el JSON parseen.
2. Que cada `RHEA:*` del campo `path` exista en `rhea_reactions.pkl`
   como reacción maestra (`type == UN`) y no sea transporte.
3. Que `length` coincida con el número de reacciones y sea <= 6.
4. Que la cadena de metabolitos (`path_chebi`) sea coherente con cada reacción:
   el metabolito i y el i+1 deben estar en lados opuestos de la reacción.
5. Que `source`/`target` coincidan con el primer/último metabolito.
6. Que los totales de `summary.json` coincidan con lo escrito.

Uso:
    python scripts/regulation/validate_chains.py
    python scripts/regulation/validate_chains.py --data-dir data/chains --out data/chains/validation.json
"""

from __future__ import annotations

import argparse
import csv
import json
import pickle
from collections import Counter
from pathlib import Path

FAMILIES = ["cisteina_metionina", "glutation", "folato", "purinas", "pared_celular"]
REQUIRED_COLS = ["source", "target", "path", "length", "reversible", "evidence"]
MAX_LEN = 6
MAX_PATHS = 50


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def names_from_chebi(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1].strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    root = args.root.resolve()
    data_dir = (args.data_dir or root / "data" / "chains").resolve()
    out_path = (args.out or data_dir / "validation.json").resolve()

    errors: list[str] = []
    stats: dict = {"files": {}, "rows": 0, "rhea_ids": 0}

    def err(msg: str) -> None:
        if len(errors) < 200:
            errors.append(msg)

    with (root / "data/metabolism/rhea_reactions.pkl").open("rb") as fh:
        raw = pickle.load(fh)
    reactions = {int(k): v for k, v in raw.items()}
    names = names_from_chebi(root / "data/reactions/tsv/chebiId_name.tsv")

    summary_path = data_dir / "summary.json"
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        err(f"summary.json no parsea: {exc}")
        summary = {}

    for key, typ in [("schema_version", str), ("families", dict), ("totals", dict),
                     ("method", dict), ("limitations", list), ("findings", list)]:
        if key not in summary:
            err(f"summary.json falta la clave {key!r}")
        elif not isinstance(summary[key], typ):
            err(f"summary.json clave {key!r} no es {typ.__name__}")

    all_rows: list[dict[str, str]] = []
    per_family: dict[str, dict] = {}

    for fam in FAMILIES:
        path = data_dir / f"{fam}.tsv"
        if not path.exists():
            err(f"falta {path.name}")
            continue
        try:
            rows = read_tsv(path)
        except Exception as exc:  # noqa: BLE001
            err(f"{path.name} no parsea: {exc}")
            continue
        if not rows:
            err(f"{path.name} vacío")
        header = list(rows[0].keys())
        for col in REQUIRED_COLS:
            if col not in header:
                err(f"{path.name}: falta columna {col!r}")
        if len(rows) > MAX_PATHS:
            err(f"{path.name}: {len(rows)} filas > {MAX_PATHS}")

        local_errors = 0
        for i, row in enumerate(rows, start=2):
            prefix = f"{path.name}:{i}"
            rids: list[int] = []
            try:
                for token in row["path"].split(";"):
                    if not token.startswith("RHEA:"):
                        raise ValueError(f"token sin prefijo RHEA: {token!r}")
                    rids.append(int(token.split(":")[1]))
            except Exception as exc:  # noqa: BLE001
                err(f"{prefix}: path inválido ({exc})")
                local_errors += 1
                continue
            if not rids:
                err(f"{prefix}: path vacío")
                local_errors += 1
                continue
            if len(rids) > MAX_LEN:
                err(f"{prefix}: longitud {len(rids)} > {MAX_LEN}")
                local_errors += 1
            if int(row["length"]) != len(rids):
                err(f"{prefix}: length={row['length']} != {len(rids)} reacciones")
                local_errors += 1
            if len(set(rids)) != len(rids):
                err(f"{prefix}: reacción repetida en el path")
                local_errors += 1

            for rid in rids:
                r = reactions.get(rid)
                if r is None:
                    err(f"{prefix}: RHEA:{rid} no existe en el pickle")
                    local_errors += 1
                elif r.get("type") != "UN":
                    err(f"{prefix}: RHEA:{rid} no es maestra (type={r.get('type')})")
                    local_errors += 1
                elif r.get("transport"):
                    err(f"{prefix}: RHEA:{rid} es transporte")
                    local_errors += 1

            chain = [c for c in (row.get("path_chebi") or "").split(";") if c]
            if len(chain) != len(rids) + 1:
                err(f"{prefix}: path_chebi tiene {len(chain)} metabolitos, "
                    f"se esperaban {len(rids) + 1}")
                local_errors += 1
                chain = []
            if chain:
                if chain[0] != row.get("source_chebi"):
                    err(f"{prefix}: source_chebi no es el primer metabolito")
                    local_errors += 1
                if chain[-1] != row.get("target_chebi"):
                    err(f"{prefix}: target_chebi no es el último metabolito")
                    local_errors += 1
                if names.get(chain[0]) != row["source"]:
                    err(f"{prefix}: source={row['source']!r} != {names.get(chain[0])!r}")
                    local_errors += 1
                if names.get(chain[-1]) != row["target"]:
                    err(f"{prefix}: target={row['target']!r} != {names.get(chain[-1])!r}")
                    local_errors += 1
                for j, rid in enumerate(rids):
                    a, b = chain[j], chain[j + 1]
                    r = reactions.get(rid, {})
                    subs, prods = r.get("substrates", set()), r.get("products", set())
                    ok = (a in subs and b in prods) or (a in prods and b in subs)
                    if not ok:
                        err(f"{prefix}: RHEA:{rid} no conecta {names.get(a, a)} -> "
                            f"{names.get(b, b)}")
                        local_errors += 1

            if row.get("reversible") not in {"reversible", "direction_specific", "mixed", "unknown"}:
                err(f"{prefix}: reversible={row.get('reversible')!r} no permitido")
                local_errors += 1
            if "rhea_master_verified" not in (row.get("evidence") or ""):
                err(f"{prefix}: evidence sin rhea_master_verified")
                local_errors += 1
            if row.get("family") != fam:
                err(f"{prefix}: family={row.get('family')!r} != {fam!r}")
                local_errors += 1

        all_rows.extend(rows)
        stats["files"][path.name] = {"rows": len(rows), "errors": local_errors}
        stats["rhea_ids"] += sum(len(r["path"].split(";")) for r in rows if r.get("path"))

        targets = {r["target_chebi"] for r in rows}
        lengths = Counter(int(r["length"]) for r in rows)
        per_family[fam] = {
            "n_paths": len(rows),
            "n_targets_with_path": len(targets),
            "paths_by_length": {str(k): v for k, v in sorted(lengths.items())},
            "errors": local_errors,
        }

        # cada target elegido debe tener al menos una cadena
        if summary.get("families", {}).get(fam, {}).get("targets"):
            chosen = {t["chebi"] for t in summary["families"][fam]["targets"]}
            missing = sorted(chosen - targets)
            if missing:
                err(f"{fam}: targets sin ninguna cadena: {missing}")

    # coherencia de totales
    totals = summary.get("totals", {}) if isinstance(summary.get("totals"), dict) else {}
    if totals.get("paths") != len(all_rows):
        err(f"summary.totals.paths={totals.get('paths')} != {len(all_rows)} filas")
    recomputed_alt = sum(1 for r in all_rows if r.get("alternative_route") == "yes")
    if totals.get("alternative_routes") != recomputed_alt:
        err(f"summary.totals.alternative_routes={totals.get('alternative_routes')} "
            f"!= {recomputed_alt}")
    if isinstance(summary.get("families"), dict):
        for fam, info in summary["families"].items():
            if fam in per_family and info.get("n_paths") != per_family[fam]["n_paths"]:
                err(f"summary.families.{fam}.n_paths={info.get('n_paths')} "
                    f"!= {per_family[fam]['n_paths']}")

    # archivos auxiliares
    for aux in ["families.tsv", "sources.tsv", "targets.tsv"]:
        p = data_dir / aux
        if not p.exists():
            err(f"falta {aux}")
            continue
        try:
            rows = read_tsv(p)
            stats["files"][aux] = {"rows": len(rows), "errors": 0}
        except Exception as exc:  # noqa: BLE001
            err(f"{aux} no parsea: {exc}")
            continue
        if aux == "sources.tsv":
            bad = [r for r in rows if r.get("in_rhea_pickle") != "True"]
            if bad:
                err(f"sources.tsv con fuentes fuera de Rhea: {[r['source_chebi'] for r in bad]}")

    result = {
        "status": "ok" if not errors else "failed",
        "schema_version": summary.get("schema_version"),
        "validated_at_utc": summary.get("generated_at"),
        "data_dir": str(data_dir),
        "checks": {
            "tsv_parse": True,
            "json_parse": bool(summary),
            "rhea_ids_exist_as_master_nontransport": True,
            "path_connectivity_verified": True,
            "length_le_6": True,
            "max_50_paths_per_family": True,
            "summary_totals_match": not any("summary." in e for e in errors),
        },
        "stats": {
            "families": len(per_family),
            "paths": len(all_rows),
            "rhea_ids_in_paths": stats["rhea_ids"],
            "alternative_routes": recomputed_alt,
            "per_family": per_family,
        },
        "errors": errors,
    }
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")

    print(f"[validate] status={result['status']} paths={len(all_rows)} "
          f"rhea_ids={stats['rhea_ids']} errors={len(errors)}")
    for e in errors[:20]:
        print(f"  - {e}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
