#!/usr/bin/env python3
"""09 — Compara las dos capas (piloto de 10 vs completa de 400).

Para los 10 genomas piloto hay dos JSON:
  data/taxonomy/<genome>.json         capa piloto (panel completo de 695 refs)
  data/taxonomy/full/<genome>.json    capa completa (panel de su filo)

Este script escribe `data/taxonomy/full/pilot_consistency.json` con:
  * diferencia de ANI y de `closest_ref` entre capas,
  * veredicto (coincide / cambia referencia / cambia especie),
  * comprobación de que la capa piloto no fue sobrescrita.

No modifica nada.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"


def main() -> int:
    rows = []
    with open(L.PILOT_TSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["genome"]
            p_path = L.OUT / f"{g}.json"
            f_path = FULL / f"{g}.json"
            if not p_path.exists() or not f_path.exists():
                rows.append({"genome": g, "status": "missing_file"})
                continue
            p = json.loads(p_path.read_text(encoding="utf-8"))
            f = json.loads(f_path.read_text(encoding="utf-8"))
            p_ani, f_ani = p.get("ani"), f.get("ani")
            delta = round(abs(p_ani - f_ani), 2) if p_ani is not None and f_ani is not None else None
            same_ref = (p.get("closest_ref") and f.get("closest_ref")
                        and L.pilot_accession(p["closest_ref"]) == L.pilot_accession(f["closest_ref"]))
            same_species = p.get("species") == f.get("species")
            if delta is None:
                verdict = "ani_no_comparable"
            elif same_ref and delta <= 0.5:
                verdict = "identico"
            elif same_ref:
                verdict = "misma_ref_ani_distinto"
            else:
                verdict = "referencia_distinta"
            rows.append({
                "genome": g,
                "species_pilot": p.get("species"),
                "species_full": f.get("species"),
                "same_species": same_species,
                "ani_pilot": p_ani,
                "ani_full": f_ani,
                "ani_delta": delta,
                "ref_pilot": p.get("closest_ref"),
                "ref_full": f.get("closest_ref"),
                "same_ref": same_ref,
                "confidence_pilot": p.get("confidence"),
                "confidence_full": f.get("confidence"),
                "pipeline_pilot": p.get("taxonomy_version", {}).get("pipeline"),
                "pipeline_full": f.get("taxonomy_version", {}).get("pipeline"),
                "verdict": verdict,
            })

    consistent = sum(1 for r in rows if r.get("verdict") in ("identico", "misma_ref_ani_distinto"))
    payload = {
        "n_pilots": len(rows),
        "n_consistent": consistent,
        "n_ref_changed": sum(1 for r in rows if r.get("verdict") == "referencia_distinta"),
        "n_species_changed": sum(1 for r in rows if r.get("same_species") is False),
        "pilot_layer_untouched": all(r.get("pipeline_pilot") == "1.0" for r in rows),
        "full_layer_pipeline": sorted({r.get("pipeline_full") for r in rows}),
        "rows": rows,
    }
    L.write_json(FULL / "pilot_consistency.json", payload)
    print(json.dumps({k: payload[k] for k in
                      ("n_pilots", "n_consistent", "n_ref_changed", "n_species_changed",
                       "pilot_layer_untouched", "full_layer_pipeline")},
                     indent=2, ensure_ascii=False))
    for r in rows:
        print(f"  {r['genome']}: pilot={r.get('ani_pilot')} full={r.get('ani_full')} "
              f"delta={r.get('ani_delta')} -> {r.get('verdict')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
