#!/usr/bin/env python3
"""02 - Selecciona los 10 genomas piloto.

Regla objetiva y reproducible: **un genoma por filo**, el de mayor carga T2
estricta del filo (gaps con reacción Rhea de transporte exacta del metabolito),
desempatando por (T2 con candidato Swiss-Prot, gaps totales, genome id). Se
eligen los 10 filos con mayor carga T2. Solo se consideran genomas con
`proteins.faa` disponible.

Entradas: data/metabolism/transport_audit.tsv, data/metabolism/gap_interpretation.tsv,
          ref/chebi_transport_rhea.tsv, ref/rhea2acc_transport.tsv
Salida:   data/metabolism/transport_t1/pilot_genomes.tsv (+ pilot_selection.json)
"""
from __future__ import annotations

import importlib
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
common = importlib.import_module("00_common")

PILOT_N = 10


def main() -> int:
    _, chebi2rhea, rhea2acc, _ = common.load_ref()

    total = defaultdict(int)
    strict = defaultdict(int)
    testable = defaultdict(int)
    for row in common.read_tsv(common.TRANSPORT_AUDIT):
        genome = row["genome"]
        total[genome] += 1
        rheas = chebi2rhea.get(row["chebi"], [])
        if rheas:
            strict[genome] += 1
            if any(rhea2acc.get(r) for r in rheas):
                testable[genome] += 1

    org, phylum = {}, {}
    for row in common.read_tsv(common.GAP_TABLE):
        org.setdefault(row["genome"], common.clean(row.get("organism", "")))
        phylum.setdefault(row["genome"], common.clean(row.get("phylum", "")))

    # --- un ganador por filo ----------------------------------------------------
    best: dict[str, str] = {}
    skipped = defaultdict(int)
    for genome in total:
        phy = phylum.get(genome, "")
        proteins_path = common.ANN / genome / "proteins.faa"
        if not phy or not proteins_path.exists():
            skipped["sin_filo_o_sin_proteoma"] += 1
            continue
        key = (strict[genome], testable[genome], total[genome], genome)
        if phy not in best:
            best[phy] = genome
            continue
        cur = best[phy]
        cur_key = (strict[cur], testable[cur], total[cur], cur)
        if key > cur_key:
            best[phy] = genome

    candidates = []
    for phy, genome in best.items():
        n_prot = len(common.fasta_ids(common.ANN / genome / "proteins.faa"))
        if n_prot == 0:
            skipped["proteoma_vacio"] += 1
            continue
        candidates.append(
            {
                "genome": genome,
                "phylum": phy,
                "organism": org.get(genome, ""),
                "gaps_total": total[genome],
                "strict_t2": strict[genome],
                "testable_t2": testable[genome],
                "proteins": n_prot,
            }
        )

    candidates.sort(
        key=lambda c: (-c["strict_t2"], -c["testable_t2"], -c["gaps_total"], c["genome"])
    )
    if len(candidates) < PILOT_N:
        raise SystemExit(
            f"ERROR: solo {len(candidates)} filos con datos; se necesitan {PILOT_N}"
        )
    pilot = candidates[:PILOT_N]

    for c in pilot:
        c["rationale"] = (
            f"mayor carga T2 del filo {c['phylum']}: {c['strict_t2']} gaps con "
            f"reacción de transporte exacta, {c['testable_t2']} con candidato Swiss-Prot"
        )

    common.write_tsv(
        common.OUT / "pilot_genomes.tsv",
        [
            "genome",
            "phylum",
            "organism",
            "gaps_total",
            "strict_t2",
            "testable_t2",
            "proteins",
            "rationale",
        ],
        (
            [
                c["genome"],
                c["phylum"],
                c["organism"],
                c["gaps_total"],
                c["strict_t2"],
                c["testable_t2"],
                c["proteins"],
                c["rationale"],
            ]
            for c in pilot
        ),
    )
    common.write_json(
        common.OUT / "pilot_selection.json",
        {
            "rule": (
                "un genoma por filo (mayor T2 estricta; desempate por T2 con "
                "candidato, gaps totales y genome id); se eligen los 10 filos "
                "con mayor carga T2 entre genomas con proteins.faa"
            ),
            "pilot_n": PILOT_N,
            "phylum_winners_total": len(best),
            "genomes_with_gaps": len(total),
            "skipped": dict(skipped),
            "selected": pilot,
            "not_selected": candidates[PILOT_N:],
        },
    )

    print(f"[02] pilotos: {len(pilot)} (de {len(candidates)} filos candidatos)")
    for c in pilot:
        print(
            f"  {c['genome']}\t{c['phylum']}\tT2={c['strict_t2']:2d} "
            f"testable={c['testable_t2']:2d} gaps={c['gaps_total']:2d} "
            f"prot={c['proteins']}"
        )
    print(
        f"[02] total T2 estricta en piloto: {sum(c['strict_t2'] for c in pilot)}; "
        f"con candidato: {sum(c['testable_t2'] for c in pilot)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
