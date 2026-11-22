#!/usr/bin/env python3
"""01 - Construye las tablas de referencia de la capa T1 de transporte.

Entradas (solo lectura):
  - data/db/transport.fasta          (18.864 proteínas Swiss-Prot)
  - data/metabolism/rhea_reactions.pkl
  - data/metabolism/rhea2uniprot_sprot.tsv
  - data/taxdump/nodes.dmp

Salidas (escritura solo en transport_t1/ref/):
  - transport_db.tsv          acc, entry, organism, taxid, superkingdom, description
  - chebi_transport_rhea.tsv  chebi, rhea  (solo transporte exacto del metabolito)
  - rhea2acc_transport.tsv    rhea, acc, superkingdom, organism
  - ref_stats.json            conteos y verificación de parseo

Verificación de parseo: el conjunto de accs de la base debe coincidir
exactamente con el conjunto de accs ligados a alguna Rhea de transporte en
`rhea2uniprot_sprot.tsv`. Si no coincide, el script aborta.
"""
from __future__ import annotations

import importlib
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
common = importlib.import_module("00_common")


def main() -> int:
    common.REF.mkdir(parents=True, exist_ok=True)

    print("[01] cargando Rhea ...", flush=True)
    rx, strict, one_side = common.load_rhea()
    transport_rheas = {int(k) for k, v in rx.items() if v.get("transport")}
    print(
        f"[01] Rhea: {len(rx)} reacciones, {len(transport_rheas)} con flag transport, "
        f"{len(strict)} metabolitos con transporte exacto",
        flush=True,
    )

    print("[01] parseando transport.fasta ...", flush=True)
    records = common.parse_transport_fasta()
    db_accs = {r["acc"] for r in records}
    print(f"[01] base: {len(records)} entradas, {len(db_accs)} accs únicas", flush=True)

    print("[01] resolviendo taxonomía NCBI ...", flush=True)
    parent = common.load_taxonomy_parents()
    for rec in records:
        rec["superkingdom"] = common.superkingdom(rec["taxid"], parent)
    sk_counts = Counter(rec["superkingdom"] for rec in records)

    print("[01] cargando rhea2uniprot_sprot.tsv ...", flush=True)
    rhea2acc = common.load_rhea2acc(common.RHEA2ACC, transport_rheas)
    linked_accs = set().union(*rhea2acc.values()) if rhea2acc else set()

    # --- verificación de parseo -------------------------------------------------
    if db_accs != linked_accs:
        faltan = sorted(linked_accs - db_accs)[:10]
        sobran = sorted(db_accs - linked_accs)[:10]
        raise SystemExit(
            "ERROR de parseo: base y rhea2uniprot no coinciden. "
            f"faltan_en_base={faltan} sobran_en_base={sobran}"
        )
    print(
        f"[01] OK: {len(db_accs)} accs de la base == accs ligados a transporte; "
        f"pares rhea-acc={sum(len(v) for v in rhea2acc.values())}",
        flush=True,
    )

    # --- salidas ---------------------------------------------------------------
    common.write_tsv(
        common.REF / "transport_db.tsv",
        [
            "acc",
            "entry",
            "organism",
            "taxid",
            "superkingdom",
            "description",
        ],
        (
            [
                r["acc"],
                r["entry"],
                common.clean(r["organism"]),
                r["taxid"],
                r["superkingdom"],
                common.clean(r["description"]),
            ]
            for r in sorted(records, key=lambda x: x["acc"])
        ),
    )

    chebi_rows = [
        [chebi, rid]
        for chebi in sorted(strict)
        for rid in sorted(strict[chebi])
    ]
    common.write_tsv(
        common.REF / "chebi_transport_rhea.tsv",
        ["chebi", "rhea"],
        chebi_rows,
    )

    acc_meta = {r["acc"]: r for r in records}
    r2a_rows = []
    for rid in sorted(rhea2acc):
        for acc in sorted(rhea2acc[rid]):
            meta = acc_meta[acc]
            r2a_rows.append([rid, acc, meta["superkingdom"], common.clean(meta["organism"])])
    common.write_tsv(
        common.REF / "rhea2acc_transport.tsv",
        ["rhea", "acc", "superkingdom", "organism"],
        r2a_rows,
    )

    rheas_with_acc = len(rhea2acc)
    stats = {
        "generated_by": "scripts/metabolism/01_build_reference.py",
        "inputs_sha256": {
            "data/db/transport.fasta": common.sha256_file(common.TRANSPORT_FASTA),
            "data/metabolism/rhea_reactions.pkl": common.sha256_file(common.RHEA_PKL),
            "data/metabolism/rhea2uniprot_sprot.tsv": common.sha256_file(common.RHEA2ACC),
            "data/taxdump/nodes.dmp": common.sha256_file(common.TAXDUMP / "nodes.dmp"),
        },
        "transport_db_entries": len(records),
        "transport_db_unique_accs": len(db_accs),
        "superkingdom": dict(sorted(sk_counts.items())),
        "rhea_total": len(rx),
        "rhea_transport_flag": len(transport_rheas),
        "rhea_transport_exact_substrate": len(
            {rid for rheas in strict.values() for rid in rheas}
        ),
        "metabolites_with_exact_transport": len(strict),
        "metabolites_one_side_only": len(one_side),
        "rhea_transport_with_swissprot": rheas_with_acc,
        "rhea_transport_without_swissprot": len(transport_rheas) - rheas_with_acc,
        "rhea_acc_pairs": sum(len(v) for v in rhea2acc.values()),
        "accs_linked_to_transport": len(linked_accs),
        "parse_check": {
            "db_accs_equal_rhea2uniprot_accs": True,
            "n_chebi_rhea_pairs": len(chebi_rows),
            "n_rhea_acc_pairs": len(r2a_rows),
        },
    }
    common.write_json(common.REF / "ref_stats.json", stats)

    print(
        f"[01] escrito: transport_db.tsv ({len(records)}), "
        f"chebi_transport_rhea.tsv ({len(chebi_rows)}), "
        f"rhea2acc_transport.tsv ({len(r2a_rows)})",
        flush=True,
    )
    print(f"[01] superkingdom: {dict(sorted(sk_counts.items()))}", flush=True)
    print(
        f"[01] Rhea de transporte con Swiss-Prot: {rheas_with_acc}/{len(transport_rheas)} "
        f"({100 * rheas_with_acc / len(transport_rheas):.1f}%)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
