"""Tests de clasificación y auditoría de cierre de brechas.

No descargan nada, no usan GPU y no tocan la red. Los fixtures son sintéticos:
el objetivo es fijar el criterio, no recalcular la instantánea real.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import pickle
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "audit_gap_closure.py"


def _load_tool():
    spec = importlib.util.spec_from_file_location("audit_gap_closure", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = _load_tool()

TRANSFER_FIELDS = [
    "complex_id", "complex_name", "n_subunits", "n_mapped", "frac_mapped",
    "identity_mean", "score", "genes_k12", "genes_genoma", "observed_n",
    "expected_n", "stoich_status", "mass_kda", "phylo_jacc", "ecs",
    "evidence_class", "eco", "limitation",
]
MODEL_FIELDS = [
    "complex_id", "name", "contig", "start", "end", "subunits", "expected",
    "missing", "status", "evidence", "confidence", "source", "ecs", "flag",
    "limitation",
]


def gap(**overrides):
    """Una fila mínima de ``gap_interpretation.tsv``."""
    row = {field: "" for field in mod.GAP_FIELDS}
    row.update({
        "genome": "g1",
        "gap": "metabolite",
        "chebi": "CHEBI:1",
        "status": "validated",
        "role": "anabolic",
    })
    row.update(overrides)
    return row


def association(**overrides):
    """Una asociación ya construida (lo que produce ``build_associations``)."""
    record = {
        "genome": "g1",
        "gap": "metabolite",
        "chebi": "CHEBI:1",
        "complex_id": "CPX-1",
        "complex_name": "test complex",
        "source": "transfer",
        "strength": "strong",
        "frac_mapped": 0.9,
        "incomplete": False,
        "stoich_status": "unknown",
        "ecs": ["1.1.1.1"],
        "rhea_overlap": False,
        "evidence_class": "documentado_inferido",
        "limitation": "no cierra brechas",
    }
    record.update(overrides)
    return record


def complex_record(**overrides):
    record = {
        "genome": "g1",
        "complex_id": "CPX-1",
        "complex_name": "test complex",
        "source": "transfer",
        "strength": "strong",
        "frac_mapped": 0.9,
        "incomplete": False,
        "stoich_status": "unknown",
        "expected_n": 2,
        "observed_n": 2,
        "ecs": ["1.1.1.1"],
        "products": {"CHEBI:1"},
        "rheas": {1},
        "evidence_class": "documentado_inferido",
        "limitation": "no cierra brechas",
    }
    record.update(overrides)
    return record


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def make_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Crea un ``data/`` mínimo con brechas, complejos y puente EC/Rhea."""
    data = tmp_path / "data"
    metabolism = data / "metabolism"
    metabolism.mkdir(parents=True)
    (data / "reactions" / "tsv").mkdir(parents=True)
    complexes = data / "complexes"
    complexes.mkdir(parents=True)

    reactions = {
        1: {"products": {"CHEBI:3"}, "substrates": {"CHEBI:30"}},
        2: {"products": {"CHEBI:4"}, "substrates": {"CHEBI:40"}},
        3: {"products": {"CHEBI:1"}, "substrates": {"CHEBI:10"}},
    }
    with (metabolism / "rhea_reactions.pkl").open("wb") as fh:
        pickle.dump(reactions, fh)
    (data / "reactions" / "tsv" / "rhea2ec.tsv").write_text(
        "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n"
        "1\tUN\t1\t1.1.1.1\n"
        "2\tUN\t2\t1.1.1.2\n"
        "3\tUN\t3\t1.1.1.3\n",
        encoding="utf-8",
    )

    gaps = [
        gap(genome="g1", gap="strong", chebi="CHEBI:1",
            mapping_validated="1", mapping_rhea="RHEA:3"),
        gap(genome="g1", gap="weak", chebi="CHEBI:2",
            hmm_evidence="PF00001"),
        gap(genome="g1", gap="complex", chebi="CHEBI:3",
            status="anomalous"),
        gap(genome="g1", gap="incomplete", chebi="CHEBI:4",
            status="anomalous"),
        gap(genome="g1", gap="direction", chebi="CHEBI:5",
            audit_status="enzymatic_missing",
            direction_status="direction_unspecified"),
        gap(genome="g1", gap="nothing", chebi="CHEBI:6",
            status="", role="", origin="", audit_status="",
            transport_status="", direction_status="", review_status="",
            hmm_evidence="", homolog_status="", mapping_validated=""),
    ]
    gap_path = metabolism / "gap_interpretation.tsv"
    write_tsv(gap_path, mod.GAP_FIELDS, gaps)

    transfer_rows = [
        {
            "complex_id": "CPX-1", "complex_name": "test complex",
            "n_subunits": 2, "n_mapped": 2, "frac_mapped": 0.9,
            "identity_mean": 70, "score": 0.9, "observed_n": 2,
            "expected_n": 2, "stoich_status": "unknown", "ecs": "1.1.1.1",
            "evidence_class": "documentado_inferido",
            "limitation": "no cierra brechas",
        },
        {
            "complex_id": "CPX-2", "complex_name": "incomplete complex",
            "n_subunits": 2, "n_mapped": 1, "frac_mapped": 0.9,
            "identity_mean": 65, "score": 0.85, "observed_n": 1,
            "expected_n": 2, "stoich_status": "incomplete", "ecs": "1.1.1.2",
            "evidence_class": "documentado_inferido",
            "limitation": "no cierra brechas",
        },
    ]
    transfer_path = complexes / "transfer_g1.tsv"
    write_tsv(transfer_path, TRANSFER_FIELDS, transfer_rows)

    write_tsv(
        complexes / "model_g1.complexes.tsv",
        MODEL_FIELDS,
        [{
            "complex_id": "cx1", "name": "model candidate",
            "status": "complex_candidate", "ecs": "", "missing": "",
            "limitation": "no cierra brechas",
        }],
    )
    return data, gap_path, transfer_path


# ---------------------------------------------------------------------------
# Evidencia directa
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "updates",
    [
        {"mapping_validated": "1", "mapping_rhea": "RHEA:1"},
        {"homolog_status": "homolog_strong", "homolog_acc": "P12345"},
        {"review_status": "domain_confirmed", "hmm_evidence": "PF00001"},
        {"review_status": "reannotated_hmm", "hmm_evidence": "PF00002"},
    ],
)
def test_strong_evidence_is_supported_not_confirmed(updates):
    out = mod.classify_gap(gap(**updates))
    assert out["state"] == "gap_evidence_supported"
    assert out["evidence_tier"] == "strong"
    assert out["closed"] is False
    assert out["state"] != "confirmed"


@pytest.mark.parametrize(
    "updates",
    [
        {"hmm_evidence": "PF00001"},
        {"homolog_status": "homolog_weak", "homolog_acc": "P12345"},
        {"homolog_status": "homolog_low", "homolog_acc": "P12345"},
        {"status": "candidate"},
        {"review_status": "evidence_high"},
        {"origin": "promiscuous", "multifunction_ec": "1.1.1.1"},
    ],
)
def test_weak_evidence_is_candidate(updates):
    out = mod.classify_gap(gap(**updates))
    assert out["state"] == "gap_candidate"
    assert out["evidence_tier"] == "weak"
    assert out["closed"] is False


def test_missing_acc_does_not_count_as_homolog_evidence():
    out = mod.classify_gap(gap(homolog_status="homolog_strong", homolog_acc=""))
    assert out["state"] == "not_evaluated"
    assert out["evidence_tier"] == "none"


# ---------------------------------------------------------------------------
# Complejos
# ---------------------------------------------------------------------------
def test_complex_strong_is_suspected_not_closed():
    out = mod.classify_gap(gap(), [association()])
    assert out["state"] == "complex_suspected"
    assert out["associated"] is True
    assert out["closed"] is False


def test_complex_incomplete_only_is_incomplete():
    out = mod.classify_gap(
        gap(),
        [association(incomplete=True, stoich_status="incomplete")],
    )
    assert out["state"] == "complex_incomplete"
    assert out["closed"] is False


def test_strong_complex_trumps_incomplete_complex():
    out = mod.classify_gap(
        gap(),
        [
            association(complex_id="CPX-STRONG"),
            association(complex_id="CPX-INCOMPLETE", incomplete=True),
        ],
    )
    assert out["state"] == "complex_suspected"
    assert out["closed"] is False


def test_weak_complex_is_only_a_candidate():
    out = mod.classify_gap(
        gap(),
        [association(strength="weak", frac_mapped=0.5)],
    )
    assert out["state"] == "gap_candidate"
    assert out["closed"] is False


def test_direct_evidence_trumps_complex_association():
    out = mod.classify_gap(gap(mapping_validated="1"), [association()])
    assert out["state"] == "gap_evidence_supported"


# ---------------------------------------------------------------------------
# Bloqueadores
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "updates,code",
    [
        ({"direction_status": "direction_unspecified"}, "direction_unresolved"),
        ({"direction_status": "thermo_unknown"}, "direction_unresolved"),
        ({"role": "catabolic"}, "catabolic"),
        ({"audit_status": "degradation"}, "catabolic"),
        ({"origin": "degradation"}, "catabolic"),
        ({"audit_status": "false_positive"}, "false_positive_direction"),
        ({"role": "entry_point"}, "entry_point"),
        ({"origin": "imported_or_entry"}, "entry_point"),
        ({"origin": "spontaneous"}, "spontaneous"),
        ({"status": "lineage_specific"}, "host_dependent"),
        ({"homolog_status": "no_homolog"}, "no_homolog"),
        ({"status": "novel"}, "novel_no_reference"),
        ({"transport_status": "transport_reaction_only"},
         "transport_reaction_only"),
        ({"transport_status": "t1_bacterial"}, "transporter_homolog"),
    ],
)
def test_blockers_are_not_closable_cpu(updates, code):
    out = mod.classify_gap(gap(**updates))
    assert out["state"] == "not_closable_cpu"
    assert code in out["blockers"]
    assert out["closed"] is False


def test_without_layers_is_not_evaluated():
    out = mod.classify_gap(gap(
        status="", role="", origin="", audit_status="", transport_status="",
        direction_status="", review_status="", hmm_evidence="",
        homolog_status="", mapping_validated="",
    ))
    assert out["state"] == "not_evaluated"
    assert out["blockers"] == []
    assert out["closed"] is False


def test_one_sided_transport_review_is_not_a_closure():
    out = mod.classify_gap(gap(
        review_status="one_sided_not_transport",
        transport_status="no_transport",
    ))
    assert out["state"] == "not_evaluated"
    assert out["closed"] is False


def test_confirmed_input_is_never_propagated():
    out = mod.classify_gap(gap(status="confirmed", review_status="confirmed"))
    assert out["state"] in mod.ALLOWED_STATES
    assert out["state"] != "confirmed"
    assert "confirmed" not in mod.ALLOWED_STATES
    assert out["closed"] is False


# ---------------------------------------------------------------------------
# Asociación por ChEBI, genoma y Rhea
# ---------------------------------------------------------------------------
def test_association_requires_same_genome():
    gaps = [gap(genome="g1", chebi="CHEBI:1")]
    complexes = {"g2": [complex_record(genome="g2", products={"CHEBI:1"})]}
    associations = mod.build_associations(gaps, complexes)
    assert associations == {}


def test_association_requires_exact_chebi():
    gaps = [gap(genome="g1", chebi="CHEBI:2")]
    complexes = {"g1": [complex_record(products={"CHEBI:1"})]}
    associations = mod.build_associations(gaps, complexes)
    assert associations == {}


def test_association_records_rhea_overlap_without_inventing_products():
    gaps = [gap(genome="g1", chebi="CHEBI:1", mapping_rhea="RHEA:10")]
    complexes = {"g1": [complex_record(products={"CHEBI:1"}, rheas={10})]}
    associations = mod.build_associations(gaps, complexes)
    assert len(associations[("g1", "CHEBI:1")]) == 1
    assert associations[("g1", "CHEBI:1")][0]["rhea_overlap"] is True

    # Solapar Rhea no basta si el producto no coincide con el ChEBI.
    complexes = {"g1": [complex_record(products={"CHEBI:999"}, rheas={10})]}
    assert mod.build_associations(gaps, complexes) == {}


# ---------------------------------------------------------------------------
# Extremo a extremo sobre un data/ sintético
# ---------------------------------------------------------------------------
def test_audit_end_to_end_is_read_only_and_never_closes(tmp_path):
    data, gap_path, transfer_path = make_fixture(tmp_path)
    gap_before = gap_path.read_bytes()
    transfer_before = transfer_path.read_bytes()

    out_md = tmp_path / "reports" / "audit.md"
    out_json = tmp_path / "reports" / "audit.json"
    result = mod.audit(data, out_md, out_json)

    assert result["totals"]["gaps"] == 6
    assert result["totals"]["states_sum"] == 6
    assert result["states"] == {
        "gap_evidence_supported": 1,
        "gap_candidate": 1,
        "complex_suspected": 1,
        "complex_incomplete": 1,
        "not_evaluated": 1,
        "not_closable_cpu": 1,
    }
    assert result["complexes"]["transfer"]["files"] == 1
    assert result["complexes"]["transfer"]["rows"] == 2
    assert result["complexes"]["transfer"]["rows_with_ecs"] == 2
    assert result["complexes"]["transfer"]["rows_with_products"] == 2
    assert result["complexes"]["transfer"]["strong"] == 2
    assert result["complexes"]["transfer"]["incomplete"] == 1
    assert result["complexes"]["model"]["rows"] == 1
    assert result["complexes"]["model"]["rows_with_ecs"] == 0
    assert result["complexes"]["associations"]["gaps_with_any"] == 2
    assert result["complexes"]["associations"]["gaps_with_strong"] == 1
    assert result["complexes"]["associations"]["gaps_with_incomplete"] == 1
    assert result["closure"]["confirmed"] == 0
    assert result["closure"]["closed"] == 0
    assert result["closure"]["associated_not_closed"] == 2
    assert len(result["not_evaluated_detail"]) == 1
    assert result["not_evaluated_detail"][0]["gap"] == "nothing"
    assert "confirmed" not in result["states"]
    assert result["states"]["complex_incomplete"] == 1

    assert out_md.exists() and out_json.exists()
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    assert payload["closure"]["confirmed"] == 0
    assert payload["states"] == result["states"]
    assert "Asociada no es cerrada" in out_md.read_text(encoding="utf-8")

    assert gap_path.read_bytes() == gap_before
    assert transfer_path.read_bytes() == transfer_before


def test_audit_rejects_missing_gap_table(tmp_path):
    empty = tmp_path / "data"
    (empty / "metabolism").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        mod.audit(empty, tmp_path / "a.md", tmp_path / "a.json")
