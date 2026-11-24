"""Pipeline headless: anota un FASTA desde cero y escribe todos los artefactos.

Uso:
    python -m webapp.pipeline genoma.fna --outdir results/
    python -m webapp.pipeline a.fna b.fna.gz --outdir results/ --threads 8
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotator import annotate, build_report  # noqa: E402
from functional import annotate as functional_annotate, write_paths  # noqa: E402
from complexes import infer as infer_complexes, enrich_gaps, enrich_rows, write as write_complexes  # noqa: E402
from complex_model import run as run_complex_model  # noqa: E402
from complex_transfer import load as load_complex_transfer  # noqa: E402
from membrane_layer import run as run_membrane  # noqa: E402
from domain_layer import run as run_domains  # noqa: E402
from structural_b import run as run_structural_b, evaluate as evaluate_structural_b  # noqa: E402
from context_layer import run as run_context  # noqa: E402
from regulation_layer import DEFAULT_MAX_GAP, run as run_regulation  # noqa: E402
import oric  # noqa: E402
from esmfold_tools import close_gaps  # noqa: E402
from structure_evidence import evidence_for_rows  # noqa: E402
from structural_annotation import scan_pdbs, summarize as summarize_enzymm  # noqa: E402
from structural_filters import run_filters  # noqa: E402
from structural_hypotheses import build as build_hypothesis, from_manifest  # noqa: E402
from foldseek_local import batch_search as batch_local_foldseek  # noqa: E402
from precomputed_hypotheses import load as load_precomputed  # noqa: E402
from metabolism import build_network  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "data"


def _num(value) -> float | None:
    """Convierte pLDDT/pTM de TSV a float sin romper si viene vacío."""
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def read_fasta_bytes(path: Path) -> bytes:
    raw = path.read_bytes()
    if path.suffix == ".gz" or raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw)
    return raw


def write_outputs(res: dict, outdir: Path, metabolic: dict | None) -> None:
    g = res["genome"]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{g}.gff3").write_text(res["gff"])
    (outdir / f"{g}.proteins.faa").write_text(res["proteins_fasta"])
    (outdir / f"{g}.candidates.tsv").write_text(res["candidates_tsv"])
    (outdir / f"{g}.report.md").write_text(res["report_md"])
    summary = {
        "name": res["name"], "genome": g, "bp": res["bp"], "gc": res["gc"],
        "coverage": res.get("coverage", 0), "annotated_bp": res.get("annotated_bp", 0),
        "features": dict(Counter(f["type"] for f in res["features"])),
        "deserts": len(res.get("deserts", [])),
        "candidates": len(res.get("candidates", [])),
    }
    if res.get("functional"):
        f = res["functional"]
        summary["functional"] = {
            "status": f.get("status"), "n_hits": f.get("n_hits"),
            "match_rate": f.get("match_rate"), "paths": f.get("paths"),
            "warnings": f.get("warnings"),
        }
    if res.get("complexes"):
        c = res["complexes"]
        summary["complexes"] = {
            "n_complete": c.get("n_complete"), "n_incomplete": c.get("n_incomplete"),
            "n_ambiguous": c.get("n_ambiguous"),
            "n_candidate": c.get("n_candidate", 0),
            "n_suspected": sum(1 for x in (c.get("complexes") or [])
                               if x.get("status") == "complex_suspected"),
            "rows": len(c.get("complexes") or []),
        }
    if res.get("context_c"):
        c = res["context_c"]
        summary["context_c"] = {"ok": c.get("ok"), "n": c.get("n", 0),
                                "reason": c.get("reason", "")}
    if res.get("regulation"):
        r = res["regulation"]
        summary["regulation"] = {
            "ok": r.get("ok"), "n_operons": r.get("n_operons", 0),
            "n_promoters": r.get("n_promoters", 0),
            "evidence": "rule_based", "status": "hypothesis_not_experimental",
            "reason": r.get("reason", ""),
        }
    if res.get("structural_b_filters"):
        bf = res["structural_b_filters"]
        summary["structural_b_filters"] = {
            "n_evaluated": bf.get("n_evaluated", 0),
            "n_no_local_db": bf.get("n_no_local_db", 0),
            "n_no_model": bf.get("n_no_model", 0),
            "network": bf.get("network", ""),
        }
    if res.get("structural_b"):
        b = res["structural_b"]
        summary["structural_b"] = {"ok": b.get("ok"),
                                   "n_available": b.get("n_available", 0),
                                   "reason": b.get("reason", "")}
    if res.get("domains"):
        d = res["domains"]
        summary["domains"] = {"ok": d.get("ok"), "n": d.get("n", 0),
                              "reason": d.get("reason", "")}
    if res.get("membrane"):
        m = res["membrane"]
        summary["membrane"] = {"ok": m.get("ok"), "rows": m.get("n_rows", 0),
                               "reason": m.get("reason", "")}
    if res.get("oric"):
        o = res["oric"]
        summary["oric"] = {k: o.get(k) for k in
                           ("status", "position", "confidence", "reason",
                            "validated_by", "chromosome", "size")}
    if metabolic:
        gaps = metabolic.get("gaps", [])
        summary["metabolic"] = {
            "source": metabolic.get("source"), "reactions": metabolic.get("reactions"),
            "proteins": metabolic.get("proteins"), "gaps": len(gaps),
            "origins": dict(Counter(x.get("origin", "") for x in gaps)),
            "transport": dict(Counter(x.get("transport_status", "") for x in gaps)),
            "direction": dict(Counter(x.get("direction_status", "") for x in gaps)),
        }
        cols = ["gap", "chebi", "support", "role", "origin", "status", "status_reason",
                "transport_status", "transport_rhea", "transporter_acc", "transport_pident",
                "transport_qcov", "direction_status", "direction_reason", "pathways",
                "homolog_status", "homolog_acc", "homolog_pident", "homolog_qcov",
                "homolog_rhea", "homolog_ec", "mapping_suggestion", "review_class", "known_exception", "exception_note", "hmm_evidence",
                "complex_id", "complex_name", "complex_status", "complex_role",
                "complex_evidence", "complex_confidence", "complex_missing",
                "complex_flag", "quaternary_status", "complex_note"]
        lines = ["\t".join(cols)]
        for x in gaps:
            lines.append("\t".join(str(x.get(c, "")) for c in cols))
        (outdir / f"{g}.gaps.tsv").write_text("\n".join(lines) + "\n")
        if res.get("esmfold_closure"):
            cols = ["genome", "gap", "query", "provider", "job_id", "status",
                    "plddt", "ptm", "pdb", "uniprot", "evidence", "limitation", "error"]
            lines = ["\t".join(cols)]
            for row in res["esmfold_closure"]:
                lines.append("\t".join(str(row.get(c, "")) for c in cols))
            (outdir / f"{g}.esmfold_closure.tsv").write_text("\n".join(lines) + "\n")
        if res.get("structural_evidence"):
            cols = ["genome", "gap", "query", "pdb", "uniprot", "plddt", "foldseek_source",
                    "foldseek_top", "foldseek_bits", "foldseek_seq_id", "alphafill_ligands",
                    "active_site_rows", "active_site_conserved", "n_sources", "verdict",
                    "closure", "error", "limitation"]
            lines = ["\t".join(cols)]
            for row in res["structural_evidence"]:
                lines.append("\t".join(str(row.get(c, "")) for c in cols))
            (outdir / f"{g}.structural_evidence.tsv").write_text("\n".join(lines) + "\n")
    if res.get("structural_hypotheses") is not None:
        cols = ["genome", "protein_id", "uniprot", "pdb", "hypothesis", "ec", "cofactor",
                "pocket", "template_mcsa", "foldseek_top", "f1_foldseek", "f2_enzymm",
                "f3_pocket", "f4_alphafill", "f5_active_site", "n_sources", "confidence",
                "review_decision", "limitation", "evidence_json"]
        lines = ["\t".join(cols)]
        for row in res["structural_hypotheses"]:
            lines.append("\t".join(str(row.get(c, "")) for c in cols))
        (outdir / f"{g}.structural_hypotheses.tsv").write_text("\n".join(lines) + "\n")
        (outdir / f"{g}.structural_hypotheses.json").write_text(
            json.dumps(res["structural_hypotheses"], indent=2, ensure_ascii=False))
    if res.get("five_filters"):
        cols = ["pdb", "genome", "uniprot", "f1_foldseek", "foldseek_top", "foldseek_bits",
                "f2_enzymm", "enzymm_hits", "enzymm_templates", "f3_pocket", "p2rank_pockets",
                "pocket_active_hits", "f4_alphafill", "alphafill_ligands", "f5_active_site",
                "active_positions", "n_pass", "priority", "source", "limitation", "error"]
        lines = ["\t".join(cols)]
        for row in res["five_filters"]:
            lines.append("\t".join(str(row.get(c, "")) for c in cols))
        (outdir / f"{g}.structural_filters.tsv").write_text("\n".join(lines) + "\n")
    if res.get("enzymm_scan"):
        cols = ["query_id", "template_pdb_id", "template_mcsa_id", "template_uniprot_id",
                "template_ec", "rmsd", "log_evalue", "matched_residues",
                "template_effective_size", "predicted_correct", "source", "limitation"]
        lines = ["\t".join(cols)]
        for row in res["enzymm_scan"]:
            lines.append("\t".join(str(row.get(c, "")) for c in cols))
        (outdir / f"{g}.enzymm_scan.tsv").write_text("\n".join(lines) + "\n")
    if res.get("esmfold_closure") is not None:
        rows = res.get("esmfold_closure") or []
        summary["esmfold"] = {
            "provider": "local" if any(r.get("provider") == "local" for r in rows) else "neuros nap",
            "n": len(rows),
            "with_pdb": sum(1 for r in rows if r.get("pdb")),
            "error": res.get("esmfold_closure_error", ""),
        }
    if res.get("structural_hypotheses") is not None:
        summary["structural_hypotheses"] = {
            "n": len(res.get("structural_hypotheses") or []),
            "supported": sum(1 for r in (res.get("structural_hypotheses") or []) if r.get("confidence") == "supported_hypothesis"),
        }
    if res.get("five_filters") is not None:
        summary["five_filters"] = {
            "n": len(res.get("five_filters") or []),
            "high": sum(1 for r in (res.get("five_filters") or []) if r.get("priority") == "high"),
            "review": sum(1 for r in (res.get("five_filters") or []) if r.get("priority") == "review"),
        }
    if res.get("enzymm_scan") is not None:
        summary["enzymm"] = {"n_hits": len(res.get("enzymm_scan") or []),
                             "source": "EnzyMM/PyJess + M-CSA"}
    if res.get("structural_evidence") is not None:
        summary["structural_evidence"] = {
            "n": len(res.get("structural_evidence") or []),
            "multi_source": sum(1 for r in (res.get("structural_evidence") or []) if r.get("verdict") == "multi_source_support"),
        }
    (outdir / f"{g}.summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> int:
    ap = argparse.ArgumentParser(description="Anota FASTA desde cero (GFF + red + brechas + informe).")
    ap.add_argument("fasta", nargs="+", type=Path)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--reference", default="", help="GCF de referencia (opcional)")
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--no-cds", action="store_true")
    ap.add_argument("--no-trna", action="store_true")
    ap.add_argument("--no-rrna", action="store_true")
    ap.add_argument("--no-ghost", action="store_true")
    ap.add_argument("--no-deserts", action="store_true")
    ap.add_argument("--operon-max-gap", type=int, default=DEFAULT_MAX_GAP,
                    help="Hueco máximo (pb) para unir CDS en un operón "
                         "putativo por reglas (default 300).")
    ap.add_argument("--esm", action="store_true", help="Puntuar ORFs con ESM-2")
    ap.add_argument("--no-metabolic", action="store_true")
    ap.add_argument("--structure-b-external", action="store_true",
                    help="Permitir P2Rank/AlphaFill/APIs en filtros B (por defecto solo local).")
    ap.add_argument("--no-membrane", action="store_true",
                    help="No ejecutar membrana/transporte (hmmscan + Pfam).")
    ap.add_argument("--no-functional", action="store_true",
                    help="No ejecutar DIAMOND/Swiss-Prot automáticamente.")
    ap.add_argument("--no-oric", action="store_true",
                    help="No ejecutar ORCA/oriC automáticamente.")
    ap.add_argument("--accession", default="",
                    help="Accession (GCF/NC) para consultar DoriC como referencia; opcional.")
    ap.add_argument("--esmfold", choices=["off", "local", "neuros nap"], default="off",
                    help="Cierre de brechas con ESMFold (local o Neurosnap API).")
    ap.add_argument("--esmfold-gaps", type=int, default=3,
                    help="Máximo de brechas a enviar/piegue (default 3).")
    ap.add_argument("--esmfold-out", type=Path, default=None,
                    help="Directorio de PDB/métricas (default OUTDIR/GENOME_esmfold).")
    ap.add_argument("--structure-evidence", action="store_true",
                    help="Cruzar Foldseek PDB100 + AlphaFill + sitios activos.")
    ap.add_argument("--structural-annotation", action="store_true",
                    help="Escanear PDBs con EnzyMM/PyJess (plantillas M-CSA).")
    ap.add_argument("--structural-filters", action="store_true",
                    help="Ejecutar los 5 filtros estructurales (F1-F5).")
    ap.add_argument("--structure-manifest", type=Path, default=None,
                    help="TSV genome/protein_id/uniprot/pdb para filtros e hipótesis.")
    ap.add_argument("--foldseek-local-db", type=Path, default=None,
                    help="Base Foldseek local (ej. data/db/foldseek_swissprot); evita la API.")
    ap.add_argument("--precomputed-hypotheses", type=Path, default=None,
                    help="TSV de hipótesis precalculadas; filtra por genoma.")
    ap.add_argument("--panel-hmm", type=Path, default=None)
    ap.add_argument("--transport-db", type=Path, default=None)
    args = ap.parse_args()

    panel = args.panel_hmm or (args.data_root / "db" / "pfam" / "panel_focus.hmm")
    tdb = args.transport_db or (args.data_root / "db" / "transport.dmnd")
    for fasta in args.fasta:
        data = read_fasta_bytes(fasta)
        options = {
            "cds": not args.no_cds, "trna": not args.no_trna, "rrna": not args.no_rrna,
            "ghost": not args.no_ghost, "deserts": not args.no_deserts,
            "esm": args.esm, "esm_top": 0 if args.esm else 0,
            "min_region_len": 150, "min_aa": 60,
        }
        res = annotate(fasta.name, data, options)
        if res.get("features"):
            try:
                res["regulation"] = run_regulation(
                    res["features"], args.outdir, res["genome"],
                    max_gap=args.operon_max_gap,
                    sequences=res.get("sequences"))
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] regulación automática no ejecutada: {exc}",
                      file=sys.stderr)
        if not args.no_functional and res.get("proteins"):
            try:
                res["functional"] = functional_annotate(
                    res["proteins"], args.data_root, args.outdir,
                    res["genome"], threads=args.threads)
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] anotación funcional no ejecutada: {exc}", file=sys.stderr)
        if res.get("functional"):
            try:
                dom = run_domains(
                    res["functional"].get("rows") or [], res.get("proteins") or [],
                    args.outdir, res["genome"], args.data_root, cpu=args.threads)
                res["domains"] = dom
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] dominios no evaluados: {exc}", file=sys.stderr)
            try:
                res["structural_b"] = run_structural_b(
                    res["functional"].get("rows") or [], args.outdir,
                    res["genome"], args.data_root)
                res["structural_b_filters"] = evaluate_structural_b(
                    res["functional"].get("rows") or [], args.outdir,
                    res["genome"], args.data_root,
                    foldseek_db=args.foldseek_local_db, cpu=args.threads,
                    external=args.structure_b_external)
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] estado estructural B no evaluado: {exc}",
                      file=sys.stderr)
            try:
                res["context_c"] = run_context(
                    res["functional"].get("rows") or [], res.get("proteins") or [],
                    args.outdir, res["genome"], cpu=args.threads)
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] contexto C no evaluado: {exc}", file=sys.stderr)
            try:
                comp = infer_complexes(res["features"],
                                       res["functional"].get("hits") or {},
                                       res["genome"])
                model = run_complex_model(
                    res["features"], res["functional"].get("hits") or {},
                    args.outdir, res["genome"], args.data_root)
                if model.get("ok"):
                    comp["complexes"] += model["complexes"]
                    comp["by_protein"].update(model["by_protein"])
                    comp["n_candidate"] = model.get("n_candidate", 0)
                transfer = load_complex_transfer(
                    args.data_root, res["genome"], res["features"])
                if transfer.get("ok"):
                    comp["complexes"] += transfer["complexes"]
                    comp["by_protein"].update(transfer["by_protein"])
                    comp["n_candidate"] = (comp.get("n_candidate", 0)
                                           + transfer.get("n_candidate", 0))
                enrich_rows(res["functional"].get("rows") or [], comp)
                write_paths(res["functional"].get("rows") or [], args.outdir,
                            res["genome"])
                comp["protein_rows"] = res["functional"].get("rows") or []
                res["complexes"] = comp
                write_complexes(comp, args.outdir)
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] complejos no evaluados: {exc}", file=sys.stderr)
        if not args.no_membrane and res.get("gff"):
            try:
                mem = run_membrane(res, args.outdir, args.data_root,
                                   cpu=args.threads)
                res["membrane"] = mem
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] membrana no ejecutada: {exc}", file=sys.stderr)
        if not args.no_oric and res.get("sequences"):
            try:
                res["oric"] = oric.run(res, args.data_root,
                                       accession=args.accession)
                oric.write(res, res["oric"], args.outdir)
            except Exception as exc:  # noqa: BLE001
                print(f"[aviso] oriC no ejecutado: {exc}", file=sys.stderr)
        res["report_md"] = build_report(res)
        metabolic = None
        if not args.no_metabolic:
            try:
                metabolic = build_network(
                    res["proteins"], args.data_root, args.reference or None,
                    threads=args.threads, panel_hmm=panel if panel.exists() else None,
                    transport_db=tdb if tdb.exists() else None,
                    producer_db=(args.data_root / "db" / "producer_family.dmnd")
                    if (args.data_root / "db" / "producer_family.dmnd").exists() else None,
                )
                if res.get("complexes"):
                    metabolic["gaps"] = enrich_gaps(
                        metabolic.get("gaps", []), res["complexes"], args.data_root)
                res["metabolic"] = metabolic
                res["report_md"] = build_report(res)
            except Exception as exc:
                print(f"[aviso] red metabólica no construida: {exc}", file=sys.stderr)
        if args.esmfold != "off" and metabolic and metabolic.get("gaps"):
            for g in metabolic["gaps"]:
                g.setdefault("genome", res["genome"])
            esm_out = args.esmfold_out or (args.outdir / f"{res['genome']}_esmfold")
            try:
                rows = close_gaps(
                    metabolic["gaps"], res.get("proteins", []), args.esmfold,
                    esm_out, limit=max(0, args.esmfold_gaps),
                    note=f"dna_poc {res['genome']}",
                )
                res["esmfold_closure"] = rows
                if args.structure_evidence and rows:
                    res["structural_evidence"] = evidence_for_rows(rows)
                if args.structural_annotation and rows:
                    pdbs = [r["pdb"] for r in rows if r.get("pdb")]
                    ann = scan_pdbs(pdbs, args.outdir / f"{res['genome']}_enzymm", jobs=args.threads)
                    res["enzymm_scan"] = summarize_enzymm(ann.get("hits") or [])
                    if ann.get("error"):
                        res["enzymm_error"] = ann["error"]
                if args.structural_filters and rows:
                    local_hits = {}
                    if args.foldseek_local_db and args.foldseek_local_db.exists():
                        pdbs = [row["pdb"] for row in rows if row.get("pdb")]
                        batch = batch_local_foldseek(pdbs, args.foldseek_local_db,
                                                     args.outdir / "foldseek_local_esm",
                                                     threads=args.threads)
                        local_hits = batch.get("hits") or {}
                        res["foldseek_local_esm"] = {k: v for k, v in batch.items() if k != "hits"}
                    out = []
                    for row in rows:
                        if not row.get("pdb"):
                            continue
                        stem = Path(row["pdb"]).stem
                        out.append(run_filters(
                            row["pdb"], genome=res["genome"],
                            uniprot=row.get("uniprot", ""),
                            foldseek_hits=(local_hits.get(stem, [])
                                           if args.foldseek_local_db else None),
                            plddt=_num(row.get("plddt")),
                            ptm=_num(row.get("ptm")),
                        ))
                    res["five_filters"] = out
                    res["structural_hypotheses"] = [
                        build_hypothesis(x, {"genome": res["genome"],
                                             "protein_id": x.get("query", x.get("protein_id", "")),
                                             "uniprot": x.get("uniprot", ""),
                                             "pdb": x.get("pdb", "")})
                        for x in out
                    ]
                if rows:
                    res["report_md"] = build_report(res)
            except Exception as exc:
                print(f"[aviso] ESMFold no ejecutado: {exc}", file=sys.stderr)
                res["esmfold_closure_error"] = str(exc)
        if args.structure_manifest and args.structure_manifest.exists():
            mrows = list(csv.DictReader(open(args.structure_manifest), delimiter="\t"))
            local_hits = {}
            if args.foldseek_local_db and args.foldseek_local_db.exists():
                pdbs = [r["pdb"] for r in mrows if r.get("pdb")]
                batch = batch_local_foldseek(pdbs, args.foldseek_local_db,
                                             args.outdir / "foldseek_local_manifest",
                                             threads=args.threads)
                local_hits = batch.get("hits") or {}
                res["foldseek_local"] = {k: v for k, v in batch.items() if k != "hits"}
            frows = []
            for r in mrows:
                if not r.get("pdb"):
                    continue
                stem = Path(r["pdb"]).stem
                frows.append(run_filters(
                    r["pdb"], genome=r.get("genome", ""),
                    uniprot=r.get("uniprot", ""),
                    foldseek_hits=(local_hits.get(stem, [])
                                   if args.foldseek_local_db else None),
                    plddt=_num(r.get("plddt")),
                    ptm=_num(r.get("ptm")),
                ))
            res["structural_hypotheses"] = from_manifest(mrows, frows)
            res["report_md"] = build_report(res)
        if args.precomputed_hypotheses:
            extra = load_precomputed(args.precomputed_hypotheses, genome=res["genome"])
            if extra:
                res["structural_hypotheses"] = (res.get("structural_hypotheses") or []) + extra
                res["report_md"] = build_report(res)
        write_outputs(res, args.outdir, metabolic)
        print(f"OK {res['genome']}: {res['bp']} pb, cobertura {res.get('coverage',0)*100:.1f}%, "
              f"brechas {len(metabolic.get('gaps', [])) if metabolic else 0} -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
