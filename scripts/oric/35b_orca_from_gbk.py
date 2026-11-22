#!/usr/bin/env python3
"""
35_orca_from_gbk.py — Predice oriC con ORCA a partir de GenBank ab initio.

API orcapy 1.0.1:
    model = joblib.load(ORCA_RFC_model_1_4_0.pkl.gz)
    orca  = ORCA.from_gbk(path, model=model)
    orca.find_oriCs(show_info=False, plot_path=None)
    -> orca.predictions  (probabilidades)
    -> orca.oriC_middles (posiciones, centro del oriC)

Multi-record: ORCA solo acepta un registro. Se procesa el cromosoma (mas
grande); los plasmidos se reportan como skipped (ORCA no modela oriC
plasmidico; usar PlasmidFinder).

v2 — Fase B del gate (separacion A/B):

  A (35a_marker_validation.py) valida los marcadores dnaA/dnaN contra la
  anotacion RefSeq ya descargada (tolerante a numeracion/guiones/mayusculas,
  unicidad por gene=, interior de contig, distancia de operon) y emite
  veredictos en data/marker_validation/<prefix>.json con status
  validated / weak / ambiguous / missing / error.

  B (este script) consume esos veredictos (--markers-dir) y SOLO emite
  predicciones confiables para genomas validated (o weak con --accept-weak).
  Sin JSON de 35A cae al gate legado por texto del gbk (substring, menos
  robusto). Tras la anotacion funcional (37) correr con
  --gbk-name annotation_func.gbk; para la ruta RefSeq-first, el gate no
  depende del gbk y basta annotation.gbk ab initio.

Salida: data/orca_predictions/<prefix>.json
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

import joblib
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("orca")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT_ROOT / "data" / "annotations_ab_initio"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "orca_predictions"
DEFAULT_MODEL = PROJECT_ROOT / "data" / "models" / "ORCA_RFC_model_1_4_0.pkl.gz"
DEFAULT_MARKERS = PROJECT_ROOT / "data" / "marker_validation"
SUMMARY_FILE = "orca_summary.json"
FINAL_STATUSES = ("ok", "low_confidence")


@dataclass
class ContigPrediction:
    contig: str
    length: int
    oric_position: int = 0
    confidence: float = 0.0
    status: str = "pending"
    reason: str = ""


@dataclass
class GenomePrediction:
    prefix: str
    timestamp: str
    model: str = ""
    chromosome: dict = field(default_factory=dict)
    plasmids: list = field(default_factory=list)
    status: str = "pending"
    error: str = ""
    gate: str = ""  # fuente del gate: 35A:validated | 35A:weak(aceptado) | gbk-text


def split_gbk(gbk_path: Path, workdir: Path) -> list[tuple[int, str, Path]]:
    """Split multi-record GenBank; returns [(length, id, path), ...]."""
    records = list(SeqIO.parse(str(gbk_path), "genbank"))
    out = []
    for i, rec in enumerate(records):
        safe = "".join(c if c.isalnum() or c in "._-" else "_"
                       for c in rec.id)[:50]
        p = workdir / f"rec_{i:02d}_{safe}.gbk"
        SeqIO.write([rec], str(p), "genbank")
        out.append((len(rec.seq), rec.id, p))
    return out


def check_marker_genes(gbk_path: Path) -> dict:
    """Return {dnaA: bool, dnaN: bool} by scanning gene qualifiers."""
    found = {"dnaA": False, "dnaN": False}
    try:
        for rec in SeqIO.parse(str(gbk_path), "genbank"):
            for feat in rec.features:
                if feat.type not in ("CDS", "gene"):
                    continue
                for key in ("gene", "locus_tag", "product", "Name"):
                    for val in feat.qualifiers.get(key, []):
                        v = val.lower()
                        if "dnaa" in v or "dnaA" in val:
                            found["dnaA"] = True
                        if "dnan" in v or "dnaN" in val:
                            found["dnaN"] = True
    except Exception:
        pass
    return found


def load_marker_verdict(markers_dir: Path, prefix: str) -> dict | None:
    """Veredicto de fase A (35a) para este genoma; None si no existe."""
    p = markers_dir / f"{prefix}.json"
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
        return d if isinstance(d, dict) and "status" in d else None
    except Exception:
        return None


def resolve_gate(chrom_path: Path, verdict: dict | None, accept_weak: bool):
    """Decide si se confia en la prediccion. El veredicto de 35A manda si existe.

    Returns (gate_open, reason, source).
    """
    if verdict is not None:
        st = verdict.get("status")
        if st == "validated":
            return True, "", "35A:validated"
        if st == "weak" and accept_weak:
            return True, "", "35A:weak(aceptado)"
        reasons = "; ".join(verdict.get("reasons", [])) or st
        return False, f"35A:{st} [{reasons[:160]}]", "35A"
    # Fallback legado: substring sobre gene/locus_tag/product/Name (menos robusto)
    markers = check_marker_genes(chrom_path)
    ok = markers["dnaA"] or markers["dnaN"]
    return (ok,
            "" if ok else "dnaA/dnaN not annotated; ORCA prediction unreliable",
            "gbk-text")


def predict_one_record(gbk_path: Path, model) -> tuple[int, float]:
    """Run ORCA on a single-record GenBank. Returns (middle_pos, confidence)."""
    from orcapy import ORCA

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        orca = ORCA.from_gbk(str(gbk_path), model=model)
        orca.find_oriCs(show_info=False, plot_path=None)

    preds = getattr(orca, "predictions", None)
    middles = getattr(orca, "oriC_middles", None)
    if preds is None or middles is None or len(preds) == 0 or len(middles) == 0:
        raise RuntimeError("ORCA no encontro oriC")

    # Soportar list, tuple, numpy array
    best_idx = 0
    best_prob = float(preds[0])
    for i, v in enumerate(preds):
        fv = float(v)
        if fv > best_prob:
            best_prob = fv
            best_idx = i
    pos = int(middles[best_idx])
    return pos, best_prob


def process_genome(
    prefix: str, gbk_path: Path, model, model_name: str,
    min_confidence: float,
    marker_verdict: dict | None = None, accept_weak: bool = False,
) -> GenomePrediction:
    ts = datetime.now().isoformat(timespec="seconds")
    if not gbk_path.exists():
        return GenomePrediction(prefix=prefix, timestamp=ts,
                                model=model_name, status="missing_gbk")

    workdir = Path(tempfile.mkdtemp(prefix=f"orca_{prefix[:20]}_"))
    try:
        records = split_gbk(gbk_path, workdir)
        if not records:
            return GenomePrediction(prefix=prefix, timestamp=ts,
                                    model=model_name, status="empty_gbk")

        records.sort(reverse=True)  # por longitud descendente
        chrom_len, chrom_id, chrom_path = records[0]

        gate_open, gate_reason, gate_src = resolve_gate(
            chrom_path, marker_verdict, accept_weak)
        try:
            pos, conf = predict_one_record(chrom_path, model)
            if not gate_open:
                status = "needs_functional_annotation"
                reason = gate_reason
            elif conf < min_confidence:
                status = "low_confidence"
                reason = ""
            else:
                status = "ok"
                reason = ""
            chrom = ContigPrediction(
                contig=chrom_id, length=chrom_len,
                oric_position=pos, confidence=conf,
                status=status, reason=reason,
            )
        except Exception as e:
            chrom = ContigPrediction(
                contig=chrom_id, length=chrom_len,
                status="failed", reason=str(e)[:300],
            )

        plasmids = [
            asdict(ContigPrediction(
                contig=pid, length=plen, status="skipped",
                reason="ORCA no modela oriC plasmidico (usar PlasmidFinder)",
            ))
            for plen, pid, _ in records[1:]
        ]

        return GenomePrediction(
            prefix=prefix, timestamp=ts, model=model_name,
            chromosome=asdict(chrom), plasmids=plasmids,
            status=chrom.status,
            error=chrom.reason if chrom.status == "failed" else "",
            gate=gate_src,
        )
    finally:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="ORCA oriC prediction from ab initio GenBank",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    p.add_argument("--gbk-name", default="annotation.gbk",
                   help="Nombre del GenBank dentro de cada directorio; usa "
                        "annotation_func.gbk tras la anotacion funcional (37)")
    p.add_argument("--markers-dir", type=Path, default=DEFAULT_MARKERS,
                   help="Veredictos de fase A (35a_marker_validation.py); "
                        "si no hay JSON por genoma, cae al gate legado del gbk")
    p.add_argument("--accept-weak", action="store_true",
                   help="Tratar status=weak de 35A como gate abierto")
    p.add_argument("--min-confidence", type=float, default=0.5,
                   help="Umbral para status=ok; debajo es low_confidence")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--only", default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)

    input_dir = args.input_dir.resolve()
    out_dir = args.output_dir.resolve()
    model_path = args.model.resolve()

    if not input_dir.is_dir():
        log.error("No existe input: %s", input_dir)
        return 2
    if not model_path.exists():
        log.error("No existe modelo ORCA: %s", model_path)
        return 2

    log.info("Cargando modelo ORCA: %s", model_path.name)
    model = joblib.load(str(model_path))
    model_name = model_path.name

    out_dir.mkdir(parents=True, exist_ok=True)

    prefixes = sorted(
        p.name for p in input_dir.iterdir()
        if p.is_dir() and (p / args.gbk_name).exists()
    )
    if args.only:
        prefixes = [p for p in prefixes if args.only in p]
    if args.limit:
        prefixes = prefixes[: args.limit]

    if not prefixes:
        log.error("No hay genomas con %s en %s", args.gbk_name, input_dir)
        return 2

    log.info("=== ORCA oriC prediction ===")
    log.info("GenBank: %s | Genomas: %d | Gate: --markers-dir=%s%s",
             args.gbk_name, len(prefixes), args.markers_dir,
             " (+--accept-weak)" if args.accept_weak else "")

    results: list[GenomePrediction] = []
    counts = {
        "ok": 0, "low_confidence": 0, "needs_functional_annotation": 0,
        "failed": 0, "skipped": 0, "other": 0,
    }

    for idx, prefix in enumerate(prefixes, 1):
        out_path = out_dir / f"{prefix}.json"
        if not args.force and out_path.exists():
            try:
                existing = json.loads(out_path.read_text())
                if existing.get("status") in FINAL_STATUSES:
                    results.append(GenomePrediction(**existing))
                    counts["skipped"] += 1
                    log.info("[%d/%d] %s — SKIP (%s, ya predicho)",
                             idx, len(prefixes), prefix, existing["status"])
                    continue
            except Exception:
                pass

        verdict = load_marker_verdict(args.markers_dir, prefix)
        result = process_genome(
            prefix, input_dir / prefix / args.gbk_name,
            model, model_name, args.min_confidence,
            marker_verdict=verdict, accept_weak=args.accept_weak,
        )
        results.append(result)
        out_path.write_text(json.dumps(asdict(result), indent=2))

        st = result.status
        if st in counts:
            counts[st] += 1
        else:
            counts["other"] += 1

        c = result.chromosome
        if st == "ok":
            log.info("[%d/%d] %s — OK oriC=%d conf=%.3f (%s, %d bp)",
                     idx, len(prefixes), prefix,
                     c.get("oric_position", 0), c.get("confidence", 0.0),
                     c.get("contig", ""), c.get("length", 0))
        elif st == "low_confidence":
            log.warning("[%d/%d] %s — LOW oriC=%d conf=%.3f",
                        idx, len(prefixes), prefix,
                        c.get("oric_position", 0), c.get("confidence", 0.0))
        elif st == "needs_functional_annotation":
            log.warning("[%d/%d] %s — NEEDS_FUNC (%s)",
                        idx, len(prefixes), prefix,
                        (c.get("reason") or "no dnaA/dnaN")[:120])
        else:
            log.warning("[%d/%d] %s — %s: %s",
                        idx, len(prefixes), prefix, st, result.error[:150])

    summary = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "model": model_name,
        "gbk_name": args.gbk_name,
        "markers_dir": str(args.markers_dir),
        "min_confidence": args.min_confidence,
        "counts": {
            "total": len(results),
            **counts,
        },
        "per_genome": [asdict(r) for r in results],
    }
    summary_path = out_dir / SUMMARY_FILE
    summary_path.write_text(json.dumps(summary, indent=2))

    log.info("=== Fin ===")
    log.info("OK: %d | Low: %d | NeedsFunc: %d | Skipped: %d | Failed: %d",
             counts["ok"], counts["low_confidence"],
             counts["needs_functional_annotation"],
             counts["skipped"], counts["failed"])
    log.info("Summary: %s", summary_path)
    reliable = counts["ok"] + counts["low_confidence"]
    return 0 if reliable > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
