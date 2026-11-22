#!/usr/bin/env python3
"""
01_compare_ncbi.py — Compare ab initio GFFs against NCBI ground truth.

For each genome:
  - Read NCBI GFF (data/annotations/<prefix>.gff.gz)
  - Read ab initio GFF (data/annotations_ab_initio/<prefix>/annotation.gff3)
  - Match features by type + strand + reciprocal overlap >= threshold
  - Compute TP/FP/FN, recall/precision/F1 per category (CDS, tRNA, rRNA)

Excluye CDS con pseudo=true del ground truth NCBI (son pseudogenes).
Mantiene partial=true (genes truncados en el contig siguen siendo reales).

Output:
  data/benchmark/<prefix>_comparison.json  (por genoma)
  data/benchmark/benchmark_summary.json    (agregado)

Uso:
    python scripts/validation/01_compare_ncbi.py
    python scripts/validation/01_compare_ncbi.py --limit 5
    python scripts/validation/01_compare_ncbi.py --min-overlap 0.6
    python scripts/validation/01_compare_ncbi.py --only GCF_000006685
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("benchmark")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_AB_INITIO = PROJECT_ROOT / "data" / "annotations_ab_initio"
DEFAULT_NCBI = PROJECT_ROOT / "data" / "annotations"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "benchmark"
SUMMARY_FILE = "benchmark_summary.json"

CATEGORIES = ("CDS", "tRNA", "rRNA")


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def _parse_attrs(attrs: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for pair in attrs.split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            out[k] = v
    return out


def parse_gff(
    path: Path,
    categories: tuple[str, ...],
    exclude_pseudo: bool = False,
) -> dict[str, dict[str, dict[str, list[tuple[int, int]]]]]:
    """
    Returns: { category: { contig: { strand: [(start, end), ...] } } }
    GFF coords are 1-based inclusive; kept as-is.
    """
    result: dict = {cat: {} for cat in categories}
    opener = gzip.open if str(path).endswith(".gz") else open

    n_features = 0
    n_excluded_pseudo = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in categories:
                continue
            if exclude_pseudo:
                attrs = _parse_attrs(cols[8])
                if attrs.get("pseudo") == "true":
                    n_excluded_pseudo += 1
                    continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            strand = cols[6]
            if strand not in ("+", "-"):
                strand = "."
            contig = cols[0]
            result[ftype].setdefault(contig, {}).setdefault(strand, []).append(
                (start, end)
            )
            n_features += 1

    if exclude_pseudo and n_excluded_pseudo:
        log.debug("  Excluidos %d CDS con pseudo=true", n_excluded_pseudo)

    return result


# ---------------------------------------------------------------------------
# Feature matching
# ---------------------------------------------------------------------------

def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _reciprocal(a: tuple[int, int], b: tuple[int, int]) -> float:
    ov = _overlap(a[0], a[1], b[0], b[1])
    if ov <= 0:
        return 0.0
    len_a = a[1] - a[0] + 1
    len_b = b[1] - b[0] + 1
    return min(ov / len_a, ov / len_b)


def match_features(
    ref: dict[str, dict[str, list[tuple[int, int]]]],
    pred: dict[str, dict[str, list[tuple[int, int]]]],
    min_overlap: float,
) -> tuple[int, int, int]:
    """
    Greedy matching by descending reciprocal overlap, scoped per (contig, strand).
    Indices are local to each (contig, strand); matched sets are also local.
    """
    common_contigs = set(ref) & set(pred)
    ref_total = sum(len(f) for c in ref for f in ref[c].values())
    pred_total = sum(len(f) for c in pred for f in pred[c].values())

    tp = 0
    for contig in common_contigs:
        for strand in set(ref[contig]) & set(pred[contig]):
            ref_feats = ref[contig][strand]
            pred_feats = pred[contig][strand]
            local_pairs: list[tuple[float, int, int]] = []
            for i, rf in enumerate(ref_feats):
                for j, pf in enumerate(pred_feats):
                    recip = _reciprocal(rf, pf)
                    if recip >= min_overlap:
                        local_pairs.append((recip, i, j))
            local_pairs.sort(reverse=True)
            matched_ref_local: set[int] = set()
            matched_pred_local: set[int] = set()
            for _, ri, pj in local_pairs:
                if ri in matched_ref_local or pj in matched_pred_local:
                    continue
                matched_ref_local.add(ri)
                matched_pred_local.add(pj)
                tp += 1

    fp = pred_total - tp
    fn = ref_total - tp
    return tp, fp, fn

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _metrics(tp: int, fp: int, fn: int) -> dict:
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return {
        "tp": tp, "fp": fp, "fn": fn,
        "recall": round(recall, 4),
        "precision": round(precision, 4),
        "f1": round(f1, 4),
    }


# ---------------------------------------------------------------------------
# Per-genome benchmark
# ---------------------------------------------------------------------------

@dataclass
class CategoryResult:
    tp: int = 0
    fp: int = 0
    fn: int = 0
    recall: float = 0.0
    precision: float = 0.0
    f1: float = 0.0
    ref_total: int = 0
    pred_total: int = 0


@dataclass
class GenomeResult:
    prefix: str
    status: str  # ok | missing_ncbi | missing_ab | no_common_contigs | parse_error
    contigs_ref: int = 0
    contigs_pred: int = 0
    contigs_common: int = 0
    categories: dict = field(default_factory=dict)
    error: str = ""


def benchmark_one(
    prefix: str,
    ab_dir: Path,
    ncbi_dir: Path,
    min_overlap: float,
) -> GenomeResult:
    ab_gff = ab_dir / prefix / "annotation.gff3"
    ncbi_gff = ncbi_dir / f"{prefix}.gff.gz"

    if not ab_gff.exists():
        return GenomeResult(prefix=prefix, status="missing_ab")
    if not ncbi_gff.exists():
        return GenomeResult(prefix=prefix, status="missing_ncbi")

    try:
        ref = parse_gff(ncbi_gff, CATEGORIES, exclude_pseudo=True)
        pred = parse_gff(ab_gff, CATEGORIES, exclude_pseudo=False)
    except Exception as e:
        return GenomeResult(
            prefix=prefix, status="parse_error", error=str(e)[:200]
        )

    contigs_ref = set().union(*[set(ref[c]) for c in CATEGORIES]) if ref else set()
    contigs_pred = set().union(*[set(pred[c]) for c in CATEGORIES]) if pred else set()
    common = contigs_ref & contigs_pred

    if not common:
        return GenomeResult(
            prefix=prefix, status="no_common_contigs",
            contigs_ref=len(contigs_ref),
            contigs_pred=len(contigs_pred),
        )

    categories_out: dict[str, dict] = {}
    for cat in CATEGORIES:
        tp, fp, fn = match_features(ref[cat], pred[cat], min_overlap)
        metrics = _metrics(tp, fp, fn)
        ref_total = sum(len(f) for c in ref[cat] for f in ref[cat][c].values())
        pred_total = sum(len(f) for c in pred[cat] for f in pred[cat][c].values())
        metrics["ref_total"] = ref_total
        metrics["pred_total"] = pred_total
        categories_out[cat] = metrics

    return GenomeResult(
        prefix=prefix,
        status="ok",
        contigs_ref=len(contigs_ref),
        contigs_pred=len(contigs_pred),
        contigs_common=len(common),
        categories=categories_out,
    )


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------

def aggregate(results: list[GenomeResult], min_overlap: float) -> dict:
    ok = [r for r in results if r.status == "ok"]

    counts = {
        "total": len(results),
        "ok": len(ok),
        "missing_ncbi": sum(1 for r in results if r.status == "missing_ncbi"),
        "missing_ab": sum(1 for r in results if r.status == "missing_ab"),
        "no_common_contigs": sum(1 for r in results if r.status == "no_common_contigs"),
        "parse_error": sum(1 for r in results if r.status == "parse_error"),
    }

    micro: dict = {}
    macro: dict = {}
    for cat in CATEGORIES:
        tp = sum(r.categories[cat]["tp"] for r in ok)
        fp = sum(r.categories[cat]["fp"] for r in ok)
        fn = sum(r.categories[cat]["fn"] for r in ok)
        ref_total = sum(r.categories[cat]["ref_total"] for r in ok)
        pred_total = sum(r.categories[cat]["pred_total"] for r in ok)
        metrics = _metrics(tp, fp, fn)
        metrics["ref_total"] = ref_total
        metrics["pred_total"] = pred_total
        micro[cat] = metrics

        if ok:
            macro[cat] = {
                "avg_recall": round(sum(r.categories[cat]["recall"] for r in ok) / len(ok), 4),
                "avg_precision": round(sum(r.categories[cat]["precision"] for r in ok) / len(ok), 4),
                "avg_f1": round(sum(r.categories[cat]["f1"] for r in ok) / len(ok), 4),
            }

    return {
        "started": datetime.now().isoformat(timespec="seconds"),
        "min_reciprocal_overlap": min_overlap,
        "exclude_pseudo": True,
        "counts": counts,
        "micro_average": micro,
        "macro_average": macro,
        "per_genome": [asdict(r) for r in results],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Benchmark ab initio GFFs against NCBI ground truth",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--ab-dir", type=Path, default=DEFAULT_AB_INITIO)
    p.add_argument("--ncbi-dir", type=Path, default=DEFAULT_NCBI)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--min-overlap", type=float, default=0.5,
                   help="Reciprocal overlap threshold [0-1]")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--only", default=None,
                   help="Solo genomas cuyo prefix contenga este substring")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)

    ab_dir = args.ab_dir.resolve()
    ncbi_dir = args.ncbi_dir.resolve()
    out_dir = args.output_dir.resolve()

    if not ab_dir.is_dir():
        log.error("No existe: %s", ab_dir)
        return 2
    if not ncbi_dir.is_dir():
        log.error("No existe: %s", ncbi_dir)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)

    prefixes = sorted(
        p.name for p in ab_dir.iterdir()
        if p.is_dir() and (p / "annotation.gff3").exists()
    )
    if args.only:
        prefixes = [p for p in prefixes if args.only in p]
    if args.limit:
        prefixes = prefixes[: args.limit]

    if not prefixes:
        log.error("No hay genomas con annotation.gff3 en %s", ab_dir)
        return 2

    log.info("=== Benchmark ab initio vs NCBI ===")
    log.info("AB dir   : %s", ab_dir)
    log.info("NCBI dir : %s", ncbi_dir)
    log.info("Output   : %s", out_dir)
    log.info("Genomas  : %d", len(prefixes))
    log.info("Min overlap: %.2f", args.min_overlap)

    results: list[GenomeResult] = []
    t0 = time.time()

    for idx, prefix in enumerate(prefixes, 1):
        result = benchmark_one(prefix, ab_dir, ncbi_dir, args.min_overlap)
        results.append(result)

        if result.status != "ok":
            log.warning("[%d/%d] %s — %s", idx, len(prefixes), prefix, result.status)
            continue

        cds = result.categories["CDS"]
        trna = result.categories["tRNA"]
        rrna = result.categories["rRNA"]
        log.info(
            "[%d/%d] %s — CDS F1=%.3f (R=%.3f P=%.3f, TP=%d FP=%d FN=%d) | "
            "tRNA F1=%.2f | rRNA F1=%.2f",
            idx, len(prefixes), prefix,
            cds["f1"], cds["recall"], cds["precision"],
            cds["tp"], cds["fp"], cds["fn"],
            trna["f1"], rrna["f1"],
        )

        if idx % 25 == 0:
            partial = aggregate(results, args.min_overlap)
            (out_dir / "benchmark_partial.json").write_text(
                json.dumps(partial, indent=2)
            )

        per_genome_path = out_dir / f"{prefix}_comparison.json"
        per_genome_path.write_text(json.dumps(asdict(result), indent=2))

    summary = aggregate(results, args.min_overlap)
    summary_path = out_dir / SUMMARY_FILE
    summary_path.write_text(json.dumps(summary, indent=2))

    elapsed = time.time() - t0
    log.info("=== Fin en %.1fs ===", elapsed)
    log.info("Counts: %s", summary["counts"])

    for cat in CATEGORIES:
        micro = summary["micro_average"][cat]
        log.info(
            "  %-4s micro: F1=%.4f R=%.4f P=%.4f (TP=%d FP=%d FN=%d)",
            cat, micro["f1"], micro["recall"], micro["precision"],
            micro["tp"], micro["fp"], micro["fn"],
        )

    log.info("Summary: %s", summary_path)
    return 0 if summary["counts"]["ok"] > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
