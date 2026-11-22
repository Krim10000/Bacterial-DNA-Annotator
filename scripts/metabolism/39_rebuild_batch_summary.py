#!/usr/bin/env python3
"""
39_rebuild_batch_summary.py — Reconstruye batch_summary.json desde los report.json.

Problema que resuelve: 33_batch_ab_initio.py reescribe batch_summary.json en
cada corrida SOLO con los items procesados en esa corrida. Si el batch se
lanzo dos veces (barrnap y despues Infernal), el summary final refleja solo
la ultima corrida. Los report.json individuales NUNCA se pierden (salvo que
--force regenere una anotacion completa, que es intencional), asi que el
summary se puede reconstruir sin re-annotar nada.

Que reporta:
  - Totales reconstruidos: genomas con report.json / sin report.json
  - Motor rRNA declarado en tool_versions de cada report.json, para saber si
    el estado en disco corresponde a la corrida con barrnap o a la de Infernal
  - Lista de genomas sin report.json o sin annotation.gbk (los huecos reales)

Uso:
    python scripts/metabolism/39_rebuild_batch_summary.py             # reconstruye y guarda
    python scripts/metabolism/39_rebuild_batch_summary.py --dry-run   # solo reporta, no escribe
    python scripts/metabolism/39_rebuild_batch_summary.py --root data/annotations_ab_initio
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("rebuild")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ROOT = PROJECT_ROOT / "data" / "annotations_ab_initio"
SUMMARY_FILE = "batch_summary.json"


@dataclass
class RebuiltItem:
    prefix: str
    input: str
    status: str  # success | missing_report | missing_gbk
    cds: int = 0
    rrna: int = 0
    trna: int = 0
    elapsed_s: float = 0.0
    genome_bp: int = 0
    contigs: int = 0
    rna_engine: str = ""
    annotated_at: str = ""


def engine_from_versions(versions) -> str:
    """Deducir el motor rRNA desde tool_versions del report.json.

    Nota: check_dependencies registra las herramientas DISPONIBLES. Con
    --rna-tool auto (default del batch), si cmsearch+indice existian se uso
    Infernal; el conteo por motor es una aproximacion honesta.
    """
    if not isinstance(versions, dict):
        return ""
    has_cm = any(k.startswith("cmsearch") for k in versions)
    has_bn = any(k.startswith("barrnap") for k in versions)
    if has_cm and has_bn:
        return "infernal (auto: cmsearch disponible)"
    if has_cm:
        return "infernal"
    if has_bn:
        return "barrnap"
    return "none"


def load_report(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def scan(root: Path) -> tuple[list[RebuiltItem], list[str]]:
    items: list[RebuiltItem] = []
    broken: list[str] = []

    if not root.is_dir():
        log.error("No existe: %s", root)
        return items, broken

    for sub in sorted(p for p in root.iterdir() if p.is_dir()):
        prefix = sub.name
        rep = load_report(sub / "report.json")
        if rep is None:
            broken.append(prefix)
            items.append(RebuiltItem(prefix=prefix, input="", status="missing_report"))
            continue

        gbk_ok = (sub / "annotation.gbk").exists()
        items.append(RebuiltItem(
            prefix=prefix,
            input=str(rep.get("input_path", "")),
            status="success" if gbk_ok else "missing_gbk",
            cds=int(rep.get("cds_features", 0)),
            rrna=int(rep.get("rrna_features", 0)),
            trna=int(rep.get("trna_features", 0)),
            elapsed_s=float(rep.get("elapsed_seconds", 0.0)),
            genome_bp=int(rep.get("genome_bp", 0)),
            contigs=int(rep.get("contigs", 0)),
            rna_engine=engine_from_versions(rep.get("tool_versions", {})),
            annotated_at=str(rep.get("timestamp", "")),
        ))
    return items, broken


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Reconstruye batch_summary.json desde los report.json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--dry-run", action="store_true",
                   help="Solo reporta; no escribe ni hace backup")
    args = p.parse_args(argv)

    root = args.root.resolve()
    items, broken = scan(root)
    if not items:
        log.error("No hay subdirectorios con report.json en %s", root)
        return 2

    ok = [i for i in items if i.status == "success"]
    n_gbk_missing = sum(1 for i in items if i.status == "missing_gbk")
    n_missing = sum(1 for i in items if i.status == "missing_report")

    engines: dict[str, int] = {}
    for i in ok:
        engines[i.rna_engine] = engines.get(i.rna_engine, 0) + 1

    log.info("=== Auditoria %s ===", root)
    log.info("Genomas con report.json+gbk : %d", len(ok))
    log.info("Con report.json pero SIN gbk: %d", n_gbk_missing)
    log.info("SIN report.json             : %d%s", n_missing,
             (" (" + ", ".join(broken[:10]) + ")") if broken else "")
    log.info("Motores rRNA declarados:")
    for eng, n in sorted(engines.items()):
        log.info("  %-40s %d", eng, n)

    total_cds = sum(i.cds for i in ok)
    total_rrna = sum(i.rrna for i in ok)
    total_trna = sum(i.trna for i in ok)
    total_bp = sum(i.genome_bp for i in ok)
    log.info("Totales: CDS=%d rRNA=%d tRNA=%d bp=%s",
             total_cds, total_rrna, total_trna, f"{total_bp:,}")

    if ok:
        ts = [i.annotated_at for i in ok if i.annotated_at]
        if ts:
            log.info("Ventana de timestamps: %s .. %s", min(ts), max(ts))

    if args.dry_run:
        log.info("Dry-run: no se escribe nada.")
        return 0

    summary_path = root / SUMMARY_FILE
    if summary_path.exists():
        bak = root / (SUMMARY_FILE + ".bak_" +
                      datetime.now().strftime("%Y%m%d_%H%M%S"))
        shutil.move(str(summary_path), str(bak))
        log.warning("Summary anterior respaldado en %s", bak.name)

    doc = {
        "rebuilt": True,
        "rebuilt_at": datetime.now().isoformat(timespec="seconds"),
        "counts": {
            "total_seen": len(items),
            "success": len(ok),
            "missing_gbk": n_gbk_missing,
            "missing_report": n_missing,
        },
        "totals": {
            "cds": total_cds,
            "rrna": total_rrna,
            "trna": total_trna,
            "bp": total_bp,
        },
        "rna_engines": engines,
        "missing_report": broken,
        "items": [asdict(i) for i in items],
    }
    summary_path.write_text(json.dumps(doc, indent=2))
    log.info("Escrito: %s", summary_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
