"""ORCA + evidencia local para predecir oriC en un genoma nuevo.

Usa ORCA sobre el cromosoma (contig más grande) y combina GC skew, cajas DnaA
e iterones. DoriC no se consulta para un FASTA arbitrario: es una verdad de
referencia para genomas con accession conocido.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

try:  # ejecutado como paquete (webapp.oric)
    from .origin_pipeline import (classify, dnaa_boxes, gc_skew_origin, iteron_scan)
except ImportError:  # ejecutado con webapp/ en sys.path
    from origin_pipeline import (classify, dnaa_boxes, gc_skew_origin, iteron_scan)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "data" / "models" / "ORCA_RFC_model_1_4_0.pkl.gz"
# Evita que orcapy/matplotlib intente escribir en ~/.config en entornos de solo lectura.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dna_mpl"))
_MODELS: dict[str, Any] = {}


def _load_model(path: Path):
    key = str(path)
    if key not in _MODELS:
        import joblib
        _MODELS[key] = joblib.load(str(path))
    return _MODELS[key]


def _largest(sequences: dict[str, str]) -> tuple[str, str]:
    name = max(sequences, key=lambda k: len(sequences[k]))
    return name, sequences[name]


def _write_gbk(sequences: dict[str, str], features: list[dict],
               hits: dict[str, dict], out: Path) -> str:
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    name, seq = _largest(sequences)
    # ORCA/BioFile espera "accession.version"; un FASTA subido no lo trae.
    acc = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)[:40]
    if "." not in acc:
        acc += ".1"
    try:  # Biopython < 1.78
        from Bio.Alphabet import generic_dna
        seq_obj = Seq(seq, generic_dna)
    except ImportError:  # Biopython >= 1.78
        seq_obj = Seq(seq)
    rec = SeqRecord(seq_obj, id=acc, name=acc[:20],
                    description="Bacterial DNA Annotator ab initio")
    rec.annotations["molecule_type"] = "DNA"
    for f in features:
        if f.get("contig") != name or f.get("type") not in ("CDS", "tRNA", "rRNA"):
            continue
        h = hits.get(f.get("id", ""), {})
        quals: dict[str, list[str]] = {"locus_tag": [f.get("id", "cds")]}
        if h.get("gene"):
            quals["gene"] = [str(h["gene"])]
        if h.get("product"):
            quals["product"] = [str(h["product"])]
        if h.get("ec"):
            quals["EC_number"] = [str(h["ec"])]
        strand = -1 if f.get("strand") == "-" else 1
        rec.features.append(SeqFeature(
            FeatureLocation(int(f["start"]) - 1, int(f["end"]), strand),
            type="CDS", qualifiers=quals))
    out.parent.mkdir(parents=True, exist_ok=True)
    from Bio import SeqIO
    SeqIO.write([rec], str(out), "genbank")
    return name


def _run_orca(gbk: Path, model_path: Path) -> tuple[int | None, float | None, str]:
    try:
        import joblib
        from orcapy import ORCA
    except Exception as exc:  # noqa: BLE001
        return None, None, f"orcapy no disponible: {exc}"
    if not model_path.exists():
        return None, None, f"falta el modelo ORCA: {model_path}"
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = _load_model(model_path)
            orca = ORCA.from_gbk(str(gbk), model=model)
            orca.find_oriCs(show_info=False, plot_path=None)
        preds, middles = orca.predictions, orca.oriC_middles
        if len(preds) == 0 or len(middles) == 0:
            return None, None, "ORCA no encontró oriC"
        best = max(range(len(preds)), key=lambda i: float(preds[i]))
        return int(middles[best]), float(preds[best]), ""
    except Exception as exc:  # noqa: BLE001
        return None, None, f"ORCA falló: {type(exc).__name__}: {exc}"


def doric_position(data_root: Path | None, accession: str) -> int | None:
    """Busca un oriC curado de DoriC por accession (solo si el usuario lo da)."""
    if not data_root or not accession:
        return None
    import csv
    d = Path(data_root) / "doric"
    if not d.exists():
        return None
    for path in sorted(d.glob("*.csv")):
        with path.open(encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                if accession not in (row.get("refseq", ""), row.get("nc", "")):
                    continue
                try:
                    orics = json.loads(row.get("orics", "") or "[]")
                    o = orics[0]
                    return (int(o["oric_start"]) + int(o["oric_end"])) // 2
                except (ValueError, KeyError, IndexError, TypeError):
                    return None
    return None


def run(result: dict[str, Any], data_root: Path | None = None,
        model_path: Path | None = None, accession: str = "") -> dict[str, Any]:
    """Devuelve estado, posición y evidencia para un resultado de annotate()."""
    sequences = result.get("sequences") or {}
    if not sequences:
        return {"status": "unavailable", "reason": "sin secuencias de contigs"}
    chrom, seq = _largest(sequences)
    features = result.get("features") or []
    hits = ((result.get("functional") or {}).get("hits") or {})
    size = len(seq)
    skew_pos, _ = gc_skew_origin(seq)
    # dnaA/dnaN por nombre funcional
    dna_a = None
    for f in features:
        h = hits.get(f.get("id", ""), {})
        if str(h.get("gene", "")).lower() == "dnaa":
            dna_a = int(f.get("start", 0))
            break
    doric = doric_position(data_root, accession)
    model = model_path or (Path(data_root) / "models" / "ORCA_RFC_model_1_4_0.pkl.gz"
                           if data_root else DEFAULT_MODEL)
    with tempfile.TemporaryDirectory(prefix="oric_") as td:
        gbk = Path(td) / "chromosome.gbk"
        _write_gbk(sequences, features, hits, gbk)
        pos, conf, err = _run_orca(gbk, Path(model))

    boxes = len(dnaa_boxes(seq, max(0, (pos or skew_pos) - 10000),
                           min(size, (pos or skew_pos) + 10000)))
    iters = len(iteron_scan(seq, pos or skew_pos))
    if pos is None and doric is not None:
        return {
            "genome": result.get("genome", ""), "status": "validated",
            "position": doric, "confidence": 1.0,
            "validated_by": "DoriC", "reason": f"DoriC para {accession}",
            "attempted": "DoriC; ORCA; GC skew",
            "chromosome": chrom, "size": size, "gc_skew": skew_pos,
            "dnaA": dna_a, "boxes": boxes, "iterons": iters,
        }
    if pos is None:
        return {
            "status": "unavailable", "reason": err, "position": None,
            "chromosome": chrom, "size": size, "gc_skew": skew_pos,
            "dnaA": dna_a, "boxes": boxes, "iterons": iters,
            "attempted": "ORCA; GC skew; cajas DnaA; iterones",
        }
    out = classify(
        genome=result.get("genome", ""), phylum="", organism="", size=size,
        marker_status="present" if dna_a is not None else "missing",
        dnaA=dna_a, orca_pos=pos, orca_conf=conf, doric=doric,
        skew=skew_pos, boxes=boxes, L=size, iterons=iters,
    )
    out.update({
        "chromosome": chrom, "size": size, "gc_skew": skew_pos,
        "dnaA": dna_a, "boxes": boxes, "iterons": iters,
        "model": Path(model).name, "doric": doric,
        "attempted": "ORCA; GC skew; cajas DnaA; iterones; DoriC" if doric else
                     "ORCA; GC skew; cajas DnaA; iterones",
    })
    return out


def gff_text(oric: dict[str, Any]) -> str:
    pos = oric.get("position")
    if not pos:
        return "##gff-version 3\n"
    return (
        f"##gff-version 3\n{oric.get('chromosome','')}\tORCA\toriC\t{pos}\t{pos}\t.\t+\t.\t"
        f"status={oric.get('status','')};confidence={oric.get('confidence','')};"
        f"validated_by={oric.get('validated_by','')}\n"
    )


def write(result: dict[str, Any], oric: dict[str, Any], outdir: Path) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    g = result.get("genome", "genome")
    json_path = outdir / f"{g}.oric.json"
    json_path.write_text(json.dumps(oric, ensure_ascii=False, indent=2),
                         encoding="utf-8")
    gff = outdir / f"{g}.oric.gff3"
    pos = oric.get("position")
    if pos:
        gff.write_text(gff_text(oric), encoding="utf-8")
    return {"json": json_path, "gff": gff}
