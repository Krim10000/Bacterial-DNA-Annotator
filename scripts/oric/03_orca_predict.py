#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_orca_predict.py — Flujo ORIC (Paso 3).
Ejecuta ORCA usando los marcadores validados y compara contra DoriC (si existe).
"""
import argparse
import json
import logging
import sys
import tempfile
import warnings
from datetime import datetime
from pathlib import Path

import joblib
from Bio import SeqIO
from orcapy import ORCA

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("orca")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
AB_DIR = PROJECT_ROOT / "data" / "annotations_ab_initio"
VAL_DIR = PROJECT_ROOT / "data" / "marker_validation"
DORIC_JSON = PROJECT_ROOT / "data" / "doric_ground_truth.json"
OUT_DIR = PROJECT_ROOT / "data" / "orca_predictions"
MODEL_PATH = PROJECT_ROOT / "data" / "models" / "ORCA_RFC_model_1_4_0.pkl.gz"

def load_doric():
    if DORIC_JSON.exists():
        return json.loads(DORIC_JSON.read_text())
    return {}

def run_orca_on_record(seq, replicon_id, markers, model):
    """Inyecta los marcadores validados y corre ORCA."""
    # Construir lista de locaciones para ORCA: [("dnaA", start, end), ...]
    gene_locs = []
    for mk in ["dnaA", "dnaN"]:
        if mk in markers and markers[mk].get("replicon") == replicon_id:
            gene_locs.append((mk, markers[mk]["start"], markers[mk]["end"]))
            
    if not gene_locs:
        raise ValueError(f"No hay marcadores validados para el replicón {replicon_id}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Usamos from_string para evitar problemas con el GenBank multi-record
        orca = ORCA.from_string(sequence=str(seq), gene_locations=gene_locs, accession=replicon_id, model=model)
        orca.find_oriCs(show_info=False, plot_path=None)

    preds = getattr(orca, "predictions", None)
    middles = getattr(orca, "oriC_middles", None)
    
    if not preds or not middles:
        raise RuntimeError("ORCA no devolvió predicciones")

    best_idx = max(range(len(preds)), key=lambda i: float(preds[i]))
    return int(middles[best_idx]), float(preds[best_idx])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not MODEL_PATH.exists():
        log.error(f"Falta el modelo ORCA: {MODEL_PATH}")
        return

    log.info("Cargando modelo ORCA...")
    model = joblib.load(MODEL_PATH)
    doric = load_doric()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Listar genomas que pasaron validación
    val_files = sorted(VAL_DIR.glob("*.json"))
    log.info(f"Procesando {len(val_files)} genomas...")

    n_ok, n_skip, n_fail = 0, 0, 0

    for idx, vf in enumerate(val_files, 1):
        prefix = vf.stem
        out_file = OUT_DIR / f"{prefix}.json"

        if out_file.exists() and not args.force:
            n_skip += 1
            continue

        val_data = json.loads(vf.read_text())
        if val_data.get("status") == "missing":
            log.warning(f"[{idx}] {prefix}: Faltan marcadores, saltando ORCA.")
            n_skip += 1
            continue

        gbk_path = AB_DIR / prefix / "annotation.gbk"
        if not gbk_path.exists():
            n_skip += 1
            continue

        try:
            # Encontrar el cromosoma principal (el registro más largo)
            records = list(SeqIO.parse(str(gbk_path), "genbank"))
            chrom = max(records, key=lambda r: len(r.seq))
            
            # Correr ORCA
            pos, conf = run_orca_on_record(chrom.seq, chrom.id, val_data.get("markers", {}), model)
            
            # Comparar con DoriC (si existe)
            doric_truth = doric.get(prefix) or doric.get(chrom.id)
            delta_bp = None
            if doric_truth and doric_truth.get("orics"):
                # Tomar el primer oriC de DoriC para este replicón
                d_start = doric_truth["orics"][0]["start"]
                # Distancia circular (el genoma bacteriano es un círculo)
                L = len(chrom.seq)
                raw_delta = abs(pos - d_start)
                delta_bp = min(raw_delta, L - raw_delta)

            result = {
                "genome": prefix,
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "replicon": chrom.id,
                "marker_status": val_data.get("status"),
                "orca_position": pos,
                "orca_confidence": round(conf, 4),
                "doric_truth": doric_truth,
                "delta_bp": delta_bp,
                "status": "ok"
            }
            out_file.write_text(json.dumps(result, indent=2))
            
            doric_str = f"(DoriC: {doric_truth['orics'][0]['start']}, Δ={delta_bp} bp)" if delta_bp is not None else "(Sin DoriC)"
            log.info(f"[{idx}/{len(val_files)}] {prefix}: pos={pos} conf={conf:.3f} {doric_str}")
            n_ok += 1

        except Exception as e:
            log.warning(f"[{idx}/{len(val_files)}] {prefix}: Falló ORCA - {str(e)[:100]}")
            n_fail += 1

    log.info(f"=== Fin ORCA ===")
    log.info(f"OK: {n_ok} | Skip: {n_skip} | Fail: {n_fail}")

if __name__ == "__main__":
    main()
