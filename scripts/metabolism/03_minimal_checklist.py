#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_minimal_checklist.py — Evalúa la completitud del metabolismo mínimo.
Usa una lista de ECs esenciales, los mapea a RHEA IDs y verifica 
qué porcentaje está presente en la red de cada genoma.
"""
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("minimal")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
NET_DIR = PROJECT_ROOT / "data" / "metabolism"
RHEA_EC = PROJECT_ROOT / "data" / "metabolism" / "rhea2ec.tsv"
OUT_FILE = PROJECT_ROOT / "data" / "metabolism" / "minimal_coverage.json"

# Lista curada de enzimas esenciales (ECs)
ESSENTIAL_ECS = {
    "Glycolysis": ["2.7.1.1", "5.3.1.1", "6.2.1.1", "2.7.1.11", "4.1.2.13", "5.4.2.11", "2.7.1.40", "2.7.2.3", "5.4.2.12", "4.2.1.11", "1.2.1.12"],
    "TCA_Cycle": ["2.3.3.1", "4.2.1.3", "1.1.1.42", "1.2.4.2", "6.2.1.4", "1.3.5.1", "4.2.1.2", "1.1.1.37"],
    "PPP_Pentose": ["1.1.1.49", "3.1.1.31", "5.3.1.6", "2.2.1.1", "4.1.2.13", "5.1.3.1"],
    "Nucleotide_Biosynthesis": ["2.4.2.14", "2.4.2.10", "6.3.4.4", "6.3.4.13", "6.3.4.5", "2.7.6.1", "3.6.1.31", "3.5.4.4", "6.3.4.4", "6.3.3.5"],
    "AminoAcid_Core": ["4.3.1.17", "1.4.1.1", "1.4.1.3", "2.6.1.1", "2.6.1.2", "2.7.2.4", "2.7.2.3", "2.7.1.39", "2.7.1.40", "2.7.1.4"],
    "Cell_Envelope": ["6.3.2.8", "1.1.1.100", "2.7.1.148", "2.7.1.49", "2.7.8.13", "2.7.8.12"]
}

def load_ec_to_rhea():
    ec_to_rhea = {}
    with open(RHEA_EC, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rhea_id = f"RHEA:{parts[0].strip()}"
                ec_raw = parts[3].strip()
                for ec in ec_raw.split(';'):
                    ec = ec.strip()
                    if ec:
                        ec_to_rhea.setdefault(ec, []).append(rhea_id)
    return ec_to_rhea

def main():
    ec_to_rhea = load_ec_to_rhea()
    
    # Mapear nuestro checklist a RHEA IDs
    checklist = {}
    for pathway, ecs in ESSENTIAL_ECS.items():
        checklist[pathway] = set()
        for ec in ecs:
            for rhea_id in ec_to_rhea.get(ec, []):
                checklist[pathway].add(rhea_id)
                
    log.info(f"Checklist definido: {sum(len(v) for v in checklist.values())} reacciones RHEA únicas.")
    
    results = {}
    net_files = list(NET_DIR.glob("*_network.json"))
    log.info(f"Evaluando {len(net_files)} genomas...")
    
    for nf in net_files:
        prefix = nf.name.replace("_network.json", "")
        data = json.load(open(nf))
        # Nodos presentes en esta bacteria
        present_rhea = {n['id'] for n in data['nodes'] if n.get('type') == 'reaction'}
        
        genome_stats = {}
        total_present = 0
        total_req = 0
        
        for pathway, req_rheas in checklist.items():
            if not req_rheas: continue
            present = len(req_rheas & present_rhea)
            total_present += present
            total_req += len(req_rheas)
            genome_stats[pathway] = round(100 * present / len(req_rheas), 1)
            
        genome_stats["OVERALL_COVERAGE"] = round(100 * total_present / total_req, 1) if total_req > 0 else 0
        results[prefix] = genome_stats

    OUT_FILE.write_text(json.dumps(results, indent=2))
    log.info(f"✅ Evaluación guardada en {OUT_FILE.name}")
    
    # Mostrar Top 3
    top = sorted(results.items(), key=lambda x: x[1]["OVERALL_COVERAGE"], reverse=True)[:3]
    for g, stats in top:
        log.info(f"  {g[:30]}: {stats['OVERALL_COVERAGE']}% completitud (Glicólisis: {stats.get('Glycolysis')}%)")

if __name__ == "__main__":
    main()
