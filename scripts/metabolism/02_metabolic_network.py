#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construye la red metabólica de un genoma a partir de sus proteínas."""
import argparse
import json
import logging
from pathlib import Path
import networkx as nx

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("metab_net")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
AB_DIR = PROJECT_ROOT / "data" / "annotations_ab_initio"
RHEA_UNIPROT = PROJECT_ROOT / "data" / "metabolism" / "rhea2uniprot_sprot.tsv"
OUT_DIR = PROJECT_ROOT / "data" / "metabolism"

def load_uniprot_to_rhea():
    """Carga el mapeo basado en la estructura real del TSV (4 columnas)."""
    uniprot_to_rhea = {}
    if not RHEA_UNIPROT.exists():
        log.error(f"No se encontró {RHEA_UNIPROT}.")
        return uniprot_to_rhea
    
    with open(RHEA_UNIPROT, 'r') as f:
        # Saltar la cabecera
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rhea_id = f"RHEA:{parts[0].strip()}"
                uniprot_id = parts[3].strip()
                
                if uniprot_id not in uniprot_to_rhea:
                    uniprot_to_rhea[uniprot_id] = []
                if rhea_id not in uniprot_to_rhea[uniprot_id]:
                    uniprot_to_rhea[uniprot_id].append(rhea_id)
                
    log.info(f"Cargadas {len(uniprot_to_rhea)} proteínas UniProt mapeadas a Rhea.")
    return uniprot_to_rhea

def extract_uniprot_ids(tsv_path):
    """Extrae IDs de UniProt de la columna 2 del output de DIAMOND (sp|Q9PLU4|...)."""
    uniprots = set()
    if not tsv_path.exists():
        return uniprots
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sseqid = parts[1]
                if "|" in sseqid:
                    up_id = sseqid.split("|")[1]
                    uniprots.add(up_id)
    return uniprots

def build_network(uniprots, up_to_rhea):
    G = nx.Graph()
    for up in uniprots:
        G.add_node(up, type="protein")
        if up in up_to_rhea:
            for rhea_id in up_to_rhea[up]:
                G.add_node(rhea_id, type="reaction")
                G.add_edge(up, rhea_id, relation="catalyzes")
    return G

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    up_to_rhea = load_uniprot_to_rhea()
    
    prefixes = sorted([p.name for p in AB_DIR.iterdir() if p.is_dir()])
    if args.limit:
        prefixes = prefixes[:args.limit]

    log.info(f"Construyendo redes para {len(prefixes)} genomas...")
    
    summary = {}
    for i, prefix in enumerate(prefixes, 1):
        tsv = AB_DIR / prefix / "functional.tsv"
        uniprots = extract_uniprot_ids(tsv)
        
        if not uniprots:
            continue
            
        G = build_network(uniprots, up_to_rhea)
        
        n_prot = len([n for n, d in G.nodes(data=True) if d['type'] == 'protein'])
        n_rx = len([n for n, d in G.nodes(data=True) if d['type'] == 'reaction'])
        n_edges = G.number_of_edges()
        
        gaps = [n for n in G.nodes() if G.degree(n) == 0 and G.nodes[n]['type'] == 'protein']
        
        summary[prefix] = {
            "prot_total": n_prot,
            "rhea_total": n_rx,
            "edges": n_edges,
            "gaps_count": len(gaps),
            "gaps_ids": gaps
        }
        
        # Evitar el FutureWarning de NetworkX usando edges="links"
        net_data = nx.node_link_data(G, edges="links")
        (OUT_DIR / f"{prefix}_network.json").write_text(json.dumps(net_data))
        
        if i % 10 == 0 or i == len(prefixes):
            log.info(f"[{i}/{len(prefixes)}] {prefix}: {n_prot} proteínas, {n_rx} reacciones, {len(gaps)} gaps.")

    (OUT_DIR / "metabolism_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("✅ Redes construidas y guardadas en data/metabolism/")

if __name__ == "__main__":
    main()
