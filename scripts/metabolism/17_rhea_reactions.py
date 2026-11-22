#!/usr/bin/env python3
"""Construye una red metabólica masiva usando TODAS las enzimas de Rhea."""
import networkx as nx, json
from pathlib import Path

def load_mappings(base_dir):
    tsv_dir = base_dir / "tsv"
    ec_to_master = {}
    master_to_lr = {}
    rhea_to_smiles = {}
    smiles_to_chebi = {}
    chebi_to_name = {}

    with open(tsv_dir / "rhea2ec.tsv", "r", errors="ignore") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                master_id = parts[2]
                ec = parts[3]
                ec_to_master.setdefault(ec, set()).add(master_id)

    with open(tsv_dir / "rhea-directions.tsv", "r", errors="ignore") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                master_to_lr[parts[0]] = parts[1]

    with open(tsv_dir / "rhea-reaction-smiles.tsv", "r", errors="ignore") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rhea_to_smiles[parts[0]] = parts[1]

    with open(tsv_dir / "rhea-chebi-smiles.tsv", "r", errors="ignore") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                smiles_to_chebi[parts[1]] = parts[0]

    with open(tsv_dir / "chebiId_name.tsv", "r", errors="ignore") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chebi_to_name[parts[0]] = parts[1].strip()

    return ec_to_master, master_to_lr, rhea_to_smiles, smiles_to_chebi, chebi_to_name

def get_human_name(smiles, smiles_to_chebi, chebi_to_name):
    chebi_id = smiles_to_chebi.get(smiles)
    if chebi_id:
        return chebi_to_name.get(chebi_id, chebi_id)
    return smiles

def main():
    root = Path.cwd()
    data_dir = root / "data" / "reactions"
    
    print("Cargando archivos de Rhea...")
    ec_to_master, master_to_lr, rhea_to_smiles, smiles_to_chebi, chebi_to_name = load_mappings(data_dir)
    print(f"Enzimas EC disponibles en Rhea: {len(ec_to_master)}")

    G = nx.DiGraph()
    
    # ESCALAMIENTO: Usamos todas las enzimas disponibles (hasta un límite seguro para PyVis)
    all_ecs = list(ec_to_master.keys())[:300] 
    
    print(f"\nConstruyendo red masiva con {len(all_ecs)} enzimas...")
    found = 0
    for ec in all_ecs:
        masters = ec_to_master.get(ec, [])
        for master in masters:
            lr_id = master_to_lr.get(master)
            if not lr_id: continue
            
            smiles_rxn = rhea_to_smiles.get(lr_id)
            if not smiles_rxn or ">>" not in smiles_rxn: continue
            
            left, right = smiles_rxn.split(">>")
            substrates = [get_human_name(smi, smiles_to_chebi, chebi_to_name) for smi in left.split(".")]
            products = [get_human_name(smi, smiles_to_chebi, chebi_to_name) for smi in right.split(".")]
            
            found += 1
            if found <= 10: # Solo imprimir las primeras 10 para no saturar la consola
                print(f"  EC {ec} -> {' + '.join(substrates)} => {' + '.join(products)}")
            
            for sub in substrates: G.add_node(sub, type="compound")
            for prod in products: G.add_node(prod, type="compound")
                
            for sub in substrates:
                for prod in products:
                    G.add_edge(sub, prod, enzyme=ec, rhea_id=lr_id)
                    
    out_file = data_dir / "reaction_network.json"
    data = nx.node_link_data(G, edges="links")
    with open(out_file, "w") as f:
        json.dump(data, f)
        
    print(f"\n✅ Se procesaron {found} reacciones reales.")
    print(f"Red masiva construida con {G.number_of_nodes()} metabolitos y {G.number_of_edges()} conexiones.")
    print(f"Guardado en {out_file}")

if __name__ == "__main__":
    main()
