#!/usr/bin/env python3
"""Busca sistemas CRISPR y Profagos en los GFF de NCBI."""
import gzip, re, json
from pathlib import Path
import sys

def analyze_defense_systems(gff_path):
    """Busca proteínas relacionadas con fagos, CRISPR y defensa."""
    results = {
        "crispr": [],
        "prophage": [],
        "defense": []
    }
    
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS": continue
                
                start, end = int(parts[3]), int(parts[4])
                attrs = parts[8]
                
                match = re.search(r'product=([^;]+)', attrs)
                if not match: continue
                
                prod = match.group(1).lower()
                
                # 1. Buscar CRISPR/Cas
                if "crispr" in prod or "cas9" in prod or "cas3" in prod or "cas1" in prod or "cas2" in prod:
                    results["crispr"].append({"start": start, "end": end, "name": prod})
                    
                # 2. Buscar Profagos (Integrasas, cápsides, colas)
                elif any(k in prod for k in ["integrase", "capsid", "terminase", "tail fiber", "phage", "prophage"]):
                    results["prophage"].append({"start": start, "end": end, "name": prod})
                    
                # 3. Sistemas de Defensa (Restricción, Toxina-Antitoxina)
                elif any(k in prod for k in ["restriction", "methyltransferase", "toxin", "antitoxin", "abortive"]):
                    results["defense"].append({"start": start, "end": end, "name": prod})
                    
    except Exception:
        pass
        
    return results

def main():
    root = Path.cwd()
    ann_dir = root / "data" / "annotations"
    out_file = root / "data" / "defense_systems.json"
    
    gff_files = list(ann_dir.glob("*.gff.gz"))
    print(f"Analizando sistemas de defensa en {len(gff_files)} genomas...")
    
    global_results = {}
    
    for i, gff in enumerate(gff_files):
        genome_id = gff.stem.replace(".gff", "")
        res = analyze_defense_systems(gff)
        
        # Solo guardar si encontramos algo
        if any(res.values()):
            global_results[genome_id] = res
            
        if (i+1) % 50 == 0:
            print(f"  [{i+1}/{len(gff_files)}] procesados...")
            
    with open(out_file, "w") as f:
        json.dump(global_results, f, indent=2)
        
    print(f"\n✅ Análisis guardado en {out_file}")
    
    total_crispr = sum(len(v["crispr"]) for v in global_results.values())
    total_phage = sum(len(v["prophage"]) for v in global_results.values())
    print(f"Se encontraron {total_crispr} proteínas CRISPR y {total_phage} proteínas de profagos en total.")

if __name__ == "__main__":
    main()
