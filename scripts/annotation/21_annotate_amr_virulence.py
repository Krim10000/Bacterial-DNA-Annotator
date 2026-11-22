#!/usr/bin/env python3
"""Busca genes de resistencia a antibióticos (AMR) y virulencia en los GFF de NCBI."""
import gzip, re, json
from pathlib import Path

# Palabras clave basadas en bases de datos como CARD, VFDB y anotaciones de NCBI
AMR_KEYWORDS = [
    "beta-lactamase", "efflux", "transporter", "resistance", "antibiotic",
    "tetracycline", "ampC", "mecA", "vanA", "aminoglycoside", "chloramphenicol"
]
VIRULENCE_KEYWORDS = [
    "toxin", "hemolysin", "adhesin", "invasion", "virulence", "pathogenicity",
    "type iii effector", "type iv secretion", "capsule", "siderophore", "pili"
]

def classify_genes(gff_path):
    """Busca palabras clave en el atributo 'product' del GFF."""
    results = {"amr": [], "virulence": []}
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS": continue
                
                attrs = parts[8]
                match = re.search(r'product=([^;]+)', attrs)
                if not match: continue
                
                prod = match.group(1).lower()
                start, end = int(parts[3]), int(parts[4])
                
                # Detectar AMR
                if any(kw in prod for kw in AMR_KEYWORDS):
                    results["amr"].append({"start": start, "end": end, "name": match.group(1)})
                # Detectar Virulencia
                elif any(kw in prod for kw in VIRULENCE_KEYWORDS):
                    results["virulence"].append({"start": start, "end": end, "name": match.group(1)})
    except Exception:
        pass
    return results

def main():
    root = Path.cwd()
    ann_dir = root / "data" / "annotations"
    out_file = root / "data" / "defense_virulence.json"
    
    gff_files = list(ann_dir.glob("*.gff.gz"))
    print(f"Analizando AMR y Virulencia en {len(gff_files)} genomas...")
    
    global_db = {}
    for i, gff in enumerate(gff_files):
        genome_id = gff.stem.replace(".gff", "")
        res = classify_genes(gff)
        
        if any(res.values()):
            global_db[genome_id] = res
            
        if (i+1) % 100 == 0:
            print(f"  [{i+1}/{len(gff_files)}] procesados...")

    with open(out_file, "w") as f:
        json.dump(global_db, f, indent=2)
        
    total_amr = sum(len(v["amr"]) for v in global_db.values())
    total_vir = sum(len(v["virulence"]) for v in global_db.values())
    print(f"\n✅ ¡Análisis completado! Guardado en {out_file}")
    print(f"Se encontraron {total_amr} genes de resistencia y {total_vir} genes de virulencia.")

if __name__ == "__main__":
    main()
