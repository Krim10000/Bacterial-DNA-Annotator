#!/usr/bin/env python3
"""Anota el oriC (documentado vs predicho) y reporta las estadísticas."""
import gzip, re, json
from pathlib import Path

def find_oric_in_gff(gff_path):
    """Busca si NCBI ya documentó el origen de replicación."""
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if "rep_origin" in line or "origin of replication" in line.lower():
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        return {"start": int(parts[3]), "end": int(parts[4]), "source": "documented_NCBI", "confidence": 1.0}
    except: pass
    return None

def predict_oric(seq):
    """Predice el oriC buscando clústeres de cajas DnaA (TTATCCACA)."""
    # Patrón con tolerancia a mutaciones biológicas (1 error)
    pattern = r"TTAT[C|T]CAC[G|A|T]" 
    boxes = [m.start() for m in re.finditer(pattern, seq)]
    
    # Buscamos clústeres de al menos 3 cajas en menos de 500 pb
    for i in range(len(boxes) - 2):
        if boxes[i+2] - boxes[i] < 500:
            num_boxes = 0
            for b in boxes:
                if boxes[i]-100 <= b <= boxes[i+2]+100: num_boxes += 1
            
            confidence = 0.7 if num_boxes == 3 else 0.9
            return {
                "start": max(0, boxes[i] - 200), 
                "end": boxes[i+2] + 200, 
                "source": "predicted_script", 
                "confidence": confidence,
                "boxes_found": num_boxes
            }
    return None

def main():
    root = Path.cwd()
    raw_dir = root / "data" / "raw"
    ann_dir = root / "data" / "annotations"
    out_file = root / "data" / "architecture_annotations.json"
    
    fasta_files = list(raw_dir.glob("*.fna.gz"))
    print(f"Anotando oriC (Documentado vs Predicho) en {len(fasta_files)} genomas...")
    
    results = {}
    documented_count = 0
    predicted_count = 0
    not_found_count = 0
    
    for i, fasta in enumerate(fasta_files):
        gid = fasta.name.replace("_genomic.fna.gz", "")
        gff_path = ann_dir / f"{gid}.gff.gz"
        
        # 1. Buscar en la documentación de NCBI
        oric = find_oric_in_gff(gff_path) if gff_path.exists() else None
        
        # 2. Si no está documentado, predecir con el script
        if not oric:
            with gzip.open(fasta, "rt") as f:
                seq = "".join([line.strip() for line in f if not line.startswith(">")]).upper()
            oric = predict_oric(seq)
            if oric:
                predicted_count += 1
            else:
                not_found_count += 1
        else:
            documented_count += 1
            
        if oric:
            results[gid] = {"oriC": oric}
            
        if (i+1) % 100 == 0:
            print(f"  [{i+1}/{len(fasta_files)}] genomas procesados...")
            
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n" + "="*50)
    print("📊 RESUMEN DEL ANÁLISIS DE oriC")
    print("="*50)
    print(f"Total de bacterias analizadas: {len(fasta_files)}")
    print(f"✅ Documentados en NCBI (Confianza 100%): {documented_count}")
    print(f"🧠 Predichos por nuestro script (Confianza 70-90%): {predicted_count}")
    print(f"❌ No encontrados (Región atípica): {not_found_count}")
    print("="*50)
    print(f"Resultados detallados guardados en {out_file}")

if __name__ == "__main__":
    main()
