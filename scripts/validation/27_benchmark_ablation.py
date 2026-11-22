"""Benchmark de ablación: mide el efecto de quitar componentes del método."""
import json, random, gzip
from pathlib import Path

def rank_genomes(meta_db, rna_db):
    """Rankea las bacterias según su completitud metabólica y de ARN."""
    scores = {}
    for gid in set(list(meta_db.keys()) + list(rna_db.keys())):
        # Puntuación: 1 punto por ruta metabólica + 0.1 por ARN
        paths = len(meta_db.get(gid, {}))
        rnas = sum(len(v) for v in rna_db.get(gid, {}).values())
        scores[gid] = paths + (rnas * 0.1)
        
    # Ordenar de mayor a menor
    return sorted(scores.items(), key=lambda x: x[1], reverse=True)

def strip_gff_annotations(gff_path, out_path, drop_rate=0.3):
    """Elimina aleatoriamente el X% de las anotaciones CDS para simular un genoma 'huérfano'."""
    dropped_genes = []
    kept_lines = []
    
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    with opener(gff_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                kept_lines.append(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9 and parts[2] == "CDS":
                if random.random() < drop_rate:
                    dropped_genes.append(line)
                    continue # Saltar esta línea (borrarla)
            kept_lines.append(line)
            
    with open(out_path, 'w') as f:
        f.writelines(kept_lines)
        
    return dropped_genes

def main():
    root = Path.cwd()
    meta_db = json.load(open(root / "data/metabolism/global_metabolism.json"))
    rna_db = json.load(open(root / "data/rna_annotations/rna_annotations.json"))
    
    print("1. Rankeando las 400 bacterias por completitud...")
    ranking = rank_genomes(meta_db, rna_db)
    
    print("\n🏆 Top 5 Bacterias mejor anotadas:")
    for gid, score in ranking[:5]:
        print(f"  - {gid} (Score: {score:.1f})")
        
    # Seleccionamos la mejor bacteria para el experimento
    best_gid = ranking[0][0]
    print(f"\n2. Seleccionada para el Benchmark: {best_gid}")
    
    # Buscar su GFF
    gff_path = list((root / "data/annotations").glob(f"{best_gid}*.gff.gz"))[0]
    stripped_gff_path = root / "data/annotations" / f"{best_gid}_STRIPPED.gff"
    
    print("\n3. Aplicando 'Ablation' (Borrando el 30% de los genes al azar)...")
    dropped = strip_gff_annotations(gff_path, stripped_gff_path, drop_rate=0.3)
    print(f"  Se borraron {len(dropped)} anotaciones de genes de la base de datos.")
    
    print("\n4. Resultados del Test (La prueba para el inversionista):")
    print(f"  - Hemos creado un genoma artificial '{stripped_gff_path.name}' al que le faltan {len(dropped)} genes.")
    print("  - Si metemos este genoma en nuestro Dashboard, el metabolismo se verá roto (brechas rojas).")
    print("  - El objetivo de nuestra IA es escanear el ADN crudo, encontrar esos genes faltantes y volver a cerrar las brechas.")
    print("\n✅ Benchmark generado con éxito. Listo para la comparación visual.")

if __name__ == "__main__":
    main()
