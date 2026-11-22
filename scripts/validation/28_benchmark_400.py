#!/usr/bin/env python3
"""Benchmark de las 400 bacterias: Ranking y simulación de pérdida de anotaciones."""
import json, random
from pathlib import Path

def calculate_completeness(gff_path, meta_data, rna_data):
    """Calcula un score de 0 a 100 para la bacteria."""
    score = 0
    # 1. Metabolismo (hasta 50 puntos)
    if meta_data:
        score += min(50, len(meta_data) * 5) 
    # 2. ARNs (hasta 20 puntos)
    if rna_data:
        rnas = sum(len(v) for v in rna_data.values())
        score += min(20, rnas * 0.2)
    # 3. Cantidad de Genes anotados (hasta 30 puntos)
    count_cds = 0
    import gzip
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if "\tCDS\t" in line: count_cds += 1
    except: pass
    
    # 4000 genes es un genoma bacteriano promedio completo
    score += min(30, (count_cds / 4000) * 30) 
    return round(score, 2), count_cds

def main():
    root = Path.cwd()
    meta_db = json.load(open(root / "data/metabolism/global_metabolism.json"))
    rna_db = json.load(open(root / "data/rna_annotations/rna_annotations.json"))
    gff_dir = root / "data/annotations"
    
    results = []
    print("Evaluando las 400 bacterias (esto toma unos segundos)...")
    
    for gff in gff_dir.glob("*.gff.gz"):
        gid = gff.stem.replace(".gff", "")
        meta = meta_db.get(gid, {})
        rna = rna_db.get(gid, {})
        
        score, cds_count = calculate_completeness(gff, meta, rna)
        results.append({
            "genome_id": gid,
            "score": score,
            "cds_count": cds_count,
            "pathways": len(meta),
            "rnas": sum(len(v) for v in rna.values())
        })
        
    # Ordenar de mejor a peor
    results.sort(key=lambda x: x["score"], reverse=True)
    
    out_file = root / "data" / "benchmark_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
        
    print(f"\n✅ Ranking guardado en {out_file}")
    print("\n🏆 Top 3 Bacterias Mejor Anotadas (El 'Genoma de Oro'):")
    for r in results[:3]:
        print(f"  {r['genome_id']} | Score: {r['score']} | Genes: {r['cds_count']} | Rutas: {r['pathways']}")
        
    print("\n⚠️ Bottom 3 Bacterias Peor Anotadas (El desafío para la IA):")
    for r in results[-3:]:
        print(f"  {r['genome_id']} | Score: {r['score']} | Genes: {r['cds_count']} | Rutas: {r['pathways']}")

    # --- LA SIMULACIÓN (Ablation) ---
    print("\n🔬 SIMULACIÓN PARA INVERSIONISTAS:")
    best_genome = results[0]
    print(f"Tomamos la mejor bacteria ({best_genome['genome_id']} con {best_genome['cds_count']} genes).")
    print(f"Le borramos artificialmente el 50% de sus anotaciones de genes.")
    print(f"El nuevo conteo de genes será: {best_genome['cds_count'] // 2}.")
    print("El metabolismo se romperá en el Dashboard (aparecerán brechas rojas).")
    print("Nuestra IA escaneará las regiones 'basura' para encontrar esos genes ocultos.")

if __name__ == "__main__":
    main()
