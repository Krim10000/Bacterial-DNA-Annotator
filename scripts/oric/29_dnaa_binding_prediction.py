#!/usr/bin/env python3
"""
Predice la secuencia de unión de DnaA basándose en su estructura 3D y la región upstream.
"""
import gzip, re, json, requests
from pathlib import Path

def get_esmfold_structure(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        res = requests.post(url, data=sequence, timeout=30)
        if res.status_code == 200:
            return res.text
    except:
        pass
    return None

def find_dnaa_gene(gff_path):
    """Busca el gen dnaA y devuelve sus coordenadas."""
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if "product=chromosomal replication initiator protein DnaA" in line or "gene=dnaA" in line:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        return int(parts[3]), int(parts[4])
    except:
        pass
    return None

def extract_upstream_region(fasta_path, dnaa_start, upstream_len=500):
    """Extrae la región antes del gen dnaA (donde está el oriC)."""
    with gzip.open(fasta_path, "rt") as f:
        seq = "".join([line.strip() for line in f if not line.startswith(">")]).upper()
    
    # Si el gen está en la hebra forward, tomamos lo de antes
    start = max(0, dnaa_start - upstream_len)
    return seq[start:dnaa_start]

def find_dnaa_boxes(seq):
    """Busca cajas DnaA con tolerancia a mutaciones."""
    # Patrón: TTAT[C|T]CAC[G|A|T] (permite variaciones)
    pattern = r"TTAT[C|T]CAC[G|A|T]"
    return [m.group(0) for m in re.finditer(pattern, seq)]

def main():
    root = Path.cwd()
    raw_dir = root / "data" / "raw"
    ann_dir = root / "data" / "annotations"
    
    # Tomamos 5 bacterias de ejemplo
    fasta_files = list(raw_dir.glob("*.fna.gz"))[:5]
    
    print("🧬 Analizando especificidad de unión de DnaA en 5 bacterias...")
    print("=" * 60)
    
    results = {}
    
    for fasta in fasta_files:
        gid = fasta.name.replace("_genomic.fna.gz", "")
        gff_path = ann_dir / f"{gid}.gff.gz"
        
        # 1. Encontrar el gen DnaA
        coords = find_dnaa_gene(gff_path)
        if not coords:
            print(f"  ❌ {gid}: No se encontró el gen dnaA.")
            continue
            
        start, end = coords
        print(f"\n📌 {gid}")
        print(f"   Gen dnaA en posición: {start} - {end}")
        
        # 2. Extraer la región upstream (donde se une)
        upstream = extract_upstream_region(fasta, start)
        
        # 3. Buscar motivos de unión
        boxes = find_dnaa_boxes(upstream)
        print(f"   Motivos de unión encontrados en la región upstream: {boxes}")
        
        results[gid] = {
            "dnaa_start": start,
            "upstream_region": upstream[:100], # Guardar solo los primeros 100 pb
            "binding_motifs": boxes
        }
        
    # Guardar resultados
    out_file = root / "data" / "dnaa_binding_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n" + "=" * 60)
    print(f"✅ Resultados guardados en {out_file}")
    print("\n💡 Interpretación:")
    print("   Si dos bacterias tienen motivos de unión similares, sus DnaA")
    print("   probablemente tienen estructuras 3D similares en el dominio de unión.")

if __name__ == "__main__":
    main()
