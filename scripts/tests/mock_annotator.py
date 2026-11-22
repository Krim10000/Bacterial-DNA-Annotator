"""Generador de anotaciones sintéticas para probar el pipeline sin datos reales."""
import re
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from lib.io import read_fasta

def find_annotations(fasta_path, output_gff):
    annotations = []
    seq_name = "unknown"
    
    for name, seq in read_fasta(fasta_path):
        seq_name = name
        # Buscar Promotores (Caja -10)
        for match in re.finditer("TATAAT", seq):
            annotations.append({"type": "promoter", "start": match.start(), "end": match.end(), "strand": "+"})
            
        # Buscar ORFs reales (mínimo 300 pb para que se vean bien)
        # Solo buscamos los primeros 500 para no saturar el gráfico
        count = 0
        for match in re.finditer("ATG([ACGT]{300,2000}?)(TAA|TAG|TGA)", seq):
            if count > 500: break
            annotations.append({"type": "CDS", "start": match.start(), "end": match.end(), "strand": "+"})
            count += 1
            
    # Escribir GFF3
    with open(output_gff, "w") as f:
        f.write("##gff-version 3\n")
        for i, ann in enumerate(annotations):
            f.write(f"{seq_name}\tmock_AI\t{ann['type']}\t{ann['start']}\t{ann['end']}\t.\t{ann['strand']}\t.\tID={ann['type']}_{i}\n")
            
    return len(annotations), len(seq)

if __name__ == "__main__":
    root = Path.cwd()
    fasta_path = list((root / "data/raw").glob("*.fna.gz"))[0]
    out_path = root / "data" / "mock_annotations.gff3"
    
    num_anns, genome_size = find_annotations(fasta_path, out_path)
    print(f"✅ {num_anns} elementos encontrados. Tamaño del genoma: {genome_size} pb.")
