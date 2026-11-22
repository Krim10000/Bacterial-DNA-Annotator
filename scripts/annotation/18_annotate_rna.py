#!/usr/bin/env python3
"""Extrae coordenadas de rRNAs y tRNAs directamente de los GFF de NCBI."""
import gzip, json
from pathlib import Path

def parse_rna_from_gff(gff_path):
    features = {"rRNA": [], "tRNA": [], "ncRNA": []}
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9: continue
                
                ftype = parts[2]
                attrs = parts[8]
                start, end = int(parts[3]), int(parts[4])
                
                if "gbkey=tRNA" in attrs or ftype == "tRNA":
                    features["tRNA"].append({"start": start, "end": end})
                elif "gbkey=rRNA" in attrs or ftype == "rRNA":
                    features["rRNA"].append({"start": start, "end": end})
                elif "gbkey=ncRNA" in attrs or ftype == "ncRNA":
                    features["ncRNA"].append({"start": start, "end": end})
    except Exception:
        pass
    return features

def main():
    root = Path.cwd()
    ann_dir = root / "data" / "annotations"
    out_dir = root / "data" / "rna_annotations"
    out_dir.mkdir(exist_ok=True)
    
    gff_files = list(ann_dir.glob("*.gff.gz"))
    print(f"Extrayendo coordenadas exactas de ARNs desde {len(gff_files)} archivos...")
    
    global_rna_db = {}
    for i, gff in enumerate(gff_files):
        genome_id = gff.stem.replace(".gff", "")
        features = parse_rna_from_gff(gff)
        if any(features.values()):
            global_rna_db[genome_id] = features
        if (i+1) % 100 == 0:
            print(f"  [{i+1}/{len(gff_files)}] procesados...")

    summary_file = out_dir / "rna_annotations.json"
    with open(summary_file, "w") as f:
        json.dump(global_rna_db, f)
    print(f"\n✅ Coordenadas guardadas en {summary_file}")

if __name__ == "__main__":
    main()
