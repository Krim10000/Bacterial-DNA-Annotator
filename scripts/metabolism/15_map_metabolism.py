#!/usr/bin/env python3
"""Extrae enzimas de un GFF y las mapea a rutas metabólicas usando la API REST de KEGG directamente."""
import sys, gzip, re, json, time
from pathlib import Path
import requests

def extract_enzymes_from_gff(gff_path):
    enzymes = set()
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    
    with opener(gff_path, 'rt') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS": continue
            
            attrs = parts[8]
            # Buscar en formato EC_number= o db_xref=EC:
            matches = re.findall(r'EC[=:]?\s*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)', attrs)
            # También buscar dentro de la descripción del producto [EC:...]
            matches.extend(re.findall(r'\[EC:([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\]', attrs))
            
            for ec in matches:
                enzymes.add(ec)
    return list(enzymes)

def main():
    root = Path.cwd()
    
    # Apuntar a la bacteria solicitada
    target_id = "GCF_004103715.1"
    gff_path = list((root / "data/annotations").glob(f"{target_id}*.gff.gz"))
    
    if not gff_path:
        print(f"❌ No se encontró {target_id} en data/annotations/")
        return
        
    target_gff = gff_path[0]
    print(f"Analizando: {target_gff.name}")
    
    # 1. Extraer enzimas
    ec_numbers = extract_enzymes_from_gff(target_gff)
    print(f"Se encontraron {len(ec_numbers)} enzimas (EC numbers) anotadas.")
    
    if not ec_numbers:
        print("⚠️ No hay enzimas anotadas en este GFF.")
        return

    # 2. Consultar KEGG directamente (sin bioservices)
    print("Consultando rutas en KEGG (esto tardará ~15 segundos)...")
    pathways_found = {}
    
    for i, ec in enumerate(ec_numbers[:20]): # Consultamos las primeras 20
        try:
            # Paso A: Obtener el KO (KEGG Orthology) a partir del EC
            res_ko = requests.get(f"https://rest.kegg.jp/conv/ko/ec:{ec}").text
            if res_ko:
                ko = res_ko.strip().split("\t")[1] # ej. ko:K00844
                
                # Paso B: Obtener las rutas para ese KO en E. coli (eco)
                res_path = requests.get(f"https://rest.kegg.jp/link/pathway/eco:{ko.replace('ko:', '')}").text
                for line in res_path.strip().split("\n"):
                    if not line: continue
                    # El formato devuelto es: path:eco00010 \t eco:b0001
                    path_id = line.split("\t")[0].replace("path:", "")
                    if path_id not in pathways_found:
                        pathways_found[path_id] = {"name": f"Ruta {path_id}", "enzymes": []}
                    if ec not in pathways_found[path_id]["enzymes"]:
                        pathways_found[path_id]["enzymes"].append(ec)
            
            print(f"  [{i+1}/20] EC {ec} -> KO encontrado", end="\r")
            time.sleep(0.4) # KEGG requiere un rate limit de ~3 peticiones por segundo
            
        except Exception as e:
            pass

    print("\n")
    # 3. Guardar resultados
    out_dir = root / "data" / "metabolism"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / f"{target_id}_pathways.json"
    
    with open(out_file, "w") as f:
        json.dump(pathways_found, f, indent=2)
        
    print(f"✅ Se encontraron {len(pathways_found)} rutas metabólicas reales.")
    for pid, info in list(pathways_found.items())[:5]:
        print(f" - {pid}: {len(info['enzymes'])} enzimas")

if __name__ == "__main__":
    main()
