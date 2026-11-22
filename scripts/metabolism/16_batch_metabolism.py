#!/usr/bin/env python3
"""Mapea proteínas (por nombre) a rutas metabólicas usando búsqueda de texto en KEGG."""
import sys, gzip, re, json, time, requests
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

def extract_products_from_gff(gff_path):
    products = set()
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS": continue
                attrs = parts[8]
                match = re.search(r'product=([^;]+)', attrs)
                if match:
                    prod = match.group(1).replace('"', '').strip()
                    # Limpiar palabras comunes que estorban en la búsqueda
                    prod = re.sub(r'(?i)(putative|probable|hypothetical|protein|subunit|domain-containing)', '', prod)
                    prod = prod.strip()
                    if len(prod) > 3:
                        products.add(prod)
    except Exception:
        pass
    return list(products)

def main():
    root = Path.cwd()
    ann_dir = root / "data" / "annotations"
    out_dir = root / "data" / "metabolism"
    out_dir.mkdir(exist_ok=True)
    
    gff_files = list(ann_dir.glob("*.gff.gz"))
    print(f"Escaneando {len(gff_files)} genomas...")
    
    # 1. Recopilar todos los nombres de proteínas únicos de las 400 bacterias
    all_products = set()
    genome_products = {}
    
    for gff in gff_files:
        gid = gff.stem.replace(".gff", "")
        prods = extract_products_from_gff(gff)
        if prods:
            genome_products[gid] = prods
            all_products.update(prods)
            
    print(f"Se encontraron {len(all_products)} proteínas ÚNICAS. Buscando en KEGG (esto tardará ~{(len(all_products)*0.4)/60:.1f} minutos)...")
    
    # 2. Buscar en KEGG (con caché para no saturar)
    cache_file = out_dir / "kegg_cache.json"
    cache = json.load(open(cache_file)) if cache_file.exists() else {}
    
    for i, prod in enumerate(list(all_products)):
        if prod in cache: continue
        
        try:
            # Paso A: Buscar la proteína en E. coli (eco) por texto
            search_url = f"https://rest.kegg.jp/find/eco/{requests.utils.quote(prod)}"
            res = requests.get(search_url).text
            if res:
                # Tomar el primer resultado (ej. eco:b0184)
                gene_id = res.strip().split("\n")[0].split("\t")[0]
                
                # Paso B: Ver en qué rutas metabólicas participa ese gen
                path_url = f"https://rest.kegg.jp/link/pathway/{gene_id}"
                path_res = requests.get(path_url).text
                paths = [l.split("\t")[0].replace("path:", "") for l in path_res.strip().split("\n") if l]
                cache[prod] = paths
            else:
                cache[prod] = []
            time.sleep(0.35) # Respetar rate limit
            if (i+1) % 50 == 0:
                print(f"   [{i+1}/{len(all_products)}] procesadas...", end="\r")
        except Exception:
            cache[prod] = []
            
    with open(cache_file, "w") as f:
        json.dump(cache, f)
        
    # 3. Construir el JSON global cruzando bacterias con rutas
    print("\nConstruyendo base de datos metabólica global...")
    global_db = {}
    for gid, prods in genome_products.items():
        global_db[gid] = {}
        for prod in prods:
            for path in cache.get(prod, []):
                if path not in global_db[gid]:
                    global_db[gid][path] = []
                if prod not in global_db[gid][path]:
                    global_db[gid][path].append(prod)
                    
    out_file = out_dir / "global_metabolism.json"
    with open(out_file, "w") as f:
        json.dump(global_db, f)
        
    print(f"✅ ¡Éxito! Datos de {len(global_db)} bacterias guardados en {out_file}")
    print(f"Tamaño del archivo: {out_file.stat().st_size / 1024 / 1024:.2f} MB")

if __name__ == "__main__":
    main()
