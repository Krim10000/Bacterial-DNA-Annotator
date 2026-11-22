#!/usr/bin/env python3
"""Compara la predicción de oriC de ORCA con la verdad de referencia de DoriC."""
import gzip, joblib, subprocess, warnings, shutil, json
from pathlib import Path
from orcapy import ORCA

def download_and_process(gid, model, temp_dir):
    # 1. Construir URL de NCBI (Usando la ruta estándar)
    # Nota: Las rutas FTP se pueden obtener de tu selected.tsv, aquí usamos un ejemplo genérico
    # Para hacerlo rápido, usaremos el archivo local si ya lo tienes en annotations, sino descargaremos
    
    # En tu pipeline real, debes iterar sobre tu selected.tsv para tener las URLs exactas
    # Este es un ejemplo de cómo se vería para una bacteria
    pass 

def main():
    root = Path.cwd()
    model_path = root / "data/models/ORCA_RFC_model_1_4_0.pkl.gz"
    model = joblib.load(model_path)
    
    results = {}
    temp_gbk = root / "data" / "temp_orca.gbff"
    
    # Leemos las rutas desde nuestro archivo selected.tsv
    sel_file = root / "data" / "selected.tsv"
    if not sel_file.exists():
        print("❌ Falta selected.tsv")
        return
        
    with open(sel_file) as f:
        f.readline() # Header
        lines = f.readlines()
        
    print(f"Procesando {len(lines)} bacterias con ORCA (esto puede tardar)...")
    
    for i, line in enumerate(lines):
        parts = line.strip().split("\t")
        if len(parts) < 4: continue
        
        ftp_path = parts[3]
        basename = ftp_path.rstrip("/").split("/")[-1]
        gid = basename.split("_")[0] + "_" + basename.split("_")[1] # Ej: GCF_004103715.1
        
        # URL del GenBank
        gbk_url = f"{ftp_path}/{basename}_genomic.gbff.gz"
        
        try:
            # Descargar y descomprimir en memoria temporal
            subprocess.run(["wget", "-q", gbk_url, "-O", str(temp_gbk) + ".gz"], check=True)
            with gzip.open(str(temp_gbk) + ".gz", 'rb') as f_in:
                with open(temp_gbk, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            # Procesar con ORCA
            orca = ORCA.from_gbk(str(temp_gbk), model=model)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                orca.find_oriCs(show_info=False, plot_path=None)
                
            if hasattr(orca, 'predictions') and orca.predictions:
                best_prob = max(orca.predictions)
                if best_prob >= 0.5: # Umbral de confianza
                    best_idx = orca.predictions.index(best_prob)
                    best_pos = orca.oriC_middles[best_idx]
                    results[gid] = {"oriC_position": int(best_pos), "confidence": float(best_prob)}
                    print(f"  [{i+1}/{len(lines)}] ✅ {gid}: {best_pos} bp (Conf: {best_prob:.2f})")
                else:
                    print(f"  [{i+1}/{len(lines)}] ⚠️ {gid}: Baja confianza ({best_prob:.2f})")
            else:
                print(f"  [{i+1}/{len(lines)}] ⚠️ {gid}: Sin predicción")
                
        except Exception as e:
            print(f"  [{i+1}/{len(lines)}] ❌ {gid}: Error - {str(e)[:50]}")
            
        finally:
            # Limpiar archivos temporales para no llenar el disco
            if Path(str(temp_gbk) + ".gz").exists(): Path(str(temp_gbk) + ".gz").unlink()
            if temp_gbk.exists(): temp_gbk.unlink()

    out_file = root / "data" / "orca_predictions.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✅ Resultados guardados en {out_file}")

if __name__ == "__main__":
    main()
