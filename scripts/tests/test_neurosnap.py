#!/usr/bin/env python3
"""Prueba de humo de la integración con la API de Neurosnap (ESMFold)."""
import json, requests, time

API_KEY = "659ab970c5e4520aa65c6b324985a414992fa2c1b6ff0f894287159c9bb930740d10ed4d726c313753d8569814d664af876184ba32768ae84424668578033b44"  # <--- PON TU CLAVE AQUÍ
SECUENCIA = "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"

HEADERS = {"X-API-KEY": API_KEY}

def enviar_trabajo():
    print("🚀 Enviando secuencia a Neurosnap (ESMFold)...")
    
    # 1. Estructura JSON para la secuencia (FFSequenceList)
    sequences = {"aa": {"proteina_prueba": SECUENCIA}}
    
    # 2. Forzar multipart/form-data. El campo se llama "Input Sequences" (con S)
    multipart_fields = [
        ("Input Sequences", (None, json.dumps(sequences)))
    ]
    
    url = "https://neurosnap.ai/api/job/submit/ESMFold?note=test_api"
    response = requests.post(url, headers=HEADERS, files=multipart_fields)
    
    if response.status_code == 200:
        job_id = response.json()
        print(f"✅ Trabajo enviado. ID: {job_id}")
        return job_id
    else:
        print(f"❌ Error al enviar: {response.status_code} - {response.text}")
        return None

def esperar_resultado(job_id):
    print("⏳ Esperando a que Neurosnap termine de plegar la proteína...")
    url = f"https://neurosnap.ai/api/job/status/{job_id}"
    
    while True:
        response = requests.get(url, headers=HEADERS)
        if response.status_code == 200:
            status = response.json()
            
            if status == "completed":
                print("✅ ¡Plegamiento completado!")
                return True
            elif status in ["failed", "cancelled", "deleted"]:
                print(f"❌ El trabajo falló con estado: {status}")
                return False
            else:
                print(f"  ... Estado: {status}. Esperando 5s...")
                time.sleep(5)
        else:
            print(f"❌ Error al consultar estado: {response.status_code}")
            return False

def descargar_resultado(job_id):
    print("📥 Buscando archivos de salida...")
    url_data = f"https://neurosnap.ai/api/job/data/{job_id}"
    response = requests.get(url_data, headers=HEADERS)
    
    if response.status_code == 200:
        data = response.json()
        archivos_salida = data.get("out", [])
        
        # Buscar el archivo .pdb en la lista de salidas
        pdb_filename = None
        for file_info in archivos_salida:
            if file_info[0].endswith(".pdb"):
                pdb_filename = file_info[0]
                break
                
        if pdb_filename:
            print(f"📄 Archivo encontrado: {pdb_filename}")
            url_file = f"https://neurosnap.ai/api/job/file/{job_id}/out/{pdb_filename}"
            file_response = requests.get(url_file, headers=HEADERS)
            
            if file_response.status_code == 200:
                with open("betalactamasa_neurosnap.pdb", "wb") as f:
                    f.write(file_response.content)
                print("✅ Archivo descargado y guardado como 'betalactamasa_neurosnap.pdb'")
            else:
                print(f"❌ Error al descargar archivo: {file_response.status_code}")
        else:
            print("⚠️ No se encontró ningún archivo .pdb en los resultados.")
    else:
        print(f"❌ Error al obtener datos del trabajo: {response.status_code}")

if __name__ == "__main__":
    if API_KEY == "TU_API_KEY_AQUI":
        print("⚠️ ERROR: Debes editar el script y poner tu API_KEY real.")
    else:
        job_id = enviar_trabajo()
        if job_id:
            if esperar_resultado(job_id):
                descargar_resultado(job_id)
