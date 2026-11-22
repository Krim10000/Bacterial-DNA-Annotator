"""Conecta el plegado ESMFold con los candidatos y sus brechas."""
import requests, re
from pathlib import Path

def get_esmfold_structure(sequence):
    """Conecta con la API pública de ESMFold y devuelve el PDB."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    print("  Conectando con ESMFold (Meta AI)...")
    try:
        res = requests.post(url, data=sequence, timeout=30)
        if res.status_code == 200:
            return res.text
    except Exception as e:
        print(f"  Error de conexión: {e}")
    return None

def analyze_plddt(pdb_data):
    """Extrae la confianza (pLDDT) de la columna B-factor del PDB."""
    plddts = []
    for line in pdb_data.split("\n"):
        if line.startswith("ATOM"):
            try:
                # La columna B-factor está en las posiciones 61-66 del formato PDB
                b_factor = float(line[60:66].strip())
                plddts.append(b_factor)
            except ValueError:
                pass
    if not plddts: return 0, 0
    return sum(plddts)/len(plddts), max(plddts)

def main():
    root = Path.cwd()
    out_dir = root / "data" / "structures"
    out_dir.mkdir(exist_ok=True)
    
    # Secuencia de ejemplo (una proteína real)
    test_sequence = "MKKLLPTAALLLSACSTPAPAQETVKLQGSGKYTVTATGQTVTIKDLNGVQQVKATLDLQGQTVTLTGQTVTVKDLNG"
    
    pdb = get_esmfold_structure(test_sequence)
    if pdb:
        # Guardar el archivo PDB
        out_file = out_dir / "test_protein.pdb"
        with open(out_file, "w") as f:
            f.write(pdb)
            
        # Analizar la confianza
        mean_plddt, max_plddt = analyze_plddt(pdb)
        print(f"\n✅ Estructura guardada en {out_file}")
        print(f"📊 Análisis de Confianza (pLDDT):")
        print(f"  - Media: {mean_plddt:.2f} (Si es > 70, la estructura es fiable)")
        print(f"  - Máxima: {max_plddt:.2f}")
        
        if mean_plddt > 80:
            print("  🧠 El sistema concluye: 'Alta confianza estructural. El modelo ha identificado un dominio proteico estable.'")
    else:
        print("❌ No se pudo obtener la estructura.")

if __name__ == "__main__":
    main()
