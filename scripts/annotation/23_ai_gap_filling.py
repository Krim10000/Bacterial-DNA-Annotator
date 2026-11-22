#!/usr/bin/env python3
"""Usa el modelo de IA para calcular la perplejidad y encontrar genes ocultos."""
import torch, numpy as np, sentencepiece as spm, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.dna_model import DNAModel

def calculate_perplexity(model, sp, seq):
    """Calcula cuán 'sorprendido' está el modelo con esta secuencia."""
    ids = sp.encode(seq, out_type=int)[:2048]
    ids = ids + [0] * (2048 - len(ids))
    input_ids = torch.tensor([ids])
    
    with torch.no_grad():
        logits = model(input_ids)
        # Calcular la pérdida (CrossEntropy) de predecir los tokens
        shift_logits = logits[..., :-1, :].contiguous()
        shift_labels = input_ids[..., 1:].contiguous()
        loss_fct = torch.nn.CrossEntropyLoss(ignore_index=0)
        loss = loss_fct(shift_logits.view(-1, shift_logits.size(-1)), shift_labels.view(-1))
        return torch.exp(loss).item() # Perplejidad

def main():
    root = Path.cwd()
    tokenizer_path = root / "data" / "tokenizer" / "dna_sp.model"
    # Busca cualquier checkpoint que tengas descargado
    checkpoints = list((root / "checkpoints_descargados").rglob("*.pt"))
    
    if not checkpoints:
        print("⚠️ No se encontró ningún modelo .pt. Descarga los checkpoints de Modal primero.")
        return
        
    model_path = checkpoints[0]
    print(f"Cargando modelo de IA: {model_path.name}")
    
    sp = spm.SentencePieceProcessor(model_file=str(tokenizer_path))
    model = DNAModel(vocab_size=sp.get_piece_size())
    model.load_state_dict(torch.load(model_path, map_location="cpu"))
    model.eval()

    # Simulamos una brecha metabólica: una región intergénica donde FALTA una enzima.
    # En un caso real, esto vendría del Dashboard.
    # Secuencia de ADN de prueba (puedes cambiarla por una región real de tu bacteria)
    test_sequences = {
        "Región Basura (Control Negativo)": "NNNNNNNNNNACGTACGTNNNNNNNNNNACGTACGTNNNNNNNNNN",
        "Gen Conocido (Control Positivo)": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAA",
        "Candidato a Gen Oculto (Brecha)": "ATGACCGTATTGCAGGCTGGTAAAGCGCTGAAAGATATCGCTAAAGAGCTGGGCGTGGAAGTTATCGAAACCGGTGAAGTGAAAGCGTAA"
    }

    print("\n--- Analizando Perplejidad (Menor es mejor) ---")
    for name, seq in test_sequences.items():
        ppl = calculate_perplexity(model, sp, seq)
        print(f"  {name}: Perplejidad = {ppl:.2f}")
        
    print("\n💡 Interpretación:")
    print("  Si 'Candidato a Gen Oculto' tiene una perplejidad baja (cercana al Gen Conocido),")
    print("  el modelo confirma que ahí hay un gen real, aunque NCBI diga que es basura.")

if __name__ == "__main__":
    main()
