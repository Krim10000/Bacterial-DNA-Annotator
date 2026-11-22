"""Utilidades compartidas para anotar y etiquetar genomas."""
import torch
import torch.nn as nn
from lib.dna_model import DNAModel

class DNAAnnotator(DNAModel):
    def __init__(self, vocab_size, num_classes=5, **kwargs):
        super().__init__(vocab_size, **kwargs)
        # Capa de clasificación multi-etiqueta
        # Usamos BCEWithLogitsLoss (no Softmax) para permitir solapamientos
        self.annotation_head = nn.Linear(kwargs.get('d_model', 768), num_classes)
        
    def forward(self, input_ids, labels=None):
        x = self.embed(input_ids)
        cos, sin = self.rope(x, x.shape[1])
        
        for block in self.blocks:
            x = block(x, cos, sin)
            
        x = self.norm(x)
        
        # Si estamos entrenando para anotar, devolvemos los logits de las etiquetas
        # en lugar de las predicciones de vocabulario.
        if labels is not None:
            logits = self.annotation_head(x)
            loss = nn.BCEWithLogitsLoss()(logits.view(-1, logits.size(-1)), labels.view(-1, labels.size(-1)))
            return loss, logits
        else:
            return self.annotation_head(x)
