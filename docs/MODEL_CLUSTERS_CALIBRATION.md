# Calibración de los 281.090 clusters del modelo K-12

Herramienta: `tools/calibrate_model_clusters.py`
Informe: `reports/model_clusters_calibration.md` · JSON: `reports/model_clusters_calibration.json`

No estima precisión: mide solapamiento y prioriza revisión.

## Números

| métrica | valor |
|---|---:|
| clusters | 281.090 |
| confianza alta | 106.547 |
| confianza media | 9.126 |
| confianza baja | 165.417 |
| clusters de 2 subunidades | 87.110 (31,0 %) |
| con al menos un gen de transferencia | 19.467 (6,9 %) |
| con intervalo colocalizado (≤50 kb) | 7.500 (2,7 %) |
| P1 (score ≥4,0) | 207 |
| P2 (3,0–4,0) | 1.350 |

## Lectura

- La confianza alta está dominada por clusters de dos subunidades.
- El solapamiento con transferencias es bajo: la mayoría de los 281.090
  clusters no tiene ancla K-12.
- El intervalo amplio es engañoso: 93,7 % toca un intervalo, pero solo 2,7 %
  está colocalizado.
- La priorización P1/P2 es operativa, no biológica.
- No hay gold independiente para las otras especies.

## Regla

`complex_candidate` = hipótesis. No demuestra interacción, estequiometría ni
función. No cierra brechas.
