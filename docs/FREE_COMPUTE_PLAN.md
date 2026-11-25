# Plan de cómputo gratuito para 2.739 proteínas

## División

- **2.570** con longitud ≤600 aa.
- **169** con longitud >600 aa.
- 10 lotes de **257** proteínas, balanceados por longitud.
- Los 169 largos quedan aparte.

## Plataformas gratuitas

| plataforma | cuota | uso |
|---|---|---|
| Foldify Open | 600 aa, 5 jobs | lotes 01–04 cuando vuelva |
| Kaggle | ~30 GPU-h/semana | lotes 05–08 |
| Colab | T4 gratuita | lote 09 |
| Modal free | $30/mes | overflow / lote 10 |
| HF ZeroGPU | cuota pequeña | casos sueltos |

## Archivos

- `data/metabolism/targets/free_batches/batch_01.tsv ... batch_10.tsv`
- `batch_XX.faa`
- `data/metabolism/targets/FREE_169_LONG.tsv`
- `notebooks/ESMFold_Free_Batches.ipynb`
- `scripts/structure/run_esmfold_batch.py`

## Ejecución en Kaggle/Colab

1. Subir `batch_XX.tsv`.
2. Abrir `notebooks/ESMFold_Free_Batches.ipynb`.
3. Activar GPU T4/P100.
4. Ajustar `BATCH`.
5. Probar con 5.
6. Subir a 257.
7. Descargar `esmfold_out.zip`.

## Runner local/GPU

```bash
python scripts/structure/run_esmfold_batch.py \
  --batch data/metabolism/targets/free_batches/batch_01.tsv \
  --outdir /tmp/esmfold_batch01 \
  --limit 5
```

Es reanudable: salta PDB ya existentes.

## Tiempo

- 2.570 × 20–30 s ≈ 14–21 GPU-h.
- Kaggle: ~30 h/semana.
- Reparto realista: **3–6 semanas**.

## Reglas

- No crear cuentas múltiples.
- No automatizar en masa donde los términos lo prohíban.
- Respetar Foldify: 5 jobs concurrentes.
- Descargar resultados de Foldify antes de 7 días.
- No pagar Neurosnap.
- No usar MetaCenter sin afiliación checa.


## Subunidades

La evaluación de complejos queda congelada hasta terminar las 2.570 estructuras. Ver `IMPORTANT_NOTE.md`.
