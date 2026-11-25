# Empezar en Kaggle o Colab

## Kaggle

1. Crear notebook.
2. Subir `batch_01.tsv` como dataset o en `/kaggle/input`.
3. Abrir `notebooks/ESMFold_Free_Batches.ipynb`.
4. Settings → Accelerator → **GPU T4 x2** o **P100**.
5. Internet → **On**.
6. Ejecutar la celda de instalación.
7. En la celda de ejecución, dejar `[:5]` para probar.
8. Verificar `metrics.tsv` y un PDB.
9. Cambiar a lote completo.
10. Descargar `esmfold_out.zip`.

## Colab

1. Subir el notebook y `batch_01.tsv` a `/content/`.
2. Runtime → Change runtime type → **T4 GPU**.
3. Ejecutar instalación.
4. `BATCH` se resuelve a `/content/batch_01.tsv`.
5. Probar 5; luego lote completo.
6. Descargar ZIP.

## Orden de lotes

```text
batch_01
batch_02
batch_03
...
batch_10
```

## Regla

- Probar 5 primero.
- No lanzar los 10 lotes a la vez.
- Kaggle ~30 h/semana.
- Descargar resultados.
- No pagar nada.
