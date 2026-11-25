# Siguiente corrida: b02_c01

## Archivo

`data/metabolism/targets/free_chunks/b02_c01.tsv`

- 25 proteínas.
- Longitudes 458–600 aa.
- ~12.910 aa totales.
- ETA estimada: **25–40 min**.

## Kaggle

1. Subir `b02_c01.tsv` al dataset `possible-enzimes`.
2. Abrir el notebook.
3. Cambiar:

```python
BATCH = "/kaggle/input/datasets/camilorojas00/possible-enzimes/b02_c01.tsv"
```

4. Ejecutar el modelo una vez.
5. Ejecutar la función.
6. `RUN_FULL = True`.
7. ZIP y descargar.

## No repetir

- No reutilizar el `metrics.tsv` del lote 01.
- Nuevo `OUT` o borrar el anterior.
- No cargar el modelo dos veces.

## Resultado esperado

```text
esmfold_out/pdb/*.pdb
esmfold_out/metrics.tsv
```
