# Trabajo en bloques pequeños

## Regla

- 25 proteínas por bloque.
- ~25–30 min por bloque.
- Descargar ZIP al terminar.
- Siguiente bloque.

## Archivos

`data/metabolism/targets/free_chunks/`

```text
b01_c01.tsv
b01_c02.tsv
...
b10_c11.tsv
```

Índice:

`data/metabolism/targets/CHUNKS_INDEX.tsv`

## Uso

En el notebook, cambia solo:

```python
BATCH = "/kaggle/input/.../b01_c01.tsv"
```

o en Colab:

```python
BATCH = "/content/b01_c01.tsv"
```

Luego:

```python
RUN_FULL = True
```

## Reanudar

El runner salta PDB existentes:

```python
if target.exists() and target.stat().st_size > 1000:
    continue
```

## Si un bloque falla

- descargar lo que haya;
- subir el ZIP a la siguiente sesión;
- descomprimir;
- repetir el mismo bloque.

## Ventaja

- máximo 25 proteínas perdidas;
- descarga frecuente;
- Kaggle no se queda 4 h sin checkpoint.
