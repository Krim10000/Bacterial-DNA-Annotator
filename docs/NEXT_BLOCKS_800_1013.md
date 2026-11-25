# Bloques restantes: 800 + 1.013

## Estado

Ya hechos:

- 257 lote 01;
- 150 b02_c03–c08;
- 300 b02_c09 en adelante.

Total cortas hechas: **757**.  
Restantes ≤600 aa: **1.813**.

## Bloque A — 800

`NEXT_800.tsv`

- 800 proteínas únicas.
- 206.766 aa.
- ETA: **5,5–7 h**.
- Incluye:
  - resto de b03_c10;
  - b03_c11;
  - b04–b06;
  - primeros 15 de b07.

## Bloque B — 1.013

`REST_1013.tsv`

- 1.013 proteínas únicas.
- 254.861 aa.
- ETA: **6,5–8 h**.
- Incluye:
  - resto de b07;
  - b08–b10.
- Al terminar: **0 cortas pendientes**.

## Después

Quedan:

- 169 proteínas >600 aa.

## Kaggle

Para 800:

```python
BATCH = "/kaggle/input/.../NEXT_800.tsv"
OUT = Path("/kaggle/working/esmfold_800")
```

Para 1.013:

```python
BATCH = "/kaggle/input/.../REST_1013.tsv"
OUT = Path("/kaggle/working/esmfold_1013")
```

## Advertencia

800 en una sola celda son ~6 h.  
Si la sesión muere, pierdes todo.  
Descarga el ZIP al terminar; si puedes, parte el bloque en 2–3 subtandas.
