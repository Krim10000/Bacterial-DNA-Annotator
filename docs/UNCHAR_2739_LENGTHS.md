# 2.739 proteínas sin AFDB ordenadas por longitud

Archivo ordenado:

`data/metabolism/targets/UNCHAR_2739_BY_LENGTH.tsv`

Columnas:

```text
uniprot  genome  protein_id  length  clusters  description
```

## Distribución

| longitud | proteínas | clusters | longitud media |
|---|---:|---:|---:|
| 0–150 | 677 | 1.463 | 105 |
| 151–300 | 1.008 | 3.609 | 223 |
| 301–500 | 750 | 3.133 | 386 |
| 501–800 | 231 | 910 | 602 |
| 801–1.024 | 45 | 204 | 889 |
| >1.024 | 28 | 106 | 1.642 |

Total: **2.739 proteínas**, **9.425 clusters**.

## Estimación ESMFold

Supuestos:

- 5–30 s por proteína de ~300 aa en A100.
- Escalado aproximado cuadrático con la longitud.
- ESMFold requiere GPU CUDA ≥16 GB.
- Las >1.024 aa necesitan fragmentación o modelo alternativo.

| grupo | GPU-h bajo | GPU-h medio | GPU-h alto |
|---|---:|---:|---:|
| ≤1.024 aa (2.711) | 4,6 | 13,7 | 27,3 |
| >1.024 aa (28) | 2,6 | 7,7 | 15,4 |
| **Total** | **7,2** | **21,4** | **42,7** |

Coste con A100 40GB:

- $2–3,5/h.
- **Bajo:** ~$14.
- **Medio:** ~$43–75.
- **Alto:** ~$85–150.

Presupuesto práctico: **$40–90** para las 2.739, con margen.

## Notas

- Las 28 proteínas >1.024 aa no se pliegan completas de una vez.
- El coste real depende de la GPU, longitud y reintentos.
- El ThinkPad con Intel UHD/Iris Xe no sirve: no hay CUDA.
