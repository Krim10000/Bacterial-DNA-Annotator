# Benchmark externo de tres herramientas

## Gold 60 kb

Gold independiente:

```text
benchmarks/external/gold/chlorobaculum_tepidum_60kb
55 CDS + 1 tRNA
```

Comando:

```bash
docker run --rm -w /work \
  -v "$PWD:/work" \
  -v "$PWD/data/db/bakta_db:/data/bakta_db:ro" \
  -v "$PWD/data/db/dfast_db:/data/dfast_db:ro" \
  dna-token-bac/benchmark-external:cpu \
  --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools all --run --threads 1 --bakta-db /data/bakta_db \
  --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db
```

## Resultado

| Herramienta | TP | FP | FN | Precisión | Recall | F1 | Tiempo | RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Bakta 1.12.2 light | 52 | 4 | 4 | 0,9286 | 0,9286 | **0,9286** | 125,6 s | 2,65 GB |
| Prokka 1.15.6 | 54 | 2 | 2 | 0,9643 | 0,9643 | **0,9643** | 45,5 s | 2,65 GB |
| DFAST 1.6.2 | 52 | 4 | 4 | 0,9286 | 0,9286 | **0,9286** | 69,8 s | 2,65 GB |

Por categoría:

| Herramienta | CDS F1 | tRNA F1 |
|---|---:|---:|
| Bakta | 0,9273 | 1,0000 |
| Prokka | 0,9636 | 1,0000 |
| DFAST | 0,9369 | 0,0000 |

DFAST no detectó el tRNA del gold.

## Segundo genoma

*C. glutamicum*:

| Herramienta | TP | FP | FN | F1 | Tiempo |
|---|---:|---:|---:|---:|---:|
| Bakta | 2.980 | 163 | 45 | 0,9663 | 2.052 s |
| Prokka | 2.980 | 163 | 45 | 0,9663 | 1.780 s |
| DFAST | 2.980 | 163 | 45 | 0,9663 | 740 s |

Coinciden en TP/FP/FN; cambian los bordes exactos.

## Lectura

En el gold pequeño **Prokka gana por 2 CDS**. En un genoma real las tres
herramientas empatan en F1. No hay evidencia para declarar superioridad
global.

## Límites

- Dos genomas no son un benchmark multi-especie.
- El gold pequeño no contiene rRNA.
- El F1 mezcla CDS/tRNA/rRNA.
- La comparación es geométrica, no funcional.
- Bakta usó DB `light`; `full` puede cambiar resultados.

## Cuatro genomas

| Genoma | Bakta F1 | Prokka F1 | DFAST F1 |
|---|---:|---:|---:|
| C. glutamicum | 0,9663 | 0,9663 | 0,9663 |
| M. tuberculosis H37Rv | 0,9485 | 0,9462 | 0,9474 |
| P. aeruginosa | 0,9848 | 0,9830 | 0,9840 |
| P. putida | 0,9809 | 0,9799 | 0,9805 |
| **Media** | **0,9701** | **0,9688** | **0,9696** |

No hay ganador único. La diferencia entre las tres es ≤0,005 F1.

Informe detallado:

```text
reports/external_benchmark/multi_v2/summary.md
```

## Held-out (3 genomas)

Genomas no usados en tuning: *S. aureus*, *C. jejuni*, *S. pneumoniae*.

| Herramienta | TP | FP | FN | F1 micro | Exactos | Tiempo total |
|---|---:|---:|---:|---:|---:|---:|
| Bakta | 6.164 | 419 | 303 | 0,9447 | 91,9 % | 1.964,7 s |
| Prokka | 6.164 | 416 | 303 | **0,9449** | 91,3 % | 1.451,7 s |
| DFAST | 6.163 | 415 | 304 | **0,9449** | 90,9 % | 852,0 s |

Diferencia entre las tres: **0,02 puntos porcentuales**. No se declara ganador.

`btheta` quedó fuera del agregado porque Bakta/DFAST renombran contigs en genomas multi-contig. Prokka dio F1 0,9639.

Informe: `reports/external_benchmark/heldout/summary.md`.
