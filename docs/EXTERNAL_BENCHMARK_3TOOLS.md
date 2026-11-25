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
