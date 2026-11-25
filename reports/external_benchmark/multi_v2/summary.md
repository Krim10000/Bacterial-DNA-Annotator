# Benchmark multi-genoma CPU

Comparación de Bakta, Prokka y DFAST contra GFF de referencia. Cuatro genomas, un hilo o cuatro según corrida; la métrica no depende de hilos.

| Genoma | Herramienta | TP | FP | FN | F1 | Exactos | Tiempo (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| cglutamicum | bakta | 2980 | 163 | 45 | 0.9663 | 2580 | 2052 |
| cglutamicum | prokka | 2980 | 163 | 45 | 0.9663 | 2534 | 1780 |
| cglutamicum | dfast | 2980 | 163 | 45 | 0.9663 | 2516 | 740 |
| mtb_h37rv | bakta | 3832 | 294 | 122 | 0.9485 | 2960 | 1200 |
| mtb_h37rv | prokka | 3822 | 303 | 132 | 0.9462 | 2921 | 961 |
| mtb_h37rv | dfast | 3827 | 298 | 127 | 0.9474 | 2923 | 521 |
| paeruginosa | bakta | 5615 | 140 | 33 | 0.9848 | 5114 | 964 |
| paeruginosa | prokka | 5605 | 151 | 43 | 0.9830 | 5104 | 1056 |
| paeruginosa | dfast | 5613 | 147 | 35 | 0.9840 | 5111 | 405 |
| pputida | bakta | 5396 | 132 | 78 | 0.9809 | 5074 | 966 |
| pputida | prokka | 5391 | 138 | 83 | 0.9799 | 5066 | 763 |
| pputida | dfast | 5396 | 137 | 78 | 0.9805 | 5052 | 399 |

## Media F1

| Herramienta | Media F1 |
|---|---:|
| bakta | 0.9701 |
| prokka | 0.9688 |
| dfast | 0.9696 |

## Lectura

- No hay ganador único: las diferencias son ≤0,005 F1.
- Bakta tiene la media ligeramente mayor en estos cuatro genomas.
- Prokka gana el gold pequeño de 60 kb por 2 CDS.
- DFAST es el más rápido, pero falló el tRNA del gold.
- La comparación es geométrica: no evalúa función, rutas ni estructura.
- Faltan más genomas y un set held-out independiente.
