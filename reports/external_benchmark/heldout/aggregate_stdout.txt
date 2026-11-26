# Benchmark held-out: Bakta, Prokka y DFAST

**Generado:** 2026-11-26T02:44:44Z  
**Alcance:** Resumen agregado de un benchmark held-out por genoma para Bakta, Prokka y DFAST contra anotaciones de referencia locales.

**Genomas held-out:** `cjejuni`, `saureus`, `spneumoniae`  
**Set de tuning excluido:** `cglutamicum`, `mtb_h37rv`, `paeruginosa`, `pputida`

**Veredicto:** Sin ganador: la dispersión de F1 micro es 0.02 puntos porcentuales (≤ 1 pp). Las diferencias no permiten declarar superioridad.

## Resumen agregado por herramienta

| Herramienta | Versión | Genomas ok | TP | FP | FN | Precisión | Recall | F1 micro | F1 medio | Exactos | Tasa exacta | Bordes distintos | Tiempo total (s) | Tiempo medio (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Bakta | 1.12.2 | 3/3 | 6164 | 419 | 303 | 0.9364 | 0.9531 | 0.9447 | 0.9453 | 5663 | 91.9% | 501 | 1964.7 | 654.9 |
| Prokka | 1.15.6 | 3/3 | 6164 | 416 | 303 | 0.9368 | 0.9531 | 0.9449 | 0.9455 | 5625 | 91.3% | 539 | 1451.7 | 483.9 |
| DFAST | 1.6.2 | 3/3 | 6163 | 415 | 304 | 0.9369 | 0.9530 | 0.9449 | 0.9455 | 5603 | 90.9% | 560 | 852.0 | 284.0 |

> TP/FP/FN y F1 micro se calculan sumando los conteos de todos los genomas. La tasa exacta es `exact_matches / TP`; un par puede ser TP con bordes distintos. El tiempo es la suma de `runtime_seconds` de cada herramienta.

## Resultado por genoma

| Genoma | Herramienta | Estado | TP | FP | FN | P | R | F1 | Exactos | Bordes distintos | Tiempo (s) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cjejuni | Bakta | `ok` | 1608 | 102 | 16 | 0.9404 | 0.9901 | 0.9646 | 1488 | 120 | 545.8 |
| cjejuni | Prokka | `ok` | 1609 | 102 | 15 | 0.9404 | 0.9908 | 0.9649 | 1487 | 122 | 375.6 |
| cjejuni | DFAST | `ok` | 1608 | 101 | 16 | 0.9409 | 0.9901 | 0.9649 | 1488 | 120 | 219.9 |
| saureus | Bakta | `ok` | 2640 | 72 | 202 | 0.9735 | 0.9289 | 0.9507 | 2341 | 299 | 862.0 |
| saureus | Prokka | `ok` | 2639 | 69 | 203 | 0.9745 | 0.9286 | 0.9510 | 2327 | 312 | 656.0 |
| saureus | DFAST | `ok` | 2639 | 69 | 203 | 0.9745 | 0.9286 | 0.9510 | 2316 | 323 | 348.5 |
| spneumoniae | Bakta | `ok` | 1916 | 245 | 85 | 0.8866 | 0.9575 | 0.9207 | 1834 | 82 | 556.9 |
| spneumoniae | Prokka | `ok` | 1916 | 245 | 85 | 0.8866 | 0.9575 | 0.9207 | 1811 | 105 | 420.1 |
| spneumoniae | DFAST | `ok` | 1916 | 245 | 85 | 0.8866 | 0.9575 | 0.9207 | 1799 | 117 | 283.6 |

## Agregado por categoría

| Herramienta | Categoría | TP | FP | FN | P | R | F1 | Exactos |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Bakta | CDS | 5967 | 416 | 303 | 0.9348 | 0.9517 | 0.9432 | 5483 |
| Bakta | tRNA | 160 | 3 | 0 | 0.9816 | 1.0000 | 0.9907 | 157 |
| Bakta | rRNA | 37 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 23 |
| Prokka | CDS | 5967 | 412 | 303 | 0.9354 | 0.9517 | 0.9435 | 5482 |
| Prokka | tRNA | 160 | 4 | 0 | 0.9756 | 1.0000 | 0.9877 | 120 |
| Prokka | rRNA | 37 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 23 |
| DFAST | CDS | 5966 | 411 | 304 | 0.9355 | 0.9515 | 0.9435 | 5483 |
| DFAST | tRNA | 160 | 4 | 0 | 0.9756 | 1.0000 | 0.9877 | 120 |
| DFAST | rRNA | 37 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |

## TP / FP / FN

- **Bakta**: TP=6164, FP=419, FN=303 (precisión 0.9364, recall 0.9531, F1 0.9447).
- **Prokka**: TP=6164, FP=416, FN=303 (precisión 0.9368, recall 0.9531, F1 0.9449).
- **DFAST**: TP=6163, FP=415, FN=304 (precisión 0.9369, recall 0.9530, F1 0.9449).

## Genomas excluidos del agregado

- `btheta`: Genoma multi-contig: Bakta renombró los contigs a contig_1/contig_2 y DFAST a sequence1/sequence2; el runner no empareja seqids genéricos con más de un contig, por lo que ambos quedaron not_comparable. Se conserva como limitación y no entra en el agregado. Estado por herramienta: Bakta=not_comparable, Prokka=ok, DFAST=not_comparable.
  - Prokka: F1 0.9639, TP/FP/FN 4725/284/70.

## Entradas y procedencia

- Imagen: `sha256:2fdb291dc37749aebde733a5dc32e9ca4d4460b3456dcc6dbb4b1d7732d6307b`.
- Hilos: [4]; timeout: [5400.0] s; umbral de solapamiento: [0.5].
- `cjejuni` — Campylobacter jejuni (GCF_000009085.1); FASTA `1a0c7e733310de29a6c715cf5e6825727f9ba3c9718b504cd5f5af03591599d9`; referencia `592fae15135b0ca61d29b40c06f21fcc0f16bb1e37dbb85da78ecb6760b8a7e1`.
- `saureus` — Staphylococcus aureus (GCF_000013425.1); FASTA `f30fbb37ecacc1c6baca951ada10065a39a29e019c81f664a5d1fca3d8818adb`; referencia `b39f4aaab06bcdfecd9080e8605039fc76e18f8f0ecb3043cba02c38b69f1571`.
- `spneumoniae` — Streptococcus pneumoniae (GCF_001457635.1); FASTA `e5e580af13c7a8ff81e795fe656615ddc0495982e6afc6c7284a7df8c08a67c6`; referencia `5724f05787121edc09f41212753e002d930258c4a35b76f0cf1843204ba607ae`.

## Limitaciones

- La referencia es el gold de esta comparación; sus criterios y errores se heredan.
- El emparejamiento es 1:1 voraz por solapamiento recíproco y no un óptimo bipartito.
- La coincidencia literal de 'product' no mide función; los nombres dependen de cada base.
- No se comparan ncRNA, tmRNA, pseudogenes (si se excluyen), dominios ni rutas.
- Una herramienta ausente no es un cero: queda como not_available sin métricas.
- La Bakta DB se detecta localmente; no se descarga y su ausencia deja Bakta en not_available.
- No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.
- Las ubicaciones compuestas de GenBank se comparan por su intervalo envolvente.
- El emparejamiento es 1:1 voraz por solapamiento recíproco al umbral indicado; no es un óptimo bipartito y otro umbral cambiaría TP/FP/FN.
- La coincidencia de 'product' no mide función; los nombres dependen de cada base.
- Los tiempos son de una sola ejecución secuencial por genoma en el mismo host; no hay réplicas ni intervalo de confianza.
- El F1 agregado es micro (TP/FP/FN sumados); el F1 medio es la media aritmética por genoma.
- El resumen no declara superioridad entre herramientas.
