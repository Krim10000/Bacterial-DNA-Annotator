# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-25T18:47:51Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 --tools all --run --threads 1 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence5 --json /work/reports/external_benchmark/benchmark5.json --md /work/reports/external_benchmark/benchmark5.md
```

Este informe no declara superioridad entre herramientas. Las métricas comparan features contra la referencia local; no evalúan función, dominios ni rutas.

## Entorno

- Python: `3.11.18` (`/opt/conda/bin/python`)
- Plataforma: `Linux 6.8.0-149-generic x86_64`
- CPU lógicas: `8`
- Umbral de solapamiento recíproco: **0.50**
- `--run`: **sí**
- Pseudogenes excluidos: **sí**
- Normalización de seqid: **sí**

### Bakta DB

- Estado: `ok` (válida: **sí**).
- Ruta: `/data/bakta_db` (origen: `argument`).
- Versión: `0` (fuente: `version.json`).
- Marcadores: `bakta.db, version.json`.

| herramienta | binario | disponible | versión | PATH / razón |
|---|---|---|---|---|
| Bakta | `bakta` | sí | 1.12.2 | /usr/local/bin/bakta |
| Prokka | `prokka` | sí | 1.15.6 | /usr/local/bin/prokka |
| DFAST | `dfast` | sí | 1.6.2 | /usr/local/bin/dfast |

- Herramientas disponibles: **bakta, prokka, dfast**.
- Herramientas ausentes: **ninguna**.

## Entradas

- FASTA: `benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna` (SHA256 `2d2dc79ed637655d…`, 1 contig, 60000 bp).
- Referencia: `benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3` (formato `gff3`, SHA256 `73d0aeadcfc2e236…`).
- Conteos de referencia: CDS **55**, tRNA **1**, rRNA **0**.
- CDS pseudo excluidos en referencia: **0**.
- SHA256 FASTA: `2d2dc79ed637655d5c12b816b248b472090cd48f395526c7adedb572416e0949`.
- SHA256 referencia: `73d0aeadcfc2e236eb30ae448f5788c45633806224a7f0450b8d19aa5caa4830`.

## Validación

- Inicio: `2026-11-25T18:43:43Z`; fin: `2026-11-25T18:47:51Z`; total: **247.80 s**.
- JSON declarado: `reports/external_benchmark/benchmark5.json`.
- Markdown declarado: `reports/external_benchmark/benchmark5.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 2d2dc79ed637655d5c12b816b248b472090cd48f395526c7adedb572416e0949 |
| reference_path | ok | referencia gff3 con 56 features |
| reference_sha256 | ok | SHA256 73d0aeadcfc2e236eb30ae448f5788c45633806224a7f0450b8d19aa5caa4830 |
| reference_features | ok | CDS=55 tRNA=1 rRNA=0 pseudo_excluidos=0 |
| tool_bakta | ok | binario ejecutable: /usr/local/bin/bakta |
| tool_prokka | ok | binario ejecutable: /usr/local/bin/prokka |
| tool_dfast | ok | binario ejecutable: /usr/local/bin/dfast |
| bakta_db | ok | estado=ok; /data/bakta_db |
| result_bakta | ok | estado=ok; sin incidencias |
| result_prokka | ok | estado=ok; sin incidencias |
| result_dfast | ok | estado=ok; sin incidencias |

## Resultados

| herramienta | fuente | estado | F1 CDS | F1 tRNA | F1 rRNA | precisión global | recall global | F1 global | segundos |
|---|---|---|---|---|---|---|---|---|---|
| Bakta | executed | ok | 0.9273 | 1.0000 | — | 0.9286 | 0.9286 | 0.9286 | 125.6 |
| Prokka | executed | ok | 0.9636 | 1.0000 | — | 0.9643 | 0.9643 | 0.9643 | 45.5 |
| DFAST | executed | ok | 0.9369 | 0.0000 | — | 0.9286 | 0.9286 | 0.9286 | 69.8 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 125.62 s, RSS máx. 2649.5 MiB (`2026-11-25T18:43:48Z` → `2026-11-25T18:45:53Z`).
- Salida: `/tmp/external_benchmark_bakta_v_sf3db5/out/benchmark.gff3` (72556 bytes, SHA256 `7064fcebc80c0b2f…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_v_sf3db5/out --prefix benchmark --threads 1 --force --skip-plot --db /data/bakta_db /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna`.
- stderr (cola): `/opt/conda/envs/bakta/lib/python3.13/site-packages/bakta/features/cds.py:51: UserWarning: sequence should be at least 100000 characters (60000 found)
  trainings_info = gene_finder.train(*seqs, translation_table=cfg.translation_table)`
- Anotación comparada: `reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence5/bakta/benchmark.gff3` (formato `gff3`, SHA256 `7064fcebc80c0b2f…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 55 | 51 | 4 | 4 | 0.9273 | 0.9273 | 0.9273 | 40 | 11 | 0.9783 |
| tRNA | 1 | 1 | 1 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 1 | 0 | 1.0000 |
| rRNA | 0 | 0 | 0 | 0 | 0 | — | — | — | 0 | 0 | — |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 51 | 51 | 0.2353 | 0.4633 | 1.0000 | 0.0000 |
| tRNA | 1 | 1 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 0 | 0 | — | — | — | — |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 55 | 4 | 7.3 | 55 | 55 | 10 | 18.2 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **4**.
- Falsos positivos (herramienta sin pareja): **4**.
- Pares solo con hebra distinta: **2**.
- Pares emparejados con algún borde distinto: **11**.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002932.3 | 3506..3958 | + | ribonuclease P protein component |
| FN | CDS | NC_002932.3 | 15585..16514 | - | outer membrane protein assembly factor BamD |
| FN | CDS | NC_002932.3 | 17554..18081 | - | F0F1 ATP synthase subunit B |
| FN | CDS | NC_002932.3 | 59083..59376 | - | FeoA family protein |
| FP | CDS | NC_002932.3 | 3635..3814 | - | hypothetical protein |
| FP | CDS | NC_002932.3 | 17605..17916 | + | hypothetical protein |
| FP | CDS | NC_002932.3 | 42293..42658 | - | hypothetical protein |
| FP | CDS | NC_002932.3 | 59137..59370 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002932.3 | 4017..4268 | -66 | +0 | 0.7925 |
| CDS | NC_002932.3 | 7605..9044 | -12 | +0 | 0.9917 |
| CDS | NC_002932.3 | 11570..12604 | +27 | +0 | 0.9739 |
| CDS | NC_002932.3 | 15001..15588 | +0 | -165 | 0.7194 |
| CDS | NC_002932.3 | 16990..17535 | +0 | -39 | 0.9286 |
| CDS | NC_002932.3 | 20271..21260 | +0 | +36 | 0.9649 |
| CDS | NC_002932.3 | 31609..32775 | +0 | +63 | 0.9488 |
| CDS | NC_002932.3 | 35348..36772 | +0 | -15 | 0.9895 |
| CDS | NC_002932.3 | 40306..40686 | +0 | +87 | 0.8141 |
| CDS | NC_002932.3 | 42801..43667 | +78 | +0 | 0.9100 |
| CDS | NC_002932.3 | 52338..54092 | +0 | -249 | 0.8581 |


## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 45.45 s, RSS máx. 2649.5 MiB (`2026-11-25T18:45:54Z` → `2026-11-25T18:46:39Z`).
- Salida: `/tmp/external_benchmark_prokka_2qu2cvpl/out/benchmark.gff` (74659 bytes, SHA256 `bd7a273938250c08…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_2qu2cvpl/out --prefix benchmark --cpus 1 --force --quiet /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna`.
- Anotación comparada: `reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence5/prokka/benchmark.gff` (formato `gff3`, SHA256 `bd7a273938250c08…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 55 | 53 | 2 | 2 | 0.9636 | 0.9636 | 0.9636 | 45 | 8 | 0.9847 |
| tRNA | 1 | 1 | 1 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 1 | 0.9865 |
| rRNA | 0 | 0 | 0 | 0 | 0 | — | — | — | 0 | 0 | — |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 53 | 53 | 0.2830 | 0.5199 | 0.8571 | 0.0000 |
| tRNA | 1 | 1 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 0 | 0 | — | — | — | — |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 55 | 4 | 7.3 | 55 | 55 | 9 | 16.4 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **2**.
- Falsos positivos (herramienta sin pareja): **2**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **9**.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002932.3 | 3285..3446 | + | 50S ribosomal protein L34 |
| FN | CDS | NC_002932.3 | 3506..3958 | + | ribonuclease P protein component |
| FP | CDS | NC_002932.3 | 3635..3814 | - | hypothetical protein |
| FP | CDS | NC_002932.3 | 42305..42838 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002932.3 | 4017..4268 | -66 | +0 | 0.7925 |
| CDS | NC_002932.3 | 7605..9044 | -12 | +0 | 0.9917 |
| CDS | NC_002932.3 | 11570..12604 | +27 | +0 | 0.9739 |
| CDS | NC_002932.3 | 15001..15588 | +0 | -165 | 0.7194 |
| CDS | NC_002932.3 | 31609..32775 | +0 | +63 | 0.9488 |
| CDS | NC_002932.3 | 35348..36772 | +0 | -15 | 0.9895 |
| CDS | NC_002932.3 | 40306..40686 | +0 | +87 | 0.8141 |
| CDS | NC_002932.3 | 41754..42221 | +0 | +21 | 0.9571 |
| tRNA | NC_002932.3 | 16805..16877 | -1 | +0 | 0.9865 |


## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 69.82 s, RSS máx. 2649.5 MiB (`2026-11-25T18:46:41Z` → `2026-11-25T18:47:51Z`).
- Salida: `/tmp/external_benchmark_dfast_z6e7xx31/out/genome.gff` (106958 bytes, SHA256 `a404f2421b2266c6…`).
- Comando: `/usr/local/bin/dfast --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna --out /tmp/external_benchmark_dfast_z6e7xx31/out --cpu 1 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence5/dfast/genome.gff` (formato `gff3`, SHA256 `a404f2421b2266c6…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 56 | 52 | 4 | 3 | 0.9286 | 0.9455 | 0.9369 | 40 | 12 | 0.9731 |
| tRNA | 1 | 0 | 0 | 0 | 1 | — | 0.0000 | 0.0000 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | 0 | 0 | — | — | — | 0 | 0 | — |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 52 | 52 | 0.8269 | 0.8714 | 1.0000 | 0.0000 |
| tRNA | 0 | 0 | — | — | — | — |
| rRNA | 0 | 0 | — | — | — | — |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 55 | 55 | 4 | 7.3 | 56 | 56 | 10 | 17.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **4**.
- Falsos positivos (herramienta sin pareja): **4**.
- Pares solo con hebra distinta: **2**.
- Pares emparejados con algún borde distinto: **12**.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002932.3 | 3506..3958 | + | ribonuclease P protein component |
| FN | CDS | NC_002932.3 | 17554..18081 | - | F0F1 ATP synthase subunit B |
| FN | CDS | NC_002932.3 | 59083..59376 | - | FeoA family protein |
| FN | tRNA | NC_002932.3 | 16805..16877 | - | tRNA-Ala |
| FP | CDS | NC_002932.3 | 3635..3814 | - | hypothetical protein |
| FP | CDS | NC_002932.3 | 17605..17916 | + | hypothetical protein |
| FP | CDS | NC_002932.3 | 42293..42658 | - | hypothetical protein |
| FP | CDS | NC_002932.3 | 59137..59370 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002932.3 | 4017..4268 | -66 | +0 | 0.7925 |
| CDS | NC_002932.3 | 7605..9044 | -12 | +0 | 0.9917 |
| CDS | NC_002932.3 | 11570..12604 | +27 | +0 | 0.9739 |
| CDS | NC_002932.3 | 15001..15588 | +0 | -165 | 0.7194 |
| CDS | NC_002932.3 | 15585..16514 | +0 | +378 | 0.7110 |
| CDS | NC_002932.3 | 16990..17535 | +0 | -39 | 0.9286 |
| CDS | NC_002932.3 | 20271..21260 | +0 | +36 | 0.9649 |
| CDS | NC_002932.3 | 31609..32775 | +0 | +63 | 0.9488 |
| CDS | NC_002932.3 | 35348..36772 | +0 | -15 | 0.9895 |
| CDS | NC_002932.3 | 40306..40686 | +0 | +87 | 0.8141 |
| CDS | NC_002932.3 | 42801..43667 | +78 | +0 | 0.9100 |
| CDS | NC_002932.3 | 52338..54092 | +0 | -249 | 0.8581 |


## Herramientas no disponibles

- Ninguna herramienta seleccionada quedó `not_available`.

## Limitaciones

- La referencia es el gold de esta comparación; sus criterios y errores se heredan.
- El emparejamiento es 1:1 voraz por solapamiento recíproco y no un óptimo bipartito.
- La coincidencia literal de 'product' no mide función; los nombres dependen de cada base.
- No se comparan ncRNA, tmRNA, pseudogenes (si se excluyen), dominios ni rutas.
- Una herramienta ausente no es un cero: queda como not_available sin métricas.
- La Bakta DB se detecta localmente; no se descarga y su ausencia deja Bakta en not_available.
- No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.
- Las ubicaciones compuestas de GenBank se comparan por su intervalo envolvente.

## Reproducción

Comando exacto con `--run` (Bakta requiere una Bakta DB local; no se descarga):

```bash
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
