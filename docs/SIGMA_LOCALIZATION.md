# Localización y fondo de la capa 1 sigma

**Fecha:** 2026-11-25  
**Alcance:** comparación de ventanas y fondos para la capa 1 (señal/localización)
del pipeline sigma en *E. coli* K-12 MG1655 y *B. subtilis* 168.  
**No se modificó** `scripts/sigma/sigma_annotator.py`, `modules/`, `webapp/` ni
`data/sigma/<genoma>/`.

## Resumen

- La implementación reproduce **exactamente** el baseline actual
  (`-80/+20`, 200 señuelos genómicos): K-12 49,48 % y S. subtilis 41,02 %
  (`delta = 0.0` contra `data/sigma/layers/<genoma>/layer_metrics.json`).
- Se barrieron 4 ventanas × 3 fondos × 2 genomas (24 combinaciones). Las
  métricas por sigma y por TSS están en `comparison.tsv`; el ranking top-1/top-3
  está en `assignment.tsv`.
- **Recomendación:** ventana **-60/+20** con fondo de **flancos locales**
  (`-500/-100` y `+100/+500`, 100 subventanas de igual largo por TSS).
  - Capa 1 macro (media de los dos genomas): **45,2 % → 53,3 % (+8,1 pp)**.
  - K-12: **49,5 % → 54,7 % (+5,2 pp)**; S. subtilis: **41,0 % → 52,0 % (+11,0 pp)**.
  - Asignación macro: top-1 **50,9 % → 51,2 % (+0,3 pp)**; top-3 cualquiera
    **87,9 % → 88,4 % (+0,55 pp)**. No hay mejora real de asignación; el cambio
    es de capa 1.
- **Alternativa conservadora:** **-60/+20 con fondo genómico**. Mantiene el
  fondo actual y sube capa 1 macro a **48,5 % (+3,2 pp)**, top-1 a **51,6 %
  (+0,7 pp)** y top-3 a **89,7 % (+1,9 pp)**.
- Los shuffles mononucleotídicos y las ventanas largas (-100/-120) **no** se
  recomiendan: el shuffle es demasiado estricto y las ventanas largas suben el
  soporte a costa de degradar top-1/top-3.

## Método

### Datos

| genoma | gold | TSS agrupados | etiquetas (TSS, sigma) | secuencia usada |
|---|---|---:|---:|---|
| `ecoli_k12` | `data/sigma/gold/ecoli_k12/tss_sigma.tsv` (RegulonDB) | 3189 | 3314 | `data/allostery/work/ecoli_k12.fna` (`NC_000913.3`) |
| `subtilis` | `data/sigma/gold/subtilis/tss_sigma.tsv` (DBTBS) | 2226 | 2226 | `data/dbtbs/genome.fna` (`basu168`) |

El gold se agrupó por `(contig, pos, strand)` con el mismo cargador de
`50_eval_layers.py`; un TSS con varios sigma se evalúa como multietiqueta.
`pos` es 1-based y el TSS es la primera base transcrita. La secuencia de K-12
usada aquí (FASTA) reproduce el baseline calculado con el GBFF, así que no hay
deriva de secuencia.

### Ventanas

`-up/+20` con `up ∈ {60, 80, 100, 120}`, en orientación de transcripción
(reverse-complemento en hebra `-`). La ventana incluye 20 nt aguas abajo del
TSS; el resto es aguas arriba.

### Fondos

| fondo | qué es | tamaño | semilla |
|---|---|---:|---|
| `genomic` | ventanas aleatorias del genoma completo, mismo largo que la ventana del TSS | 200 por (sigma, largo) | 17 + `crc32(sigma) % 1000` (idéntica a `50_eval_layers.py`) |
| `shuffle` | shuffles mononucleotídicos de **la ventana del TSS** | 100 por TSS | `crc32(genoma:tss:hebra:ventana:salt)` |
| `flank` | subventanas de igual largo muestreadas de los flancos locales `[-500,-100)` y `[+100,+500)` en orientación de transcripción | 100 por TSS (≈50 de cada flanco) | `crc32(genoma:tss:hebra:ventana:salt)` |

- En `flank` se muestrean inicios sin reemplazo cuando el flanco lo permite; si
  un flanco no alcanza, se rellena con el otro. En esta corrida todos los TSS
  tuvieron 100 subventanas (mínimo 100, sin fallback a shuffle).
- Para `-60/+20` la ventana central y los flancos son disjuntos. Para `-120/+20`
  el flanco aguas arriba se solapa 20 nt con el extremo distal de la ventana
  (no con las cajas -35/-10); se documenta como límite.

### Métrica

Para el sigma verdadero de cada TSS, con `score` = mejor score de
`score_sigma` en la ventana y `null` = distribución del mismo score en el fondo:

- `pvalue = (1 + #{null >= score}) / (1 + n_null)`
- `pct = #{null < score} / n_null`
- `z = (score - media(null)) / sd(null)` (`ddof=1`; `sd=0` usa 1)
- `soporte capa 1` = `pct >= 0,95` **y** `z >= 1,0` (definición actual de
  `layer1_localization`).
- `p<=0,05` y `pct>=0,95` se reportan por separado.
- **Global por TSS:** un TSS cuenta como soportado si **algún** sigma verdadero
  cumple la regla. El `z` global de la tabla es el máximo `z` entre los sigma
  verdaderos del TSS.
- **Top-1/top-3:** se reutiliza el ranking del evaluador
  `(z, -pvalue, score)` sobre todos los sigma; top-3 = primeros tres.

Los porcentajes de fondos locales son **condicionales a cada TSS** y no son
directamente intercambiables con los del fondo genómico global.

### Reproducción

```bash
# Barrido completo (24 combinaciones)
.venv/bin/python scripts/sigma/70_local_background.py --all --jobs 8

# Recalcular recomendación/salidas desde un comparison.json existente
.venv/bin/python scripts/sigma/70_local_background.py \
  --from-json data/sigma/localization/comparison.json

# Pruebas
.venv/bin/python -m pytest tests/test_sigma_local.py -q
.venv/bin/python scripts/sigma/70_local_background.py --selftest
```

## Resultados globales

`n` = TSS evaluados. `p<=0.05`, `pct>=0.95` y `soporte` son fracciones de TSS;
`z mediano` es la mediana del máximo `z` verdadero por TSS. `top-1` y `top-3 any`
son sobre el ranking de todos los sigma.

### E. coli K-12

| ventana | fondo | n | p<=0.05 | pct>=0.95 | soporte capa 1 | z mediano | top-1 | top-3 any |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| -60/+20 | genomic | 3189 | 51,1 % | 52,3 % | 52,3 % | 1,74 | 57,6 % | 91,3 % |
| -60/+20 | shuffle | 3189 | 34,2 % | 37,9 % | 37,9 % | 1,31 | 45,6 % | 81,0 % |
| -60/+20 | flank | 3189 | 53,2 % | 54,7 % | 54,7 % | 1,84 | 55,1 % | 89,2 % |
| -80/+20 | genomic | 3189 | 47,4 % | 49,5 % | 49,5 % | 1,68 | 57,0 % | 90,2 % |
| -80/+20 | shuffle | 3189 | 30,4 % | 33,8 % | 33,8 % | 1,19 | 42,9 % | 78,4 % |
| -80/+20 | flank | 3189 | 53,8 % | 54,5 % | 54,4 % | 1,75 | 50,7 % | 88,3 % |
| -100/+20 | genomic | 3189 | 43,6 % | 45,9 % | 45,9 % | 1,65 | 55,2 % | 89,5 % |
| -100/+20 | shuffle | 3189 | 27,9 % | 31,3 % | 31,3 % | 1,13 | 40,8 % | 76,4 % |
| -100/+20 | flank | 3189 | 55,1 % | 56,0 % | 55,8 % | 1,70 | 48,9 % | 86,8 % |
| -120/+20 | genomic | 3189 | 41,4 % | 42,1 % | 42,1 % | 1,57 | 53,3 % | 89,2 % |
| -120/+20 | shuffle | 3189 | 26,3 % | 29,2 % | 29,2 % | 1,04 | 39,3 % | 74,9 % |
| -120/+20 | flank | 3189 | 57,0 % | 57,5 % | 57,3 % | 1,71 | 46,2 % | 84,9 % |

### B. subtilis

| ventana | fondo | n | p<=0.05 | pct>=0.95 | soporte capa 1 | z mediano | top-1 | top-3 any |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| -60/+20 | genomic | 2226 | 42,0 % | 44,6 % | 44,6 % | 1,41 | 45,6 % | 88,2 % |
| -60/+20 | shuffle | 2226 | 28,2 % | 31,5 % | 31,5 % | 1,03 | 39,8 % | 74,3 % |
| -60/+20 | flank | 2226 | 50,4 % | 52,0 % | 52,0 % | 1,70 | 47,3 % | 87,6 % |
| -80/+20 | genomic | 2226 | 39,9 % | 41,0 % | 41,0 % | 1,29 | 44,8 % | 85,5 % |
| -80/+20 | shuffle | 2226 | 25,2 % | 28,1 % | 28,1 % | 0,92 | 37,5 % | 72,7 % |
| -80/+20 | flank | 2226 | 51,4 % | 52,3 % | 52,2 % | 1,63 | 45,5 % | 85,8 % |
| -100/+20 | genomic | 2226 | 40,3 % | 40,7 % | 40,7 % | 1,22 | 41,5 % | 83,9 % |
| -100/+20 | shuffle | 2226 | 24,3 % | 26,8 % | 26,8 % | 0,86 | 36,1 % | 72,1 % |
| -100/+20 | flank | 2226 | 52,9 % | 53,8 % | 53,5 % | 1,60 | 43,5 % | 84,7 % |
| -120/+20 | genomic | 2226 | 32,5 % | 33,2 % | 33,2 % | 1,13 | 37,7 % | 82,3 % |
| -120/+20 | shuffle | 2226 | 22,2 % | 24,8 % | 24,8 % | 0,82 | 35,2 % | 70,7 % |
| -120/+20 | flank | 2226 | 53,4 % | 54,4 % | 54,2 % | 1,57 | 42,5 % | 82,7 % |

## Comparación macro y regla de recomendación

La regla implementada maximiza la media macro de `soporte capa 1` entre los dos
genomas, pero **rechaza** combinaciones que empeoren top-1 o top-3 macro más de
1 punto porcentual respecto de la referencia actual (`-80/+20`, `genomic`).

| ventana | fondo | soporte macro | p<=0.05 macro | top-1 macro | top-3 macro | ¿pasa guardarraíl? |
|---|---|---:|---:|---:|---:|---|
| -80/+20 | genomic | 45,2 % | 43,7 % | 50,9 % | 87,9 % | sí (referencia) |
| -60/+20 | genomic | 48,5 % | 46,5 % | 51,6 % | 89,7 % | sí |
| **-60/+20** | **flank** | **53,3 %** | **51,8 %** | **51,2 %** | **88,4 %** | **sí (elegida)** |
| -80/+20 | flank | 53,3 % | 52,6 % | 48,1 % | 87,0 % | no |
| -100/+20 | flank | 54,7 % | 54,0 % | 46,2 % | 85,7 % | no |
| -120/+20 | flank | 55,8 % | 55,2 % | 44,3 % | 83,8 % | no |
| -120/+20 | genomic | 37,6 % | 37,0 % | 45,5 % | 85,7 % | no |

La combinación `-60/+20 flank` y `-80/+20 flank` empatan prácticamente en
soporte macro (53,3 %), pero `-80/+20` baja top-1 a 48,1 % y queda fuera. Entre
las aceptadas, `-60/+20 flank` tiene el mayor soporte.

## Detalle por sigma de la recomendación (-60/+20, flank)

### E. coli K-12

| sigma | n | p<=0.05 | pct>=0.95 | soporte | z mediano | precision | recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| sigma70 | 1997 | 50,7 % | 52,4 % | 52,4 % | 1,74 | 0,869 | 0,536 |
| sigma54 | 99 | 29,3 % | 29,3 % | 29,3 % | 0,80 | 0,159 | 0,394 |
| sigma32 | 312 | 64,1 % | 65,1 % | 65,1 % | 2,22 | 0,483 | 0,663 |
| sigma38 | 241 | 44,8 % | 46,1 % | 46,1 % | 1,64 | 0,238 | 0,664 |
| sigma28 | 143 | 35,0 % | 37,8 % | 37,8 % | 1,44 | 0,277 | 0,497 |
| sigma24 | 522 | 60,5 % | 61,5 % | 61,5 % | 1,95 | 0,552 | 0,638 |

### B. subtilis

| sigma | n | p<=0.05 | pct>=0.95 | soporte | z mediano | precision | recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| sigma70 | 1949 | 55,7 % | 57,6 % | 57,5 % | 1,86 | 0,916 | 0,514 |
| sigma54 | 11 | 54,5 % | 54,5 % | 54,5 % | 1,96 | 0,040 | 0,455 |
| sigma32 | 13 | 7,7 % | 7,7 % | 7,7 % | 0,13 | 0,024 | 0,231 |
| sigma38 | 182 | 9,9 % | 9,9 % | 9,9 % | 0,26 | 0,050 | 0,159 |
| sigma28 | 71 | 15,5 % | 15,5 % | 15,5 % | 0,36 | 0,119 | 0,197 |
| sigma24 | 0 | — | — | — | — | — | — |

### Cambio por sigma respecto del baseline

| sigma | K-12 -80/genomic | K-12 -60/flank | Δ K-12 | S. sub -80/genomic | S. sub -60/flank | Δ S. sub |
|---|---:|---:|---:|---:|---:|---:|
| sigma70 | 48,1 % | 52,4 % | +4,3 pp | 45,9 % | 57,5 % | +11,6 pp |
| sigma54 | 28,3 % | 29,3 % | +1,0 pp | 36,4 % | 54,5 % | +18,2 pp |
| sigma32 | 50,0 % | 65,1 % | +15,1 pp | 7,7 % | 7,7 % | +0,0 pp |
| sigma38 | 49,4 % | 46,1 % | -3,3 pp | 3,8 % | 9,9 % | +6,0 pp |
| sigma28 | 38,5 % | 37,8 % | -0,7 pp | 9,9 % | 15,5 % | +5,6 pp |
| sigma24 | 52,3 % | 61,5 % | +9,2 pp | — | — | — |

La mejora global está traccionada por `sigma70` en ambos genomas, `sigma24` y
`sigma32` en K-12, y `sigma70`/`sigma54`/`sigma38` en S. subtilis. En K-12,
`sigma38` y `sigma28` bajan levemente.

### Soporte por sigma y ventana con fondo flank

| sigma | genoma | n | -60/+20 | -80/+20 | -100/+20 | -120/+20 |
|---|---|---:|---:|---:|---:|---:|
| sigma70 | K-12 | 1997 | 52,4 % | 52,5 % | 54,0 % | 55,3 % |
| sigma54 | K-12 | 99 | 29,3 % | 31,3 % | 32,3 % | 34,3 % |
| sigma32 | K-12 | 312 | 65,1 % | 64,1 % | 64,1 % | 64,7 % |
| sigma38 | K-12 | 241 | 46,1 % | 45,2 % | 48,5 % | 50,6 % |
| sigma28 | K-12 | 143 | 37,8 % | 34,3 % | 36,4 % | 39,2 % |
| sigma24 | K-12 | 522 | 61,5 % | 61,1 % | 61,3 % | 62,8 % |
| sigma70 | S. sub | 1949 | 57,5 % | 57,7 % | 58,7 % | 59,2 % |
| sigma54 | S. sub | 11 | 54,5 % | 63,6 % | 63,6 % | 54,5 % |
| sigma32 | S. sub | 13 | 7,7 % | 7,7 % | 7,7 % | 7,7 % |
| sigma38 | S. sub | 182 | 9,9 % | 10,4 % | 13,7 % | 17,0 % |
| sigma28 | S. sub | 71 | 15,5 % | 16,9 % | 19,7 % | 21,1 % |

El soporte sube con ventanas más largas, pero la asignación cae (ver tabla
macro). Ese aumento no es una mejora de localización: `score_sigma` busca la
mejor geometría en **toda** la ventana, así que agregar 40-60 nt aguas arriba
aumenta la chance de encontrar una geometría espuria.

## Asignación

| genoma | configuración | top-1 | top-3 cualquiera | top-3 todos | n |
|---|---|---:|---:|---:|---:|
| K-12 | referencia -80/genomic | 57,0 % | 90,2 % | 88,9 % | 3189 |
| K-12 | -60/genomic | 57,6 % | 91,3 % | 90,2 % | 3189 |
| K-12 | recomendada -60/flank | 55,1 % | 89,2 % | 88,0 % | 3189 |
| S. subtilis | referencia -80/genomic | 44,8 % | 85,5 % | 85,5 % | 2226 |
| S. subtilis | -60/genomic | 45,6 % | 88,2 % | 88,2 % | 2226 |
| S. subtilis | recomendada -60/flank | 47,3 % | 87,6 % | 87,6 % | 2226 |

La recomendación **no** mejora la asignación de forma apreciable: a nivel macro
top-1 queda en +0,3 pp y top-3 en +0,55 pp, con un intercambio entre genomas
(K-12 -1,9 pp en top-1; S. subtilis +2,5 pp). Si el objetivo fuera top-1, la
alternativa `-60/genomic` es mejor para K-12 (57,6 %).

## Qué no funcionó

1. **Shuffles mononucleotídicos como fondo principal.** Bajan el soporte a
   37,9 % (K-12) y 31,5 % (S. subtilis) en `-60/+20`, y top-1 a 45,6 % / 39,8 %.
   El shuffle conserva composición pero destruye estructura dinucleotídica y la
   gramática del promotor; para un par de PWM de 6-11 pb el nulo es demasiado
   exigente. Sirve como control negativo, no como fondo recomendado.
2. **Ventanas largas (-100/-120) con flancos.** Suben el soporte a 54,7 %/55,8 %
   macro, pero top-1 cae a 46,2 %/44,3 % y top-3 a 85,7 %/83,8 %. No pasan el
   guardarraíl. El problema es que `score_sigma` no restringe la geometría a la
   distancia biológica del TSS.
3. **Ventanas largas con fondo genómico.** Al contrario, el soporte cae:
   K-12 45,9 %/42,1 % y S. subtilis 40,7 %/33,2 %. El nulo genómico de mayor
   largo y la búsqueda en más posiciones no compensan.
4. **-80/+20 con flancos.** Empata en soporte macro con -60/+20 (53,3 %) pero
   top-1 baja a 48,1 %; el guardarraíl lo rechaza. El óptimo de capa 1 no es
   monótono en el largo de la ventana.
5. **S. subtilis sigma32 y sigma38.** Siguen muy bajos con cualquier ventana y
   fondo (7,7 % y 9,9 % con la recomendación; baseline 7,7 % y 3,8 %). No se
   resuelve cambiando ventana/fondo; probablemente hay desajuste entre etiqueta
   y PWM disponible, o una arquitectura de promotor distinta. No se declara
   mejora.
6. **Asignación top-1.** Ningún fondo local la mejora de forma clara; el cambio
   recomendado es de capa 1, no de discriminación entre sigma.

## Límites

- **Capa 1 no es un llamador de TSS.** Mide si el sigma verdadero tiene mejor
  señal que un fondo; no prueba que el TSS sea correcto ni descarta TSS nuevos.
- **Geometría no restringida.** `score_sigma` maximiza el score en toda la
  ventana y para cualquier espaciador permitido; no exige que las cajas estén a
  una distancia biológica del TSS. Esto favorece ventanas largas y explica por
  qué -120 no se recomienda.
- **Fondos no intercambiables.** El fondo genómico es global y compartido por
  (sigma, largo); los fondos locales son condicionales a cada TSS. Los
  porcentajes entre fondos no son directamente comparables; la comparación
  válida es dentro de cada fondo y con la asignación como control.
- **Resolución estadística.** Con 100 decoys locales, el p mínimo es 1/101 ≈
  0,0099 y el percentil tiene granularidad de 0,01. No hay corrección por
  múltiples TSS/sigma. Los empates se tratan como en el evaluador: `pct` cuenta
  estrictamente menores y `pvalue` incluye el empate.
- **Gold heterogéneo.** K-12 es RegulonDB (3314 etiquetas/3189 TSS, 124 TSS
  multietiqueta) y S. subtilis es DBTBS (2226/2226, single-label). No son
  colecciones exhaustivas; S. subtilis tiene n=13 para sigma32 y n=11 para
  sigma54, así que sus porcentajes son inestables.
- **Flancos.** Se muestrean 100 subventanas de 400 pb por flanco; pueden
  solaparse entre sí y no se excluyen promotores internos. Para `up > 100`, el
  flanco aguas arriba `[-500,-100)` se solapa con el extremo distal de la
  ventana (20 pb para -120); no toca las cajas -35/-10, pero es un límite.
- **Recomendación discrecional.** El guardarraíl de 1 pp es una decisión de
  diseño, no un test estadístico. La combinación elegida maximiza capa 1 entre
  las que no degradan top-1/top-3 macro más de 1 pp.
- **No integrado al pipeline.** `50_eval_layers.py` no se modificó; usar la
  recomendación requeriría agregar el nulo de flancos por TSS al evaluador.

## Archivos

- `scripts/sigma/70_local_background.py` — barrido, fondos, métricas, regla de
  recomendación y `--from-json`.
- `tests/test_sigma_local.py` — pruebas de ventanas, flancos, shuffles,
  estadística, scorer vectorizado, recomendación y smoke test.
- `data/sigma/localization/comparison.json` — resultados completos, baseline,
  recomendación y límites.
- `data/sigma/localization/comparison.tsv` — tabla larga por
  (genoma, ventana, fondo, scope) con n, p<=0,05, pct>=0,95, z mediano y soporte.
- `data/sigma/localization/assignment.tsv` — top-1/top-3 y precisión/recall por
  sigma para cada combinación.

## Recomendación final

**Usar `-60/+20` con fondo de flancos locales (`-500/-100`, `+100/+500`,
100 subventanas por TSS).** Es la combinación con mayor soporte de capa 1 entre
las que mantienen la asignación macro dentro de 1 pp del baseline: 53,3 % vs
45,2 % macro (+8,1 pp), con top-1 51,2 % vs 50,9 % y top-3 88,4 % vs 87,9 %.
La ganancia es de sensibilidad de capa 1, no de asignación.

Si se prefiere no cambiar el fondo, la alternativa mínima es **`-60/+20` con
fondo genómico**: 48,5 % macro (+3,2 pp), top-1 51,6 % (+0,7 pp) y top-3 89,7 %
(+1,9 pp), con la misma definición de nulo que el evaluador actual.
