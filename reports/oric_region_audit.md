# Auditoría de regiones oriC (punto ORCA vs intervalo DoriC)

**Generado:** 2026-11-26T01:26:36+00:00

**Alcance:** 400 genomas con estado oriC; solo datos locales. No se modificaron JSON de predicción, GFF3, `browser.py` ni `app.py`.

**Regla:** DoriC aporta el intervalo; ORCA solo el punto. Sin región DoriC se deja `span_source=ORCA (punto)` y `start=end=pos`: no se estima un ancho.

## Cobertura

| categoría | n | % |
|---|---:|---:|
| genomas | 400 | 100,0 |
| con región DoriC | 321 | 80.2 |
| solo punto ORCA | 79 | 19.8 |

### Lectura actual vs índice completo

| lectura | región | solo punto | sin lectura |
|---|---:|---:|---:|
| `load_oric` actual (fix-first) | 0 | 95 | 305 |
| índice de regiones | 321 | 79 | 0 |
| cambio | +321 | -16 | -305 |

**Hallazgo bloqueante:** `load_oric` exige `chosen.pos`. Solo los 95 JSON de `orca_predictions_fix` tienen `chosen`, y ninguno trae `doric_truth`; los 288 JSON de `orca_predictions` usan `orca_position` y `load_oric` los descarta. Resultado actual: 0 regiones, 95 puntos y 305 genomas sin ninguna lectura de oriC.

El índice recupera 321 regiones: 67 de genomas con fix que estaban en punto, 239 de predicciones sin fix y 15 de genomas sin predicción ORCA que sí tienen DoriC.

### Solo punto: por qué

| motivo | n |
|---|---:|
| punto de fix sin región DoriC | 28 |
| predicción ORCA sin DoriC | 49 |
| sin ORCA (posición de respaldo GC skew) | 2 |

De los 17 genomas sin predicción ORCA, 15 tienen región DoriC y 2 quedan solo punto. Sus posiciones de respaldo son redondas (600.000, 580.000, …) y confianza 0,3: no son un punto ORCA medido.

### Estado oriC

| estado | total | con región | solo punto |
|---|---:|---:|---:|
| validated | 271 | 271 | 0 |
| consensus | 77 | 0 | 77 |
| lineage_specific | 10 | 10 | 0 |
| anomalous | 42 | 40 | 2 |

## Validez de las 321 regiones DoriC

- No cruzan el origen: 257 (incluye 23 con `end == 0` resueltas por longitud).
- Cruzan el origen (`end > 0` y `end < start`): 64.
- `end == 0` es un centinela en 23 regiones; se reconstruye con la longitud declarada y no se asume cruce.
- Se normalizan a segmentos 1-based inclusivos; no se recorta ni se extiende el intervalo publicado.
- Diferencia entre longitud declarada y cobertura de los endpoints: {'0': 149, '1': 172}. ±1 pb corresponde a convenciones 0/1-based de las fuentes DoriC y se reporta, no se corrige.
- Filas DoriC cromosómicas locales: 434 para 321 genomas. 96 genomas tienen 2–4 filas.

### Discrepancia punto vs región

- El ancla que muestra el navegador cae dentro de la región en 70/321 genomas; mediana de distancia 7,067 pb.
- El punto ORCA **crudo** cae dentro de la región en 22/321 genomas; mediana de distancia 16,544 pb.

No se fuerza la coincidencia: si el punto cae fuera, se muestran ambos y se declara la distancia. El ancla de un genoma con fix puede ser DoriC/dnaA/dnaN, no ORCA; por eso se reportan las dos.

## Hallazgos de la auditoría

1. **`oriC_all.gff3` no contiene regiones DoriC.** Sus filas son ventanas de punto: 399 filas de 21 pb, 1 fila de 17 pb. No se usan como intervalo.
2. **`doric_ground_truth.json` colapsa filas.** Hay 797.145 claves de accesión; para estos 400 genomas, 321 tienen entrada. En los 96 genomas con varias filas cromosómicas, el JSON conserva la última; en 36 de ellos otra fila está más cerca del ancla (solo genomas con ancla ORCA/fix). El índice conserva la verdad canónica y lista las alternativas; no reemplaza la elección sin revisión. Además, en 5 genomas la GT apunta a una fila de plásmido: se descarta y se usa la fila cromosómica local (o solo punto).
3. **Regiones que cruzan el origen.** `load_oric` devuelve `start > end` en esos casos; un consumidor que asuma `start <= end` no dibuja ni extrae nada. El índice añade `segments`, `wrap` y `end_effective`.
4. **`end == 0` no es una coordenada.** En 23 filas DoriC es un centinela; se usa `oric_length` (igual a `len(oric_seq)`) para reconstruir el final sin inventar un ancho. No siempre implica cruzar el origen.
5. **El punto puede no ser ORCA.** En los 95 JSON de `orca_predictions_fix`, `chosen.cfg` es DoriC en 48 casos y dnaA/dnaN/sin marcador en 47; el punto mostrado es el ancla elegida, no necesariamente ORCA. El índice separa `point.pos`, `point.kind` y `raw_orca_pos`.
6. **Coordenadas DoriC heterogéneas.** Según la fuente, `oric_start` puede ser 0-based o 1-based; `oric_length`/`oric_seq` son la referencia más estable. El intervalo se conserva y la diferencia de ±1 pb se declara.

## Archivos

- `data/oric_regions/oriC_regions.json` — índice (entrada por genoma, región, punto, alternativas, discrepancia).
- `data/oric_regions/doric_ground_truth_400.json` — subconjunto compacto de DoriC (generado por `42_extract_doric_subset.py`).
- `reports/oric_region_inventory.tsv` — inventario de los 400 genomas.
- `webapp/oric_regions.py` — normalización, segmentos, extracción y selección pura.
- `tests/test_oric_regions.py` — pruebas puras de región/secuencia.

## Reproducir

```bash
python scripts/oric/42_extract_doric_subset.py
python scripts/oric/41_build_region_index.py
pytest tests/test_oric_regions.py -q
```

El builder no modifica `data/orca_predictions*` ni `doric_ground_truth.json`; solo lee y escribe los archivos nuevos (índice, subconjunto DoriC, inventario y reporte).
