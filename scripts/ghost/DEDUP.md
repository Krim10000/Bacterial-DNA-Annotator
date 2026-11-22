# Experimento: umbral de deduplicación entre genomas

**Pregunta:** ¿a qué identidad agrupamos las cadenas candidatas de los 400
genomas para no procesar duplicados, sin fusionar genes distintos?

## Método

1. Se tomaron **20.000 proteínas ya anotadas** de los 400 genomas, con hit
   Swiss-Prot de **longitud completa** (alineamiento ≥70% de la proteína
   consulta y del sujeto). Esto evita que un fragmento herede el nombre de una
   proteína grande.
2. Verdad de referencia: mismo gen = mismo `GN=` del hit Swiss-Prot; si no hay,
   mismo acceso.
3. Se agruparon con `diamond cluster` a 80, 90, 95, 97 y 100% de identidad.
4. Para el recall se corrió un **all-vs-all** de DIAMOND (46.719 pares con
   identidad y cobertura ≥80%) y se definió "debería agruparse" como: mismo gen
   **y** identidad real ≥ umbral.
5. Métricas: precisión (de los pares agrupados, ¿son el mismo gen?) y recall
   (de los pares que deberían agruparse, ¿se agruparon?).

## Resultados

| umbral | pares agrupados | precisión | recall | F1 |
|---|---|---|---|---|
| 80% | 299 | 0,910 | 0,941 | 0,925 |
| 90% | 100 | 0,890 | 0,986 | 0,936 |
| **95%** | **56** | **0,964** | **1,000** | **0,982** |
| 97% | 44 | 1,000 | 1,000 | 1,000 |
| 100% | 22 | 1,000 | 1,000 | 1,000 |

## Decisión

**Umbral = 95%.**

- Es el menor umbral con precisión ≥0,95 (no fusiona genes distintos de forma
  apreciable) y recall 1,0.
- A 90% la precisión cae a 0,89: empieza a mezclar genes distintos.
- 97% es la alternativa conservadora (precisión 1,0, F1 1,0) a cambio de un
  21% menos de compresión.

## Salvedades

- La muestra son proteínas anotadas y está dominada por ortólogos lejanos; por
  eso la compresión observada es pequeña (0,26% de clusters). En el conjunto de
  ORFs candidatos, donde abundan las copias casi idénticas entre cepas, la
  compresión será mayor.
- La estimación de precisión se apoya en 56 pares; en la corrida de 400 genomas
  conviene recalcularla sobre una muestra mayor.
- Para péptidos cortos (60–150 aa) la identidad es más ruidosa: si aparecen
  muchas fusiones dudosas, subir a 97%.

## Artefactos

- `out/dedup/dedup_thresholds.tsv` — métricas por umbral (clusters, pureza).
- `out/dedup/dedup_pairs.tsv` — precisión/recall reales por umbral.
- `out/dedup/dedup_pairs_choice.json` — umbral elegido y criterio.
- `out/dedup/sample.faa`, `clusters_*.tsv`, `allvall.tsv` — datos del experimento.
