# Recalibración del PWM de promotor para *B. subtilis*

**Pregunta**: ¿el F1 de 0,516 en *B. subtilis* con todas las reglas se debe a que
la matriz de promotor es de *E. coli*? ¿Vuelve a ~0,75 al recalibrarla?

**Respuesta corta**: sí, el problema **era la matriz** (0,516 → 0,750), pero la
regla de promotor **no aporta nada** ni siquiera con la matriz correcta
(0,750 frente a 0,752 sin la regla). La matriz estaba mal; la regla, además, es
redundante.

---

## 1 · Qué se construyó

`scripts/validation/build_promoter_pwm_subtilis.py` → `data/dbtbs/promoter_pwms_subtilis.tsv`

- **1.949 TSS de σA** (de 3.390 totales; σA es el factor housekeeping, análogo
  de σ70). Todos con coordenada válida.
- Ventana de 45 pb aguas arriba, corregida por hebra.
- **No hay caja −35**: al buscar hasta 90 pb aguas arriba, el mejor 6-mero es
  poli-T (`TTTTTT`, 0,57 bits), no `TTGACA`. Los promotores σA de este set son
  en su mayoría *extended −10*.
- Matriz resultante: **una sola caja de 12 pb**, `TATAATTAAAAT`, en −12..−1
  (caja −10 + región de inicio), 1,71 bits.

La señal existe pero es débil: `TATAAT` exacto aparece en **30 de 1.949** TSS en
la posición −12, frente a **1** en el control desplazado 200 pb (30× de
enriquecimiento). Es un motivo real, no un artefacto.

### Comparación con la matriz de *E. coli*

| | *E. coli* (σ70) | *B. subtilis* (σA) |
|---|---|---|
| cajas | −35 y −10, 6 pb cada una | una sola, 12 pb |
| consenso −10 | `TATAAT` | `TATAAT` (−12..−7) |
| consenso −35 | `TTGACA` | **no hay** |
| información | 4,83 + 3,67 bits | 1,71 bits |

Correlación de la parte −10 (6 pb alineados): **r = 0,577**. La caja −10 se
parece; la diferencia estructural es la **ausencia de −35**.

## 2 · Cambios de código

| fichero | cambio |
|---|---|
| `scripts/regulation/find_operons.py` | `promoter_scan` acepta una matriz de **una sola caja** (si falta −35 o −10). El camino de dos cajas queda idéntico. |
| `scripts/regulation/rule_operons.py` | nueva opción `--pwm <ruta>`; `--out` y `--gbff` se resuelven a absoluto; la comparación contra RegulonDB solo se hace si `--genome ecoli_k12`. |

Ninguno de los tres cambios altera el resultado de K-12 con su matriz: verificado
en §4.

## 3 · Resultado en *B. subtilis*

GT: 2.252 operones de BSGatlas, 1.999 pares evaluables.

| configuración | disparos R3 | P | R | **F1** |
|---|---|---|---|---|
| sin R3 | 0 | 0,863 | 0,666 | **0,752** |
| matriz *E. coli*, umbral 8 (el de antes) | 1.779 | 0,886 | 0,364 | **0,516** |
| matriz nueva, umbral 3 | 1.286 | 0,878 | 0,468 | 0,610 |
| matriz nueva, umbral 4 | 297 | 0,865 | 0,636 | 0,733 |
| matriz nueva, umbral 4,2 | 203 | 0,863 | 0,647 | 0,740 |
| matriz nueva, umbral 4,4 | 103 | 0,863 | 0,657 | 0,746 |
| matriz nueva, umbral 4,6 | 23 | 0,863 | 0,663 | **0,750** |
| matriz nueva, umbral 4,8 | 0 | 0,863 | 0,666 | 0,752 |

Lectura: con la matriz nueva y el umbral bien puesto, R3 deja de destruir el F1
(0,516 → 0,750). Pero **en el mejor punto con la regla activa (4,6) el F1 sigue
por debajo de no usarla** (0,750 vs 0,752), y la precisión no mejora: los 23
cortes que añade no traen ni un par correcto nuevo.

## 4 · Control en K-12: la regla tampoco ayuda donde se calibró

| configuración | disparos R3 | P | R | **F1** |
|---|---|---|---|---|
| sin R3 | 0 | 0,863 | 0,670 | **0,755** |
| umbral 8 (por defecto) | 873 | 0,873 | 0,448 | 0,592 |
| umbral 10 | 313 | 0,869 | 0,617 | 0,722 |
| umbral 12 | 59 | 0,866 | 0,663 | 0,751 |
| umbral 16 | 0 | 0,863 | 0,670 | 0,755 |

El patrón es el mismo que en *B. subtilis*: R3 solo alcanza la paridad cuando
**deja de disparar**. El 0,755 de referencia se obtiene con R3 efectivamente
apagada. El valor por defecto (`PROM_MIN = 8.0`) **penaliza K-12** en 0,16 de F1
si alguien lo deja tal cual.

## 5 · Interpretación

El diagnóstico era correcto a medias:

1. **La matriz sí era el problema.** Con la de *E. coli*, 1.779 cortes; con la
   propia, 23. El 0,516 era un artefacto de usar un motivo ajeno.
2. **La regla no sirve como frontera.** Ni con la matriz correcta: un promotor
   es una señal de *dónde empieza la transcripción*, y muchas veces eso ocurre
   **dentro** de un operón (promotores internos, TSS secundarios). Cortar ahí
   fragmenta operones verdaderos. En *B. subtilis*, 1.199 de los 2.252 operones
   tienen 2+ TSS: la mitad del problema.
3. **La distancia hace todo el trabajo.** El Otsu calibrado sin etiquetas da
   40,0 pb en *B. subtilis* frente a 42,5 en K-12, y con la regla de distancia
   sola se llega a 0,752.

## 6 · Límites del experimento

- **Techo del 95,3%**: 309 genes del GT no son CDS en BSGatlas; las cifras son
  sobre 1.999 de 2.097 pares.
- **Umbral elegido sobre el mismo GT** que se evalúa. Es una fuga leve: por eso
  importa que incluso el mejor punto (4,6) no supere a no usar la regla. Aun
  así, un umbral elegido en otra especie podría ser distinto.
- **Una sola especie nueva**: no permite separar "σA de *B. subtilis* es débil"
  de "los TSS de BSGatlas son heterogéneos".
- **σA agregado**: 1.949 TSS pueden mezclar promotores primarios y secundarios;
  el 1,5% con `TATAAT` exacto sugiere eso.

## 7 · Recomendación

1. **Dejar R3 apagada por defecto** (o exigir `--prom` ≥ 12) hasta que haya una
   señal que distinga promotores de frontera de promotores internos.
2. **Usar los promotores como anotación, no como regla**: ya están en la GUI
   (componentes del operón) y en el navegador. Ahí aportan contexto sin cortar
   operones.
3. Si se quiere recuperar R3, el siguiente experimento es **condicionar el
   promotor a que no haya un TSS interno conocido** (o a que la distancia al gen
   previo lo haga plausible como inicio de operón), no volver a calibrar la
   matriz.

## Reproducir

```bash
python scripts/validation/build_promoter_pwm_subtilis.py
python scripts/regulation/rule_operons.py \
    --gbff data/dbtbs/subtilis.gbff --genome subtilis \
    --hairpins data/operon_finder/subtilis_hairpins.tsv \
    --pwm data/dbtbs/promoter_pwms_subtilis.tsv --prom 4.6 \
    --out data/operon_finder/sub_new_p4.6
python scripts/validation/benchmark_subtilis.py \
    --pred "nuevo:data/operon_finder/sub_new_p4.6/subtilis_operons_rules.tsv"
```
