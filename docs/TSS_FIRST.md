# El promotor manda: resultado central del módulo de operones

## Resumen

Usar los sitios de inicio de transcripción (TSS) como frontera **bate al motor
de cinco reglas** en las dos especies medidas, y con una sola regla:

| especie | motor de 5 reglas | regla del TSS |
|---|---:|---:|
| *E. coli* K-12 | 0,732 | **0,785** |
| *B. subtilis* 168 | 0,752 | **0,816** |

Y uniendo las dos fuentes de TSS disponibles, K-12 sube a **0,799** con la
precisión más alta de todas las configuraciones (0,834).

## Por qué

Un promotor marca dónde **empieza** un transcrito. Si no hay TSS entre dos genes,
se copian juntos. La distancia intergénica es solo un *proxy* de eso, y por eso
falla en los casos que vimos: corta con huecos de 88–128 pb donde no había
ningún promotor.

Medido en *B. subtilis*: de 1.779 cortes falsos con la matriz de *E. coli* a 23
con la propia, y de F1 0,516 a 0,750. El umbral de distancia, en cambio, apenas
cambia entre especies (42,5 pb en K-12, 40,0 en *B. subtilis*): esa señal sí es
biología compartida.

## Las dos fuentes de TSS, y por qué unirlas

| fuente | TSS | F1 operones | P | R | F1 TU |
|---|---:|---:|---:|---:|---:|
| curados (`PromoterSet`) | 3.964 | 0,785 | 0,743 | 0,833 | 0,760 |
| alto rendimiento (dRNA-seq, apoyo ≥2) | 13.405 | 0,756 | 0,721 | 0,793 | 0,736 |
| **unión** | **15.115** | **0,798** | **0,832** | 0,767 | **0,793** |

Los curados aportan precisión; los de alto rendimiento, cobertura. La unión gana
a los dos por separado. Se construye con
`scripts/validation/build_tss_union.py`.

## La verificación de circularidad (importante)

El 0,785 de K-12 usaba TSS de RegulonDB para predecir operones de RegulonDB: la
curaduría se confirmaba a sí misma. Con TSS **independientes** (dRNA-seq de 15
estudios) el resultado es **0,756** — sigue batiendo al motor (0,732), pero el
margen real es +0,024, no +0,053.

**Cifra honesta del sistema: 0,798 (unión), con la advertencia de que parte de
la mejora sobre 0,756 viene de datos de la misma casa que el ground truth.**

## El techo, medido

De 1.118 fronteras reales de operón en K-12:

| evidencia local | fronteras explicadas |
|---|---:|
| TSS curados + terminadores | 776 (69%) |
| + TSS de alto rendimiento | **945 (85%)** |
| sin evidencia de ningún tipo | 173 (15%) |

El 15% restante no es un fallo del método: son fronteras sin promotor ni
terminador anotado. Para recuperarlas hace falta más dRNA-seq o evidencia de
regulación, no mejor código.

## Cómo se usa

```bash
bash run_pipeline.sh tss     # une curados + alto rendimiento
bash run_pipeline.sh rules   # usa la unión si existe
```

`rule_operons.py --tss <fichero>` corta donde hay un TSS de la misma hebra
(R6_tss). Si se pasa `--tss`, la distancia se desactiva con `--gap 99999`:
mezclar ambas empeora (0,726 frente a 0,785), porque la distancia vuelve a cortar
donde el promotor ya decidió.

## Límites

- **Un TSS delimita una unidad de transcripción, no un operón.** Un operón puede
  tener varias TUs (21% en K-12, 38% en *B. subtilis*). Contra el GT de TUs la
  regla del TSS da 0,793; mezclar ambos niveles en una métrica es un error.
- **Sin TSS no hay regla.** Para un genoma nuevo sin datos de transcripción, la
  caída es a distancia, con la confianza declarada.
- **La familia Rho de terminadores no es detectable por secuencia** (AUC 0,67);
  los TTS anotados son escasos e incompletos.

## Validación en *B. subtilis*: la unión NO generaliza tal cual

Se repitió el experimento con las dos fuentes de *B. subtilis*: **BSGatlas**
(3.390 TSS) y **DBTBS** (865 TSS, extraídos de sus páginas de promotor con
posición absoluta y TSS estimado a ±2 pb).

| ground truth | TSS BSGatlas | TSS DBTBS | **unión** |
|---|---:|---:|---:|
| operones BSGatlas (2.252) | **0,799** | 0,659 | 0,795 |
| operones DBTBS (996) | **0,902** | 0,812 | 0,895 |

**La unión no mejora; empeora un poco.** Y la causa es medible: **43% de los TSS
de DBTBS ya están en BSGatlas** (±10 pb), porque BSGatlas integra DBTBS entre sus
fuentes. Los 497 restantes son nuevos pero ruidosos (TSS estimado desde la región
−42:+4, no medido).

### La regla correcta

La unión **no es universal**: funciona cuando las fuentes son **independientes**
(K-12: curados vs dRNA-seq de 15 estudios) y no cuando una contiene a la otra
(*B. subtilis*: BSGatlas ⊃ DBTBS).

Antes de unir hay que medir la independencia: si >40% de la fuente nueva ya está
presente, aporta ruido más que cobertura.

El hallazgo de K-12 se mantiene, pero con su condición explícita: **unir fuentes
independientes**, no "unir todo".
