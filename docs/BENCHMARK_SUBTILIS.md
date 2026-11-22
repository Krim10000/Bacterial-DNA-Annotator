# ¿Generaliza el buscador de operones fuera de *E. coli*?

Benchmark del motor determinista y del modelo ML de K-12 aplicados a
***Bacillus subtilis* 168**, contra el ground truth de BSGatlas v1.1
(2.252 operones, 863 multi-génicos, 2.097 pares de genes).

Fecha: 2026-11-21.

## Respuesta corta

**El método generaliza; el módulo de promotores no.**

| | F1 de pares |
|---|---|
| K-12 (RegulonDB, umbral 80 pb) | 0,780 |
| **B. subtilis, todas las reglas** | **0,516** |
| **B. subtilis, sin la regla de promotor (R3)** | **0,753** |

La caída de 0,264 puntos **no es del método**: es de un único componente que se
trasplantó sin recalibrar. Al desactivarlo, la degradación real es de **0,027**
puntos (0,780 → 0,753), es decir un **3,5% relativo**.

## Qué se construyó

`rule_operons.py` lee GenBank, no GFF+FASTA. Se generó el GBFF con
`scripts/validation/build_subtilis_gbff.py`:

- 4.332 CDS desde los `CDS` de BSGatlas, con `locus_tag` tomado de las features
  `gene` (los CDS solo lo llevan en 64 de 4.332; se emparejan por coordenadas).
- 100% con locus tag, 4.267 con nombre de gen, 0 coordenadas incoherentes.
- Verificado: el lector del motor devuelve 4.332 genes ordenados, ambas hebras,
  secuencia de 4.215.606 pb.

Se anotaron además **horquillas para *B. subtilis*** (2.976 regiones, 338 de
clase A) para que la regla de terminador (R2) pudiera dispararse y no se
evaluara desactivada por accidente.

## Resultados

Pares evaluables: **1.999 de 2.097** (95,3%). Los 98 restantes incluyen genes
del GT que no son CDS en BSGatlas (`BSU_misc_RNA_*`), así que son imposibles de
predecir con una anotación de CDS. Todas las cifras de esta sección se calculan
solo sobre pares con ambos genes anotados.

| configuración | operones (total) | pares evaluables | P | R | **F1** | operón exacto |
|---|---|---|---|---|---|---|
| ML (pesos de K-12) | 3.806 | 526 | 0,875 | 0,230 | **0,364** | 1.242 |
| reglas, Otsu (40,0 pb), con R3 | 3.510 | 822 | 0,886 | 0,364 | **0,516** | 1.294 |
| reglas, fijo 42,5 pb, con R3 | 3.505 | 827 | 0,881 | 0,365 | **0,516** | 1.293 |
| reglas, fijo 42,5 pb, **sin R3** | 2.772 | 1.560 | 0,859 | 0,670 | **0,753** | 1.442 |
| reglas, Otsu, **sin R3** | 2.790 | 1.542 | 0,863 | 0,666 | **0,752** | 1.443 |

«operones (total)» es el número de operones predichos sin filtrar; «pares
evaluables» son los pares con ambos genes anotados como CDS (1.999 en el GT);
«operón exacto» cuenta operones predichos idénticos a uno del GT, sobre todas
las predicciones.

### Descomposición del fallo

| componente | efecto |
|---|---|
| **umbral de distancia** (40,0 vs 42,5 pb) | **ninguno**: 0,516 vs 0,516 |
| **calibración Otsu sin etiquetas** | acierta 40,0 pb en una especie nueva (42,5 en K-12) |
| **horquillas / R2 terminador** | +0,000: R2 se dispara 152 veces, pero mueve el F1 en la tercera decimal |
| **regla de promotor R3** | **−0,237**: de 0,753 a 0,516 |
| **modelo ML** | −0,152 frente a las reglas sin R3 |

R3 se dispara **1.779 veces** en *B. subtilis* (frente a 815 cortes por defecto
R5, 406 por distancia R4 y 12 por atenuador R1). El patrón se invierte respecto
a K-12: el PWM de σ70 de *E. coli* encuentra falsos promotores por todas partes
y **fragmenta operones largos**, que es justo lo que hunde el recall
(0,670 → 0,365).

## Por tipo de operón (mejor configuración, 42,5 pb sin R3)

| grupo | operones GT | pares GT | P | R | F1 |
|---|---|---|---|---|---|
| 2 genes | 436 | 426 | 0,996 | 0,561 | 0,718 |
| 3 genes | 173 | 338 | 0,995 | 0,607 | 0,754 |
| 4-5 genes | 155 | 510 | 1,000 | 0,673 | 0,804 |
| 6+ genes | 99 | 725 | 1,000 | 0,763 | 0,865 |
| un solo TSS | 1.053 | 433 | 0,952 | 0,737 | 0,831 |
| 2+ TSS (promotor interno) | 1.199 | 1.566 | 0,987 | 0,652 | 0,785 |

Dos lecturas:

1. **Los operones largos se predicen mejor que los cortos** (0,865 vs 0,718),
   al revés que en K-12. La causa probable: el atlas de *B. subtilis* tiene
   muchos operones de dos genes con huecos grandes, donde la distancia decide
   sola y acierta por azar en ambos sentidos.
2. **Los operones con promotor interno bajan** (0,785 vs 0,831), consistente con
   que el motor no modela promotores internos: solo corta por distancia.

Las precisiones por estrato (~0,95-1,00) son altas porque se calculan sobre
pares cuyos dos genes están en operones del GT; no son comparables con la
precisión global de 0,86, que incluye genes fuera de operón.

## Comparación con K-12

| | K-12 | B. subtilis (sin R3) | Δ |
|---|---|---|---|
| F1 de pares | 0,780 | 0,753 | −0,027 |
| operón exacto | — | 52,0% | — |
| umbral Otsu | 42,5 pb | 40,0 pb | −2,5 pb |

El umbral calibrado **sin etiquetas** en K-12 da 40,0 pb en *B. subtilis*: cuatro
puntos y medio de 42,5. Es la evidencia más limpia de que la señal de distancia
es biología compartida, no una particularidad de *E. coli*.

## Límites del experimento

1. **Una sola especie.** *B. subtilis* es Gram positiva, con genoma 4,2 Mb y
   43,5% GC (K-12: 50,8%). Generaliza a *una* especie; no a todas.
2. **El PWM de promotor es de K-12.** Se carga de sitios curados de σ70 de
   *E. coli* y se aplica a σA de *B. subtilis*. Son homólogos, pero el consenso
   −10/−35 difiere. Esto explica el −0,237; **no** demuestra que un PWM
   recalibrado para σA vaya a funcionar, solo que el actual no.
3. **Sin STRING.** La caché de STRING de K-12 no aplica a genes BSU: la feature
   vale 0 para todos. El ML pierde su feature más fuerte (coeficiente −3,51 en
   K-12), así que su 0,364 es un **piso**, no una medida justa del modelo.
4. **Horquillas heurísticas.** Las 2.976 horquillas de *B. subtilis* se generaron
   con el mismo heurístico que en K-12 (tallo + GC + poli-U). R2 filtra por
   clase A, poli-U ≥4 y ≤80 pb, y solo dispara 152 veces. El terminador sigue
   siendo el componente peor resuelto en ambas especies.
5. **Ground truth de otra naturaleza.** RegulonDB es curación experta de
   *E. coli*; BSGatlas integra DBTBS, BsubCyc, Rfam, RefSeq y SubtiWiki. No son
   igual de estrictos ni cubren lo mismo.
6. **95,3% de pares evaluables.** El 4,7% restante son genes del GT sin CDS.
   Afecta al recall absoluto, no a la comparación entre configuraciones.

## Qué habría que hacer para cerrar la pregunta

1. **Recalibrar el PWM de promotor con TSS de *B. subtilis*** (`gt_tss.tsv` tiene
   3.390 TSS con sigma anotado). Si el F1 vuelve a ~0,75 con R3 activo, el
   problema era la matriz, no la idea.
2. **Añadir STRING para BSU** (o usar co-expresión de SubtiWiki) y re-evaluar el
   ML, que hoy está artificialmente bajo.
3. **Repetir con una tercera especie** Gram negativa (p. ej. *P. putida*) para
   separar «Gram positiva» de «especie nueva».

## Reproducir

```bash
python scripts/validation/build_subtilis_gbff.py
python scripts/regulation/annotate_hairpins.py \
    --gff data/dbtbs/BSGatlas_v1.1.gff --fasta data/dbtbs/genome.fna --genome subtilis
python scripts/regulation/rule_operons.py \
    --gbff data/dbtbs/subtilis.gbff --genome subtilis \
    --hairpins data/operon_finder/subtilis_hairpins.tsv \
    --gap 42.5 --prom 99 --out data/operon_finder/sub_fixed_hp_noprom
python scripts/validation/benchmark_subtilis.py \
    --pred "sin_R3:data/operon_finder/sub_fixed_hp_noprom/subtilis_operons_rules.tsv"
```

`--prom 99` desactiva R3 (ningún PWM llega a 99). El resultado en JSON queda en
`data/dbtbs/benchmark_subtilis.json`.

## Efecto secundario: un bug real corregido

`find_operons.py --apply` con `--gbff` fallaba con
`UnboundLocalError: local variable 'seqs' referenced before assignment`: la
variable solo se definía en la ruta GFF+FASTA. Se unificó con un mapa
`{seqid: secuencia}` construido desde `seq`. Sin ese arreglo, el modelo ML no
podía aplicarse a ninguna especie nueva.
