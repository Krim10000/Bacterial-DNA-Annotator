# Buscador de operones para genomas no anotados

Método **independiente** de la capa L2 (`scripts/regulation/`): no lee sus
salidas, predice todo desde la secuencia.

## Por qué existe

En un genoma sin anotar no hay componentes que leer. Solo *E. coli* K-12 tiene
componentes curados (RegulonDB), y se usa para dos cosas:

1. **Calibrar los detectores** de componentes (terminador, promotor).
2. **Entrenar** el clasificador de fronteras.

La aplicación es al genoma 401.

## Componentes de novo

| componente | método | validación |
|---|---|---|
| Terminador intrínseco | horquilla GC (≤1 error, bucle 2–10) + tracto poli-T, sin exigir adyacencia | **AUC 0.936** contra los 375 terminadores curados |
| Promotor | PWM σ70 −10/−35 con espaciador 14–21 | consensos canónicos (`TATAAT`/`TTGACA`) |
| RBS | motivo Shine-Dalgarno en ±20 pb del inicio | — |

Lección del desarrollo: la primera versión exigía coincidencia exacta del tallo
**y** poli-T inmediatamente detrás. Su AUC era 0.445 (peor que azar). Separar las
dos señales y tolerar un error en el tallo la llevó a 0.936.

## Fronteras

Regresión logística sobre `[log distancia, solapamiento, terminador, promotor,
RBS, misma EC]`, entrenada en K-12 (2.897 fronteras etiquetadas: 1.118 cortes).
Se segmenta el cromosoma cortando donde P(frontera) ≥ umbral.

Coeficientes en K-12: `log_dist +1.31`, `same_ec −2.41`, `overlap +3.17`,
`promoter +0.11`, `terminator +0.04`.

## Resultado en K-12 (contra RegulonDB, 827 operones)

| versión | P | R | F1 | exact-set |
|---|---|---|---|---|
| L2 (`regulon-builder`) | 0.844 | 0.761 | 0.808 | 0.629 |
| `find_operons` solo componentes | 0.819 | 0.751 | 0.784 | 0.550 |
| **+ score funcional STRING** | **0.867** | **0.769** | **0.815** | **0.617** |

El salto viene de un único feature. STRING entra con el coeficiente más fuerte
del modelo (**−3.51**; score alto → no hay corte) y absorbe parte del peso de
`log_dist` (1.18 → 0.82) y de `same_ec` (−3.99 → −1.25). Coincide con la
literatura: Taboada et al. 2010 llega a 94,6% en E. coli con distancia + score
funcional, sin terminadores ni TSS.

En los casos **documentados (S)**: recall de pares 0.524 → **0.754**, exacto
43% → 55%.

### Punto de operación de alta precisión

El barrido muestra el compromiso: con umbral 0.2, **precisión 0.957** cubriendo
49,7% de los pares verdaderos. Ese es el material de la capa 1 de la cascada
(`cascade.py`, pendiente de la tabla por capa).

Dos métodos construidos por caminos distintos convergen en ~0,78–0,81. Eso es
información: sugiere que el techo de la evidencia derivable de la secuencia está
ahí, y que para subir hace falta otro tipo de dato (TSS/TTS experimentales,
expresión), no otro algoritmo.

Por tamaño, este método rinde peor que L2 en 4–6 y 7+ (0.486 y 0.429 vs 0.577 y
0.442): es más conservador y sub-une. El umbral elegido por barrido fue 0.7
(F1 0.785); con 0.5 daba 0.764.

## Uso

```bash
# validar los componentes
python scripts/regulation/find_operons.py --validate

# entrenar en K-12 (una vez)
python scripts/regulation/find_operons.py \
    --gbff data/test_ecoli.gbff --genome ecoli_k12 --train --sweep

# aplicar a un genoma no anotado
python scripts/regulation/find_operons.py \
    --gff data/annotations_ab_initio/<GCF>/annotation.gff3 \
    --fasta data/raw/<GCF>_genomic.fna.gz --genome <GCF> --apply
```

Salidas en `data/operon_finder/`: `boundary_model.pkl` y
`<genoma>_operons.tsv` (formato compatible con
`scripts/validation/benchmark_operons.py`).

## Límites

- Entrenado en **un** genoma. El umbral y los coeficientes son de *E. coli*; en
  otros linajes la distancia intergénica y la composición cambian.
- No usa sintenía entre genomas (existe, pero en K-12 cubre ~10%).
- El 7% de CDS de K-12 sin entrada en RegulonDB no se usa para entrenar.
- `overlap` con coeficiente positivo es contraintuitivo (el solapamiento suele
  indicar acoplamiento traduccional); conviene revisarlo con más genomas.

## ¿Y ViennaRNA para el fin del operón? (probado)

Sí es el estándar: ARNold usa ΔG de RNAfold con corte −10 kcal/mol más tracto
poli-U, y hay toda una familia (TransTermHP, RNIE, WebGeSter, BATTER-TPE,
RhoTermPredict para los dependientes de Rho).

Se probó aquí (`compare_structure.py`, ViennaRNA 2.7.2) sobre **intergénicas
reales** de K-12, no sobre secuencias ya sabidas terminadoras:

### Corrección importante: esos AUC eran artefacto de longitud

La primera medición (proxy 0.852, ViennaRNA 0.829) usaba negativos
**emparejados por longitud** con los terminadores. La puntuación crece con el
tamaño de la ventana, así que ese esquema la infla.

Medición corregida (`robustness.py`), con negativos que son **todas las
intergénicas reales**:

| esquema | proxy | ViennaRNA |
|---|---|---|
| bootstrap sobre intergénicas (IC95%) | **0.649** [0.608, 0.687] | — |
| ventana fija 60 nt, máximo deslizante | 0.644 | **0.426** |
| negativos barajados (mononucleótido) | 0.923 | — |

ViennaRNA queda **peor que el azar** con ventanas fijas: el ΔG de una
intergénica real no es menos estable que el de un terminador curado. Y el corte
≥6 marca el **100% de las intergénicas** en otros genomas: el umbral no porta.

Conclusión corregida: **ninguna de las dos señales discrimina bien a escala
genómica**. La afirmación anterior ("el ΔG no supera al proxy") era cierta, pero
por una razón más fuerte: ninguna de las dos sirve como feature potente. Eso es
coherente con el coeficiente del terminador en el modelo (~0.05, despreciable) y
con la cobertura experimental del 16% en el lado 3'.

Conclusión: la estructura secundaria no cierra el problema del fin de operón.
Los propios trabajos de la familia lo dicen — no todos los terminadores
predichos terminan la transcripción in vivo, y los dependientes de Rho no tienen
horquilla. Para fronteras reales hace falta **mapeo experimental de extremos 3'**
(Term-seq, dRNA-seq) o co-expresión por RNA-seq.

## Extremos experimentales (TTS/TSS): descargados y probados

`scripts/validation/fetch_ht_termini.py` trae de RegulonDB los datasets de
alto rendimiento:

| dataset | filas | únicos |
|---|---|---|
| TUS (unidades de transcripción) | 9.829 | 8.358 |
| TTS (extremos 3') | 3.099 | 2.281 |
| TSS (extremos 5') | 60.248 | ~28.000 |

Cobertura real de las fronteras de K-12 (±150 pb):

| evidencia | fronteras cubiertas |
|---|---|
| TSS cerca del inicio del gen siguiente | **61%** |
| TTS cerca del fin del gen previo | **16%** |
| ambos | 13% |

El lado 5' está razonablemente mapeado; el 3' casi no. Eso explica por qué el
"fin del operón" es el problema duro.

**Como feature, no mejora**: F1 0.807 con TTS/TSS frente a 0.815 solo con
STRING. Los coeficientes son positivos (`tss_near` +0.47, `tts_near` +0.32) pero
un TSS cercano **no implica frontera de operón** — los promotores internos
existen. Queda como opción (`--with-termini`), apagada por defecto.

## Marcador final del método

| configuración | F1 |
|---|---|
| distancia + componentes | 0.784 |
| + posición del promotor | 0.783 |
| **+ STRING** | **0.815** |
| + TSS/TTS experimentales | 0.807 |

Solo el score funcional de STRING movió la aguja, y coincide con la literatura
(Taboada et al. 2010: 94,6% con distancia + funcional). Las features
estructurales y experimentales no superan el techo de ~0,81.

## Desagregación del terminador por clase estructural

`terminator_classes.py` separa los 375 terminadores curados en clases por
componentes medibles (tallo, GC del tallo, bucle, poli-U) y los contrasta con
600 ventanas de 48 nt de intergénicas reales.

| clase | regla | % de curados | enriquecimiento | lectura |
|---|---|---|---|---|
| A fuerte | stem≥6, GC≥0.60, poliU≥4 | 30% | **3.9×** | se identifica bien |
| B media | stem≥5, poliU≥3 | 52% | **1.0×** | indistinguible del fondo |
| C débil | resto | 18% | 0.4× | no se identifica |

Datos clave:

- **El loop no discrimina**: la clase B tiene tallos más largos (11 pb) que la A
  (10 pb) y aun así no se distingue del fondo. Lo que separa es **GC del tallo**
  (0.73 en A vs 0.50 en B) y **longitud del poli-U** (6 vs 5; en C, 2).
- **10% del catálogo es dependiente de Rho** (37 de 375): sin horquilla, invisible
  a cualquier método estructural.
- Aceptando solo clase A: 112 verdaderos de 375 y 46 falsos de 600 (~0.71 de
  precisión, 30% de recall).
- **Ambigüedad conformacional** (`pf_fold`): los terminadores curados ganan menos
  energía en el ensemble que las intergénicas (0.61 vs 0.83 kcal/mol), es decir
  pliegan de forma más determinista. La diferencia existe pero es pequeña.

## ¿Toda horquilla termina la transcripción? (no)

Contra los 827 operones multigénicos curados, posición de los 375 terminadores:

| situación | n |
|---|---|
| en el extremo 3' de su operón (≤80 pb) | **9** |
| **dentro del operón, lejos del final** | **70** |
| fuera de todo operón multigénico | 296 |
| dependientes de Rho | 37, de los cuales **8 con horquilla clase A** |

Dos hechos duros:

1. **70 terminadores curados están dentro de un operón**, no en su final: la
   horquilla termina el transcrito y el operón sigue. Es la firma empírica de la
   terminación condicional (atenuación) — la misma secuencia pliega a terminador
   o a antiterminador según el ribosoma y el metabolito.
2. **8 terminadores dependientes de Rho tienen horquilla clase A**: la presencia
   de una horquilla fuerte no dice qué mecanismo actúa.

Más: las horquillas clase A aparecen en ~4% de las ventanas intergénicas, muy por
encima de los 375 terminadores curados. Es decir, la horquilla **no es necesaria**
(10% son Rho) ni **suficiente** (mayoría de horquillas no terminan, y hay
terminación interna).

## Motor determinista por reglas (sin modelo)

`rule_operons.py` reemplaza el clasificador logístico por una tabla de decisión
con motivo auditable por frontera (salida `*_boundary_rules.tsv`):

| regla | cuándo | decisión |
|---|---|---|
| R1 atenuador | el gen previo es péptido líder (≤35 aa) | UNIR |
| R1b atenuador funcional | horquilla clase A y STRING ≥0.70 | UNIR |
| R2 terminador | horquilla A, poli-U ≥4, ≤80 pb del gen previo | CORTAR |
| R3 promotor | σ70 con score alto en posición plausible | CORTAR |
| R4 distancia | hueco ≥ umbral | CORTAR |
| R5 defecto | nada de lo anterior | UNIR |

**R1 y R1b son la traducción directa de "una horquilla dentro de un operón es un
atenuador"**: 7 operones curados empiezan con péptido líder (`pheL`, `thrL`,
`hisL`, `ilvL`, `tnaC`, `rhoL`, `tyrT`) y 70 terminadores están dentro de su
operón. La regla no corta ahí.

Rendimiento (K-12, umbral de hueco 80 pb, promotor ≥14):

| método | F1 de pares | tiempo |
|---|---|---|
| modelo logístico (con STRING) | 0.815 | ~10 min |
| reglas deterministas | **0.780** | **8 s** |
| L2 del otro agente | 0.808 | — |

Recall por tipo de operón (reglas):

| tipo | ops | pares | recall |
|---|---|---|---|
| multigénico sigma alternativo | 40 | 220 | **0.927** |
| multigénico simple | 191 | 878 | **0.861** |
| multigénico sin evidencia | 234 | 873 | 0.806 |
| multigénico con promotores internos | 362 | 2.728 | **0.719** |

El tipo con promotores internos concentra el 58% de los pares y es el cuello de
botella, igual que en el modelo. El determinista cuesta ~0.035 de F1 pero gana
auditabilidad total y 75× de velocidad.

Cautela: el umbral de 80 pb salió de un barrido sobre K-12; es calibración, no
un valor universal. En otro genoma debería recalibrarse con su propia
distribución de huecos.

## Lote completo: 311 genomas

`annotate_hairpins.py` corrió sobre los 311 genomas con FASTA y GFF: **808.830
horquillas**. El resumen está en `data/operon_finder/summary.json`.

Primer resultado del lote: la clase A dependía del **GC del genoma** (r=0,48
sobre 311 genomas). Con el umbral fijo de *E. coli* (`GC tallo ≥ 0,60`), un
genoma de 27% GC daba 1% de clase A; con umbral relativo (`GC genoma + 0,10`)
da ~12%. El total apenas cambia (17,3% → 15,9%), pero se reparte entre especies.

Los TSV se re-etiquetaron: `class` es la normalizada y `class_fixed` conserva la
calibrada en K-12. Para recalcular: `summarize_hairpins.py --rewrite`.

Nota de proceso: **no editar `annotate_all.sh` mientras corre** — bash lee los
scripts incrementalmente y termina con errores de sintaxis en la cola.

## Método vigente (2026-11-21): el promotor decide

Tras medir todas las piezas, el camino por defecto es:

| pieza | estado | aporte medido |
|---|---|---|
| **TSS (unión curados + alto rendimiento)** | **principal** | 0,785 → **0,798** F1 de pares |
| distancia (Otsu) | respaldo cuando no hay TSS | 0,732 sola |
| R2 terminador | apagada con TSS | 0 o negativo |
| R3 promotor (PWM) | apagada | 0 o negativo |
| segunda pasada EC/sintenia | **diagnóstico**, no defecto | +0,003 |
| sintenia densa (RBH) | diagnóstico | 0,000 |

```bash
bash run_pipeline.sh tss     # construye la unión de TSS
bash run_pipeline.sh rules   # usa la unión si existe (--tss ... --gap 99999)
```

**No mezclar TSS y distancia**: cortar por ambas empeora (0,726 frente a 0,798),
porque la distancia vuelve a cortar donde el promotor ya decidió.

Herramientas de diagnóstico que no están en el camino por defecto:
`refine_with_synteny.py`, `boundary_prior.py`, `boundary_classifier.py`.

El detalle y las advertencias (circularidad, techo del 85%) están en
[`docs/TSS_FIRST.md`](../../docs/TSS_FIRST.md).
