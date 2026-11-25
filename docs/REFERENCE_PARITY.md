# Paridad con genomas de referencia

Este documento cubre dos informes complementarios sobre el panel
`data/reference_genomes.json`:

1. **paridad de conteos** (`scripts/reference/reference_parity.py`): compara
   cuántos CDS/tRNA/rRNA anota el sistema frente a la referencia;
2. **paridad por coordenadas** (`scripts/reference/coordinate_parity.py`):
   compara contig, hebra y coordenadas de cada feature, con emparejamiento 1:1
   y tolerancia configurable. Se documenta en la sección final
   [Paridad por coordenadas](#paridad-por-coordenadas).

## Objetivo (conteos)

Medir, sin evaluar función ni coordenadas, si el número de features que anota el
sistema se acerca al de la anotación curada de los genomas del panel
`data/reference_genomes.json`.

La comparación de esta primera parte es deliberadamente de **conteos por
clase**:

- CDS;
- tRNA;
- rRNA.

No mide productos, función, exactitud de coordenadas, identidad de secuencia ni
calidad de las predicciones. Un conteo parecido no implica que las features sean
las mismas. Para la comparación geométrica existe el informe de coordenadas
descrito más abajo.

## Entradas y selección

- Manifiesto: `data/reference_genomes.json`.
- Resumen del sistema: campo `system_summary` de cada registro, normalmente
  `results/<id>/<id>.summary.json`.
- Referencia: campo `annotation_path`, normalmente
  `data/reference_genomes/<id>/<id>_genomic.gff.gz`.

Se incluye un genoma cuando cumple las dos condiciones:

1. `status` es `annotated` o `available`;
2. declara `system_summary`.

Los registros con status elegible pero sin `system_summary` no se inventan ni se
completan con otros archivos: quedan listados como omitidos en el informe. Los
registros con otros status (`downloaded`, etc.) también quedan fuera.

## Reglas de conteo

- Se leen GFF3 planos o `.gz`; se ignoran comentarios y la sección `##FASTA`.
- Solo se usa la columna de tipo: `CDS`, `tRNA` y `rRNA`.
- Un CDS multiparte (varias filas con el mismo `ID` en el mismo contig) cuenta
  como **un** feature; las filas sin `ID` cuentan una a una.
- No se filtran pseudogenes por atributo (`pseudo=true`). Los registros
  `pseudogene` no son `CDS` y por tanto no entran en el conteo de CDS; si un CDS
  lleva `pseudo=true`, se cuenta como CDS porque la comparación es por tipo.
- No se cuentan `ncRNA`, `tmRNA`, `SRP_RNA`, `RNase_P_RNA`, `gene`, `exon` ni
  otros tipos, aunque estén en el GFF de referencia.
- Los conteos del sistema se toman del bloque `features` del resumen; no se
  recalculan desde el GFF3 del sistema.

## Cálculo y estados

Para cada genoma y clase:

- `delta = sistema - referencia`;
- `% sistema/referencia = 100 * sistema / referencia` (vacío si la referencia
  no tiene esa clase).

Estados:

| estado | significado |
|---|---|
| `ok` | CDS, tRNA y rRNA coinciden exactamente y la referencia tiene las tres clases. |
| `partial` | La comparación fue posible, pero al menos un delta es distinto de cero, o la referencia carece de alguna clase. |
| `not_available` | Falta o es ilegible el resumen o el GFF, el manifiesto no declara la ruta, o el GFF no contiene ninguna de las tres clases. |

El estado `partial` es el esperado cuando hay diferencias reales de conteo; no
se usa un umbral arbitrario de tolerancia. El agregado "micro" del informe suma
sistema y referencia por clase sobre los genomas comparables.

## Salidas y uso

```bash
.venv/bin/python scripts/reference/reference_parity.py
```

Genera:

- `reports/reference_parity.tsv`: una fila por genoma, con conteos, delta,
  porcentaje, estado, rutas y motivo. Es la salida para consumo automático.
- `reports/reference_parity.md`: informe legible con cobertura, tabla por
  genoma, motivos, entradas, agregado y límites.

Opciones útiles:

```bash
.venv/bin/python scripts/reference/reference_parity.py \
  --root . \
  --manifest data/reference_genomes.json \
  --tsv reports/reference_parity.tsv \
  --md reports/reference_parity.md
```

El script no usa red, GPU ni binarios externos. No modifica `results/`, los GFF
ni el manifiesto.

## Resultados actuales

Instantánea del **2026-11-25** (para cifras vigentes, ver
`reports/reference_parity.md`). En ese momento había 42 registros en el
manifiesto: 7 evaluados, 2 omitidos con status elegible pero sin
`system_summary` (`ecoli_k12`, `subtilis`) y 33 con otros status.

| genoma | estado | CDS sis | CDS ref | Δ CDS | % CDS | tRNA sis | tRNA ref | Δ tRNA | % tRNA | rRNA sis | rRNA ref | Δ rRNA | % rRNA |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mtb_h37rv` | partial | 3986 | 3906 | +80 | 102.05 | 45 | 45 | 0 | 100.00 | 2 | 3 | -1 | 66.67 |
| `cglutamicum` | partial | 3041 | 2999 | +42 | 101.40 | 60 | 59 | +1 | 101.69 | 12 | 18 | -6 | 66.67 |
| `pputida` | partial | 5507 | 5471 | +36 | 100.66 | 75 | 75 | 0 | 100.00 | 15 | 22 | -7 | 68.18 |
| `paeruginosa` | partial | 5715 | 5572 | +143 | 102.57 | 64 | 63 | +1 | 101.59 | 8 | 13 | -5 | 61.54 |
| `crescentus` | partial | 3738 | 3886 | -148 | 96.19 | 52 | 51 | +1 | 101.96 | 4 | 6 | -2 | 66.67 |
| `synechocystis` | partial | 3780 | 3698 | +82 | 102.22 | 41 | 41 | 0 | 100.00 | 2 | 6 | -4 | 33.33 |
| `saureus` | partial | 2675 | 2767 | -92 | 96.68 | 61 | 59 | +2 | 103.39 | 11 | 16 | -5 | 68.75 |

Agregado micro de esos 7 genomas:

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 28442 | 28299 | +143 | 100.51 |
| tRNA | 398 | 393 | +5 | 101.27 |
| rRNA | 54 | 84 | -30 | 64.29 |

Lectura rápida: los conteos de CDS y tRNA están cerca de la referencia en
términos agregados; rRNA queda sistemáticamente por debajo. Este informe no
atribuye esa brecha a una causa concreta, pero es compatible con que el pipeline
no llame todas las copias o tipos de rRNA (por ejemplo 23S/5S). Confirmarlo
requiere analizar las coordenadas y el caller de rRNA, fuera del alcance de esta
versión.

## Límites

- **Solo conteos**: no evalúa función, productos, coordenadas, ni identidad de
  secuencia.
- **El gold no es verdad absoluta**: hereda criterios y posibles errores de la
  anotación RefSeq/NCBI usada como referencia.
- **rRNA**: las diferencias pueden deberse a tipos no detectados (23S/5S),
  copias parciales o diferencias de criterio; no se diagnostican aquí.
- **ncRNA y otros tipos** quedan fuera, aunque estén en el GFF.
- **Cobertura del panel**: solo se comparan genomas ya anotados y con
  `system_summary`; los demás no se estiman ni se imputan.
- **Sin causalidad**: un delta no indica por sí solo un error del sistema ni de
  la referencia.
- **Sin uso para entrenamiento**: este panel es de evaluación; no se usa como
  gold de entrenamiento ni para declarar superioridad.

## Pruebas (conteos)

```bash
.venv/bin/python -m pytest tests/test_reference_parity.py -q
```

Los tests usan manifiesto, resumen y GFF sintéticos; cubren paridad exacta,
deltas, porcentajes, CDS multiparte, gzip, ausencia de resumen/GFF, resumen
ilegible, GFF sin clases objetivo y filtrado del manifiesto.

---

## Paridad por coordenadas

### Objetivo

Comparar, feature a feature, los **CDS, tRNA y rRNA** que el sistema escribe en
su GFF (`results/<id>/<id>.gff3`) contra los del GFF de referencia. A diferencia
del informe de conteos, aquí importan:

- el contig y la hebra;
- el inicio y el fin (GFF 1-based, inclusivo);
- cuántos features del sistema no tienen pareja (FP) y cuántos de la referencia
  no fueron encontrados (FN);
- la coincidencia exacta de bordes y la deriva media/mediana de inicio y fin.

El objetivo es responder si el sistema anota **las mismas features en los
mismos lugares**, no solo un número parecido.

### Entradas y selección

- Manifiesto: `data/reference_genomes.json`.
- GFF del sistema: si el registro declara `system_result`, se busca
  `<system_result>/<id>.gff3`; si no, `results/<id>/<id>.gff3`. También se
  aceptan `.gff3.gz`, `.gff` y `annotation.gff3(.gz)`.
- GFF de referencia: campo `annotation_path`, normalmente
  `data/reference_genomes/<id>/<id>_genomic.gff.gz`.
- Se incluyen los registros con `status` `annotated` o `available`. Si falta
  cualquiera de los dos GFF, el manifiesto no declara `annotation_path` o la
  referencia no contiene CDS/tRNA/rRNA, el genoma queda como `not_available`:
  no se imputan coordenadas ni se completa con otros archivos.

### Reglas de emparejamiento

1. Se leen GFF3 planos o `.gz`; se ignoran comentarios y la sección `##FASTA`.
2. Las filas multiparte con el mismo `ID`, contig, tipo y hebra se colapsan a
   un intervalo `[min(inicio), max(fin)]`. Las filas sin `ID` cuentan una a una.
3. Los GFF del sistema pueden traer la hebra negativa con `inicio > fin`; el
   intervalo se normaliza a `[menor, mayor]` antes de comparar.
4. El emparejamiento es **1:1** y exige **mismo contig, mismo tipo y misma
   hebra**. No se normalizan nombres de contig ni se cruzan hebras.
5. Un par es candidato si cumple al menos una condición:
   - **solape recíproco**: la fracción de solape respecto de cada feature es
     mayor o igual a `--min-reciprocal-overlap` (0.5 por defecto);
   - **bordes próximos**: `|Δinicio| <= tolerancia` y `|Δfin| <= tolerancia`
     (5 pb por defecto).
6. La asignación se resuelve con un greedy por calidad más caminos de aumento,
   de modo que sea de **cardinalidad máxima** y prefiera, en este orden, bordes
   exactos, menor deriva máxima y mayor solape. Así se conservan primero los
   pares exactos sin perder pares posibles.

La tolerancia de bordes rescata features pequeñas o desplazadas donde el solape
recíproco por sí solo no llegaría al 50 %. No relaja contig, tipo ni hebra.

### Métricas y estados

Por genoma y por clase (CDS, tRNA, rRNA), y también agregadas:

- **TP**: pares 1:1 asignados;
- **FP**: features del sistema sin pareja;
- **FN**: features de la referencia sin pareja;
- **precisión** = TP / (TP + FP); **recall** = TP / (TP + FN);
  **F1** = 2·P·R / (P + R);
- **coincidencia exacta**: pares con inicio y fin idénticos, en número y como
  porcentaje de los TP;
- **Δinicio** y **Δfin**: `sistema − referencia` en pb; media y mediana sobre
  los pares emparejados (signadas; positivas = el sistema extiende hacia la
  derecha).

Estados:

| estado | significado |
|---|---|
| `ok` | Todas las clases con datos tienen `FP=0`, `FN=0`, al menos un par 1:1 y deriva máxima de bordes dentro de la tolerancia. |
| `partial` | La comparación fue posible, pero hay FP/FN o algún par con deriva mayor que la tolerancia. |
| `not_available` | Falta o es ilegible el GFF del sistema o el de referencia, el manifiesto no declara `annotation_path`, o la referencia no contiene CDS/tRNA/rRNA. No se imputa nada. |

A nivel de clase, una clase ausente en ambas anotaciones se marca
`not_available` y no impide que el genoma alcance `ok` si las clases con datos
son perfectas. Un genoma con `ok` es paridad geométrica **dentro de la
tolerancia**, no necesariamente identidad exacta; el porcentaje de pares
exactos se informa aparte.

### Salidas y uso

```bash
.venv/bin/python scripts/reference/coordinate_parity.py
```

Genera:

- `reports/reference_coordinate_parity.tsv`: una fila por genoma y clase
  (`CDS`, `tRNA`, `rRNA`, `ALL`), con TP/FP/FN, tasas, bordes, tolerancia,
  rutas y motivo. Es la salida para consumo automático.
- `reports/reference_coordinate_parity.md`: informe legible con cobertura,
  tablas por genoma y clase, agregado micro, motivos, entradas y límites.

Opciones útiles:

```bash
.venv/bin/python scripts/reference/coordinate_parity.py \
  --root . \
  --manifest data/reference_genomes.json \
  --tolerance 5 \
  --min-reciprocal-overlap 0.5 \
  --tsv reports/reference_coordinate_parity.tsv \
  --md reports/reference_coordinate_parity.md
```

El script no usa red, GPU ni binarios externos. No modifica `results/`, los GFF
ni el manifiesto.

### Resultados actuales

Instantánea del **2026-11-25** (para cifras vigentes, ver
`reports/reference_coordinate_parity.md`). Con el manifiesto de ese momento:
25 registros elegibles, **24 comparables** (`ok` 0, `partial` 24), 1
`not_available` (`ecoli_k12`, sin GFF del sistema ni `annotation_path`) y 17
omitidos por status. Todos los genomas comparables tienen al menos un FP, un FN
o una deriva de bordes superior a 5 pb; por eso ningún genoma alcanza `ok`.
El panel se sigue anotando en paralelo, así que estas cifras son una instantánea
y cambian al regenerar el informe.

Agregado micro de los 24 genomas comparables (suma de TP/FP/FN):

| clase | TP | FP | FN | precisión | recall | F1 | exactos | % exactos | Δinicio medio | Δfin medio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CDS | 81213 | 2664 | 2164 | 0.9682 | 0.9740 | 0.9711 | 70649 | 86.99 % | +0.30 | -0.56 |
| tRNA | 1451 | 29 | 0 | 0.9804 | 1.0000 | 0.9901 | 1314 | 90.56 % | +0.07 | +0.04 |
| rRNA | 226 | 0 | 127 | 1.0000 | 0.6402 | 0.7806 | 0 | 0.00 % | +3.07 | -3.12 |
| ALL | 82890 | 2693 | 2291 | 0.9685 | 0.9731 | 0.9708 | 71963 | 86.82 % | +0.30 | -0.56 |

Lectura rápida:

- **CDS**: precisión 96.8 %, recall 97.4 %, F1 0.971; 87.0 % de los pares
  tiene bordes exactos y la mediana de deriva es 0 pb. Los FP/FN son
  compatibles con diferencias de calling de inicio/fin o de genes presentes en
  una sola anotación, no con un corrimiento global.
- **tRNA**: F1 0.990 y recall 1.000; el sistema no deja tRNAs de la referencia
  sin pareja. Hay 29 FP (tRNAs extra o duplicados respecto del gold).
- **rRNA**: precisión 1.000 pero recall 0.640; el sistema llama un subconjunto
  de las copias/tipos de rRNA de la referencia. Además, **ningún** par de rRNA
  tiene ambos bordes exactos: la mediana de Δinicio es +3 pb y la de Δfin
  -3 pb, compatible con un recorte sistemático de bordes. Esta es la principal
  fuente de estados `partial` y coincide con la brecha de conteos de rRNA del
  informe anterior.
- La deriva máxima reportada por genoma puede ser grande (cientos o miles de
  pb) en CDS largos con inicio/fin alternativos que igualmente superan el 50 %
  de solape recíproco. El informe la muestra explícitamente; un TP no implica
  bordes exactos.

### Límites (coordenadas)

- **Depende del GFF del sistema**: se compara lo que el pipeline escribió, no
  el resumen de conteos.
- **Solo tres clases**: `ncRNA`, `tmRNA`, `pseudogene`, `gene`, `exon` y otros
  tipos quedan fuera aunque estén en el GFF.
- **Emparejamiento geométrico**: un par 1:1 no demuestra que sea el mismo gen
  ortólogo ni que la secuencia traducida coincida.
- **Contig y hebra exactos**: diferencias de nomenclatura de contig o hebras no
  normalizadas producen FP/FN masivos; no se intenta resolverlas por homología
  ni por orden.
- **Multiparte**: se colapsa a `[min(inicio), max(fin)]`; no se evalúan
  intrones, fases ni exones individuales.
- **Solape recíproco**: el umbral del 50 % permite emparejar CDS largos con
  cambios grandes de inicio/fin. La deriva se reporta y `ok` exige que no
  supere la tolerancia.
- **El gold no es verdad absoluta**: hereda criterios y posibles errores de la
  anotación RefSeq/NCBI o de la fuente usada.
- **Cobertura del panel**: los genomas sin GFF del sistema o sin GFF de
  referencia quedan como `not_available`; no se estiman ni se imputan.
- **Sin causalidad ni uso para entrenamiento**: un FP/FN no indica por sí solo
  un error del sistema ni de la referencia.

### Pruebas (coordenadas)

```bash
.venv/bin/python -m pytest tests/test_coordinate_parity.py -q
```

Los tests usan manifiesto, GFF del sistema y GFF de referencia sintéticos;
cubren paridad exacta y métricas, tolerancia de bordes, umbral de solape
recíproco, TP/FP/FN y tasas, hebras/contigs distintos, CDS multiparte,
coordenadas invertidas en hebra negativa, gzip, ausencia de GFF, manifiesto sin
`annotation_path`, referencia sin clases objetivo, GFF vacío, `##FASTA` y
líneas inválidas, filtrado de status, agregado micro, escritura de TSV/Markdown
y no modificación de las entradas.
