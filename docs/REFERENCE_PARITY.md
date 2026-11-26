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

Instantánea del **2026-11-25** (cifras vigentes en
`reports/reference_parity.md`). El manifiesto tiene 42 registros: 40 evaluados
(todos `annotated`/`available` con `system_summary`) y 2 omitidos por no
declarar `system_summary` (`ecoli_k12`, `subtilis`).

Agregado micro de los 40 genomas:

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 143681 | 142980 | +701 | 100.49 |
| tRNA | 2423 | 2334 | +89 | 103.81 |
| rRNA | 514 | 515 | -1 | 99.81 |

El rRNA pasa de **335/515 (65.05 %)** antes de la corrección a **514/515
(99.81 %)** después. El único delta restante es `paeruginosa` (12/13: falta una
copia de 5S que el modelo Rfam local no recupera; no se imputa). Los deltas de
CDS y tRNA son los de la anotación general y quedan fuera de esta corrección.

## Diagnóstico de la brecha rRNA (P0)

Antes de la corrección, el agregado de coordenadas daba rRNA con precisión
1.000, recall 0.6514, F1 0.7889 y **0 pares con bordes exactos**. La brecha no
venía del emparejamiento ni de las coordenadas:

- el parser GFF ya es 1-based inclusivo y normaliza la hebra negativa;
- los `seqid` y las hebras coincidían (precisión 1.000, sin FP);
- no había features fusionados: 16S/23S/5S son filas separadas en la
  referencia y en el sistema.

La causa era la **llamada**: el pipeline usaba solo `barrnap` 0.9.

1. **Falta el modelo 23S en la base de barrnap.** En esta máquina
   `/usr/share/barrnap/db/bac.hmm` contiene solo `16S_rRNA` (RF00177) y
   `5S_rRNA` (RF00001): no contiene `LSU_rRNA_bacteria`/RF02541. Por eso el
   sistema escribía **0 copias de 23S** en todos los genomas.
2. **El E-value por defecto descarta copias de 5S.** `barrnap` usa `-E 1e-6`;
   hay copias de 5S con E=0.023 (p. ej. *H. pylori*) que quedan fuera.
3. La referencia RefSeq/NCBI anota los tres tipos con Infernal/Rfam. El repo ya
   incluye `data/rfam_models/bac_rRNA.cm` con exactamente RF00177 (16S),
   RF02541 (23S) y RF00001 (5S).

Evidencia por genoma (conteo 16S/23S/5S):

| genoma | antes | después (recall) | referencia |
|---|---|---|---|
| `cglutamicum` | 6/0/6 = 12 | 6/6/6 = 18 | 6/6/6 = 18 |
| `hpylori` | 2/0/0 = 2 | 2/2/2 = 6 | 2/2/2 = 6 |
| `synechocystis` | 2/0/0 = 2 | 2/2/2 = 6 | 2/2/2 = 6 |
| `dradiodurans` | 3/0/0 = 3 | 3/3/3 = 9 | 3/3/3 = 9 |
| `subtilis` | 10/0/10 = 20 | 10/10/10 = 30 | 10/10/10 = 30 |

La corrección está en `webapp/rrna_calling.py` (parsers, deduplicación y
subconjuntos de modelos), `webapp/annotator.py` (`run_rrna` prefiere
`cmsearch` y cae a `barrnap --evalue 10`) y
`scripts/annotation/00_annotate_single.py` (preferencia Infernal y raíz de
datos correcta). La reanotación del panel se hace con
`scripts/reference/reannotate_rrna.py`, que **no usa el GFF de referencia como
entrada**: llama rRNA sobre el FASTA y solo reescribe las filas `rRNA` del GFF
del sistema y `features.rRNA` del resumen.

Estrategias de `rrna_calling`:

- `exact`: CM completo para 16S/23S/5S; reproduce los bordes de Infernal/Rfam.
- `recall`: CM completo para 5S y HMM para 16S/23S. Es la usada para el panel
  completo: recupera todas las copias y usa bordes de CM en 5S (la mayoría de
  los 167 pares exactos del agregado).
- `balanced`: CM completo para 16S/5S y HMM para 23S.
- `fast`: HMM para los tres (solo localización).

Validación en modo `exact`: `cglutamicum` (18/18), `hpylori` (6/6),
`tpallidum` (6/6), `bburgdorferi` (5/5) y `mgenitalium` (3/3) coinciden
**base a base** con la referencia. `ctrachomatis` coincide en conteo pero no en
bordes (el 23S de la referencia difiere hasta 68 pb del alineamiento local), lo
que se reporta como límite del modelo/versión de Rfam, no como un fallo de
emparejamiento.

## Límites

- **Solo conteos**: no evalúa función, productos, coordenadas, ni identidad de
  secuencia.
- **El gold no es verdad absoluta**: hereda criterios y posibles errores de la
  anotación RefSeq/NCBI usada como referencia.
- **rRNA**: la brecha de 23S/5S se diagnosticó y corrigió (ver arriba); el
  panel quedó en 514/515. Queda una copia de 5S de `paeruginosa` sin llamar y
  no se imputa.
- **`subtilis`**: su GFF se reanotó, pero el manifiesto no declara
  `system_summary`; por eso queda fuera de conteos y sí entra en coordenadas.
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
ilegible, GFF sin clases objetivo y filtrado del manifiesto. La llamada de rRNA
y la reescritura del GFF tienen tests propios en
`tests/test_rrna_calling.py`.

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

Instantánea del **2026-11-25** (cifras vigentes en
`reports/reference_coordinate_parity.md`). El panel tiene 42 registros
elegibles; 41 comparables (incluye `subtilis`, cuyo GFF no tiene
`system_summary`) y 1 `not_available` (`ecoli_k12`, sin GFF del sistema).

Antes/después de reanotar solo rRNA, con los mismos GFF de referencia y el
mismo emparejamiento 1:1:

| clase | TP | FP | FN | precisión | recall | F1 | exactos | % exactos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rRNA antes | 355 | 0 | 190 | 1.0000 | 0.6514 | 0.7889 | 0 | 0.00 % |
| rRNA después | 544 | 0 | 1 | 1.0000 | 0.9982 | 0.9991 | 167 | 30.70 % |

El único FN restante es una copia de 5S de `paeruginosa` (13 en la referencia,
12 en el sistema). No hay FP: la reanotación no añadió features fuera de la
referencia. Los 167 pares exactos corresponden sobre todo a 5S (CM completo);
16S y 23S se localizaron con HMM en la estrategia `recall`, así que sus bordes
quedan a pocos pb (Δinicio medio +1.01 pb, Δfin medio -1.59 pb).

Genomas con paridad rRNA 100 % exacta tras aplicar `--strategy exact`:
`cglutamicum`, `hpylori`, `tpallidum`, `bburgdorferi` y `mgenitalium` (deriva 0
pb). `ctrachomatis` tiene los 6 rRNA emparejados pero bordes distintos en
16S/23S/5S (deriva máxima 68 pb), por diferencias entre la referencia y el
modelo Rfam local.

CDS y tRNA no cambian en esta corrección (se listan en el informe para
contexto): CDS F1 0.9695, tRNA F1 0.9819. El agregado `ALL` pasa de F1 0.9691 a
0.9698; el cambio es pequeño porque rRNA es una fracción minoritaria del total.

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
- **rRNA con `recall`**: 5S usa CM completo (bordes exactos en la mayoría de
  las copias) y 16S/23S usan HMM (localización aproximada). Para bordes exactos en los tres tipos hay que
  usar `--strategy exact`, más lento. La reanotación se aplicó en modo `recall`
  al panel completo y en modo `exact` a seis genomas de validación.
- **Una copia de 5S sin llamar**: `paeruginosa` queda con 12/13; no se imputa
  ni se copia de la referencia.
- **`ctrachomatis`**: conteo correcto pero bordes distintos (hasta 68 pb); se
  atribuye a diferencias de versión/criterio del modelo Rfam, no al
  emparejamiento.

### Reanotación aplicada

```bash
# panel completo (recall; 5S con CM, 16S/23S localizados con HMM)
.venv/bin/python scripts/reference/reannotate_rrna.py \
  --strategy recall --apply --jobs 3 --cpu 2
# validación con CM completo
.venv/bin/python scripts/reference/reannotate_rrna.py \
  --only cglutamicum --only hpylori --strategy exact --apply
.venv/bin/python scripts/reference/reannotate_rrna.py \
  --only bburgdorferi --only mgenitalium --only ctrachomatis \
  --only tpallidum --strategy exact --apply
```

La reanotación **solo** reescribe filas `rRNA` del GFF del sistema y
`features.rRNA` del resumen; no toca CDS/tRNA ni `data/reference_genomes/`.
Los informes antes/después quedan en `reports/rrna_reannotation.*`,
`reports/reference_parity_before.*` y
`reports/reference_coordinate_parity_before.*`. Los `*_before.*` se generaron
sobre una reconstrucción **solo-barrnap** de los GFF del sistema en una carpeta
temporal de trabajo, porque la reanotación se aplicó in situ; sus conteos
coinciden con la columna `antes` de `reports/rrna_reannotation.tsv`.

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
y no modificación de las entradas. La llamada de rRNA (parsers `cmsearch`/
`barrnap`, deduplicación, subconjuntos de modelos y reescritura selectiva del
GFF) se prueba en `tests/test_rrna_calling.py`.
