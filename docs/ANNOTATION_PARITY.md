# Paridad de anotación base

Herramienta: [`tools/annotation_parity.py`](../tools/annotation_parity.py).
Informe actual: [`reports/annotation_parity.md`](../reports/annotation_parity.md)
y [`reports/annotation_parity.json`](../reports/annotation_parity.json).
Tests: [`tests/test_annotation_parity.py`](../tests/test_annotation_parity.py).

## Para qué sirve

Compara **nuestra** anotación base (GFF3 de
`data/annotations_ab_initio/<genoma>/annotation.gff3`, o un GFF de `out_fl`)
contra un **GFF3 de referencia local** (NCBI/RefSeq, Bakta, Prokka, etc.) y
describe diferencias que se pueden revisar: genes ausentes o extra,
fragmentación, fusiones, cambios de hebra, desplazamientos de borde, ORFs
pequeños, ncRNA y productos/hipotéticas.

No es un benchmark de superioridad entre anotadores. Para eso existe
[`tools/benchmark_external.py`](../tools/benchmark_external.py), que compara
varias herramientas y calcula precisión/recall contra una referencia. Aquí el
objetivo es **diagnosticar nuestra anotación** y decir qué revisar.

## Uso

Sin argumentos, el genoma de ejemplo y la referencia se autodescubren:

```bash
.venv/bin/python tools/annotation_parity.py \
  --genome GCF_000006685.1_ASM668v1
```

Salidas por defecto:

```text
reports/annotation_parity.md
reports/annotation_parity.json
```

Referencia explícita (NCBI, Bakta o Prokka):

```bash
.venv/bin/python tools/annotation_parity.py \
  --ours data/annotations_ab_initio/GCF_000006685.1_ASM668v1/annotation.gff3 \
  --reference /ruta/al/referencia.gff3 \
  --md reports/annotation_parity.md \
  --json reports/annotation_parity.json
```

GFF de `out_fl`:

```bash
.venv/bin/python tools/annotation_parity.py \
  --ours out_fl/annotation/demo.gff3 \
  --reference data/annotations/GCF_000006685.1_ASM668v1.gff.gz
```

Si los nombres de contig del GFF de `out_fl` no coinciden con la referencia, el
informe dirá `contig_mismatch` y `coordinates.status = not_evaluable`; no es un
fallo del comando, es una comparación de coordenadas que no se puede hacer.

Sin referencia:

```bash
.venv/bin/python tools/annotation_parity.py --ours /ruta/anotacion.gff3
```

El estado será `not_evaluable`: no se calculan precision/recall ni diferencias,
y no se inventan números. El informe solo describe nuestra anotación.

Opciones útiles:

| opción | defecto | efecto |
|---|---|---|
| `--ours` | autodescubre | GFF propio o directorio |
| `--reference` | autodescubre | GFF de referencia o directorio |
| `--genome` | `GCF_000006685.1_ASM668v1` sin args | nombre para autodescubrir |
| `--ours-products` | `auto` | `auto`, `none` o TSV `id<TAB>...<TAB>producto` |
| `--overlap-threshold` | `0.5` | solapamiento recíproco mínimo para emparejar |
| `--boundary-alert` | `30` | nt de diferencia de borde que genera alerta |
| `--coord-tolerance` | `0` | nt de tolerancia para coordenadas equivalentes |
| `--small-orf-max-aa` | `150` | umbral estricto de ORF pequeño |
| `--max-details` | `200` | máximo de ítems detallados por bloque en JSON |
| `--exclude-pseudo/--no-exclude-pseudo` | excluir | excluir CDS con `pseudo=true` |
| `--root` | raíz del repo | raíz para autodescubrimiento |

## Autodescubrimiento

Nuestra anotación, en orden:

1. `data/annotations_ab_initio/<genoma>/annotation.gff3`
2. `data/annotations_ab_initio/<genoma>/annotation.gff`
3. `data/annotations_ab_initio/<genoma>/<genoma>.gff3`
4. `out_fl/annotation/<genoma>.gff3`
5. `out_fl/<genoma>.gff3`
6. `out_fl/<genoma>/annotation.gff3`

Referencia, en orden:

1. `data/annotations/<genoma>.gff.gz`
2. `data/annotations/<genoma>.gff3.gz`
3. `data/annotations/<genoma>.gff`
4. `data/annotations/<genoma>.gff3`
5. las mismas variantes con `_genomic`
6. `data/references/<genoma>.*`

Se leen GFF normales y gzip (por extensión o por magic bytes). No hay red ni
descargas.

## Qué mide

### Clases

- **CDS**: se excluyen pseudogenes por defecto (`pseudo=true`,
  `pseudogene=true`). Los CDS multiparte con el mismo `ID` se agrupan en una
  sola feature; se conserva el sobre mínimo–máximo y la suma de segmentos.
- **tRNA** y **rRNA**.
- **ncRNA**: `ncRNA`, `misc_RNA`, `tmRNA`, `RNase_P_RNA`, `SRP_RNA`, etc.

Las coordenadas GFF3 se conservan 1-based inclusive.

### Conteos y emparejamiento

- Conteos por clase, delta y features dentro de contigs compartidos.
- Emparejamiento 1:1 voraz por **solapamiento recíproco** mínimo
  (`--overlap-threshold`, 0.5 por defecto), dentro de
  `(clase, contig, hebra)`.
- `recall = emparejados / features de referencia en alcance`.
- `precision = emparejados / features nuestras en alcance`.
- `exactos`: mismo `start`, `end` y hebra.
- `--coord-tolerance` cuenta además pares equivalentes dentro de N nt.
- No se alinean secuencias ni se hace liftover.

### Diferencias accionables

El informe separa dos bloques:

**Accionables** (`differences.actionable`): cambios que conviene revisar.

| categoría | significado |
|---|---|
| `missing_cds` | CDS de la referencia ausente en nuestra anotación |
| `extra_cds` | CDS nuestro sin contraparte en la referencia |
| `fragmented_in_ours` | un gen de referencia cubierto por varias features nuestras |
| `merged_in_ours` | una feature nuestra que abarca varios genes de referencia |
| `strand_conflict` | mismo locus anotado en hebra opuesta |
| `boundary_shift` | borde 5'/3' desplazado más de `--boundary-alert` nt |
| `length_discordance` | largo muy distinto entre features emparejadas |
| `multipart_structure_mismatch` | mismo sobre, segmentos multiparte distintos |
| `weak_overlap` | solapamiento recíproco entre 0.1 y el umbral |
| `missing_trna`, `missing_rrna`, `missing_ncrna` | clase ausente en nuestra anotación |
| `extra_trna`, `extra_rrna`, `extra_ncrna` | clase extra en nuestra anotación |
| `hypothetical_disagreement` | un lado dice hipotética y el otro no |
| `contig_mismatch` | no hay ningún contig compartido |
| `contig_only_reference`, `contig_only_ours` | contigs presentes solo en un lado |

**Revisión** (`differences.review`): diferencias de nomenclatura o
limitaciones de la referencia, no errores de anotación.

| categoría | significado |
|---|---|
| `product_mismatch` | dos productos anotados con poco solapamiento de tokens |
| `pseudo_in_reference` | nuestro CDS corresponde a un pseudogén de la referencia |
| `pseudo_in_ours` | la referencia tiene CDS donde nosotros marcamos pseudogén |
| `class_absent_in_reference` | la referencia no anota esa clase; no se evalúa |
| `ncrna_type_mismatch` | subtipo ncRNA distinto en un par emparejado |

Cada ítem incluye `category`, `severity` (`high`/`medium`/`low`), `feature_class`,
`seqid`, coordenadas de referencia y nuestras, producto, métricas y un mensaje.
El JSON guarda hasta `--max-details` por bloque; los conteos por categoría son
completos.

### Productos e hipotéticas

- Primero se usa el atributo `product` del GFF.
- Si no hay, se busca un `functional.tsv` vecino (p. ej. el Swiss-Prot del
  pipeline) o un TSV indicado con `--ours-products`.
- `% hypothetical` se calcula sobre los CDS **con** producto, con un
  clasificador conservador (`hypothetical`, `uncharacterized`, `unknown
  protein`, `DUF####`, `UPF####`). No usa `putative` a secas.
- Si no hay producto en un lado, el porcentaje es `null`/`—`; no se infiere.
- En pares con producto en ambos lados se reportan coincidencia literal,
  Jaccard de tokens y desacuerdos. Los nombres distintos entre bases se
  separan como revisión, no como error.

### ORFs pequeños

`--small-orf-max-aa` (150 por defecto) usa una **estimación**:
`span_nt // 3 - 1`, salvo que `partial` indique que falta el extremo 3'. No se
traduce la secuencia. Se cuentan ORFs pequeños en cada lado, cuántos quedaron
emparejados y cuántos son ausentes/extra accionables.

### ncRNA

Se comparan totales y subtipos (`RNase_P_RNA`, `SRP_RNA`, `tmRNA`, etc.). Si
la referencia no trae ncRNA, el estado es `reference_absent` y no se penaliza
a nuestra anotación. Si la referencia trae ncRNA y nosotros no, se reporta
`missing_ncrna`.

### Contigs

Solo se comparan features de contigs compartidos. Si no hay ninguno:

- los conteos siguen siendo descriptivos;
- `coordinates.status = not_evaluable`;
- se emite `contig_mismatch` y no se enumeran cientos de ausentes/extra falsos.

## Estados

- `evaluated`: hay referencia utilizable.
- `not_evaluable`: no hay referencia local, la ruta no existe o el GFF no tiene
  features comparables.
- Si hay referencia pero no hay contigs compartidos, el informe queda
  `evaluated` (conteos descriptivos) pero `coordinates.status = not_evaluable`.
- Por clase: `evaluated`, `not_present`, `reference_absent`, `ours_absent`,
  `not_evaluable_contigs`.
- Los valores no calculados son `null` en JSON y `—` en Markdown.

## Salida JSON

Bloques principales:

```text
status, reason, genome, generated_at, command
inputs.{ours,reference}
parameters
contigs
summary.{CDS,tRNA,rRNA,ncRNA}
coordinates
products.{ours,reference,paired}
small_orfs
ncrna
differences.{actionable,review}
limitations
```

`inputs.*.sha256` permite verificar que el informe corresponde a los archivos
citados. `command` reproduce la ejecución.

## Ejemplo de lectura

En el caso actual (`GCF_000006685.1_ASM668v1` contra el GFF de RefSeq local):

- CDS nuestros **902** vs referencia **891** (899 líneas CDS menos 7
  pseudogenes excluidos y 1 CDS multiparte agrupado);
- **884** CDS emparejados, **803** con coordenadas exactas;
- tRNA **37/37** y rRNA **6/6** exactos;
- ncRNA de referencia **3** (`RNase_P_RNA`, `SRP_RNA`, `tmRNA`) sin contraparte
  en nuestra anotación;
- diferencias accionables separadas en bordes, largos, hipotéticas,
  ausentes/extra y ncRNA;
- productos nuestros desde `functional.tsv` (Swiss-Prot), no desde el GFF.

El informe es la fuente de detalle y de las cifras completas; este documento
solo resume el orden de magnitud para no quedar desactualizado.

## Límites

- La comparación es de anotación, no de función, homología ni identidad.
- El emparejamiento voraz no es un óptimo bipartito.
- Un GFF de referencia puede omitir tRNA/rRNA/ncRNA o incluir pseudogenes.
- Los nombres de producto dependen de la base y la convención de cada
  herramienta; por eso se separan en revisión.
- Nuestros productos pueden venir de Swiss-Prot (`functional.tsv`), no del GFF;
  la comparación con RefSeq es orientativa.
- El largo en aa es una estimación desde el ADN.
- No se descargan referencias, no se ejecutan Bakta/Prokka/DFAST y no se usan
  GPU ni red.
- Un resultado `not_evaluable` no es una medida de calidad: es ausencia de
  referencia o de contigs comparables.

## Tests

```bash
.venv/bin/python -m pytest tests/test_annotation_parity.py -q
```

Los tests usan GFF sintéticos y cubren paridad exacta, tolerancia de
coordenadas, ausentes/extra, ORFs pequeños, bordes, largo, multiparte,
fragmentación, fusión, hebra, productos, ncRNA, pseudogenes, gzip,
autodescubrimiento, `out_fl`, truncado de detalles, contigs incompatibles y el
caso sin referencia (`not_evaluable`).
