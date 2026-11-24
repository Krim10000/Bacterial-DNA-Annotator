# Benchmark reproducible por capas

**Generado:** 2026-11-24T01:17:29Z

**Raíz:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`

**Alcance:** Solo datos locales; sin red, sin reentrenar y sin modificar pipeline/GUI/docs.

**Comando:** `python tools/benchmark_layers.py --md reports/layer_benchmark.md --json reports/layer_benchmark.json`

Este informe no cierra brechas ni inventa gold: cada capa declara qué se pudo medir y qué queda como `not_evaluable`.

## Resumen

| capa | estado | métrica principal | no evaluable |
|---|---|---|---|
| A/B/C | parcial | 1.337.745 proteínas; A 42,4 %, B 9,2 %, C 48,4 % | precisión de función/dominios |
| Membrana | parcial | 319.715 filas; high+medium 33,6 % | precisión de localización |
| Complejos | interno (no ciego) | 86/273 gold 1:1 a J≥0,8; 86/960 clusters con mejor J≥0,8 | precisión predictiva ciega |
| oriC | parcial | validated 271, consensus 77, anomalous 42; DoriC 306 | exactitud independiente |

## 1. A/B/C

**Estado:** parcial. Se miden conteos y cobertura descriptiva; no hay gold local de dominios/función.

### Conteos

| etiqueta | proteínas | % |
|---|---|---|
| A · función nombrada por homología | 566.957 | 42,4 % |
| B · homólogo no caracterizado / DUF | 123.200 | 9,2 % |
| C · sin hit Swiss-Prot | 647.588 | 48,4 % |

- Genomas: **400**.
- Proteínas: **1.337.745**.
- `batch_summary.totals.cds`: **1.337.745**; coincide con las filas: **True**.

### Etiqueta × fuerza de evidencia

| etiqueta | strong | moderate | weak | none | unknown | total |
|---|---|---|---|---|---|---|
| A · función nombrada por homología | 76.124 | 490.833 | 0 | 0 | 0 | 566.957 |
| B · homólogo no caracterizado / DUF | 0 | 13.047 | 110.153 | 0 | 0 | 123.200 |
| C · sin hit Swiss-Prot | 0 | 0 | 0 | 647.588 | 0 | 647.588 |

### Cobertura descriptiva

- Proteínas por genoma: mínimo **531**, mediana **2.863,0**, máximo **11.514**.
- Porcentaje A por genoma: mínimo **6,3 %**, mediana **45,6 %**, máximo **65,5 %**.
- Con UniProt: **690.157** (51,6 %).
- Con gen: **690.157** (51,6 %).
- `subunit_hint=named_subunit`: **80.839** (6,0 %).

### No evaluable

- **precisión/recall de función nombrada (A) contra gold**: No existe una tabla local de gold de función para los 400 genomas. 'A' significa hit Swiss-Prot con nombre, no función confirmada experimentalmente.
- **precisión/recall de dominios (B/C) contra gold**: No hay gold local de Pfam/InterPro para estos genomas. El panel Pfam de membrana es un panel de búsqueda, no verdad de referencia.

### Limitaciones

- A/B/C es una ruta por homología, no una llamada de función confirmada.
- 'C' significa sin hit Swiss-Prot; no significa que la proteína no tenga función.
- Los conteos dependen de la versión de Swiss-Prot/DIAMOND usada por el pipeline.
- No se evalúa la identidad de la función (solo si hay un nombre asignado).
- Existe anotación NCBI local (data/annotations) con product, pero se usa para estructura/llamada de genes; no se trata como gold de dominios/función.

## 2. Membrana

**Estado:** parcial. La validación local es estructural; no hay gold biológico de localización.

### Composición

| tipo | filas | % |
|---|---|---|
| tm | 104.828 | 32,8 % |
| signal_peptide | 128.208 | 40,1 % |
| lipoprotein | 18.560 | 5,8 % |
| transporter | 58.069 | 18,2 % |
| secretion | 10.050 | 3,1 % |

- Genomas: **400**; filas: **319.715**; proteínas únicas: **319.715**.

### Confianza y transportadores

| confianza | filas | % |
|---|---|---|
| high | 54.951 | 17,2 % |
| medium | 52.485 | 16,4 % |
| low | 212.279 | 66,4 % |

- high+medium: **107.436** (33,6 %).
- Transportadores: **58.069**; con sustrato explícito **5.456** (9,4 %).

### Validación estructural local

- Fuente: `data/membrane/validation_full.json` (checked_at `2026-11-21T03:54:32Z`).
- ok: **True**; processed_ok: **True**.
- genomas: **400**; filas: **319.715**; features GFF: **319.715**.
- errores: **0**; pendientes: **0**; genomas_ok: **400**.
- Consistente con los TSV leídos: **True**.

### No evaluable

- **precisión de localización (TM/SP/lipoproteína/secreción)**: La validación local (05_validate_outputs.py) comprueba cabecera, coordenadas, protein_id y correspondencia TSV-GFF; no hay gold biológico local (TMHMM/SignalP/LipoP/UniProt) para medir aciertos.

### Limitaciones

- TM/SP/lipobox son heurísticos; no se usaron TMHMM/SignalP/LipoP.
- El panel Pfam es de 255 modelos; algunas familias están subrepresentadas.
- 66,4 % de las filas son confianza 'low'; 'low' es hipótesis, no localización confirmada.
- La validación estructural no demuestra que la asignación sea correcta.

## 3. Complejos

**Estado:** benchmark interno y **no ciego**. El modelo se ajustó con el mismo gold K-12; el solapamiento no es precisión predictiva.

### Conjuntos

| conjunto | n | genes | detalle |
|---|---|---|---|
| gold K-12 | 273 | 731 | evidencia: experimental 133, inferida 140 |
| clusters K-12 | 960 | 2.934 | confianza: alta 498, baja 233, media 229 |

### Concordancia con gold por umbral de Jaccard

| Jaccard | gold cubierto (cualquier cluster) | pares 1:1 | recall gold | concordancia clusters | recall genes | pureza genes |
|---|---|---|---|---|---|---|
| 0,8 | 88 (32,2 %) | 86 | 31,5 % | 9,0 % | 36,5 % | 100,0 % |
| 0,5 | 198 (72,5 %) | 167 | 61,2 % | 17,4 % | 68,8 % | 84,2 % |

`concordancia clusters` = pares 1:1 / total de clusters. Se evita llamarla precisión porque los clusters no catalogados no son falsos positivos conocidos.

### Cobertura del gold a Jaccard ≥ 0,8

| grupo | gold | cubiertos | % |
|---|---|---|---|
| evidencia: experimental | 133 | 35 | 26,3 % |
| evidencia: inferida | 140 | 53 | 37,9 % |
| tamaño: 2 | 123 | 23 | 18,7 % |
| tamaño: 3-4 | 112 | 57 | 50,9 % |
| tamaño: 5+ | 38 | 8 | 21,1 % |

`cubiertos` = al menos un cluster con Jaccard ≥ 0,8; no es el emparejamiento 1:1.

### Etiquetas existentes y artefactos

- Etiquetas `evidencia` de los clusters: **documentado 86, nuevo 757, parcial 84, subcomplejo 33**.
- Artefactos: **400** `transfer_*.tsv` y **400** `model_*.complexes.tsv`; no hay gold para medirlos.

### Advertencia de circularidad

- El modelo de pares se ajustó con gold_k12.tsv y los clusters K-12 se generaron aplicándolo al mismo genoma. El solapamiento con el gold no es una validación ciega.
- No mide: Precisión predictiva en genomas no vistos ni complejos ausentes del catálogo ComplexPortal.

### No evaluable

- **precisión predictiva ciega / externa**: El gold K-12 se usó para ajustar el modelo y para etiquetar los clusters; no hay un conjunto K-12 held-out ni otro genoma con gold de complejos.
- **falsos positivos en clusters 'nuevo'**: No hay verdad para complejos fuera de ComplexPortal; un cluster no catalogado no puede contarse como falso positivo.

### Limitaciones

- Benchmark interno y no ciego: el solapamiento con gold está favorecido por construcción.
- Gold solo de E. coli K-12 (ComplexPortal), 273 complejos con >=2 genes.
- Complejos grandes (50S, 30S, complejo I) quedan casi todos fuera de Jaccard >=0,8.
- No se evalúa estequiometría ni interacción física; 'documentado' es coincidencia de subunidades.
- Los clusters 'nuevo'/'parcial' no son necesariamente falsos: el mundo está abierto.

## 4. oriC

**Estado:** parcial. DoriC se separa del consenso; `validated` usa DoriC como ancla y no se presenta como exactitud independiente.

### Estados

| estado | genomas | % |
|---|---|---|
| anomalous | 42 | 10,5 % |
| consensus | 77 | 19,2 % |
| lineage_specific | 10 | 2,5 % |
| validated | 271 | 67,8 % |

- DoriC disponible: **306** (76,5 %).
- DoriC local: release **DoriC 14.0**, CSV `DoriC14_release_20261012.csv`, filas **530.178**.
- Sin DoriC: **94** (23,5 %); estados: anomalous 7, consensus 77, lineage_specific 10.

### `validated_by`

| validated_by | genomas | % |
|---|---|---|
| DoriC + GC skew | 6 | 1,5 % |
| GC skew | 7 | 1,8 % |
| GC skew (linaje sin dnaA) | 10 | 2,5 % |
| ORCA | 35 | 8,8 % |
| ORCA + DoriC + dnaA | 265 | 66,2 % |
| ORCA + proximidad dnaA | 77 | 19,2 % |

### Distancia a DoriC

| grupo | n | mediana pb | ≤1 kb | ≤10 kb | ≤50 kb |
|---|---|---|---|---|---|
| validated (anclado a DoriC) | 271 | 4.755,0 | 87 (32,1 %) | 176 (64,9 %) | 267 (98,5 %) |
| anomalous con DoriC | 35 | 257.863,0 | 0 (0,0 %) | 0 (0,0 %) | 0 (0,0 %) |
| todos con DoriC | 306 | 6.654,5 | 87 (28,4 %) | 176 (57,5 %) | 267 (87,2 %) |
| ORCA-only vs DoriC | 306 | 16.544,0 | 39 (12,8 %) | 126 (41,2 %) | 211 (69,0 %) |

- `validated` con DoriC: **271**; usan DoriC: **271**; fuera de 50 kb: **4**.
- ORCA-only vs DoriC: `delta_bp` recomputado coincide: **True** (mismatches=0).
- `consensus`: **77**; con DoriC: **0**; distancia a verdad externa: **not_evaluable**.
- `lineage_specific`: **10**; con DoriC: **0**; distancia a verdad externa: **not_evaluable**.

Ejemplos de `validated` fuera de 50 kb: `GCF_000484535.1_ASM48453v1` (1.304.960 pb), `GCF_055387105.1_ASM5538710v1` (1.042.460 pb), `GCF_055390195.1_ASM5539019v1` (681.318 pb), `GCF_021057185.1_ASM2105718v1` (267.813 pb)

Estos 4 casos sugieren un desajuste de coordenadas/replicón; no se corrigen en este benchmark y no deben presentarse como aciertos.

### No evaluable

- **exactitud independiente del pipeline oriC**: Los 271 estados validated usan DoriC (271 con DoriC explícito en validated_by); no hay ningún validated sin DoriC. La distancia mide consistencia con el ancla, no acierto independiente.
- **distancia a verdad externa para consensus**: Los 77 consensus no tienen DoriC (0 con DoriC). El estado combina ORCA y proximidad a dnaA, no una verdad externa.
- **distancia para lineage_specific y anomalous sin DoriC**: lineage_specific=10, anomalous sin DoriC=7; no hay verdad externa para medir distancia.

### Limitaciones

- DoriC es una referencia curada, no una verdad universal; 94/400 genomas no tienen DoriC.
- DoriC puede cruzar el origen (start > end); se usa start como ancla, igual que 03_orca_predict.py.
- ORCA asume linealidad; no se rota el contig en este benchmark.
- La tabla de estados no trae replicón: para la distancia se usa el contig más largo como cromosoma.
- El estado 'validated' no es accuracy: 271/271 usan DoriC como ancla.
- El estado 'anomalous' se define por discrepancia; su distancia no mide desempeño global.

## Fuentes

| capa | archivo |
|---|---|
| abc | `data/annotations_ab_initio` |
| abc | `data/annotations_ab_initio/batch_summary.json` |
| membrane | `data/membrane/full` |
| membrane | `data/membrane/gff/full` |
| membrane | `scripts/membrane/05_validate_outputs.py` |
| membrane | `data/membrane/validation_full.json` |
| complexes | `data/complexes/gold_k12.tsv` |
| complexes | `data/complexes/clusters_ecoli_k12_evidencia.tsv` |
| complexes | `data/complexes/clusters_ecoli_k12.tsv` |
| complexes | `scripts/complexes/annotate_evidence.py` |
| complexes | `data/complexes/transfer_*.tsv` |
| complexes | `data/complexes/model_*.complexes.tsv` |
| oric | `data/orca_predictions_fix/oriC_status.tsv` |
| oric | `data/orca_predictions` |
| oric | `data/annotations_ab_initio` |
| oric | `data/doric` |
| oric | `scripts/oric/03_orca_predict.py` |

## Limitaciones globales

- Solo se usaron datos locales; no hay descargas, entrenamiento ni gold externo.
- Las métricas describen estos 400 genomas y este snapshot de datos, no una población.
- Un `not_evaluable` no es un cero: significa que falta la verdad de referencia.
- No se modificó el pipeline, la GUI ni la documentación.
