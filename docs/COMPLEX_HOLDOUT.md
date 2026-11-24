# Holdout de complejos contra ComplexPortal

Este documento describe `tools/complex_holdout.py`, su contrato de evaluación y
por qué, con los datos locales actuales, el holdout independiente es
`not_evaluable`.

## Propósito

Medir cobertura y, cuando sea posible, precisión/recall de los clusters y
transferencias de complejos contra catálogos ComplexPortal locales. La regla
central es **no confundir entrenamiento con generalización**:

- si el catálogo es del genoma de entrenamiento, la evaluación es `in_sample`;
- si el catálogo está en el gold usado para entrenar, no es un holdout;
- si no hay catálogo independiente, el estado es `not_evaluable` y no se
  reporta precisión inventada.

## Por qué K-12 no es un test independiente

`data/complexes/model_k12.genomic.json` declara:

```text
trained_genome      = ecoli_k12
gold                = data/complexes/gold_k12.tsv
gold_complexes      = 273
n_pos               = 1773
n_neg               = 265042
```

Los 273 complejos de `gold_k12.tsv` se usaron para entrenar. Además, los
catálogos crudos disponibles son solo de *E. coli* K-12 (taxón 83333). Por eso:

- evaluar contra `gold_k12.tsv` es **in-sample**;
- evaluar contra los complejos de `gold_k12.tsv` no es holdout;
- no hay ningún genoma fuera del entrenamiento con catálogo local.

## Insumos

| insumo | formato | uso |
|---|---|---|
| `data/complexes/gold_*.tsv` | `complex_id`, `locus_tags`, `n_genes`, ... | catálogo curado; puede ser gold de entrenamiento |
| `data/complexes/raw/complexportal_*.tsv` | `#Complex ac`, `Expanded participant list`, `Taxonomy identifier`, ... | catálogo crudo; se mapea UniProt→locus_tag con `k12_genes.tsv` |
| `data/complexes/k12_genes.tsv` | `locus_tag`, `uniprot`, ... | mapeo de participantes crudos a genes K-12 |
| `data/complexes/model_k12.genomic.json` | JSON | procedencia: genoma y gold de entrenamiento |
| `data/complexes/clusters_ecoli_k12_evidencia.tsv` | `cluster_id`, `genes`, `confidence`, `evidencia`, ... | predicciones K-12 |
| `data/complexes/transfer_*.tsv` | `complex_id`, `genes_k12`, ... | cobertura derivada del catálogo K-12 |

## Estrategias de holdout

### Por genoma

Se buscan catálogos cuyo genoma no esté en `genomas_entrenamiento`. Si hay al
menos uno, se buscan predicciones en:

```text
clusters_<genoma>.tsv
model_<genoma>.complexes.tsv
<genoma>.complexes.tsv
```

El catálogo y las predicciones deben compartir el mismo namespace de genes. Si
no, el estado es `not_evaluable` con el motivo correspondiente.

### Por complejo

Se toman complejos del catálogo crudo que **no** están en el gold de
entrenamiento (`raw_only`). Es un holdout por complejo, pero no por genoma: el
genoma K-12 sigue siendo el genoma de entrenamiento. Se reporta como
`exploratorio_no_independiente`.

Solo se evalúan complejos con ≥2 genes mapeados, porque los clusters del modelo
tienen ≥2 subunidades. Los homómeros y los complejos con un único gen mapeado se
excluyen y se cuentan.

## Estados

| estado | significado |
|---|---|
| `evaluable` | hay catálogo independiente y predicciones mapeadas; se puede calcular precisión/recall |
| `not_evaluable` | no hay catálogo independiente, o hay catálogo pero faltan predicciones; **no** se inventa precisión |
| `in_sample_k12_no_holdout` | evaluación contra el gold de entrenamiento; solo contexto, no generalización |
| `exploratorio_no_independiente` | holdout por complejo con fuga de genoma; solo recall exploratorio |
| `cobertura_derivada_no_holdout` | cobertura de transferencias derivadas del mismo catálogo; no es precisión |

## Métricas

El emparejamiento usa **conjuntos de genes**, no `curado_id`:

- `Jaccard = |P ∩ C| / |P ∪ C|`;
- una predicción es TP si su mejor Jaccard ≥ umbral;
- un complejo del catálogo se recupera si alguna predicción alcanza el umbral;
- `precision = TP_predicciones / predicciones`;
- `recall = complejos_recuperados / complejos_elegibles`;
- `F1 = 2·P·R / (P+R)`.

Umbrales por defecto: `0.3, 0.5, 0.7, 0.9`; umbral principal `0.5`.
El catálogo elegible excluye complejos con <2 genes mapeados.

Estas métricas no evalúan estequiometría, isoformas, interfaces ni interacción
física.

## Resultado actual (instantánea)

Con los datos locales:

| dato | valor |
|---|---:|
| catálogos locales | 2 (`gold_k12`, `complexportal_83333`) |
| complejos crudos | 324 |
| complejos gold | 273 |
| complejos raw-only | 51 |
| raw-only con ≥2 genes mapeados | 10 |
| genomas de entrenamiento | `ecoli_k12` |
| catálogos independientes | 0 |
| **holdout independiente** | **`not_evaluable`** |

El holdout por complejo (exploratorio) da recall **0,600** a Jaccard ≥ 0,5
(6/10), pero no es generalización: el genoma es K-12 y la selección raw-only no
es aleatoria.

El contexto in-sample K-12 (contra los 273 complejos gold) da, a Jaccard ≥ 0,5:
precisión **0,174**, recall **0,725**, F1 **0,281**. Son cifras in-sample y no
deben citarse como rendimiento del modelo.

La cobertura de transferencias es **258/273** complejos gold (94,5 %) y
**701/731** genes gold (95,9 %), con 255 complejos completos, 3 parciales y 15
sin cobertura. Es cobertura derivada, no precisión.

### Qué se puede y qué no se puede concluir

- **No** se puede concluir precisión/recall de generalización: no hay catálogo
  independiente.
- **Sí** se puede reportar cobertura del catálogo, el recall exploratorio de los
  10 complejos raw-only y el contexto in-sample K-12.
- En el desglose in-sample por evidencia, `nuevo` tiene precisión 0 **por
  construcción** (no está en el gold); no es un error del modelo.
- El recall del holdout por complejo no es una estimación insesgada: la
  selección raw-only excluye homómeros y complejos con un solo gen mapeado.

Los números completos están en `reports/complex_holdout.json` y
`reports/complex_holdout.md`.

## Cómo ejecutar

```bash
.venv/bin/python -m pytest tests/test_complex_holdout.py -q
.venv/bin/python tools/complex_holdout.py \
    --out-json reports/complex_holdout.json \
    --out-md reports/complex_holdout.md
```

Opciones útiles:

```bash
--complexes-dir data/complexes
--raw-dir data/complexes/raw
--model-json data/complexes/model_k12.genomic.json
--clusters data/complexes/clusters_ecoli_k12_evidencia.tsv
--thresholds 0.3,0.5,0.7,0.9
--main-threshold 0.5
```

El script no usa red y solo escribe los reportes indicados.

## Cómo habilitar un holdout independiente

1. Elegir un taxón con entrada en ComplexPortal y con genoma en los 400.
2. Descargar y mapear UniProt→locus_tag:

   ```bash
   python scripts/complexes/fetch_complexportal.py \
     --taxid <taxid> \
     --functional data/annotations_ab_initio/<genoma>/functional.tsv \
     --genes data/regulation/<genoma>/genes.tsv
   ```

3. Guardar `data/complexes/gold_<taxid>.tsv` (o `gold_<genoma>.tsv`) y
   `data/complexes/raw/complexportal_<taxid>.tsv`.
4. Asegurar que `<genoma>` no es `ecoli_k12` ni otro genoma de entrenamiento.
5. Añadir predicciones en el mismo namespace (`clusters_<genoma>.tsv` o
   `model_<genoma>.complexes.tsv`). Si el catálogo se llama
   `gold_<taxid>.tsv`, el evaluador usa `<taxid>` como etiqueta de genoma; o
   renombrar el catálogo a `gold_<genoma>.tsv` o nombrar las predicciones
   `clusters_<taxid>.tsv`.
6. Reejecutar. El estado pasará a `evaluable` solo si el catálogo es
   independiente y las predicciones están mapeadas.

## Límites

- Solo hay catálogo curado de K-12; no hay gold para las otras especies.
- Los clusters `model_<genoma>.complexes.tsv` de los 400 genomas no se evalúan
  aquí: no hay catálogo ComplexPortal de esas especies.
- K-12 y los 273 complejos gold son entrenamiento; la evaluación K-12 es
  in-sample.
- El holdout por complejo raw-only no es aleatorio: excluye homómeros y
  complejos con un solo gen mapeado.
- El mapeo crudo UniProt→locus_tag depende de `k12_genes.tsv`; las isoformas
  se colapsan al acceso base.
- El emparejamiento por Jaccard no evalúa estequiometría ni interfaces.
- Las transferencias RBH no tienen cobertura de alineamiento y no demuestran
  interacción física.
- No hay validación experimental en otras especies.

## Referencias locales

- `docs/COMPLEX_CATALOGS.md`
- `docs/COMPLEX_TRANSFER.md`
- `docs/COMPLEX_MODEL_400.md`
- `docs/COMPLEX_COVERAGE.md`
