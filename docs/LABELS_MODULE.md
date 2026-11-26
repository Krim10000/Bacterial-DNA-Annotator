# Módulo de etiquetas anidadas

## Principio

El módulo `labels` no anota genes desde cero. **Consume** las salidas de
anotación, metabolismo, regulación y estructura, y construye una red de
etiquetas con procedencia.

```text
gen
→ reacción / EC
→ ruta / módulo
→ proceso / condición
→ regulón
```

## Orden en el pipeline

```text
10  annotation
20  metabolism
30  regulation
32  labels
35  sigma
40  structure
```

`labels` va antes de `sigma` y después de las capas que ya producen función,
reacciones, operones y regulones.

## Contrato del módulo

```python
id = "labels"
order = 32
requires = ("annotation",)
```

No requiere red. No descarga. No reanota.

### Entradas

- `data/regulation/<g>/genes.tsv`
- `data/regulation/<g>/operons.tsv`
- `data/regulation/<g>/operon_genes.tsv`
- `data/regulation/<g>/promoters.tsv` / regulón
- salidas de `metabolism` (EC, Rhea, reacciones)
- `functional.tsv` o equivalente
- dominios/estructura, si existen

### Salidas

```text
data/labels/<g>/label_graph.tsv
data/labels/<g>/gene_labels.tsv
data/labels/<g>/operon_labels.tsv
data/labels/<g>/label_enrichment.tsv
data/labels/<g>/label_summary.json
```

## Esquema

### `label_graph.tsv`

```text
label_id
label_name
parent_id
level
ontology
description
```

Debe ser un DAG. No una lista plana.

### `gene_labels.tsv`

```text
genome
gene_id
label_id
label_name
level
evidence_type: experimental | curated | inferred | predicted | text
source
confidence
value
```

### `operon_labels.tsv`

```text
genome
operon_id
label_id
label_name
score
n_genes
n_genome
enrichment
fdr
```

### `label_enrichment.tsv`

```text
genome
operon_id
label_id
odds_ratio
pvalue
fdr
background
method
```

## Etiquetas mínimas

```text
proceso
├── respuesta a estrés
│   ├── shock térmico
│   ├── estrés de envoltura
│   ├── estrés oxidativo
│   └── limitación de hierro
├── metabolismo
│   ├── nitrógeno
│   ├── azufre
│   └── carbohidratos
├── motilidad
│   └── flagelo
└── diferenciación
    └── esporulación
```

## Reglas de construcción

### 1. No duplicar

- No volver a llamar EC desde sigma.
- No volver a predecir operones.
- No volver a hacer DIAMOND funcional.
- No usar el gold de sigma como etiqueta.

### 2. Evidencia

Cada etiqueta declara:

```text
fuente
tipo
confianza
independencia
```

### 3. Agregación

Gen → operón:

- no doble conteo;
- no sumar parálogos;
- usar media/máximo robusto;
- declarar método.

### 4. Enriquecimiento

Comparar contra:

- genoma;
- operones;
- clado.

Siempre con FDR.

### 5. Propagación

Si un operón es “shock térmico”, también es “respuesta a estrés”.

Con decaimiento y sin inflar todos los padres.

## API para sigma

```python
labels_for_gene(gene_id)
labels_for_operon(operon_id)
enrichment_for_operon(operon_id)
labels_for_sigma(sigma_id)
evidence_path(sigma_id, label_id)
```

`sigma` consume `operon_labels.tsv` y `gene_labels.tsv`.

Salida:

```text
sigma_candidate
context_labels
context_evidence
context_confidence
```

No:

```text
sigma = label
```

## Validación

- DAG sin ciclos;
- labels con fuente;
- FDR;
- ablación;
- provenance holdout;
- K-12/B. subtilis contra regulones conocidos;
- no usar gold de sigma como feature.

## MVP

1. Cargar genes/operones existentes.
2. Construir etiquetas de proceso desde GO/KEGG/EC.
3. Agregar a operón.
4. Calcular enriquecimiento.
5. Exportar TSV/JSON.
6. Integrar como insumo opcional de sigma.

## Lo que no es

- no es un clasificador de sigma;
- no es una tabla reacción→sigma;
- no es una prueba de función;
- no reemplaza RegulonDB/DBTBS.

## Estado

Diseño documentado. Implementación pendiente.
