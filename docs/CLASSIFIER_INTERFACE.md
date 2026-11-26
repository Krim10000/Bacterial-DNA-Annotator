# Interfaz con el clasificador de proteínas

## Entrada obligatoria

```text
data/labels/gene_evidence.tsv
```

Columnas:

```text
genome
gene_id
protein_id
source
evidence
aspect
term_id
term_name
confidence
with_from
coverage
notes
```

## Entrada opcional

```text
data/labels/context/gene_context.tsv
```

Columnas:

```text
genome
gene_id
operon_id
position
n_members
n_annotated
neighbor_genes
context_terms
source
evidence
confidence
```

## Salida esperada

```text
data/classifier/protein_classes.tsv
```

Columnas:

```text
genome
gene_id
protein_id
protein_class
class_probability
class_reason
route
source
evidence
with_from
```

## Reglas

- No reanotar GO, EC, ortología ni dominios.
- No usar el gold de sigma como feature.
- Separar train/test por genoma o filo.
- Declarar `no_data` cuando no haya evidencia.
- La probabilidad no es confirmación.

## Coordinación

El bus es:

```text
docs/AGENT_BUS.md
```

Cada mensaje debe firmarse con el nombre del agente, no con `lead`.
