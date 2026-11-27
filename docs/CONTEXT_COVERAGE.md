# Cobertura por contexto de operón

## Pregunta

¿Cuántos genes `no_data` están dentro de un operón con vecinos anotados?

## Resultado actual

```text
no_data:                       859.533
con contexto:                  205.801
  internos del operón:          83.793
  en borde del operón:         122.008
sin archivo de operones:       152.719
operón sin vecino anotado:     476.531
sin contexto:                   24.482
```

## Interpretación

- **205.801** genes tienen un contexto que puede priorizar su análisis.
- **83.793** están en el medio de un operón con función conocida.
- No se les asigna GO directo.
- Se guarda:

```text
source = operon_context
evidence = IGC
confidence = low
context_terms = términos de vecinos anotados
```

## Archivo

```text
data/labels/context/gene_context.tsv
data/labels/context/context_summary.json
```

Columnas:

```text
genome
gene_id
operon_id
position          internal | edge
n_members
n_annotated
neighbor_genes
context_terms
source
evidence
confidence
```

## Regla

Un vecino anotado **no demuestra** que el gen huérfano tenga esa función.

Sirve para:

- priorizar;
- agrupar;
- proponer una hipótesis;
- alimentar al clasificador.

No sirve para:

- asignar GO fuerte;
- cerrar una brecha;
- confirmar una ruta.

## Siguiente

- usar contexto como feature del clasificador;
- combinarlo con dominio/estructura;
- no fusionarlo como GO directo.

## Contexto conservado por cluster

Además del vecino local, se buscó el mismo cluster en otro genoma:

```text
no_data:                         859.532
con contexto por cluster:         64.330
sin cluster:                     177.319
cluster sin donante anotado:     617.883
```

Archivo:

```text
data/labels/context/gene_cluster_context.tsv
```

Columnas:

```text
genome
gene_id
protein_id
cluster
n_donors
context_terms
source = cluster_context
evidence = IGC
confidence = low
```

No es GO directo; es contexto conservado.
