# Interfaz L2 → anotación/integración (`lead`)

Contrato de consumo de la capa de regulación que produce `regulon-builder`.
Estado: **v0**. Sujeto a cambio durante el sprint de L2; el gate está pendiente
(ver `data/regulation/GATE.md` cuando exista).

Regla: `lead` **lee** estas salidas, no las edita. Si hace falta un campo, se
pide al `regulon-builder`.

## Rutas

    data/regulation/<genoma>/genes.tsv       # tabla génica estable
    data/regulation/<genoma>/operons.tsv     # unidades transcripcionales
    data/regulation/<genoma>/promoters.tsv   # promotor por operón + sigma
    data/regulation/<genoma>/regulon.json    # TF -> operón (bloqueado por gate)
    data/regulation/GATE.md                  # veredicto del gate L2

Para K-12 el `<genoma>` es `ecoli_k12`, no el GCF.

## `genes.tsv` — una fila por gen

Columnas observadas (23):

    genome contig start end strand gene_id locus_tag gene_name product ec acc
    fam partial rbs_motif rbs_source start_type conf rna gene_type source
    status evidence attempted confidence

Notas de uso:

- **`strand` es `1`/`-1`**, no `+`/`-`. No asumir signo textual.
- `locus_tag` es la clave de unión con RegulonDB para K-12 (`b0001`…), y con
  `data/parsed_annotations/<gen>_labels.json` para el resto.
- `ec` viene vacío en buena parte del set de 400; el enriquecimiento por
  `/EC_number` del GenBank lo hace `scripts/annotation/14_parse_gff.py`
  (~32% de los CDS), no esta capa.
- `status`/`evidence`/`attempted`/`confidence` usan el vocabulario de
  `scripts/ghost/ghost/status.py`.

## `operons.tsv` — una fila por operón

Columnas (19):

    operon_id label contig start end strand n_genes genes gene_names gene_ids
    n_cds span_bp min_gap max_gap n_links_evidence status evidence attempted
    confidence

- `operon_id` = `<genoma>:<contig>:<start>:<strand>`.
- `genes` = b-numbers separados por coma; `gene_names` y `gene_ids` alineados
  por posición con `genes`.
- `n_genes` cuenta CDS; `n_cds` coincide salvo genes de ARN dentro del operón.
- Incluye operones de un solo gen: filtrar `n_genes >= 2` para análisis de
  co-transcripción.

## `promoters.tsv` — un promotor por operón (o `unresolved`)

    operon_id label contig strand box35 mm35 box10 mm10 spacer ext10 score
    dist_10_to_start n_sites sigma sigma_status status evidence attempted
    confidence

- `sigma_status=unresolved` cuando no hay asignación respaldada (hoy la
  mayoría). **No interpretar `unresolved` como "sin promotor"**: significa
  "no hay evidencia suficiente con los modelos actuales".
- `score` 0 con cajas vacías = no se encontró motivo, no "promotor débil".

## Qué consume el `lead` y para qué

1. **Contexto de genes fantasma / desiertos** (`scripts/ghost/`): si un
   candidato cae dentro de un operón anotado, su contexto de co-transcripción
   es evidencia adicional. Join por `contig` + coordenadas + `strand`.
2. **Brechas metabólicas**: un gen mapeado a reacción que pertenece a un
   operón con `status` alto sube la confianza de que la enzima se exprese.
3. **Informe final**: contar operones y promotores con evidencia por genoma.

## Límites que hay que arrastrar al informe

- Gate L2 **no superado** (pares F1 0,793; exacto 0,606). El número de
  operones es correcto; el error está en las fronteras de los largos.
- Cobertura de promotores ~24%; sigma resuelto solo por nombre de producto.
- L3 bloqueado hasta `GATE.md`.
- La interfaz puede cambiar: verificar con `head -1` antes de un join en
  producción.
