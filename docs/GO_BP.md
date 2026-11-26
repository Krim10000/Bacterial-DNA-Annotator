# GO Biological Process (BP)

## Regla oficial

Los tres aspectos GO son **independientes a nivel de anotación**:

```text
MF no implica BP
CC no implica BP
```

La única propagación automática válida es:

```text
MF --part_of--> BP
```

según:

```text
GO_REF:0000108
GORULE:0000023
```

y solo si el enlace `part_of` existe en **`go.obo`**. `go-basic.obo` elimina
los enlaces inter-aspecto; por eso este módulo usa `go.obo`.

## Lo que NO se hace

```text
CC → BP          no existe relación
EC → BP          EC solo mapea a MF
Rhea → BP        Rhea solo mapea a MF
KEGG reaction → BP  no hay mapping oficial reacción→BP
ruta → BP        no es automático; requiere IGC y ruta completa
```

Una enzima del ciclo de Krebs **no** justifica `GO:0006099` por sí sola.
`GO:0004333 fumarate hydratase activity` no tiene `part_of` a
`GO:0006099 tricarboxylic acid cycle`.

## Evidencia

| Inferencia | Código | Condición |
|---|---|---|
| MF→BP por `part_of` | IEA / GO_REF:0000108 | enlace en `go.obo` |
| MF→BP manual | IC | curador, `with/from` |
| ortología 1:1 | ISO | fuente experimental |
| similitud | ISS | proteína similar |
| ruta completa | IGC | `with/from` ruta/operón |
| contexto | IGC | solo hipótesis |
| CC→BP | — | prohibido |
| EC/Rhea→BP | — | no soportado |

## Implementación

```text
webapp/go_bp.py
modules/go_bp/module.py
order = 33
```

Insumos:

```text
data/labels/go_mf_cc/<genoma>/gene_go.tsv
data/labels/external/go/go.obo
```

Salida:

```text
data/labels/go_bp/<genoma>/gene_bp.tsv
data/labels/go_bp/<genoma>/bp_summary.json
```

Cada fila:

```text
genome
gene_id
acc
bp_go_id
bp_go_name
mf_go_id
mf_go_name
evidence = IEA
source = GO_REF:0000108
with_from = MF GO ID
confidence = automatic
```

## Cobertura

```text
442 genomas
81.683 filas BP
64.836 genes con BP
```

K-12: 844 filas, 628 genes.
B. subtilis: 719 filas, 612 genes.

## Lo que falta

- BP por ruta completa con IGC (no implementado).
- BP por ortología (ISO/ISS).
- BP por ncRNA: Rfam2GO directo, no inferido.
- Validación contra procesos curados.

## Citas

- GO annotations: https://geneontology.org/docs/go-annotations/
- FAQ: https://geneontology.org/docs/faq/
- GO-CAM: https://geneontology.org/docs/gocam-overview/
- GO_REF:0000108: https://geneontology.org/GO_REF/0000108
- GORULE:0000023: https://github.com/geneontology/go-site/blob/master/metadata/rules/gorule-0000023.md
- GO_REF:0000003 (EC→MF): https://geneontology.org/GO_REF/0000003
- GO_REF:0000116 (Rhea→MF): https://geneontology.org/GO_REF/0000116
- GAF 2.2: https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
