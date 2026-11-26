# Contrato `gene_evidence.tsv`

## Objetivo

Una sola tabla de evidencia por gen para:

- clasificador de proteínas;
- sigma;
- circuitos de regulación.

## Columnas

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

## Semántica

| columna | significado |
|---|---|
| `source` | archivo o capa que produjo el dato |
| `evidence` | código GO/ECO o `no_data` |
| `aspect` | `F`, `C`, `P` |
| `term_id` | GO/EC/RHEA/pathway |
| `confidence` | high/medium/low/unknown |
| `with_from` | donante o MF de origen |
| `coverage` | direct/inferred/none |
| `notes` | detalle de procedencia |

## Reglas

- No duplicar términos con la misma fuente/evidencia.
- Si no hay evidencia: una fila `no_data`.
- No usar `confidence` como prueba.
- No transferir IEA como experimental.

## Estado

Piloto validado:

```text
K-12 + B. subtilis
13.606 genes
36.622 filas
data/labels/gene_evidence_k12_bsu.tsv
```

El archivo completo de 442 genomas se cortó antes de terminar y quedó como
`gene_evidence.partial.tsv`; **no debe consumirse**.

Script:

```bash
python scripts/labels/21_build_gene_evidence.py \
  --genomes ecoli_k12 subtilis \
  --out data/labels/gene_evidence_k12_bsu.tsv
```

## Pendiente

- ejecutar por lotes de genomas;
- añadir dominios/ortología;
- validar duplicados;
- publicar `gene_evidence.tsv` completo.

## Cobertura completa

```text
442 genomas
1.522.305 genes
3.502.471 filas
1.241.855 filas no_data
280.450 genes con evidencia
```

El archivo completo es:

```text
data/labels/gene_evidence.tsv
```

La ausencia de una fila de término ya no es ambigua: hay una fila `no_data`
por cada gen sin evidencia.
