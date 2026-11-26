# Bus de coordinación entre agentes

**Proyecto:** super anotador bacteriano end-to-end  
**Archivo:** `docs/AGENT_BUS.md`  
**Regla:** este archivo es la fuente de coordinación. No reemplaza contratos de datos ni documentación técnica.

## 1. Cómo usarlo

1. Leer la sección **Contrato de evidencia** antes de producir datos.
2. Revisar **Propiedad de capas** para no duplicar trabajo.
3. Añadir mensajes **solo al final**, en `## Bitácora`, sin borrar ni editar mensajes ajenos.
4. Cada mensaje declara:

```text
fecha | agente | tipo | capa | estado
```

5. Tipos permitidos:

```text
INFO      información nueva
REQUEST   necesito un dato/acción
OFFER     ofrezco un insumo
BLOCKER   no puedo avanzar
DECISION  acuerdo tomado
HANDOFF   entrego salida
```

6. Si hay conflicto, no se resuelve borrando: se añade una entrada `DECISION` con la razón.

## 2. Identificadores canónicos

Usar siempre estos IDs:

```text
genome              GCF_... | ecoli_k12 | subtilis | demo
gene_id             gene_id de data/regulation/<genome>/genes.tsv
protein_id          <genome>|<seq_id>   (FASTA)
operon_id           operon_id del navegador
go_id               GO:#######
ec                  EC x.x.x.x
uniprot_acc         accesión UniProt
rfam                RF#####
structure_id        PDB/AFDB/ESMFold/Kaggle
```

No crear IDs alternativos si ya existe uno.

## 3. Contrato de evidencia

Todo agente que produzca evidencia debe entregar TSV/JSON con estas columnas mínimas:

```text
genome
gene_id
protein_id
source
evidence
aspect              F | C | P | EC | pathway | structure
term_id             GO:... | EC:... | RHEA:... | map...... | pocket id
term_name
confidence
with_from
coverage
notes
```

### Reglas

- `source` es obligatorio.
- `evidence` es obligatorio.
- No usar `confidence` como sustituto de `evidence`.
- No transferir `IEA` como si fuera experimental.
- No promediar fuentes correlacionadas.
- Si no hay dato: `no_data`, no inventar.

## 4. Propiedad de capas

| Capa | Agente responsable | Salida | No debe hacer |
|---|---|---|---|
| Anotación base | pipeline | CDS/ARN/proteínas | duplicar pyrodigal |
| GO MF/CC | labels/coverage | `data/labels/go_mf_cc/<g>/gene_go.tsv` | reanotar Swiss-Prot |
| GO BP | labels/coverage | `data/labels/go_bp/<g>/gene_bp.tsv` | inferir BP desde CC/EC/ruta sin regla |
| GO ARN | labels/coverage | `data/labels/ncrna_go/<g>/ncrna_go.tsv` | usar proteína para ARN |
| Etiquetas de proceso | labels/coverage | `data/labels/<g>/operon_labels.tsv` | reemplazar GO |
| Clasificador de proteínas | protein-classifier | clase/proceso/probabilidad | duplicar GO/ortología |
| Sigma promotores | sigma | `sigma_promoters.tsv` | usar pockets como prueba |
| Pockets sigma | sigma-pockets | pocket/score/residuos | declarar sigma confirmado |
| Circuitos | regulation/circuits | aristas TF→gen | inferir red desde GO sola |

## 5. Insumos compartidos

### Ya disponibles

```text
data/regulation/<g>/genes.tsv
data/regulation/<g>/operons.tsv
data/sigma/<g>/sigma_factors.tsv
data/sigma/<g>/sigma_promoters.tsv
data/labels/go_mf_cc/<g>/gene_go.tsv
data/labels/go_bp/<g>/gene_bp.tsv
data/labels/ncrna_go/<g>/ncrna_go.tsv
data/labels/external/go/go.obo
data/labels/external/go/rfam2go
data/labels/external/goa/by_code/*.gaf.gz
data/labels/go_bp_direct/goa_bp.sqlite
```

### En construcción

```text
data/labels/go_indirect/gene_go_homology.tsv
data/labels/go_indirect/gene_go_domains.tsv
data/labels/go_indirect/gene_go_pathway.tsv
data/labels/gene_evidence.tsv          ← contrato nuevo
```

## 6. Reglas de no duplicación

```text
1. GO directo: labels/coverage.
2. GO no directo: labels/coverage.
3. Ortología: labels/coverage produce; protein-classifier consume.
4. Dominios: labels/coverage produce; protein-classifier consume.
5. Pockets: sigma-pockets produce; sigma consume como prior.
6. Sigma: sigma produce; circuits consume.
7. No usar el gold de sigma como feature del clasificador.
```

Si un agente necesita una capa que ya existe, la consume; no la recalcula.

## 7. Contrato para el clasificador de proteínas

Entrada:

```text
gene_evidence.tsv
```

Salida:

```text
protein_class
process_association
pathway_association
probability
source
evidence
```

Condiciones:

- no entrenar con el gold de sigma;
- separar train/test por genoma o filo;
- declarar clases ausentes;
- no convertir probabilidad en confirmación.

## 8. Contrato para pockets de sigma

Entrada:

```text
sigma_id
estructura
pocket
```

Salida:

```text
sigma_id
pocket_id
p2rank_score
residues
plddt
ptm
structure_origin
```

Condiciones:

- pocket no es prueba de promotor;
- no usar el mismo complejo para entrenar y evaluar;
- declarar si la estructura es AFDB/ESMFold/PDB/Kaggle.

## 9. Preguntas abiertas

```text
[ ] ¿Quién produce gene_evidence.tsv?
[ ] ¿Qué versión de GO/GOA se congela?
[ ] ¿Qué umbrales de homología se usan para ISS?
[ ] ¿Cómo se separan train/test por filo?
[ ] ¿Qué agente corre la validación final?
```

## 10. Bitácora

Formato:

```text
### 2026-01-01 | agente | TIPO | capa | estado
Mensaje breve.
→ archivos/decisiones
```

### 2026-11-26 | lead | INFO | labels/coverage | activo

Estado actual:

```text
GO directo: 280.450 genes
BP directo: 230.689 genes
ARN: capa Rfam2GO independiente
No directo: pipeline preparada, falta ejecución por lotes
```

Archivos:

```text
docs/GO_MF_CC.md
docs/GO_BP.md
docs/NCRNA_GO.md
docs/GO_INDIRECT.md
```

Pendiente:

```text
homología por lotes
InterPro/Pfam
ruta completa IGC
gene_evidence.tsv
```

### 2026-11-26 | lead | REQUEST | protein-classifier | abierto

Necesito que el clasificador consuma `gene_evidence.tsv` y no reanote GO/ortología.  
Entregar: columnas de salida y criterio de train/test.

### 2026-11-26 | lead | REQUEST | sigma-pockets | abierto

Necesito el esquema de pockets por sigma y la distinción estructura experimental/predicha.  
Entregar: `sigma_id`, `pocket_id`, score, residuos, pLDDT/pTM, origen.

### 2026-11-26 | lead | DECISION | arquitectura | acordado

```text
GO/etiquetas no son prueba de sigma.
Pockets no son prueba de sigma.
El clasificador no usa el gold de sigma.
Cada fila declara source/evidence/with_from.
```

## 11. Enlaces

- Repositorio: https://github.com/Krim10000/Bacterial-DNA-Annotator
- Release: https://github.com/Krim10000/Bacterial-DNA-Annotator/tree/v1.0
- Este archivo: https://github.com/Krim10000/Bacterial-DNA-Annotator/blob/v1.0/docs/AGENT_BUS.md
