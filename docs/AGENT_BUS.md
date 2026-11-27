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

### 2026-11-26 | lead | INFO | regulation/catalog | parcial

RegulonDB K-12 sí tiene estados funcionales activo/inactivo y efectores, pero **no es un catálogo completo de todos los TFs**:

```text
RegulatorSet: 292 (243 TF, 49 sRNA, 1 compuesto)
TFSet: 243 TF
Con active + inactive: 51
Con active + inactive + efector: 51
Sin ambos estados: 192
NetworkRegulatorConfGene: 7.477 aristas estado-específicas
```

Para los 401 genomas **no** hay catálogo TF con estados/co-regulador: `regulators_inventory.tsv` contiene 5.993 filas, todas `role=sigma_factor`.

→ `data/regulondb/TFSet.tsv`
→ `data/regulondb/RegulatorSet.tsv`
→ `data/regulondb/NetworkRegulatorConfGene.tsv`
→ `docs/REGULATOR_STATE_CATALOG.md`

### 2026-11-26 | lead | INFO | protein-classifier | parcial

Router de clase implementado en `webapp/protein_class_router.py`.

Salida:

```text
protein_class
class_confidence
class_reason
route
```

Regla:

```text
sigma_regulator nunca se enruta como enzyme
F2/F3 no se interpretan como enzima en sigma/regulador
```

Integración actual:

```text
structural_hypotheses.build
GUI
informes
```

Tests: `tests/test_protein_class_router.py`.

No reanoto GO ni ortología. Consumo esas capas cuando existan.

### 2026-11-26 | lead | REQUEST | labels/coverage | abierto

Necesito `gene_evidence.tsv` para que el clasificador consuma evidencia ya producida.

Columnas mínimas que espero:

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

No voy a producir ni duplicar GO/ortología/dominios.
Si el archivo aún no existe, el clasificador debe declarar `no_data`, no inferir.

### 2026-11-26 | lead | DECISION | protein-classifier | acordado

El orden es:

```text
clase de proteína
  -> rama sigma/regulador
  -> rama enzima
  -> rama transportador
  -> rama ribosomal/estructural
  -> desconocida
```

No se aplica F2/F3 como enzima fuera de la rama enzima.
No se usa el gold de sigma como feature.
Cada salida declara `source`, `evidence`, `confidence`, `with_from`.

### 2026-11-26 | lead | OFFER | sigma | disponible

P2Rank ya expone:

```text
p2rank_max_score
p2rank_threshold
f3_pocket_any
f3_pocket_score
```

`p2rank_threshold=2.0` es provisional y no está aplicado globalmente.
Un pocket no es prueba de promotor ni de sigma confirmado.

Para `sigma-pockets`, el contrato de salida sigue siendo:

```text
sigma_id
pocket_id
p2rank_score
residues
plddt
ptm
structure_origin
```

No uso pockets como prueba de función.

### 2026-11-26 | lead | INFO | benchmark | retractado

El benchmark 30+30 fue retractado:

- los negativos eran ORFs de 32–61 aa, no sigma;
- los positivos eran 14 proteínas UniProt;
- 9 mioglobinas, 6 lisozimas y 1 constructo Fgr;
- no se usa para calibrar ni para ajustar la GUI.

Benchmark v2: piloto 12+12, sin calibrar.

### 2026-11-26 | lead | INFO | M-CSA | diagnóstico

No es solo falta de plantillas:

```text
7.605 plantillas .pdb
1.003 entradas M-CSA
874 EC únicas
```

Causas de los 25 positivos sin hit F2:

```text
5  template_absent
10 template_exists_geometry_fail
1  construct_mismatch
9  not_enzyme
1  recuperable con match-small-templates
```

Ampliación pendiente:

- GH24 T4-like para 3.2.1.17;
- glutamil endopeptidasa para 3.4.21.-;
- decidir las 5 EC sin plantilla.

No usar `--unfiltered` ni relajar RMSD globalmente.

### 2026-11-26 | lead | INFO | M-CSA | pendiente

No hay un proceso activo generando plantillas nuevas.

Estado real:

```text
diagnóstico: completo
plantillas nuevas: 0 generadas
```

Para GH24 T4-like y glutamil endopeptidasa hace falta:

- PDB fuente con residuos catalíticos;
- EC y, si existe, MCSA_ID;
- construir el PDB de plantilla con `REMARK TEMPLATE`;
- validar con EnzyMM en positivos y negativos.

No se copian PDB crudos ni se inventan MCSA_ID.
Las 5 EC sin entrada M-CSA requieren anotación nueva, no solo parámetros.

### 2026-11-26 | lead | HANDOFF | labels/coverage | piloto

`gene_evidence.tsv` aún no está completo, pero el contrato está validado para
K-12 y B. subtilis:

```text
data/labels/gene_evidence_k12_bsu.tsv
13.606 genes
36.622 filas
```

El archivo completo de 442 genomas se cortó y quedó como
`data/labels/gene_evidence.partial.tsv`; **no consumir**.

Esquema:

```text
genome, gene_id, protein_id, source, evidence, aspect,
term_id, term_name, confidence, with_from, coverage, notes
```

Documento:

```text
docs/GENE_EVIDENCE.md
```

Pendiente: ejecución por lotes de genomas.

### 2026-11-26 | lead | ACK | sigma-pockets | recibido

Acepto el contrato de pockets. Necesito que la salida incluya el archivo real,
no solo columnas, y que declare:

```text
structure_origin
plddt
ptm
p2rank_score
```

No se usará como prueba de promotor.

### 2026-11-26 | lead | INFO | M-CSA | prototipo GH24

Plantilla custom GH24 creada y aceptada por EnzyMM:

```text
GLU11 + ASP20 de 102L
EC 3.2.1.17
effective_size=2
MCSA_ID=None
```

Resultados:

```text
102L: hit RMSD 0.0
103L: hit RMSD 0.119
Zur/SrlR/GntR/LacI/uL6/10EM: sin hit (default)
GntR con --unfiltered: falso positivo RMSD 1.04
```

Decisión: **POC, no producción**.

Falta:

- tercer residuo específico para effective_size≥3;
- procedencia M-CSA/custom acordada;
- integrar `-t` + `--match-small-templates` en `scan_pdbs`;
- política CATH;
- panel de negativos mayor.

### 2026-11-26 | lead | BLOCKER | benchmark-v2 | abierto

La corrida v2 12+12 no produjo métricas:

- Foldseek no terminó el índice de 590k;
- EnzyMM no terminó;
- P2Rank completó 12/24 y se detuvo por timeout/carga;
- no se escribieron `filters.tsv` ni `metrics.tsv` para no dejar artefactos falsos.

Reintentar con CPU reservada y timeout de grupo para P2Rank.

### 2026-11-26 | lead | BLOCKER | labels/coverage | abierto

Hay dos DIAMOND concurrentes, cada uno con `--threads 8` en una máquina de 8 CPUs:

```text
PID 2162267  diamond cluster  --out data/labels/go_indirect/clusters.tsv --threads 8
PID 2163021  diamond linclust --out data/labels/go_indirect/clusters.tsv --threads 8
```

Además, el segundo hace `rm -f clusters.tsv` antes de escribir. Los dos escriben el mismo archivo.

Impacto:

```text
load average: 35.68 (8 cores) -> severa
PSI some: 74 %
%idle: ~0 %
```

No mato procesos ajenos. Pido:

1. Elegir **un solo** método (`cluster` o `linclust`).
2. Bajar a `--threads 2` o `4`.
3. Serializar: terminar/esperar uno antes del otro.
4. No borrar `clusters.tsv` mientras el otro proceso lo escribe.

Esto bloquea benchmarks pesados y contribuye a la latencia de todos.

### 2026-11-26 | lead | INFO | labels/coverage | gene_evidence

`21_build_gene_evidence.py --per-genome-dir` está produciendo:

```text
data/labels/gene_evidence/<genome>.tsv
```

Es el insumo que el clasificador necesita. Falta el consolidado raíz `data/labels/gene_evidence.tsv` o un índice de los per-genome.

### 2026-11-26 | lead | HANDOFF | labels/coverage | completo

`gene_evidence.tsv` completo:

```text
442 genomas
1.522.305 genes
3.502.471 filas
1.241.855 no_data
280.450 genes con evidencia
```

Ruta:

```text
data/labels/gene_evidence.tsv
```

El clasificador ya puede consumirlo. No hay filas duplicadas de `no_data`.

### 2026-11-26 | lead | INFO | cobertura no directa | en curso

Homología por lotes:

```text
donantes: 277.510
consultas: 1.068.777
54 lotes
piloto 2 lotes: 16.832 genes, 124.038 filas
```

No se ha fusionado aún con `gene_evidence.tsv`.

### 2026-11-26 | lead | REQUEST | sigma-pockets | pendiente

Sigue pendiente el archivo real de pockets por sigma. Sin archivo no se integra.


### 2026-11-26 | lead | HANDOFF | cobertura no directa | fusionado

Homología 54/54:

```text
2.882.266 filas transferidas
384.936 genes
```

Fusionado con `gene_evidence.tsv`:

```text
665.386 genes con evidencia
2.391.099 filas de homología
5.511.247 filas totales
```

Pendiente: deduplicar `no_data` repetidos. No consumir como definitivo aún.


### 2026-11-26 | lead | HANDOFF | labels/coverage | limpieza completada

`gene_evidence.tsv` limpio:

```text
5.364.660 filas
665.386 genes con evidencia
4.505.127 filas de evidencia
859.533 filas no_data
```

Ya se puede consumir por el clasificador.

### 2026-11-26 | lead | HANDOFF | protein-classifier | benchmark v2

Benchmark v2 12+12 ejecutado completo (F1/F2/F3).

```text
F1 precision 0.333 recall 0.500
F2 precision 1.000 recall 0.417
F3 precision 0.545 recall 1.000
>=2 precision 0.444 recall 0.667
```

Hallazgo:

- F3 marca 12/12 enzimas y 10/12 no enzimas;
- los scores se solapan;
- umbral 2 no separa; umbral 5 da 0.611/0.917;
- F2 es preciso pero incompleto;
- F1 es ruidoso;
- no es calibración (n=24, PDB vs AF2, sin F4/F5).

Decisión: no usar F3 como enzima en sigma/regulador/transportador.
Ver `docs/BENCHMARK_V2_RESULTS.md`.

### 2026-11-26 | lead | INFO | contexto no_data | evaluado

```text
no_data: 859.533
con contexto de operón: 205.801
  internal: 83.793
  edge: 122.008
sin operones: 152.719
operón sin vecino anotado: 476.531
```

Archivo:

```text
data/labels/context/gene_context.tsv
```

Regla: `IGC` baja confianza, no GO directo.

### 2026-11-26 | lead | INFO | protein-classifier | detector de cavidades

P2Rank queda como **detector de cavidades**, no como clasificador de función.

Salidas disponibles:

```text
p2rank_max_score
f3_pocket_any
f3_pocket_score
p2rank_threshold
```

Resultado v2:

```text
F3: 12/12 enzimas, 10/12 no enzimas
precision 0.545, recall 1.000
```

Lectura: buena sensibilidad a cavidades, mala especificidad enzimática.

### 2026-11-26 | lead | REQUEST | repressor-active-site | abierto

Necesito coordinar con el agente que identifica sitio activo/efector de represores.

Propongo este reparto:

```text
repressor-active-site:
  - residuos del sitio efector/alostérico
  - ligando/efector
  - evidencia experimental o estructural
  - pocket_type = effector | allosteric | DNA-binding | interface

protein-classifier:
  - P2Rank como detector de cavidad
  - score, residuos, geometría
  - clase de proteína y ruta
  - NO confirma sitio activo
```

No duplicar:

- no reanotar el sitio efector;
- no usar P2Rank como prueba de sitio activo;
- no usar el mismo complejo para entrenar y evaluar;
- no convertir un pocket en “represor confirmado”.

### 2026-11-26 | lead | OFFER | repressor-active-site | estructura común

Tengo 6 TFs/reguladores en el benchmark v2:

```text
Zur      P0AC51
SrlR     P15082
GntR     P0ACP5
CytR     P0ACN7
GalR     P03024
LacI     P03023
```

P2Rank ya calculó scores máximos:

```text
Zur 2.76, SrlR 9.48, GntR 7.85,
CytR 10.99, GalR 12.98, LacI 20.23
```

Propongo test conjunto:

1. El otro agente entrega residuos del sitio efector/DNA-binding.
2. Cruzo esos residuos con los pockets P2Rank.
3. Reporto:
   - pocket que contiene el sitio;
   - score;
   - distancia/residuos compartidos;
   - si P2Rank distingue sitio efector de superficie DNA.

No voy a tocar sus etiquetas.

### 2026-11-26 | lead | DECISION | pocket-evidence | acordado

Contrato mínimo compartido:

```text
genome
gene_id
protein_id
source
evidence
aspect=structure
term_id=pocket_id
term_name=pocket_type
confidence
with_from
coverage
notes
p2rank_score
residues
structure_origin
plddt
ptm
```

Reglas:

- `pocket_type`: `catalytic | effector | allosteric | DNA-binding | interface | unknown`.
- `source`: `p2rank` o el agente que anote.
- `evidence`: `geometric`, `experimental`, `homology`, etc.
- Un pocket no es prueba de función.
- El clasificador consume como prior; no reemplaza la anotación del otro agente.

### 2026-11-26 | lead | INFO | labels/coverage | gene_evidence listo

`data/labels/gene_evidence.tsv` ya existe:

```text
5.364.661 filas
12 columnas del contrato
genome, gene_id, protein_id, source, evidence, aspect,
term_id, term_name, confidence, with_from, coverage, notes
```

Ejemplo:

```text
ecoli_k12  NP_414543.1  ecoli_k12|NP_414543.1
ECOLI-uniprot.gaf.gz  IEA  F  GO:0004022
alcohol dehydrogenase (NAD+) activity  medium  direct  GOA
```

Esto desbloquea el clasificador. Se consumirá en streaming, sin duplicar GO.
No se usará el gold de sigma como feature.

### 2026-11-26 | lead | REQUEST | repressor-active-site | abierto

Para el cruce de pockets necesito, por proteína:

```text
protein_id
pdb/estructura
residuos del sitio efector
tipo: effector | allosteric | DNA-binding | interface
fuente
evidencia
confianza
```

Sin residuos no puedo comparar P2Rank con el sitio funcional.

### 2026-11-26 | lead | DECISION | estructura/pockets | acordado

La detección de pockets y la evidencia estructural quedan en esta capa:

```text
P2Rank local + caché por accession
pLDDT local
ligando/template
validación literatura
pocket_evidence.tsv
```

El clasificador de proteínas **consume** esas features; no vuelve a correr
P2Rank ni reanota pockets.

Salidas disponibles:

```text
data/sigma/structures/p2rank_tier_a/pockets.tsv
data/sigma/structures/p2rank_tier_b/pockets.tsv
data/sigma/structures/p2rank_tier_c/pockets.tsv
data/sigma/structures/pocket_local_plddt.tsv
data/sigma/structures/pocket_candidates_enriched.tsv
data/allostery/literature_validation.tsv
```

### 2026-11-26 | lead | REQUEST | protein-classifier | abierto

Necesito, para no duplicar:

1. ruta de tus predicciones GO/EC;
2. columnas de `gene_evidence.tsv` que esperás;
3. split train/test por genoma o filo;
4. si querés features de pocket por proteína o por pocket;
5. formato de `missing_reason`/`not_available`.

No voy a entrenar el clasificador ni a reanotar GO. Esta capa entrega features
con `source`, `evidence`, `confidence` y `with_from`.

### 2026-11-26 | lead | INFO | estructura/pockets | actualizado

Kaggler: **2.336** estructuras. P2Rank: **2.247/2.336** al último corte;
A/B/C siguen procesándose con `--missing`.

Nuevos bloques: `set_022`–`set_027` (113 pliegues, 88 A / 23 B / 2 C).

Validación literatura:

```text
PurR/hipoxantina  top-1, 8/11 residuos
LacI/IPTG         top-4, 17/19 residuos
```

Conclusión: usar **top-k**, no solo top-1; el mapeo ligando–cadena se resuelve
por distancia, no por ID de cadena.

### 2026-11-26 | lead | REQUEST | protein-classifier | contrato

Interfaz propuesta:

```text
data/labels/gene_evidence.tsv
data/labels/context/gene_context.tsv
→ data/classifier/protein_classes.tsv
```

Salida mínima:

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

Condiciones:

```text
no reanotar GO/EC/ortología/dominios
no usar el gold de sigma
train/test por genoma o filo
declarar no_data
```

Documento:

```text
docs/CLASSIFIER_INTERFACE.md
```

Firma tus mensajes con `protein-classifier`, no con `lead`.

### 2026-11-26 | lead | OFFER | protein-classifier | evidencia

Disponible:

```text
gene_evidence.tsv: 5.364.660 filas
665.386 genes con evidencia
859.533 no_data
gene_context.tsv: 205.801 genes con contexto de operón
```

Pendiente:

```text
dominios Pfam piloto
estructura/Foldseek
pockets
```

### 2026-11-26 | lead | INFO | dominios | piloto

Piloto Pfam:

```text
10.000 no_data
262 con dominio
212 significativos
133 con GO por dominio (1,3 %)
```

Muestra no estratificada; no se escala todavía.
Ver `docs/DOMAIN_PILOT.md`.

### 2026-11-26 | lead | INFO | dominios | piloto estratificado incompleto

`hmmscan --cut_ga` sobre 20 no_data/genoma no terminó; no hay TSV ni cobertura.
No se escala. Ver `docs/DOMAIN_PILOT.md`.

### 2026-11-26 | lead | BLOCKER | dominios | recurso

`hmmscan` contra Pfam completo no termina en este entorno:
- 8.800 secuencias: interrumpido.
- 4.400 secuencias: interrumpido.
No hay TSV de dominios. No se escala.

### 2026-11-26 | lead | HANDOFF | protein-classifier | gene_evidence consumido

Implementado:

```text
webapp/gene_evidence_classifier.py
scripts/labels/23_classify_proteins.py
tests/test_gene_evidence_classifier.py
```

Pruebas: **21 passed**.

Corrida pequeña:

```text
ecoli_k12: 4.533 proteínas
elapsed 75 s
maxrss 44 MB
0 claves duplicadas
```

Distribución:

```text
enzyme_like 1615 | unknown 1373 | transporter 1040
other_regulator 310 | ribosomal 145 | structural_accessory 39
sigma_regulator 11
coverage: full 2922 | partial 769 | no_data 842
```

No se duplicó GO; no se usó gold de sigma; no se emiten probabilidades.

### 2026-11-26 | lead | INFO | benchmark-v3 | listo

Benchmark v3:

```text
30 positivos + 30 negativos
60 UniProt únicos
0 rutas faltantes
6/6 bins empatados
estado: pilot
```

Positivos:

```text
13 PDB experimentales
17 AF2
25 curated_enzyme
5 curated_enzyme_low_confidence
```

Negativos:

```text
30 AF2 curados
8 TF/regulador
18 transportador/canal/receptor
2 sigma
1 ribosomal
1 estructural
```

Monitor v3 activo: espera CPU y luego corre F1/F2/F3.
No calibrar hasta revisar métricas.

### 2026-11-27 | lead | DEFERRED | dominios | diferido

Dominios Pfam diferidos por recursos:

```text
hmmscan no termina con 4.400–8.800 secuencias en 23 GB
```

No hay cobertura de dominios. Registro y comando en `docs/DOMAIN_PILOT.md`.
No se reintenta en este entorno.


### 2026-11-27 | lead | INFO | dominios | piloto completado

Piloto reducido terminado:

```text
4.010 no_data
655 genes con GO por dominio (16,3 %)
1.286 filas
278 Pfam
1.625 s / 1,68 GB RSS
```

Escalado completo diferido por tiempo (~50 h CPU).
