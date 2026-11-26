# Meta-clasificación de genes/operones para sigma

## Principio

Ninguna fuente basta. El objetivo no es una etiqueta única, sino una **red de
evidencia**:

```text
gen
→ reacción / EC
→ ruta / módulo
→ proceso / condición
→ regulón
→ sigma candidato
```

Cada arista tiene:

```text
fuente
tipo de evidencia
confianza
independencia
```

## Regla

Una reacción química no “es” de shock térmico. Puede estar **asociada** a una
condición en un contexto concreto.

Por eso:

```text
GO:0009408 response to heat
```

es una etiqueta de proceso, no una propiedad de la enzima.

## Ontología anidada

Usar DAG, no una lista plana:

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

Fuentes:

- GO + GO slims;
- KEGG modules/pathways;
- MetaCyc;
- SEED;
- eggNOG/COG;
- InterPro;
- EcoCyc/RegulonDB;
- SubtiWiki;
- MicroScope.

## Modelo de datos

### Etiqueta

```text
label_id
label_name
parent_id
level
ontology
```

### Gen-etiqueta

```text
gene_id
label_id
evidence_type: experimental | curated | inferred | predicted | text
source
confidence
```

### Operón-etiqueta

```text
operon_id
label_id
score
enrichment
background
n_genes
```

### Sigma-etiqueta

```text
sigma_id
label_id
direction: supported | absent | conflicting
weight
source
```

## Cómo se calcula

### 1. Gen

Cada gen recibe etiquetas de:

- EC/reacción;
- ruta;
- dominio;
- localización;
- GO;
- literatura.

### 2. Operón

Agregar sin doble conteo:

```text
score_operon(label) = f(genes del operón)
```

No sumar todos los genes si están duplicados o son parálogos.

### 3. Enriquecimiento

Comparar contra:

- genoma;
- background de operones;
- clado filogenético.

Salida:

```text
label
n_operon
n_genoma
OR / Fisher
FDR
```

### 4. Propagación jerárquica

Si un operón es “shock térmico”, también es “respuesta a estrés”.

La propagación debe tener decaimiento y no inflar todos los padres.

### 5. Sigma

El prior de sigma se construye desde:

```text
regulón conocido
+ promotores predichos
+ etiquetas de los genes diana
+ historia evolutiva
+ contexto
```

No desde una sola etiqueta.

## Historia evolutiva

El árbol aporta:

- orthogroup;
- ganancia/pérdida;
- expansión linaje-específica;
- conservación del regulón;
- clado para validación.

No aporta:

- función por sí solo;
- promotor;
- condición.

Uso correcto:

```text
sigma N1 pertenece al orthogroup OGxxxx
OGxxxx está conservado en el clado X
el clado X tiene genes de shock térmico
→ hipótesis, no prueba
```

## Reglas de evidencia

### Fuerte

- TSS dependiente de sigma;
- regulón curado;
- unión física;
- expresión condicional.

### Media

- promotor predicho;
- enriquecimiento de ruta;
- conservación de sintenia;
- ortología.

### Débil

- vecindad;
- co-ocurrencia;
- texto;
- dominio aislado.

### Prohibida

- usar la misma anotación como feature y gold;
- asumir que una ruta implica un sigma;
- etiquetar un gen por su vecino;
- propagar sin FDR;
- llamar “confirmado” a una hipótesis.

## Salida accionable

Para cada sigma:

```text
sigma_id
native_name
family
orthogroup
evidence:
  inventory
  promoter
  regulon
  labels
  evolution
  structure
hypothesis: shock térmico
confidence: media
rationale: ...
```

Para cada operón:

```text
operon_id
labels: [shock térmico, chaperonas, respuesta a estrés]
sigma_candidates: [sigma32, sigma70]
```

## Consultas

```python
labels_for_gene(gene_id)
labels_for_operon(operon_id)
enrichment_for_operon(operon_id)
labels_for_sigma(sigma_id)
evidence_path(sigma_id, label_id)
```

## Validación

- enriquecimiento contra background;
- FDR;
- ablación de etiquetas;
- provenance holdout;
- comparar K-12/B. subtilis;
- no usar el gold de sigma como etiqueta de proceso.

## Piloto

1. K-12 y B. subtilis.
2. GO/KEGG/EC.
3. Etiquetas de operón.
4. Enriquecimiento de regulones conocidos.
5. Ver si sigma32 recupera shock térmico, sigma54 nitrógeno, sigma28 flagelo.
6. Solo después, 400.

## Límites

- GO/KEGG tienen cobertura desigual;
- una reacción puede estar en varias condiciones;
- la anotación funcional puede ser incompleta;
- el contexto no prueba sigma;
- la literatura no tiene un benchmark estándar para esto.
