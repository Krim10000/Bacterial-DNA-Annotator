# Red de control regulatorio

## Objetivo

Describir lazos de control en las 400 especies:

```text
TF / sigma
→ promotor / operón
→ gen(es)
→ proceso
```

Y detectar motivos:

```text
autoregulación
feed-forward
feedback
co-regulación
```

## Lo que ya existe

### K-12

```text
RegulonDB
243 TF
3.444 interacciones TF–promotor con secuencia
1.485 pares represor + evidencia experimental
51 TF con estado activo/inactivo + efector
```

### 400 genomas

```text
1.555 represores anotados
14.850 TF con familia declarada
sigma: 402 genomas con salida
operones/promotores: 400 genomas
GO/contexto: capas propias
```

Pero:

```text
no hay TF–promotor curado para los 400
no hay estados/efectores para los 400
```

## Esquema de aristas

```text
edge_id
genome
source_type      TF | sigma | operon | gene
source_id
target_type      operon | gene | sigma | TF
target_id
sign             + | - | ?
mechanism        binding | promoter | orthology | expression | context
evidence         EXP | IMP | IDA | ISS | ISO | IGC | IEA
source           RegulonDB | DBTBS | PWM | ortholog | expression | context
with_from
confidence
coverage
notes
```

## Fuentes de aristas

### 1. Sigma → operón

Ya existe:

```text
sigma_promoters.tsv
```

Evidencia:

```text
PWM
k-mer
inventario
contexto
```

### 2. TF → promotor

Jerarquía:

```text
1. RegulonDB/DBTBS experimental
2. PWM conocido + promotor
3. ortología de sitio/regulón
4. expresión/perturbación
5. contexto
```

### 3. TF → TF / sigma → sigma

Para lazos:

```text
TF regula sigma
sigma regula TF
TF regula TF
sigma regula sigma
```

## Signo

No inferir de GO.

Fuentes:

```text
rol curado (represor/activador)
RegulonDB/DBTBS
familia TF con rol conocido
expresión
```

Si no hay evidencia:

```text
sign = ?
```

## Detección de lazos

### Autoregulación

```text
TF → su propio operón
sigma → su propio operón
```

### Feed-forward

```text
TF → sigma
TF → operón
sigma → operón
```

### Feedback

```text
TF → ... → TF
sigma → ... → sigma
```

### Co-regulación

```text
TF A y TF B → mismo operón
sigma + TF → mismo promotor
```

## Fases

### Fase 1: K-12

```text
grafo curado
evaluar contra RegulonDB
línea base real
```

### Fase 2: B. subtilis

```text
DBTBS/BSGatlas
sigmulón derivado
grafo descriptivo
```

### Fase 3: 400

```text
transferencia por ortología
PWM/sitios conservados
sigma + TF + operón
```

### Fase 4: validación

```text
holdout por genoma/filo
precisión por tipo de arista
cobertura
signo
```

## Lo que no se debe hacer

```text
TF vecino ≠ TF que regula
pocket ≠ unión a ADN
GO de proceso ≠ dirección
correlación ≠ causalidad
```

## Métricas

```text
precisión por tipo de arista
recall de TF–promotor
exactitud de signo
número de lazos
cobertura por especie
```

## Prioridad inmediata

```text
1. Construir grafo K-12 desde RegulonDB.
2. Validar contra el mismo RegulonDB con holdout.
3. Añadir sigma–operón.
4. Añadir B. subtilis.
5. Transferir a 400 con ISS/IGC.
```

## Archivos propuestos

```text
data/regulatory/edges.tsv
data/regulatory/nodes.tsv
data/regulatory/loops.tsv
data/regulatory/summary.json
```
