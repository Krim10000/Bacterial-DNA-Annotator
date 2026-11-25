# Pipeline para familias uncharacterized/DUF

Diagrama: `docs/PIPELINE_DIAGRAM.png`

## Punto de partida

- 1.337.745 proteínas.
- 690.157 con función nombrada.
- 647.588 sin hit Swiss-Prot.
- 54.031 uncharacterized/DUF con UniProt.
- 25.058 clusters al 50%.
- 3.931 accesiones únicas.

## Fase actual

- 1.192 con AFDB.
- 2.570 cortas en ESMFold.
- 169 largas aparte.
- 587 cortas completadas.
- QC pLDDT/pTM.

## Después de las estructuras

### 1. Filtros estructurales

- F1 Foldseek.
- F2 EnzyMM/M-CSA.
- F3 P2Rank.
- F4 AlphaFill.
- F5 solo mapa propio.
- Propagación a clusters.

### 2. Familias conservadas

- Clusterizar por secuencia (30–90%).
- Clusterizar por estructura (Foldseek all-vs-all).
- Añadir dominio/familia.
- Elegir representante.
- Asignar `FAM_xxxx`.

No se llama función: se define familia.

### 3. Historia evolutiva con 400 especies

1. Orthogroups con MMseqs2/OrthoFinder.
2. Árbol de especies (GTDB/16S).
3. Árbol de familia.
4. Presencia/ausencia por especie.
5. Ganancia/pérdida/duplicación.
6. Sintenía y conservación de vecindad.
7. Filotipo: core, variable, linaje-específico.

Cuidado:

- no inferir HGT solo por topología;
- no inferir función por conservación;
- la historia explica distribución, no bioquímica.

### 4. Definición negativa

No se puede probar “no es enzima” por ausencia.

Se registra:

```text
no evidencia de dominio catalítico
no evidencia de TM
no evidencia de unión a ADN
no evidencia de homólogo funcional
```

También:

- calidad del ORF;
- start/RBS;
- longitud;
- conservación;
- solapamiento;
- bajo pLDDT.

Esto elimina ORFs dudosos, no proteínas reales.

### 5. Vecindad génica

- operón;
- genes conservados alrededor;
- coocurrencia en 400 genomas;
- sintenía;
- genes conocidos de la misma ruta;
- sistemas de defensa/móviles.

La vecindad es hipótesis, no prueba.

### 6. Complejos

Solo después de todos los monómeros:

- pares curados por operón;
- AlphaFold-Multimer/ColabFold;
- ipTM/PAE;
- no correr 1.823 pares;
- no confirmar interacción in vivo.

### 7. Catálogo

Una fila por familia:

```text
FAM_id
representante
miembros
clusters
taxonomía
dominio
estructura
F1-F5
contexto
posible complejo
confianza
límites
```

## Regla

- `FAM` no es función.
- `DUF` no es enzima.
- `no evidencia` no es `no existe`.
- `complejo predicho` no es interacción confirmada.
