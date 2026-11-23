# Tres caminos de anotación de proteínas

## Estado de implementación (2026-11-23)

Ahora el pipeline de un FASTA nuevo:

- ejecuta DIAMOND/Swiss-Prot automáticamente;
- escribe `*.functional.tsv`, `*.functional_annotation.json` y
  `*.protein_paths.tsv`;
- clasifica cada proteína en A, B o C;
- conserva hit, identidad, e-value y producto.

Todavía no ejecuta automáticamente los pasos estructurales de B/C. Para B/C el
pipeline deja el trabajo preparado (`next_step`) y los filtros estructurales
siguen siendo opcionales (`--structure-manifest`, `--structural-filters`).
El camino C sigue siendo un plan CPU para los 400; no se debe presentar como
completo.

## Resumen

| camino | proteínas | % total | dato clave |
|---|---:|---:|---|
| A. Función nombrada | **690.157** | 51,6% | Swiss-Prot con función |
| B. UniProt uncharacterized/DUF | **54.031** | 4,0% | subconjunto de A; tiene accesión/AFDB |
| C. Sin hit Swiss-Prot | **647.588** | 48,4% | sin función curada; 446.927 clusters |

A + C = 1.337.745.  
B está contenido en A; se trata aparte por su falta de función específica.

## Camino A — función nombrada

- 690.157 proteínas.
- Hit Swiss-Prot con función.
- Se usa la anotación existente.
- No requiere estructura para una función básica.

## Camino B — uncharacterized/DUF con UniProt

- 54.031 proteínas.
- 25.058 clusters al 50%.
- 3.931 accesiones únicas.
- 1.192 con AFDB.
- 2.739 sin AFDB (2.570 ≤600 aa; 169 >600).
- Se aplican los 5 filtros y se propagan hipótesis.

Salidas:

- `UNCHAR_*`;
- `structural_hypotheses`;
- `CANDIDATES_508_REVIEW`.

## Camino C — sin hit Swiss-Prot

- 647.588 proteínas.
- 446.927 clusters al 50%.
- 360.741 clusters singleton (80,7% de clusters; 55,7% de proteínas).
- No hay AFDB directo para ORFs arbitrarios.
- No se pliegan todos.

### Singleton

Un cluster singleton tiene **una sola secuencia** en el conjunto a ese umbral.

No significa:

- que no tenga homólogo en UniRef/PDB;
- que sea una proteína rara;
- que no tenga dominio.

Sí significa:

- no se reduce por clustering dentro de este conjunto;
- hay que anotarlo por dominio/homología/contexto;
- puede ser un ORF dudoso y requerir control de calidad.

### Qué se hace en C

1. Reclustering/representantes.
2. MMseqs2/DIAMOND vs UniRef90/TrEMBL/PDB/eggNOG.
3. Pfam/InterPro en los que sigan sin señal.
4. TMHMM/Phobius/señal.
5. Contexto génico/operón.
6. Taxonomía/conservación.
7. Estructura solo para el top.

## Caminos en el pipeline

La documentación del pipeline debe distinguir:

```text
A: anotación existente
B: hipótesis estructurales sobre UniProt/AFDB
C: reducción CPU de no-hit + estructura selectiva
```

El camino B está implementado con `--structure-manifest`,
`--structural-filters` y `--precomputed-hypotheses`.

El camino C todavía es un plan CPU; no está automatizado por completo.

## Regla

- Un singleton no es “basura”.
- Un no-hit no es “sin homología”.
- Una hipótesis no es `confirmed`.
- La anotación completa significa evidencia completa, no función inventada.
