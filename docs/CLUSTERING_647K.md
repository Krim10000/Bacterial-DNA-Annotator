# Clustering de 647.588 proteínas sin hit

## Resultado

| métrica | valor |
|---|---:|
| Secuencias | 647.588 |
| Clusters al 50% | 446.927 |
| Reducción | 31,0% |
| Clusters con >1 | 86.186 |
| Clusters con >5 | 9.521 |
| Singletons | 360.741 |
| Tamaño máximo | 296 |

## Lectura

- El 50% deja muchísimos singletons.
- No basta para reducir a un conjunto manejable.
- Los 647k siguen sin estructura AFDB directa.

## Comparación

| conjunto | secuencias | clusters 50% | reducción |
|---|---:|---:|---:|
| Sin hit | 647.588 | 446.927 | 31% |
| Uncharacterized con UniProt | 54.031 | 25.058 | 54% |

## Implicación

- 446k clusters siguen siendo demasiados para los 5 filtros.
- Los filtros estructurales requieren PDB.
- Para 647k habría que plegar o usar otra señal.

## Ruta viable

1. Trabajar sobre los **54.031** con AFDB.
2. Esos ya son **25.058 clusters**.
3. Aplicar los 5 filtros a esos representantes.
4. Para 647k: clustering más agresivo, dominio/familia y selección estratificada.
