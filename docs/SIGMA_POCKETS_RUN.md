# P2Rank sobre sigma

## Qué se está corriendo

```text
sigma con estructura local: 169
origen: AFDB
P2Rank 2.5.1
threads: 2
nice: 19
```

No se corren los 5.993 sigma: solo los que tienen estructura.

## Salida

```text
data/sigma/structures/p2rank_sigma/<sigma_id>/
data/sigma/structures/sigma_pockets.tsv
data/sigma/structures/p2rank_sigma_summary.json
```

Columnas:

```text
sigma_id
genome
gene_id
accession
structure_origin
structure_path
pocket_id
pocket_rank
score
probability
sas_points
surf_atoms
residues
plddt
```

## Regla

```text
pocket del sigma ≠ promotor
pocket del sigma ≠ sigma confirmado
```

Es un prior estructural, no una prueba.

## Estado

Corrida lanzada en background, reanudable: si existe la predicción, no se repite.

## Clasificador

Corrida completa de 442 terminada:

```text
1.524.918 filas
full:    557.763
partial: 107.623
no_data: 859.532
```
