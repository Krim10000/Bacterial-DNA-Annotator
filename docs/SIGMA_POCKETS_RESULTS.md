# P2Rank sobre sigma: resultados

## Corrida

```text
sigma con estructura AFDB: 169
P2Rank 2.5.1
2 threads
nice 19
```

## Resultado

```text
sigma con al menos un pocket: 125 / 169
sigma sin pocket:              44
pockets:                       833
```

Salida:

```text
data/sigma/structures/sigma_pockets.tsv
data/sigma/structures/p2rank_sigma_summary.json
data/sigma/structures/p2rank_sigma/<sigma_id>/
```

Columnas:

```text
sigma_id
genome
gene_id
accession
structure_origin = AFDB
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

## Lectura

- 125 sigma tienen cavidad predicha.
- 44 no tienen pocket: `no_data`.
- El `plddt` es del modelo AFDB, no del pocket.
- No hay pTM ni volumen Å³ en esta salida.
- No hay complejo sigma–ADN.

## Regla

```text
pocket del sigma ≠ promotor
pocket del sigma ≠ sigma confirmado
```

Es un prior estructural para el clasificador/sigma, no una prueba.

## Pendiente

- mapear tipo de sigma en la tabla de pockets;
- integrar como feature;
- validar contra sitios conocidos;
- no usarlo como gold.
