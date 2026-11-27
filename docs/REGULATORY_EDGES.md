# Aristas TF–promotor

## Fuentes descargadas

### RegulonDB K-12

```text
data/regulondb/TF-RISet.tsv
5.785 líneas
```

### PRODORIC

```text
409 matrices
3.215 sitios de unión
2.076 sitios mapeados a nuestras especies
```

Mapeo:

```text
ecoli_k12       1.071
subtilis          687
paeruginosa       318
```

## Tabla común

```text
data/regulatory/evidence_edges.tsv
7.861 aristas
```

Columnas:

```text
genome
source_db
source_version
source_organism
tf_name
matrix_id
target
target_type
regulation
binding_sequence
evidence
confidence
match_type
source_url
notes
```

## Reglas

- `match_type = curated_same_species` para RegulonDB/PRODORIC directos.
- No hay inferencia por ortología en esta tabla.
- `confidence` no es probabilidad.
- No mezclar fuente curada con inferencia.

## Pipeline

```bash
python scripts/regulatory/01_fetch_prodoric.py
python scripts/regulatory/02_build_edges.py
```

## GUI

Modo:

```text
🧩 Regulación
```

Filtros por genoma, fuente y confianza; descarga TSV.

## Bases no disponibles

```text
CollecTF     → el endpoint público devuelve 404 / app sin API clara
RegPrecise   → sin volcado directo localizado
P2TF         → el sitio devuelve una app que no expone el volcado
```

No se inventaron datos.

## Compresión

```text
data/regulatory/external/prodoric/matrices.tar.gz
data/regulatory/external/prodoric/matrices.json.gz
```

Las matrices crudas ya no se usan para la tabla final.

## ENTRAF

Catálogo de TFs añadido:

```text
668 entradas
361 mapeadas
307 no_match
```

Ver `docs/ENTRAF.md` y la pestaña `Catálogo ENTRAF` en la GUI.
