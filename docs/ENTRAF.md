# ENTRAF

## Fuente

**ENTRAF**  
*An Encyclopedia of well-annotated DNA-binding TRanscription factors in Bacteria and Archaea.*

Computational Biology Lab  
IIMAS-UNAM Sede Mérida  
Km. 5.5 Carr. Sierra Papacal - Chuburna Pto.  
Tablaje 31257, Sierra Papacal, 97302 Mérida, Yucatán  
+52 (999) 406 00 03 ext. 7610  
ernesto.perez@iimas.unam.mx  
twitter.com/eperez09  
Computational Biology Lab on FB

URL:

http://web.pcyt.unam.mx/EntrafDB/

## Agradecimiento

Agradecemos a ENTRAF y al Computational Biology Lab, IIMAS-UNAM Sede Mérida,
por la curación y publicación de esta base.

## Cobertura

```text
668 entradas
668 con PDB
255 con PubMed
361 mapeadas a nuestras especies
307 no_match
```

Mapeo:

```text
ecoli_k12                        230
subtilis                         127
paeruginosa                        2
GCF_000091545.1_ASM9154v1          2
```

## Campos

```text
source_db
source_version
entraf_id
gene_name
uniprot_id
organism
sequence
length
function
role
family
pdb
pubmed
evidence
genome_match
match_type
hit_type
confidence
inference
source_url
citation
acknowledgement
retrieved_at
```

## Tabla

```text
data/regulatory/entraf_entries.tsv
```

`hit_type`:

```text
direct_species
no_match
```

`confidence`:

```text
high    PDB + PubMed
medium  PDB
low     sin PDB
```

`inference = none`: no se infirió ortología en esta tabla.

## Pipeline

```bash
python scripts/regulatory/03_scrape_entraf.py
python scripts/regulatory/04_build_entraf_table.py
```

## GUI

```text
🧩 Regulación → Catálogo ENTRAF
```

## Límite

`no_match` no significa que el TF no exista en esa especie; significa que el
organismo de ENTRAF no coincide con nuestro catálogo. Se puede usar después
para transferencia por ortología, marcada como inferencia.

## Compresión

```text
raw_html.tar.gz
```

El HTML crudo se conserva comprimido; la tabla activa queda sin comprimir.
