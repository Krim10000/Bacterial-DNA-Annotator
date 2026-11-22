# Datos del proyecto

Este directorio está **vacío en la versión publicada** a propósito.

Los resultados del proyecto (400 genomas anotados, redes metabólicas, tablas de
regulación y estructura) pesan más de 20 GB y, sobre todo, **incluyen datos de
terceros que no se pueden redistribuir**:

- RegulonDB (uso académico / no comercial),
- KEGG (suscripción, sin derecho de redistribución),
- DoriC (uso académico, licencia no explícita),
- STRING, DBTBS y BATTER (acuerdos académicos o sin licencia declarada).

Ver [`../DATA_NOTICE.md`](../DATA_NOTICE.md).

## Qué se puede hacer sin datos

- **Anotar FASTA**: funciona con los ejemplos de `webapp/samples/` o con tus
  propios genomas.
- **Navegador**: el genoma `demo` viene incluido en `webapp/demo/` (150 kb
  anotados desde cero). Es ilustrativo, no evidencia experimental.
- **Descargar** GFF3, proteínas, candidatos y desiertos.
- Instalar y probar la interfaz, el navegador y los créditos.

## Qué necesita datos

| modo | necesita |
|---|---|
| Rutas metabólicas | tablas en `data/metabolism/`, `data/regulation/`, etc. |
| Navegador de genoma | `data/regulation/<genoma>/genes.tsv` y pistas asociadas |
| Móviles, membrana, taxonomía, cadenas | sus carpetas en `data/` |

Si faltan, la app lo dice y muestra el comando. No inventa resultados.

## Cómo generarlos

```bash
# 1. Descargar las bases una vez
python -m webapp.setup_data --data-root data          # Swiss-Prot + Rhea
python -m webapp.setup_data --data-root data --pfam   # + Pfam (~600 MB)

# 2. Anotar los genomas del proyecto
python scripts/00_data_prep/00_download_catalog.py
python scripts/00_data_prep/04_download_fasta.py
python scripts/01_annotation/01_annotate_batch.py

# 3. Capas derivadas
bash run_pipeline.sh components
bash run_pipeline.sh rules
bash run_pipeline.sh hairpins
bash run_pipeline.sh closure
```

Cada etapa escribe su propia carpeta y es idempotente. Los tiempos y las
dependencias están en [`../scripts/README.md`](../scripts/README.md).

## Regenerar los datos de RegulonDB

```bash
python scripts/11_validation/regulondb_gt.py
```

Requiere aceptar las condiciones de RegulonDB y **no** se debe publicar el
resultado en un repositorio abierto. Es para uso local/revisión.

## Convención de carpetas

```
data/<módulo>/<archivo>.tsv
data/<módulo>/<genoma>.json
data/<módulo>/README.md      # si hay reglas propias
```

Nada de datos se escribe fuera de `data/`, `logs/`, `out*/` o `reports/`.
