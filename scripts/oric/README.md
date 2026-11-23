# Flujo de Trabajo: Orígenes de Replicación (oriC)

Este módulo se encarga de predecir la posición del origen de replicación (oriC) 
en genomas bacterianos, resolviendo los problemas de ambigüedad de marcadores.

## Prerrequisitos
- Tener el entorno virtual activo.
- Haber corrido el pipeline de anotación (`scripts/annotation/`) para tener los GFF y FASTA crudos.
- Haber descargado y descomprimido la base de datos DoriC en `data/doric/`.

## Pasos del Flujo

### Paso 1: Validación de Marcadores (`01_validate_markers.py`)
Busca los genes `dnaA` y `dnaN` en las anotaciones de NCBI. Verifica unicidad, 
posición (no pegados a los bordes) y coherencia del operón. Resuelve casos 
ambiguos generando un reporte de candidatos.
* **Uso:** `python scripts/oric/01_validate_markers.py`
* **Salida:** `data/marker_validation/validated_genomes.txt`

### Paso 2: Cruce con DoriC (`02_doric_lookup.py`)
Lee la base de datos DoriC (curada experimentalmente) y extrae la posición 
real del oriC para los genomas que la tengan anotada. Esto sirve como 
"Ground Truth" para medir la precisión de ORCA.
* **Uso:** `python scripts/oric/02_doric_lookup.py`
* **Salida:** `data/doric_ground_truth.json`

### Paso 3: Predicción con ORCA (`03_orca_predict.py`)
Ejecuta ORCA (Machine Learning) sobre los genomas validados. 
Si el genoma estaba en DoriC, compara la predicción contra la verdad absoluta. 
Si era ambiguo, usa DoriC o heurísticas para elegir el mejor candidato.
* **Uso:** `python scripts/oric/03_orca_predict.py`
* **Salida:** `data/orca_predictions/{genoma}.json`

Ver estado, límites y plan de integración en
[`docs/INTEGRATIONS.md`](../../docs/INTEGRATIONS.md).
