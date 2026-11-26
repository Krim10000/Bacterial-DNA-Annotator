# ARN y oriC en el navegador

## tRNA, rRNA y ncRNA

Al seleccionar un gen de tipo ARN:

- se muestra la secuencia en orientación 5′→3′;
- se descarga como FASTA de ARN (`.rna.fna`);
- **no se traduce** a proteína;
- no se muestra un visor 3D vacío ni se inventa una estructura.

No hay modelos 3D de ARN locales en el proyecto. La ficha lo declara.

## Horquillas

Al seleccionar una horquilla:

- coordenadas de la región;
- hebra;
- clase/confianza;
- ΔG, tallo, bucle, poli-U y genes vecinos;
- secuencia de ARN (T→U), si el genoma está disponible;
- aviso de que el plegamiento no se dibuja.

No se dibuja un tallo-bucle inventado: no hay estructura 3D de ARN local.

## oriC

`load_oric` conserva:

```text
pos       punto ORCA
start/end región DoriC si existe
source    ORCA/DoriC
confidence
span_source
```

La pista dibuja una **banda** para la región y una línea para el punto ORCA.
La ficha muestra:

- inicio, fin y longitud;
- punto ORCA;
- confianza;
- fuente;
- secuencia de la región, si existe;
- aviso de que oriC no tiene hebra codificante: la dirección no define un
  transcrito.

Si solo hay un punto, se declara `ORCA (punto)` y no se inventa un intervalo.

## Estructura secundaria integrada

La ficha de ARN y la ficha de horquilla llaman a `webapp/rna_structure.py`:

- ViennaRNA `RNA.fold` si está disponible;
- Nussinov como respaldo;
- `not_available` si no hay backend.

Se muestra el dot-bracket y, si existe, la energía MFE. La interfaz declara que
es una predicción, no una estructura 3D ni una validación experimental.

### Región DoriC vs punto: auditoría y límites

La evidencia se separa en dos piezas que no son intercambiables:

- **punto**: ORCA o el ancla elegida (DoriC/dnaA/dnaN);
- **región**: solo si DoriC publica un intervalo.

`webapp/oric_regions.py` normaliza la región a segmentos 1-based inclusivos,
resuelve cruces de origen (`end < start`) y el centinela `end == 0`, y extrae
la secuencia sin inventar endpoints. El índice compacto vive en
`data/oric_regions/oriC_regions.json`; el generador es
`scripts/oric/41_build_region_index.py` y el subconjunto DoriC se crea con
`scripts/oric/42_extract_doric_subset.py`.

Cobertura auditada de los 400 genomas:

| lectura | región DoriC | solo punto | sin lectura |
|---|---:|---:|---:|
| `load_oric` actual | 0 | 95 | 305 |
| índice de regiones | 321 | 79 | 0 |

El bloqueo actual es concreto: `load_oric` exige `chosen.pos`; solo los 95
JSON de `orca_predictions_fix` lo tienen y ninguno incluye `doric_truth`. Los
288 JSON de `orca_predictions` usan `orca_position` y se descartan. Por eso hoy
la interfaz solo muestra 95 puntos y ningún intervalo. El índice recupera 321
regiones: 67 de genomas con fix que estaban en punto, 239 de predicciones sin
fix y 15 de genomas sin ORCA que sí tienen DoriC.

Reglas que no se negocian:

- sin región DoriC, `start == end == pos` y `span_source = "ORCA (punto)"`;
- no se estima un ancho a partir del punto;
- `oriC_all.gff3` son ventanas de 21 pb (una de 17 pb por borde de contig), no
  regiones DoriC; no se usan como intervalo;
- los endpoints publicados no se reescriben: si la longitud declarada difiere
  en ±1 pb por convención 0/1-based, se reporta como `length_delta`;
- `doric_ground_truth.json` colapsa varias filas DoriC en una; el índice
  conserva la verdad canónica solo si coincide con una fila **cromosómica**
  local. En 5 genomas la GT apunta a una fila de plásmido: se descarta y se usa
  la fila cromosómica (o solo punto). Las alternativas quedan con su distancia.

Integración pendiente: `load_oric` puede devolver
`webapp.oric_regions.load_region(data_root, genome)`. Esta auditoría **no**
modificó `browser.py`; mientras no se haga ese cambio, la interfaz sigue
leyendo el JSON de predicción y muestra solo los 95 puntos.

Pruebas:

```bash
python -m pytest tests/test_oric_regions.py -q
```
