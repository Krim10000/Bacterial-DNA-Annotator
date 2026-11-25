# Capa TSS/TTS (`webapp/tss_layer.py`)

Capa de solo lectura para **cargar, normalizar y enlazar** inicios (TSS) y
finales (TTS) de transcripción ya calculados. No usa red, no usa GPU y no
predice posiciones: si no hay una tabla local para el genoma, devuelve
`not_available`.

## 1 · Qué es un TSS y un TTS

- **TSS** (*transcription start site*): coordenada 1-based del primer
  nucleótido del transcrito primario. Es el extremo 5' del ARN, no el sitio de
  unión del factor sigma ni la caja -10/-35.
- **TTS** (*transcription termination site*): coordenada del extremo 3' del
  transcrito. En las tablas de este repositorio suele ser el inicio de la
  horquilla de terminación (o el extremo izquierdo de un intervalo), no
  necesariamente el último nucleótido incorporado.

Un TSS no es automáticamente un promotor funcional: puede ser un inicio
leaderless, un promotor interno, un transcrito antisentido o un artefacto de
mapeo. Un TTS no siempre termina el transcrito del gen vecino: puede ser
antisentido o corresponder a un terminador dependiente de Rho. Por eso la capa
conserva `source`, `class`, `location` y `mechanism`, y no fusiona evidencias
distintas.

Convención de coordenadas: **1-based, inclusiva, por contig y hebra**. Una
posición 0 se rechaza en lugar de convertirla silenciosamente a 1: sumar 1 a un
BED sería inventar una coordenada.

## 2 · Qué hace y qué no hace

Hace:

1. Descubre tablas precalculadas para el genoma pedido.
2. Normaliza `strand`, `pos`, `source`, `class`, además de `location` y
   `mechanism`.
3. Detecta y elimina duplicados exactos.
4. Enlaza el TSS más cercano aguas arriba y el TTS más cercano aguas abajo de
   cada feature (misma hebra y contig, ventana configurable).
5. Resume por tipo, clase, ubicación, fuente, hebra y contig.

No hace:

- No predice TSS ni TTS a partir de secuencia.
- No transfiere coordenadas entre ensamblajes ni especies.
- No convierte promotores putativos de reglas (`promoters.tsv`) en TSS.
- No consulta DBTBS, RegulonDB, BSGatlas ni ninguna API.
- No ejecuta modelos, GPU ni entrenamiento.

## 3 · Fuentes y nivel de evidencia

| `source` canónico | Ejemplo local | Qué es | Límite |
|---|---|---|---|
| `DBTBS` | `data/dbtbs/dbtbs_tss.tsv` | TSS curados de *B. subtilis* con sigma y literatura | DBTBS no es descarga automatizable; puede estar incompleto |
| `BSGatlas` | `data/dbtbs/gt_tss.tsv`, `gt_terminators.tsv`, `subtilis_tss_tts.tsv` | Integración publicada que incluye DBTBS/BsubCyc y otros; TSS con sigma y terminadores con energía | No todos los terminadores son experimentales; muchos son predicciones estructurales |
| `RegulonDB` | `data/regulondb/tss_k12_union.tsv`, `TerminatorSet.tsv` | TSS curados y terminadores de *E. coli* K-12 | Es una especie; la unión mezcla fuentes y la tabla no trae contig |
| `RegulonDB alto rendimiento` | filas `HT xN` de la unión K-12 | TSS de estudios de alto rendimiento | Puede incluir ruido y promotores no canónicos |
| `transferido K-12` | `data/dbtbs/bsu_transferred_tss.tsv` | TSS transferidos por ortología desde K-12 | Inferencia, no medición en *B. subtilis*; confianza menor |
| `rule_based` | no se consume aquí | Promotores por PWM/reglas (`webapp/regulation_layer.py`) | Es una hipótesis declarada, no un TSS; no se carga como TSS |
| `precalculado` | tabla del usuario | Archivo local sin etiqueta de fuente reconocida | La procedencia real es responsabilidad de quien lo generó |

La evidencia no es equivalente entre filas. El código no la colapsa: dos
fuentes en la misma coordenada se conservan como corroboración.

## 4 · Archivos que busca

Orden de prioridad:

1. **Específico del genoma**:
   `data/regulation/<genoma>/`
   - `tss_tts.tsv`, `tss_tts.json`, `transcription.tsv` (tabla combinada), o
   - `tss.tsv` + `tts.tsv`, o
   - `terminators.tsv`, `tss.json`, `tts.json`, etc.
   - También mira **un** nivel de subdirectorio.
2. **Tabla combinada global**:
   `data/dbtbs/*_tss_tts.tsv` o `*.json`. El prefijo se empareja con el nombre
   del genoma y sus alias (`subtilis`, `bsu`, `Bacillus subtilis 168`…).
3. **Respaldos curados del repositorio**, solo si no hay tabla específica:
   - *B. subtilis*: `data/dbtbs/gt_tss.tsv` + `gt_terminators.tsv`.
   - *E. coli* K-12: `data/regulondb/tss_k12_union.tsv` + `TerminatorSet.tsv`.

Se ignoran a propósito `promoters.tsv`, `promoters_meta.json`, `benchmark.json`,
`regulon.json`, `expression_rules.tsv` y demás tablas que no son TSS/TTS. Un
JSON de resumen (solo conteos, como `subtilis_tss_tts.json`) no aporta
posiciones: si se selecciona (por ejemplo, en `data/dbtbs/`), `load` devuelve
`not_available` y añade una nota de que probablemente sea un resumen. En
`data/regulation/<genoma>/` los archivos `*_meta.json` se ignoran cuando hay
tablas con posiciones, y si son lo único presente se sigue buscando en las
otras rutas.

### Formatos aceptados

**TSV**: cabecera con nombres. Se reconocen, entre otras, estas columnas:

| campo canónico | columnas leídas |
|---|---|
| tipo | `type`, `kind`, `feature_type`, `tipo` (si falta, se infiere del nombre del archivo) |
| posición | `pos`, `position`, `tss`, `start`, `leftEndPos`, `2)leftEndPos`, `middle`, `site`… |
| fin de intervalo | `end`, `stop`, `rightEndPos`, `3)rightEndPos`, `fin`… |
| contig | `contig`, `chrom`, `chromosome`, `seqid`, `replicon`, `locus`… |
| hebra | `strand`, `hebra`, `orientation`, `4)strand`… |
| clase | `class`, `clase`, `sigma`, `strength`, `fuerza`… |
| ubicación | `location`, `ubicacion`, `anomaly`, `position_class`… |
| fuente | `source`, `fuente`, `provenance`, `origin`… |
| mecanismo | `mechanism`, `terminator_type`, `7)type`… |

**JSON**: una lista de registros, o un objeto con listas en `tss`/`tts`, o
`rows`/`records`. Cada registro puede ser un objeto con los campos anteriores o,
si la sección ya dice el tipo, un entero/string de posición.

## 5 · Normalización

### `strand`

Se reduce a `+`, `-` o `.` (desconocida):

| entrada | salida |
|---|---|
| `+`, `1`, `plus`, `forward`, `fwd`, `directo` | `+` |
| `-`, `-1`, `minus`, `reverse`, `rev`, `antisentido` | `-` |
| vacío, `?`, `unknown`, otro | `.` |

Una hebra desconocida no se enlaza a un gen: sería una asignación a ciegas.

### `pos`

Entero 1-based. Se aceptan `"1,234"`, `"250.0"` o `250`; se rechazan `0`,
negativos, decimales no enteros y texto. Las filas inválidas se cuentan en
`n_rejected`, no se corrigen.

### `source`

Se mapean alias conocidos (`bsgatlas`→`BSGatlas`, `dbtbs`→`DBTBS`,
`transferido_K12`→`transferido K-12`, `alto rendimiento`→`RegulonDB alto
rendimiento`, etc.). Una etiqueta no reconocida se conserva tal cual para no
perder procedencia.

### `class`

Depende del tipo:

- **TSS**: factor sigma canónico, p. ej. `A`/`SigA` → `sigmaA`,
  `sigma70` → `sigma70`. `?`, vacío o una clase que no es sigma → `unknown`.
- **TTS**: fuerza estructural `A_fuerte`→`strong`, `B_media`→`medium`,
  `C_debil`→`weak`, `batter`→`batter`. Sin etiqueta de fuerza → `unknown`.

`class` **no** es la ubicación. La ubicación (`intergenic`, `antisense`,
`inside_gene`, `gene_edge`, `contig_edge`) se guarda en `location`. Así no se
mezcla «sigma E» con «antisentido», que son ejes distintos.

`mechanism` guarda `rho_independent`/`rho_dependent` cuando la fuente lo trae;
es importante porque una horquilla local solo puede explicar terminadores
Rho-independientes.

### Duplicados

Clave exacta: `(type, contig, pos, strand, source)`. La primera aparición se
conserva y las demás se cuentan en `n_duplicates`. Dos fuentes distintas en la
misma coordenada **no** son duplicado: se mantienen y se cuentan en
`n_shared_positions`.

## 6 · API

```python
from webapp import tss_layer

loaded = tss_layer.load("subtilis", "data")
if loaded["status"] == "available":
    features = tss_layer.attach(features, loaded["tss"], loaded["tts"])
    print(tss_layer.summary(loaded))
```

### `load(genome, data_root)`

`genome` puede ser `"subtilis"`, `"bsu"`, `"Bacillus subtilis 168"`, el nombre
de un directorio de `data/regulation/` o el nombre de un FASTA (se le quitan
extensiones). `data_root` puede ser `data/` o la raíz del repositorio.

Devuelve:

```python
{
  "genome": "subtilis",
  "data_root": ".../data",
  "status": "available",        # available | partial | not_available
  "available": True,
  "reason": "",
  "tss": [...],
  "tts": [...],
  "n_tss": 3390,
  "n_tts": 3009,
  "n_duplicates": 0,
  "n_rejected": 0,
  "n_shared_positions": 0,
  "sources": [{"path": "...", "format": "tsv", "hint": "both",
               "source": "BSGatlas", "n_tss": 3390, "n_tts": 3009,
               "n_rejected": 0, "notes": []}],
  "searched": [".../data/dbtbs/*_tss_tts.tsv|json"],
  "limitation": "...",
}
```

Cada registro TSS/TTS tiene:

```python
{"type": "TSS", "genome": "subtilis", "contig": "basu168",
 "pos": 157, "end": None, "strand": "+", "source": "BSGatlas",
 "class": "sigmaA", "location": "intergenic", "mechanism": "",
 "detail": "evidence=..."}
```

Los `n_tss`/`n_tts` de `sources[]` son filas leídas del archivo **antes** de
quitar duplicados; los `n_tss`/`n_tts` de nivel superior son los registros
conservados. La diferencia queda visible en `n_duplicates`.

`status="partial"` significa que hay un solo tipo (solo TSS o solo TTS); no se
rellena el otro con una predicción. `status="not_available"` incluye el motivo y
las rutas buscadas.

### `attach(features, tss, tts, *, max_distance=600)`

Enlaza, para cada feature con `contig`, `start`, `end` y hebra `+`/`-`:

- `tss_upstream`: TSS más cercano **antes** del feature en la dirección de
  transcripción (`pos < start` en `+`; `pos > end` en `-`), a
  `<= max_distance`.
- `tts_downstream`: TTS más cercano **después** (`pos > end` en `+`;
  `pos < start` en `-`), a `<= max_distance`.
- `tss_upstream_distance` / `tts_downstream_distance`: pb hasta ese registro.
- `n_tss_inside` / `n_tts_inside`: marcas dentro de `[start, end]` en la misma
  hebra.

Si no hay marca cercana deja `None`; no fabrica una. La ventana de 600 pb es un
límite conservador de enlace, no un umbral biológico. La función muta los dicts
y devuelve la misma lista.

### `summary(loaded)` o `summary(tss, tts)`

Devuelve conteos (`n_tss`, `n_tts`, `n_total`, `n_duplicates`, `n_rejected`,
`n_shared_positions`), `by_strand`, `tss_classes`, `tts_classes`,
`tss_locations`, `tts_locations`, `tss_sources`, `tts_sources`, `mechanisms`,
`contigs`, `pos_min`, `pos_max`, `status` y la `limitation`.

## 7 · Límites conocidos

1. **Genoma nuevo sin tabla ⇒ `not_available`.** Es el comportamiento correcto:
   no se heredan coordenadas de otra especie ni se aproximan por homología.
2. **No hay transferencia automática entre ensamblajes.** Un TSS de
   `GCF_...` no se puede usar en otro assembly sin un mapeo explícito; la capa no
   lo hace.
3. **Algunas tablas históricas no traen `contig`.** `tss_k12_union.tsv` y
   `gt_tss.tsv` son de cromosoma único y dejan `contig=""`; `attach` no las
   enlaza con features de un contig con nombre. Es una limitación conservadora,
   no un fallo silencioso.
4. **`promoters.tsv` no es TSS.** Son promotores putativos por reglas/PWM de
   `regulation_layer.py`; esta capa los ignora para no presentar una hipótesis
   como medición.
5. **TTS en intervalo.** Para tablas con `start`/`end`, `pos` es el inicio del
   intervalo y `end` el otro extremo. No se elige un único nucleótido 3' sin
   conocer la fuente; el consumidor puede usar `end`/`start` según la hebra.
6. **Terminadores dependientes de Rho.** Una horquilla local no los detecta; si
   la fuente no los trae, no existen en la tabla.
7. **TSS transferidos.** `transferido K-12` es ortología, no medición; el
   `source` lo deja explícito.
8. **La clase no es una probabilidad.** `sigmaA` o `strong` describen la
   anotación de la fuente, no una confianza calibrada.

## 8 · Por qué no se puede predecir sin un modelo

Predecir un TSS a partir de la secuencia no es identificar un motivo corto:

- Las cajas -10/-35 son degeneradas y su espaciador varía; el mismo genoma usa
  varios factores sigma con consensos distintos.
- La posición del TSS depende de la polimerasa, factores sigma alternativos,
  proteínas accesorias, superenrollamiento y contexto; la secuencia local no
  alcanza.
- Hay promotores internos, leaderless, antisentido y UTRs largos; una regla
  posición-fija falla en todos ellos.
- En este repositorio, el PWM de promotor da ~19% de sensibilidad con ~3% de
  falsos positivos (`scripts/regulation/README.tss_tts.md`): insuficiente para
  afirmar una coordenada.
- Un terminador dependiente de Rho no tiene una estructura local que una
  horquilla pueda ver; predecirlo exige un modelo entrenado con datos de
  terminación, no una regla de secuencia.

Por eso esta capa no predice: devolver una posición inventada sería peor que
devolver `not_available`. Predecir requeriría un modelo entrenado y calibrado
(por especie o con validación cruzada), datos de entrenamiento y una evaluación
honesta de sensibilidad/especificidad; no es el alcance de esta capa CPU.

## 9 · Pruebas

```bash
.venv/bin/python -m pytest tests/test_tss_layer.py -q
```

Cubren:

- TSV combinado sintético: normalización de `strand`, `pos`, `source`, `class`,
  duplicados exactos y filas inválidas.
- JSON sintético con registros y JSON de resumen sin posiciones.
- Archivos separados `tss.tsv`/`tts.tsv` e ignorado de `promoters.tsv`.
- Ausencia total y aislamiento entre especies (un genoma nuevo no hereda
  *B. subtilis*).
- `attach` en hebras `+`/`-`, ventana, contig distinto y hebra desconocida.
- `summary` desde la carga y desde listas.
- Tabla real de *B. subtilis* si está presente: 3.390 TSS y 3.009 TTS.

## 10 · Integración

La capa no modifica `app.py`, `browser.py` ni `pipeline.py`. Para integrarla,
el llamador puede hacer `loaded = tss_layer.load(genome, data_root)` y pasar
`loaded["tss"]`/`loaded["tts"]` al navegador, o usar `attach` antes de mostrar
un gen. Si `status` es `not_available`, la UI debe mostrar la ausencia de datos
en lugar de una marca vacía o una posición aproximada.
