# PAE (Predicted Aligned Error)

Estado: lector local en `webapp/pae_utils.py`, con pruebas en
`tests/test_pae_utils.py`. No descarga archivos, no ejecuta AlphaFold y no
modifica la GUI.

## Qué es

El **PAE** es una matriz que estima el error de la **posición relativa** entre
pares de residuos de un modelo de AlphaFold. Para la celda `(i, j)`:

> si se alinea la estructura sobre el residuo `i`, ¿cuánto se espera que se
> desvíe la posición del residuo `j`?

- Se mide en **Ångström (Å)**.
- **Valores bajos = mejor**: hay más confianza en la posición relativa.
- La diagonal es 0 (un residuo respecto de sí mismo).
- **No es simétrica**: `PAE[i][j]` y `PAE[j][i]` pueden diferir porque la
  alineación se hace sobre un residuo u otro.

La matriz es `NxN`, con una fila/columna por residuo (o por token en modelos
que incluyen ácidos nucleicos u otras entidades). En AFDB el JSON trae:

```json
{
  "predicted_aligned_error": [[0, 12, 3], [12, 0, 4], [3, 4, 0]],
  "max_predicted_aligned_error": 31
}
```

`max_predicted_aligned_error` es el **máximo posible** de la escala (no el
máximo observado en esa matriz). AFDB redondea los valores a enteros y usa un
corte visual habitual de 31 Å. AlphaFold 3, AlphaFold Server y ColabFold suelen
usar la clave `pae` y añadir `token_chain_ids` para saber a qué cadena
pertenece cada fila/columna.

## Qué NO es

- **No es una distancia** entre residuos.
- **No es un mapa de contactos** ni una probabilidad de contacto.
- **No mide la calidad local** del plegamiento; para eso está pLDDT.
- **No es ipTM** ni una probabilidad de que exista un complejo.
- **No es evidencia experimental**.
- **No distingue** por sí solo entre desorden, flexibilidad, una orientación
  alternativa de dominios o un error de plegamiento.

## Cómo lo lee el proyecto

`load_matrix(path)` y `load_pae(path)` aceptan:

| Formato | Ejemplo de claves | Origen típico |
|---|---|---|
| AFDB | `predicted_aligned_error`, `max_predicted_aligned_error` | AlphaFold DB |
| AlphaFold 3 / Server | `pae`, `token_chain_ids` | AlphaFold Server / AF3 |
| ColabFold | `pae`, `plddt`, `ptm`, `iptm` | salida de ColabFold |
| Variantes | `pae_matrix`, `aligned_error`, `matrix` | conversiones locales |
| Matriz directa | `[[...], ...]` | JSON mínimo |
| Lista de modelos | `[{"pae": [...]}, ...]` | ColabFold; se usa el primer modelo con PAE |

### API

```python
from webapp.pae_utils import (
    interface_pae,
    load_matrix,
    load_pae,
    mean_pae,
    max_pae,
    pair_summary,
    residue_summary,
    summarize,
)
```

| Función | Devuelve |
|---|---|
| `load_matrix(path)` | `list[list[float]]` NxN; lanza `FileNotFoundError` si no existe y `ValueError` si el JSON o la matriz son inválidos. |
| `load_pae(path)` | Diccionario con `matrix`, `n`, `mean`, `max`, `chains`, `chain_order`, `chain_reason`, `plddt`, `ptm`, `iptm`, `max_predicted_aligned_error`, `source_key`, `model_index`, `limitation`. |
| `summarize(path)` | Solo `available`, `ok`, `n`, `mean`, `max`, `reason`, `limitation`. No lanza por archivo ausente o inválido. |
| `mean_pae(matrix)` | Media de toda la matriz, incluida la diagonal. |
| `max_pae(matrix)` | Máximo de toda la matriz, incluida la diagonal. |
| `residue_summary(matrix)` | Lista por residuo: `residue` (1-based), `mean_pae`, `max_pae`, `min_pae`, excluyendo el par consigo mismo. |
| `pair_summary(matrix)` | Resumen de pares ordenados fuera de la diagonal, más `highest` y `lowest` con `i`/`j` 1-based. |
| `interface_pae(source, chain_a, chain_b, chains=None)` | PAE de la interfaz direccional `chain_a -> chain_b`, con `n_a`, `n_b`, `n_pairs`, `mean`, `max`, `median` (y alias `mean_pae`/`max_pae`), `ok` y `reason`. |

Convenciones:

- Las matrices se manejan con índices **0-based**; los resúmenes por residuo y
  por par muestran numeración **1-based**.
- `summarize.mean` y `summarize.max` incluyen la diagonal. Para el resumen por
  residuo y por par se excluye el par `i == i`.
- `pair_summary` considera pares **ordenados**: `(i, j)` no es lo mismo que
  `(j, i)`.
- `interface_pae` usa filas de `chain_a` y columnas de `chain_b`; es
  direccional y conviene mirar ambas direcciones.
- `chain_a`/`chain_b` aceptan el ID de cadena (`"A"`) o su índice de orden
  0-based (`0`). Si el JSON no declara cadenas, no se inventan límites: se
  devuelve `ok=False` con el motivo.
- El argumento opcional `chains` permite aportar cadenas externas, por
  ejemplo `{"A": [0, 1], "B": [2, 3]}` o `{"A": {"start": 0, "end": 1}}`.
  Esos índices y rangos son **0-based** e inclusivos.

### Semántica de `available` y `ok`

- `available=True`: el archivo existe.
- `ok=True`: además, el JSON contiene una matriz PAE válida y `n`, `mean` y
  `max` son reales.
- Archivo ausente: `available=False`, `ok=False`, `n=0`, `mean=None`,
  `max=None`.
- Archivo corrupto o no NxN: `available=True` (existe), `ok=False`,
  `n=0`, `mean=None`, `max=None`, `reason` explica el problema.

`reason` es siempre explícito y `limitation` recuerda el alcance del PAE.
Cuando no hay PAE no se rellena con pLDDT ni con valores estimados.

## pLDDT vs PAE vs ipTM

| Métrica | Qué mide | Escala | Nivel | Uso típico | Límite |
|---|---|---|---|---|---|
| **pLDDT** | Confianza local en la posición de cada residuo | 0–100 (mayor mejor) | Por residuo | Calidad local, regiones desordenadas | No dice si dos dominios están bien colocados entre sí |
| **PAE** | Error esperado de la posición relativa de un residuo al alinear sobre otro | Å (menor mejor) | Par de residuos / dominios | Empaquetamiento y orientación de dominios, interfaz | No es distancia ni contacto; no prueba interacción |
| **ipTM** | Confianza en la interfaz de un complejo | 0–1 (mayor mejor) | Interfaz / complejo | Priorizar interfaces en multímeros | Depende de la predicción; no es afinidad ni estequiometría |
| **pTM** | Confianza global del modelo (tipo TM-score) | 0–1 (mayor mejor) | Modelo completo | Calidad global comparativa | Puede ocultar un dominio mal orientado |

En AFDB el **pLDDT** va en la columna B-factor del PDB/mmCIF. El PAE va en un
JSON aparte. En complejos, **ipTM** y **pTM** suelen venir en el JSON de
confianza de AlphaFold-Multimer/AlphaFold 3, no en el PAE de AFDB de un
monómero. Este lector guarda `ptm`/`iptm` si están presentes, pero no los
calcula.

### Umbrales orientativos

No hay un corte universal; estos rangos son convenciones de lectura:

| pLDDT | Lectura |
|---|---|
| > 90 | Muy alta confianza local |
| 70–90 | Confianza alta |
| 50–70 | Baja confianza |
| < 50 | Muy baja; probablemente desordenado o mal modelado |

| PAE | Lectura |
|---|---|
| < 5 Å | Posición relativa confiable |
| 5–10 Å | Confianza intermedia |
| 10–15 Å | Poco confiable |
| > 15 Å | Orientación/posición relativa incierta |

Para ipTM se usan a menudo cortes alrededor de 0.6–0.8, pero dependen del
sistema y del objetivo. El umbral nunca debe presentarse como una probabilidad
calibrada de que la interfaz sea real.

## Cómo interpretarlo

1. **Bloques en la diagonal.** Un bloque de PAE bajo sugiere un dominio o una
   unidad rígida con buena confianza interna.
2. **Bloques fuera de la diagonal.** PAE bajo entre dos bloques sugiere que su
   orientación relativa está bien definida. PAE alto entre bloques, con PAE
   bajo dentro de cada uno, sugiere que los dominios están bien plegados pero
   su empaquetamiento relativo es incierto (o que hay flexibilidad real).
3. **Interfaz.** Comparar `A -> B` y `B -> A`. PAE bajo y simétrico en ambas
   direcciones refuerza la confianza en la interfaz; un valor alto no demuestra
   que no haya interacción, solo que el modelo no la resuelve con confianza.
4. **Combinar señales.** PAE bajo en la interfaz + ipTM alto + pLDDT razonable
   es más informativo que cualquiera de los tres por separado. PAE bajo con
   pLDDT bajo en la región merece cautela.
5. **Declarar la procedencia.** Un PAE de AFDB, AlphaFold Server, ColabFold o
   una conversión local no son intercambiables sin registrar versión, origen y
   mapeo de residuos.

Ejemplo de lectura:

```text
Dominio 1 (residuos 1-120)   PAE interno bajo
Dominio 2 (residuos 121-300) PAE interno bajo
PAE 1-120 vs 121-300         alto y asimétrico
```

Interpretación: el modelo confía en cada dominio por separado, pero no en cómo
están colocados uno respecto del otro. No implica que la proteína esté mal
plegada ni que los dominios no interactúen.

## Límites

- **Expectativa, no medida.** El PAE es una salida del modelo; no reemplaza
  cristalografía, cryo-EM, SAXS, crosslinking ni validación funcional.
- **Una sola conformación.** No representa ensambles, cambios alostéricos,
  induced fit ni desorden intrínseco.
- **Depende de la alineación.** El significado de `(i, j)` es condicional a
  alinear sobre `i`; la matriz no es simétrica.
- **No identifica residuos de interfaz.** `interface_pae` resume un bloque;
  no dice qué pares concretos forman la interfaz.
- **No da afinidad, cinética, estequiometría ni función.**
- **Escala y redondeo.** AFDB redondea a enteros y satura en un máximo
  (`max_predicted_aligned_error`); no confundir ese máximo posible con el
  máximo observado.
- **Costo cuadrático.** Una matriz `NxN` ocupa O(N²); para proteínas muy
  grandes el JSON y el cálculo de resúmenes pueden ser costosos. La lectura se
  hace localmente y sin streaming.
- **Numeración y fragmentos.** En modelos de fragmentos, los índices del PAE
  pueden no coincidir con la numeración de la proteína completa. Verificar
  `token_res_ids`/metadatos antes de cruzar con anotaciones.
- **Cadenas.** Sin `token_chain_ids`, `chain_ids`, `chains` o un mapeo externo,
  no se puede separar una interfaz de forma segura. El módulo prefiere fallar
  con un motivo antes que suponer límites.
- **Entidades no proteicas.** En AF3 la matriz es por token; ligandos, iones,
  agua o ácidos nucleicos pueden ocupar filas/columnas. No todos los flujos de
  AFDB clásico incluyen esas entidades.
- **No sustituye pLDDT/ipTM.** Cada métrica responde una pregunta distinta.

## Por qué no hay PAE para todos

El PAE no está dentro del PDB/mmCIF: es un archivo adicional. Por eso puede
faltar aunque la estructura esté disponible.

- **AFDB**: el JSON de PAE se descarga aparte. La API expone `paeDocUrl`, pero
  no todas las copias locales lo tienen; las descargas masivas por FTP pueden
  incluir solo coordenadas y confianza agregada. Conviene comprobar la URL
  antes de asumir que existe.
- **ESMFold**: produce pLDDT y coordenadas, pero no un PAE comparable al de
  AlphaFold. El pipeline local puede generar un modelo ESMFold para una
  proteína sin entrada AFDB; en ese caso no hay PAE.
- **PDB experimental**: no trae PAE predicho. Se podría calcular con un modelo
  propio, pero no es parte de este lector.
- **Modelos por homología, Foldseek, AlphaFill o plantillas**: no traen PAE de
  AlphaFold.
- **Complejos**: AlphaFold-Multimer/Server/ColabFold pueden producir PAE, pero
  la interfaz solo es calculable si el archivo conserva la asignación de
  cadena por residuo/token. Un PAE de monómero AFDB no la tiene.
- **Fragmentos y versiones**: el PAE puede existir para un fragmento y no para
  la proteína completa, o corresponder a otra versión del modelo. Hay que
  verificar el identificador y la numeración.

En el proyecto, cuando falta el PAE:

- `summarize(path)` devuelve `available=False`/`ok=False` y un `reason`;
- `interface_pae(...)` devuelve `ok=False` y explica si faltan cadenas;
- no se rellena el hueco con pLDDT, ipTM ni un valor inventado.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_pae_utils.py -q
```

Las pruebas usan JSON sintéticos y no requieren red ni GPU. Cubren formato
AFDB, AF3/Server/ColabFold, matriz directa, lista de modelos, archivo ausente,
JSON corrupto, matrices no NxN, valores negativos/NaN, resúmenes por
residuo/par, interfaz direccional, cadenas explícitas, cadenas solapadas y
metadatos de cadena inválidos.

## Referencias

- EMBL-EBI Training, *Downloading files* (campos del JSON PAE de AFDB):
  <https://www.ebi.ac.uk/training/online/courses/navigating-alphafold-database/understanding-the-structure-prediction-page/downloading-files/>
- AlphaFold Protein Structure Database, FAQ y API:
  <https://alphafold.ebi.ac.uk/faq> · <https://alphafold.ebi.ac.uk/api-docs>
- Jumper et al., *Highly accurate protein structure prediction with
  AlphaFold*, Nature (2021):
  <https://doi.org/10.1038/s41586-021-03819-2>
- Varadi et al., *AlphaFold Protein Structure Database in 2024*, Nucleic Acids
  Research (2024):
  <https://doi.org/10.1093/nar/gkad1011>
- Evans et al., *Protein complex prediction with AlphaFold-Multimer* (2022):
  <https://doi.org/10.1101/2021.10.04.463034>
- Abramson et al., *Accurate structure prediction of biomolecular interactions
  with AlphaFold 3*, Nature (2024):
  <https://doi.org/10.1038/s41586-024-07487-w>
