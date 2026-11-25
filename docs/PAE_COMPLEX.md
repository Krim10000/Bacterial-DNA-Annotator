# PAE de interfaz para complejos

Estado: dos capas nuevas, construidas sobre el lector local
`webapp/pae_utils.py`.

- `webapp/pae_complex.py`: cálculo de la interfaz por par de subunidades.
  Pruebas en `tests/test_pae_complex.py`.
- `webapp/complex_pae_ui.py`: bloque para la ficha de complejos en Streamlit,
  con búsqueda local del PAE y tabla `A→B`/`B→A`. Pruebas en
  `tests/test_complex_pae_ui.py`.

Ninguna capa descarga archivos ni ejecuta AlphaFold. Este cambio **no edita**
`app.py`, `browser.py` ni `gene_card.py`: la integración se documenta abajo con
un snippet para la ficha existente.

Este documento explica **cómo se calcula la interfaz de un complejo**, cómo se
leen los estados y por qué el resultado **no demuestra interacción ni
estequiometría**. Para la definición general del PAE, formatos de archivo y
resúmenes por residuo, ver `docs/PAE.md`.

## Qué calcula

Dado un PAE y un mapeo de cadenas a subunidades:

1. Agrupa las cadenas que pertenecen a cada subunidad (una subunidad puede
   tener varias cadenas, por ejemplo un homómero).
2. Para cada par de subunidades `A`/`B`, calcula dos bloques direccionales:
   - `A -> B`: filas de las cadenas de `A`, columnas de las cadenas de `B`;
   - `B -> A`: filas de las cadenas de `B`, columnas de las cadenas de `A`.
3. Publica media y máximo por dirección, media combinada ponderada y máximo del
   par.
4. Agrega media y máximo de todo el complejo y cuenta los pares `ok`,
   `not_available` e `invalid`.

El PAE es direccional: `A -> B` y `B -> A` no son intercambiables. Un valor bajo
en una dirección y alto en la otra describe una confianza asimétrica en la
posición relativa, no una prueba de contacto.

## API

```python
from webapp.pae_complex import (
    LIMITATION,
    STATUSES,
    analyze_complex,
    interface_pair,
    subunits_from_chain_subunits,
)
```

| Función | Devuelve |
|---|---|
| `analyze_complex(source, subunits, chains=None)` | Análisis de todos los pares de subunidades: `status`, `pairs`, `summary`, `limitation`, etc. |
| `interface_pair(source, subunit_a, subunit_b, subunits, chains=None)` | Solo el par indicado, con `a_to_b` y `b_to_a`. |
| `subunits_from_chain_subunits({cadena: subunidad})` | Invierte el mapeo a `{subunidad: [cadenas]}`; agrupa homómeros. Lanza `ValueError` si el mapeo es inválido. |

`source` puede ser:

- una ruta local a un `*.pae.json` (AFDB, AlphaFold 3/Server, ColabFold, matriz
  directa o lista de modelos, según `pae_utils`);
- una matriz PAE `NxN` ya cargada.

`chains` es el mapeo opcional de cadenas a índices de residuo que acepta
`pae_utils.interface_pae`:

```python
chain_indices = {
    "A": [0, 1, 2],
    "B": {"start": 3, "end": 5},
    "C": [6, 7],
}
```

Si se aporta `chains`, tiene prioridad sobre los metadatos de cadena del JSON.
Los índices y rangos son **0-based** e inclusivos. Es la vía para usar un PAE
AFDB de monómero o una matriz sin `token_chain_ids`, siempre que el usuario
declare los límites.

### Forma del resultado de `analyze_complex`

```python
{
    "status": "ok",              # ok | not_available | invalid
    "ok": True,
    "available": True,           # el archivo existe / la matriz se aportó
    "reason": "ok",
    "source": "complejo.pae.json",
    "source_kind": "file",       # file | matrix
    "subunits": {"RpoA": ["A"], "RpoB": ["B"], "RpoC": ["C"]},
    "n_subunits": 3,
    "pairs": [
        {
            "subunit_a": "RpoA",
            "subunit_b": "RpoB",
            "status": "ok",
            "ok": True,
            "mean_pae": 5.0,     # ponderada por pares de residuos
            "max_pae": 8.0,
            "n_pairs": 8,        # pares de residuos sumando ambas direcciones
            "a_to_b": {...},
            "b_to_a": {...},
            "reason": "ok",
        },
        # ... un elemento por combinación de subunidades
    ],
    "summary": {
        "n_pairs": 3,            # pares de subunidades
        "n_ok": 3,
        "n_not_available": 0,
        "n_invalid": 0,
        "n_residue_pairs": 24,
        "mean_pae": 5.0,
        "max_pae": 8.0,
    },
    "limitation": "<texto de límites>",
}
```

Cada dirección (`a_to_b`/`b_to_a`) contiene:

- `status`, `ok`, `available`;
- `chains_a`, `chains_b`;
- `expected_chain_pairs`, `n_chain_pairs`, `n_chain_pairs_ok`;
- `n_pairs`, `mean_pae`, `max_pae`;
- `chain_pairs`: detalle por cada combinación de cadenas, con
  `requested_chain_a`, `resolved_chain_a`, `n_a`, `n_b`, `n_pairs`,
  `mean_pae`, `max_pae`, `median_pae` y `reason`.

## Mapeo de subunidades

Formato principal: `{subunidad: cadena | [cadenas]}`.

```python
# Heterómero con una cadena por subunidad.
subunits = {
    "RpoA": ["A"],
    "RpoB": ["B"],
    "RpoC": ["C"],
}

# Homómero: dos cadenas para una misma subunidad.
subunits = {
    "RpoA": ["A", "C"],
    "RpoB": ["B"],
}
```

Reglas:

- Se exigen **al menos dos subunidades**; con una sola no hay par de interfaz.
- Cada cadena debe aparecer **una sola vez** en todo el mapeo. Una cadena en
  dos subunidades es `invalid` porque la media no tendría una asignación única.
- Una lista vacía, un ID vacío, un booleano o un tipo no admitido es `invalid`.
- Las etiquetas de subunidad se recortan; dos etiquetas que colisionan tras el
  recorte son `invalid`.
- El orden de las subunidades define el orden de las combinaciones en `pairs`.
- No se infieren subunidades desde los IDs de cadena. Un PAE que declara
  cadenas `A`/`B` no dice si son dos protómeros, dos subunidades distintas o
  una cadena con un ligando. Sin mapeo, el estado es `not_available`.
- Si el mapeo disponible viene por cadena, invertirlo:

```python
grouped = subunits_from_chain_subunits({
    "A": "RpoA",
    "B": "RpoB",
    "C": "RpoA",
})
# {"RpoA": ["A", "C"], "RpoB": ["B"]}
```

### Cadenas y ordinales

`pae_utils` acepta un ID de cadena o un índice de orden 0-based. En un mapeo de
subunidades, un ordinal podría asignar una cadena equivocada en silencio. Por
eso `pae_complex` exige que el ID solicitado coincida con la cadena resuelta:
si `"0"` termina resolviendo a `"A"`, el par se marca `invalid` con un motivo
explícito. **No se inventa una cadena para completar el mapeo.**

## Estados

| Estado | Cuándo | Valores numéricos |
|---|---|---|
| `ok` | PAE válido, cadenas declaradas y todos los pares de cadenas calculables. | `mean_pae`/`max_pae` reales. |
| `not_available` | No existe el archivo PAE; el PAE no declara cadenas y no se aportó `chains`; falta el mapeo de subunidades. | Todo `None`/0; nunca se rellena con pLDDT, ipTM ni estimaciones. |
| `invalid` | El PAE existe pero es ilegible, no es NxN o tiene valores inválidos; el mapeo de subunidades es contradictorio; una cadena no existe, está repetida o los índices se solapan. | Todo `None`/0. |

Los estados aparecen en el resultado global, en cada par, en cada dirección y
en cada par de cadenas. La precedencia al agregar es:

```text
invalid > not_available > ok
```

Es decir, un par con una dirección `invalid` es `invalid` aunque la otra sea
`ok`; el detalle de la dirección válida se conserva en `a_to_b`/`b_to_a`, pero
no se publica una media combinada parcial.

`available` es distinto de `ok`: indica que el archivo existe (o que se aportó
una matriz), aunque el contenido no sea utilizable. Un archivo ausente produce
`available=False`; un JSON corrupto produce `available=True` y
`status="invalid"`.

## Agregación

Para una dirección con varias cadenas por subunidad se calculan todos los
bloques cruzados. Por ejemplo, si `A` tiene las cadenas `A1`, `A2` y `B` tiene
`B1`:

```text
A -> B = bloques A1->B1 y A2->B1
```

- La media de la dirección es la media ponderada por número de pares de
  residuos de cada bloque.
- El máximo es el máximo de los bloques.
- Si **algún** bloque de la dirección no es `ok`, la dirección no publica media
  ni máximo agregados: solo el detalle en `chain_pairs`. Es una regla
  conservadora para no presentar una media incompleta como si fuera la
  interfaz completa.

La media de un par `A`/`B` combina ambas direcciones con la misma ponderación:

```text
mean_par = (mean(A->B) * n(A->B) + mean(B->A) * n(B->A))
           / (n(A->B) + n(B->A))
```

El `summary` del complejo aplica la misma ponderación a todos los pares `ok`.
`n_residue_pairs` es la suma de pares de residuos de esos pares. No se agrega
la mediana porque no es combinable de forma exacta a partir de resúmenes; el
detalle por bloque sí incluye `median_pae`.

## Interpretación

1. **PAE bajo en ambas direcciones.** La posición relativa entre las
   subunidades está bien definida en el modelo. Es compatible con una interfaz
   resuelta, pero no la demuestra.
2. **PAE asimétrico.** Una dirección puede ser más confiable que la otra. Hay
   que mirar ambas; promediarlas sin mostrar el detalle oculta la asimetría.
3. **PAE alto.** No demuestra que no haya interacción: indica que el modelo no
   resuelve con confianza la posición relativa. Puede deberse a desorden,
   flexibilidad, una orientación alternativa o un error de plegamiento.
4. **Combinar señales.** PAE de interfaz bajo + pLDDT razonable + ipTM alto es
   más informativo que cualquiera de los tres por separado. `pae_complex` no
   calcula pLDDT ni ipTM; solo lee PAE y, si el JSON los trae, no los usa como
   sustituto.
5. **Revisar el mapeo.** Una media baja con un mapeo equivocado sigue siendo un
   número bajo. Registrar qué cadena corresponde a qué subunidad y con qué
   numeración de residuos es parte del resultado.
6. **Estequiometría.** El número de cadenas del mapeo es una declaración del
   usuario. El PAE no informa cuántas copias de cada subunidad existen en el
   ensamblado real ni en qué proporción.

Ejemplo de lectura:

```text
Subunidad A (cadenas A, C)  PAE interno bajo
Subunidad B (cadena B)      PAE interno bajo
A -> B                      media 9 Å
B -> A                      media 1 Å
```

Interpretación: el modelo confía en `B` mirando desde `A` mucho más que en `A`
mirando desde `B`. El par es asimétrico. No se puede concluir que `A` y `B`
interactúen ni que la interfaz tenga una estequiometría concreta.

## Límites

- **El PAE no demuestra interacción.** Es una expectativa de error relativo
  del modelo (Å), no una medida de contacto, distancia, afinidad ni función.
- **El PAE no demuestra estequiometría.** Tampoco la composición del complejo,
  el orden de ensamblaje ni la existencia de un complejo biológico.
- **El mapeo es una declaración.** `pae_complex` no valida que la cadena `A`
  sea realmente `RpoA`; solo comprueba que la cadena exista y no esté repetida.
  Un mapeo biológicamente incorrecto produce un resultado incorrecto.
- **No identifica residuos de interfaz.** Resume bloques; no dice qué pares de
  residuos se tocan ni si hay puentes salinos, hidrofóbicos o de otro tipo.
- **No es un mapa de contactos.** Un PAE bajo no equivale a contacto y un PAE
  alto no equivale a ausencia de contacto.
- **Una sola conformación.** No representa ensambles, alostería, *induced fit*
  ni desorden intrínseco.
- **Dependencia del modelo y de la alineación.** El valor `(i, j)` es
  condicional a alinear sobre `i`; la matriz no es simétrica.
- **Escala y redondeo.** AFDB redondea a enteros y satura en
  `max_predicted_aligned_error`. Ese máximo posible no es el máximo observado.
- **Fragmentos y numeración.** En modelos de fragmentos, los índices del PAE
  pueden no coincidir con la numeración de la proteína completa. Verificar
  `token_res_ids`/metadatos antes de cruzar con anotaciones.
- **Entidades no proteicas.** En AF3 la matriz es por token; ligandos, iones,
  agua o ácidos nucleicos pueden ocupar filas/columnas. El mapeo de subunidades
  debe declararlos si forman parte del cálculo.
- **Agregación conservadora.** Si una subunidad tiene varias cadenas y una
  combinación falla, no se publica la media de la dirección. Revisar
  `chain_pairs` para ver qué bloque falló.
- **Costo.** La matriz es O(N²) y se calculan todos los pares de cadenas de
  cada par de subunidades. En complejos grandes el costo puede ser alto; la
  lectura es local y no usa streaming.
- **Sin descarga ni cálculo de PAE.** Si el archivo no existe, no se genera un
  PAE sintético ni se reutiliza el de otro modelo. Tampoco se infiere un PAE
  desde ESMFold, PDB experimental, Foldseek o homología.

## Ficha de complejos (`complex_pae_ui`)

`webapp/complex_pae_ui.py` hace visible el PAE en la ficha de un complejo. No
descarga archivos, no ejecuta AlphaFold y **no inventa un mapeo
cadena-subunidad**. Si no hay PAE local o no hay mapeo declarado, muestra
`not_available`; no rellena la tabla con pLDDT, ipTM ni estimaciones.

### API

```python
from webapp.complex_pae_ui import (
    CAVEAT,
    complex_pae_view,
    find_complex_pae,
    render_complex_pae,
    table_rows,
)
```

| Función | Devuelve |
|---|---|
| `render_complex_pae(complex_row, data_root, key_prefix="")` | Dibuja el bloque en Streamlit y devuelve la vista. |
| `complex_pae_view(complex_row, data_root)` | La misma vista sin tocar Streamlit (para pruebas u otras integraciones). |
| `find_complex_pae(complex_row, data_root)` | Resuelve la ruta local del PAE: `status`, `path`, `reason`, `declared`, `candidates`. |
| `table_rows(view)` | Filas listas para `st.dataframe`, con media/máximo `A→B` y `B→A`. |

`complex_row` puede ser un `dict` o una fila de pandas (`Series`). Se leen
`complex_id`/`id`, `subunits` y, si existe, `pae` (o
`pae_path`/`pae_json`/`pae_file`). `data_root` es el directorio `data/` que ya
usa la app. `key_prefix` namespacia la clave del dataframe cuando la ficha se
dibuja más de una vez.

La vista siempre tiene la misma forma:

```python
{
    "status": "ok",              # ok | not_available | invalid
    "ok": True,
    "reason": "ok",
    "complex_id": "ATP_SYNTHASE:...",
    "subunits": ["atpa", "atpb"],
    "pae_path": "/.../complejo.pae.json",  # None si no se encontró
    "pae_found": True,
    "pae_declared": False,
    "declared_chains": ["A", "B"],  # si hay PAE pero falta el mapeo
    "mapping": {"RpoA": ["A"], "RpoB": ["B"]},
    "mapping_source": "complex_row",  # complex_row | sidecar | pae_json
    "analysis": {...},           # salida de webapp.pae_complex.analyze_complex
    "rows": [...],               # una fila por par de subunidades
    "summary": {...},            # resumen del complejo
    "limitation": "<texto de límites>",
    "caveat": "<aviso breve que se muestra siempre>",
}
```

### Búsqueda local del PAE

Orden:

1. **Ruta declarada** en la fila, por cualquiera de las claves `pae`,
   `pae_path`, `pae_json` o `pae_file`. Acepta un string, una lista de
   strings o un dict con `{"path": "...", "chain_subunits": {...}}`.
   - Una ruta absoluta se usa tal cual.
   - Una ruta relativa se prueba bajo `data_root`, bajo su padre (la raíz del
     proyecto) y bajo el directorio actual.
   - Una ruta remota (`http://`, `https://`, `ftp://`, `s3://`, `gs://`, etc.)
     se rechaza con `not_available`: **no se descarga**.
2. Si no hay ruta declarada (o no existe), se busca recursivamente en
   `data_root/metabolism/targets` por archivos `*.pae.json`.
   - Prioridad 1: el nombre del archivo coincide con `complex_id` o `id`
     (también se prueba la variante con `:` reemplazado por `_`).
   - Prioridad 2: el nombre contiene **todas** las subunidades declaradas, por
     ejemplo `P37494__P59653.pae.json` para `subunits="P37494;P59653"`.
   - Varios archivos con la misma prioridad: `invalid`, con los nombres en el
     motivo. **No se elige uno al azar.**
   - Un archivo que solo coincide con una subunidad no se usa para completar
     un complejo de dos o más: no se rellena la evidencia que falta.

La búsqueda automática no sale de `data_root` ni usa la red. Una ruta
declarada, en cambio, se usa tal cual (o resuelta contra las raíces indicadas
arriba), incluso si queda fuera de `data_root`. Si la ruta declarada no existe,
se informa en `declared_missing` y se intenta la búsqueda local; el estado final
sigue siendo explícito.

### Mapeo cadena-subunidad

El mapeo es **una declaración**, no una inferencia. Precedencia:

1. `complex_row`:
   - si `pae` es un dict, sus claves de mapeo;
   - `chain_subunits`, `chain_to_subunit` o `pae_chain_subunits`
     (`{cadena: subunidad}`);
   - `subunit_chains`, `subunit_to_chains` o `pae_subunit_chains`
     (`{subunidad: cadena | [cadenas]}`).
   - si `subunits` es un objeto `{subunidad: cadena | [cadenas]}`, se usa
     como mapeo explícito; sus claves son las etiquetas de subunidad.
2. Un sidecar junto al PAE. Para `complejo.pae.json` se admiten
   `complejo.chains.json`, `complejo.chain_subunits.json`,
   `complejo.subunit_chains.json`, `complejo.map.json`,
   `complejo.pae.json.chains.json` y `complejo.pae.json.map.json`. El
   contenido usa las mismas claves:
   ```json
   {"chain_subunits": {"A": "RpoA", "B": "RpoB"}}
   ```
3. Claves de mapeo en el nivel superior del propio JSON del PAE (por ejemplo,
   un PAE de AlphaFold 3 con metadatos añadidos por el pipeline).

Si una fuente presente declara un mapeo mal formado, el estado es `invalid` y
no se prueba la siguiente fuente. **No se deriva el mapeo del orden de las
cadenas, del nombre de la subunidad ni de que un ID de cadena coincida con una
etiqueta.** Un PAE que declara `A`/`B` y una fila que declara subunidades `A`/`B`
sigue siendo `not_available` sin mapeo explícito.

### Estados y mensajes

| Estado | Cuándo | Qué dibuja |
|---|---|---|
| `ok` | PAE local, cadenas y mapeo explícito; todos los pares calculables. | Resumen, tabla y aviso. |
| `not_available` | No hay PAE local; la ruta remota no se descarga; falta el mapeo; el PAE no declara cadenas. | `st.info` con el literal `not_available` y el motivo. |
| `invalid` | Varios PAE candidatos; mapeo mal formado; PAE corrupto o no NxN; cadena desconocida. | `st.warning` con el literal `invalid`; tabla solo si hay algún valor calculable. |

En `not_available`, `available` puede ser `True` (el archivo existe) pero el
contenido no permite separar subunidades. En `invalid`, no se publica una media
parcial como si fuera la interfaz completa.

### Tabla

Una fila por par de subunidades, con columnas:

```text
Par · Subunidad A · Cadenas A · Subunidad B · Cadenas B
A→B media (Å) · A→B máx (Å) · B→A media (Å) · B→A máx (Å)
Pares (residuos) · Estado · Motivo
```

Las medias son las de `pae_complex` (ponderadas por pares de residuos) y el
máximo es el máximo del bloque direccional. El resumen del complejo muestra
pares evaluados, media, máximo y número de pares de residuos. La tabla se
construye con `table_rows(view)` y no depende de pandas.

El bloque muestra siempre `CAVEAT`:

> El PAE es error relativo del modelo (Å), no una distancia ni evidencia
> experimental. Un PAE de interfaz bajo no demuestra interacción física ni
> estequiometría; la asignación de cadenas a subunidades es una declaración y
> no se infiere del PAE.

### Cómo integrarlo en la ficha existente

En `webapp/browser.py`, la ficha de complejo vive en la pestaña `t_cx`. Basta
añadir el import (el archivo usa imports planos) y una llamada después de los
metadatos actuales:

```python
from complex_pae_ui import render_complex_pae

# ... dentro de `with t_cx:`, después de mostrar la ficha del complejo:
render_complex_pae(
    r.to_dict(),
    data_root,
    key_prefix=f"cx_{genome}_{r.get('complex_id', '')}",
)
```

Si se importa como paquete, usar
`from webapp.complex_pae_ui import render_complex_pae`. La función no modifica
el estado de Streamlit salvo dibujar el bloque; devuelve la vista por si la
ficha quiere mostrar o descargar algo adicional.

Alternativa para pruebas o para otra vista:

```python
from webapp.complex_pae_ui import complex_pae_view, table_rows

view = complex_pae_view(row, data_root)
if view["status"] == "ok":
    for row in table_rows(view):
        print(row["Par"], row["A→B media (Å)"], row["B→A media (Å)"])
else:
    print(view["status"], view["reason"])
```

## Ejemplo

```python
from webapp.pae_complex import analyze_complex

result = analyze_complex(
    "data/complejo.pae.json",
    {
        "RpoA": ["A"],
        "RpoB": ["B"],
        "RpoC": ["C"],
    },
)

if result["status"] == "ok":
    for pair in result["pairs"]:
        print(
            pair["subunit_a"],
            pair["subunit_b"],
            "A->B",
            pair["a_to_b"]["mean_pae"],
            "B->A",
            pair["b_to_a"]["mean_pae"],
        )
    print("media del complejo:", result["summary"]["mean_pae"])
else:
    print(result["status"], result["reason"])
```

Para analizar un PAE sin metadatos de cadena hay que declarar los índices:

```python
result = analyze_complex(
    "AF-P1.pae.json",
    {"A": ["A"], "B": ["B"]},
    chains={"A": [0, 1], "B": [2, 3]},
)
```

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_complex_pae_ui.py tests/test_pae_complex.py -q
```

Las pruebas usan matrices, JSON y un Streamlit falso, sin red ni GPU. Cubren:

- matrices asimétricas de 3 subunidades y medias/máximos por par;
- varias cadenas por subunidad y agregación ponderada;
- ambas direcciones (`a_to_b`/`b_to_a`);
- archivo con `token_chain_ids`, AFDB sin cadenas y `chains` explícito;
- archivo ausente, JSON inválido y matriz no NxN;
- mapeos vacíos, de una subunidad, con cadenas repetidas o tipos inválidos;
- cadenas desconocidas, índices solapados y rechazo de la resolución ordinal;
- resultado conservador cuando falla un bloque de una dirección;
- `interface_pair`, `subunits_from_chain_subunits` y el texto de límites;
- búsqueda del PAE por ID, por `:` reemplazado por `_` y por todas las
  subunidades; un PAE de una sola subunidad no se acepta para un complejo;
- ruta declarada absoluta/relativa, ruta ausente y ruta remota rechazada;
- ambigüedad de candidatos como `invalid`, sin elegir al azar;
- mapeo desde la fila (incluido `subunits` como objeto), desde un sidecar y
  desde el JSON del PAE;
- ausencia de mapeo como `not_available`, incluso si los IDs de cadena y de
  subunidad coinciden;
- tabla con media/máximo `A→B` y `B→A`, homómeros con dos cadenas por
  subunidad y aviso de interacción/estequiometría;
- render con Streamlit falso para `ok`, `not_available` e `invalid`, y
  `key_prefix` en la clave del dataframe.

`tests/test_complex_pae_ui.py` no necesita un runtime de Streamlit: sustituye
`complex_pae_ui.st` por un objeto que registra llamadas y permite comprobar que
no se dibuja una tabla cuando no hay valores.

## Relación con `docs/PAE.md`

`docs/PAE.md` describe el formato, la escala, la interpretación general y los
resúmenes por residuo/par de `pae_utils`. Este documento añade solo la capa de
complejos: mapeo de subunidades, dos direcciones por par y agregación. Las
limitaciones de `PAE.md` siguen aplicando; las de esta página las extienden al
caso de subunidades y estequiometría.
