# PAE de interfaz para complejos

Estado: capa nueva `webapp/pae_complex.py`, construida sobre el lector local
`webapp/pae_utils.py`. Pruebas en `tests/test_pae_complex.py`. No descarga
archivos, no ejecuta AlphaFold, no modifica la GUI y no toca `app.py`,
`browser.py` ni `gene_card.py`.

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
.venv/bin/python -m pytest tests/test_pae_complex.py -q
```

Las pruebas usan matrices y JSON sintéticos, sin red ni GPU. Cubren:

- matrices asimétricas de 3 subunidades y medias/máximos por par;
- varias cadenas por subunidad y agregación ponderada;
- ambas direcciones (`a_to_b`/`b_to_a`);
- archivo con `token_chain_ids`, AFDB sin cadenas y `chains` explícito;
- archivo ausente, JSON inválido y matriz no NxN;
- mapeos vacíos, de una subunidad, con cadenas repetidas o tipos inválidos;
- cadenas desconocidas, índices solapados y rechazo de la resolución ordinal;
- resultado conservador cuando falla un bloque de una dirección;
- `interface_pair`, `subunits_from_chain_subunits` y el texto de límites.

## Relación con `docs/PAE.md`

`docs/PAE.md` describe el formato, la escala, la interpretación general y los
resúmenes por residuo/par de `pae_utils`. Este documento añade solo la capa de
complejos: mapeo de subunidades, dos direcciones por par y agregación. Las
limitaciones de `PAE.md` siguen aplicando; las de esta página las extienden al
caso de subunidades y estequiometría.
