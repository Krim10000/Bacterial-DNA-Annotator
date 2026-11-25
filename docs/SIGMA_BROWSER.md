# Sigma en el navegador y vista Sigma

> **Estado:** integración de solo lectura sobre las salidas del módulo `sigma`.
> No recalcula anotaciones, no escribe datos y no importa Streamlit desde
> `modules/`. El navegador (`webapp/browser.py`) suma una pista opcional; la
> vista dedicada (`webapp/sigma_view.py`) reúne métricas, filtros y descargas.

---

## 1. Qué se agregó

### 1.1 Pista `sigma` en el navegador

`browser.render()` ahora lee `data/sigma/<genoma>/sigma_promoters.tsv` con
`load_sigma_promoters()` y, si hay filas, enciende una casilla nueva:

```text
Genes · Operones · Promotores · Horquillas · Sospecha metabólica ·
TSS / TTS · Complejos · σ promotores
```

La pista dibuja, en una banda propia por hebra:

- una **marca vertical** en el `start` del operón (posición aguas arriba de la
  unidad transcripcional);
- una **banda tenue** sobre `start`–`end` del operón, para que se vea qué
  unidad anota esa fila;
- **color por tipo de σ** (`sigma70`, `sigma54`, `sigma38`, `sigma32`,
  `sigma28`, `sigma24`, `sigmaE`);
- **opacidad y trazo por confianza** (`confidence`): alta = sólida y opaca,
  media = atenuada, baja = línea discontinua, desconocida = gris tenue.

La banda `sigma` está separada **28 px** de `complejos` y de `sospecha`; su
banda mide 18 px, así que quedan ~11 px de aire a cada lado. En la hebra `−`
esto evita el solapamiento que había con la banda de complejos.

La misma pista aparece en:

- la **leyenda** del navegador (colores y significado de la opacidad);
- la tabla **«Todos los elementos en pantalla»** (`kind="sigma"`, una fila por
  operón, con `id` = `operon_id`);
- la **ficha al hacer clic** (en el SVG o en la tabla), con σ, tipo, confianza,
  `evidence`, `detected_by`, `alternatives`, `provenance` y los soportes de
  ruta/expresión/inventario/localización;
- los enlaces compartibles `?sel=sigma:<operon_id>`.

### 1.2 Vista dedicada `sigma_view.render(data_root, genome)`

`webapp/sigma_view.py` pinta:

- **métricas**: total de operones anotados, `n asignados`
  (`sigma_status=assigned`), ambigüos, sin resolver y confianza
  fuerte/media/baja;
- **resumen por tipo de σ** (total, asignados y confianza por tipo);
- **filtros** por tipo de σ y por `confidence_label`, y **búsqueda** por operón
  o gen;
- **tabla** con σ, clase, arquitectura, confianza, estado, score, p-valor,
  evidencia y `detected_by`;
- **detalle de fila** con `evidence`, `detected_by`, `alternatives` y
  `provenance` (JSON legible);
- **descargas TSV y JSON** de exactamente la selección filtrada.

---

## 2. Cómo se resuelven las coordenadas (decisión documentada)

`sigma_promoters.tsv` **no trae coordenadas propias**: una fila por operón.
`load_sigma_promoters()` resuelve, en este orden:

1. el operón completo en `data/regulation/<genoma>/operons.tsv` (mismo
   `operon_id`): se usan `start`, `end`, `contig` y `strand`;
2. el `start` codificado en el propio `operon_id`
   (`genoma:contig:start:hebra`), con `end = start`;
3. el primer gen (`first_gene`) si se puede identificar en `genes.tsv`.

**La marca se dibuja en el `start` del operón**, no en un TSS: el módulo sigma
escanea una ventana aguas arriba y el TSV no conserva la posición exacta del
promotor. Es una aproximación explícita, no una medición. Si no se puede
resolver ninguna coordenada, la fila se conserva en la tabla sigma pero no se
dibuja ni entra en «elementos en pantalla».

---

## 3. Carga defensiva

`load_sigma_promoters(data_root, genome)` nunca lanza por datos faltantes:

| Situación | Resultado |
|---|---|
| No existe `data/sigma/<genoma>/sigma_promoters.tsv` | DataFrame vacío con columnas estables + aviso en la UI |
| Archivo vacío o solo cabecera | DataFrame vacío + aviso |
| TSV con columnas faltantes | rellena las columnas del esquema con `""` |
| `operons.tsv`/`genes.tsv` ausentes | usa el fallback del `operon_id` o deja la fila sin coordenadas |

Los metadatos quedan en `df.attrs`:

```python
df.attrs["sigma_available"]  # bool
df.attrs["sigma_source"]     # ruta del TSV
df.attrs["sigma_warning"]    # texto para la UI ("" si todo bien)
```

La pista **no rompe** el resto del navegador: genes, operones, promotores,
horquillas, TSS/TTS, complejos, sospechas, solapamientos y la vista circular
siguen funcionando igual.

`load_tables()` deja además `tables["sigma"]` con esa misma tabla (vacía si no
hay archivo), para que el SVG, la tabla de pantalla y la ficha compartan una
sola fuente.

---

## 4. Integración final en `app.py`

La vista quedó como **modo propio**, no como una sección debajo del navegador:

```text
Menú lateral → 🧭 Explorar genoma → 🧪 Sigmas
```

En `app.py`:

- `import sigma_view` junto a los demás imports de `webapp/`;
- `MODE_LABELS["sigma"] = "mode_sigma"`;
- `MENU_GROUPS` incluye `"sigma"` en el bloque `explorar`;
- `?modo=sigmas` y `?modo=sigma` resuelven a `sigma`;
- el despacho llama `sigma_view.render(DATA_ROOT, elegido)` con la especie
  activa.

El navegador conserva su pista **σ promotores**; la vista Sigma es la página
de análisis, filtros y descargas. No hace falta tocar `run_pipeline.sh`,
`modules/sigma/`, `scripts/sigma/` ni `data/sigma/`.

---

## 5. Tests

```bash
.venv/bin/python -m pytest \
  tests/test_sigma_browser.py \
  tests/test_browser_strands.py \
  tests/test_browser_ux.py -q
```

`tests/test_sigma_browser.py` cubre, sin levantar Streamlit:

- carga de un TSV sintético y resolución por operón/`operon_id`;
- archivos ausentes, vacíos y con columnas mínimas;
- forma estable de la tabla y helpers de confianza/color/JSON;
- SVG: filtrado por contig/ventana, color por tipo, opacidad/trazo por
  confianza y `data-kind="sigma"`;
- integración en «elementos en pantalla», `element_span` y la ficha;
- métricas, filtros y tabla de `sigma_view` como funciones puras;
- comprobación por código fuente de que `render` y la ficha están cableados.

Los tests existentes del navegador siguen pasando sin cambios.

---

## 6. Límites

- **No confirma nada.** `sospecha_fuerte` es una hipótesis priorizada, no un
  promotor medido. La pista muestra lo que el módulo anotó, con su confianza.
- **Gold K-12:** top-1 57,0 % frente a una mayoría trivial de 62,6 %;
  top-3 90,2 %. La vista muestra este aviso; no es una validación.
- **La posición no es un TSS.** Es el `start` del operón; el TSV no conserva la
  coordenada exacta del promotor.
- **Una fila por operón.** No hay marcas por sitio individual ni por alternativa
  de σ; las alternativas se ven en la ficha (`alternatives`).
- **Datos precalculados.** El navegador y la vista no recalculan PWMs, fondo
  nulo, inventario ni expresión. Si falta el TSV, se avisa y no se inventa.
- **`expression_support=not_available`** significa que no había tabla de
  expresión, no que el gen no se exprese.
- **Genomas sin salida sigma** no muestran la pista ni la vista; el resto del
  navegador sigue disponible.
