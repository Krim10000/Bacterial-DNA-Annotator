# Arquitectura de la aplicación

## Vista de pájaro

```
┌──────────────────────────────────────────────────────────────┐
│  Navegador (Streamlit)          webapp/app.py                │
│  ├─ i18n (es/en/zh)             webapp/i18n.py               │
│  ├─ textos y créditos           textos.py · credits.py       │
│  ├─ navegador clicable          browser.py + components/     │
│  └─ requisitos/tiempos/GPU      requirements_ui.py           │
└───────────────┬──────────────────────────────────────────────┘
                │
        ┌───────▼────────┐        ┌──────────────────────────┐
        │  Anotación     │        │  Datos ya calculados     │
        │  annotator.py  │        │  data/<módulo>/*.tsv     │
        │  pipeline.py   │        │  (400 genomas)           │
        └───────┬────────┘        └──────────────────────────┘
                │
     pyrodigal · tRNAscan-SE · barrnap · DIAMOND · ghost
                │
        GFF3 · proteínas · TSV · informe Markdown · JSON
```

La GUI **no recalcula** los 400 genomas: lee tablas ya generadas y las muestra.
Solo el modo *Anotar FASTA* ejecuta herramientas sobre datos nuevos. Esa
separación es lo que permite explorar sin esperar y mantener el pipeline
reproducible.

## Capa de módulos (`modules/`)

Cada dominio se declara en `modules/<id>/module.py` con un `MODULE`:

- `options`: qué controles dibuja la GUI (bool, int, float, text, select).
- `status(data_root)`: disponible/incompleto y por qué.
- `estimate(config)`: tiempo aproximado.
- `selftest()`: prueba rápida sin descargas.
- `run(ctx, inputs)`: ejecución opcional con salidas y métricas.

`modules/registry.py` descubre las carpetas con `pkgutil`; `webapp/module_ui.py`
construye la vista **🧩 Módulos** sin listas centrales. Las pruebas
`tests/test_modules.py` validan descubrimiento, opciones y una ejecución real.
El contrato está en [`MODULES.md`](MODULES.md).

Los módulos no importan Streamlit. `modules/_bridge.py` añade `webapp/` al
`sys.path` para reutilizar `annotator`, `metabolism`, etc.

## Módulos de `webapp/`

| archivo | responsabilidad |
|---|---|
| `app.py` | composición de la interfaz, modos y navegación. **No** contiene lógica científica. |
| `annotator.py` | anotación de un FASTA: CDS, tRNA, rRNA, ghost, desiertos, comparación. |
| `pipeline.py` | misma ruta que la GUI, sin Streamlit (`dna-annotate`). |
| `browser.py` | navegador de genoma: pistas, búsqueda, fichas al hacer clic. |
| `components/genome_track/` | componente HTML/JS sin build que dibuja el SVG y devuelve el clic. |
| `i18n.py` + `locales/` | traducción por texto español como clave, selector y parche de Streamlit. |
| `textos.py` | textos explicativos largos, resueltos en el idioma activo. |
| `credits.py` | autoría, herramientas, licencias y citas. |
| `requirements_ui.py` | requisitos, tiempos, glosario y estado de GPU. |
| `pathways.py`, `metabolism.py` | catálogo, KEGG y red metabólica. |
| `operons.py`, `closures.py` | operones/regulación y cierres de brecha. |
| `structure_evidence.py`, `structural_filters.py` | filtros estructurales F1–F5. |
| `esmfold_ui.py`, `esmfold_tools.py` | estructura 3D por API o GPU local. |
| `hypotheses_ui.py`, `structural_hypotheses.py` | hipótesis revisables. |
| `family_logic.py`, `gap_logic.py`, `exception_logic.py` | reglas de explicación (sin GUI). |
| `transport_layer.py`, `producer_homology.py`, `active_sites.py` | capas de evidencia. |

Regla: **un módulo de lógica no importa Streamlit** salvo los que terminan en
`_ui` o `browser`/`app`. Así se pueden probar sin interfaz.

## i18n

- El **texto en español es la clave**. `tr("texto")` busca la traducción exacta
  y, si no existe, deja el original. Nunca falla.
- `T["clave"]` de `textos.py` traduce al vuelo; el español sigue siendo la
  fuente.
- `install_streamlit_i18n()` envuelve los métodos de `DeltaGenerator` para
  traducir etiquetas y textos sin tocar las ~500 llamadas de la GUI. No traduce
  valores de datos ni opciones de filtros.
- Para añadir cadenas: escribe en español → `python tools/extract_ui_strings.py`
  → traduce las claves nuevas en `locales/en.json` y `locales/zh.json`.
  `tests/test_i18n.py` falla si falta una traducción o se rompe un placeholder.
- Para añadir un idioma: copia `en.json`, traduce y agrega el código en
  `LANGUAGES` (`webapp/i18n.py`).

## Demo incluida

`webapp/demo/` contiene 150 kb de *Chlorobaculum tepidum* anotados desde cero
(137 genes, 34 agrupaciones por proximidad). `browser.load_tables()` la usa
cuando el genoma es `demo` y no hay datos locales. Las filas se marcan
`status=demo` y `confidence=low`; no se mezclan con las métricas del proyecto.
Se regenera con `python tools/build_demo_data.py`.

## Navegador clicable

1. `browser.svg()` genera un SVG con pistas y grupos
   `<g class="hit" data-kind="gene" data-id="b0001">`.
2. El componente `genome_track` recibe el SVG por `args`, lo inserta y escucha
   clics. Al pulsar, envía `setComponentValue({kind, id, ts})`.
3. `browser.render()` recibe el valor y dibuja una ficha con los datos de la
   tabla correspondiente.
4. **Respaldo:** si Streamlit no procesa el valor (ocurre con algunos SVG
   grandes), el componente espera ~1,2 s y recarga la app con
   `?sel=kind:id&pos=...`, que `render()` interpreta igual.

El componente es HTML+JS puro: no hay `npm`, `build` ni dependencias nuevas.
Protocolo documentado en `components/genome_track/README.md`.

## Datos y estados

- `data/` se organiza por módulo (`metabolism`, `regulation`, `operon_finder`,
  `membrane`, …). Las salidas grandes están fuera de git; ver `DATA_NOTICE.md`.
- Toda tabla declara `evidence`, `confidence` y `status` cuando aplica.
- `confirmed` no se usa sin validación experimental/EC/contexto.
- El pipeline es **idempotente**: cada etapa se puede volver a ejecutar sin
  duplicar resultados.

## Pruebas

```bash
python -m pytest tests -q                 # app, i18n, SVG, instalador, demo
cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests
```

- `tests/test_i18n.py` — cobertura y formato de traducciones.
- `tests/test_browser_svg.py` — XML válido, escape y elementos clicables.
- `tests/test_installer.py` — sintaxis del instalador y protocolo del componente.
- `tests/test_annotate_demo.py` — anotación real de `webapp/samples/demo.fna`.
- `tests/test_imports.py` — los módulos cargan sin Streamlit en ejecución.

La CI (`.github/workflows/ci.yml`) corre estas pruebas en Python 3.10, 3.11 y
3.12, además de un arranque real de Streamlit.
