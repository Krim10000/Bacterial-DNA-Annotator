🌐 [Español](ARCHITECTURE.md) · **English**

# Application architecture

## Bird's-eye view

```text
Browser (Streamlit) — webapp/app.py
  ├─ i18n (es/en/zh) — webapp/i18n.py
  ├─ texts and credits — textos.py · credits.py
  ├─ clickable genome browser — browser.py + components/
  └─ requirements/times/GPU — requirements_ui.py
        │
        ├── Annotation — annotator.py · pipeline.py
        └── Already computed data — data/<module>/*.tsv (400 genomes)
                 │
      pyrodigal · tRNAscan-SE · barrnap · DIAMOND · ghost
                 │
      GFF3 · proteins · TSV · Markdown report · JSON
```


The GUI **does not recompute** the 400 genomes: it reads already generated
tables and displays them. Only the *Annotate FASTA* mode runs tools on new
data. That separation is what makes it possible to explore without waiting and
to keep the pipeline reproducible.

## Module layer (`modules/`)

Each domain is declared in `modules/<id>/module.py` with a `MODULE`:

- `options`: which controls the GUI draws (bool, int, float, text, select).
- `status(data_root)`: available/incomplete and why.
- `estimate(config)`: approximate time.
- `selftest()`: quick test without downloads.
- `run(ctx, inputs)`: optional execution with outputs and metrics.

`modules/registry.py` discovers the folders with `pkgutil`; `webapp/module_ui.py`
builds the **🧩 Modules** view without central lists. The tests
`tests/test_modules.py` validate discovery, options and a real execution.
The contract is in [`MODULES.en.md`](MODULES.en.md).

Modules do not import Streamlit. `modules/_bridge.py` adds `webapp/` to
`sys.path` to reuse `annotator`, `metabolism`, etc.

## `webapp/` modules

| file | responsibility |
|---|---|
| `app.py` | interface composition, modes and navigation. It does **not** contain scientific logic. |
| `annotator.py` | annotation of a FASTA: CDS, tRNA, rRNA, ghost, deserts, comparison. |
| `pipeline.py` | same path as the GUI, without Streamlit (`dna-annotate`). |
| `browser.py` | genome browser: tracks, search, cards on click. |
| `components/genome_track/` | build-free HTML/JS component that draws the SVG and returns the click. |
| `i18n.py` + `locales/` | translation keyed by Spanish text, selector and Streamlit patch. |
| `textos.py` | long explanatory texts, resolved in the active language. |
| `credits.py` | authorship, tools, licenses and citations. |
| `requirements_ui.py` | requirements, times, glossary and GPU status. |
| `pathways.py`, `metabolism.py` | catalog, KEGG and metabolic network. |
| `operons.py`, `closures.py` | operons/regulation and gap closures. |
| `structure_evidence.py`, `structural_filters.py` | structural filters F1–F5. |
| `esmfold_ui.py`, `esmfold_tools.py` | 3D structure via API or local GPU. |
| `hypotheses_ui.py`, `structural_hypotheses.py` | reviewable hypotheses. |
| `family_logic.py`, `gap_logic.py`, `exception_logic.py` | explanation rules (no GUI). |
| `transport_layer.py`, `producer_homology.py`, `active_sites.py` | evidence layers. |

Rule: **a logic module does not import Streamlit** except those ending in
`_ui` or `browser`/`app`. That way they can be tested without an interface.

## i18n

- **The Spanish text is the key**. `tr("texto")` looks up the exact translation
  and, if it does not exist, leaves the original. It never fails.
- `T["clave"]` in `textos.py` translates on the fly; Spanish remains the
  source.
- `install_streamlit_i18n()` wraps the `DeltaGenerator` methods to translate
  labels and texts without touching the ~500 GUI calls. It does not translate
  data values or filter options.
- To add strings: write in Spanish → `python tools/extract_ui_strings.py`
  → translate the new keys in `locales/en.json` and `locales/zh.json`.
  `tests/test_i18n.py` fails if a translation is missing or a placeholder
  breaks.
- To add a language: copy `en.json`, translate and add the code in
  `LANGUAGES` (`webapp/i18n.py`).

## Included demo

`webapp/demo/` contains 150 kb of *Chlorobaculum tepidum* annotated from
scratch (137 genes, 34 clusters by proximity). `browser.load_tables()` uses it
when the genome is `demo` and there is no local data. Rows are marked
`status=demo` and `confidence=low`; they are not mixed with the project
metrics. It is regenerated with `python tools/build_demo_data.py`.

## Clickable browser

1. `browser.svg()` generates an SVG with tracks and groups
   `<g class="hit" data-kind="gene" data-id="b0001">`.
2. The `genome_track` component receives the SVG through `args`, inserts it
   and listens for clicks. On press, it sends `setComponentValue({kind, id, ts})`.
3. `browser.render()` receives the value and draws a card with the data from
   the corresponding table.
4. **Fallback:** if Streamlit does not process the value (this happens with
   some large SVGs), the component waits ~1.2 s and reloads the app with
   `?sel=kind:id&pos=...`, which `render()` interprets the same way.

The component is pure HTML+JS: there is no `npm`, `build` or new dependencies.
Protocol documented in `components/genome_track/README.md`.

## Data and states

- `data/` is organized by module (`metabolism`, `regulation`, `operon_finder`,
  `membrane`, …). Large outputs are outside git; see `DATA_NOTICE.md`.
- Every table declares `evidence`, `confidence` and `status` when applicable.
- `confirmed` is not used without experimental/EC/context validation.
- The pipeline is **idempotent**: each stage can be re-run without duplicating
  results.

## Tests

```bash
python -m pytest tests -q                 # app, i18n, SVG, instalador, demo
cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests
```

- `tests/test_i18n.py` — coverage and format of translations.
- `tests/test_browser_svg.py` — valid XML, escaping and clickable elements.
- `tests/test_installer.py` — installer syntax and component protocol.
- `tests/test_annotate_demo.py` — real annotation of `webapp/samples/demo.fna`.
- `tests/test_imports.py` — the modules load without Streamlit at runtime.

CI (`.github/workflows/ci.yml`) runs these tests on Python 3.10, 3.11 and
3.12, plus a real Streamlit startup.
