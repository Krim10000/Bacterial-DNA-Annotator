🌐 [Español](MODULES.md) · **English**

# Self-discovered modules

The application has two layers:

1. **Core** (`webapp/`): annotation, browser, data, i18n.
2. **Modules** (`modules/`): each domain declares its options, requirements,
   time, status, self-test and — if applicable — a run function.

The GUI **has no hand-written tab per module**. The
**🧩 Modules** view discovers the folders and builds on its own:

- the module selector,
- the checkboxes, numbers, texts and selectors,
- the status (`✅ available` / `⚠️ incomplete`) and the reason,
- the estimated time,
- the **🧪 Test module** button,
- the **▶️ Run module** button if the module implements it.

## Adding a module

```text
modules/
  mi_modulo/
    __init__.py
    module.py     # exporta MODULE
```

Minimal `module.py`:

```python
from modules.base import Module, ModuleStatus, Option

def _selftest() -> str:
    return "carga OK"

MODULE = Module(
    id="mi_modulo",
    name="Mi módulo",
    summary="Qué hace, en una frase.",
    order=70,
    category="genómica",
    options=(
        Option("activo", "Activar análisis", "bool", True,
               help="Qué hace esta opción."),
        Option("threads", "Hilos", "int", 4, min=1, max=64,
               advanced=True),
    ),
    selftest_fn=_selftest,
)
```

On saving the folder:

```bash
python -m modules --list
python -m modules --selftest mi_modulo
```

It appears in the GUI without touching `app.py`.

## Contract

### `Option`

| field | use |
|---|---|
| `key` | unique identifier within the module |
| `label` | Spanish text; the GUI translates it |
| `kind` | `bool`, `int`, `float`, `text`, `select`, `multiselect` |
| `default` | initial value; must be in `choices` if `select` |
| `help` | short explanation for the user |
| `choices` | `select`/`multiselect` options |
| `min`, `max`, `step` | number range |
| `advanced` | shown under "Advanced options" |
| `depends_on` | only used if another option is true |
| `unit` | visible unit (`nt`, `aa`, s) |

`Module.config(values)` validates types, ranges and unknown keys.

### `Module`

| field | use |
|---|---|
| `id` | unique, lowercase and `_` |
| `name`, `summary` | texts for the GUI |
| `order` | order in the selector |
| `requires` | binaries that must exist (`pyrodigal`, `diamond`…) |
| `gpu` | whether it needs a GPU |
| `experimental` | visible notice |
| `options` | GUI options |
| `run_fn` | `(RunContext, inputs) -> ModuleResult`; optional |
| `selftest_fn` | `() -> str`, fast, no downloads |
| `estimate_fn` | `(config) -> str` |
| `status_fn` | `(data_root) -> ModuleStatus` |

### `ModuleResult`

```python
ModuleResult(
    module="mi_modulo",
    ok=True,
    outputs={"gff": Path("...")},
    metrics={"genes": 137},
    warnings=["..."],
    error=None,
)
```

The GUI shows metrics, allows downloading small outputs and saves the
result in `st.session_state` so that another module can reuse it (for example,
Metabolism consumes the proteins from Annotation).

## Design rules

1. **No Streamlit in `modules/`.** The interface is generated from the contract.
2. **No downloads in `selftest()`.** It must run in CI in seconds.
3. **Explicit status.** If a binary or a datum is missing, `status()` says so.
4. **No effects from importing.** Importing a module must not compute anything.
5. **One module, one folder.** There is no central list to maintain.
6. **Declared outputs.** Every file it produces is returned in `outputs`.
7. **Honest estimated time.** If it depends on the genome, say so.

## Included modules

| id | name | options | GUI run |
|---|---|---|---|
| `annotation` | FASTA annotation | 10 | yes |
| `metabolism` | Metabolic network and gaps | 5 | yes (requires data) |
| `regulation` | Regulation and operons | 4 | ✅ per-genome summary |
| `structure` | Structure and filters | 6 | ✅ global summary |
| `mobile` | Mobile elements and defense | 4 | ✅ per-genome summary |
| `ghost` | Ghost genes | 3 | ✅ runs on FASTA |
| `taxonomy` | Taxonomy and ANI | 4 | ✅ per-genome summary |
| `membrane` | Membrane and transport | 4 | ✅ per-genome summary |
| `targets` | Essential targets | 3 | ✅ global summary |
| `validation` | Validation and benchmarks | 4 | ✅ global summary |

The summary modules read already-computed tables/JSON and return verifiable
counts; they **do not recompute** the pipeline. `ghost` does run the ORF engine
on an input FASTA. All of them can be tested with
`python -m modules --selftest <id>`.

## CLI

```bash
python -m modules --list
python -m modules --describe annotation
python -m modules --selftest annotation
python -m modules --selftest all
```

`--list` shows the status and whether it uses a GPU. `--describe` returns the
contract in JSON, useful for tests and for generating documentation.

## Tests

```bash
python -m pytest tests/test_modules.py -q
```

It checks:

- discovery and uniqueness of ids;
- unique options and valid defaults;
- range validation;
- real execution of Annotation on `webapp/samples/demo.fna`;
- fast self-test.

## Relationship with the `scripts/` reorganization

The module layer is the **stable contract**; scripts can be reorganized
behind it without touching the GUI. The folder plan is in
[`SCRIPTS_REORG_PLAN.md`](SCRIPTS_REORG_PLAN.md). The current modules call
`webapp.annotator` and `webapp.metabolism`; those that do not yet run from the
GUI show their status and their self-test, and will keep working when the
internal paths change.
