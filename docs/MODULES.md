# Módulos autodescubiertos

🌐 **Español** · [English](MODULES.en.md)

La aplicación tiene dos capas:

1. **Núcleo** (`webapp/`): anotación, navegador, datos, i18n.
2. **Módulos** (`modules/`): cada dominio declara sus opciones, requisitos,
   tiempo, estado, autoprueba y —si aplica— una función de ejecución.

La GUI **no tiene una pestaña escrita a mano por módulo**. La vista
**🧩 Módulos** descubre las carpetas y construye sola:

- el selector de módulo,
- las casillas, números, textos y selectores,
- el estado (`✅ disponible` / `⚠️ incompleto`) y el motivo,
- el tiempo estimado,
- el botón **🧪 Probar módulo**,
- el botón **▶️ Ejecutar módulo** si el módulo lo implementa.

## Añadir un módulo

```text
modules/
  mi_modulo/
    __init__.py
    module.py     # exporta MODULE
```

`module.py` mínimo:

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

Al guardar la carpeta:

```bash
python -m modules --list
python -m modules --selftest mi_modulo
```

Aparece en la GUI sin tocar `app.py`.

## Contrato

### `Option`

| campo | uso |
|---|---|
| `key` | identificador único dentro del módulo |
| `label` | texto en español; la GUI lo traduce |
| `kind` | `bool`, `int`, `float`, `text`, `select`, `multiselect` |
| `default` | valor inicial; debe estar en `choices` si es `select` |
| `help` | explicación breve para el usuario |
| `choices` | opciones de `select`/`multiselect` |
| `min`, `max`, `step` | rango de números |
| `advanced` | se muestra en «Opciones avanzadas» |
| `depends_on` | solo se usa si otra opción es verdadera |
| `unit` | unidad visible (`nt`, `aa`, s) |

`Module.config(values)` valida tipos, rangos y claves desconocidas.

### `Module`

| campo | uso |
|---|---|
| `id` | único, minúsculas y `_` |
| `name`, `summary` | textos para la GUI |
| `order` | orden en el selector |
| `requires` | binarios que deben existir (`pyrodigal`, `diamond`…) |
| `gpu` | si necesita GPU |
| `experimental` | aviso visible |
| `options` | opciones de la GUI |
| `run_fn` | `(RunContext, inputs) -> ModuleResult`; opcional |
| `selftest_fn` | `() -> str`, rápido, sin descargas |
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

La GUI muestra métricas, permite descargar salidas pequeñas y guarda el
resultado en `st.session_state` para que otro módulo lo reutilice (por ejemplo,
Metabolismo consume las proteínas de Anotación).

## Reglas de diseño

1. **Nada de Streamlit en `modules/`.** La interfaz se genera desde el contrato.
2. **Sin descargas en `selftest()`.** Debe correr en CI en segundos.
3. **Estado explícito.** Si falta un binario o un dato, `status()` lo dice.
4. **Sin efectos por importar.** Importar un módulo no debe calcular nada.
5. **Un módulo, una carpeta.** No hay lista central que mantener.
6. **Salidas declaradas.** Todo archivo que produce se devuelve en `outputs`.
7. **Tiempo estimado honesto.** Si depende del genoma, dilo.

## Módulos incluidos

| id | nombre | opciones | ejecución GUI |
|---|---|---|---|
| `annotation` | Anotación de FASTA | 10 | sí |
| `metabolism` | Red metabólica y brechas | 5 | sí (requiere datos) |
| `regulation` | Regulación y operones | 4 | ✅ resumen por genoma |
| `structure` | Estructura y filtros | 6 | ✅ resumen global |
| `mobile` | Móviles y defensa | 4 | ✅ resumen por genoma |
| `ghost` | Genes fantasma | 3 | ✅ ejecuta sobre FASTA |
| `taxonomy` | Taxonomía y ANI | 4 | ✅ resumen por genoma |
| `membrane` | Membrana y transporte | 4 | ✅ resumen por genoma |
| `targets` | Targets esenciales | 3 | ✅ resumen global |
| `validation` | Validación y benchmarks | 4 | ✅ resumen global |

Los módulos de resumen leen tablas/JSON ya calculados y devuelven conteos
verificables; **no recalculan** el pipeline. `ghost` sí ejecuta el motor de ORFs
sobre un FASTA de entrada. Todos se pueden probar con
`python -m modules --selftest <id>`.

## CLI

```bash
python -m modules --list
python -m modules --describe annotation
python -m modules --selftest annotation
python -m modules --selftest all
```

`--list` muestra el estado y si usa GPU. `--describe` devuelve el contrato en
JSON, útil para tests y para generar documentación.

## Pruebas

```bash
python -m pytest tests/test_modules.py -q
```

Comprueba:

- descubrimiento y unicidad de ids;
- opciones únicas y defaults válidos;
- validación de rangos;
- ejecución real de Anotación sobre `webapp/samples/demo.fna`;
- autoprueba rápida.

## Relación con la reorganización de `scripts/`

La capa de módulos es el **contrato estable**; los scripts pueden reorganizarse
detrás sin tocar la GUI. El plan de carpetas está en
[`SCRIPTS_REORG_PLAN.md`](SCRIPTS_REORG_PLAN.md). Los módulos actuales llaman a
`webapp.annotator` y `webapp.metabolism`; los que aún no ejecutan desde la GUI
muestran su estado y su autoprueba, y seguirán funcionando cuando cambien las
rutas internas.
