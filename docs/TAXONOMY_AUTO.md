# Taxonomía automática local

**Estado 2026-11-24:** capa implementada en `webapp/taxonomy_layer.py` y
probada en `tests/test_taxonomy_layer.py`. Normaliza la taxonomía
**precalculada** de un genoma a un TSV de cinco columnas. No consulta
NCBI/GTDB por red, no usa GPU y no ejecuta fastANI ni GTDB-Tk.

## Alcance

- Lee `taxonomy/full/<genome>.json` y, como respaldo,
  `taxonomy/<genome>.json`.
- Escribe `<genome>.taxonomy.tsv` con
  `taxon`, `rank`, `confidence`, `method`, `source`.
- Si no hay JSON local para el genoma, devuelve `status="not_available"` con
  una razón. **No inventa especie, género ni linaje.**
- No edita `webapp/pipeline.py`, `webapp/app.py` ni ningún artefacto del
  pipeline.

Para un FASTA nuevo sin JSON precalculado, esta capa no puede asignar
taxonomía: `not_available` es ausencia de dato, **no** evidencia de que el
genoma no tenga parientes. Automatizarlo exigiría fastANI/GTDB-Tk locales o
un panel de referencia, que no forman parte de esta tarea.

## Datos locales que detecta

`available(data_root)` devuelve `(True, "")` si existe al menos un JSON
taxonómico bajo `<data_root>/taxonomy/`. También acepta la raíz del proyecto
y resuelve `<raíz>/data/taxonomy/`.

`inventory(data_root)` informa además, sin crear ni modificar nada:

| grupo | rutas detectadas | papel |
|---|---|---|
| taxonomía principal | `taxonomy/full/<genome>.json` (400) | fuente del TSV |
| taxonomía piloto | `taxonomy/<genome>.json` (10) | respaldo |
| ANI | `taxonomy/ani_hits.tsv`, `taxonomy/full/ani_hits_full.tsv`, `taxonomy/full/ani_hits_full.jsonl`, `taxonomy/full/ani_groups.json`, `taxonomy/full/scope_plan.tsv` | inventario; no se recalcula |
| genomas/referencias | `selected.tsv`, `taxonomy/refpanel.tsv`, `taxonomy/full/genome_meta.tsv`, `taxonomy/full/panel_full.tsv`, `mobile_typing/pilot_genomes.tsv` | inventario; no se usan para inventar taxonomía |

Los mapas de ANI y de genomas se detectan y se reportan, pero no son
obligatorios para normalizar un JSON. Un JSON válido sin `ani` devuelve
`ani=null` en el resultado; no se convierte en 0 % ni en una especie nueva.

## API

```python
from pathlib import Path
from webapp.taxonomy_layer import available, inventory, run

ok, why = available(Path("data"))
info = inventory(Path("data"))
resultado = run(
    "GCF_000006685.1_ASM668v1",
    Path("data"),
    Path("out"),
)
```

### `available(data_root) -> tuple[bool, str]`

- `True, ""`: hay al menos un JSON taxonómico local.
- `False, razón`: no existe `taxonomy/`, no hay JSON o `data_root` no existe.

### `inventory(data_root) -> dict`

Devuelve `available`, `reason`, `data_root`, `taxonomy_dir`, `full_dir`,
`n_full`, `n_pilot`, `n_taxonomy`, `ani_maps` y `genome_maps`. Es un
inventario de solo lectura, pensado para auditar qué datos locales hay.

### `run(genome, data_root, outdir) -> dict`

Busca `taxonomy/full/<genome>.json`; si no existe, prueba
`taxonomy/<genome>.json`. Requiere que `genome` sea un nombre simple (sin
`/` ni `\`).

En éxito:

| campo | contenido |
|---|---|
| `status` | `"ok"` |
| `ok` | `True` |
| `rows` | filas normalizadas |
| `n_rows` | número de filas |
| `tsv` | ruta de `<outdir>/<genome>.taxonomy.tsv` |
| `source` / `source_rel` | JSON leído (absoluto / relativo a `data_root`) |
| `taxonomy_status` | `assigned`, `provisional` u otro estado del JSON |
| `confidence` | confianza declarada por el JSON (`unknown` si falta) |
| `species` | campo `species` del JSON; `None` si no existe |
| `closest_ref`, `ani` | copiados del JSON; pueden ser `null` |
| `method` | método del primer rango del linaje |
| `conflicts` | `evidence.conflicts` del JSON |
| `limitation` | límite interpretativo de la capa |

En fallo:

```python
{
    "status": "not_available",
    "ok": False,
    "reason": "...",
    "rows": [],
    "n_rows": 0,
    "tsv": "",
    "source": "",
}
```

No se escribe TSV cuando falta el genoma, el JSON es inválido/ilegible, no es
un objeto o no contiene ningún rango/especie utilizable.

### `normalize_taxonomy(data, source="") -> list[dict]`

Función pura para probar la normalización sin tocar disco. Cada fila tiene
exactamente los cinco campos del contrato.

## Normalización

### Rangos

Se aceptan los prefijos de las cadenas GTDB/NCBI:

| prefijo | `rank` |
|---|---|
| `d__` | `domain` |
| `p__` | `phylum` |
| `c__` | `class` |
| `o__` | `order` |
| `f__` | `family` |
| `g__` | `genus` |
| `s__` | `species` |
| `t__` | `strain` |

Un token vacío, sin prefijo, con prefijo desconocido o con valor `unknown` se
ignora. **No se adivina el rango por la posición** ni se rellenan huecos: un
rango ausente no aparece en el TSV.

### Preferencia de linaje

1. `gtdb_lineage` (o `evidence.gtdb.lineage`), si no es `unknown`.
2. `ncbi_lineage` (o `evidence.ncbi.lineage`).
3. `evidence.ncbi.lineage_ranks` / `evidence.ncbi.ranks` para JSON antiguos.

Se elige **un** linaje. No se fusionan NCBI y GTDB en una misma fila. Las
diferencias entre ambos siguen en el JSON original y se copian en
`resultado["conflicts"]`; el TSV no las oculta ni las resuelve.

### Especie

Una fila `species` solo se escribe si:

- el linaje elegido ya la incluye, o
- el JSON declara `species` (o `species_ncbi`/`species_gtdb`) y **el linaje
  elegido no trae especie**.

En el segundo caso `method="species_json"`: se normaliza un campo declarado,
no se deduce del género ni se consulta una base externa. Si no hay especie
declarada, no se escribe especie.

### Confianza y procedencia

- `confidence` se copia del JSON. La capa **no recalcula** high/medium/low;
  si falta, usa `unknown`.
- `method` vale `gtdb`, `ncbi` o `species_json`.
- `source` es la ruta del JSON relativa a `data_root`, por ejemplo
  `taxonomy/full/GCF_000006685.1_ASM668v1.json`.
- Las filas se ordenan `domain` → `phylum` → ... → `species` → `strain`.

## Salida

`<outdir>/<genome>.taxonomy.tsv`, UTF-8, separado por tabuladores:

```text
taxon	rank	confidence	method	source
Bacteria	domain	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydiota	phylum	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydiia	class	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydiales	order	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydiaceae	family	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydia	genus	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
Chlamydia muridarum	species	high	gtdb	taxonomy/full/GCF_000006685.1_ASM668v1.json
```

## Cobertura real (solo lectura, 2026-11-24)

- 410 registros taxonómicos locales: 400 en `full/` y 10 piloto.
- 400/400 genomas `full` normalizan con `status="ok"`.
- Linaje elegido: 349 GTDB r232; 51 sin GTDB, con NCBI + `species_json`.
- Filas por genoma: 391 con 7 rangos; 7 con 6; 1 con 5; 1 con 4. Los rangos
  faltantes no se rellenan (p. ej. un JSON NCBI sin `order`/`family` produce
  una tabla más corta).
- Estados del JSON: 382 `assigned`, 18 `provisional`.
- Confianza: 74 `high`, 308 `medium`, 18 `low`.

## Límites

1. **Taxonomía precalculada.** No hay ANI ni GTDB-Tk en esta capa. Un FASTA
   nuevo sin JSON local queda `not_available`; no se ejecuta nada por red.
2. **`not_available` no es evidencia negativa.** Solo dice que falta el JSON
   local.
3. **`ani=null` no es ANI 0 %.** Significa que no hubo comparación o
   referencia; la capa no lo interpreta.
4. **No valida la especie.** Copia la confianza declarada; no comprueba el
   umbral 95 % ni la simetría del ANI.
5. **Un solo linaje en el TSV.** La otra clasificación (NCBI o GTDB) no se
   fusiona; se conserva en el JSON y en `conflicts`.
6. **Versiones locales.** GTDB r232 y el taxdump/assembly summary usados por
   el JSON pueden no incluir especies recientes. La capa no actualiza bases.
7. **No infiere fenotipo.** No deduce patogenicidad, Gram, hospedador ni
   metabolismo a partir del nombre taxonómico.
8. **Sin red, sin GPU, sin binarios externos.** Solo `csv`, `json` y
   `pathlib`.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_taxonomy_layer.py -q
```

Cubren: detección de mapas locales, JSON válido y campos exactos, genoma
ausente, JSON inválido, especie declarada como respaldo y ausencia de
invención de especie cuando el JSON no la trae.
