# Capa sigma → operón → regulón

> **Estado:** v1.0.0 · capa derivada · CPU · sin red · sin GPU.
>
> **Qué hace:** integra la asignación sigma por operón
> (`data/sigma/<genoma>/sigma_promoters.tsv`) con la estructura de operones y,
> cuando existe, con las aristas TF→operón de
> `data/regulation/<genoma>/regulon.json`. Produce una tabla por genoma y una
> tabla global opcional.
>
> **Qué no hace:** no valida experimentalmente nada, no agrega dirección
> (activador/represor), no calcula fuerza de unión ni infiere TF que no estén
> en `regulon.json`. `tf_overlap` es un cruce de anotaciones, no una prueba de
> regulación.

---

## 1. Por qué existe

Los pasos anteriores producen tres capas separadas:

- `sigma_promoters.tsv`: una asignación de sigma por operón, con evidencia PWM,
  confianza y estado;
- `operons.tsv` / `operon_genes.tsv` / `genes.tsv`: qué genes forman cada
  operón y qué producto/EC tiene cada gen;
- `regulon.json`: aristas TF→operón ya calculadas por el constructor de
  regulones.

Faltaba una vista única que respondiera: **¿qué sigma tiene este operón, qué
genes/funciones arrastra y qué TFs del regulón lo tocan?** Esta capa resuelve
ese join sin re-calcular ni inventar evidencia.

---

## 2. Uso

### 2.1 Un genoma

```bash
.venv/bin/python scripts/sigma/40_build_sigma_regulon.py --genome ecoli_k12
```

Salida:

```text
data/sigma/regulon/ecoli_k12/sigma_regulon.tsv
data/sigma/regulon/ecoli_k12/sigma_regulon_summary.json
```

### 2.2 Todos los genomas procesables

```bash
.venv/bin/python scripts/sigma/40_build_sigma_regulon.py --all --jobs 8
```

Salida adicional:

```text
data/sigma/regulon/sigma_regulon_all.tsv
data/sigma/regulon/sigma_regulon_all_summary.json
```

`--all` descubre genomas que tienen `sigma_promoters.tsv` **y** las tres tablas
de regulación (`genes.tsv`, `operons.tsv`, `operon_genes.tsv`). En este checkout
son **401**: los 400 genomas GCF más `ecoli_k12`. Quedan fuera `subtilis` (no
tiene promotores sigma) y `ecoli_k12_test` (no tiene tablas de regulación).

### 2.3 Opciones útiles

| Opción | Uso |
|---|---|
| `--out DIR` | Directorio de salida. Default: `data/sigma/regulon`. |
| `--data-root DIR` | Raíz de datos de entrada. Default: `data`. |
| `--jobs N` | Procesos en paralelo para `--all`. Default: `1`. |
| `--limit N` / `--offset N` | Procesa un subconjunto de genomas. |
| `--selftest` | Autoprueba sin datos reales. |

Ejemplo de lote acotado para pruebas:

```bash
.venv/bin/python scripts/sigma/40_build_sigma_regulon.py \
  --all --limit 5 --out data/sigma/regulon/probe
```

### 2.4 Autoprueba y tests

```bash
.venv/bin/python scripts/sigma/40_build_sigma_regulon.py --selftest
.venv/bin/python -m pytest tests/test_sigma_regulon.py -q
```

---

## 3. Entradas

| Entrada | Ruta | Uso | ¿Obligatoria? |
|---|---|---|---|
| Asignación sigma | `data/sigma/<genoma>/sigma_promoters.tsv` | Una fila por operón con sigma, estado, confianza, evidencia y contexto | Sí |
| Operones | `data/regulation/<genoma>/operons.tsv` | `label`, `gene_names` de respaldo y total de operones | Sí |
| Genes del operón | `data/regulation/<genoma>/operon_genes.tsv` | Orden de genes y composición del operón | Sí |
| Genes | `data/regulation/<genoma>/genes.tsv` | Producto, EC y nombres de gen | Sí |
| Regulón | `data/regulation/<genoma>/regulon.json` | Aristas TF→operón para `tf_overlap`/`tf_names` | No |
| Gene map | `data/metabolism/genemaps/<genoma>_gene_map.tsv` | Rhea por gen (enriquecimiento opcional de `target_functions`) | No |

Si falta `regulon.json`, la fila queda con
`tf_overlap=not_available` y `tf_names=not_available`. Si el archivo existe pero
está corrupto, el genoma falla con un error explícito; no se degrada en
silencio a `not_available`.

### 3.1 Cómo se unen las tablas

1. `sigma_promoters.tsv` se deduplica por `operon_id`. Si hubiera más de una
   fila por operón, se conserva la de mayor prioridad:
   `assigned` > `ambiguous` > `unresolved`, luego `confidence_score` y `score`.
   En los datos actuales no hay duplicados.
2. `operon_genes.tsv` se ordena por `position` para obtener los genes del
   operón en orden.
3. `genes.tsv` aporta `product` y `ec` por `gene_id`.
4. `regulon.json` se indexa por `target_operon`. Si una arista no trae
   `target_operon`, se intenta mapear `target_genes` por `gene_id`, nombre de
   gen o locus tag, pero **solo si el alias apunta a un único operón**. Los
   alias ambiguos se descartan y se cuentan en `tf_stats`.
5. El gene map opcional se une por `gene_id` y, si no coincide, por
   `(contig, start, end)`. Los IDs de gene map de los genomas GCF suelen ser
   proteínas WP_ que no coinciden con los `gene_id` de `operon_genes.tsv`; por
   eso el respaldo por coordenadas es el que más aporta.

---

## 4. Esquema de `sigma_regulon.tsv`

Una fila por operón con asignación sigma. Todos los campos son texto en el TSV
(sin comillas CSV; los tabs/saltos de línea se normalizan a espacios).

| Columna | Tipo / valores | Descripción |
|---|---|---|
| `genome` | texto | Genoma procesado. |
| `operon_id` | texto | ID de operón usado en `operons.tsv`/`sigma_promoters.tsv`. |
| `operon_label` | texto | Etiqueta del operón; respaldo: `label` de `operons.tsv` o primer gen. |
| `sigma` | `sigma70`, `sigma54`, `sigma32`, `sigma38`, `sigma28`, `sigma24`, … | Sigma asignado por `sigma_promoters.tsv`. |
| `sigma_class` | `housekeeping`, `alternative`, `stress`, `stationary`, `flagellar`, `ecf`, `unknown`, … | Clase del sigma. Se copia de la entrada; el fallback solo se usa si falta. |
| `confidence_label` | `sospecha_fuerte`, `sospecha_media`, `sospecha_baja`, `sin_sospecha` | Confianza de la asignación sigma. |
| `sigma_status` | `assigned`, `ambiguous`, `unresolved` | Estado de la asignación sigma. |
| `target_genes` | texto | Nombres de gen del operón separados por coma, en orden de `position`. Respaldo: `gene_names` de `operons.tsv`. |
| `target_functions` | texto | `gen:producto [EC x] [RHEA:...]` separados por `; `. Máximo 8 entradas; si se recorta, agrega `(+N más)`. Si no hay producto/EC/Rhea, queda el nombre del gen. Si no hay ningún gen, `not_available`. |
| `tf_overlap` | entero ≥ 0 o `not_available` | **Número de TFs distintos** del regulón que regulan ese operón (no número de aristas). `0` = hay `regulon.json` pero ningún TF; `not_available` = no hay `regulon.json`. |
| `tf_names` | texto | Nombres de TF distintos separados por `; `, ordenados alfabéticamente. Vacío si `tf_overlap=0`; `not_available` si no hay regulón. |
| `pathway_context` | `yes: término; término`, `yes`, `no`, `not_available` | Contexto de `pathway_support`/`pathway_hits` de `sigma_promoters.tsv`. `yes:` incluye los términos; `not_available` si no hay dato. |
| `evidence` | texto | Evidencia original de `sigma_promoters.tsv` (PWM, margen, inventario, contexto, expresión). |
| `detected_by` | texto | Detectores originales (`promoter_pwm:...`, `inventory:...`, etc.). |
| `provenance` | JSON (texto) | Provenance original de `sigma_promoters.tsv` (contig, ventana, seed, versión). |

### 4.1 Ejemplo de fila

```text
ecoli_k12	ecoli_k12:NC_000913.3:1020410:-	sulA	sigma70	housekeeping	sospecha_baja	ambiguous	sulA	sulA:cell division inhibitor SulA	1	LexA	no	PWM sigma70 z=1.546 ...	promoter_pwm:sigma70;...	{"contig": "NC_000913.3", ...}
```

### 4.2 Ejemplo con contexto y varios genes

```text
ecoli_k12	ecoli_k12:NC_000913.3:1032139:+	hyaABCDEF	sigma54	alternative	sospecha_baja	ambiguous	hyaA,hyaB,hyaC,hyaD,hyaE,hyaF	hyaA:hydrogenase 1 small subunit [EC 1.12.5.1]; ...	0		yes: hyaA:hydrogenase; hyaB:hydrogenase; ...	PWM sigma54 ...	...
```

---

## 5. Esquema de `sigma_regulon_summary.json`

| Clave | Descripción |
|---|---|
| `genome` | Genoma. |
| `version` | Versión del constructor. |
| `n_edges` | Filas de la tabla = asignaciones sigma→operón emitidas. |
| `n_operons_total` | Operones en `operons.tsv`. |
| `n_operons_with_sigma` | Operones con fila en la tabla. |
| `n_operons_without_sigma` | `n_operons_total - n_operons_with_sigma`. |
| `n_sigma_rows_input` | Filas leídas de `sigma_promoters.tsv` (antes de deduplicar). |
| `n_duplicates_resolved` | Filas descartadas al deduplicar por operón. |
| `sigma_counts` | Conteo por valor de `sigma`. |
| `sigma_class_counts` | Conteo por `sigma_class`. |
| `confidence_counts` | Conteo por `confidence_label`. |
| `sigma_status_counts` | Conteo por `sigma_status`. |
| `pathway_context_counts` | `yes` / `no` / `not_available`. |
| `tf_available` | `true` si existe `regulon.json`. |
| `n_tfs` / `tfs` | TFs distintos que tocan al menos un operón de la tabla. |
| `n_operons_with_tf_overlap` | Filas con `tf_overlap > 0`. |
| `n_operons_sin_tf` | Filas con `tf_overlap == 0` (hay regulón, pero sin TF). |
| `n_operons_tf_not_available` | Filas con `tf_overlap == not_available`. |
| `tf_stats` | Aristas de entrada, mapeadas, sin mapear, ambiguas y sin TF. |
| `inputs` | SHA-256 de cada insumo, o `not_available` si falta. |

`sigma_regulon_all_summary.json` agrega los conteos de todos los genomas del
lote y agrega `n_genomes` y `failures`.

---

## 6. Resultados K-12 (verificados)

Comando ejecutado:

```bash
.venv/bin/python scripts/sigma/40_build_sigma_regulon.py --genome ecoli_k12
```

Salida:

```text
ecoli_k12: 2544 edges · 60 con TF · 2484 sin TF
```

Números reales del `sigma_regulon_summary.json`:

| Métrica | Valor |
|---|---:|
| `n_edges` | 2544 |
| `n_operons_total` | 2544 |
| `n_operons_with_sigma` | 2544 |
| `n_operons_without_sigma` | 0 |
| `n_duplicates_resolved` | 0 |
| `sigma38` | 745 |
| `sigma70` | 608 |
| `sigma28` | 384 |
| `sigma24` | 350 |
| `sigma32` | 244 |
| `sigma54` | 213 |
| `sospecha_baja` | 1098 |
| `sin_sospecha` | 992 |
| `sospecha_media` | 454 |
| `assigned` | 454 |
| `ambiguous` | 2033 |
| `unresolved` | 57 |
| `pathway_context=yes` | 165 |
| `pathway_context=no` | 2379 |
| `n_operons_with_tf_overlap` | 60 |
| `n_operons_sin_tf` | 2484 |
| `n_operons_tf_not_available` | 0 |
| `n_tfs` | 7 (`ArcA`, `CRP`, `FNR`, `IHF`, `LacI`, `LexA`, `TrpR`) |

Detalle del cruce con `regulon.json`:

- 65 aristas TF→operón de entrada;
- 63 operones destino distintos;
- 60 de esos operones tienen fila en `sigma_regulon.tsv`;
- 62 aristas caen en operones con fila: 58 operones con 1 TF y 2 operones con
  2 TFs;
- 3 aristas apuntan a operones que no tienen asignación sigma y no aparecen en
  la tabla;
- 0 aristas ambiguas o sin mapear.

`ecoli_k12` no tiene `gene_map.tsv`, así que `target_functions` usa solo
producto/EC de `genes.tsv` (sin Rhea). El campo `inputs.gene_map.tsv` queda
`not_available`.

---

## 7. Límites

1. **Una asignación por operón.** La tabla hereda la decisión de
   `sigma_promoters.tsv`: un solo sigma por operón. No representa promotores
   alternativos ni sigma secundarios del mismo operón.
2. **No hay dirección ni signo.** `regulon.json` puede tener `sign`, pero esta
   capa no lo propaga: `tf_names` solo indica que el TF aparece como regulador,
   no si activa o reprime.
3. **No hay validación experimental.** `confidence_label` y `sigma_status` son
   la confianza del scanner PWM, no evidencia de laboratorio.
4. **`regulon.json` es disperso.** En K-12 cubre 65 aristas y 7 TFs; por eso
   solo 60 de 2544 operones tienen `tf_overlap > 0`. No debe leerse como que el
   resto no está regulado.
5. **No todos los genomas tienen regulón.** Si no existe `regulon.json`, el
   valor es `not_available`; no se completa con ortología ni con RegulonDB.
6. **El gene map es opcional y no siempre une por ID.** En genomas GCF los
   `gene_id` del gene map (WP_) no coinciden con `operon_genes.tsv`; se usa el
   respaldo por coordenadas. Si no hay match, no se agrega Rhea.
7. **`target_functions` se acota.** Máximo 8 entradas y máximo 2 IDs Rhea por
   gen para que la tabla global de 401 genomas sea manejable. El recorte se
   declara con `(+N más)`.
8. **Sin expresión.** Esta capa no incorpora `expression_rules.tsv` ni tablas de
   expresión; si se necesita, es un join posterior.
9. **`n_edges` no es el número de aristas del regulón.** Es el número de filas
   sigma→operón. El número de aristas TF se reporta en `tf_stats`.
10. **No reemplaza a `regulon.json`.** Es una vista integrada; la fuente de
    verdad de TF→operón sigue siendo `data/regulation/<genoma>/regulon.json`.

---

## 8. Reproducibilidad

- El script no usa red, no descarga datos y no escribe fuera de `--out`.
- `sigma_regulon_summary.json` guarda SHA-256 de cada insumo usado.
- El orden de filas por genoma sigue el orden de `sigma_promoters.tsv`; la tabla
  global ordena por nombre de genoma.
- Los tests puros están en `tests/test_sigma_regulon.py` y no dependen de los
  400 genomas.

---

## 9. Relación con otros documentos

- [`SIGMA_MODULE.md`](SIGMA_MODULE.md): cómo se calcula la asignación sigma.
- [`SIGMA_GROUPS.md`](SIGMA_GROUPS.md): agrupación de factores sigma por
  ortólogo.
- [`PLAN_REGULON_BUILDER.md`](../PLAN_REGULON_BUILDER.md): construcción de
  `regulon.json`.
