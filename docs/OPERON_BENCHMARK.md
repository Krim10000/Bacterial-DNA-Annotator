# Benchmark de operones por reglas

Herramienta: `tools/benchmark_operons.py`
Tests: `tests/test_benchmark_operons.py`
Informe: `reports/operon_benchmark.md` · JSON: `reports/operon_benchmark.json`

Este benchmark compara operones predichos por reglas contra un gold local
tipo RegulonDB. **No usa red, no descarga datos y no inventa gold**: si no
encuentra un gold aplicable, devuelve `not_evaluable` con la razón y las rutas
buscadas.

## Ejecución

```bash
# caso real local: K-12 (gold RegulonDB) contra reglas deterministas
.venv/bin/python tools/benchmark_operons.py \
  --gold data/regulondb/gt/gt_operons.tsv \
  --pred data/operon_finder/ecoli_k12_operons_rules.tsv \
  --genome ecoli_k12 \
  --json reports/operon_benchmark.json \
  --md reports/operon_benchmark.md

# autodetección: busca <genoma>.operons_auto.tsv, GFF y gold específico
.venv/bin/python tools/benchmark_operons.py --genome ecoli_k12

# tests sintéticos
.venv/bin/python -m pytest tests/test_benchmark_operons.py -q
```

## Entradas y autodetección

`--pred` acepta:

- TSV de operones: `operons_auto`, reglas, etc.
- GFF/GFF3 con features `operon` (o `operon_region`, `transcription_unit`,
  `regulon`) y `genes=`/`gene_ids=`; si faltan, se derivan de hijos con
  `Parent` al operón.
- GFF sin features de operón pero con `operon_id`/`tu_id` en los CDS.

`--gold` acepta los mismos formatos. En TSV se detectan alias de columnas:
`genes`/`gene_ids`/`gene_list`/`members`, `strand`/`sign`,
`start`/`begin`, `end`/`stop`, `contig`/`seqid`. Los genes se separan por
coma, punto y coma, barra o espacios. Si existe `mapping_ok`, solo se evalúan
las filas `yes` (convención de `gt_operons.tsv`).

Autodetección local (sin red), en este orden:

- Predicción: `results/<genoma>.operons_auto.tsv`,
  `out/<genoma>.operons_auto.tsv`, `out_fl/<genoma>.operons_auto.tsv`,
  `<genoma>.gff3`, `data/operon_finder/<genoma>_operons_rules.tsv`,
  `data/regulation/<genoma>/operons.tsv`.
- Gold: `data/regulation/<genoma>/gold_operons.tsv`,
  `data/benchmarks/<genoma>_gold_operons.tsv`,
  `data/gold/<genoma>_gold_operons.tsv` y, **solo si el genoma es K-12**,
  `data/regulondb/gt/gt_operons.tsv`.

Un gold genérico de K-12 no se aplica a otro genoma. Para un genoma nuevo hay
que aportar `--gold`; si no existe, el resultado es `not_evaluable`.

## Métricas

Se calculan dos vistas porque un operón es a la vez un conjunto de genes y un
intervalo genómico. Cada vista reporta precisión, recall, F1 y exact-match,
global y separada por hebra `+` / `-`.

### Nivel gen

- **Pares de genes**: todos los pares dentro de un operón predicado frente a
  los pares dentro de un operón gold. Es la métrica estándar de la literatura:
  un operón partido en dos pierde un par, no todos. `precisión = TP/pares
  predichos`, `recall = TP/pares gold`, `F1` armónica.
- **Conjuntos exactos (exact-set)**: el conjunto de genes del operón coincide
  exactamente con uno de gold. Se reporta P/R/F1 y la tasa de recuperación
  (`n_exact / n_gt`).

Esta vista requiere IDs de gen comparables. Si `operons_auto` usa `cds_N` y el
gold usa b-numbers, no se inventa una correspondencia: las métricas por gen
quedan `not_evaluable` y solo se reportan las de intervalo. Si se dispone de
`gt_gene_aliases.tsv` (columna `name`, `bnumber`), se aplica automáticamente
cuando está junto al gold; también se puede pasar `--alias`.

### Nivel intervalo

- **Solapamiento recíproco**: un intervalo predicho y uno gold se emparejan
  1:1 si ambas fracciones cubiertas son `>= --min-overlap` (0.5 por defecto):
  `min(intersección/pred, intersección/gold)`. Requiere misma hebra y contig
  compatible; un contig vacío en el gold se trata como comodín.
- **Intervalo exacto (exact-match)**: mismos `start`/`end` (y contig/hebra
  compatibles). Es muy sensible a los bordes y se informa aparte.

La coincidencia es 1:1 y determinista (mayor solapamiento, luego IDs). Un
mismo operón gold no puede "salvar" a dos predicciones.

### Filtros y conteos

- `--min-genes 2` (default): se evalúan operones con al menos 2 genes. Con
  `--min-genes 1` se incluyen los de un gen.
- `--min-overlap 0.5` (default): umbral de solapamiento recíproco.
- `--genome`: nombre para autodetección y reporte.
- `--json` / `--md`: rutas de salida (default `reports/operon_benchmark.*`).
- `--root`: raíz para autodetección (útil en tests; default el repo).

El JSON incluye `status`, `evaluation_level`, rutas, parámetros, conteos,
`by_strand` (`all`, `+`, `-`, y `unknown` si aparece), avisos y un bloque
`promoters` explícitamente `not_evaluable` en esta iteración.

## Resultado local actual (K-12)

Gold: `data/regulondb/gt/gt_operons.tsv` (2605 filas; 2581 con
`mapping_ok=yes`; 827 operones con ≥2 genes tras el filtro).
Predicción: `data/operon_finder/ecoli_k12_operons_rules.tsv` (2490 filas; 847
operones con ≥2 genes).

| hebra | pares de genes P | R | F1 | exact-set | intervalo P | R | F1 | intervalo exacto |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 0.828 | 0.770 | 0.798 | 0.565 | 0.761 | 0.780 | 0.771 | 0.552 |
| `+` | 0.868 | 0.740 | 0.799 | 0.569 | 0.754 | 0.795 | 0.774 | 0.561 |
| `-` | 0.796 | 0.799 | 0.797 | 0.561 | 0.769 | 0.765 | 0.767 | 0.543 |

Lectura: el nivel gen y el de intervalo coinciden en F1 ≈ 0.77–0.80. El
exact-match de conjuntos/intervalos es más duro (≈0.55): muchas predicciones
se solapan bien pero no reproducen exactamente los bordes o el conjunto
completo. La hebra `+` tiene mejor precisión por gen y la `-` mejor recall;
no se observa un sesgo de hebra extremo.

Avisos del caso real: 21 filas gold con `mapping_ok` distinto de `yes`
excluidas; 3 filas sin genes; 66 genes de la predicción fuera del universo
gold (excluidos de métricas por gen); el gold no declara contig y se trata como
comodín.

## Cuándo es `not_evaluable`

- No hay gold local aplicable (`gold_not_found`).
- El gold indicado no existe o no se puede leer (`gold_unreadable`).
- No hay predicción local (`prediction_not_found`).
- Gold/predicción no tienen operones tras `--min-genes`
  (`gold_without_evaluable_operons`, `prediction_without_operons`).
- Nivel gen sin IDs compartidos: el bucket queda `not_evaluable`, pero el
  nivel intervalo puede seguir evaluado (`evaluation_level="interval_only"`).

Un `not_evaluable` no es un cero: indica que faltan datos comparables.

## Qué NO mide

- **No es un holdout ciego**: las reglas de K-12 pueden haberse calibrado con
  el mismo gold. Mide consistencia contra el gold, no generalización a un
  genoma nuevo. Para generalización hace falta un holdout externo.
- **No evalúa promotores/TSS**: `promoters_auto.tsv` anota una ventana upstream
  de 50 pb, no un TSS experimental. Compararlo con `gt_promoters.tsv` requiere
  definir un emparejamiento TSS→ventana que aquí no se inventa; por eso el
  bloque `promoters` queda `not_evaluable`.
- **No corrige nomenclatura**: IDs no mapeados se excluyen de la vista por gen
  y eso puede favorecer la precisión; se reporta la cobertura de mapeo.
- **No valida expresión, terminadores, reguladores ni cotranscripción**.
- **No usa red ni bases externas**: solo archivos locales.

## Vista de intervalo y `operons_auto`

El benchmark es útil para un FASTA nuevo precisamente porque la vista de
intervalo no necesita IDs de gen: aunque `operons_auto` nombre `cds_1`,
`cds_2`, etc., las coordenadas se pueden comparar con un gold que use
b-numbers. La vista por gen solo se reporta si hay IDs/alias compartidos. Si no
hay gold local para ese genoma, el resultado correcto es `not_evaluable`, no
un número inventado.

## Regla

Un benchmark interno no es una métrica de producción. Reportar F1 sin declarar
gold, filtros, versión de reglas y si fue ciego es peor que no reportarlo.
