# Handoff: integración del anotador de sigma al pipeline

> **Para el agente que mantiene el pipeline.** Este documento dice qué copiar,
> qué mergear, qué no tocar y cómo verificar que quedó integrado.
>
> Módulo: **`sigma`** · versión **1.0.0** · experimental · CPU · sin GPU · sin red.

---

## 1. Qué se integró

| Pieza | Ruta | Acción |
|---|---|---|
| Núcleo | `scripts/sigma/sigma_annotator.py` | Copiar |
| CLI anotación | `scripts/sigma/10_annotate_sigma.py` | Copiar |
| CLI agrupación | `scripts/sigma/20_group_sigma.py` | Copiar |
| Módulo GUI | `modules/sigma/module.py` | Copiar (autodescubierto) |
| Tests | `tests/test_sigma_annotator.py` | Copiar |
| Documentación | `docs/SIGMA_MODULE.md`, `docs/SIGMA_GROUPS.md`, `docs/SIGMA_AUDIT.md` | Copiar |
| Preview GUI | `webapp/module_ui.py` | **Mergear**, no sobrescribir |
| Pista navegador | `webapp/browser.py` (pista `sigma`) | **Mergear**, no sobrescribir |
| Vista Sigma | `webapp/sigma_view.py` | Copiar |
| Modo/atajo GUI | `webapp/app.py` (`sigma` + `_go_to_sigma`) | **Mergear**, no sobrescribir |
| Etapa CLI | `run_pipeline.sh` (`sigma`) | **Mergear** |
| Tabla de módulos | `docs/MODULES.md` | Añadir fila |
| Traducciones | `webapp/locales/_source.json`, `en.json`, `zh.json` | Añadir claves nuevas |

---

## 2. Contrato de datos

El módulo **lee**:

```text
data/regulation/<genoma>/genes.tsv
data/regulation/<genoma>/operons.tsv
data/regulation/<genoma>/operon_genes.tsv
data/regulation/<genoma>/regulators_inventory.tsv
data/regulondb/gt/promoter_pwms.tsv
data/annotations_ab_initio/<genoma>/proteins.faa
data/expression/<genoma>.tsv        # opcional
```

El módulo **escribe**:

```text
data/sigma/<genoma>/sigma_factors.tsv
data/sigma/<genoma>/sigma_promoters.tsv
data/sigma/<genoma>/sigma_summary.json
data/sigma/<genoma>/sigma_report.md
data/sigma/sigma_orthologs.tsv      # global, tras 20_group_sigma.py
data/sigma/summary_all.tsv          # global, tras --all
```

**No editar a mano** `data/sigma/`; se regenera. **No tocar**
`data/regulation/`; el módulo es de solo lectura ahí.

---

## 3. Campos que no se pueden perder

### `sigma_promoters.tsv`

- `sigma`, `sigma_class`, `sigma_arch`, `spacer`
- `score`, `zscore`, `pvalue`, `score_margin`
- `sigma_status` (`assigned` / `ambiguous` / `unresolved`)
- `confidence` (`high` / `medium` / `low`)
- `confidence_label` (`sospecha_fuerte` / `sospecha_media` / `sospecha_baja` / `sin_sospecha`)
- `evidence`
- `detected_by`
- `pipeline_step`
- `alternatives`
- `pathway_support`, `pathway_hits`
- `expression_support`, `expression_detail`
- `inventory_support`, `inventory_detail`
- `location_support`
- `provenance`

### `sigma_factors.tsv`

- `sigma_type`, `sigma_class`
- `status`, `confidence`, `confidence_label`
- `evidence`, `detected_by`, `pipeline_step`
- `countable`, `quality_flag`
- `ortholog_group`, `copy_number`, `is_expansion`
- `provenance`

**Regla:** si se vuelve a correr `10_annotate_sigma.py`, no se debe perder
`ortholog_group`; el núcleo ya lo preserva. Si se cambia el orden de columnas,
actualizar `OUTPUT_FACTOR_COLS` y `OUTPUT_PROMOTER_COLS`.

---

## 4. GUI

### 4.1 Dónde aparece

No va en **“1 · Herramientas disponibles”**: esa lista es de binarios externos.
El anotador es un **módulo** y aparece en:

```text
🧭 Explorar genoma → 🧪 Sigmas          # ver resultados
🧩 Módulos → Sigmas y promotores        # ejecutar el módulo
```

Enlace directo: `?modo=sigmas&genoma=ecoli_k12`.

Además se agregó un atajo:

```text
Más → 🧪 Anotador de sigmas
```

### 4.2 Qué debe mostrar

- Pista **σ promotores** en el navegador (marca + banda, color por tipo,
  opacidad por confianza).
- Vista **🧭 Explorar genoma → 🧪 Sigmas** con métricas, filtros, tabla,
  evidencia y descargas.
- Tabla de `sigma_promoters.tsv`.
- Selector de fila.
- JSON completo de la fila (evidencia).
- Descargas TSV/JSON/Markdown.

El preview se implementó en `webapp/module_ui.py` con la convención:

```python
outputs["preview"] = ruta_tsv
```

Si el pipeline tiene otra UI de módulos, replicar esa convención: tabla + fila
seleccionada + descarga, sin duplicar el archivo `preview` en la lista de
descargas.

---

## 5. Integración en el pipeline

`run_pipeline.sh` debe tener:

```bash
sigma)
  "$PY" scripts/sigma/10_annotate_sigma.py --all --jobs "${SIGMA_JOBS:-4}"
  "$PY" scripts/sigma/20_group_sigma.py --all --threads "${SIGMA_JOBS:-4}"
  ;;
```

Uso:

```bash
bash run_pipeline.sh sigma
```

No agregar este paso a `all` sin decisión explícita: el lote completo tarda
~20–30 min y escribe ~700 MB.

---

## 6. Verificación

```bash
# 1. El módulo se descubre
python -m modules --list | grep sigma

# 2. Autoprueba
python -m modules --selftest sigma
python scripts/sigma/10_annotate_sigma.py --selftest
python scripts/sigma/20_group_sigma.py --selftest

# 3. Un genoma
python scripts/sigma/10_annotate_sigma.py --genome ecoli_k12 --null-windows 80

# 4. Agrupación global
python scripts/sigma/20_group_sigma.py --all --threads 8

# 5. Tests
python -m pytest tests/test_sigma_annotator.py tests/test_modules.py \
  tests/test_menu_nav.py tests/test_i18n.py -q
```

### Benchmark y capas (2026-11-25)

Evaluación fina: [`SIGMA_LAYERED_EVAL.md`](SIGMA_LAYERED_EVAL.md).

| Métrica | K-12 | B. subtilis |
|---|---:|---:|
| Top-1 | 57,0 % | 44,8 % |
| Mayoría trivial | 62,6 % | 87,6 % |
| Top-3 cualquiera | 90,2 % | 85,5 % |
| TSS multi-sigma | 124 | 0 |

**No presentar el módulo como validado.** El top-1 no supera la mayoría
trivial; la confianza ordena hipótesis pero no está calibrada. Ver
`docs/SIGMA_BENCHMARK.md`.

### Valores de referencia del lote (2026-11-25)

| Métrica | Valor |
|---|---:|
| Genomas | 401 |
| Factores | 5.993 |
| Contables | 2.135 |
| `sigma_unresolved` | 3.765 |
| Grupos ortólogos | 2.432 |
| Grupos contables | 663 |
| Expansiones | 213 |

Los números pueden cambiar con `window`, `null_windows`, `max_p`, `min_z` o
`min_margin`; el reporte debe declarar los parámetros usados.

---

## 7. Lo que NO debe hacer el agente

- No etiquetar `sigma70` por defecto.
- No contar `sigma_unresolved` como familia.
- No borrar `ortholog_group` al reanotar.
- No tocar `data/regulation/`.
- No meter el módulo en “Herramientas disponibles” como si fuera un binario.
- No prometer dirección/magnitud de efector.
- No llamar `confirmed` a una asignación.
- No commitear `data/sigma/` (está bajo `data/**`, ignorado por git).

---

## 8. Mensaje corto para el agente

> Integra el módulo `sigma` desde `modules/sigma/` y `scripts/sigma/`. Añade la
> etapa `sigma` a `run_pipeline.sh`. Mergea el preview de `webapp/module_ui.py`
> y el atajo `_go_to_sigma` de `webapp/app.py`; no sobrescribas esos archivos.
> Verifica con `python -m modules --selftest sigma` y
> `pytest tests/test_sigma_annotator.py`. No toques `data/regulation/`. El
> módulo debe conservar `evidence`, `detected_by`, `pipeline_step`, `provenance`,
> `countable`, `quality_flag` y `ortholog_group` en cada fila. No lo pongas en
> “Herramientas disponibles”: es un módulo, va en 🧩 Módulos.
