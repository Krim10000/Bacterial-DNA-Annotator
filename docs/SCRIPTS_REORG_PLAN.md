# Plan de reorganización profunda de `scripts/`

**Estado:** **ejecutada** el 2026-11-22 en el árbol de trabajo. Este documento
se conserva como registro del plan y de las decisiones.

> Resultado: 27 → 17 carpetas, 23 movimientos/fusiones, **127 elementos movidos**
> y **548 referencias reescritas en 206 archivos**. La auditoría
> `tools/audit_script_paths.py` pasa con 615 referencias comprobadas; 37 tests
> pytest y 50 tests de `ghost` en verde. Se corrigieron además rutas que no
> detectaba el plan: `webapp/annotator.py` al paquete `ghost`, las rutas
> partidas `ROOT / "scripts" / "<carpeta>"` y dos `ANNOTATOR` que apuntaban a
> archivos inexistentes.

> La nueva capa `modules/` (ver [`MODULES.md`](MODULES.md)) desacopla la GUI de
> las rutas de `scripts/`. Eso reduce el riesgo de esta reorganización: los
> módulos seguirán llamando a `webapp.annotator`/`webapp.metabolism`, y las
> rutas solo cambian detrás de esa interfaz.

**Objetivo:** pasar de 27 carpetas de etapa a 17 carpetas por dominio, sin
romper importaciones, rutas de datos ni el pipeline, y dejando una migración
reversible.

> Regla que gobierna todo: **cada script debe seguir en
> `scripts/<una-carpeta>/<archivo>.py`**. Hay **130 scripts** que calculan la
> raíz con `parents[2]`. Añadir un nivel (por ejemplo
> `scripts/regulation/x.py`) rompería esas rutas. Por eso la
> reorganización renombra y fusiona carpetas, pero **no anida**.

## 1 · Punto de partida medido

| métrica | valor |
|---|---|
| scripts Python totales | 340 |
| scripts publicados (sin legado) | 301 |
| carpetas de etapa | 27 (más `__pycache__`, `_review`, `scripts_orginales`) |
| scripts con `parents[2]` | 130 |
| scripts con `parents[1]` | 29 |
| orquestadores `orchestrator.py` | 7 |
| scripts `run_*.sh` | 14 |
| referencias únicas a rutas `scripts/...` | 596 |
| referencias rotas detectadas | **12** (ver §2) |
| documentos con rutas `scripts/...` | 29 |
| consumidores de `12_operon_finder` por `sys.path` | 11 |
| destinos propuestos | 17 carpetas |

La auditoría se ejecuta con:

```bash
python tools/audit_script_paths.py
```

Hoy falla con 12 rutas. **No se debe empezar a mover hasta que pase**, porque si
no es imposible distinguir un fallo nuevo de uno viejo.

## 2 · Referencias rotas que hay que arreglar antes de mover

| referencia rota | dónde | corrección |
|---|---|---|
| `scripts/00_data_prep/14_parse_gff.py` | `PLAN_LEAD_ANNOTATION.md`, comentario en `01_annotation/14_parse_gff.py` | `scripts/01_annotation/14_parse_gff.py` |
| `scripts/03_oric/36_doric_lookup.py` | `DATA_NOTICE.md` | `scripts/03_oric/02_doric_lookup.py` |
| `scripts/15_allostery/00_build_ground_truth.py` | `docs/PLAN_ALOSTERIA.md` | nombre real del script o quitar la referencia |
| `scripts/32_ab_initio_annotator.py` | `01_annotation/00_annotate_single.py`, `01_annotate_batch.py` | `scripts/01_annotation/00_annotate_single.py` |
| `scripts/33_batch_ab_initio.py` | `01_annotation/01_annotate_batch.py` | `scripts/01_annotation/01_annotate_batch.py` |
| `scripts/34_benchmark_ab_initio.py` | `02_benchmark/01_compare_ncbi.py` | `scripts/02_benchmark/01_compare_ncbi.py` |
| `scripts/35a_marker_validation.py` | `03_oric/01_validate_markers.py` | `scripts/03_oric/01_validate_markers.py` |
| `scripts/36_doric_lookup.py` | `_review/`, `scripts_orginales/` | legado: se elimina, no se corrige |
| `scripts/37_functional_annotation.py` | comentarios de `05_metabolism/37*`, `38*`, legado | anteponer `scripts/05_metabolism/` |
| `scripts/38_batch_functional_annotation.py` | ídem | ídem |
| `scripts/39_rebuild_batch_summary.py` | ídem | ídem |
| `scripts/db/dl_pfam.sh` | `DATA_NOTICE.md` | ruta real (`scripts/db/` no existe) |

## 3 · Mapa destino

| origen | destino | scripts | acción | riesgo |
|---|---|---:|---|---|
| `00_data_prep/` | `data_prep/` | 9 | renombrar | bajo |
| `01_annotation/` | `annotation/` | 8 | renombrar + fusionar | medio |
| `04_defense/` | `annotation/` | 4 | fusionar | medio |
| `03_oric/` | `oric/` | 8 | renombrar | bajo |
| `05_metabolism/` | `metabolism/` | 29 | renombrar + fusionar | medio |
| `12_transport_t1/` | `metabolism/` | 6 | fusionar | bajo |
| `06_mobile/` | `mobile/` | 13 | renombrar | bajo |
| `07_membrane/` | `membrane/` | 7 | renombrar | bajo |
| `08_taxonomy/` | `taxonomy/` | 13 | renombrar | bajo |
| `06_structure_3d/` | `structure/` | 3 | fusionar | medio |
| `15_structure/` | `structure/` | 14 | renombrar + fusionar | medio |
| `15_allostery/` | `allostery/` | 30 | renombrar | medio |
| `10_regulation/` | `regulation/` | 7 | renombrar + fusionar | alto |
| `13_tss_tts/` | `regulation/` | 4 | fusionar | alto |
| `10_chains/` | `regulation/` | 2 | fusionar | bajo |
| `12_operon_finder/` | `regulation/` | 13 | **fusionar** | alto |
| `11_evidence/` | `evidence/` | 15 | renombrar + fusionar | medio |
| `13_reannotation/` | `evidence/` | 3 | fusionar | bajo |
| `11_validation/` | `validation/` | 25 | renombrar + fusionar | alto |
| `02_benchmark/` | `validation/` | 4 | fusionar | medio |
| `12_interpretability/` | `interpretability/` | 6 | renombrar | bajo |
| `09_ghost/` | `ghost/` | 4 | renombrar (conserva `ghost/` y `tests/`) | alto |
| `14_targets/` | `targets/` | 9 | renombrar | bajo |
| `lib/` | `lib/` | 5 | conservar | bajo |
| `scripts/tests/` | `tests/` | 2 | conservar | bajo |

Resultado: **17 carpetas** (15 de dominio, `lib` y `tests`). Se eliminan
`scripts_orginales/` y `_review/` del árbol de trabajo (ya están fuera del repo
publicado).

La simulación (`python tools/reorg_dry_run.py`) confirma:

- **23 movimientos/fusiones** (los otros cuatro orígenes son `lib`, `tests` y
  los dos de ML, que se conservan).
- **483 ocurrencias** de rutas a reescribir en **184 archivos**.
- **una sola colisión**: `annotation/orchestrator.py` (se resuelve renombrando
  el de defensa).

### Decisiones del mapa

- **`12_operon_finder` se fusiona con `regulation`.** Once scripts de
  `13_tss_tts` y `11_validation` inyectan su ruta en `sys.path`; se actualizan
  las once rutas y se valida que sigan importando. Es el punto de mayor riesgo
  de la migración y por eso tiene una puerta de validación propia.
- **`regulation` absorbe `12_operon_finder`, `13_tss_tts` y `10_chains`.** No
  hay colisiones de nombres de archivo; la única colisión es
  `orchestrator.py` en `annotation`, que se resuelve aparte.
- **`04_defense` se fusiona con `annotation`.** El `orchestrator.py` de defensa
  pasa a `orchestrator_defense.py`.
- **`06_structure_3d` y `15_structure` se fusionan.** `15_structure` no tiene
  `orchestrator.py`; no hay colisión.
- **El ML custom queda fuera de alcance de v1.0.** No se mueven `07_ml/`,
  `08_eval/`, los scripts ML de la raíz (`06_train_tokenizer_v2/v3.py`,
  `07_precompute_v2.py`, `09_train_model*.py`, `12_eval_probe*.py`,
  `13_comprehensive_eval.py`, `14_robust_eval.py`) ni `modal_train*.py`. Se
  quedan como están y se documentan como experimentales. La reorganización no
  los toca para no mezclar dos cambios de riesgo distinto.

## 4 · Fases, con puertas de validación

### Fase 0 — Congelar y medir (0,5 día)

```bash
git worktree add ../DNA_TOKEN_BAC_reorg -b reorg/scripts-v2
git tag pre-reorg-scripts
python tools/audit_script_paths.py > .scratch/audit_before.txt || true
python -m pytest tests -q > .scratch/tests_before.txt
```

Criterio: la rama nueva compila y los tests pasan **antes** de tocar nada.

### Fase 1 — Arreglar las 12 referencias rotas (0,5 día)

Se corrigen en la rama nueva, sin mover carpetas. Criterio:

```bash
python tools/audit_script_paths.py   # debe salir 0
```

### Fase 2 — Mover carpetas (0,5 día)

- `git mv scripts/<origen> scripts/<destino>` para las carpetas que no se
  fusionan.
- Para las fusiones, `git mv` archivo por archivo y resolver la colisión
  `orchestrator.py`.
- No tocar `data/`, `webapp/`, `run_pipeline.sh` todavía: primero se comprueba
  que los archivos existen en su nuevo sitio.

Criterio parcial: `find scripts -maxdepth 2 -name '*.py'` muestra 17 carpetas y
ningún archivo quedó en una ruta inesperada.

### Fase 3 — Reescribir referencias (0,5–1 día)

Script de migración (no manual) que reemplace **prefijos ordenados de más largo
a más corto**:

```text
scripts/00_data_prep/      -> scripts/data_prep/
scripts/01_annotation/     -> scripts/annotation/
scripts/04_defense/        -> scripts/annotation/
...
scripts/12_operon_finder/  -> scripts/regulation/
```

Incluye las formas:

- `scripts/<viejo>/...` en `.py`, `.sh`, `.md`, `.yml`, `Dockerfile`.
- `ROOT / "scripts" / "<viejo>"` en Python.
- `HERE="scripts/<viejo>"` en shell.
- `parents[2]` **no se toca**: la profundidad no cambia.

Archivos que hay que revisar a mano:

| archivo | qué cambia |
|---|---|
| `run_pipeline.sh` | ~20 rutas de etapa |
| `scripts/*/orchestrator.py` (7) | rutas internas y colisión de `orchestrator.py` |
| `scripts/06_mobile/run_all.sh` | `HERE` y 9 comandos |
| `scripts/08_taxonomy/run_all.sh`, `run_full.sh` | `HERE` |
| `scripts/09_ghost/run_one.sh`, `run_batch.sh`, `run_closure_domain.sh` | `SCRIPT_DIR` y rutas |
| `scripts/11_evidence/run_*.sh` (4) | `S="scripts/11_evidence"` |
| `scripts/regulation/run_rules_all.sh` | ruta del script |
| `scripts/14_targets/run_all.sh` | 9 rutas |
| `scripts/15_structure/run_all.sh` | 5 rutas |
| `scripts/regulation/*.py` (antes `13_tss_tts/`) | 4 `sys.path` a `12_operon_finder`, ahora `regulation` |
| `scripts/validation/*.py` | 7 `sys.path` a `12_operon_finder`, ahora `regulation` |
| `scripts/regulation/_common.py`, `14_validate.py` | `09_ghost` y `11_validation` |
| `scripts/regulation/*` | imports de `12_operon_finder`, ahora en la misma carpeta |
| `webapp/Dockerfile`, `.github/workflows/*` | `scripts/09_ghost` |
| `docs/*` (29 archivos), `README*`, `scripts/README.md` | rutas citadas |

Criterio:

```bash
python tools/audit_script_paths.py
python -m compileall -q scripts webapp
```

### Fase 4 — Validación funcional (1 día)

Puertas obligatorias, en este orden:

1. `python -m pytest tests -q` → 32/32.
2. `cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests` → 50/50.
3. `bash run_pipeline.sh doctor`.
4. `python -c` que importe `scripts/lib` y `scripts/ghost/ghost`.
5. Smoke de `--help` en cada entrypoint:

   ```bash
   python scripts/annotation/00_annotate_single.py --help
   python scripts/regulation/rule_operons.py --help
   python scripts/metabolism/04_network_trace.py --help
   python scripts/evidence/07_build_evidence.py --help
   python scripts/validation/regulondb_gt.py --help
   python scripts/ghost/ghost/run.py --help
   ```

6. **Equivalencia de salida** con la demo:

   ```bash
   # antes del cambio, en main
   python -m webapp.pipeline webapp/samples/demo.fna --outdir /tmp/before
   # después del cambio, en la rama nueva
   python -m webapp.pipeline webapp/samples/demo.fna --outdir /tmp/after
   diff <(sort /tmp/before/demo.gff3) <(sort /tmp/after/demo.gff3)
   diff /tmp/before/demo.summary.json /tmp/after/demo.summary.json
   ```

   Los GFF deben ser idénticos salvo rutas internas; el `summary.json` debe
   coincidir en conteos.

7. GUI: `python tools/gui_smoke.py http://127.0.0.1:8501` (demo, 3 idiomas,
   clic en gen).
8. CI verde en la rama.

### Fase 5 — Limpieza y documentación (0,5 día)

- Actualizar `scripts/README.md` con el mapa nuevo.
- Actualizar `docs/ARCHITECTURE.md`, `docs/INSTALL*`, `README*`.
- Borrar `scripts/_review/` y `scripts/scripts_orginales/` del árbol (el bundle
  de respaldo los conserva).
- Revisar `14_targets/`: no tiene `__main__`; decidir si se convierte en CLI o
  se marca como biblioteca.
- Revisar `05_metabolism/41_essential_pathways.py` (stub deprecado).
- Fusionar `11_evidence/09_confirm70.py` y `10_confirm_medium.py` (casi
  duplicados).

### Fase 6 — ML custom: fuera de alcance

No se toca en esta reorganización. `07_ml/`, `08_eval/`, los scripts ML de la
raíz de `scripts/` y `modal_train*.py` se quedan como están y se documentan como
experimentales. Si en el futuro se integran al producto, tendrán su propia
migración (cambian `parents[1]`, `sys.path` y rutas de checkpoints).

## 5 · Matriz de riesgo

| riesgo | probabilidad | impacto | mitigación |
|---|---|---|---|
| Romper `parents[2]` al anidar | media | alto | prohibido anidar; validar con `audit_script_paths.py` |
| Colisión `orchestrator.py` | alta | bajo | renombrar el de defensa |
| `12_operon_finder` en 11 `sys.path` | alta | alto | fusionar en `regulation` y actualizar las 11 rutas en el mismo commit |
| `09_ghost` en CI/Docker | media | alto | actualizar Dockerfile y workflows |
| Referencias en docs no detectadas | media | medio | auditoría + revisión de 29 docs |
| Salidas distintas por imports | baja | alto | comparación GFF/summary antes/después |
| Pérdida de scripts en `git mv` | baja | alto | rama + tag + bundle + `git status` |

## 6 · Rollback

- Rama `reorg/scripts-v2` y tag `pre-reorg-scripts`.
- Bundle de respaldo ya existe
  (`DNA_TOKEN_BAC_pre_v1_20261122_1328.bundle`).
- Si falla la validación: no se fusiona la rama; `main` queda intacto.
- Los datos (`data/`, `logs/`, `out*/`) no se mueven, así que el rollback no
  afecta a resultados.

## 7 · Estimación

| fase | tiempo |
|---|---|
| 0. Congelar y medir | 0,5 día |
| 1. Arreglar 12 referencias | 0,5 día |
| 2. Mover carpetas | 0,5 día |
| 3. Reescribir referencias | 0,5–1 día |
| 4. Validación funcional | 1 día |
| 5. Limpieza y documentación | 0,5 día |
| **Total** | **3–4 días** |

## 8 · Criterios de aceptación

- [ ] `tools/audit_script_paths.py` sale 0.
- [ ] `compileall` sin errores.
- [ ] 32 tests pytest y 50 unittest de ghost.
- [ ] `run_pipeline.sh doctor` funciona.
- [ ] `--help` de los entrypoints principales.
- [ ] GFF/summary de la demo equivalentes.
- [ ] GUI smoke con demo y clic.
- [ ] CI verde.
- [ ] `scripts/README.md` y docs actualizados.
- [ ] Sin `scripts_orginales/` ni `_review/` en el árbol de trabajo.
- [ ] El repo público sigue sin datos, modelos ni archivos >50 MB.

## 9 · Decisiones cerradas

1. **Nombres:** nombres de dominio (`annotation`, `metabolism`, `regulation`…).
   Los prefijos numéricos actuales no reflejan el orden real y se documentan en
   `scripts/README.md`.
2. **`12_operon_finder`:** se fusiona con `regulation`. Se asume el coste de
   actualizar las 11 rutas de `sys.path` y se valida con una puerta específica.
3. **ML custom:** fuera de alcance de v1.0. No se mueve ni se toca; queda
   documentado como experimental.
4. **Legado:** se borra `scripts_orginales/` y `_review/` del working tree. El
   bundle de respaldo los conserva.
5. **Compatibilidad:** sin symlinks. Un único cambio limpio; los symlinks
   complican Docker/Windows y la auditoría.

La migración se ejecutó sobre el árbol de trabajo con respaldo previo en
`.scratch/reorg_backup/` y un informe en `.scratch/reorg_applied.json`. No se
requiere una segunda ejecución.
