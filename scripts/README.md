# Scripts del pipeline

El árbol se organizó por **dominio**, no por número de carpeta. La regla que
gobierna todo: cada script vive en `scripts/<una-carpeta>/<archivo>.py`, a dos
niveles de la raíz. Hay **130 scripts** que calculan la raíz con `parents[2]`;
añadir un nivel los rompería.

> Reorganización ejecutada en noviembre de 2026: 27 → 17 carpetas, 23
> movimientos/fusiones y 548 referencias reescritas. Ver
> [`docs/SCRIPTS_REORG_PLAN.md`](../docs/SCRIPTS_REORG_PLAN.md) para el detalle
> y el informe de la simulación.

## Mapa actual

| carpeta | contenido | entrypoints principales | tiempo | GPU |
|---|---|---|---|---|
| `data_prep/` | catálogo NCBI, filtrado, muestreo, corpus | `00_download_catalog.py`, `05_build_corpus.py` | 2–5 min | no |
| `annotation/` | anotación ab initio, ARN, arquitectura, defensa | `00_annotate_single.py`, `01_annotate_batch.py`, `14_parse_gff.py` | ~3 min/genoma | no |
| `oric/` | origen de replicación (ORCA, DoriC) | `35b_orca_from_gbk.py`, `40_triage_anomalous.py` | ~1 min/genoma | no |
| `metabolism/` | red metabólica, Rhea, brechas, transporte T1 | `04_network_trace.py`, `09_make_gene_map.py` | ~2 min/genoma | no |
| `mobile/` | elementos móviles, CRISPR, AMR, profagos | `01_extract_candidates.py`, `05_classify.py` | ~25 min/barrido | no |
| `membrane/` | membrana, transporte y secreción | `03_annotate_membrane.py`, `06_run_full.py` | ~5 min/genoma | no |
| `taxonomy/` | ANI, GTDB y especie | `07_assign_taxonomy_full.py` | 1–5 min/genoma | no |
| `structure/` | afdb, Foldseek, AlphaFill, ESMFold, filtros F1–F5 | `06_map_active_sites.py`, `09_consensus.py` | 1–5 min/proteína | API/opcional |
| `allostery/` | alostería, Boltz-2, ANM, OpenMM | `18_genome_pipeline.py`, `22_campaign.py` | horas | sí (Modal) |
| `regulation/` | operones, promotores, TSS/TTS, cadenas, reglas, horquillas | `rule_operons.py`, `11_predict_operons.py`, `pipeline.py` | **8 s** (reglas) | no |
| `evidence/` | evidencia por brecha, re-anotación | `07_build_evidence.py`, `build_reannotation.py` | 10–60 min | no |
| `validation/` | RegulonDB, DBTBS, benchmarks, sintenía | `regulondb_gt.py`, `benchmark_operons.py` | 2–10 min | no |
| `interpretability/` | embeddings y firmas | `01_extract_embeddings.py` | 10–40 min | recomendada |
| `ghost/` | genes fantasma, desiertos, cierres por dominio | `ghost/run.py`, `scan_domain_closure.py` | ~5 s/genoma | no |
| `targets/` | targets esenciales P1/P2/P3 | `00_prepare.py` (sin `__main__`) | n/a | no |
| `lib/` | utilidades compartidas | `io.py`, `dna_v2.py` | — | — |
| `tests/` | pruebas sueltas | — | — | — |

Fuera del alcance de la reorganización de v1.0: `07_ml/`, `08_eval/` y los
scripts ML de la raíz de `scripts/` (modelo custom, experimental). No se movieron
para no mezclar dos cambios de riesgo distinto.

## Cómo ejecutar

```bash
python scripts/00_doctor.py
bash run_pipeline.sh components
python scripts/regulation/rule_operons.py --help
python scripts/metabolism/04_network_trace.py --help
```

Cada etapa escribe en `data/<módulo>/` y `logs/`, nunca en el código. Los
datos grandes están fuera de git y se regeneran.

## Convenciones

1. **Un entrypoint por carpeta** con `if __name__ == "__main__"`. Si un script
   no lo tiene, es una biblioteca.
2. **Rutas relativas a la raíz** con `parents[2]`; nunca rutas absolutas.
3. **Argumentos por CLI** (`--data-root`, `--out`, `--threads`), no constantes
   globales.
4. **Salidas con `evidence`, `confidence` y `status`.** Nada se marca
   `confirmed` sin EC, contexto o experimento.
5. **No aplicar correcciones automáticamente**: las sugerencias se marcan
   `suggested`.
6. **No cambiar la profundidad** de los scripts. Si hace falta una subcarpeta,
   primero hay que migrar los `parents[2]`; no se hace en v1.0.

## Validación

```bash
python tools/audit_script_paths.py
python -m compileall -q scripts webapp modules
python -m pytest tests -q
cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests
bash run_pipeline.sh doctor
```

`audit_script_paths.py` revisa tanto `scripts/<ruta>` literales como las rutas
construidas con `"scripts" / "<carpeta>"`.

## Scripts de legado

Las copias anteriores a la reorganización no forman parte del árbol actual ni
de la versión publicada. Se conservan solo en el respaldo local del autor
(`.scratch/reorg_backup/`); el bundle de git también las contiene.

## Orden recomendado para entender el código

1. `ghost/ghost/run.py` — el flujo más pequeño y completo.
2. `regulation/rule_operons.py` — reglas deterministas.
3. `metabolism/04_network_trace.py` — núcleo de brechas.
4. `evidence/07_build_evidence.py` — contrato de evidencia.
5. `webapp/annotator.py` — lo que ve el usuario al subir un FASTA.
