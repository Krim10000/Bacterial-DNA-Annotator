# Cobertura real de móviles, defensa, CRISPR, AMR y virulencia

Pregunta: de los supuestos 400 genomas del proyecto, ¿cuántos tienen realmente
datos de móviles/defensa/AMR, qué tipos hay y qué falta para llegar a 400?

La respuesta no se basa en el número 400 escrito a mano: `tools/mobile_coverage.py`
lee el inventario de disco y audita cada artefacto. El resultado reproducible
está en [`reports/mobile_coverage.json`](../reports/mobile_coverage.json) y
[`reports/mobile_coverage.md`](../reports/mobile_coverage.md).

## Respuesta corta

- El inventario real declara **400 genomas** (`data/mobile_typing/full_genomes.tsv`);
  la auditoría no asume ese número, lo lee del archivo.
- **400/400** tienen TSV tipado, GFF3, archivo de espaciadores y capa previa.
  Son **44.832 filas** tipadas.
- **310/400** tienen FASTA y por tanto CRISPR de novo evaluable. Los **90**
  restantes no tienen FASTA: sus filas de móviles/defensa/AMR/virulencia existen,
  pero CRISPR queda `not_available` (no evaluado). Por eso la cobertura
  **completa es 310/400** y **parcial 90/400**, no 400/400.
- No falta ningún TSV ni GFF del inventario. Lo que falta para que los 400 sean
  evaluables de punta a punta es el **FASTA de 90 genomas**.
- No se ejecutó ninguna base externa ni se creó/rellenó ninguna anotación.

## Fuentes auditadas

| fuente | qué aporta | cobertura observada |
|---|---|---:|
| `data/mobile_typing/full_genomes.tsv` | inventario (filo, organismo, contigs, FASTA) | 400 |
| `data/mobile_typing/full/*.tsv` | capa tipada (10 columnas) | 400 archivos · 44.832 filas |
| `data/mobile_typing/gff/full/*.gff3` | GFF3 por contig | 400 |
| `data/mobile_typing/crispr_spacers/full/*.tsv` | espaciadores | 400 archivos · 245 con filas |
| `data/mobile_typing/_work/full/` | intermedios del pipeline | 400 genomas |
| `data/mobile_defense/features.tsv` | capa previa (sin contig/hebra) | 400 genomas · 161.688 filas |
| `data/annotations_ab_initio/*/annotation.gff3` | entrada ab initio | 400 |
| `data/raw/*_genomic.fna.gz` | FASTA para CRISPR | 310 |
| `data/mobile_typing/summary_full.json` | resumen declarado | 400 procesados · 44.832 filas |
| `data/mobile_typing/validation_full.json` | validación declarada | `pass`, 0 errores |

El TSV del piloto (`data/mobile_typing/*.tsv`, 10 genomas) sigue existiendo y es
un subconjunto del barrido completo; la auditoría prioriza `full/`.

## Tipos agregados

Mapeo aplicado, sin inventar tipos: `transposon` → `mobile_element` y
`cas` → `crispr`. El resto se conserva tal cual.

| categoría | tipos incluidos | filas tipadas | genomas tipados | filas previas | genomas previos |
|---|---|---:|---:|---:|---:|
| `mobile_element` | `transposon`, `mobile_element` | 25.142 | 400 | 22.866 | 400 |
| `prophage` | `prophage` | 858 | 265 | 6.981 | 392 |
| `defense` | `defense` | 2.939 | 385 | 33.754 | 400 |
| `crispr` | `crispr`, `cas` | 5.037 | 296 | 3.279 | 247 |
| `amr` | `amr` | 1.051 | 271 | 81.542 | 400 |
| `virulence` | `virulence` | 9.805 | 394 | 13.266 | 398 |

Desglose CRISPR: **1.749 arrays** en 245 genomas + **3.288 filas `cas`** en 245
genomas; **20.687 espaciadores**. La unión de arrays y `cas` cubre 296 genomas.

La capa previa y la tipada **no miden lo mismo**. `features.tsv` viene de los
JSON previos por palabras clave y no tiene contig ni hebra; el TSV tipado filtra
esos ítems con nombre/Pfam/contexto y añade coordenadas. Los 81.542 ítems AMR
previos no son 81.542 determinantes: el barrido documentado en
`data/mobile_typing/README_full.md` los reduce a 1.051 tras descartar
transportadores genéricos. Esa diferencia es metodológica, no cobertura perdida.

## Pipeline que generó los datos

La atribución es estática: se deduce de `scripts/mobile/22_run_full.sh`,
`scripts/mobile/paths.py` y los nombres de salida. Los TSV no guardan metadatos
del generador y no se ejecutó ningún script para esta auditoría.

| etapa | script | salida | encontrado |
|---|---|---|---:|
| inventario | `scripts/mobile/19_prepare_full.py` | `mobile_typing/full_genomes.tsv` | 1 |
| extracción/mapeo | `scripts/mobile/01_extract_candidates.py --full` | `_work/full/*/candidates.tsv` | 400 |
| Pfam | `scripts/mobile/20_pfam_pool.py` | `_work/full/*/pfam_hits.tsv` | 400 |
| CRISPR | `scripts/mobile/21_crispr_pool.py` | `_work/full/*/crispr_arrays.tsv` | 310 |
| tipado | `scripts/mobile/05_classify.py --full` | `_work/full/*/elements.tsv` | 400 |
| escritura | `scripts/mobile/06_write_outputs.py --full` | `full/*.tsv`, `gff/full/*.gff3`, `crispr_spacers/full/*.tsv` | 400 / 400 / 400 |
| validación | `scripts/mobile/07_validate.py --full` | `validation_full.json` | 1 |

Orquestador: `bash scripts/mobile/22_run_full.sh all` (reanudable por etapas).
El piloto de 10 usa `scripts/mobile/run_all.sh`.

Los JSON previos tienen sus propios generadores, **fuera** de `scripts/mobile/`:

| sección | script | JSON |
|---|---|---|
| móviles | `scripts/annotation/22_annotate_transposons.py` | `data/mobile_elements.json` |
| defensa/CRISPR/profago | `scripts/annotation/20_annotate_defense.py` | `data/defense_systems.json` |
| AMR/virulencia | `scripts/annotation/21_annotate_amr_virulence.py` | `data/defense_virulence.json` |

No se encontró en el repositorio el script que aplana esos tres JSON en
`data/mobile_defense/features.tsv`; la auditoría usa el archivo existente tal
como está.

## Qué falta para 400

| brecha | n | lectura |
|---|---:|---|
| TSV tipado ausente | 0 | el inventario está cubierto |
| GFF3 ausente | 0 | idem |
| Espaciadores ausentes con FASTA | 0 | idem |
| **FASTA ausente** | **90** | CRISPR de novo no evaluable |
| CRISPR no evaluado | 90 | `not_available`, no “sin arrays” |
| Sin capa previa | 0 | `features.tsv` cubre 400 |
| Sin ninguna salida | 0 | todos tienen al menos TSV/GFF/previo |

Con los datos actuales, “llegar a 400” no requiere recalcular los 400: requiere
conseguir el FASTA de los 90 genomas y volver a correr las etapas dependientes de
secuencia (`22_run_full.sh crispr` y luego `build`). Eso no se hizo aquí: la
tarea prohíbe red y ejecución de bases/herramientas externas. La lista completa
de los 90 está en `reports/mobile_coverage.json` (`gaps.fasta_missing`).

### Brechas de método (no de archivos)

Aunque los 400 tuvieran FASTA, la capa actual sigue sin:

- DefenseFinder/PADLOC: sistemas de defensa sin nombre detectable.
- AMRFinderPlus/ResFinder: alelos y mutaciones puntuales de resistencia.
- ISEScan/ISfinder: extremos de IS y elemento completo (hoy se anota la transposasa).
- CRISPRCasFinder/minced: subtipo CRISPR inferido del array.
- MacSyFinder/VFDB: virulencia más allá del nombre del gen.

Estas herramientas no están instaladas/ejecutadas en este frente CPU. Son
límites de método, no genomas marcados como ausentes.

Deuda ya documentada en `data/mobile_typing/README_full.md`: 32 ítems previos
fuera de contig excluidos, 4.358 filas en coordenadas NCBI-only (sin CDS ab
initio) y 16 ítems CRISPR previos sin tipar. No se “arreglaron” ni se
reubicaron.

## `not_available`: qué significa exactamente

El JSON tiene una entrada por genoma con `status` (`available`/`not_available`),
`coverage` (`full`/`partial`/`none`) y un mapa `not_available` con la capa y el
motivo. Capas auditadas:

| capa | `not_available` cuando… |
|---|---|
| `typed_mobile` | no hay TSV o el TSV está vacío |
| `gff` | no hay GFF3 o está vacío |
| `crispr_scan` | no hay FASTA; o hay FASTA pero falta el TSV de espaciadores |
| `prior_features` | no hay filas para ese genoma en `features.tsv` |
| `fasta` | no existe `data/raw/<genoma>_genomic.fna.gz` |

`not_available` significa **no hay dato en disco**, no que la biología esté
ausente. El caso importante es CRISPR: en los 90 genomas sin FASTA no se puede
decir “sin arrays”; se dice “no evaluado”. Un genoma con FASTA y 0 espaciadores
sí aparece como `crispr_scan: available` con `rows: 0`.

`prior_features` es una capa auxiliar: si no hay filas para un genoma se reporta
en `not_available`, pero una cobertura con TSV + GFF + CRISPR evaluable sigue
siendo `full`. Las capas que determinan `full`/`partial` son `typed_mobile`,
`gff` y `crispr_scan`.

## Cobertura por genoma

En los 400 genomas del inventario:

- **310** tienen cobertura `full` (TSV + GFF + CRISPR evaluable).
- **90** tienen cobertura `partial`: TSV/GFF/previo presentes, CRISPR no
  evaluado por falta de FASTA; `not_available = {fasta, crispr_scan}`.
- **0** quedan `none` (sin TSV tipado).

El detalle por genoma y por filo está en
[`reports/mobile_coverage.md`](../reports/mobile_coverage.md); la versión
máquina-legible, en [`reports/mobile_coverage.json`](../reports/mobile_coverage.json).

## Reproducir

```bash
# Auditoría (solo lee; escribe únicamente los dos informes)
.venv/bin/python tools/mobile_coverage.py

# Tests con datos sintéticos y caso vacío
.venv/bin/python -m pytest tests/test_mobile_coverage.py -q
```

No se ejecutó hmmscan, DefenseFinder, CRISPRCasFinder, AMRFinderPlus ni ninguna
descarga. La auditoría contrasta, además, el `summary_full.json` declarado con el
recuento real: 400 genomas y 44.832 filas coinciden; `validation_full.json`
declara `pass` con 0 errores.

## Limitaciones de esta auditoría

- Verifica **presencia, tamaño y recuentos**, no la corrección biológica de cada
  fila; para eso está `scripts/mobile/07_validate.py` y la revisión manual.
- No comprueba frescura ni hashes contra las entradas; un artefacto podría ser
  más antiguo que el GFF/FASTA actual.
- La atribución a scripts es por convención de rutas, no por metadatos embebidos.
- Los recuentos de la capa previa incluyen ruido por palabras clave y no deben
  citarse como determinantes confirmados.
- No se infiere ausencia biológica a partir de `not_available`.
