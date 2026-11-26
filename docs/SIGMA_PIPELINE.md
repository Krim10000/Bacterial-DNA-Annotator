# Pipeline de sigma: orden canónico y responsabilidades

## Regla

Ninguna capa vuelve a hacer el trabajo de otra. Cada paso consume la salida
anterior y deja un archivo con contrato explícito.

## Orden

| # | Capa | Script | Entrada | Salida | No repetir |
|---:|---|---|---|---|---|
| 0 | Anotación base | `scripts/annotation/00_annotate_single.py` o pipeline | FASTA | GFF, proteínas, genes | No re-predice CDS |
| 1 | Operones/promotores de operón | `scripts/regulation/find_operons.py` / pipeline | GFF | `operons.tsv`, `operon_genes.tsv`, `promoters.tsv` | No re-llama operones |
| 2 | Inventario sigma | `scripts/regulation/12_promoters_sigma.py` | genes/proteoma + HMMER/DIAMOND | `regulators_inventory.tsv` | No usar su `promoters.tsv` como capa final |
| 3 | Ortología global | `scripts/sigma/20_group_sigma.py` | inventario | `ortholog_group` | No recalcular por genoma |
| 4 | Escaneo sigma | `scripts/sigma/10_annotate_sigma.py` / `annotate_genome` | operones + secuencia + inventario + PWM | `data/sigma/<g>/sigma_factors.tsv`, `sigma_promoters.tsv` | No reescanear con `12` |
| 5 | Calibración | `63_train_kmer_calibrator.py` (offline) y política `hybrid` | gold de una especie | modelo JSON | No calibrar con gold en producción |
| 6 | Contexto | `51_eval_context.py` (offline/opcional) | PWM + features | contexto/ablación | No reanotar genes/operones |
| 7 | Dossier | `scripts/sigma/02_build_dossier.py` | salidas 2–4 | `reports/sigma_dossiers/*.json` | No reescanear |
| 8 | Evaluación | `scripts/sigma/50_eval_layers.py`, `51_eval_context.py` | gold | métricas | No usar en producción |
| 9 | GUI | `webapp/browser.py`, `sigma_view.py` | `sigma_promoters.tsv` | visualización | No recalcula |

## Duplicaciones prohibidas

- `12_promoters_sigma.py` y `sigma_annotator.py` **no** deben usarse ambos para
  el mismo promotor. El primero aporta inventario; el segundo es la capa final
  de promotores.
- `12` escribe `promoters.tsv` por compatibilidad histórica; no es la fuente
  de `sigma_promoters.tsv`.
- `webapp/regulation_layer.py` produce `promoters_auto` de inicio de operón;
  no es una llamada de sigma.
- TSS/TTS se consumen de su capa; no se re-predicen.
- Estructura/Foldseek/pockets no deciden promotor; solo identidad/QC.
- El k-mer de *B. subtilis* no se aplica a otras especies.
- El contexto no se usa como prueba; solo como prior con ablación.

## Contratos

### `data/regulation/<g>/regulators_inventory.tsv`

```text
genome, gene_id, gene_name, product, family, role,
pfam_checked, pfam_hit, status, evidence, attempted, confidence
```

### `data/sigma/<g>/sigma_factors.tsv`

```text
genome, gene_id, locus_tag, gene_name, product, sigma_type,
sigma_class, status, confidence, promoter_hits, ortholog_group,
copy_number, evidence, provenance
```

### `data/sigma/<g>/sigma_promoters.tsv`

```text
genome, operon_id, contig, strand, first_gene, sigma,
architecture, spacer, score, zscore, pvalue, margin,
sigma_status, confidence, alternatives, pathway_support, provenance
```

### `reports/sigma_dossiers/<g>.json`

```text
genome, n_factors, n_known, n_lineage_or_partial,
n_unknown_family, n_promoters, factors[]
```

## Dossier actual

Ejecutado sobre 402 genomas:

```text
factores:                 6.001
con ortólogo conocido:    2.141
linaje/parcial:              94
familia desconocida:      3.766
promotores candidatos:   788.933
promotores con ruta:      29.427
genomas sin factores:           0
```

Los 3.766 `unknown_family` son el siguiente cuello: no son necesariamente
sigma nuevos; son inventario sin familia asignable. Hay que resolverlos con
dominios/ortología, no con más escaneo de promotores.

## Estado

- Pipeline base y sigma ya existen para 402 genomas.
- El dossier ya se puede regenerar.
- Lo que falta es mejorar inventario/ortología, contexto y validación,
  no volver a ejecutar la anotación base.
