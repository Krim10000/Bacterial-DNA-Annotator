# Pipeline headless: anotar un FASTA desde cero

```bash
python -m webapp.pipeline genoma.fna --outdir results/ --threads 8
python -m webapp.pipeline a.fna b.fna.gz --outdir results/ --esm
python -m webapp.pipeline genoma.fna --outdir results/ \
  --esmfold "neuros nap" --esmfold-gaps 3 \
  --structure-evidence --structural-annotation --structural-filters \
  --structure-manifest manifest.tsv \
  --foldseek-local-db data/db/foldseek_swissprot
```

## Salidas por genoma

- `<genoma>.gff3`: CDS, tRNA, rRNA.
- `<genoma>.proteins.faa`: proteínas traducidas.
- `<genoma>.candidates.tsv`: ORFs de novo.
- `<genoma>.gaps.tsv`: brechas con rol, origen, transporte, dirección y HMM.
- `<genoma>.report.md`: informe con cobertura, desiertos, candidatos y brechas.
- `<genoma>.summary.json`: resumen legible por máquina.
- `<genoma>.esmfold_closure.tsv`: cierre por estructura (si `--esmfold`).
- `<genoma>.structural_evidence.tsv`: evidencia multi-fuente (si `--structure-evidence`).
- `<genoma>.enzymm_scan.tsv`: plantillas M-CSA (si `--structural-annotation`).
- `<genoma>.structural_filters.tsv`: 5 filtros (si `--structural-filters`).
- `<genoma>.structural_hypotheses.tsv/.json`: hipótesis revisables (si `--structure-manifest` o filtros).

## Opciones

- `--no-cds`, `--no-trna`, `--no-rrna`, `--no-ghost`, `--no-deserts`.
- `--esm`: puntuar ORFs con ESM-2.
- `--no-metabolic`: omitir red y brechas.
- `--reference GCF`: reutilizar red de referencia.
- `--panel-hmm`, `--transport-db`: rutas a recursos.
- `--esmfold off|local|"neuros nap"`: cierre de brechas con ESMFold.
- `--esmfold-gaps N`: máximo de brechas a plegar (default 3).
- `--esmfold-out DIR`: PDB/métricas.
- `--esmfold local` exige GPU/CUDA/esm/openfold; si no, el pipeline avisa
  y no simula el resultado.
- `--structure-evidence`: Foldseek PDB100 + AlphaFill + sitios activos.
- `--structural-annotation`: EnzyMM/PyJess contra plantillas M-CSA.
- `--structural-filters`: F1 Foldseek, F2 EnzyMM, F3 P2Rank, F4 AlphaFill, F5 sitio activo.
- `--structure-manifest`: TSV `genome protein_id uniprot pdb` para aplicar filtros a proteínas con estructura.
- `--foldseek-local-db`: base Foldseek local; evita la API y el 429.
- `--precomputed-hypotheses`: TSV de hipótesis precalculadas; filtra por genoma.

## Garantías

- Por defecto no sube datos a la nube.
- `--esmfold "neuros nap"` **sí** envía la secuencia a Neurosnap; requiere autorización y créditos.
- `--esmfold local` no envía secuencias, pero exige GPU y dependencias.
- Brechas con `direction_unspecified` si no hay ΔG.
- Transporte solo T1 bacteriano; lo demás queda en revisión.
- Informe declara fuentes y límites.

## Explicaciones por familia

El informe incluye una tabla por familia KEGG con:

- brechas,
- homólogo presente,
- transporte,
- explicación,
- acción.

Se calcula con `webapp/family_logic.py` y aplica a cualquier FASTA.

## Sugerencias de mapeo (subidos)

Para cada brecha con homólogo de enzima productora, el pipeline adjunta:

- `homolog_rhea`: reacción propuesta.
- `homolog_ec`: EC de la reacción si existe.
- `mapping_suggestion`: `suggested` (fuerte) o `review` (débil).

El informe incluye una tabla **“Sugerencias de mapeo (revisables)”**. No se
aplican automáticamente; el usuario decide.

## Clases de revisión

Cada brecha de un FASTA nuevo recibe `review_class`:

- `transport`: hay reacción de transporte (T1/T2).
- `annotation_gap`: homólogo fuerte/débil de la enzima productora.
- `unknown`: falta evidencia.

Aparece en `gaps.tsv`, en la GUI y en el informe.
