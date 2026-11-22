# Auditoría del pipeline `scripts/ghost`

Estado a la fecha. Este documento es la fuente de verdad de **qué hay, qué
hace, qué entra y qué sale** antes de añadir características nuevas.

---

## 1. Inventario de módulos

| módulo | función | entrada | salida |
|---|---|---|---|
| `seqio.py` | FASTA/GFF, traducción, IUPAC | FASTA, GFF3 | secuencias, features |
| `orfs.py` | regiones y ORFs en 6 marcos | secuencia + intervalos | ORFs con ambigüedad |
| `scoring.py` | codones, RBS, composición, ranking | ORF | `codon_llr`, `rbs`, `rank` |
| `pipeline.py` | orquesta candidatos T1 | FASTA + GFF | candidatos + stats |
| `run.py` | CLI de candidatos | flags | TSV/JSON/FAA |
| `protein_lm.py` | ESM-2 8M (CPU) | péptidos | log-verosimilitud |
| `validate.py` | control positivo/negativo | genoma | recall, FPR |
| `gap_conservation.py` | ¿el filo puede producir el metabolito? | traza + redes | `candidate_loss` / `phylum_absent` |
| `gap_demand.py` | rareza global de la brecha | 400 trazas | `specific` / `common` |
| `biomass.py` | ancestros de precursores de biomasa | Rhea RDF | prioridad (exploratorio) |
| `gap_targets.py` | enzimas objetivo y cierres | traza + Swiss-Prot | `*_gap_targets.tsv` |
| `gap_link.py` | candidato → donante → Rhea → brecha | candidatos + trazas | `*_gap_links.tsv` |
| `structure.py` | ESMFold + Foldseek | top de candidatos | PDB, pLDDT, hits |
| `labels.py` | esquema de etiqueta con fuente/evidencia/confianza | cualquier fila | campos `label_*` |
| `empty_regions.py` | desiertos sin CDS/RNA/oriC/móviles | GFF + JSONs | `empty_regions.tsv` |
| `gap_refine.py` | brechas biosintéticas con demanda | trazas + Rhea | `refined_gaps.tsv` |
| `gap_closure.py` | cierre cis candidato↔enzima | candidatos + brechas | `closures.tsv` |
| `tiers.py` | Tier 2: hypothetical/DUF como candidatas | `functional.tsv` + `proteins.faa` | catálogo Tier 2 |
| `tier3.py` | Tier 3: modelos alternativos de genes imprecisos | GFF + FASTA | catálogo Tier 3 |
| `eval_known.py` | DIAMOND + Pfam sobre representantes | `representatives.faa` | `eval_known.tsv` |
| `dedup.py` | dedup global de candidatos (95%) | `*_ghost_candidates.json` | clusters + representantes |
| `dedup_eval.py` | experimento de umbral de dedup | proteínas anotadas | métricas por umbral |
| `dedup_pairs.py` | precisión/recall reales por umbral | all-vs-all | `dedup_pairs.tsv` |
| `report.py` | reporte Markdown | salidas | `*_report.md` |
| `analyze.py` | **orquestador de un genoma** | flags | todo lo anterior + manifiesto |

---

## 2. Puntos de entrada

### Instalado en el proyecto
```bash
cd ~/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc

# un genoma, pipeline completo
bash scripts/ghost/run_one.sh GCF_000006985.1_ASM698v1

# orquestador explícito (candidatos + brechas + reporte)
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.analyze \
  --genome <GCF> --data-root data --out out

# lote (reanudable)
bash scripts/ghost/run_batch.sh        # 400
bash scripts/ghost/run_batch.sh 20     # prueba
```

### Pruebas
```bash
cd scripts/ghost
PYTHONPATH=. ../../.venv/bin/python -m unittest discover -s tests -v
```

---

## 3. Contrato de datos

### Entradas esperadas
```
data/raw/<GCF>_genomic.fna.gz            (o data/raw_excluded/)
data/annotations_ab_initio/<GCF>/annotation.gff3
data/metabolism/<GCF>_network_trace.json
data/metabolism/<GCF>_network.json
data/metabolism/rhea2uniprot_sprot.tsv
data/species_phylum.tsv
data/db/proteomes.dmnd
data/db/swissprot.dmnd
```

### Salidas por genoma
```
out/<GCF>_ghost_candidates.tsv/json
out/<GCF>_ghost_orfs.faa
out/<GCF>_diamond_proteomes.tsv
out/<GCF>_gap_conservation.tsv
out/<GCF>_gap_targets.tsv
out/<GCF>_report.md
out/<GCF>_analysis_manifest.json
out/<GCF>_analyze.log
```

---

## 4. Entorno

- `.venv` del proyecto: `torch 2.14.1+cpu`, `numpy 1.26.4`, `biopython`,
  `pyrodigal`, `sentencepiece`, `fair-esm`.
- Sin GPU. ESM-2 8M corre en CPU (~0,3 s/péptido).
- `TORCH_HOME` debe apuntar a `dna_poc/.cache/torch` (los lanzadores lo hacen).
- `diamond 2.1.11`, `hmmer`, `foldseek` (binario en `dna_poc/foldseek/bin`).

---

## 5. Evidencia de validación

| prueba | resultado |
|---|---|
| Unitarias + sintéticas + dedup + tiers + labels (`tests/`) | **41/41 OK** |
| End-to-end sintético (ADN ambiguo) | ORF con X/variantes, CDS conocido filtrado |
| Recuperación de 100 CDS ocultados, genoma 1 | recall **0,99**, identidad 0,99 |
| Ídem, genoma 2 | recall **0,99**, identidad 0,99 |
| Control negativo (secuencias barajadas) | FPR **0,000** en umbrales 0,75 y 0,55 |
| ESM-2 8M reales vs barajadas | AUC **0,987** |
| Foldseek (beta-lactamasa) | AF-P62593 al **99,6%**, e-value 2·10⁻⁵³ |
| Orquestador `analyze` (genoma 1) | ver manifiesto |

Calibración: CDS reales 0,62 de LLR; ORFs intergénicos −0,24.

---

## 6. Problemas conocidos y deuda técnica

1. **`partial` ORFs** (sin codón de parada, corren hasta el borde). Medido:
   excluirlos reduce candidatos 637 -> 401 (-37%) pero baja el recall de
   0,99 a 0,97 en el genoma 1 (0,99 en el genoma 2). Decisión: se marcan con
   `partial` y se **penalizan x0,7** en el ranking; `--exclude-partial` los
   elimina si se prioriza precisión.
2. **`biomass.py` demasiado permisivo**: 13.074 ancestros (43% de Rhea). Se
   conserva como exploratorio; el filtro que se usa es `gap_conservation`.
3. **`gap_targets.py` desincronizado** entre workspace y proyecto (se corrige
   en esta auditoría).
4. **ESMFold API**: 504 para 65–120 aa; `structure.py` reintenta y degrada.
5. **Salidas planas** en `out/`: no hay subdirectorio por genoma.
6. **`run_batch.sh`** no genera reporte global al final.
7. **`enumerate_regions`** en `orfs.py` duplica lógica de `pipeline.py` y no
   se usa; candidato a eliminar.
8. **Sin Pfam, T2 ni T3**: el dedup global ya está implementado y el umbral
   elegido es **95%** (`out/dedup/REPORT.md`). Falta correrlo sobre el catálogo
   de los 400 y añadir dominios/T2/T3.
9. **`report.py`** depende de `species_phylum.tsv` para inferir donantes.
10. **GUI**: cubre CDS, tRNA, rRNA, ORFs fantasma y ESM. Faltan ORCA,
    transposones, defensa y BGCs (dependen de bases que no se distribuyen).

---

## 7. Checklist de reproducibilidad

- [x] Entorno documentado y aislado (`.venv` del proyecto).
- [x] Bases de datos presentes (`proteomes.dmnd`, `swissprot.dmnd`,
      `foldseek_swissprot`).
- [x] Pesos ESM en `.cache/torch`.
- [x] Semilla fija en validación (`--seed 42`).
- [x] Pruebas automáticas (20).
- [x] Manifiesto por corrida (`*_analysis_manifest.json`).
- [ ] Versionado de los umbrales (`min_aa`, `pad`, `min_region_len`) en config.
- [ ] Snapshot de las bases (hash) por corrida.

---

## 8. Qué sigue (después de esta auditoría)

1. **Dedup** con `diamond cluster` y elección empírica del umbral contra genes
   conocidos (90/95/97%).
2. **Pfam-A** + `hmmscan` para dominios.
3. **T2** (53.596 hypothetical) y **T3** (anotación imprecisa).
4. **`gap_closure`**: empalme candidato ↔ enzima objetivo con evidencia.
5. **Lote 400** solo cuando 1–4 estén validados.
