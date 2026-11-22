# Plan de regulación L2/L3 (actualizado con el estado real)

## Inventario verificado

- `data/test_ecoli.gbff` (11 MB) y RegulonDB (`data/regulondb/gt`): 827 operones,
  4.699 pares génicos de verdad.
- `genes.tsv` de E. coli K-12: 4.318 CDS, 102 tokens sin resolver.
- `operons.tsv`: 821 operones; 59,5% de CDS en operones ≥2; media 3,0 genes.
- `promoters.tsv`: 631/2.599 operones con promotor; 1.968 sin resolver.
- `synteny_pairs.tsv`: 39.407 pares en ≥2 genomas.
- ARN: `rna_annotations.json` **sin hebra** (solo start/end).
- 10 checkpoints borrados del disco pero presentes en git (`D`).

## L2: resultado del gate

| métrica | valor |
|---|---|
| Pares génicos: precisión / recall / F1 | 0,886 / 0,717 / **0,793** |
| Conjunto exacto: precisión / recall / F1 | 0,608 / 0,603 / **0,606** |
| Cobertura de genes | 93,3% |
| Operones 2–3 genes: F1 | 0,637 |
| Operones 4–6 genes: F1 | 0,549 |
| Operones 7+ genes: F1 | **0,395** |
| Operones con promotor | 24% |

**Veredicto:** el gate estricto (≥15/20 operones curados con los mismos genes y
orden) **no está demostrado**: el recall exacto es 0,60 (~12/20). El predictor
es útil a nivel de pares (F1 0,79), pero **no pasa el gate** para escribir L3.

## Estado de higiene (actualizado)

- [x] **ARN con hebra**: recomputado desde NCBI GFF; 20.323 tRNA + 2.998 rRNA,
      100% con `strand`, `product`, `locus_tag`.
- [x] **Parser GFF**: 400/400 `parsed_annotations` con `locus_tag`, `gene`,
      `product`, `EC`; `promoter_reason` documentado.
- [x] **Checkpoints**: los 10 restaurados (`git status` sin `D`).
- [x] **Reporte de evaluación (v1)**: `reports/eval_probe_v1.json`,
      AUROC 0,667 (muestra 65 ventanas; v1 débil en promotores).
- [ ] **Evaluación comprensiva v2**: 13 falla por desajuste v1/v2; falta
      checkpoint v2 o adaptar el script.
- [ ] **Operones 7+ / sigma**: pendiente del `regulon-builder`.

## Qué se corrigió (histórico)

1. **Operones largos (7+)**: recall 0,33; `max_genes=15` y el umbral pueden
   truncar. Añadir evidencia funcional/sintenía y revisar el tope.
2. **Fronteras**: exact-set 0,61; calibrar mejor distancia + RBS + sintenía, y
   evitar sobre-fragmentación.
3. **Promotor/sigma**: 631/2.599. Faltan modelos Pfam de sigma
   (`PF00140`, `PF04542`, `PF04539`, `PF04545`, `PF04997`); el inventario por
   nombre de producto no basta. Descargar/procesar esos HMM y reescanear.
4. **ARN sin hebra**: recomputar `rna_annotations` desde el GFF de NCBI; sin
   hebra no se puede decidir si un tRNA corta un operón.
5. **Parser GFF**: `14_parse_gff.py` descarta locus_tag/product/gene/EC y busca
   promotores donde NCBI no los emite; conservar atributos y no reportar 0 como
   “no hay”.
6. **genemaps**: 466 filas de mediana vs ~2.500 CDS (19%); no es tabla génica.
   La tabla génica completa ahora es `genes.tsv`.
7. **Checkpoints**: restaurar los 10 borrados (`git checkout -- checkpoints_descargados/...`).
8. **Evaluación del modelo**: no hay reportes de 11/12/13; generarlos.

## L3: regulón (solo si L2 pasa)

- Nodos: operones + TFs; aristas TF→operón con signo, score, posición y
  distancia.
- Motivos: tabla curada y pequeña (CRP, FNR, ArcA, NarL, Lrp, ArgR, TrpR, LexA,
  Fur, OxyR, SoxS, Fis, IHF, H-NS, RpoS), con cita; **no inventar motivos**.
- Regla lógica por operón: `lacZYA: ON IF (CRP AND NOT LacI)`.
- **Honestidad**: `logic_status=assumed` salvo evidencia documental; el reporte
  cuenta solo las reglas con respaldo.
- Validación: SOS/LexA, lac, trp, ara, mal.

## L4: FBA (aplazado)

- `cobra 0.32.1` instalado, pero la cobertura gen→reacción es ~19%; un FBA hoy
  dejaría 4/5 del metabolismo sin restringir. No mezclar con L3.

## Integración con el plan general

- L2/L3 es una **vía paralela** a brechas/genes fantasma; comparte `genes.tsv`,
  `functional.tsv`, `rbs_motif` y Pfam.
- El gate de L2 es bloqueante para L3; no se escribe L3 con F1 exacto 0,61.
- La higiene de ARN/parser/checkpoints es P0 porque contamina etiquetas y
  validaciones.

## Próximo

1. L3: regulón TF→operón y validación (el otro agente).
2. Mejorar 7+ y fronteras (F1 0,442 aún bajo).
3. Mantener la higiene base y no tocar `scripts/regulation/`.
