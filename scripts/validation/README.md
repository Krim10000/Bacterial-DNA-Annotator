# Validación de la capa regulatoria (L2/L3)

Verdad de referencia y *scorer* para medir el predictor de operones (L2) y el
regulón (L3) contra datos curados, en vez de contra una lista hecha a mano.

## Por qué

RegulonDB v13+ dejó de publicar ficheros planos: los datasets se sirven por
GraphQL (`getDataOfFile`). `regulondb_gt.py` los descarga una vez, los
normaliza y deja un manifiesto con versión, fecha y sha256.

El punto clave: los genes vienen en b-numbers (`b0001`…), que son exactamente
los `locus_tag` de NC_000913.3 — el mismo genoma que `data/test_ecoli.gbff`.
No hay que mapear identificadores.

## Regenerar

```bash
python scripts/validation/regulondb_gt.py            # usa caché si existe
python scripts/validation/regulondb_gt.py --fetch    # fuerza descarga
python scripts/validation/regulondb_gt.py --offline  # prohíbe red
```

Salida en `data/regulondb/` (ignorado por git: uso académico/no comercial):

| fichero | filas | contenido |
|---|---|---|
| `gt_operons.tsv` | 2.605 | operón, coords, hebra, genes, confianza, evidencia |
| `gt_tus.tsv` | 3.767 | unidad de transcripción + promotor + sigma + TFs |
| `gt_promoters.tsv` | 4.057 | TSS, sigma, cajas −10/−35 con posición y secuencia |
| `gt_tfs.tsv` | 243 | TF, gen, b-number, efector (metabolito) |
| `gt_gene_aliases.tsv` | — | nombre de gen → b-number |
| `manifest.json` | — | versión 14.5.1, fecha, recuentos, sha256, cita |

Solo 827 operones son puntuables con exactitud (≥2 genes y genes mapeados).
Los 24 restantes son elementos de inserción cuyos nombres no resuelven
(`insAB4` vs `insA4|insB4`); se excluyen con `mapping_ok=no`.

## Puntuar un predictor

```bash
python scripts/validation/benchmark_operons.py \
    --pred data/regulation/ecoli_k12/operons.tsv --json out.json
```

El TSV predicho necesita cabecera y, como mínimo: `genes` (separador coma,
barra o punto y coma; b-numbers o nombres) y `strand` (`+`/`-`). Se aceptan
alias de columna (`genes|gene_ids|gene_list`, `strand|sign`).

Métricas:

- **pair-level** (principal): de los pares de genes que la predicción declara
  co-transcritos, cuántos lo están de verdad. Robusta a diferencias de borde.
- **exact-set**: coincidencia exacta del conjunto de genes. Dura; sirve para
  distinguir "casi bien" de "bien".
- Cobertura de genes, genes no resueltos y genes en más de un operón.

Control de sanidad: puntuar `gt_operons.tsv` contra sí mismo debe dar
P = R = F1 = 1.000.

## Límites que hay que declarar

1. Es **un** genoma (K-12 MG1655) de los 400 del proyecto. La calibración del
   umbral de distancia debe ser por genoma; K-12 no es representativo.
2. RegulonDB está curado sobre literatura de *E. coli*: sobrerrepresenta genes
   estudiados. Un F1 alto en K-12 no garantiza transferencia a otros filos.
3. Los operones de un solo gen no aportan a la métrica de pares y no deben
   usarse para inflar el recuento.
4. Uso académico/no comercial. Los datos no se versionan; se re-descarga con el
   manifiesto como referencia de versión.

## No editar

`scripts/validation/` y `data/regulondb/` son el contrato de validación. Si
hace falta otra vista del ground truth, se pide; no se edita desde el
predictor.

## PWMs desde sitios reales

```bash
python scripts/validation/build_tf_pwms.py --test --promoters
```

Genera, a partir de RegulonDB:

| salida | contenido |
|---|---|
| `data/regulondb/gt/tf_pwms.tsv` | 82 TF con ≥8 sitios (3.790 sitios), log-odds + consenso |
| `data/regulondb/gt/promoter_pwms.tsv` | cajas −10/−35 por factor sigma |

Control de cordura: los consensos de σ70 salen `TATAAT` (−10) y `TTGACA` (−35),
y σ32/σ54 reproducen sus elementos canónicos. AUC self-test mediano de las PWM
de TF: 0.985.

Para qué: la capa L2 usa hoy un patrón con umbral de mismatches y su precisión
posicional a ±10 nt es P=0.31 / R=0.19. Una PWM elige el máximo local, que es
lo que hace falta para decidir si un operador solapa el −10 (de ahí sale el
signo de una arista de L3).
