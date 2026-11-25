# Paridad con genomas de referencia

## Objetivo

Medir, sin evaluar función ni coordenadas, si el número de features que anota el
sistema se acerca al de la anotación curada de los genomas del panel
`data/reference_genomes.json`.

La comparación es deliberadamente de **conteos por clase**:

- CDS;
- tRNA;
- rRNA.

No mide productos, función, exactitud de coordenadas, identidad de secuencia ni
calidad de las predicciones. Un conteo parecido no implica que las features sean
las mismas.

## Entradas y selección

- Manifiesto: `data/reference_genomes.json`.
- Resumen del sistema: campo `system_summary` de cada registro, normalmente
  `results/<id>/<id>.summary.json`.
- Referencia: campo `annotation_path`, normalmente
  `data/reference_genomes/<id>/<id>_genomic.gff.gz`.

Se incluye un genoma cuando cumple las dos condiciones:

1. `status` es `annotated` o `available`;
2. declara `system_summary`.

Los registros con status elegible pero sin `system_summary` no se inventan ni se
completan con otros archivos: quedan listados como omitidos en el informe. Los
registros con otros status (`downloaded`, etc.) también quedan fuera.

## Reglas de conteo

- Se leen GFF3 planos o `.gz`; se ignoran comentarios y la sección `##FASTA`.
- Solo se usa la columna de tipo: `CDS`, `tRNA` y `rRNA`.
- Un CDS multiparte (varias filas con el mismo `ID` en el mismo contig) cuenta
  como **un** feature; las filas sin `ID` cuentan una a una.
- No se filtran pseudogenes por atributo (`pseudo=true`). Los registros
  `pseudogene` no son `CDS` y por tanto no entran en el conteo de CDS; si un CDS
  lleva `pseudo=true`, se cuenta como CDS porque la comparación es por tipo.
- No se cuentan `ncRNA`, `tmRNA`, `SRP_RNA`, `RNase_P_RNA`, `gene`, `exon` ni
  otros tipos, aunque estén en el GFF de referencia.
- Los conteos del sistema se toman del bloque `features` del resumen; no se
  recalculan desde el GFF3 del sistema.

## Cálculo y estados

Para cada genoma y clase:

- `delta = sistema - referencia`;
- `% sistema/referencia = 100 * sistema / referencia` (vacío si la referencia
  no tiene esa clase).

Estados:

| estado | significado |
|---|---|
| `ok` | CDS, tRNA y rRNA coinciden exactamente y la referencia tiene las tres clases. |
| `partial` | La comparación fue posible, pero al menos un delta es distinto de cero, o la referencia carece de alguna clase. |
| `not_available` | Falta o es ilegible el resumen o el GFF, el manifiesto no declara la ruta, o el GFF no contiene ninguna de las tres clases. |

El estado `partial` es el esperado cuando hay diferencias reales de conteo; no
se usa un umbral arbitrario de tolerancia. El agregado "micro" del informe suma
sistema y referencia por clase sobre los genomas comparables.

## Salidas y uso

```bash
.venv/bin/python scripts/reference/reference_parity.py
```

Genera:

- `reports/reference_parity.tsv`: una fila por genoma, con conteos, delta,
  porcentaje, estado, rutas y motivo. Es la salida para consumo automático.
- `reports/reference_parity.md`: informe legible con cobertura, tabla por
  genoma, motivos, entradas, agregado y límites.

Opciones útiles:

```bash
.venv/bin/python scripts/reference/reference_parity.py \
  --root . \
  --manifest data/reference_genomes.json \
  --tsv reports/reference_parity.tsv \
  --md reports/reference_parity.md
```

El script no usa red, GPU ni binarios externos. No modifica `results/`, los GFF
ni el manifiesto.

## Resultados actuales

Instantánea del **2026-11-25** (para cifras vigentes, ver
`reports/reference_parity.md`). En ese momento había 42 registros en el
manifiesto: 7 evaluados, 2 omitidos con status elegible pero sin
`system_summary` (`ecoli_k12`, `subtilis`) y 33 con otros status.

| genoma | estado | CDS sis | CDS ref | Δ CDS | % CDS | tRNA sis | tRNA ref | Δ tRNA | % tRNA | rRNA sis | rRNA ref | Δ rRNA | % rRNA |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mtb_h37rv` | partial | 3986 | 3906 | +80 | 102.05 | 45 | 45 | 0 | 100.00 | 2 | 3 | -1 | 66.67 |
| `cglutamicum` | partial | 3041 | 2999 | +42 | 101.40 | 60 | 59 | +1 | 101.69 | 12 | 18 | -6 | 66.67 |
| `pputida` | partial | 5507 | 5471 | +36 | 100.66 | 75 | 75 | 0 | 100.00 | 15 | 22 | -7 | 68.18 |
| `paeruginosa` | partial | 5715 | 5572 | +143 | 102.57 | 64 | 63 | +1 | 101.59 | 8 | 13 | -5 | 61.54 |
| `crescentus` | partial | 3738 | 3886 | -148 | 96.19 | 52 | 51 | +1 | 101.96 | 4 | 6 | -2 | 66.67 |
| `synechocystis` | partial | 3780 | 3698 | +82 | 102.22 | 41 | 41 | 0 | 100.00 | 2 | 6 | -4 | 33.33 |
| `saureus` | partial | 2675 | 2767 | -92 | 96.68 | 61 | 59 | +2 | 103.39 | 11 | 16 | -5 | 68.75 |

Agregado micro de esos 7 genomas:

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 28442 | 28299 | +143 | 100.51 |
| tRNA | 398 | 393 | +5 | 101.27 |
| rRNA | 54 | 84 | -30 | 64.29 |

Lectura rápida: los conteos de CDS y tRNA están cerca de la referencia en
términos agregados; rRNA queda sistemáticamente por debajo. Este informe no
atribuye esa brecha a una causa concreta, pero es compatible con que el pipeline
no llame todas las copias o tipos de rRNA (por ejemplo 23S/5S). Confirmarlo
requiere analizar las coordenadas y el caller de rRNA, fuera del alcance de esta
versión.

## Límites

- **Solo conteos**: no evalúa función, productos, coordenadas, ni identidad de
  secuencia.
- **El gold no es verdad absoluta**: hereda criterios y posibles errores de la
  anotación RefSeq/NCBI usada como referencia.
- **rRNA**: las diferencias pueden deberse a tipos no detectados (23S/5S),
  copias parciales o diferencias de criterio; no se diagnostican aquí.
- **ncRNA y otros tipos** quedan fuera, aunque estén en el GFF.
- **Cobertura del panel**: solo se comparan genomas ya anotados y con
  `system_summary`; los demás no se estiman ni se imputan.
- **Sin causalidad**: un delta no indica por sí solo un error del sistema ni de
  la referencia.
- **Sin uso para entrenamiento**: este panel es de evaluación; no se usa como
  gold de entrenamiento ni para declarar superioridad.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_reference_parity.py -q
```

Los tests usan manifiesto, resumen y GFF sintéticos; cubren paridad exacta,
deltas, porcentajes, CDS multiparte, gzip, ausencia de resumen/GFF, resumen
ilegible, GFF sin clases objetivo y filtrado del manifiesto.
