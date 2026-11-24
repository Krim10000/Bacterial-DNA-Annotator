# Tres caminos de anotación de proteínas

## Estado de implementación (2026-11-23)

Ahora el pipeline de un FASTA nuevo:

- ejecuta DIAMOND/Swiss-Prot automáticamente;
- escribe `*.functional.tsv`, `*.functional_annotation.json` y
  `*.protein_paths.tsv`;
- clasifica cada proteína en A, B o C;
- conserva hit, identidad, e-value y producto.

Para los 400, `functional.tsv` ya existía, pero `protein_paths.tsv` no estaba
persistido. Se reconstruyó con `tools/rebuild_protein_paths.py` sin volver a
correr DIAMOND: **400/400 genomas** tienen ahora `protein_paths.tsv`.

Todavía no ejecuta automáticamente los pasos estructurales de B/C. Para B/C el
pipeline deja el trabajo preparado (`next_step`) y los filtros estructurales
siguen siendo opcionales (`--structure-manifest`, `--structural-filters`).
El camino C sigue siendo un plan CPU para los 400; no se debe presentar como
completo.

## Etiquetas interpretables

`protein_paths.tsv` ahora incluye:

| columna | valores |
|---|---|
| `path` | A / B / C |
| `label` | `A_named_function`, `B_uncharacterized_homolog`, `C_no_swissprot_hit` |
| `evidence_strength` | `strong`, `moderate`, `weak`, `none` |
| `evidence_toward` | `named_function`, `uncharacterized_family`, `unknown` |
| `interpretation` | explicación legible con pident y cobertura |
| `subunit_hint` | `named_subunit` / `none`; solo si el producto contiene *subunit* o *component* |
| `subunit_evidence` | `product_name_only`; no es prueba de pertenencia a complejo |
| `structural_status` | `afdb_model_available` / `no_afdb_model` (solo B) |
| `domain` | mejor dominio Pfam del panel para C; hipótesis |
| `domain_evalue` | e-value del dominio; no confirma función |
| `tm_segments` | segmentos hidrofóbicos heurísticos (C); no TMHMM |
| `signal_peptide` | `yes/no/uncertain` (C); no SignalP |
| `lipobox` | motivo [LVI][ASTM][AG]C (C) |
| `context_hypothesis` | `membrane_like/secreted_like/soluble_like/unknown`; hipótesis |
| `tmhmm_status` | `ok/not_installed/error/no_sequence/no_prediction`; opcional |
| `filter_status` | B: `evaluated_local_no_external` / `not_evaluated_no_local_db` |
| `filter_verdict` | B: `local_only_partial` cuando no se consultan F3/F4 |

Reglas:

- **A** puede ser `strong` si pident ≥70 % y cobertura ≥80 %; aun así es
  homología, no función confirmada.
- **B** nunca pasa de `moderate`; una DUF no resuelve la función.
- **C** es `none`: no prueba ausencia de función.

## Resumen

| camino | proteínas | % total | dato clave |
|---|---:|---:|---|
| A. Función nombrada | **566.957** | 42,4% | Swiss-Prot con función |
| B. UniProt uncharacterized/DUF | **123.200** | 9,2% | homólogo no caracterizado/DUF |
| C. Sin hit Swiss-Prot | **647.588** | 48,4% | sin función curada; 446.927 clusters |

Total: **1.337.745** proteínas. A + B = 690.157 (el antiguo conteo “A” incluía
las B). B ya no se presenta como subconjunto de A: es una etiqueta separada.

`subunit_hint=named_subunit`: 80.839 proteínas (6,0 %). Es solo nombre, no
pertenencia a complejo.

## Camino A — función nombrada

- 566.957 proteínas.
- Hit Swiss-Prot con función.
- Se usa la anotación existente.
- No requiere estructura para una función básica.

## Camino B — uncharacterized/DUF con UniProt

- 123.200 proteínas.
- 25.058 clusters al 50%.
- 3.931 accesiones únicas.
- 1.192 con AFDB.
- 2.739 sin AFDB (2.570 ≤600 aa; 169 >600).
- Se aplican los 5 filtros y se propagan hipótesis.

Salidas:

- `UNCHAR_*`;
- `structural_hypotheses`;
- `CANDIDATES_508_REVIEW`.

## Camino C — sin hit Swiss-Prot

- 647.588 proteínas.
- 446.927 clusters al 50%.
- 360.741 clusters singleton (80,7% de clusters; 55,7% de proteínas).
- No hay AFDB directo para ORFs arbitrarios.
- No se pliegan todos.

### Singleton

Un cluster singleton tiene **una sola secuencia** en el conjunto a ese umbral.

No significa:

- que no tenga homólogo en UniRef/PDB;
- que sea una proteína rara;
- que no tenga dominio.

Sí significa:

- no se reduce por clustering dentro de este conjunto;
- hay que anotarlo por dominio/homología/contexto;
- puede ser un ORF dudoso y requerir control de calidad.

### Qué se hace en C

1. Reclustering/representantes.
2. MMseqs2/DIAMOND vs UniRef90/TrEMBL/PDB/eggNOG.
3. Pfam/InterPro en los que sigan sin señal.
4. TMHMM/Phobius/señal.
5. Contexto génico/operón.
6. Taxonomía/conservación.
7. Estructura solo para el top.

## Caminos en el pipeline

La documentación del pipeline debe distinguir:

```text
A: anotación existente
B: hipótesis estructurales sobre UniProt/AFDB
C: reducción CPU de no-hit + estructura selectiva
```

El camino B está implementado con `--structure-manifest`,
`--structural-filters` y `--precomputed-hypotheses`.

El camino C todavía es un plan CPU; no está automatizado por completo.

La etiqueta A/B/C **no es una llamada de complejo**. B estructural solo
registra si hay modelo AFDB local y, opcionalmente, filtros; C dominio/contexto
solo aporta familia Pfam y localización putativa. `subunit_hint` es un eje
independiente y no cambia la etiqueta; ver [`SUBUNIT_NAMES.md`](SUBUNIT_NAMES.md).

## Complejos

Una subunidad de un complejo obligado puede caer en A/B/C; el pipeline no
ensambla el complejo ni exige estequiometría. Ver [`COMPLEX_GAPS.md`](COMPLEX_GAPS.md).
No confundir "subunidad anotada" con "complejo funcional".

`subunit_hint` es solo el nombre del producto. En K-12, 20,5 % de los nombres
con *subunit/component* quedan fuera del catálogo curado, pero varios son
subunidades reales no catalogadas; 44,2 % de las subunidades curadas no llevan
esa palabra. Ver [`SUBUNIT_NAMES.md`](SUBUNIT_NAMES.md).

## Regla

- Un singleton no es “basura”.
- Un no-hit no es “sin homología”.
- Una hipótesis no es `confirmed`.
- La anotación completa significa evidencia completa, no función inventada.
