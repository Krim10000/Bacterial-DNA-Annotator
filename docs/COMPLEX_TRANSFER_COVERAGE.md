# Cobertura de transferencia K-12 sin catálogo ComplexPortal independiente

Pregunta: en los 400 genomas del panel, ¿qué cubre la transferencia de
complejos K-12 cuando **no** hay un catálogo ComplexPortal independiente contra
el cual compararla?

La respuesta corta: cubre **19.085 complejos transferidos** y **56.222 de
95.278 subunidades** (59.01 %), pero **no es una validación**. Ningún genoma
del panel tiene catálogo ComplexPortal independiente local: los 400 quedan en
`not_available`. El único catálogo local es el de K-12, que es el prior de la
transferencia y no se usa como validación independiente.

## Qué se auditó

Script: [`scripts/complexes/transfer_coverage_audit.py`](../scripts/complexes/transfer_coverage_audit.py).

Insumos, solo lectura:

- `data/complexes/transfer_<genoma>.tsv` — transferencia por ortología K-12;
- `data/complexes/model_<genoma>.complexes.tsv` — candidatos del modelo K-12;
- `data/complexes/gold_k12.tsv` — catálogo K-12, solo para mapear b-numbers a
  nombres de gen y comparar con las firmas curadas;
- `data/complexes/catalog_<genoma|taxid>.tsv` o
  `data/complexes/gold_<genoma|taxid>.tsv` — catálogo independiente, si existe;
- `data/taxonomy/full/<genoma>.json` — opcional, solo para reconocer un
  catálogo independiente nombrado por taxid.

No se descargó ComplexPortal, no se consultó la red y no se modificó ningún
`transfer_*.tsv`, `model_*.tsv` ni catálogo.

## Resultado en los 400

Generado en [`reports/complex_transfer_coverage.tsv`](../reports/complex_transfer_coverage.tsv)
y [`reports/complex_transfer_coverage.md`](../reports/complex_transfer_coverage.md).

| métrica | valor |
|---|---:|
| genomas auditados | **400** |
| con `transfer_*.tsv` con datos | **399** |
| con `transfer_*.tsv` solo cabecera | **1** (`GCF_021057185.1_ASM2105718v1`) |
| con `model_*.complexes.tsv` | **400** |
| con catálogo ComplexPortal independiente | **0** |
| sin catálogo independiente (`not_available`) | **400** |
| complejos transferidos | **19.085** |
| transferidos con `frac_mapped >= 0.8` (`transfer_suspected`) | **8.408** (44.06 %) |
| transferidos con `frac_mapped < 0.8` (`transfer_candidate`) | **10.677** |
| subunidades esperadas (suma `n_subunits`) | **95.278** |
| subunidades mapeadas (suma `n_mapped`) | **56.222** |
| cobertura global de subunidades | **0.5901** (59.01 %) |
| filas con `n_mapped == n_subunits` | **7.398** |
| estequiometría `complete` / `incomplete` / `unknown` | **0 / 1.785 / 17.300** |
| coincidencias con firmas curadas | **1.180** |
| firmas `complete` / `incomplete` | **247 / 933** |
| genomas con al menos una firma | **382** |
| firmas distintas con match | **13 de 18** |
| candidatos del modelo K-12 | **281.090** |
| complejos K-12 distintos transferidos (IDs únicos en los 400 TSV) | **258 de 273** |

La media no ponderada de las medias de `frac_mapped` por genoma es **0.7041**;
la cobertura global (subunidades mapeadas / esperadas) es **0.5901**. Son
métricas distintas: la primera da el mismo peso a cada genoma, la segunda pesa
por subunidades.

## Estados

`transfer_suspected` usa el mismo umbral que el pipeline
(`complex_suspected`, `frac_mapped >= 0.8`). `transfer_candidate` es la
transferencia con menos subunidades mapeadas. Ninguno de los dos es
`confirmed`.

La estequiometría se declara solo cuando `expected_n` está presente:

- `complete`: `observed_n >= expected_n`; **0 filas** en el panel actual;
- `incomplete`: `expected_n` presente y `observed_n < expected_n`; **1.785**;
- `unknown`: sin `expected_n`; **17.300**.

`n_mapped == n_subunits` (7.398 filas) es otra cosa: todas las subunidades del
catálogo K-12 están mapeadas, aunque `expected_n` no permita evaluar
estequiometría. El catálogo K-12 tampoco es un patrón de estequiometría para
otra especie.

## Firmas curadas

Se reutilizan las firmas de [`webapp/complexes.py`](../webapp/complexes.py)
(`SIGNATURES`) y se cruzan con los nombres de gen K-12 de `gold_k12.tsv`:

- `signature_complete`: están todos los genes esenciales de la firma;
- `signature_incomplete`: hay al menos 2 genes de la firma, pero falta algún
  esencial.

Una firma completa **no** confirma el complejo: solo dice que la transferencia
contiene los genes esenciales de esa firma. La lista es mínima y no es un
catálogo exhaustivo.

## Modelo K-12

Los 400 genomas tienen `model_*.complexes.tsv`. Sus **281.090** filas son
`complex_candidate`: sospechas por vecindad/operón, sin estructura ni
estequiometría. El modelo no es un catálogo y no se suma a la transferencia
como si fueran dos validaciones independientes.

## Limitaciones

- **No hay catálogo independiente.** `not_available` en las 400 filas significa
  ausencia local, no ausencia remota en ComplexPortal. La comprobación remota
  está en [`COMPLEXPORTAL_400.md`](COMPLEXPORTAL_400.md); este script no
  consulta la red.
- **La transferencia no es validación.** Es RBH K-12 con identidad >= 30 % y
  evalue <= 1e-5, sin cobertura de alineamiento en los TSV actuales.
- **Cobertura baja no es ausencia.** `frac_mapped` bajo puede deberse a
  complejo ausente, anotación incompleta, fragmentación de complejos grandes u
  ortología no detectada.
- **Estequiometría incompleta.** `expected_n` falta en la mayoría de las filas;
  no hay ninguna fila `complete` en el panel.
- **Firmas mínimas.** Solo usan nombres de gen K-12 y genes esenciales; no
  prueban interacción física ni estequiometría.
- **Modelo candidato.** Cada fila del modelo K-12 es una hipótesis; no cierra
  brechas ni demuestra interacción física.
- **Un genoma sin transferencia.** `GCF_021057185.1_ASM2105718v1` tiene el TSV
  de transferencia solo con cabecera; su fila queda con 0 complejos
  transferidos y sí tiene modelo K-12.

## Reproducir

```bash
.venv/bin/python scripts/complexes/transfer_coverage_audit.py
.venv/bin/python -m pytest tests/test_transfer_coverage_audit.py -q
```

Salidas:

- [`reports/complex_transfer_coverage.tsv`](../reports/complex_transfer_coverage.tsv):
  una fila por genoma, con catálogo, transferencia, estados, subunidades,
  firmas, modelo y limitaciones;
- [`reports/complex_transfer_coverage.md`](../reports/complex_transfer_coverage.md):
  resumen, estados, cobertura, rankings y tabla completa.

Tests: `tests/test_transfer_coverage_audit.py`, con TSV sintéticos. No leen los
400 genomas ni descargan ComplexPortal. Comprueban catálogo independiente
`not_available`, transferencia como inferencia, conteos, firmas, TSV/Markdown,
taxid y ausencia de modificación de los insumos.
