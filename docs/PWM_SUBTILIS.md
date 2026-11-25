# PWMs específicas de *B. subtilis*

## Qué se construyó

`scripts/validation/build_sigma_pwms_subtilis.py`

Salida:

```text
data/dbtbs/promoter_pwms_subtilis_all.tsv
```

Se deriva de:

```text
data/dbtbs/gt_tss.tsv
data/dbtbs/genome.fna
```

No se transfieren matrices de *E. coli*.

## Matrices

| Etiqueta | Sigma canónico | TSS | Caja | Consenso | Info |
|---|---|---:|---|---|---:|
| A | sigma70 | 1.949 | -10 | TATAATTAAAAT | 1,71 |
| B | sigma38 | 182 | -10 | ATAAAAAAAAAA | 1,42 |
| D | sigma28 | 71 | -10 | TAAAAATAAATA | 2,11 |
| I | sigma32 | 13 | -10 | TCAAAAACAAAA | 2,68 |
| L | sigma54 | 11 | -10 | CTGGAATATTAA | 4,01 |

Para σA no aparece caja −35 estable; los promotores son mayoritariamente
*extended −10*. Las matrices B/D/I/L se construyen con pocos TSS y son
hipótesis de baja cobertura, no equivalentes funcionales de los sigma de K-12.

## Uso

En el módulo sigma, para `subtilis` se usa automáticamente:

```text
data/dbtbs/promoter_pwms_subtilis_all.tsv
```

CLI:

```bash
python scripts/sigma/10_annotate_sigma.py \
  --genome subtilis \
  --pwm data/dbtbs/promoter_pwms_subtilis_all.tsv
```

Evaluación:

```bash
python scripts/sigma/50_eval_layers.py \
  --genome subtilis --source gold \
  --pwm data/dbtbs/promoter_pwms_subtilis_all.tsv
```

## Resultado

En la evaluación con el gold de B. subtilis:

```text
top-1: 44,79 %
top-3: 85,53 %
mayoría trivial: 87,56 %
```

La matriz propia no supera la mayoría trivial. El problema no era solo la
matriz: la arquitectura de promotores de B. subtilis y la mezcla de TSS
primarios/secundarios siguen limitando el top-1.

## Límites

- No es una validación experimental.
- B/D/I/L tienen muy pocos TSS.
- No se inventa caja −35.
- σ54 se escanea con una sola caja; su geometría ideal es −24/−12.
- El gold es DBTBS/BSGatlas.
