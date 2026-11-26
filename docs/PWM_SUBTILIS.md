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

> **Corrección (2026-11-25).** El 44,79 % que figuraba aquí **no era con la
> matriz propia**: `50_eval_layers.py` parseaba `--pwm` pero lo ignoraba y
> seguía usando `data/regulondb/gt/promoter_pwms.tsv` (E. coli). Con la
> matriz propia, 200 señuelos y ventana −80/+20, el resultado es:
>
> | Configuración | Top-1 | Top-3 |
> |---|---:|---:|
> | E. coli PWM (bug `--pwm`) | 44,79 % | 85,53 % |
> | **PWM propia de B. subtilis** | **29,69 %** | **80,68 %** |
> | PWM propia + capa de contexto | 43,40 % | 86,88 % |
>
> La matriz propia **no supera la mayoría trivial** (87,56 %). El contexto
> recupera parte de la diferencia, pero tampoco la supera. Ver
> [`SIGMA_CONTEXT.md`](SIGMA_CONTEXT.md) para el protocolo, la ablación y
> los límites.

En la evaluación original con el gold de B. subtilis (matriz de E. coli por
el bug de `--pwm`) se reportaba:

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
