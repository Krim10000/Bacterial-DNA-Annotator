# Calibrador k-mer de sigma

## Qué es

Un modelo **específico de especie** que estima, para cada sigma del inventario,
una probabilidad a partir de n-gramos de la ventana del promotor.

No sustituye a la PWM: la política `hybrid` solo usa el k-mer cuando su
probabilidad supera un umbral; si no, cae a `two_stage`.

## Entrenamiento

```bash
python scripts/sigma/63_train_kmer_calibrator.py \
  --genome subtilis \
  --gold data/sigma/gold/subtilis/tss_sigma.tsv \
  --source subtilis \
  --up 40 --down 5 --ngram 4 \
  --out data/sigma/policy/subtilis_kmer.json
```

El JSON contiene vocabulario, escala, coeficientes e interceptos. No se
guardan secuencias ni etiquetas individuales.

## Uso

```bash
python scripts/sigma/10_annotate_sigma.py \
  --genome subtilis \
  --pwm data/dbtbs/promoter_pwms_subtilis_all.tsv \
  --kmer-model data/sigma/policy/subtilis_kmer.json \
  --policy hybrid --kmer-threshold 0.6
```

El módulo sigma usa automáticamente el modelo incluido para `subtilis` y la
política `hybrid`.

## Resultado

GroupKFold 5, ventana −40/+5, 4-mer:

```text
top-1:     0,9003
macro-F1:  0,5184
mayoría:   0,8756
```

## Límites

- Se entrena con el gold de esa especie; no es universal.
- No es validación independiente.
- El k-mer puede aprender sesgos de anotación.
- No hay negativos reales.
- Si falta el modelo, no se inventan probabilidades.
