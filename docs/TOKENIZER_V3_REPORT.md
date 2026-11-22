# Tokenizer v3 — entrenamiento y evaluación

## Objetivo

Corregir los problemas del tokenizer v2 sin perder resolución biológica:

- vocabulario más fino que 4,2 nt/token;
- IUPAC completo (v2 manda `B,D,H,V` a `<unk>`);
- `<pad>`, `<unk>`, `<bos>`, `<eos>`, `<mask>` separados;
- sin el prefijo `▁` de SentencePiece;
- motivos reguladores protegidos.

## Restricción combinatoria

Con alfabeto de 4 letras, el número de piezas posibles de longitud ≤ L es
`4^1 + … + 4^L`. Por eso:

| max piece | piezas ACGT posibles |
|---:|---:|
| 3 | 84 |
| 4 | 340 |
| 5 | 1.364 |
| 6 | 5.460 |
| 7 | 21.844 |

**Un vocabulario de 8k con max piece 3 es imposible.** El diseño inicial
(8k/3) se corrigió a un barrido de puntos factibles.

## Barrido entrenado y evaluado

Corpus: 168.373 ventanas de 2.048 pb (345 MB). Evaluación: 6 genomas
(4 no vistos + 2 vistos), 2 Mb por genoma.

| tokenizer | vocab/max | nt/token | tokens 1.386 Gbp | IUPAC | motivos 1 token | `▁` | veredicto |
|---|---:|---:|---:|---:|---:|---:|---|
| v2 | 4.096/6 | 4,00 | 346M | ❌ | 7/7 | 2 | baseline |
| v3a | 1.024/5 | 3,43 | 404M | ✅ | 7/7 | 0 | bueno |
| v3b | 512/5 | 3,28 | 423M | ✅ | 7/7 | 0 | bueno |
| **v3c** | **256/4** | **2,97** | **467M** | ✅ | 7/7 | 0 | **cumple 2–3** |

`motivos 1 token` ignora el `▁` de v2; mide si el motivo se parte o no.

## Tokenizer elegido: v3c 256/4

Archivos (canónico = copia de `_256`):

- `data/tokenizer_v3/dna_sp_v3.model`
- `data/tokenizer_v3/tokenizer_meta_v3.json`
- `data/tokenizer_v3/EVAL.md`
- variantes del barrido: `data/tokenizer_v3_1024/`, `data/tokenizer_v3_512/`, `data/tokenizer_v3_256/`

Métricas:

| métrica | v2 | v3c |
|---|---:|---:|
| vocab | 4.096 | 256 |
| nt/token | 4,004 | **2,969** |
| tokens por ventana de 1 kb | 252,5 | 337,9 |
| tokens para 1,386 Gbp | 346.134.714 | 466.778.782 |
| unk rate en genomas | 0 | 0 |
| roundtrip | 120/120 | 120/120 |
| IUPAC completo | ❌ | ✅ |
| motivos 1 token | 7/7 | 7/7 |
| `▁` | 2 | 0 |
| `<pad>/<unk>/<bos>/<eos>/<mask>` | -1/0/1/2/0 | **0/1/2/3/4** |

### Motivos en contexto

| motivo | v2 | v3c |
|---|---:|---:|
| TATAAT | 1,000 | 1,000 |
| TTGACA | 1,000 | 1,000 |
| AGGAGG | 0,696 | 0,703 |
| ATG | 0,989 | 0,993 |
| TAA | 1,000 | 1,000 |
| TAG | 1,000 | 1,000 |
| TGA | 0,978 | 0,658 |

**Límite de solapamiento:** los símbolos atómicos no resuelven motivos que
comparten bases.

- `TGA` comparte la última base con `ATG`, `GTG` y `TTG`; en `ATGA` gana
  `ATG` y `TGA` queda partido.
- `AGGAGG` contiene `GGAGG`; algunas apariciones se resuelven como
  `A` + `GGAGG`.
- `TATAAT` y `TTGACA` no se solapan con otros motivos protegidos: 1,00.

Para el modelo nuevo: definir una **prioridad de solapamiento** o un
pretokenizador por límites de codón. No basta con `user_defined_symbols`.

## Veredicto

| gate | resultado |
|---|---|
| 2–3 nt/token | ✅ 2,97 |
| IUPAC completo | ✅ |
| 7/7 motivos atómicos | ✅ |
| sin `▁` | ✅ |
| tokens especiales distintos | ✅ |
| roundtrip exacto | ✅ 120/120 |

**v3c 256/4 es el tokenizer recomendado para el próximo modelo.** Si el
vocabulario de 256 resulta limitante, usar v3b 512/5 (3,28 nt/token) como
fallback.

## Reproducir

```bash
cd dna_poc
source .venv/bin/activate

# Barrido
python scripts/06_train_tokenizer_v3.py --vocab-size 1024 --max-piece 5 \
    --out data/tokenizer_v3_1024
python scripts/06_train_tokenizer_v3.py --vocab-size 512 --max-piece 5 \
    --out data/tokenizer_v3_512
python scripts/06_train_tokenizer_v3.py --vocab-size 256 --max-piece 4 \
    --out data/tokenizer_v3_256

# Evaluación contra v2
python scripts/06_eval_tokenizer_v3.py \
    --v3 data/tokenizer_v3_256/dna_sp_v3.model \
    --out data/tokenizer_v3_256
```

## Límites

- Es una evaluación **intrínseca**. El impacto en el modelo requiere entrenar
  un modelo con v3c y comparar firmas, PPL y promotores.
- 256 piezas es un vocabulario pequeño; si el modelo subajusta, 512/5 es el
  fallback.
- La tasa de TGA/AGGAGG en contexto no es 1,0 por el solapamiento descrito.
- El corpus de entrenamiento son 310 genomas de vida libre; la evaluación usa
  4 genomas no vistos (simbiontes) y 2 vistos.
