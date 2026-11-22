# Diseño del próximo modelo de ADN — dictamen basado en 20 genomas

> Documento de trabajo. Resume qué aprende el modelo v2 (epoch 9), qué no
> aprende, y qué cambiar en el próximo entrenamiento. Los números de firmas
> provienen de `data/embeddings/v2_epoch9_20g_256/` (20 genomas, 9.849
> ventanas de 256 pb, 384 dims, ajuste GC + genoma).

## 0. Revisión 2026-11-20: benchmark de modelos pretrained

Antes de entrenar desde cero se evaluaron modelos existentes con las mismas
ventanas y etiquetas. Resultado (600 ventanas, 20 genomas, LOGO por genoma,
ajuste GC+genoma):

| modelo | balanced acc | AUC CDS | AUC tRNA | AUC rRNA | AUC intergénico | promotor |
|---|---:|---:|---:|---:|---:|---:|
| **ProkBERT-mini-long** (27M, 4 kb) | **0,825** | **0,981** | **0,992** | **0,967** | **0,909** | **0,988** |
| NT-v2 50M (12 kb) | 0,510 | 0,818 | 0,840 | 0,663 | 0,626 | 0,908 |
| v2 epoch9 (propio) | 0,376 | 0,669 | 0,708 | 0,701 | 0,528 | 0,767 |
| 3-mer baseline | 0,502 | 0,719 | 0,732 | 0,842 | 0,638 | 0,850* |
| GC baseline | 0,361 | 0,503 | 0,652 | 0,699 | 0,629 | 0,455 |

\* conteo de TATAAT para promotores.

**Consecuencia:** el próximo paso ya no es “entrenar el conv-transformer desde
cero”. Es:

1. **Fine-tunear ProkBERT** para componentes/promotor local (27M, bacterial,
   4 kb); es el mejor por un margen grande.
2. **Usar NT-v2 12 kb** (o un modelo long-context en GPU) para operones.
3. **Destilar** ProkBERT + NT-v2 en un encoder 12–24 kb con tokenizer v3c
   256/4, si se necesita un solo modelo.
4. **V3 desde cero solo como fallback**, si el fine-tuning de ProkBERT no
   generaliza a genomas no vistos.

El benchmark completo, capacidad del equipo y caveats están en
[`PRETRAINED_BENCHMARK.md`](PRETRAINED_BENCHMARK.md).

## 1. Qué se midió

- 20 genomas de 20 filos: 16 de vida libre + 4 simbiontes/parásitos
  (`raw_excluded`).
- 9.849 ventanas de 256 pb: 1.828 CDS, 900 tRNA, 143 rRNA, 2.140 intergénico,
  1.810 upstream, 555 transposón, 211 profago, 179 CRISPR.
- Embeddings de 384 dims del checkpoint `dna_model_v2_epoch9.pt`.
- Firma por componente: centroide, `d` de Cohen vs resto, FDR por ventana,
  meta-análisis por genoma, consistencia de signo. Ajuste por GC/GC² y por
  genoma (efectos fijos).

## 2. Qué aprende hoy (evidencia)

| componente | dims de firma | `|d|` máx | `q_meta` mín | consistencia |
|---|---:|---:|---:|---:|
| tRNA | 249 | 1,00 | 2,5e−10 | 20/20 |
| rRNA | 63 | 1,30 | 1,8e−3 | 16-18/20 |
| CDS | 32 | 0,42 | 7,0e−4 | 16-17/20 |
| transposón | 25 | 0,33 | 3,6e−2 | 16-17/20 |
| intergénico | 0 | — | — | — |
| upstream | 0 | — | — | — |
| profago | 0 | — | — | — |
| CRISPR | 0 | — | — | — |

Lectura:

- **tRNA es el componente mejor codificado**: 249 dims, `|d|≈1`, consistencia
  20/20, `q_meta` 1e−10. Es una señal robusta, no un artefacto de GC.
- **rRNA tiene efecto grande** (`|d|≈1,3`) pero solo 143 ventanas; 63 dims
  sobreviven y la consistencia es 0,8–0,9.
- **CDS tiene una firma débil pero reproducible** (`|d|≈0,4`, 32 dims,
  `q_meta` 1e−3). El modelo distingue CDS del resto, pero con poco margen.
- **Intergénico no tiene firma tras controlar GC y genoma.** Antes del ajuste
  aparecían dims; eran composición/linaje. El modelo no aprendió un “código de
  no-gen”.
- **Upstream/promotor no tiene firma.** Ni con ventana local de 256 pb. Es
  coherente con el AUROC ~0,57 en genomas no vistos.
- **Profago/CRISPR no tienen firma estable.** CRISPR solo aparece en 3 genomas;
  profago tiene 77 ventanas y ningún dim estable.

Esto responde a la pregunta de diseño: el modelo actual **sí** tiene un espacio
latente útil para tRNA/rRNA/CDS, pero **no** para regulación, intergénico ni
elementos móviles. El próximo modelo debe conservar lo primero y atacar lo
segundo.

### 2.1 Estabilidad entre genomas (split-half)

20 particiones aleatorias de 10 vs 10 genomas, misma firma en cada mitad:

| componente | Jaccard medio | corr(d) media | dims estables | lectura |
|---|---:|---:|---:|---|
| tRNA | **0,775** | **0,943** | 263 | firma robusta |
| rRNA | 0,146 | 0,765 | 127 | patrón estable, selección dispersa por pocas ventanas |
| CDS | 0,149 | 0,366 | 52 | firma débil/inestable |
| intergénico | 0,001 | **−0,478** | 0 | no hay firma; las mitades se contradicen |

- **tRNA es la única firma que replica casi intacta del piloto de 6 genomas al
  set de 20** (Jaccard 6g vs 20g = 0,88).
- rRNA tiene un patrón global estable (`corr(d)=0,76`) pero las dims
  seleccionadas cambian entre mitades: hay pocas ventanas de rRNA (143) y la
  selección es inestable. Con más rRNA, la firma se estabilizaría.
- CDS tiene `corr(d)=0,37`: el modelo distingue CDS, pero la firma no es
  reproducible entre linajes. Es el punto más débil de la representación.
- Intergénico tiene correlación negativa: no existe un patrón consistente de
  “no gen”.

### 2.2 Clasificación LOGO (20 genomas)

| métrica | valor |
|---|---:|
| Accuracy | 0,394 |
| Balanced accuracy | 0,407 |
| Azar (4 clases) | 0,250 |
| AUC tRNA | 0,748 |
| AUC rRNA | 0,785 |
| AUC CDS | 0,605 |
| AUC intergénico | 0,474 |

El patrón de tRNA/rRNA generaliza entre genomas; CDS es débil; intergénico está
por debajo del azar. Esto no es un problema de classifier: es el techo del
embedding actual.

### 2.3 Composición

En los 20 genomas, `entropy` tiene `|r parcial|≈0,71` y `gc_skew≈0,44`; el GC
bruto correlaciona 0,785. `is_intergenic` correlaciona 0,286 en bruto pero solo
0,084 tras controlar GC: su “señal” es composición/contexto. `is_upstream` cae
a 0,043. La mayor parte del espacio latente sigue siendo composición, no
anotación.

## 3. Por qué el entrenamiento se estancó/osciló en epoch 9

No hay logs de entrenamiento en el repo, pero el setup explica la oscilación:

1. **Batch efectivo 2** con 32k posiciones: el gradiente por paso tiene muy
   alta varianza; el loss por batch oscila aunque la tendencia sea plana.
2. **LR constante 1e−4, sin warmup, sin scheduler y sin gradient clipping.**
   AdamW con batch 2 es propenso a rebotes; no hay nada que amortigüe.
3. **76 % de padding.** Cada muestra de 32k posiciones tiene ~7,6k tokens
   reales; el loss se calcula solo sobre los enmascarados, así que el número
   efectivo de tokens por paso es minúsculo (~2×7,6k×15 % ≈ 2.300).
4. **MLM con máscara aleatoria distinta cada época**, lo que añade ruido al
   loss de una época a otra.
5. **Máscara de atención densa** de 32k×32k: no reduce cómputo; el paso es
   caro y el gradiente puede volverse inestable.

Detener en epoch 9 con val plana fue razonable, pero la causa no es “el modelo
ya aprendió todo”: es un problema de **optimización + empaquetado de datos**.
El próximo entrenamiento debe cambiar esas tres cosas antes de aumentar
capacidad.

Evidencia de los checkpoints: entre epoch 6 y epoch 9 los pesos cambian
(`max|Δ|≈0,035`, media 0,016) pero el loss sobre una secuencia fija es
prácticamente idéntico (474,6 vs 474,2). El modelo se mueve sin mejorar: es una
meseta con oscilación, no un óptimo. La perplejidad de 474 es peor que un
predictor de frecuencia de tokens (213), así que tampoco es un problema de
“ya está bien”.

## 4. Qué cambiar en el modelo nuevo

### 4.1 Tokenizer

Problemas del actual:

- 4,2 nt/token (max piece 6): colapsa motivos de 4–6 pb; TATAAT/TTGACA son
  tokens únicos, pero muchos k-mers reguladores no.
- `B`, `D`, `H`, `V` van a `<unk>`; solo `R,Y,K,M,S,W,N` están protegidos.
- `mask` = token 1 (`<s>`, BOS) y `pad` = token 0 (`<unk>`), no hay
  `<mask>`/`<pad>` dedicados.
- Prefijo `▁` en cada secuencia, un artefacto de SentencePiece.

Resultado medido (barrido real, 6 genomas, 4 no vistos):

| tokenizer | vocab/max | nt/token | tokens 1,386 Gbp | IUPAC | motivos | ▁ | veredicto |
|---|---:|---:|---:|---:|---:|---:|---|
| v2 | 4.096/6 | 4,00 | 346M | ❌ | 7/7 | 2 | baseline |
| v3a | 1.024/5 | 3,43 | 404M | ✅ | 7/7 | 0 | bueno |
| v3b | 512/5 | 3,28 | 423M | ✅ | 7/7 | 0 | bueno |
| **v3c** | **256/4** | **2,97** | **467M** | ✅ | 7/7 | 0 | **cumple 2–3** |

- **Recomendación: v3c 256/4.** Es el único que cumple 2–3 nt/token y mantiene
  IUPAC completo, motivos atómicos, tokens especiales y sin `▁`.
- Si el vocabulario de 256 resulta limitante en el modelo, usar **512/5**
  (3,28 nt/token) como fallback.
- **8k/3 es imposible**: con alfabeto de 4 letras solo hay
  `4^1+4^2+4^3 = 84` piezas de ≤3 nt. Un vocabulario de 8k requiere
  `max_piece ≥ 7` (4^7 = 16.384).
- **Límite de solapamiento:** proteger motivos que comparten bases no garantiza
  que ambos queden atómicos. Medido en contexto: `TGA` 0,658 (choca con
  `ATG`/`GTG`/`TTG`), `AGGAGG` 0,703 (contiene `GGAGG`), `TATAAT` y `TTGACA`
  1,00. El modelo nuevo necesita una **prioridad de solapamiento** o un
  pretokenizador por límites de codón; no basta `user_defined_symbols`.
- `character_coverage=1.0`, `byte_fallback=False`; proteger los 15 IUPAC y los
  motivos (`TATAAT`, `TTGACA`, `AGGAGG`, `ATG`, `TAA`, `TAG`, `TGA`).
- Tokens especiales separados: `<pad>`, `<mask>`, `<unk>`, `<bos>`, `<eos>`.
- Objetivo: 2–3 nt/token. Con 1,386 Gbp → ~467M tokens para 20–30M parámetros.
- Evaluar el tokenizer por **cobertura de motivos**, no solo por compresión.

### 4.2 Contexto y empaquetado

- El contexto real actual es **32 kb** (32k nucleótidos → 7,6k tokens), no
  100 kb.
- Empaquetar secuencias a **8k tokens** (~24 kb con v3c 256/4) o **16k tokens**
  (~48 kb) sin padding; separar contigs con `<eos>`.
- 8k tokens cubre la mayoría de operones; 16k para los largos. No hace falta
  32k tokens para el 99 % de las regiones.
- Si se usa 16k+, atención local eficiente (block-sparse/FlexAttention), no una
  máscara densa.

### 4.3 Arquitectura

- 20–30M parámetros, `d_model=512`, 8 capas, 8 cabezas, RoPE.
- Stem convolucional 1D (kernel 5–7) antes del transformer para motivos
  locales; o un token-level CNN + transformer.
- **Attention pooling** para el embedding de ventana, en lugar de mean pooling
  (hoy el mean diluye la señal local; las firmas de 256 pb son más fuertes que
  las de 1 kb).
- Cabeza token-level y cabeza window-level; no solo MLM.
- Reverse-complement augmentation y, si se puede, consistencia RC del
  embedding.

### 4.4 Objetivo

- MLM 15 % con **80/10/10** y enmascarado de tokens completos; añadir
  **span masking** de 3–5 tokens para motivos/operones.
- **Multi-task auxiliar** (esto es lo que falta para intergénico/promotor):
  - clasificación CDS / tRNA / rRNA / intergénico / upstream;
  - predicción de hebra;
  - oriC;
  - promotor (cuando haya etiquetas; RegulonDB para E. coli).
- Pérdida contrastiva por componente: misma clase en genomas distintos cerca,
  clases distintas lejos. Con 20 genomas ya se ve que tRNA/rRNA/CDS son
  separables; el contraste debe empujar intergénico/upstream.
- Alternativa más eficiente que MLM: ELECTRA (replaced-token detection), que
  usa todos los tokens y no solo el 15 %.

### 4.5 Optimización

- Batch efectivo 16–32 con gradient accumulation (batch físico 2–4).
- Warmup 5–10 %, cosine decay, pico 3e−4–1e−3, grad clip 1,0, weight decay
  0,01, bf16.
- EMA de pesos; checkpoint por val PPL y por sonda de componentes, no por
  train loss.
- Logging de loss suavizado (EMA) y de val cada N pasos; la oscilación de train
  con batch 2 no debe ser el criterio de parada.
- Early stopping con `balanced accuracy` de componentes y AUROC de promotor.

### 4.6 Datos y evaluación

- Usar los 400 genomas, no solo 310; separar **por genoma/especie**, no por
  contig. El split actual por contig filtra información entre train/val.
- Incluir simbiontes para diversidad, pero hold-out por filo.
- Evaluación obligatoria:
  - firmas por componente y su estabilidad entre mitades de genomas;
  - AUROC de promotor con negativos **barajados por dinucleótidos**;
  - PPL de CDS vs intergénico;
  - sondas lineales por componente con split por genoma;
  - comparación contra baseline GC/3-mer.
- No usar un único AUROC como criterio.

## 5. Gates de éxito para el modelo nuevo

| gate | umbral |
|---|---|
| tRNA/rRNA/CDS: balanced accuracy LOGO | ≥ 0,70 |
| intergénico: balanced accuracy LOGO | ≥ 0,60 |
| AUROC promotor (negativos barajados) | ≥ 0,85 |
| PPL CDS < PPL intergénico | ≥ 10 % |
| Firmas estables entre mitades (Jaccard / corr d) | ≥ 0,6 / ≥ 0,8 |
| Superar baseline GC/3-mer | en todas las tareas |

## 6. Primer experimento recomendado (revisado)

1. **Fine-tunear ProkBERT** (27M, 4 kb): LoRA o full fine-tune, 400 genomas,
   hold-out por género/especie; evaluar componentes/promotor con los gates.
2. **Probar NT-v2 250M y Evo2 1B en Modal A10G** para contexto 12 kb/8 kb;
   comparar con ProkBERT en operones.
3. **Stack o destilación**: ProkBERT para ventanas locales + NT-v2 para
   contexto; o destilar ambos en un encoder 12–24 kb con v3c 256/4.
4. **Solo si lo anterior falla**, entrenar el conv-transformer V3 desde cero
   con v3c 256/4, 8k tokens, atención completa, attention pooling y cabezas
   auxiliares.
5. Medir con los mismos gates y con baselines de posición/GC/3-mer.

## 7. Por qué esto es distinto de “entrenar más”

El epoch 9 no está limitado por épocas: está limitado por **qué ve y cómo se
optimiza**. Con 76 % de padding, batch 2, LR constante y una máscara densa,
más épocas solo repiten la misma oscilación. El cambio de tokenizer, el
empaquetado y las cabezas auxiliares atacan la causa, no el síntoma.
